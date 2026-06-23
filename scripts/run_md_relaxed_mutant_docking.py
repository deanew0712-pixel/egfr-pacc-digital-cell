from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
from openmm import CustomExternalForce, LangevinMiddleIntegrator, Platform, unit
from openmm.app import HBonds, Modeller, NoCutoff, PDBFile, Simulation

from run_basic_pacc_docking import (
    DRUGS,
    MUTATION_SETS,
    ROOT,
    GridBox,
    prepare_receptor_pdbqt_obabel,
    run_vina,
)


BASE_DOCKING_DIR = ROOT / "reports" / "docking_pacc_core"
OUT_DIR = ROOT / "reports" / "docking_md_relaxed"


def read_box() -> GridBox:
    meta = json.loads((BASE_DOCKING_DIR / "docking_box_metadata.json").read_text(encoding="utf-8"))
    center = meta["box_center"]
    size = meta["box_size"]
    return GridBox(center[0], center[1], center[2], size[0], size[1], size[2])


def parse_vina_best_score(log_path: Path) -> float:
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            return float(parts[1])
    raise RuntimeError(f"Could not parse Vina score from {log_path}")


def backbone_restraint(system, topology, positions, restraint_k_kj_mol_nm2: float) -> int:
    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", restraint_k_kj_mol_nm2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    restrained = 0
    for atom in topology.atoms():
        if atom.name in {"N", "CA", "C", "O"}:
            pos = positions[atom.index].value_in_unit(unit.nanometer)
            force.addParticle(atom.index, pos)
            restrained += 1
    system.addForce(force)
    return restrained


def heavy_atom_rmsd(initial_pdb: Path, relaxed_pdb: Path) -> float:
    initial = PDBFile(str(initial_pdb))
    relaxed = PDBFile(str(relaxed_pdb))
    pairs = []
    for atom_i, atom_r in zip(initial.topology.atoms(), relaxed.topology.atoms()):
        if atom_i.element is not None and atom_i.element.symbol != "H":
            pairs.append((atom_i.index, atom_r.index))
    if not pairs:
        return math.nan
    initial_pos = initial.positions
    relaxed_pos = relaxed.positions
    total = 0.0
    for idx_i, idx_r in pairs:
        a = initial_pos[idx_i].value_in_unit(unit.angstrom)
        b = relaxed_pos[idx_r].value_in_unit(unit.angstrom)
        total += sum((a[j] - b[j]) ** 2 for j in range(3))
    return math.sqrt(total / len(pairs))


def relax_receptor(
    input_pdb: Path,
    relaxed_pdb: Path,
    *,
    max_iterations: int,
    md_steps: int,
    restraint_k_kj_mol_nm2: float,
) -> dict[str, Any]:
    from openmm.app import ForceField

    pdb = PDBFile(str(input_pdb))
    forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield, pH=7.4)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
    )
    restrained_atoms = backbone_restraint(system, modeller.topology, modeller.positions, restraint_k_kj_mol_nm2)
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    platform = Platform.getPlatformByName("CPU")
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    state0 = simulation.context.getState(getEnergy=True)
    energy0 = state0.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    simulation.minimizeEnergy(maxIterations=max_iterations)
    if md_steps > 0:
        simulation.step(md_steps)
        simulation.minimizeEnergy(maxIterations=max_iterations)
    state1 = simulation.context.getState(getEnergy=True, getPositions=True)
    energy1 = state1.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    relaxed_pdb.parent.mkdir(parents=True, exist_ok=True)
    with relaxed_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(simulation.topology, state1.getPositions(), handle, keepIds=True)
    return {
        "initial_energy_kcal_mol": energy0,
        "relaxed_energy_kcal_mol": energy1,
        "energy_delta_kcal_mol": energy1 - energy0,
        "backbone_restrained_atoms": restrained_atoms,
        "heavy_atom_rmsd_A": heavy_atom_rmsd(input_pdb, relaxed_pdb),
    }


def run_relaxed_docking(exhaustiveness: int, seed: int, max_iterations: int, md_steps: int) -> pd.DataFrame:
    box = read_box()
    receptors_dir = OUT_DIR / "receptors"
    poses_dir = OUT_DIR / "poses"
    logs_dir = OUT_DIR / "logs"
    receptors_dir.mkdir(parents=True, exist_ok=True)
    poses_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    relax_rows: list[dict[str, Any]] = []
    score_rows: list[dict[str, Any]] = []
    for variant in MUTATION_SETS:
        base_pdb = BASE_DOCKING_DIR / "receptors" / f"EGFR_4LRM_chainA_{variant}.pdb"
        relaxed_pdb = receptors_dir / f"EGFR_4LRM_chainA_{variant}_openmm_relaxed.pdb"
        relaxed_pdbqt = receptors_dir / f"EGFR_4LRM_chainA_{variant}_openmm_relaxed.pdbqt"
        if relaxed_pdb.exists():
            metrics = {
                "initial_energy_kcal_mol": math.nan,
                "relaxed_energy_kcal_mol": math.nan,
                "energy_delta_kcal_mol": math.nan,
                "backbone_restrained_atoms": math.nan,
                "heavy_atom_rmsd_A": heavy_atom_rmsd(base_pdb, relaxed_pdb),
            }
        else:
            metrics = relax_receptor(
                base_pdb,
                relaxed_pdb,
                max_iterations=max_iterations,
                md_steps=md_steps,
                restraint_k_kj_mol_nm2=1000.0,
            )
        prepare_receptor_pdbqt_obabel(relaxed_pdb, relaxed_pdbqt)
        relax_rows.append(
            {
                "variant_model": variant,
                "input_pdb": str(base_pdb),
                "relaxed_pdb": str(relaxed_pdb),
                "relaxed_pdbqt": str(relaxed_pdbqt),
                "relaxation_method": "OpenMM amber14/GBn2 pH7.4 hydrogens, backbone-restrained minimization"
                + (f" plus {md_steps} restrained MD steps" if md_steps > 0 else ""),
                **metrics,
            }
        )
        for drug in DRUGS:
            ligand_pdbqt = BASE_DOCKING_DIR / "ligands" / f"{drug}.pdbqt"
            pose_path = poses_dir / f"{variant}__{drug}__openmm_relaxed_vina_pose.pdbqt"
            log_path = logs_dir / f"{variant}__{drug}__openmm_relaxed_vina.log"
            if not (pose_path.exists() and log_path.exists() and log_path.stat().st_size > 0):
                run_vina(relaxed_pdbqt, ligand_pdbqt, pose_path, log_path, box, exhaustiveness, seed)
            score_rows.append(
                {
                    "variant_model": variant,
                    "drug": drug,
                    "vina_score_kcal_mol": parse_vina_best_score(log_path),
                    "receptor_model": "openmm_relaxed",
                    "pose_pdbqt": str(pose_path),
                    "log_file": str(log_path),
                }
            )
    relax_df = pd.DataFrame(relax_rows)
    relax_df.to_csv(OUT_DIR / "openmm_relaxed_receptor_metrics.csv", index=False)
    scores = pd.DataFrame(score_rows)
    scores["rank_within_variant"] = scores.groupby("variant_model")["vina_score_kcal_mol"].rank(
        method="first", ascending=True
    )
    scores.to_csv(OUT_DIR / "openmm_relaxed_docking_scores.csv", index=False)
    return scores


def summarize_against_unrelaxed(scores: pd.DataFrame) -> pd.DataFrame:
    base = pd.read_csv(BASE_DOCKING_DIR / "docking_scores.csv")
    base = base.rename(
        columns={
            "vina_score_kcal_mol": "unrelaxed_score_kcal_mol",
            "rank_within_variant": "unrelaxed_rank",
        }
    )[["variant_model", "drug", "unrelaxed_score_kcal_mol", "unrelaxed_rank"]]
    relaxed = scores.rename(
        columns={
            "vina_score_kcal_mol": "relaxed_score_kcal_mol",
            "rank_within_variant": "relaxed_rank",
        }
    )[["variant_model", "drug", "relaxed_score_kcal_mol", "relaxed_rank"]]
    merged = base.merge(relaxed, on=["variant_model", "drug"], how="inner")
    merged["score_delta_relaxed_minus_unrelaxed"] = (
        merged["relaxed_score_kcal_mol"] - merged["unrelaxed_score_kcal_mol"]
    )
    merged["rank_changed"] = merged["relaxed_rank"] != merged["unrelaxed_rank"]
    merged.to_csv(OUT_DIR / "openmm_relaxed_vs_unrelaxed_comparison.csv", index=False)
    return merged


def plot_relaxed_outputs(comparison: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for variant, group in comparison.groupby("variant_model"):
        axes[0].scatter(group["unrelaxed_score_kcal_mol"], group["relaxed_score_kcal_mol"], label=variant, s=45)
    lo = min(comparison["unrelaxed_score_kcal_mol"].min(), comparison["relaxed_score_kcal_mol"].min()) - 0.1
    hi = max(comparison["unrelaxed_score_kcal_mol"].max(), comparison["relaxed_score_kcal_mol"].max()) + 0.1
    axes[0].plot([lo, hi], [lo, hi], color="#6b7280", linewidth=1)
    axes[0].set_xlabel("Unrelaxed Vina score (kcal/mol)")
    axes[0].set_ylabel("OpenMM-relaxed Vina score (kcal/mol)")
    axes[0].set_title("Docking score sensitivity", fontsize=10)
    axes[0].legend(frameon=False, fontsize=7)

    pivot = comparison.pivot_table(
        index="variant_model",
        columns="drug",
        values="score_delta_relaxed_minus_unrelaxed",
        aggfunc="first",
    ).loc[[name for name in MUTATION_SETS if name in comparison["variant_model"].values], DRUGS]
    image = axes[1].imshow(pivot.values, cmap="coolwarm", aspect="auto")
    axes[1].set_xticks(range(len(pivot.columns)), pivot.columns, rotation=35, ha="right")
    axes[1].set_yticks(range(len(pivot.index)), pivot.index)
    axes[1].set_title("Relaxed - unrelaxed score delta", fontsize=10)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            axes[1].text(j, i, f"{pivot.values[i, j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=axes[1], fraction=0.046, pad=0.04, label="kcal/mol")
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT_DIR / f"openmm_relaxed_docking_sensitivity.{ext}", dpi=300)
    plt.close(fig)


def frame_to_markdown(df: pd.DataFrame) -> str:
    display = df.copy()
    for col in display.columns:
        if pd.api.types.is_float_dtype(display[col]):
            display[col] = display[col].map(lambda value: "" if pd.isna(value) else f"{value:.3f}")
    lines = [
        "| " + " | ".join(display.columns.astype(str)) + " |",
        "| " + " | ".join(["---"] * len(display.columns)) + " |",
    ]
    for _, row in display.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in display.columns) + " |")
    return "\n".join(lines)


def write_report(scores: pd.DataFrame, comparison: pd.DataFrame) -> None:
    best = scores.sort_values(["variant_model", "rank_within_variant"]).groupby("variant_model", as_index=False).first()
    rank_changes = comparison[comparison["rank_changed"]].copy()
    lines = [
        "# OpenMM-Relaxed Mutant Model Docking",
        "",
        "## Scope",
        "",
        "WT, G719S, L861Q, and S768I receptor models were relaxed with OpenMM before repeating the same Vina docking screen.",
        "This is a restrained local relaxation sensitivity analysis, not a production MD ensemble or FEP calculation.",
        "",
        "## Best Drug by Relaxed Receptor Model",
        "",
        frame_to_markdown(best[["variant_model", "drug", "vina_score_kcal_mol", "rank_within_variant"]]),
        "",
        "## Relaxed vs Unrelaxed Comparison",
        "",
        frame_to_markdown(
            comparison[
                [
                    "variant_model",
                    "drug",
                    "unrelaxed_score_kcal_mol",
                    "relaxed_score_kcal_mol",
                    "score_delta_relaxed_minus_unrelaxed",
                    "unrelaxed_rank",
                    "relaxed_rank",
                    "rank_changed",
                ]
            ]
        ),
        "",
        "## Interpretation",
        "",
        "- Relaxed receptor docking is used here as a sensitivity check for the coarse structural supplement.",
        "- Rank changes after relaxation should be interpreted as receptor-conformation sensitivity, not definitive drug superiority.",
        "- The preferred manuscript language is 'OpenMM-relaxed receptor sensitivity analysis' rather than 'production MD validation'.",
    ]
    if not rank_changes.empty:
        lines.extend(["", "Rank changes were observed in:", "", frame_to_markdown(rank_changes[["variant_model", "drug", "unrelaxed_rank", "relaxed_rank"]])])
    (OUT_DIR / "openmm_relaxed_docking_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260620)
    parser.add_argument("--max-iterations", type=int, default=500)
    parser.add_argument("--md-steps", type=int, default=0)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    scores = run_relaxed_docking(args.exhaustiveness, args.seed, args.max_iterations, args.md_steps)
    comparison = summarize_against_unrelaxed(scores)
    plot_relaxed_outputs(comparison)
    write_report(scores, comparison)
    print(f"Wrote OpenMM-relaxed docking outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
