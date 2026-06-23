from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from openmm import CustomExternalForce, LangevinMiddleIntegrator, Platform, unit
from openmm.app import HBonds, Modeller, NoCutoff, PDBFile, Simulation


DRUGS = ["afatinib", "osimertinib", "furmonertinib"]
VARIANTS = ["WT", "G719S", "L861Q", "S768I"]


@dataclass(frozen=True)
class GridBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float


def run_command(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=str(cwd),
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


def read_box(input_dir: Path) -> GridBox:
    meta = json.loads((input_dir / "docking_box_metadata.json").read_text(encoding="utf-8"))
    center = meta["box_center"]
    size = meta["box_size"]
    return GridBox(center[0], center[1], center[2], size[0], size[1], size[2])


def locate_vina(vina_arg: str | None) -> str:
    if vina_arg:
        return vina_arg
    found = shutil.which("vina")
    if found:
        return found
    raise RuntimeError("AutoDock Vina executable not found. Install vina or pass --vina /path/to/vina.")



def available_openmm_platforms() -> list[str]:
    return [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]


def resolve_openmm_platform(requested: str) -> str:
    available = available_openmm_platforms()
    if requested.lower() != "auto":
        if requested in available:
            return requested
        raise RuntimeError(
            f"Requested OpenMM platform '{requested}' is not available. "
            f"Available platforms: {available}. Use --platform auto to avoid stopping."
        )
    for candidate in ["CUDA", "OpenCL", "CPU", "Reference"]:
        if candidate in available:
            return candidate
    raise RuntimeError(f"No usable OpenMM platform found. Available platforms: {available}")
def prepare_receptor_pdbqt_obabel(input_pdb: Path, output_pdbqt: Path) -> None:
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["obabel", "-ipdb", str(input_pdb), "-opdbqt", "-O", str(output_pdbqt), "--partialcharge", "gasteiger", "-xr"]
    result = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    output_pdbqt.with_suffix(".obabel_receptor.log").write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
    if result.returncode != 0 or not output_pdbqt.exists() or output_pdbqt.stat().st_size == 0:
        raise RuntimeError(f"OpenBabel receptor PDBQT conversion failed for {input_pdb}: {result.stderr}")


def add_backbone_restraints(system, topology, positions, restraint_k: float) -> int:
    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", restraint_k)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    count = 0
    for atom in topology.atoms():
        if atom.name in {"N", "CA", "C", "O"}:
            pos = positions[atom.index].value_in_unit(unit.nanometer)
            force.addParticle(atom.index, pos)
            count += 1
    system.addForce(force)
    return count


def heavy_atom_rmsd(initial_pdb: Path, relaxed_pdb: Path) -> float:
    initial = PDBFile(str(initial_pdb))
    relaxed = PDBFile(str(relaxed_pdb))
    total = 0.0
    count = 0
    for atom_i, atom_r in zip(initial.topology.atoms(), relaxed.topology.atoms()):
        if atom_i.element is None or atom_i.element.symbol == "H":
            continue
        a = initial.positions[atom_i.index].value_in_unit(unit.angstrom)
        b = relaxed.positions[atom_r.index].value_in_unit(unit.angstrom)
        total += sum((a[j] - b[j]) ** 2 for j in range(3))
        count += 1
    return math.sqrt(total / count) if count else math.nan


def relax_receptor(
    input_pdb: Path,
    output_pdb: Path,
    *,
    platform_name: str,
    max_iterations: int,
    md_steps: int,
    restraint_k_kj_mol_nm2: float,
) -> dict[str, float | int | str]:
    from openmm.app import ForceField

    pdb = PDBFile(str(input_pdb))
    forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield, pH=7.4)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    restrained_atoms = add_backbone_restraints(system, modeller.topology, modeller.positions, restraint_k_kj_mol_nm2)
    integrator = LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = Platform.getPlatformByName(platform_name)
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    energy0 = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    simulation.minimizeEnergy(maxIterations=max_iterations)
    if md_steps > 0:
        simulation.step(md_steps)
        simulation.minimizeEnergy(maxIterations=max_iterations)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy1 = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with output_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(simulation.topology, state.getPositions(), handle, keepIds=True)
    return {
        "platform": platform_name,
        "initial_energy_kcal_mol": energy0,
        "relaxed_energy_kcal_mol": energy1,
        "energy_delta_kcal_mol": energy1 - energy0,
        "backbone_restrained_atoms": restrained_atoms,
        "heavy_atom_rmsd_A": heavy_atom_rmsd(input_pdb, output_pdb),
    }


def run_vina(
    vina: str,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    output_pose: Path,
    log_path: Path,
    box: GridBox,
    exhaustiveness: int,
    seed: int,
) -> float:
    output_pose.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        vina,
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        f"{box.center_x:.3f}",
        "--center_y",
        f"{box.center_y:.3f}",
        "--center_z",
        f"{box.center_z:.3f}",
        "--size_x",
        f"{box.size_x:.1f}",
        "--size_y",
        f"{box.size_y:.1f}",
        "--size_z",
        f"{box.size_z:.1f}",
        "--exhaustiveness",
        str(exhaustiveness),
        "--seed",
        str(seed),
        "--out",
        str(output_pose),
    ]
    result = run_command(cmd, cwd=output_pose.parent)
    log_path.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
    return parse_vina_score(log_path)


def parse_vina_score(log_path: Path) -> float:
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            return float(parts[1])
    raise RuntimeError(f"Could not parse Vina score from {log_path}")


def compare_with_unrelaxed(input_dir: Path, relaxed_scores: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    base = pd.read_csv(input_dir / "docking_scores.csv")
    base = base.rename(columns={"vina_score_kcal_mol": "unrelaxed_score_kcal_mol", "rank_within_variant": "unrelaxed_rank"})
    base = base[["variant_model", "drug", "unrelaxed_score_kcal_mol", "unrelaxed_rank"]]
    relaxed = relaxed_scores.rename(columns={"vina_score_kcal_mol": "relaxed_score_kcal_mol", "rank_within_variant": "relaxed_rank"})
    relaxed = relaxed[["variant_model", "drug", "relaxed_score_kcal_mol", "relaxed_rank"]]
    comparison = base.merge(relaxed, on=["variant_model", "drug"], how="inner")
    comparison["score_delta_relaxed_minus_unrelaxed"] = comparison["relaxed_score_kcal_mol"] - comparison["unrelaxed_score_kcal_mol"]
    comparison["rank_changed"] = comparison["relaxed_rank"] != comparison["unrelaxed_rank"]
    comparison.to_csv(output_dir / "openmm_gpu_relaxed_vs_unrelaxed_comparison.csv", index=False)
    return comparison


def plot_outputs(comparison: pd.DataFrame, output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for variant, group in comparison.groupby("variant_model"):
        axes[0].scatter(group["unrelaxed_score_kcal_mol"], group["relaxed_score_kcal_mol"], label=variant, s=45)
    lo = min(comparison["unrelaxed_score_kcal_mol"].min(), comparison["relaxed_score_kcal_mol"].min()) - 0.1
    hi = max(comparison["unrelaxed_score_kcal_mol"].max(), comparison["relaxed_score_kcal_mol"].max()) + 0.1
    axes[0].plot([lo, hi], [lo, hi], color="#6b7280", linewidth=1)
    axes[0].set_xlabel("Unrelaxed Vina score (kcal/mol)")
    axes[0].set_ylabel("OpenMM-relaxed Vina score (kcal/mol)")
    axes[0].set_title("GPU-relaxed receptor sensitivity", fontsize=10)
    axes[0].legend(frameon=False, fontsize=7)
    pivot = comparison.pivot_table(
        index="variant_model",
        columns="drug",
        values="score_delta_relaxed_minus_unrelaxed",
        aggfunc="first",
    ).loc[VARIANTS, DRUGS]
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
        fig.savefig(output_dir / f"openmm_gpu_relaxed_docking_sensitivity.{ext}", dpi=300)
    plt.close(fig)


def write_report(relax_metrics: pd.DataFrame, scores: pd.DataFrame, comparison: pd.DataFrame, output_dir: Path) -> None:
    best = scores.sort_values(["variant_model", "rank_within_variant"]).groupby("variant_model", as_index=False).first()
    lines = [
        "# OpenMM GPU-Relaxed EGFR-PACC Docking",
        "",
        "## Scope",
        "",
        "EGFR kinase-domain receptor models were relaxed by OpenMM restrained minimization on GPU/accelerated runtime, followed by the same AutoDock Vina screen.",
        "This is a receptor-conformation sensitivity analysis, not production MD, FEP, or experimental affinity measurement.",
        "",
        "## Relaxation Metrics",
        "",
        relax_metrics.to_markdown(index=False),
        "",
        "## Best Drug By Relaxed Receptor",
        "",
        best[["variant_model", "drug", "vina_score_kcal_mol", "rank_within_variant"]].to_markdown(index=False),
        "",
        "## Relaxed vs Unrelaxed",
        "",
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
        ].to_markdown(index=False),
        "",
        "## Manuscript Boundary",
        "",
        "- Report as restrained receptor-relaxation sensitivity.",
        "- Do not report as production MD or covalent binding free energy.",
        "- Use rank stability rather than absolute Vina score magnitude as the main interpretation.",
    ]
    (output_dir / "openmm_gpu_relaxed_docking_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("colab_input"))
    parser.add_argument("--output-dir", type=Path, default=Path("colab_output/openmm_gpu_relaxed"))
    parser.add_argument("--vina", default=None)
    parser.add_argument("--platform", default="auto", help="OpenMM platform: auto, CUDA, OpenCL, CPU, or Reference.")
    parser.add_argument("--max-iterations", type=int, default=500)
    parser.add_argument("--md-steps", type=int, default=1000)
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260620)
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    vina = locate_vina(args.vina)
    platform_name = resolve_openmm_platform(args.platform)
    print(f"OpenMM available platforms: {available_openmm_platforms()}")
    print(f"Using OpenMM platform: {platform_name}")
    box = read_box(input_dir)

    relax_rows = []
    score_rows = []
    for variant in VARIANTS:
        receptor_pdb = input_dir / "receptors" / f"EGFR_4LRM_chainA_{variant}.pdb"
        relaxed_pdb = output_dir / "receptors" / f"EGFR_4LRM_chainA_{variant}_openmm_gpu_relaxed.pdb"
        relaxed_pdbqt = output_dir / "receptors" / f"EGFR_4LRM_chainA_{variant}_openmm_gpu_relaxed.pdbqt"
        metrics = relax_receptor(
            receptor_pdb,
            relaxed_pdb,
            platform_name=platform_name,
            max_iterations=args.max_iterations,
            md_steps=args.md_steps,
            restraint_k_kj_mol_nm2=1000.0,
        )
        prepare_receptor_pdbqt_obabel(relaxed_pdb, relaxed_pdbqt)
        relax_rows.append({"variant_model": variant, "input_pdb": str(receptor_pdb), "relaxed_pdb": str(relaxed_pdb), **metrics})
        for drug in DRUGS:
            ligand_pdbqt = input_dir / "ligands" / f"{drug}.pdbqt"
            pose_path = output_dir / "poses" / f"{variant}__{drug}__openmm_gpu_relaxed_vina_pose.pdbqt"
            log_path = output_dir / "logs" / f"{variant}__{drug}__openmm_gpu_relaxed_vina.log"
            score = run_vina(vina, relaxed_pdbqt, ligand_pdbqt, pose_path, log_path, box, args.exhaustiveness, args.seed)
            score_rows.append(
                {
                    "variant_model": variant,
                    "drug": drug,
                    "vina_score_kcal_mol": score,
                    "receptor_model": "openmm_gpu_relaxed",
                    "pose_pdbqt": str(pose_path),
                    "log_file": str(log_path),
                }
            )

    relax_metrics = pd.DataFrame(relax_rows)
    relax_metrics.to_csv(output_dir / "openmm_gpu_relaxed_receptor_metrics.csv", index=False)
    scores = pd.DataFrame(score_rows)
    scores["rank_within_variant"] = scores.groupby("variant_model")["vina_score_kcal_mol"].rank(method="first", ascending=True)
    scores.to_csv(output_dir / "openmm_gpu_relaxed_docking_scores.csv", index=False)
    comparison = compare_with_unrelaxed(input_dir, scores, output_dir)
    plot_outputs(comparison, output_dir)
    write_report(relax_metrics, scores, comparison, output_dir)
    print(f"Wrote GPU-relaxed docking outputs to {output_dir}")


if __name__ == "__main__":
    main()



