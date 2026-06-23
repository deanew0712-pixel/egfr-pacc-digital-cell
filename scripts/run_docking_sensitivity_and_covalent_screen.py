from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem

from run_basic_pacc_docking import (
    DRUGS,
    MUTATION_SETS,
    OBABEL_EXE,
    ROOT,
    TEMPLATE_CHAIN,
    VINA_EXE,
    GridBox,
    prepare_receptor_pdbqt_obabel,
    run_vina,
)


BASE_DOCKING_DIR = ROOT / "reports" / "docking_pacc_core"
OUT_DIR = ROOT / "reports" / "docking_workflow_sensitivity"


@dataclass(frozen=True)
class ReceptorWorkflow:
    name: str
    description: str
    ph: float | None = None
    gasteiger: bool = False


WORKFLOWS = [
    ReceptorWorkflow(
        "obabel_default",
        "OpenBabel PDBQT conversion from heavy-atom PDB; default OpenBabel charges.",
    ),
    ReceptorWorkflow(
        "obabel_gasteiger",
        "OpenBabel PDBQT conversion with Gasteiger partial charges.",
        gasteiger=True,
    ),
    ReceptorWorkflow(
        "pdbfixer_ph7p4_obabel_gasteiger",
        "PDBFixer protonation at pH 7.4, nonpolar hydrogens removed, followed by OpenBabel Gasteiger PDBQT.",
        ph=7.4,
        gasteiger=True,
    ),
]


def read_box() -> GridBox:
    meta = json.loads((BASE_DOCKING_DIR / "docking_box_metadata.json").read_text(encoding="utf-8"))
    center = meta["box_center"]
    size = meta["box_size"]
    return GridBox(center[0], center[1], center[2], size[0], size[1], size[2])


def add_hydrogens_with_pdbfixer(input_pdb: Path, out_pdb: Path, ph: float) -> None:
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(input_pdb))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=ph)
    with out_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)


def obabel_prepare_receptor(
    input_pdb: Path,
    out_pdbqt: Path,
    *,
    preserve_hydrogens: bool = False,
    polar_hydrogens_only: bool = False,
    gasteiger: bool = False,
) -> None:
    cmd = [
        str(OBABEL_EXE),
        "-ipdb",
        str(input_pdb),
        "-opdbqt",
        "-O",
        str(out_pdbqt),
    ]
    if preserve_hydrogens:
        cmd.append("-xh")
    if polar_hydrogens_only:
        cmd.append("--DelNonPolarH")
    if gasteiger:
        cmd.extend(["--partialcharge", "gasteiger"])
    cmd.append("-xr")
    result = subprocess.run(
        cmd,
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    log_path = out_pdbqt.with_suffix(".obabel_receptor.log")
    log_path.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
    if result.returncode != 0 or not out_pdbqt.exists() or out_pdbqt.stat().st_size == 0:
        raise RuntimeError(f"OpenBabel receptor preparation failed for {input_pdb}: {result.stderr}")


def prepare_workflow_receptors() -> dict[tuple[str, str], Path]:
    receptors: dict[tuple[str, str], Path] = {}
    receptor_dir = OUT_DIR / "receptors"
    receptor_dir.mkdir(parents=True, exist_ok=True)
    for variant in MUTATION_SETS:
        base_pdb = BASE_DOCKING_DIR / "receptors" / f"EGFR_4LRM_chainA_{variant}.pdb"
        for workflow in WORKFLOWS:
            workflow_dir = receptor_dir / workflow.name
            workflow_dir.mkdir(parents=True, exist_ok=True)
            input_pdb = base_pdb
            if workflow.ph is not None:
                input_pdb = workflow_dir / f"EGFR_4LRM_chainA_{variant}_{workflow.name}.pdb"
                add_hydrogens_with_pdbfixer(base_pdb, input_pdb, workflow.ph)
            out_pdbqt = workflow_dir / f"EGFR_4LRM_chainA_{variant}_{workflow.name}.pdbqt"
            if workflow.name == "obabel_default":
                prepare_receptor_pdbqt_obabel(input_pdb, out_pdbqt)
            else:
                obabel_prepare_receptor(
                    input_pdb,
                    out_pdbqt,
                    preserve_hydrogens=workflow.ph is not None,
                    polar_hydrogens_only=workflow.ph is not None,
                    gasteiger=workflow.gasteiger,
                )
            receptors[(variant, workflow.name)] = out_pdbqt
    return receptors


def parse_vina_best_score(log_path: Path) -> float:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        match = re.match(r"\s*1\s+(-?\d+(?:\.\d+)?)\s+", line)
        if match:
            return float(match.group(1))
    for line in text.splitlines():
        if "REMARK VINA RESULT:" in line:
            return float(line.split()[3])
    raise RuntimeError(f"Could not parse Vina score from {log_path}")


def run_workflow_sensitivity(exhaustiveness: int, seed: int) -> pd.DataFrame:
    box = read_box()
    receptors = prepare_workflow_receptors()
    rows: list[dict[str, Any]] = []
    poses_dir = OUT_DIR / "poses"
    logs_dir = OUT_DIR / "logs"
    poses_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    for workflow in WORKFLOWS:
        for variant in MUTATION_SETS:
            receptor_pdbqt = receptors[(variant, workflow.name)]
            for drug in DRUGS:
                ligand_pdbqt = BASE_DOCKING_DIR / "ligands" / f"{drug}.pdbqt"
                pose_path = poses_dir / f"{workflow.name}__{variant}__{drug}__vina_pose.pdbqt"
                log_path = logs_dir / f"{workflow.name}__{variant}__{drug}__vina.log"
                if not (pose_path.exists() and log_path.exists() and log_path.stat().st_size > 0):
                    run_vina(receptor_pdbqt, ligand_pdbqt, pose_path, log_path, box, exhaustiveness, seed)
                rows.append(
                    {
                        "workflow": workflow.name,
                        "workflow_description": workflow.description,
                        "variant_model": variant,
                        "drug": drug,
                        "vina_score_kcal_mol": parse_vina_best_score(log_path),
                        "pose_pdbqt": str(pose_path),
                        "log_file": str(log_path),
                    }
                )
    scores = pd.DataFrame(rows)
    scores["rank_within_variant_workflow"] = scores.groupby(["workflow", "variant_model"])[
        "vina_score_kcal_mol"
    ].rank(method="first", ascending=True)
    scores.to_csv(OUT_DIR / "receptor_workflow_sensitivity_scores.csv", index=False)
    return scores


def summarize_workflow_sensitivity(scores: pd.DataFrame) -> pd.DataFrame:
    summary = (
        scores.groupby(["variant_model", "drug"])
        .agg(
            mean_score=("vina_score_kcal_mol", "mean"),
            sd_score=("vina_score_kcal_mol", "std"),
            min_score=("vina_score_kcal_mol", "min"),
            max_score=("vina_score_kcal_mol", "max"),
            median_rank=("rank_within_variant_workflow", "median"),
            best_rank_count=("rank_within_variant_workflow", lambda x: int((x == 1).sum())),
        )
        .reset_index()
    )
    summary["score_range_kcal_mol"] = summary["max_score"] - summary["min_score"]
    summary.to_csv(OUT_DIR / "receptor_workflow_sensitivity_summary.csv", index=False)
    return summary


def parse_pdbqt_models(path: Path) -> list[dict[str, Any]]:
    models: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.replace("\x00", "")
        if line.startswith("MODEL"):
            current = {"model": int(line.split()[1]), "atoms": {}, "score": math.nan, "smiles_idx_pairs": []}
        elif line.startswith("ENDMDL"):
            if current is not None:
                models.append(current)
            current = None
        elif current is not None and line.startswith("REMARK VINA RESULT:"):
            current["score"] = float(line.split()[3])
        elif current is not None and line.startswith("REMARK SMILES IDX"):
            nums = [int(x) for x in line.split()[3:]]
            current["smiles_idx_pairs"].extend(zip(nums[0::2], nums[1::2]))
        elif current is not None and line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            current["atoms"][serial] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
    if not models:
        models.append({"model": 1, "atoms": {}, "score": math.nan, "smiles_idx_pairs": []})
    return models


def atom_distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def atom_angle(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> float:
    ba = [a[i] - b[i] for i in range(3)]
    bc = [c[i] - b[i] for i in range(3)]
    dot = sum(ba[i] * bc[i] for i in range(3))
    norm = math.sqrt(sum(x * x for x in ba)) * math.sqrt(sum(x * x for x in bc))
    if norm == 0:
        return math.nan
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / norm))))


def find_covalent_cys_sg(receptor_pdb: Path) -> dict[str, Any]:
    cys_atoms = []
    box = read_box()
    center = (box.center_x, box.center_y, box.center_z)
    for line in receptor_pdb.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM  "):
            continue
        if line[17:20].strip() == "CYS" and line[12:16].strip() == "SG":
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            cys_atoms.append(
                {
                    "chain": line[21].strip() or TEMPLATE_CHAIN,
                    "residue_number": int(line[22:26]),
                    "atom_name": "SG",
                    "x": xyz[0],
                    "y": xyz[1],
                    "z": xyz[2],
                    "distance_to_box_center": atom_distance(xyz, center),
                }
            )
    if not cys_atoms:
        raise RuntimeError(f"No CYS SG atom found in {receptor_pdb}")
    return min(cys_atoms, key=lambda row: row["distance_to_box_center"])


def smiles_to_pdbqt_atom_map(model: dict[str, Any]) -> dict[int, int]:
    return {smiles_idx: pdbqt_serial for smiles_idx, pdbqt_serial in model["smiles_idx_pairs"]}


def warhead_smiles_indices(drug: str) -> dict[str, int]:
    ligands = pd.read_csv(ROOT / "data" / "processed" / "chembl_ligands.csv")
    smiles = str(ligands[ligands["drug"].eq(drug)].iloc[0]["canonical_smiles"])
    mol = Chem.MolFromSmiles(smiles)
    mol = max(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True), key=lambda frag: frag.GetNumHeavyAtoms())
    match = mol.GetSubstructMatch(Chem.MolFromSmarts("C=CC(=O)N"))
    if not match:
        raise RuntimeError(f"No acrylamide warhead found for {drug}")
    beta_c, alpha_c, carbonyl_c, oxygen, nitrogen = [idx + 1 for idx in match]
    return {
        "beta_c_smiles_idx": beta_c,
        "alpha_c_smiles_idx": alpha_c,
        "carbonyl_c_smiles_idx": carbonyl_c,
        "oxygen_smiles_idx": oxygen,
        "nitrogen_smiles_idx": nitrogen,
    }


def run_covalent_geometry_screen() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for variant in MUTATION_SETS:
        receptor_pdb = BASE_DOCKING_DIR / "receptors" / f"EGFR_4LRM_chainA_{variant}.pdb"
        cys = find_covalent_cys_sg(receptor_pdb)
        sg = (cys["x"], cys["y"], cys["z"])
        for drug in DRUGS:
            warhead = warhead_smiles_indices(drug)
            pose_path = BASE_DOCKING_DIR / "poses" / f"{variant}__{drug}__vina_pose.pdbqt"
            for model in parse_pdbqt_models(pose_path):
                atom_map = smiles_to_pdbqt_atom_map(model)
                beta_serial = atom_map.get(warhead["beta_c_smiles_idx"])
                alpha_serial = atom_map.get(warhead["alpha_c_smiles_idx"])
                if beta_serial is None or alpha_serial is None:
                    continue
                beta = model["atoms"].get(beta_serial)
                alpha = model["atoms"].get(alpha_serial)
                if beta is None or alpha is None:
                    continue
                distance = atom_distance(sg, beta)
                angle = atom_angle(sg, beta, alpha)
                geometry_penalty = abs(distance - 3.5) / 1.5 + abs(angle - 107.0) / 45.0
                rows.append(
                    {
                        "variant_model": variant,
                        "drug": drug,
                        "pose_model": model["model"],
                        "vina_score_kcal_mol": model["score"],
                        "template_covalent_cys": f"{cys['chain']}:CYS{cys['residue_number']}:SG",
                        "cys_sg_to_acrylamide_beta_c_distance_A": distance,
                        "sg_beta_alpha_angle_deg": angle,
                        "geometry_penalty_lower_is_better": geometry_penalty,
                        "near_attack_pose": bool(distance <= 5.0 and 70.0 <= angle <= 150.0),
                        "pose_pdbqt": str(pose_path),
                    }
                )
    df = pd.DataFrame(rows)
    df = df.sort_values(
        ["variant_model", "drug", "near_attack_pose", "geometry_penalty_lower_is_better", "vina_score_kcal_mol"],
        ascending=[True, True, False, True, True],
    )
    df.to_csv(OUT_DIR / "covalent_warhead_geometry_screen.csv", index=False)
    best = df.groupby(["variant_model", "drug"], as_index=False).first()
    best.to_csv(OUT_DIR / "covalent_warhead_geometry_best_poses.csv", index=False)
    return df


def plot_outputs(scores: pd.DataFrame, covalent: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    pivot = scores.pivot_table(
        index=["variant_model", "drug"],
        columns="workflow",
        values="vina_score_kcal_mol",
        aggfunc="first",
    )
    axes[0].boxplot([pivot[col].dropna() for col in pivot.columns], tick_labels=pivot.columns, showfliers=False)
    axes[0].set_ylabel("Vina score (kcal/mol)")
    axes[0].set_title("Receptor preparation sensitivity", fontsize=10)
    axes[0].tick_params(axis="x", rotation=45, labelsize=7)

    best = covalent.groupby(["variant_model", "drug"], as_index=False).first()
    for drug, group in best.groupby("drug"):
        axes[1].scatter(
            group["cys_sg_to_acrylamide_beta_c_distance_A"],
            group["sg_beta_alpha_angle_deg"],
            label=drug,
            s=38,
        )
    axes[1].axvspan(0, 5.0, color="#e0f2fe", alpha=0.45)
    axes[1].axhspan(70, 150, color="#dcfce7", alpha=0.25)
    axes[1].set_xlabel("Cys SG to acrylamide beta-C distance (A)")
    axes[1].set_ylabel("SG-betaC-alphaC angle (deg)")
    axes[1].set_title("Pre-covalent warhead geometry", fontsize=10)
    axes[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT_DIR / f"docking_workflow_and_covalent_sensitivity.{ext}", dpi=300)
    plt.close(fig)


def frame_to_markdown(df: pd.DataFrame) -> str:
    display = df.copy()
    for col in display.columns:
        if pd.api.types.is_float_dtype(display[col]):
            display[col] = display[col].map(lambda value: "" if pd.isna(value) else f"{value:.3f}")
    headers = [str(col) for col in display.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in display.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in display.columns) + " |")
    return "\n".join(lines)


def write_report(scores: pd.DataFrame, summary: pd.DataFrame, covalent: pd.DataFrame, exhaustiveness: int) -> None:
    workflow_lines = [
        "# Receptor Workflow Sensitivity and Covalent Geometry Screen",
        "",
        "## Scope",
        "",
        "This analysis extends the basic EGFR-PACC docking screen with receptor preparation sensitivity and a pre-covalent warhead-geometry screen.",
        "It is intended as supplementary structural evidence, not as a formal covalent binding-energy calculation.",
        "",
        "## Receptor Workflow Sensitivity",
        "",
        f"AutoDock Vina was rerun with exhaustiveness {exhaustiveness} across {len(WORKFLOWS)} receptor-preparation workflows:",
        "",
    ]
    for workflow in WORKFLOWS:
        workflow_lines.append(f"- `{workflow.name}`: {workflow.description}")
    top_rank = (
        scores[scores["rank_within_variant_workflow"].eq(1)]
        .groupby(["variant_model", "drug"])
        .size()
        .reset_index(name=f"top_rank_count_across_{len(WORKFLOWS)}_workflows")
        .sort_values(["variant_model", f"top_rank_count_across_{len(WORKFLOWS)}_workflows"], ascending=[True, False])
    )
    workflow_lines.extend(
        [
            "",
            f"Top-rank counts across {len(WORKFLOWS)} workflows:",
            "",
            frame_to_markdown(top_rank),
            "",
            "Score-range summary:",
            "",
            frame_to_markdown(summary),
            "",
            "## Covalent Geometry Screen",
            "",
            "The acrylamide warhead was detected from each ligand SMILES using the `C=CC(=O)N` SMARTS pattern.",
            "For each Vina pose, the distance from the ATP-pocket cysteine sulfur to the acrylamide beta-carbon and the SG-betaC-alphaC angle were calculated.",
            "The ATP-pocket cysteine in the 4LRM-derived template is reported by the structure as CYS A800; this is the local template residue used for the geometry screen.",
            "",
            "Best geometry-ranked poses by receptor-drug pair:",
            "",
            frame_to_markdown(
                covalent.groupby(["variant_model", "drug"], as_index=False)
            .first()[
                [
                    "variant_model",
                    "drug",
                    "pose_model",
                    "vina_score_kcal_mol",
                    "template_covalent_cys",
                    "cys_sg_to_acrylamide_beta_c_distance_A",
                    "sg_beta_alpha_angle_deg",
                    "near_attack_pose",
                ]
                ]
            ),
            "",
            "## Interpretation Boundary",
            "",
            "- Receptor protonation/workflow sensitivity is a robustness check for the non-covalent Vina score ranking.",
            "- The covalent analysis is a pre-reaction geometry screen, not a tethered covalent AutoDock4 or FEP calculation.",
            "- A strict covalent docking upgrade should use a validated flexible-side-chain/tethered protocol and should benchmark reaction-state ligand preparation.",
        ]
    )
    (OUT_DIR / "workflow_sensitivity_and_covalent_screen_report.md").write_text(
        "\n".join(workflow_lines) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260620)
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    scores = run_workflow_sensitivity(args.exhaustiveness, args.seed)
    summary = summarize_workflow_sensitivity(scores)
    covalent = run_covalent_geometry_screen()
    plot_outputs(scores, covalent)
    write_report(scores, summary, covalent, args.exhaustiveness)
    print(f"Wrote receptor workflow sensitivity and covalent geometry outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
