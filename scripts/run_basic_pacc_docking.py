from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.PDB import MMCIFParser, PDBIO, Select
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]


PYTHON_RUNTIME = Path(
    r"C:\Users\Administrator\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"
)
SCRIPT_DIR = PYTHON_RUNTIME.parent / "Scripts"
MEEKO_RECEPTOR = SCRIPT_DIR / "mk_prepare_receptor.exe"
MEEKO_LIGAND = SCRIPT_DIR / "mk_prepare_ligand.exe"
OBABEL_EXE = SCRIPT_DIR / "obabel.exe"
VINA_EXE = ROOT / "tools" / "vina" / "vina_1.2.7_win.exe"

TEMPLATE_CIF = ROOT / "data" / "raw" / "pdb" / "4LRM.cif"
TEMPLATE_CHAIN = "A"
TEMPLATE_LIGAND = "YUN"

DRUGS = ["afatinib", "osimertinib", "furmonertinib"]
MUTATION_SETS = {
    "WT": [],
    "G719S": ["GLY-719-SER"],
    "L861Q": ["LEU-861-GLN"],
    "S768I": ["SER-768-ILE"],
}


class ChainProteinSelect(Select):
    def accept_chain(self, chain):
        return 1 if chain.id == TEMPLATE_CHAIN else 0

    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0


@dataclass(frozen=True)
class GridBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float


def run_command(cmd: list[str], cwd: Path = ROOT) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=str(cwd),
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


def extract_template_protein_and_box(work_dir: Path) -> tuple[Path, GridBox, dict[str, Any]]:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("4LRM", TEMPLATE_CIF)
    protein_pdb = work_dir / "4LRM_chainA_protein_only.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(protein_pdb), ChainProteinSelect())

    ligand_coords = []
    for residue in structure[0][TEMPLATE_CHAIN]:
        if residue.id[0].strip() and residue.resname == TEMPLATE_LIGAND:
            ligand_coords.extend([atom.coord for atom in residue])
    if not ligand_coords:
        raise RuntimeError(f"Could not find ligand {TEMPLATE_LIGAND} in chain {TEMPLATE_CHAIN} of {TEMPLATE_CIF}.")
    xs = [float(coord[0]) for coord in ligand_coords]
    ys = [float(coord[1]) for coord in ligand_coords]
    zs = [float(coord[2]) for coord in ligand_coords]
    center = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
    # A 24 A cube is intentionally broad for cross-drug, coarse supplementary docking.
    box = GridBox(center[0], center[1], center[2], 24.0, 24.0, 24.0)
    meta = {
        "template_pdb": "4LRM",
        "chain": TEMPLATE_CHAIN,
        "box_reference_ligand": TEMPLATE_LIGAND,
        "box_center": [round(value, 3) for value in center],
        "box_size": [box.size_x, box.size_y, box.size_z],
        "note": "Search box centered on 4LRM chain A co-crystal ligand YUN.",
    }
    return protein_pdb, box, meta


def pdbfixer_prepare(input_pdb: Path, out_pdb: Path, mutations: list[str]) -> None:
    # PDBFixer is used to add missing side-chain atoms after point mutation.
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(input_pdb))
    if mutations:
        fixer.applyMutations(mutations, TEMPLATE_CHAIN)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    with out_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)


def prepare_receptor_pdbqt(pdb_path: Path, pdbqt_path: Path, box: GridBox) -> None:
    if OBABEL_EXE.exists():
        prepare_receptor_pdbqt_obabel(pdb_path, pdbqt_path)
        return
    prepare_receptor_pdbqt_meeko(pdb_path, pdbqt_path, box)


def prepare_receptor_pdbqt_obabel(pdb_path: Path, pdbqt_path: Path) -> None:
    cmd = [
        str(OBABEL_EXE),
        "-ipdb",
        str(pdb_path),
        "-opdbqt",
        "-O",
        str(pdbqt_path),
        "-xr",
    ]
    result = subprocess.run(
        cmd,
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    log_path = pdbqt_path.with_suffix(".obabel_receptor.log")
    log_path.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
    if not pdbqt_path.exists() or pdbqt_path.stat().st_size == 0:
        raise RuntimeError(f"OpenBabel receptor PDBQT preparation failed for {pdb_path}: {result.stderr}")


def prepare_receptor_pdbqt_meeko(pdb_path: Path, pdbqt_path: Path, box: GridBox) -> None:
    cmd = [
        str(MEEKO_RECEPTOR),
        "--read_pdb",
        str(pdb_path),
        "-o",
        str(pdbqt_path.with_suffix("")),
        "-p",
        str(pdbqt_path),
        "-a",
        "--box_center",
        f"{box.center_x:.3f}",
        f"{box.center_y:.3f}",
        f"{box.center_z:.3f}",
        "--box_size",
        f"{box.size_x:.1f}",
        f"{box.size_y:.1f}",
        f"{box.size_z:.1f}",
    ]
    result = subprocess.run(
        cmd,
        cwd=str(ROOT),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    log_path = pdbqt_path.with_suffix(".meeko_receptor.log")
    log_path.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
    if result.returncode == 0 and pdbqt_path.exists():
        return
    if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
        return
    write_simple_receptor_pdbqt(pdb_path, pdbqt_path)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(
            "\nFALLBACK: Meeko receptor preparation failed; wrote simplified heavy-atom PDBQT with zero charges.\n"
        )


def autodock_atom_type(element: str, atom_name: str) -> str:
    element = (element or "").strip().upper()
    atom_name = atom_name.strip().upper()
    if not element:
        element = re.sub(r"[^A-Z]", "", atom_name)[:1]
    mapping = {
        "C": "C",
        "N": "N",
        "O": "O",
        "S": "S",
        "P": "P",
        "F": "F",
        "CL": "Cl",
        "BR": "Br",
        "I": "I",
        "MG": "Mg",
        "ZN": "Zn",
        "FE": "Fe",
        "MN": "Mn",
        "CA": "Ca",
    }
    return mapping.get(element, element[:1] if element else "C")


def write_simple_receptor_pdbqt(pdb_path: Path, pdbqt_path: Path) -> None:
    lines = [
        "REMARK Simplified receptor PDBQT generated because Meeko receptor preparation failed.",
        "REMARK Heavy atoms only; zero charges; basic AutoDock atom types.",
    ]
    serial = 1
    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_name = line[12:16]
        element = line[76:78].strip()
        if (element or atom_name.strip()[:1]).upper().startswith("H"):
            continue
        ad_type = autodock_atom_type(element, atom_name)
        # PDBQT columns: PDB atom record + partial charge + AutoDock atom type.
        base = f"{line[:6]}{serial:5d}{line[11:66]}"
        lines.append(f"{base:66s}{0.0:10.3f} {ad_type:>2s}")
        serial += 1
    lines.append("END")
    pdbqt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return mol
    return max(frags, key=lambda frag: frag.GetNumHeavyAtoms())


def prepare_ligand_sdf_and_pdbqt(drug: str, out_dir: Path, seed: int = 20260620) -> tuple[Path, Path]:
    ligands = pd.read_csv(ROOT / "data" / "processed" / "chembl_ligands.csv")
    row = ligands[ligands["drug"].eq(drug)].iloc[0]
    smiles = str(row["canonical_smiles"])
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError(f"Could not parse SMILES for {drug}: {smiles}")
    mol = largest_fragment(mol)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        raise RuntimeError(f"RDKit embedding failed for {drug} with status {status}")
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    mol.SetProp("_Name", drug)
    sdf_path = out_dir / f"{drug}_desalted_3d.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()

    pdbqt_path = out_dir / f"{drug}.pdbqt"
    run_command([str(MEEKO_LIGAND), "-i", str(sdf_path), "-o", str(pdbqt_path)])
    return sdf_path, pdbqt_path


def parse_vina_score(stdout: str, pose_path: Path) -> float:
    for pattern in [r"^\s*1\s+(-?\d+\.\d+)\s+", r"REMARK VINA RESULT:\s+(-?\d+\.\d+)"]:
        for line in stdout.splitlines():
            match = re.search(pattern, line)
            if match:
                return float(match.group(1))
    if pose_path.exists():
        for line in pose_path.read_text(encoding="utf-8", errors="ignore").splitlines():
            match = re.search(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", line)
            if match:
                return float(match.group(1))
    return math.nan


def run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    out_pose: Path,
    log_path: Path,
    box: GridBox,
    exhaustiveness: int,
    seed: int,
) -> float:
    cmd = [
        str(VINA_EXE),
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
        "--num_modes",
        "9",
        "--energy_range",
        "4",
        "--seed",
        str(seed),
        "--out",
        str(out_pose),
    ]
    result = run_command(cmd)
    log_path.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
    return parse_vina_score(result.stdout, out_pose)


def make_heatmap(scores: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib.pyplot as plt

    pivot = scores.pivot(index="variant_model", columns="drug", values="vina_score_kcal_mol")
    pivot = pivot.loc[[name for name in MUTATION_SETS if name in pivot.index], DRUGS]
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    im = ax.imshow(pivot.values, cmap="viridis_r", aspect="auto")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_title("Basic Vina docking score matrix")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = pivot.iloc[i, j]
            ax.text(j, i, f"{value:.1f}", ha="center", va="center", color="white" if value < -7 else "black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Vina score (kcal/mol)")
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out_dir / f"basic_pacc_docking_heatmap.{ext}", dpi=300)
    plt.close(fig)


def parse_pdbqt_first_model(path: Path) -> pd.DataFrame:
    rows = []
    in_model = False
    model_seen = False
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("MODEL"):
            if model_seen:
                break
            in_model = True
            model_seen = True
            continue
        if line.startswith("ENDMDL") and in_model:
            break
        if line.startswith(("ATOM", "HETATM")) and (in_model or not model_seen):
            try:
                rows.append(
                    {
                        "atom_name": line[12:16].strip(),
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                        "atom_type": line.split()[-1],
                    }
                )
            except Exception:
                continue
    return pd.DataFrame(rows)


def parse_receptor_pdb_atoms(path: Path) -> pd.DataFrame:
    rows = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        element = line[76:78].strip() or re.sub(r"[^A-Za-z]", "", line[12:16]).strip()[:1]
        if element.upper().startswith("H"):
            continue
        try:
            rows.append(
                {
                    "residue": f"{line[17:20].strip()}{line[22:26].strip()}",
                    "atom_name": line[12:16].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element": element,
                }
            )
        except Exception:
            continue
    return pd.DataFrame(rows)


def atoms_within(receptor: pd.DataFrame, ligand: pd.DataFrame, cutoff: float = 6.0) -> pd.DataFrame:
    if receptor.empty or ligand.empty:
        return receptor.iloc[0:0].copy()
    keep = []
    lig = ligand[["x", "y", "z"]].to_numpy()
    for idx, row in receptor.iterrows():
        xyz = row[["x", "y", "z"]].to_numpy(dtype=float)
        d2 = ((lig - xyz) ** 2).sum(axis=1)
        if float(d2.min()) <= cutoff**2:
            keep.append(idx)
    return receptor.loc[keep].copy()


def make_pose_overview(scores: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib.pyplot as plt

    best = scores.sort_values(["variant_model", "vina_score_kcal_mol"]).groupby("variant_model").head(1)
    best = best.set_index("variant_model").loc[[name for name in MUTATION_SETS if name in best["variant_model"].values]].reset_index()
    if best.empty:
        return
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    colors = {
        "afatinib": "#2563EB",
        "osimertinib": "#0F766E",
        "furmonertinib": "#B45309",
    }
    fig = plt.figure(figsize=(7.0, 6.0))
    for panel, (_, row) in enumerate(best.iterrows(), start=1):
        variant = row["variant_model"]
        drug = row["drug"]
        receptor_pdb = out_dir / "receptors" / f"EGFR_4LRM_chainA_{variant}.pdb"
        pose_pdbqt = Path(row["pose_pdbqt"])
        receptor = parse_receptor_pdb_atoms(receptor_pdb)
        ligand = parse_pdbqt_first_model(pose_pdbqt)
        pocket = atoms_within(receptor, ligand, cutoff=6.0)
        ax = fig.add_subplot(2, 2, panel, projection="3d")
        if not pocket.empty:
            ax.scatter(pocket["x"], pocket["y"], pocket["z"], s=5, c="#9CA3AF", alpha=0.22, depthshade=False)
        if not ligand.empty:
            ax.scatter(
                ligand["x"],
                ligand["y"],
                ligand["z"],
                s=22,
                c=colors.get(drug, "#111827"),
                alpha=0.95,
                depthshade=True,
            )
            # Draw approximate ligand bonds by distance to improve pose readability.
            xyz = ligand[["x", "y", "z"]].to_numpy(dtype=float)
            for i in range(len(xyz)):
                for j in range(i + 1, len(xyz)):
                    dist = float(((xyz[i] - xyz[j]) ** 2).sum() ** 0.5)
                    if 0.9 <= dist <= 1.9:
                        ax.plot(
                            [xyz[i, 0], xyz[j, 0]],
                            [xyz[i, 1], xyz[j, 1]],
                            [xyz[i, 2], xyz[j, 2]],
                            color=colors.get(drug, "#111827"),
                            linewidth=0.9,
                            alpha=0.75,
                        )
            center = xyz.mean(axis=0)
            radius = 8.0
            ax.set_xlim(center[0] - radius, center[0] + radius)
            ax.set_ylim(center[1] - radius, center[1] + radius)
            ax.set_zlim(center[2] - radius, center[2] + radius)
        ax.set_title(f"{variant}: {drug}, {row['vina_score_kcal_mol']:.2f} kcal/mol")
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.view_init(elev=22, azim=-58)
    fig.suptitle("Representative best-scoring Vina poses in EGFR ATP pocket", y=0.98)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out_dir / f"basic_pacc_docking_pose_overview.{ext}", dpi=300)
    plt.close(fig)


def markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "None."
    display = frame.copy()
    for col in display.columns:
        display[col] = display[col].map(lambda value: "" if pd.isna(value) else str(value))
    header = "| " + " | ".join(display.columns) + " |"
    sep = "| " + " | ".join(["---"] * len(display.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in display.astype(str).values.tolist()]
    return "\n".join([header, sep, *rows])


def write_report(scores: pd.DataFrame, box_meta: dict[str, Any], out_dir: Path) -> None:
    best = scores.sort_values(["variant_model", "vina_score_kcal_mol"], ascending=[True, True])
    lines = [
        "# Basic EGFR-PACC Docking Report",
        "",
        "## Scope",
        "",
        "This is a coarse AutoDock Vina 1.2.7 docking screen intended as a CSBJ-oriented supplementary structural analysis.",
        "It is not a covalent docking, MD, FEP, or binding free-energy calculation.",
        "",
        "## Inputs",
        "",
        f"- Receptor template: {box_meta['template_pdb']} chain {box_meta['chain']}.",
        f"- Docking box reference ligand: {box_meta['box_reference_ligand']}.",
        f"- Box center: {box_meta['box_center']}.",
        f"- Box size: {box_meta['box_size']} Angstrom.",
        "- Receptor variants: WT, G719S, L861Q, S768I.",
        "- Ligands: afatinib, osimertinib, furmonertinib.",
        "- Receptor point mutants were generated with PDBFixer and were not MD-relaxed.",
        "- Receptor PDBQT files were generated with OpenBabel. This route avoids the RDKit residue-valence parsing failure observed with Meeko receptor preparation on the current Windows/Python runtime.",
        "- Ligands were desalted, embedded in 3D with RDKit, minimized, and converted to PDBQT with Meeko.",
        "",
        "## Score Summary",
        "",
        markdown_table(best),
        "",
        "## Interpretation Boundary",
        "",
        "- EGFR TKIs such as afatinib, osimertinib, and furmonertinib can form covalent interactions; this Vina run models only non-covalent poses.",
        "- PDBFixer-generated point-mutant side chains are not MD-relaxed.",
        "- Receptor PDBQT files were prepared with OpenBabel; future higher-stringency studies should compare against alternative receptor preparation and protonation workflows.",
        "- Scores should be interpreted as coarse pose/compatibility evidence and used only as supplementary support for the structure-informed ranking framework.",
        "- The manuscript should report these values as approximate Vina docking scores, not as experimental affinity or DeltaG_bind.",
        "",
        "## Outputs",
        "",
        "- `docking_scores.csv`",
        "- `basic_pacc_docking_heatmap.png/pdf/svg`",
        "- `basic_pacc_docking_pose_overview.png/pdf/svg`",
        "- `poses/*.pdbqt`",
        "- `receptors/*.pdb` and `receptors/*.pdbqt`",
        "- `ligands/*.sdf` and `ligands/*.pdbqt`",
    ]
    (out_dir / "basic_pacc_docking_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run coarse Vina docking for core EGFR-PACC mutations.")
    parser.add_argument("--out-dir", default=str(ROOT / "reports" / "docking_pacc_core"))
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260620)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    receptors_dir = out_dir / "receptors"
    ligands_dir = out_dir / "ligands"
    poses_dir = out_dir / "poses"
    logs_dir = out_dir / "logs"
    for path in [out_dir, receptors_dir, ligands_dir, poses_dir, logs_dir]:
        path.mkdir(parents=True, exist_ok=True)

    if not VINA_EXE.exists():
        raise FileNotFoundError(f"Vina executable not found: {VINA_EXE}")
    protein_pdb, box, box_meta = extract_template_protein_and_box(out_dir)
    (out_dir / "docking_box_metadata.json").write_text(json.dumps(box_meta, indent=2), encoding="utf-8")

    receptor_pdbqts: dict[str, Path] = {}
    for variant, mutations in MUTATION_SETS.items():
        prepared_pdb = receptors_dir / f"EGFR_4LRM_chainA_{variant}.pdb"
        prepared_pdbqt = receptors_dir / f"EGFR_4LRM_chainA_{variant}.pdbqt"
        pdbfixer_prepare(protein_pdb, prepared_pdb, mutations)
        prepare_receptor_pdbqt(prepared_pdb, prepared_pdbqt, box)
        receptor_pdbqts[variant] = prepared_pdbqt

    ligand_pdbqts: dict[str, Path] = {}
    for drug in DRUGS:
        _, ligand_pdbqt = prepare_ligand_sdf_and_pdbqt(drug, ligands_dir, seed=args.seed)
        ligand_pdbqts[drug] = ligand_pdbqt

    rows: list[dict[str, Any]] = []
    for variant, receptor_pdbqt in receptor_pdbqts.items():
        for drug, ligand_pdbqt in ligand_pdbqts.items():
            pose = poses_dir / f"{variant}__{drug}__vina_pose.pdbqt"
            log = logs_dir / f"{variant}__{drug}__vina.log"
            score = run_vina(
                receptor_pdbqt,
                ligand_pdbqt,
                pose,
                log,
                box,
                exhaustiveness=args.exhaustiveness,
                seed=args.seed,
            )
            rows.append(
                {
                    "variant_model": variant,
                    "drug": drug,
                    "vina_score_kcal_mol": score,
                    "rank_within_variant": math.nan,
                    "receptor_template": "4LRM_chainA",
                    "mutation_model": ";".join(MUTATION_SETS[variant]) if MUTATION_SETS[variant] else "WT",
                    "pose_pdbqt": str(pose),
                    "log_file": str(log),
                    "method_note": "Coarse non-covalent Vina docking; mutant side chains generated with PDBFixer and not MD-relaxed; receptor PDBQT prepared with OpenBabel.",
                }
            )

    scores = pd.DataFrame(rows)
    scores["rank_within_variant"] = scores.groupby("variant_model")["vina_score_kcal_mol"].rank(method="first")
    scores = scores.sort_values(["variant_model", "rank_within_variant"])
    scores.to_csv(out_dir / "docking_scores.csv", index=False)
    make_heatmap(scores, out_dir)
    make_pose_overview(scores, out_dir)
    write_report(scores, box_meta, out_dir)
    print(
        {
            "scores": str((out_dir / "docking_scores.csv").resolve()),
            "report": str((out_dir / "basic_pacc_docking_report.md").resolve()),
            "heatmap": str((out_dir / "basic_pacc_docking_heatmap.png").resolve()),
        }
    )
    print(scores.to_string(index=False))


if __name__ == "__main__":
    main()
