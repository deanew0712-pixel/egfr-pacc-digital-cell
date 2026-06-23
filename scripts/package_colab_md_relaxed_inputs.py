from __future__ import annotations

import shutil
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "reports" / "docking_pacc_core"
PACKAGE_DIR = ROOT / "reports" / "colab_md_relaxed_package"
INPUT_DIR = PACKAGE_DIR / "colab_input"
ZIP_PATH = PACKAGE_DIR / "colab_egfr_pacc_md_input.zip"


def copy_file(src: Path, dst: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(src)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def main() -> None:
    if PACKAGE_DIR.exists():
        shutil.rmtree(PACKAGE_DIR)
    INPUT_DIR.mkdir(parents=True, exist_ok=True)

    copy_file(SOURCE / "docking_box_metadata.json", INPUT_DIR / "docking_box_metadata.json")
    copy_file(SOURCE / "docking_scores.csv", INPUT_DIR / "docking_scores.csv")
    for variant in ["WT", "G719S", "L861Q", "S768I"]:
        copy_file(
            SOURCE / "receptors" / f"EGFR_4LRM_chainA_{variant}.pdb",
            INPUT_DIR / "receptors" / f"EGFR_4LRM_chainA_{variant}.pdb",
        )
    for drug in ["afatinib", "osimertinib", "furmonertinib"]:
        copy_file(SOURCE / "ligands" / f"{drug}.pdbqt", INPUT_DIR / "ligands" / f"{drug}.pdbqt")

    copy_file(ROOT / "scripts" / "colab_md_relaxed_egfr_pacc.py", PACKAGE_DIR / "colab_md_relaxed_egfr_pacc.py")
    copy_file(ROOT / "docs" / "COLAB_GPU_MD_RELAXED_RUNBOOK_CN.md", PACKAGE_DIR / "COLAB_GPU_MD_RELAXED_RUNBOOK_CN.md")
    copy_file(ROOT / "colab_egfr_pacc_md_relaxed_runner.ipynb", PACKAGE_DIR / "colab_egfr_pacc_md_relaxed_runner.ipynb")

    with ZipFile(ZIP_PATH, "w", compression=ZIP_DEFLATED) as zf:
        for path in PACKAGE_DIR.rglob("*"):
            if path == ZIP_PATH or path.is_dir():
                continue
            zf.write(path, path.relative_to(PACKAGE_DIR))

    manifest = sorted(str(path.relative_to(PACKAGE_DIR)) for path in PACKAGE_DIR.rglob("*") if path.is_file())
    (PACKAGE_DIR / "MANIFEST.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(f"Wrote {ZIP_PATH}")


if __name__ == "__main__":
    main()
