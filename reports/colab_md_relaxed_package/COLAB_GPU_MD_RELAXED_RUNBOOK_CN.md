# Google Colab 免费 GPU 跑 EGFR-PACC MD-Relaxed Docking

日期：2026-06-21

## 目的

在 Google Colab T4 GPU 或小型云 GPU 实例上完成：

1. WT/G719S/L861Q/S768I EGFR kinase-domain receptor restrained minimization；
2. relaxed receptor PDBQT 生成；
3. relaxed receptor × afatinib/osimertinib/furmonertinib AutoDock Vina 重跑；
4. relaxed vs unrelaxed docking rank/score sensitivity 对比；
5. 可作为 Major Revision 阶段补充材料的表格和图。

该分析应写作：

> OpenMM restrained receptor-relaxation sensitivity analysis.

不应写作：

> production MD, FEP, or definitive covalent binding free-energy validation.

## 本地已准备文件

运行脚本：

- `scripts/colab_md_relaxed_egfr_pacc.py`

需要上传到 Colab 的输入文件：

- `reports/docking_pacc_core/docking_box_metadata.json`
- `reports/docking_pacc_core/docking_scores.csv`
- `reports/docking_pacc_core/receptors/EGFR_4LRM_chainA_WT.pdb`
- `reports/docking_pacc_core/receptors/EGFR_4LRM_chainA_G719S.pdb`
- `reports/docking_pacc_core/receptors/EGFR_4LRM_chainA_L861Q.pdb`
- `reports/docking_pacc_core/receptors/EGFR_4LRM_chainA_S768I.pdb`
- `reports/docking_pacc_core/ligands/afatinib.pdbqt`
- `reports/docking_pacc_core/ligands/osimertinib.pdbqt`
- `reports/docking_pacc_core/ligands/furmonertinib.pdbqt`

建议在 Colab 中整理为：

```text
colab_input/
  docking_box_metadata.json
  docking_scores.csv
  receptors/
    EGFR_4LRM_chainA_WT.pdb
    EGFR_4LRM_chainA_G719S.pdb
    EGFR_4LRM_chainA_L861Q.pdb
    EGFR_4LRM_chainA_S768I.pdb
  ligands/
    afatinib.pdbqt
    osimertinib.pdbqt
    furmonertinib.pdbqt
```

## Colab 环境设置

在 Colab 菜单选择：

```text
Runtime -> Change runtime type -> T4 GPU
```

然后执行：

```bash
!nvidia-smi
```

安装依赖：

```bash
!pip -q install openmm pandas matplotlib tabulate
!apt-get -qq update
!apt-get -qq install -y openbabel
```

安装 AutoDock Vina。优先用 conda/micromamba 或下载官方 Linux binary；如果已有 `vina` 可执行文件，确保：

```bash
!vina --help
```

能正常输出帮助信息。

## 上传输入

方式一：手动上传 zip 后解压。

```python
from google.colab import files
uploaded = files.upload()
```

然后：

```bash
!unzip -q colab_egfr_pacc_md_input.zip -d .
```

方式二：挂载 Google Drive。

```python
from google.colab import drive
drive.mount('/content/drive')
```

再复制输入目录到 `/content/colab_input`。

## 运行

上传或复制脚本到 `/content/colab_md_relaxed_egfr_pacc.py` 后运行：

```bash
!python colab_md_relaxed_egfr_pacc.py \
  --input-dir colab_input \
  --output-dir colab_output/openmm_gpu_relaxed \
  --platform CUDA \
  --max-iterations 500 \
  --md-steps 1000 \
  --exhaustiveness 8
```

如果 T4 时间紧张，可以先跑快速版：

```bash
!python colab_md_relaxed_egfr_pacc.py \
  --input-dir colab_input \
  --output-dir colab_output/openmm_gpu_relaxed_quick \
  --platform CUDA \
  --max-iterations 200 \
  --md-steps 0 \
  --exhaustiveness 4
```

正式补充材料建议使用：

- `--max-iterations 500`
- `--md-steps 1000`
- `--exhaustiveness 8`

## 输出文件

主要输出：

- `openmm_gpu_relaxed_receptor_metrics.csv`
- `openmm_gpu_relaxed_docking_scores.csv`
- `openmm_gpu_relaxed_vs_unrelaxed_comparison.csv`
- `openmm_gpu_relaxed_docking_report.md`
- `openmm_gpu_relaxed_docking_sensitivity.png/pdf/svg`
- `receptors/*_openmm_gpu_relaxed.pdb`
- `receptors/*_openmm_gpu_relaxed.pdbqt`
- `poses/*_openmm_gpu_relaxed_vina_pose.pdbqt`

下载输出：

```bash
!zip -qr openmm_gpu_relaxed_outputs.zip colab_output/openmm_gpu_relaxed
```

```python
from google.colab import files
files.download("openmm_gpu_relaxed_outputs.zip")
```

## 结果回填到本项目

下载后解压到：

```text
reports/docking_md_relaxed_gpu/
```

建议补充到投稿图目录：

```text
reports/figures_publication/figureS_openmm_gpu_relaxed_docking_sensitivity.png
reports/figures_publication/figureS_openmm_gpu_relaxed_docking_sensitivity.pdf
reports/figures_publication/figureS_openmm_gpu_relaxed_docking_sensitivity.svg
```

## 稿件写法

推荐句子：

> As a receptor-conformation sensitivity analysis, WT and PACC-mutant kinase-domain structures were subjected to OpenMM restrained minimization followed by repeated Vina docking. Rank stability after restrained relaxation was assessed against the unrelaxed receptor docking results.

限制句子：

> This analysis was not intended to estimate covalent reaction energetics or absolute binding free energies, and should be interpreted as a receptor-conformation sensitivity check.

## 如果 Colab 报错

1. `No platform named CUDA`

确认 Colab 已选择 GPU：

```bash
!nvidia-smi
```

2. `vina not found`

安装或上传 Vina binary，然后运行时加入：

```bash
--vina /content/vina
```

3. OpenBabel PDBQT 转换失败

确认：

```bash
!obabel -V
```

4. 时间不够

先用 quick 参数：

```text
--max-iterations 200 --md-steps 0 --exhaustiveness 4
```

拿到结果后再决定是否跑正式版。
