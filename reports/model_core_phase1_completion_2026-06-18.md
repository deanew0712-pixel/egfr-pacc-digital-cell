# 第一阶段模型核心完成记录

生成时间：2026-06-18 Asia/Singapore

## 完成范围

已完成 Layer 4 三个核心子模块的可运行代码：

1. `src/egfr_pacc/classifier.py`
   - 新增 `classify_robichaux()`。
   - 新增 `resolve_compound()`。
   - 保留 `classify_variants()` 作为旧脚本兼容入口。
   - Robichaux Fig.5 四分类工作规则写入 `configs/mutation_class_rules.yaml`：
     `pacc`、`classical_like`、`t790m_like`、`ex20ins`。
   - 复合突变消解逻辑写入 YAML：
     `pacc_plus_classical_like`、`pacc_plus_t790m_like`、`any_plus_ex20ins`。

2. `src/egfr_pacc/bayes_model.py`
   - 新增 `load_locked_prior_matrix()`。
   - 新增 `deterministic_partial_pooling()`。
   - 新增 `build_pymc_model()`。
   - 新增 `class_drug_prior_table()`。
   - PyMC 模型只读取 `configs/priors.yaml` 的锁定先验，不读取本地临床结局。
   - 当前默认运行使用 deterministic locked-prior partial pooling；项目 `.venv`
     已安装并验证 PyMC 采样路径，可通过 `rank_drugs(..., use_pymc=True)` 显式启用。
   - PyMC 结构为 class-level partial pooling：`global_logit`、`class_offset`、
     `sigma_pair` 与非中心化 pair-level `theta_raw/theta_logit`。

3. `src/egfr_pacc/epdri.py`
   - 新增三层 collapse 版 EPDRI：
     - Layer 1：class-level hierarchical prior。
     - Layer 2：CNS、bypass/co-mutation、compound conflict、treatment line modifiers。
     - Layer 3：uncertainty discount + evidence-level grade cap。
   - 输出 `epdri_score`、`recommendation_grade`、`uncertainty_tier`、
     `evidence_cap_grade`、三层 layer score 和 modifier notes。

## 已接入 demo pipeline

`src/egfr_pacc/ranking.py` 已改为：

- 从 `bayes_model.class_drug_prior_table()` 读取 pooled class-drug prior。
- 从 `epdri.compute_epdri()` 计算最终 EPDRI。
- 输出 layer score 与 evidence cap，便于后续论文图表和误差分析。

已验证命令：

```bash
.venv/bin/python -m py_compile src/egfr_pacc/*.py scripts/*.py
.venv/bin/python scripts/build_preliminary_digital_cell.py
.venv/bin/python scripts/run_demo_case.py
.venv/bin/python scripts/run_batch_cases.py
```

PyMC 采样路径已验证：

```bash
PYTHONPATH=src .venv/bin/python - <<'PY'
from egfr_pacc.bayes_model import class_drug_prior_table
class_drug_prior_table(use_pymc=True, draws=500, tune=500)
PY
```

正式 posterior 已重跑：

- chains：4
- tune：2000
- draws：2000
- target_accept：0.995
- divergences：0
- max R-hat：1.0017
- min bulk ESS：1567.4825
- min tail ESS：1665.7512

生成文件：

- `reports/pymc_phase1_prior_posterior.csv`
- `reports/pymc_phase1_prior_posterior_4chains_2000draws.csv`
- `reports/pymc_phase1_prior_diagnostics_4chains_2000draws.csv`
- `reports/demo_case_ranking_pymc.csv`
- `reports/demo_case_ranking_pymc_formal_4chains_2000draws.csv`
- `reports/pymc_environment_2026-06-18.txt`

## Demo case 输出

Demo case：

- EGFR variants：G719S + L861Q。
- CNS disease：yes。
- TP53 mutation：yes。

分类结果：

- classes：`pacc` + `classical_like`。
- primary class：`pacc`。
- compound rule：`preserve_both_classes_with_pacc_bias`。
- uncertainty delta：1。

Top ranking snapshot：

| Drug | EPDRI | Grade | Uncertainty |
| --- | ---: | --- | --- |
| afatinib | 0.552 | C | high |
| osimertinib | 0.487 | C | high |
| furmonertinib | 0.473 | D | high |
| lazertinib | 0.466 | D | high |
| aumolertinib | 0.459 | D | high |

解释：

- Demo case 带有 compound conflict、TP53 bypass risk 和 CNS disease。
- 因 uncertainty high，最终推荐等级受到 evidence/uncertainty cap 约束。
- 伏美替尼目前仍是 locked prior D + emerging evidence，不直接覆盖 locked prior，
  因此分数接近 osimertinib，但 grade 仍被 evidence cap 限制为 D。

## Figure 草图

已生成论文草图基础图：

- `reports/figures/figure1_framework_flow.png`
- `reports/figures/figure2_demo_case_ranking.png`

Figure 1：EGFR-PACC digital cell MVP framework。

Figure 2：demo case collapsed EPDRI ranking output。

## 当前注意事项

- PyMC 已安装在项目本地虚拟环境 `.venv` 中。新电脑复现时运行：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip setuptools wheel
.venv/bin/python -m pip install -r requirements.txt -r requirements-pymc.txt
```

- 当前已验证版本见 `reports/pymc_environment_2026-06-18.txt`。
- 默认 demo/batch pipeline 仍不依赖 PyMC，避免每个病例重复采样；需要 PyMC
  posterior prior 时显式调用 `use_pymc=True`。
- PyMC 模型不应导入 30 例本地 outcome 进行估计；本地 outcome 仅用于后续锁表后的校准与验证。
- 正式 posterior 表已使用 4 chains、2000 tune、2000 draws 重跑，并达到 0 divergence。
