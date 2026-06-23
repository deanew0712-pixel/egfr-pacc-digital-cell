from __future__ import annotations

import csv
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data" / "processed" / "reference"
CONFIG_OUT = ROOT / "configs" / "open_fallback_digital_cell.yaml"
REPORT_OUT = ROOT / "reports" / "preliminary_digital_cell_build_report.md"

GENERATED_AT = datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def grouped_counts(rows: list[dict]) -> dict:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[row["mutation"]].append(row)
    summary = {}
    for mutation, items in sorted(grouped.items()):
        source_counts = {}
        for item in items:
            source = item["source"]
            try:
                count = float(item.get("record_count") or 0)
            except ValueError:
                count = 0
            source_counts[source] = source_counts.get(source, 0) + count
        summary[mutation] = {
            "sources": sorted({item["source"] for item in items}),
            "source_record_counts": source_counts,
            "evidence_rows": len(items),
        }
    return summary


def frequency_summary(rows: list[dict]) -> dict:
    summary: dict[str, dict] = defaultdict(dict)
    for row in rows:
        cohort = row["cohort"]
        metric = row["metric"]
        value = row.get("frequency", "")
        summary[cohort][metric] = {
            "frequency": float(value) if value not in {"", None} else None,
            "sample_count": int(float(row["sample_count"])) if row.get("sample_count") else 0,
            "mutation_count": int(float(row["mutation_count"])) if row.get("mutation_count") else 0,
            "source": row["source"],
            "ancestry_race_available": str(row.get("ancestry_race_available")).lower() == "true",
        }
    return dict(sorted(summary.items()))


def cns_summary(rows: list[dict]) -> dict:
    out = {}
    for row in rows:
        out[row["drug"]] = {
            "cns_activity": row["cns_activity"],
            "cns_modifier_score": float(row["cns_modifier_score"]),
            "pubchem_cid": row.get("pubchem_cid", ""),
            "chembl_id": row.get("chembl_id", ""),
            "openfda_label_found": str(row.get("openfda_label_found")).lower() == "true",
            "fda_label_evidence_strength": row.get("fda_label_evidence_strength", ""),
            "fda_label_effective_time": row.get("fda_label_effective_time", ""),
            "source": row["source"],
            "source_url_or_id": row["source_url_or_id"],
            "license_note": row["license_note"],
        }
    return out


def source_priority_summary(rows: list[dict]) -> dict:
    out: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        out[row["output_table"]].append(
            {
                "field": row["field"],
                "priority_rank": int(row["priority_rank"]),
                "source": row["source"],
                "dependency_status": row["dependency_status"],
                "source_url_or_id": row["source_url_or_id"],
            }
        )
    return {key: sorted(value, key=lambda item: item["priority_rank"]) for key, value in sorted(out.items())}


def completion_checks(paths: dict[str, Path], config: dict) -> list[dict]:
    checks = []
    for name, path in paths.items():
        rows = read_csv(path)
        checks.append({"item": name, "status": "complete" if rows else "empty", "rows": len(rows), "path": str(path)})
    for optional_name in ["cosmic", "genie", "drugbank"]:
        path = ROOT / "data" / "raw" / optional_name
        has_files = path.exists() and any(p for p in path.iterdir() if not p.name.startswith("._"))
        checks.append(
            {
                "item": f"optional_enhanced_{optional_name}",
                "status": "available" if has_files else "warning_missing_optional",
                "rows": "",
                "path": str(path),
            }
        )
    checks.append(
        {
            "item": "model_config",
            "status": "complete" if config.get("mutation_catalog") and config.get("cns_modifier_table") else "incomplete",
            "rows": "",
            "path": str(CONFIG_OUT),
        }
    )
    return checks


def write_report(config: dict, checks: list[dict]) -> None:
    lines = [
        "# EGFR-PACC 初步数字细胞构建报告",
        "",
        f"生成时间：{GENERATED_AT}",
        "",
        "## 构建策略",
        "",
        "本轮已将 COSMIC、GENIE、DrugBank 从必需依赖改为 optional enhanced dependency。初步数字细胞使用开放替代源构建：cBioPortal、TCGA/GDC、CIViC、ClinVar、Cancer Hotspots、ChEMBL、PubChem、DailyMed、openFDA 和 manual literature seed。",
        "",
        "## 输出文件",
        "",
        f"- 机器可读模型快照：`{CONFIG_OUT.relative_to(ROOT)}`",
        "- mutation catalog：`data/processed/reference/egfr_open_mutation_catalog.csv`",
        "- enhanced mutation catalog：`data/processed/reference/egfr_enhanced_mutation_catalog.csv`（如 COSMIC 已下载并解析）",
        "- frequency prior：`data/processed/reference/open_egfr_pacc_frequency.csv`",
        "- CNS modifier：`data/processed/reference/egfr_tki_cns_modifier_table.csv`",
        "- source priority：`data/processed/reference/source_priority_map.csv`",
        "",
        "## 完成度检查",
        "",
        "| 项目 | 状态 | 行数 | 路径 |",
        "| --- | --- | ---: | --- |",
    ]
    for check in checks:
        lines.append(f"| {check['item']} | {check['status']} | {check['rows']} | `{check['path']}` |")

    lines.extend(
        [
            "",
            "## 当前可用模型组件",
            "",
            f"- EGFR 重点突变目录：{len(config['mutation_catalog'])} 个 mutation bucket。",
            f"- EGFR 增强突变目录：{len(config.get('enhanced_mutation_catalog', {}))} 个 mutation bucket。",
            f"- 开放队列频率：{len(config['frequency_priors'])} 个 cohort。",
            f"- EGFR 药物 CNS modifier：{len(config['cns_modifier_table'])} 个药物。",
            "",
            "## 关键限制",
            "",
            "- CNS modifier 目前包含 manual literature seed，尚未完成正式文献锁表。",
            "- OncoKB、COSMIC、GENIE、DrugBank 均可增强模型，但缺失不会阻断 MVP。",
            "- 当前模型可用于候选排序与流程验证，不能视为临床验证模型。",
            "",
            "## 建议下一步",
            "",
            "1. 锁定 Robichaux、LUX-Lung、ACHILLES/TORG1834、UNICORN 的 mutation x drug ORR/PFS 表。",
            "2. 用锁表更新 `configs/priors.yaml`，保留 `source_priority_map.csv` 的字段级溯源。",
            "3. 将 37 例临床数据先导入 feature-only 表，确认不含 PFS/ORR/OS。",
            "4. 冻结 priors/baselines 后再打开 outcome，做 EPDRI vs Baseline-1 的探索性验证。",
            "5. 后续有 OncoKB/COSMIC/GENIE/DrugBank 时作为 enhanced layer 增量比较，不回改已冻结 open baseline。",
        ]
    )
    REPORT_OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    paths = {
        "egfr_open_mutation_catalog": REFERENCE / "egfr_open_mutation_catalog.csv",
        "egfr_enhanced_mutation_catalog": REFERENCE / "egfr_enhanced_mutation_catalog.csv",
        "open_egfr_pacc_frequency": REFERENCE / "open_egfr_pacc_frequency.csv",
        "egfr_tki_cns_modifier_table": REFERENCE / "egfr_tki_cns_modifier_table.csv",
        "source_priority_map": REFERENCE / "source_priority_map.csv",
    }
    mutation_rows = read_csv(paths["egfr_open_mutation_catalog"])
    enhanced_mutation_rows = read_csv(paths["egfr_enhanced_mutation_catalog"]) if paths["egfr_enhanced_mutation_catalog"].exists() else []
    frequency_rows = read_csv(paths["open_egfr_pacc_frequency"])
    cns_rows = read_csv(paths["egfr_tki_cns_modifier_table"])
    priority_rows = read_csv(paths["source_priority_map"])

    config = {
        "metadata": {
            "generated_at": GENERATED_AT,
            "status": "preliminary_open_fallback_digital_cell",
            "clinical_use": "not_for_clinical_decision",
            "optional_enhanced_dependencies": ["COSMIC", "GENIE", "DrugBank", "OncoKB"],
            "required_before_outcome_review": ["manual_literature_orr_pfs_lock_table"],
        },
        "mutation_catalog": grouped_counts(mutation_rows),
        "enhanced_mutation_catalog": grouped_counts(enhanced_mutation_rows) if enhanced_mutation_rows else {},
        "frequency_priors": frequency_summary(frequency_rows),
        "cns_modifier_table": cns_summary(cns_rows),
        "source_priority": source_priority_summary(priority_rows),
    }
    CONFIG_OUT.write_text(yaml.safe_dump(config, sort_keys=False, allow_unicode=True), encoding="utf-8")
    checks = completion_checks(paths, config)
    write_report(config, checks)
    print({"config": str(CONFIG_OUT), "report": str(REPORT_OUT), "checks": len(checks)})


if __name__ == "__main__":
    main()
