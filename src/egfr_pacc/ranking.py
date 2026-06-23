from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from .bayes_model import class_drug_prior_table
from .classifier import ClassificationResult
from .epdri import EpDriInput, compute_epdri
from .io import config_path, load_yaml


@dataclass(frozen=True)
class CaseFeatures:
    patient_id: str
    egfr_variants: list[str]
    brain_metastasis: bool = False
    leptomeningeal_metastasis: bool = False
    treatment_line: int = 1
    tp53_mut: bool = False
    rb1_loss: bool = False
    met_amplification: bool = False
    pik3ca_mut: bool = False


def has_bypass_risk(case: CaseFeatures) -> bool:
    return any([case.tp53_mut, case.rb1_loss, case.met_amplification, case.pik3ca_mut])


def rank_drugs(
    case: CaseFeatures,
    classification: ClassificationResult,
    use_pymc: bool = False,
) -> pd.DataFrame:
    priors = load_yaml(config_path("priors.yaml"))
    drug_props = load_yaml(config_path("drug_properties.yaml"))["drugs"]
    baselines = load_yaml(config_path("baselines.yaml"))
    modifiers = priors.get("modifiers", {})
    pooled_priors = class_drug_prior_table(config_path("priors.yaml"), use_pymc=use_pymc)
    class_priors = pooled_priors[
        pooled_priors["primary_class"].eq(classification.primary_class)
    ].copy()

    baseline0 = baselines["baseline_0"]["recommendations"].get(classification.primary_class, ["specialist_review"])
    cns_present = case.brain_metastasis or case.leptomeningeal_metastasis
    bypass = has_bypass_risk(case)

    rows = []
    for _, prior in class_priors.iterrows():
        drug = prior["drug"]
        props = drug_props.get(drug, {})
        epdri = compute_epdri(
            EpDriInput(
                base_prior_score=float(prior.get("pooled_prior_score", prior.get("locked_prior_score", 0.5))),
                evidence_level=str(prior.get("evidence_level", "D")),
                cns_activity=str(props.get("cns_activity", "unknown")),
                cns_present=cns_present,
                bypass_risk=bypass,
                compound_conflict=classification.compound_conflict_flag,
                treatment_line=case.treatment_line,
                classification_uncertainty_delta=classification.uncertainty_delta,
                modifiers=modifiers,
            )
        )
        rows.append(
            {
                "patient_id": case.patient_id,
                "drug": drug,
                "epdri_score": epdri.epdri_score,
                "recommendation_grade": epdri.recommendation_grade,
                "uncertainty_tier": epdri.uncertainty_tier,
                "primary_class": classification.primary_class,
                "baseline_0_recommended": drug in baseline0,
                "evidence_level": epdri.evidence_level,
                "evidence_cap_grade": epdri.evidence_cap_grade,
                "locked_prior_score": round(float(prior.get("locked_prior_score", 0.5)), 3),
                "pooled_prior_score": round(float(prior.get("pooled_prior_score", 0.5)), 3),
                "layer1_hierarchical_prior": epdri.layer_scores["layer1_hierarchical_prior"],
                "layer2_clinical_modified": epdri.layer_scores["layer2_clinical_modified"],
                "layer3_uncertainty_discounted": epdri.layer_scores["layer3_uncertainty_discounted"],
                "uncertainty_discount": epdri.uncertainty_discount,
                "pooling_method": prior.get("pooling_method", ""),
                "prior_source": prior.get("prior_source", ""),
                "modifier_notes": "; ".join(epdri.modifier_notes) if epdri.modifier_notes else "none",
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["epdri_score", "baseline_0_recommended"], ascending=[False, False]
    )
