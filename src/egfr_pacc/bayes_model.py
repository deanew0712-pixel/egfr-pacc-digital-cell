from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .io import config_path, load_yaml


def _clip_probability(value: float) -> float:
    return max(0.01, min(0.99, float(value)))


def logit(value: float) -> float:
    p = _clip_probability(value)
    return math.log(p / (1.0 - p))


def inv_logit(value: float) -> float:
    return 1.0 / (1.0 + math.exp(-float(value)))


@dataclass(frozen=True)
class PriorMatrix:
    rows: pd.DataFrame
    classes: list[str]
    drugs: list[str]


def load_locked_prior_matrix(priors_path: str | Path | None = None) -> PriorMatrix:
    priors = load_yaml(priors_path or config_path("priors.yaml"))
    rows: list[dict[str, Any]] = []
    for class_name, drug_map in priors.get("class_drug_priors", {}).items():
        for drug, entry in drug_map.items():
            score = _clip_probability(float(entry.get("score", 0.5)))
            rows.append(
                {
                    "primary_class": class_name,
                    "drug": drug,
                    "locked_prior_score": score,
                    "locked_prior_logit": logit(score),
                    "evidence_level": entry.get("evidence_level", "D"),
                    "prior_source": entry.get("source", ""),
                }
            )

    frame = pd.DataFrame(rows)
    classes = sorted(frame["primary_class"].unique().tolist()) if not frame.empty else []
    drugs = sorted(frame["drug"].unique().tolist()) if not frame.empty else []
    return PriorMatrix(rows=frame, classes=classes, drugs=drugs)


def deterministic_partial_pooling(
    priors_path: str | Path | None = None,
    class_weight: float = 0.2,
    global_weight: float = 0.1,
) -> pd.DataFrame:
    """Class-level hierarchical shrinkage without using local outcomes.

    This is the deterministic fallback used when PyMC is unavailable. It reads
    only locked priors, shrinks each class-drug logit toward the class mean and
    global mean, and returns pooled probabilities for EPDRI.
    """
    matrix = load_locked_prior_matrix(priors_path)
    frame = matrix.rows.copy()
    if frame.empty:
        return frame

    global_mean = frame["locked_prior_logit"].mean()
    class_means = frame.groupby("primary_class")["locked_prior_logit"].transform("mean")
    direct_weight = max(0.0, 1.0 - class_weight - global_weight)
    pooled_logit = (
        direct_weight * frame["locked_prior_logit"]
        + class_weight * class_means
        + global_weight * global_mean
    )
    frame["pooled_logit_mean"] = pooled_logit
    frame["pooled_prior_score"] = pooled_logit.map(inv_logit).clip(0.01, 0.99)
    frame["pooling_method"] = "deterministic_locked_prior_partial_pooling"
    frame["pooling_note"] = "No local outcomes used; PyMC not required for this fallback."
    return frame


def build_pymc_model(priors_path: str | Path | None = None, prior_sigma: float = 0.35):
    """Build a compact PyMC hierarchical prior model.

    Parameters are intentionally minimal:

    - global_logit
    - class_offset
    - sigma_pair
    - pooled class-drug theta_logit

    The locked prior logits from `priors.yaml` enter as pseudo-observed prior
    anchors. No local clinical outcome data are read or estimated here.
    """
    try:
        import pymc as pm
    except ImportError as exc:
        raise RuntimeError(
            "PyMC is not installed. Install `pymc` to sample this model, or use "
            "`deterministic_partial_pooling()` for the runnable fallback."
        ) from exc

    matrix = load_locked_prior_matrix(priors_path)
    frame = matrix.rows.reset_index(drop=True)
    if frame.empty:
        raise ValueError("No class-drug priors found in priors.yaml.")

    class_index = {name: idx for idx, name in enumerate(matrix.classes)}
    class_codes = frame["primary_class"].map(class_index).to_numpy()
    prior_logits = frame["locked_prior_logit"].to_numpy()
    global_center = float(frame["locked_prior_logit"].mean())

    coords = {
        "class": matrix.classes,
        "pair": list(range(len(frame))),
    }
    with pm.Model(coords=coords) as model:
        global_logit = pm.Normal("global_logit", mu=global_center, sigma=1.0)
        sigma_class = pm.HalfNormal("sigma_class", sigma=0.75)
        sigma_pair = pm.HalfNormal("sigma_pair", sigma=0.75)
        class_offset_raw = pm.Normal("class_offset_raw", mu=0.0, sigma=1.0, dims="class")
        class_offset = pm.Deterministic(
            "class_offset",
            class_offset_raw * sigma_class,
            dims="class",
        )
        class_mean_logit = pm.Deterministic(
            "class_mean_logit",
            global_logit + class_offset,
            dims="class",
        )
        theta_raw = pm.Normal(
            "theta_raw",
            mu=0.0,
            sigma=1.0,
            dims="pair",
        )
        theta_logit = pm.Deterministic(
            "theta_logit",
            class_mean_logit[class_codes] + theta_raw * sigma_pair,
            dims="pair",
        )
        pm.Normal(
            "locked_prior_logit",
            mu=theta_logit,
            sigma=prior_sigma,
            observed=prior_logits,
            dims="pair",
        )
        pm.Deterministic("theta_probability", pm.math.sigmoid(theta_logit), dims="pair")
    return model, frame


def class_drug_prior_table(
    priors_path: str | Path | None = None,
    use_pymc: bool = False,
    draws: int = 500,
    tune: int = 500,
    chains: int = 2,
    random_seed: int = 20260618,
    target_accept: float = 0.95,
) -> pd.DataFrame:
    """Return class-drug prior probabilities for ranking.

    `use_pymc=False` is the default so the project remains runnable on a clean
    open-data environment. Set `use_pymc=True` after installing PyMC to sample
    the compact hierarchical model.
    """
    if not use_pymc:
        return deterministic_partial_pooling(priors_path)

    try:
        import pymc as pm
    except ImportError:
        return deterministic_partial_pooling(priors_path)

    model, frame = build_pymc_model(priors_path)
    with model:
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            cores=1,
            random_seed=random_seed,
            target_accept=target_accept,
            progressbar=False,
        )
    theta = idata.posterior["theta_probability"].mean(dim=("chain", "draw")).to_numpy()
    out = frame.copy()
    out["pooled_prior_score"] = theta
    out["pooled_logit_mean"] = [logit(x) for x in theta]
    out["pooling_method"] = "pymc_hierarchical_locked_prior_partial_pooling"
    out["pooling_note"] = "Sampled from locked priors only; no local outcomes used."
    return out
