"""Cross-cohort calibration arm.

Inputs:
    - Genomics TIME profiles (Arm 2 output, ~50 in `time_profiles/tcga_hnsc/`)
    - IHC TIME profiles (Arm 1 output, ~5 in `time_profiles/deepliif_sample/`)

Method (Approach B from `~/Downloads/AI/P4-IHC-Genomics-TIME-Plan.md` §4):
    1. For each IHC profile, find the K-nearest TCGA-HNSC neighbors on
       subsite metadata (or by overall profile similarity when subsite
       is not informative).
    2. Per cell-type, fit a linear calibration mapping
       `ihc_density = a * genomics_density + b`
       across the (IHC, mean-of-NN-genomics) paired set.
    3. Hold-out validation: leave-one-IHC-out, refit, report MAE +
       immune-phenotype agreement.
    4. Surface `predict_time_from_genomics(genomics_profile) -> TIMEProfile`
       that applies the calibration and writes the calibrated prediction
       to `out_dir/calibration/calibrated_<id>.json`.

Caveats and scope (v0.3):
    With n=5 IHC and n=~50 genomics, the calibration is a *demonstration*
    of the integration pattern, not a statistical claim. The README climax
    table reports the held-out validation numbers transparently — the
    intercept-only baseline is included so a reader can see how much of
    the signal the calibration is actually adding on top of "predict the
    cohort mean."
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from hnscc_time.time_schema import (
    CELL_TYPES,
    Provenance,
    RegionProfile,
    TIMEProfile,
    derive_immune_phenotype,
    derive_til_score,
)

CALIBRATION_METHOD = "nearest_neighbor_subsite + linear_per_cell_type"
CALIBRATION_VERSION = "v0.3.0"


def load_profiles(profiles_dir: Path) -> dict[str, list[TIMEProfile]]:
    """Load TIME profiles from disk grouped by cohort.

    Returns {"tcga_hnsc": [...], "deepliif_sample": [...]}
    """
    out: dict[str, list[TIMEProfile]] = {}
    for cohort_subdir in profiles_dir.iterdir():
        if not cohort_subdir.is_dir():
            continue
        cohort_name = cohort_subdir.name
        profiles: list[TIMEProfile] = []
        for json_path in sorted(cohort_subdir.glob("*.json")):
            data = json.loads(json_path.read_text())
            # The TIMEProfile pydantic model handles nested RegionProfile / Provenance.
            profiles.append(TIMEProfile.model_validate(data))
        out[cohort_name] = profiles
    return out


def _profile_to_vector(profile: TIMEProfile) -> np.ndarray:
    """Per-cell-type density vector from tumor_core region."""
    core = profile.regions.get("tumor_core")
    if core is None:
        return np.zeros(len(CELL_TYPES))
    return np.array([
        core.CD3_density,
        core.CD8_density,
        core.FoxP3_density,
        core.PanCK_density,
    ])


def _nearest_neighbors(
    target_vector: np.ndarray,
    reference_vectors: np.ndarray,
    k: int = 5,
) -> np.ndarray:
    """Indices of k-NN in reference set (Euclidean over profile vector)."""
    dists = np.linalg.norm(reference_vectors - target_vector, axis=1)
    return np.argsort(dists)[:k]


def fit_calibration(
    ihc_profiles: list[TIMEProfile],
    genomics_profiles: list[TIMEProfile],
    k: int = 5,
) -> dict[str, dict[str, float]]:
    """Fit per-cell-type linear calibration mapping.

    Returns:
        {cell_type: {"a": slope, "b": intercept, "n_pairs": int}}

    Per cell type, builds (x, y) pairs where:
        x = mean k-NN genomics density for the IHC sample's profile shape
        y = IHC density
    Then fits `y = a*x + b` via least squares.
    """
    if not ihc_profiles or not genomics_profiles:
        raise ValueError(
            f"calibration needs >=1 IHC and >=1 genomics profile; "
            f"got {len(ihc_profiles)} IHC, {len(genomics_profiles)} genomics"
        )

    ihc_mat = np.array([_profile_to_vector(p) for p in ihc_profiles])
    gen_mat = np.array([_profile_to_vector(p) for p in genomics_profiles])

    coefs: dict[str, dict[str, float]] = {}
    for ci, cell in enumerate(CELL_TYPES):
        xs: list[float] = []
        ys: list[float] = []
        for ihc_i, ihc_vec in enumerate(ihc_mat):
            nn_idx = _nearest_neighbors(ihc_vec, gen_mat, k=min(k, len(gen_mat)))
            xs.append(float(gen_mat[nn_idx, ci].mean()))
            ys.append(float(ihc_mat[ihc_i, ci]))
        x = np.array(xs)
        y = np.array(ys)
        # Least-squares fit y = a*x + b (closed form).
        if len(x) < 2 or x.std() == 0:
            a, b = 0.0, float(y.mean())
        else:
            a = float(np.cov(x, y, bias=True)[0, 1] / (x.var() + 1e-12))
            b = float(y.mean() - a * x.mean())
        coefs[cell] = {"a": a, "b": b, "n_pairs": int(len(x))}
    return coefs


def loo_validate(
    ihc_profiles: list[TIMEProfile],
    genomics_profiles: list[TIMEProfile],
    k: int = 5,
) -> dict[str, dict[str, float]]:
    """Leave-one-IHC-out validation. Per cell type returns MAE."""
    if len(ihc_profiles) < 2:
        return {cell: {"loo_mae": float("nan"), "intercept_mae": float("nan")} for cell in CELL_TYPES}

    ihc_mat = np.array([_profile_to_vector(p) for p in ihc_profiles])
    gen_mat = np.array([_profile_to_vector(p) for p in genomics_profiles])

    results: dict[str, dict[str, float]] = {}
    for ci, cell in enumerate(CELL_TYPES):
        cal_residuals: list[float] = []
        intercept_residuals: list[float] = []
        for held_out_i in range(len(ihc_profiles)):
            kept = [i for i in range(len(ihc_profiles)) if i != held_out_i]
            kept_ihc = [ihc_profiles[i] for i in kept]
            coefs = fit_calibration(kept_ihc, genomics_profiles, k=k)
            cal = coefs[cell]
            # Predict held-out by mean-NN genomics * a + b.
            held_vec = ihc_mat[held_out_i]
            nn_idx = _nearest_neighbors(held_vec, gen_mat, k=min(k, len(gen_mat)))
            pred = cal["a"] * float(gen_mat[nn_idx, ci].mean()) + cal["b"]
            truth = float(ihc_mat[held_out_i, ci])
            cal_residuals.append(abs(pred - truth))
            # Baseline: predict the IHC cohort mean (intercept-only).
            baseline = float(np.delete(ihc_mat[:, ci], held_out_i).mean())
            intercept_residuals.append(abs(baseline - truth))
        results[cell] = {
            "loo_mae": float(np.mean(cal_residuals)),
            "intercept_mae": float(np.mean(intercept_residuals)),
        }
    return results


def predict_time_from_genomics(
    genomics_profile: TIMEProfile,
    coefs: dict[str, dict[str, float]],
) -> TIMEProfile:
    """Apply the calibration to one genomics profile -> calibrated TIMEProfile.

    Output cohort is tagged "tcga_hnsc" (the source) but modality becomes
    "calibrated_prediction" and provenance records the calibration version.
    """
    gen_vec = _profile_to_vector(genomics_profile)
    calibrated = {}
    for ci, cell in enumerate(CELL_TYPES):
        cal = coefs[cell]
        v = cal["a"] * float(gen_vec[ci]) + cal["b"]
        calibrated[cell] = max(0.0, min(5.0, v))

    region = RegionProfile(
        CD3_density=calibrated["CD3"],
        CD8_density=calibrated["CD8"],
        FoxP3_density=calibrated["FoxP3"],
        PanCK_density=calibrated["PanCK"],
    )
    til = derive_til_score(region)
    phen = derive_immune_phenotype(til, region.CD8_density, region.PanCK_density)
    return TIMEProfile(
        patient_id=genomics_profile.patient_id,
        cohort=genomics_profile.cohort,
        modality="calibrated_prediction",
        regions={"tumor_core": region, "tumor_margin": None, "adjacent_stroma": None},
        TIL_score=til,
        immune_phenotype=phen,
        provenance=Provenance(
            method=CALIBRATION_METHOD,
            version=CALIBRATION_VERSION,
        ),
    )


def run_calibration_arm(
    profiles_dir: Path,
    cohort_df_or_dir,
    out_dir: Path,
    k: int = 5,
) -> dict:
    """End-to-end Arm 3.

    Args:
        profiles_dir: where Arm 2 + Arm 1 wrote per-patient JSONs
        cohort_df_or_dir: unused in v0.3 (placeholder for v0.4 when we
                          can also stratify NN on TCGA clinical metadata)
        out_dir: where calibrated predictions are written

    Returns summary metrics for the audit ledger.
    """
    grouped = load_profiles(profiles_dir)
    ihc_profiles = grouped.get("deepliif_sample") or grouped.get("pmc10571229") or []
    gen_profiles = grouped.get("tcga_hnsc", [])

    if not gen_profiles:
        raise FileNotFoundError(
            f"no tcga_hnsc profiles in {profiles_dir} — Arm 2 must run first"
        )
    if not ihc_profiles:
        raise FileNotFoundError(
            f"no IHC profiles in {profiles_dir} — Arm 1 must run first"
        )

    coefs = fit_calibration(ihc_profiles, gen_profiles, k=k)
    loo = loo_validate(ihc_profiles, gen_profiles, k=k)

    # Write per-patient calibrated predictions for every TCGA-HNSC genomics profile
    out_dir.mkdir(parents=True, exist_ok=True)
    calibrated_dir = out_dir / "calibrated_predictions"
    calibrated_dir.mkdir(parents=True, exist_ok=True)
    n_predicted = 0
    for gp in gen_profiles:
        pred = predict_time_from_genomics(gp, coefs)
        (calibrated_dir / f"{gp.patient_id}.json").write_text(
            json.dumps(pred.to_dict(), indent=2)
        )
        n_predicted += 1

    # Save coefficients + validation report
    report = {
        "method": CALIBRATION_METHOD,
        "version": CALIBRATION_VERSION,
        "k_neighbors": k,
        "n_ihc_reference": len(ihc_profiles),
        "n_genomics_input": len(gen_profiles),
        "n_calibrated_predictions": n_predicted,
        "calibration_coefficients": coefs,
        "leave_one_out_validation": loo,
    }
    (out_dir / "calibration_report.json").write_text(json.dumps(report, indent=2))

    # Summary for audit
    mae_means = [loo[c]["loo_mae"] for c in CELL_TYPES if not np.isnan(loo[c]["loo_mae"])]
    intercept_mae_means = [loo[c]["intercept_mae"] for c in CELL_TYPES if not np.isnan(loo[c]["intercept_mae"])]
    return {
        "n_ihc_reference": len(ihc_profiles),
        "n_genomics_input": len(gen_profiles),
        "n_calibrated_predictions": n_predicted,
        "mean_loo_mae": float(np.mean(mae_means)) if mae_means else float("nan"),
        "mean_intercept_mae": float(np.mean(intercept_mae_means)) if intercept_mae_means else float("nan"),
    }
