"""Arm 4 — HPV± overall-survival stratification for TCGA-HNSC.

HPV status is the single strongest clinical stratifier in head-and-neck squamous
cell carcinoma: HPV-positive (largely oropharyngeal) disease has a markedly
better prognosis. This arm reproduces that signal honestly on the TCGA-HNSC
cases that carry a molecular HPV test result, using overall survival.

Data (pinned in ``data/manifest.yaml`` under ``hpv:``, fetched into
``data/tcga_hnsc/hpv_status.tsv``) is a tidy TSV with one row per case::

    submitter_id    hpv_status    os_days    os_event

``hpv_status`` is ``positive`` / ``negative`` (mapped from the GDC molecular test
``test_result``: Positive/Amplified -> positive, Negative/Not Amplified ->
negative); ``os_event`` is 1 for death, 0 for censored.

The deliverable is descriptive: at the size of the HPV-tested subset the protective
direction is expected but may not reach significance, and HPV's prognostic effect
is strongest in the oropharynx — both are reported rather than hidden.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

HPV_TSV = Path("data/tcga_hnsc/hpv_status.tsv")
REQUIRED_COLUMNS = ("submitter_id", "hpv_status", "os_days", "os_event")


def load_hpv_survival(path: Path = HPV_TSV) -> pd.DataFrame:
    """Load the pinned HPV-status + survival table.

    Raises ``FileNotFoundError`` if the table is absent so the pipeline can skip
    the arm gracefully (the same contract the other arms use).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run `make data` (or scripts/download_tcga_hnsc.sh) first."
        )
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"hpv_status.tsv missing columns: {missing}")
    df = df[df["hpv_status"].isin(["positive", "negative"])].copy()
    df["os_days"] = df["os_days"].astype(float)
    df["os_event"] = df["os_event"].astype(int)
    return df.sort_values("submitter_id").reset_index(drop=True)


def hpv_survival_summary(df: pd.DataFrame) -> dict[str, Any]:
    """KM medians + log-rank + univariate Cox HR for HPV+ vs HPV-.

    Returns a JSON-able dict. ``cox_hr`` < 1 means HPV+ is protective (lower
    hazard). ``interpretation`` states the direction and whether it clears the
    conventional p<0.05 bar at this cohort size — descriptive, not confirmatory.
    """
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test

    pos = df[df["hpv_status"] == "positive"]
    neg = df[df["hpv_status"] == "negative"]

    def _km_median(sub: pd.DataFrame) -> float | None:
        if sub.empty:
            return None
        kmf = KaplanMeierFitter().fit(sub["os_days"], sub["os_event"])
        med = kmf.median_survival_time_
        return None if (med is None or np.isinf(med)) else float(med)

    lr = logrank_test(pos["os_days"], neg["os_days"], pos["os_event"], neg["os_event"])

    cox = CoxPHFitter().fit(
        pd.DataFrame(
            {
                "os_days": df["os_days"],
                "os_event": df["os_event"],
                "hpv_pos": (df["hpv_status"] == "positive").astype(int),
            }
        ),
        duration_col="os_days",
        event_col="os_event",
    )
    hr = float(np.exp(cox.params_["hpv_pos"]))
    ci_low, ci_high = (float(x) for x in np.exp(cox.confidence_intervals_.loc["hpv_pos"].values))
    p = float(cox.summary.loc["hpv_pos", "p"])
    logrank_p = float(lr.p_value)

    protective = hr < 1.0
    significant = logrank_p < 0.05
    interpretation = (
        f"HPV+ shows the {'expected protective' if protective else 'unexpected adverse'} "
        f"direction (Cox HR {hr:.2f} vs HPV-), "
        + (
            f"significant at logrank p={logrank_p:.3f}."
            if significant
            else f"but at trend level only (logrank p={logrank_p:.3f}) in this HPV-tested "
            "subset; HPV's prognostic effect concentrates in oropharyngeal disease."
        )
    )

    return {
        "n_hpv_positive": int(len(pos)),
        "n_hpv_negative": int(len(neg)),
        "events_hpv_positive": int(pos["os_event"].sum()),
        "events_hpv_negative": int(neg["os_event"].sum()),
        "km_median_os_days_hpv_positive": _km_median(pos),
        "km_median_os_days_hpv_negative": _km_median(neg),
        "logrank_p": logrank_p,
        "cox_hr_hpv_pos_vs_neg": hr,
        "cox_hr_95ci": [ci_low, ci_high],
        "cox_p": p,
        "direction_protective": protective,
        "significant_p05": significant,
        "interpretation": interpretation,
    }


def make_hpv_km_plot(df: pd.DataFrame, out_path: Path) -> Path:
    """Write a Kaplan-Meier plot (HPV+ vs HPV-) to ``out_path``."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    fig, ax = plt.subplots(figsize=(6, 4))
    for label, sub in (("HPV+", df[df.hpv_status == "positive"]),
                       ("HPV-", df[df.hpv_status == "negative"])):
        if not sub.empty:
            KaplanMeierFitter().fit(
                sub["os_days"], sub["os_event"], label=f"{label} (n={len(sub)})"
            ).plot_survival_function(ax=ax, ci_show=False)
    ax.set_xlabel("Overall survival (days)")
    ax.set_ylabel("Survival probability")
    ax.set_title("TCGA-HNSC overall survival by HPV status")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)
    return out_path


def run_hpv_arm(data_dir: Path, out_dir: Path) -> dict[str, Any]:
    """Arm 4 entry point: load HPV survival table, summarize, write artifacts."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = load_hpv_survival(Path(data_dir) / "tcga_hnsc" / "hpv_status.tsv")
    summary = hpv_survival_summary(df)

    import json

    (out_dir / "hpv_survival.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    try:
        summary["km_plot"] = str(make_hpv_km_plot(df, out_dir / "hpv-km.png"))
    except Exception as exc:  # noqa: BLE001 — plot is optional; summary is the deliverable
        summary["km_plot_skipped"] = f"{type(exc).__name__}: {exc}"
    return summary
