"""Arm 5 — Python <-> R bridge for a cross-method deconvolution check.

The genomics arm (Arm 2) scores immune cell types with a Python ssGSEA-style
**mean-rank** of marker genes. This arm runs an *independent* estimator in **R**
(a marker **z-score** deconvolution, different normalization) over the same TPM
matrix, brings the scores back across the language boundary, and checks the two
methods **agree** (Spearman per cell type). High agreement is a cross-method
robustness signal; the bridge itself (a deterministic `Rscript` subprocess) is the
integration capability being demonstrated.

Design notes:
- **Self-contained R**: `scripts/deconv_markers.R` uses base R only (no CRAN /
  Bioconductor install), so it runs wherever R is present.
- **Graceful skip**: if `Rscript` is not on PATH the arm raises
  :class:`RDeconvUnavailable`, which the pipeline records as a skip — exactly like
  the IHC arm skipping when the `ihc` extra is absent. Nothing here runs in CI or
  the sandbox unless R is installed.
- **No result is claimed without a run**: the agreement number is produced only
  when the arm actually executes on a host with R.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd

from hnscc_time.genomics import IMMUNE_SIGNATURES, cohort_score_matrix, load_tpm
from hnscc_time.time_schema import CELL_TYPES

R_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "deconv_markers.R"


class RDeconvUnavailable(RuntimeError):
    """Raised when Rscript (or the R script file) is unavailable."""


def r_available() -> bool:
    """True if an `Rscript` interpreter is on PATH and the R script exists."""
    return shutil.which("Rscript") is not None and R_SCRIPT.exists()


def _tpm_matrix(cohort_df: pd.DataFrame) -> pd.DataFrame:
    """Build a gene_name x submitter_id TPM matrix from the cohort's STAR files."""
    cols: dict[str, pd.Series] = {}
    for sid, row in cohort_df.iterrows():
        path = row.get("star_counts_path")
        if path is None or pd.isna(path):
            continue
        path = Path(path)
        if path.exists():
            cols[str(sid)] = load_tpm(path)
    if not cols:
        raise FileNotFoundError("no STAR counts files found for the cohort")
    return pd.DataFrame(cols).fillna(0.0)


def run_r_deconvolution(tpm_matrix: pd.DataFrame, timeout: float = 300.0) -> pd.DataFrame:
    """Run the R marker-zscore deconvolution over a TPM matrix.

    Writes the matrix to a temp TSV, invokes ``Rscript deconv_markers.R in out``,
    and reads the per-sample x cell-type scores back. Raises
    :class:`RDeconvUnavailable` if R is not present.
    """
    if not r_available():
        raise RDeconvUnavailable(
            f"Rscript not on PATH or missing {R_SCRIPT.name}; install R to run Arm 5."
        )
    with tempfile.TemporaryDirectory() as td:
        in_tsv = Path(td) / "tpm.tsv"
        out_tsv = Path(td) / "scores.tsv"
        # index_label keeps R's column-1 = gene_name explicit + unambiguous.
        tpm_matrix.to_csv(in_tsv, sep="\t", index_label="gene_name")
        proc = subprocess.run(  # noqa: S603 — fixed args, no shell
            ["Rscript", "--vanilla", str(R_SCRIPT), str(in_tsv), str(out_tsv)],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"Rscript failed (code {proc.returncode}): {proc.stderr[-500:]}")
        return pd.read_csv(out_tsv, sep="\t", index_col=0)


def cross_method_agreement(
    py_scores: pd.DataFrame, r_scores: pd.DataFrame
) -> dict[str, Any]:
    """Spearman agreement between the Python (rank) and R (z-score) scores.

    Aligns on shared samples + cell types and reports a per-cell-type rho plus the
    mean. Pure function — unit-testable offline with synthetic frames.
    """
    shared = py_scores.index.intersection(r_scores.index)
    cells = [c for c in CELL_TYPES if c in py_scores.columns and c in r_scores.columns]
    per_cell: dict[str, float | None] = {}
    for c in cells:
        a = py_scores.loc[shared, c]
        b = r_scores.loc[shared, c]
        # Spearman = Pearson on ranks; guard zero-variance.
        ar, br = a.rank(), b.rank()
        if ar.std() == 0 or br.std() == 0 or len(shared) < 3:
            per_cell[c] = None
        else:
            per_cell[c] = float(ar.corr(br))
    valid = [v for v in per_cell.values() if v is not None]
    return {
        "n_samples": int(len(shared)),
        "cell_types": cells,
        "spearman_per_cell": per_cell,
        "mean_spearman": (sum(valid) / len(valid)) if valid else None,
    }


def run_deconv_arm(cohort_df: pd.DataFrame, out_dir: Path) -> dict[str, Any]:
    """Arm 5 entry point: R deconvolution + cross-method agreement vs Arm 2."""
    if not r_available():
        raise RDeconvUnavailable("Rscript unavailable")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    tpm = _tpm_matrix(cohort_df)
    r_scores = run_r_deconvolution(tpm)
    py_scores = cohort_score_matrix(cohort_df)
    agreement = cross_method_agreement(py_scores, r_scores)

    import json

    r_scores.to_csv(out_dir / "r_deconv_scores.tsv", sep="\t")
    (out_dir / "cross_method_agreement.json").write_text(json.dumps(agreement, indent=2))
    return {
        "n_samples": agreement["n_samples"],
        "n_markers": int(sum(len(v) for v in IMMUNE_SIGNATURES.values())),
        "mean_spearman": agreement["mean_spearman"],
        "spearman_per_cell": agreement["spearman_per_cell"],
    }
