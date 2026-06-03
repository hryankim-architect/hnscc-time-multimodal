"""End-to-end pipeline entry point.

This module provides the pipeline entry-point pattern used across repos. Each repo
replaces the body of ``run_pipeline`` with the actual bioinformatics work
(e.g. P3's VCF→HRD score, P1's Nextflow orchestration, P2's QC classifier,
P4's IHC + genomics calibration), but keeps the surrounding shape::

    audit_start  →  tracking_start  →  body  →  tracking_end  →  audit_end

The body must be deterministic enough that the canary smoke test exercises
the same code path with a fixture input.
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import click
import yaml

from hnscc_time import audit, tracking


def _run_id(name: str) -> str:
    stamp = datetime.now(UTC).strftime("%Y%m%d-%H%M%S")
    return f"{name}-{stamp}"


def _checksum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_manifest(manifest_path: Path, out_dir: Path) -> dict[str, Any]:
    """Download every entry in the manifest; verify SHA-256 checksums."""
    out_dir.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("r", encoding="utf-8") as fh:
        manifest = yaml.safe_load(fh) or {}

    results: list[dict[str, Any]] = []
    for entry in manifest.get("inputs", []):
        url = entry["url"]
        rel = entry["path"]
        expected = entry.get("sha256")
        size_mb = entry.get("size_mb")
        dest = out_dir / rel
        dest.parent.mkdir(parents=True, exist_ok=True)

        if dest.exists() and expected and _checksum(dest) == expected:
            results.append({"path": str(dest), "status": "cached"})
            continue

        urllib.request.urlretrieve(url, dest)
        actual = _checksum(dest)
        if expected and actual != expected:
            results.append({
                "path": str(dest),
                "status": "checksum_mismatch",
                "expected": expected,
                "actual": actual,
            })
            continue
        results.append({
            "path": str(dest),
            "status": "downloaded",
            "sha256": actual,
            "size_mb": size_mb,
        })

    return {"inputs": results}


def run_pipeline(run_name: str, out_dir: Path, data_dir: Path | None = None) -> dict[str, Any]:
    """Three-arm multimodal integration pipeline.

    Arm 2 (Genomics) runs on real TCGA-HNSC subset data when present;
    Arm 1 (IHC) runs on the 5 DeepLIIF Sample_Large_Tissues ROIs when the
    `ihc` extra is installed; Arm 3 (Cross-cohort calibration) joins them.

    Each arm is independently optional — `run_pipeline` skips arms whose
    inputs are missing and records the skip in the audit ledger rather
    than failing the run.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    job_id = _run_id(run_name)
    data_dir = data_dir or Path("data")

    audit.emit(
        action="pipeline_start",
        job_id=job_id,
        fields={"out_dir": str(out_dir), "data_dir": str(data_dir)},
    )

    metrics: dict[str, float] = {}
    arm_summaries: dict[str, dict[str, Any]] = {}

    with tracking.run(name=job_id, experiment="hnscc_time"):
        tracking.log_params({"run_name": run_name})
        t_pipeline = time.time()

        # ---- Arm 2 (Genomics) -----------------------------------------
        arm2_t0 = time.time()
        try:
            from hnscc_time import cohort as cohort_mod
            from hnscc_time import genomics as genomics_mod

            cohort_df = cohort_mod.load_cohort(data_dir)
            audit.emit(
                action="cohort.tcga_hnsc.assembled",
                job_id=job_id,
                fields={"n_patients": int(len(cohort_df))},
            )
            arm2 = genomics_mod.run_genomics_arm(cohort_df, out_dir / "time_profiles")
            audit.emit(
                action="genomics.time_profiles.computed",
                job_id=job_id,
                fields=arm2,
            )
            for k, v in arm2.items():
                if isinstance(v, int | float):
                    metrics[f"arm2_{k}"] = float(v)
            arm_summaries["arm2_genomics"] = arm2
        except FileNotFoundError as exc:
            audit.emit(
                action="arm2_skipped_missing_data",
                job_id=job_id,
                fields={"reason": str(exc)},
            )
            arm_summaries["arm2_genomics"] = {"skipped": str(exc)}
        metrics["arm2_elapsed_ms"] = (time.time() - arm2_t0) * 1000.0

        # ---- Arm 1 (IHC) ---------------------------------------------
        # Broad exception catch is intentional: an arm-internal crash
        # (e.g. cellpose API drift, malformed image) must NOT kill the
        # subsequent arms. The audit ledger records the exception type
        # + message so the failure is debuggable from the chain alone.
        arm1_t0 = time.time()
        try:
            from hnscc_time import ihc as ihc_mod

            arm1 = ihc_mod.run_ihc_arm(data_dir, out_dir / "time_profiles")
            audit.emit(
                action="ihc.time_profiles.computed",
                job_id=job_id,
                fields=arm1,
            )
            for k, v in arm1.items():
                if isinstance(v, int | float):
                    metrics[f"arm1_{k}"] = float(v)
            arm_summaries["arm1_ihc"] = arm1
        except ImportError as exc:
            audit.emit(
                action="arm1_skipped_no_ihc_extra",
                job_id=job_id,
                fields={"reason": str(exc)},
            )
            arm_summaries["arm1_ihc"] = {"skipped": "ihc extra not installed"}
        except FileNotFoundError as exc:
            audit.emit(
                action="arm1_skipped_missing_data",
                job_id=job_id,
                fields={"reason": str(exc)},
            )
            arm_summaries["arm1_ihc"] = {"skipped": str(exc)}
        except Exception as exc:  # noqa: BLE001 — intentional arm isolation
            audit.emit(
                action="arm1_failed",
                job_id=job_id,
                fields={"reason": f"{type(exc).__name__}: {exc}"},
            )
            arm_summaries["arm1_ihc"] = {"failed": f"{type(exc).__name__}: {exc}"}
        metrics["arm1_elapsed_ms"] = (time.time() - arm1_t0) * 1000.0

        # ---- Arm 3 (Cross-cohort calibration) ------------------------
        # Broad except for the same reason: a calibration crash should
        # log + degrade, not corrupt the per-arm summary contract.
        arm3_t0 = time.time()
        try:
            from hnscc_time import calibrate as calibrate_mod

            arm3 = calibrate_mod.run_calibration_arm(
                profiles_dir=out_dir / "time_profiles",
                cohort_df_or_dir=data_dir,
                out_dir=out_dir / "calibration",
            )
            audit.emit(
                action="calibration.trained",
                job_id=job_id,
                fields=arm3,
            )
            for k, v in arm3.items():
                if isinstance(v, int | float):
                    metrics[f"arm3_{k}"] = float(v)
            arm_summaries["arm3_calibration"] = arm3
        except FileNotFoundError as exc:
            audit.emit(
                action="arm3_skipped_missing_profiles",
                job_id=job_id,
                fields={"reason": str(exc)},
            )
            arm_summaries["arm3_calibration"] = {"skipped": str(exc)}
        except ValueError as exc:
            audit.emit(
                action="arm3_skipped_insufficient_pairs",
                job_id=job_id,
                fields={"reason": str(exc)},
            )
            arm_summaries["arm3_calibration"] = {"skipped": str(exc)}
        except Exception as exc:  # noqa: BLE001 — intentional arm isolation
            audit.emit(
                action="arm3_failed",
                job_id=job_id,
                fields={"reason": f"{type(exc).__name__}: {exc}"},
            )
            arm_summaries["arm3_calibration"] = {"failed": f"{type(exc).__name__}: {exc}"}
        metrics["arm3_elapsed_ms"] = (time.time() - arm3_t0) * 1000.0

        metrics["pipeline_elapsed_ms"] = (time.time() - t_pipeline) * 1000.0
        tracking.log_metrics(metrics)

    # Write the per-run summary artifact (used by Makefile + tests)
    artifact_path = out_dir / f"{run_name}.json"
    with artifact_path.open("w", encoding="utf-8") as fh:
        json.dump(
            {"job_id": job_id, "metrics": metrics, "arms": arm_summaries},
            fh, indent=2, sort_keys=True,
        )

    audit.emit(
        action="pipeline_end",
        job_id=job_id,
        fields={"metrics": metrics, "artifact_path": str(artifact_path)},
    )

    return {
        "job_id": job_id,
        "metrics": metrics,
        "arms": arm_summaries,
        "artifact_path": str(artifact_path),
    }


@click.group()
def cli() -> None:
    """hnscc_time demonstration pipeline."""


@cli.command()
@click.option(
    "--manifest",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("data/manifest.yaml"),
)
@click.option(
    "--out",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("data"),
)
def fetch(manifest: Path, out: Path) -> None:
    """Download public inputs declared in the manifest."""
    result = fetch_manifest(manifest, out)
    click.echo(json.dumps(result, indent=2))


@cli.command()
@click.option("--name", default="demo")
@click.option(
    "--out",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("artifacts"),
)
def run(name: str, out: Path) -> None:
    """Run the end-to-end pipeline."""
    result = run_pipeline(name, out)
    click.echo(json.dumps(result, indent=2))


if __name__ == "__main__":
    cli()
