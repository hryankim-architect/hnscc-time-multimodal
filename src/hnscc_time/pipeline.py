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
import urllib.error
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


def _download(url: str, dest: Path, *, timeout: float = 60.0, retries: int = 4) -> None:
    """Download ``url`` to ``dest`` with a per-attempt timeout + retries."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                dest.write_bytes(resp.read())
            return
        except (urllib.error.URLError, TimeoutError, OSError) as exc:  # noqa: PERF203
            last = exc
            time.sleep(min(2**attempt, 8))
    raise RuntimeError(f"failed to download {url} after {retries} attempts: {last}")


def _star_meta(rel: str, url: str) -> tuple[str, str, str]:
    """(file_id, case_submitter_id, file_name) from a STAR input url + path.

    Path layout is ``tcga_hnsc/star_counts/<case>__<file_name>``; url ends in the
    GDC ``file_id``.
    """
    name = Path(rel).name
    case, _, file_name = name.partition("__")
    file_id = url.rsplit("/", 1)[1]
    return file_id, case, file_name


def write_subset_manifest(rows: list[tuple[str, str, str]], out_dir: Path) -> Path:
    """Write ``tcga_hnsc/_subset_manifest.tsv`` (file_id/case/file_name).

    This is the cohort index ``cohort.load_cohort`` requires; building it from the
    manifest's STAR ``inputs`` is what lets ``make data`` (not just the shell
    download script) produce the layout ``make run`` reads.
    """
    sm = out_dir / "tcga_hnsc" / "_subset_manifest.tsv"
    sm.parent.mkdir(parents=True, exist_ok=True)
    body = "file_id\tcase_submitter_id\tfile_name\n" + "".join(
        f"{fid}\t{case}\t{fn}\n" for fid, case, fn in sorted(rows)
    )
    sm.write_text(body, encoding="utf-8")
    return sm


def fetch_clinical(clinical: dict[str, Any], case_ids: list[str], out_dir: Path) -> dict[str, Any]:
    """Fetch + canonicalize + verify the clinical TSV for ``case_ids``.

    POSTs the GDC ``/cases`` query in the manifest's ``clinical`` block, then
    canonicalizes to ``header + rows sorted`` (byte-stable) before writing, so the
    on-disk sha256 matches the pinned value across runs.
    """
    rel = clinical["path"]
    dest = out_dir / rel
    expected = clinical.get("sha256")
    if dest.exists() and expected and _checksum(dest) == expected:
        return {"path": str(dest), "status": "cached"}

    payload = {
        "filters": {"op": "in", "content": {"field": "submitter_id", "value": sorted(case_ids)}},
        "fields": ",".join(clinical.get("fields", [])),
        "format": "TSV",
        "size": str(len(case_ids) + 10),
    }
    req = urllib.request.Request(
        clinical["endpoint"],
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        raw = resp.read()
    lines = raw.decode("utf-8").rstrip("\n").split("\n")
    canon = ("\n".join([lines[0]] + sorted(lines[1:])) + "\n").encode("utf-8")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(canon)
    actual = _checksum(dest)
    status = "downloaded" if (not expected or actual == expected) else "checksum_mismatch"
    return {"path": str(dest), "status": status, "sha256": actual}


def _hpv_status(case: dict[str, Any]) -> str:
    """Map a case's molecular HPV test(s) to positive / negative / unknown."""
    results: set[str] = set()
    for fu in case.get("follow_ups") or []:
        for mt in fu.get("molecular_tests") or []:
            if str(mt.get("laboratory_test", "")).lower() == "human papillomavirus":
                r = str(mt.get("test_result", "")).lower()
                if r in ("positive", "amplified"):
                    results.add("positive")
                elif r in ("negative", "not amplified"):
                    results.add("negative")
    if "positive" in results:
        return "positive"
    if "negative" in results:
        return "negative"
    return "unknown"


def _case_survival(case: dict[str, Any]) -> tuple[float, int] | tuple[None, None]:
    dem = case.get("demographic") or {}
    dx = (case.get("diagnoses") or [{}])[0]
    vs = dem.get("vital_status")
    if vs == "Dead" and dem.get("days_to_death") is not None:
        return float(dem["days_to_death"]), 1
    if vs == "Alive" and dx.get("days_to_last_follow_up") is not None:
        return float(dx["days_to_last_follow_up"]), 0
    return None, None


def fetch_hpv(hpv: dict[str, Any], out_dir: Path) -> dict[str, Any]:
    """Fetch + canonicalize + verify the HPV-status + survival table (Arm 4).

    Queries GDC for the project's molecular HPV tests, joins overall survival,
    keeps cases with a known HPV result + usable survival, and writes a tidy TSV
    (rows sorted by submitter_id) so the on-disk sha256 matches the pinned value.
    """
    rel = hpv["path"]
    dest = out_dir / rel
    expected = hpv.get("sha256")
    if dest.exists() and expected and _checksum(dest) == expected:
        return {"path": str(dest), "status": "cached"}

    payload = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": "cases.project.project_id",
                                         "value": [hpv["project"]]}},
                {"op": "in", "content": {"field": "follow_ups.molecular_tests.laboratory_test",
                                         "value": [hpv["laboratory_test"]]}},
            ],
        },
        "size": "600",
        "expand": "follow_ups.molecular_tests,demographic,diagnoses",
        "fields": "submitter_id",
    }
    req = urllib.request.Request(
        hpv["endpoint"],
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=90) as resp:
        hits = json.loads(resp.read())["data"]["hits"]

    rows: list[tuple[str, str, str, str]] = []
    for case in hits:
        status = _hpv_status(case)
        t, e = _case_survival(case)
        if status in ("positive", "negative") and t is not None:
            rows.append((case["submitter_id"], status, f"{t:.1f}", str(e)))
    rows.sort()
    body = "submitter_id\thpv_status\tos_days\tos_event\n" + "".join(
        "\t".join(r) + "\n" for r in rows
    )
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(body, encoding="utf-8")
    actual = _checksum(dest)
    status = "downloaded" if (not expected or actual == expected) else "checksum_mismatch"
    return {"path": str(dest), "status": status, "sha256": actual, "n_cases": len(rows)}


def fetch_manifest(manifest_path: Path, out_dir: Path) -> dict[str, Any]:
    """Fetch every public input and reproduce the layout ``make run`` reads.

    Downloads the STAR ``inputs`` (sha256-verified), derives
    ``tcga_hnsc/_subset_manifest.tsv`` from them, and fetches the ``clinical``
    block — so ``make data && make run`` is reproducible from this manifest
    alone, no shell download script required.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("r", encoding="utf-8") as fh:
        manifest = yaml.safe_load(fh) or {}

    results: list[dict[str, Any]] = []
    rows: list[tuple[str, str, str]] = []
    for entry in manifest.get("inputs", []):
        url = entry["url"]
        rel = entry["path"]
        expected = entry.get("sha256")
        dest = out_dir / rel

        if dest.exists() and expected and _checksum(dest) == expected:
            results.append({"path": str(dest), "status": "cached"})
        else:
            _download(url, dest)
            actual = _checksum(dest)
            if expected and actual != expected:
                results.append({
                    "path": str(dest),
                    "status": "checksum_mismatch",
                    "expected": expected,
                    "actual": actual,
                })
            else:
                results.append({"path": str(dest), "status": "downloaded", "sha256": actual})
        rows.append(_star_meta(rel, url))

    out: dict[str, Any] = {"inputs": results}

    if rows:
        out["subset_manifest"] = str(write_subset_manifest(rows, out_dir))

    clinical = manifest.get("clinical")
    if clinical:
        case_ids = [case for _, case, _ in rows]
        out["clinical"] = fetch_clinical(clinical, case_ids, out_dir)

    hpv = manifest.get("hpv")
    if hpv:
        out["hpv"] = fetch_hpv(hpv, out_dir)

    return out


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

        # ---- Arm 4 (HPV± overall survival) ---------------------------
        # Clinical-survival stratifier; independent of the genomics/IHC arms.
        # Skips gracefully if hpv_status.tsv was not fetched.
        arm4_t0 = time.time()
        try:
            from hnscc_time import hpv as hpv_mod

            arm4 = hpv_mod.run_hpv_arm(data_dir, out_dir / "hpv")
            audit.emit(
                action="hpv.survival.computed",
                job_id=job_id,
                fields={
                    k: arm4[k]
                    for k in ("n_hpv_positive", "n_hpv_negative", "cox_hr_hpv_pos_vs_neg",
                              "logrank_p", "direction_protective", "significant_p05")
                    if k in arm4
                },
            )
            for k, v in arm4.items():
                if isinstance(v, int | float | bool):
                    metrics[f"arm4_{k}"] = float(v)
            arm_summaries["arm4_hpv_survival"] = arm4
        except FileNotFoundError as exc:
            audit.emit(
                action="arm4_skipped_missing_data",
                job_id=job_id,
                fields={"reason": str(exc)},
            )
            arm_summaries["arm4_hpv_survival"] = {"skipped": str(exc)}
        except Exception as exc:  # noqa: BLE001 — intentional arm isolation
            audit.emit(
                action="arm4_failed",
                job_id=job_id,
                fields={"reason": f"{type(exc).__name__}: {exc}"},
            )
            arm_summaries["arm4_hpv_survival"] = {"failed": f"{type(exc).__name__}: {exc}"}
        metrics["arm4_elapsed_ms"] = (time.time() - arm4_t0) * 1000.0

        # ---- Arm 5 (R-bridge deconvolution cross-check) --------------
        # Optional Python<->R bridge: skips cleanly when R is absent (CI/sandbox)
        # or when the cohort has no STAR data.
        arm5_t0 = time.time()
        try:
            from hnscc_time import cohort as cohort_mod
            from hnscc_time import deconv_r as deconv_mod

            cohort_df5 = cohort_mod.load_cohort(data_dir)
            arm5 = deconv_mod.run_deconv_arm(cohort_df5, out_dir / "deconv_r")
            audit.emit(
                action="deconv.r_bridge.cross_method",
                job_id=job_id,
                fields={k: arm5[k] for k in ("n_samples", "mean_spearman") if k in arm5},
            )
            if isinstance(arm5.get("mean_spearman"), int | float):
                metrics["arm5_mean_spearman"] = float(arm5["mean_spearman"])
            arm_summaries["arm5_deconv_r"] = arm5
        except deconv_mod.RDeconvUnavailable as exc:
            audit.emit(
                action="arm5_skipped_no_R",
                job_id=job_id,
                fields={"reason": str(exc)},
            )
            arm_summaries["arm5_deconv_r"] = {"skipped": "R not available"}
        except FileNotFoundError as exc:
            audit.emit(
                action="arm5_skipped_missing_data",
                job_id=job_id,
                fields={"reason": str(exc)},
            )
            arm_summaries["arm5_deconv_r"] = {"skipped": str(exc)}
        except Exception as exc:  # noqa: BLE001 — intentional arm isolation
            audit.emit(
                action="arm5_failed",
                job_id=job_id,
                fields={"reason": f"{type(exc).__name__}: {exc}"},
            )
            arm_summaries["arm5_deconv_r"] = {"failed": f"{type(exc).__name__}: {exc}"}
        metrics["arm5_elapsed_ms"] = (time.time() - arm5_t0) * 1000.0

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
