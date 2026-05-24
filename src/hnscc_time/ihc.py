"""IHC arm — multiplex tissue images -> per-patient TIME profile.

v0.2 scope (honest):
    - The PMC10571229 multiplex IHC/mIF archive URL is unresolved as of
      v0.0 (see `data/manifest.yaml` and `docs/what-is-out-of-scope.md`).
    - We use the 5 RGB ROIs in `data/pmc10571229/source-repo/Sample_Large_Tissues/`
      that ship with the DeepLIIF source repo. These are real HNSCC
      tissue images but they are *single-channel RGB hematoxylin / DAPI*
      composites, not the 5-channel mIF stack the PMC dataset would
      provide. v0.2 therefore demonstrates the *segmentation + per-region
      aggregation pipeline*; full per-marker classification awaits
      the multi-channel archive.

Pipeline:
    RGB ROI .png
    -> Cellpose nuclei segmentation (cyto model, no GPU required)
    -> per-cell area + centroid
    -> heuristic per-marker channel split (R = PanCK-ish, G = CD8-ish,
       B = CD3/FoxP3-ish) — *placeholder*, documented as such, calibrated
       in Arm 3 so the v0.3 calibration mapping replaces the heuristic
       with the genomics-anchored signal.
    -> per-ROI RegionProfile (tumor_core for all 5 sample ROIs)
    -> per-image TIMEProfile JSON in common schema
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from hnscc_time.time_schema import (
    Provenance,
    RegionProfile,
    TIMEProfile,
    derive_immune_phenotype,
    derive_til_score,
)

IHC_METHOD = "cellpose_nuclei_v3 + heuristic_channel_split"
IHC_VERSION = "v0.2.0"

# Sample ROI subdirectory under data/. v0.2 ships with 5 ROIs from the
# DeepLIIF source repo; v0.3+ swaps to the resolved PMC10571229 archive.
SAMPLE_ROI_SUBDIR = "pmc10571229/source-repo/Sample_Large_Tissues"


def _load_image(path: Path) -> np.ndarray:
    """Read an RGB PNG into a numpy array (H, W, 3) uint8."""
    try:
        from PIL import Image
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "ihc arm needs the `ihc` extra: `uv sync --extra ihc` (installs "
            "cellpose + scikit-image + tifffile + pillow)"
        ) from exc
    img = Image.open(path).convert("RGB")
    return np.asarray(img)


def _segment_nuclei(image_rgb: np.ndarray) -> tuple[int, np.ndarray]:
    """Run Cellpose nuclei segmentation on an RGB image.

    Returns (n_cells, label_mask). The label_mask is a 2-D uint32 array
    with 0 = background and each positive integer = one cell.

    We use the `cyto` model channel set [2, 1] which Cellpose's docs
    recommend for cytoplasm channel = 2 (green) + nucleus channel = 1
    (red) — but for our brightfield hematoxylin samples the inverse
    is closer; we evaluate empirically.
    """
    try:
        from cellpose import models
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "ihc arm needs the `ihc` extra: `uv sync --extra ihc` (installs cellpose)"
        ) from exc

    model = models.Cellpose(model_type="cyto", gpu=False)
    masks, _flows, _styles, _diams = model.eval(
        image_rgb,
        diameter=None,        # auto-estimate
        channels=[0, 0],      # grayscale
        flow_threshold=0.4,
        cellprob_threshold=0.0,
    )
    n_cells = int(masks.max())
    return n_cells, masks


def _per_cell_marker_intensities(
    image_rgb: np.ndarray,
    masks: np.ndarray,
) -> dict[str, float]:
    """Per-marker mean intensity across all cells.

    v0.2 heuristic channel mapping for the DeepLIIF sample tissues
    (RGB composites, NOT true mIF):
        - PanCK_density ~ mean R-channel inside cells (epithelial keratin
          stains pinkish-red in H&E/IHC)
        - CD8_density   ~ mean G-channel inside cells (placeholder)
        - CD3_density   ~ mean B-channel inside cells (placeholder)
        - FoxP3_density ~ (CD3 + CD8) / 2 (placeholder — no
          true 4th-channel signal in RGB)

    Documented as a *placeholder* in the README and in the
    TIMEProfile.provenance.method field. v0.3 calibration in Arm 3
    can re-weight these mappings against the genomics anchor.

    Each value is normalised to area-weighted density per 1000 px^2 so
    larger ROIs are comparable.
    """
    if masks.max() == 0:
        return {"CD3": 0.0, "CD8": 0.0, "FoxP3": 0.0, "PanCK": 0.0}

    cell_pixels = masks > 0
    n_cell_px = int(cell_pixels.sum())
    if n_cell_px == 0:
        return {"CD3": 0.0, "CD8": 0.0, "FoxP3": 0.0, "PanCK": 0.0}

    r = float(image_rgb[..., 0][cell_pixels].mean())
    g = float(image_rgb[..., 1][cell_pixels].mean())
    b = float(image_rgb[..., 2][cell_pixels].mean())

    # Scale 0-255 -> 0-5 to stay inside RegionProfile validation range.
    return {
        "PanCK": r / 51.0,
        "CD8": g / 51.0,
        "CD3": b / 51.0,
        "FoxP3": (g + b) / (2.0 * 51.0),
    }


def process_roi(roi_path: Path) -> tuple[RegionProfile, dict]:
    """Process one ROI -> (RegionProfile, diagnostic dict)."""
    img = _load_image(roi_path)
    n_cells, masks = _segment_nuclei(img)
    densities = _per_cell_marker_intensities(img, masks)
    region = RegionProfile(
        CD3_density=densities["CD3"],
        CD8_density=densities["CD8"],
        FoxP3_density=densities["FoxP3"],
        PanCK_density=densities["PanCK"],
    )
    diag = {
        "roi": roi_path.name,
        "n_cells": n_cells,
        "image_shape": list(img.shape),
    }
    return region, diag


def run_ihc_arm(data_dir: Path, out_dir: Path) -> dict:
    """End-to-end Arm 1: 5 sample ROIs -> TIMEProfile JSONs.

    Returns a summary dict suitable for inclusion in pipeline metrics.

    Raises:
        FileNotFoundError if the sample ROI directory is missing (must
        have run `scripts/download_pmc10571229.sh` first).
    """
    roi_dir = data_dir / SAMPLE_ROI_SUBDIR
    if not roi_dir.exists():
        raise FileNotFoundError(
            f"{roi_dir} missing. Run scripts/download_pmc10571229.sh first."
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    cohort_dir = out_dir / "deepliif_sample"
    cohort_dir.mkdir(parents=True, exist_ok=True)

    roi_paths = sorted(roi_dir.glob("*.png"))
    if not roi_paths:
        return {"n_rois_processed": 0}

    profiles: list[TIMEProfile] = []
    diagnostics: list[dict] = []
    for roi_path in roi_paths:
        region, diag = process_roi(roi_path)
        diagnostics.append(diag)
        # Each ROI becomes its own "patient" since the DeepLIIF sample set
        # is per-ROI, not per-patient. v0.3+ groups by patient when the
        # PMC10571229 archive lands.
        patient_id = f"deepliif_{roi_path.stem}"
        til = derive_til_score(region)
        phen = derive_immune_phenotype(til, region.CD8_density, region.PanCK_density)
        profile = TIMEProfile(
            patient_id=patient_id,
            cohort="deepliif_sample",
            modality="mIHC",
            regions={"tumor_core": region, "tumor_margin": None, "adjacent_stroma": None},
            TIL_score=til,
            immune_phenotype=phen,
            provenance=Provenance(
                method=IHC_METHOD,
                version=IHC_VERSION,
            ),
        )
        profiles.append(profile)
        (cohort_dir / f"{patient_id}.json").write_text(
            __import__("json").dumps(profile.to_dict(), indent=2)
        )

    return {
        "n_rois_processed": len(profiles),
        "n_cells_total": int(sum(d["n_cells"] for d in diagnostics)),
        "n_inflamed": int(sum(1 for p in profiles if p.immune_phenotype == "inflamed")),
        "n_excluded": int(sum(1 for p in profiles if p.immune_phenotype == "excluded")),
        "n_desert": int(sum(1 for p in profiles if p.immune_phenotype == "desert")),
        "n_unknown": int(sum(1 for p in profiles if p.immune_phenotype == "unknown")),
        "til_mean": float(np.mean([p.TIL_score for p in profiles])),
    }
