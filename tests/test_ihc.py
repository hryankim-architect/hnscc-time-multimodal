"""Unit tests for src/hnscc_time/ihc.py.

The IHC arm depends on the `ihc` extra (cellpose + scikit-image + pillow).
Tests `importorskip` so CI without the extra still passes.
"""

from __future__ import annotations

import numpy as np
import pytest

# These are heavy optional deps; skip the whole module if missing.
pytest.importorskip("PIL")
pytest.importorskip("cellpose")

from hnscc_time import ihc
from hnscc_time.time_schema import CELL_TYPES, RegionProfile


def test_per_cell_marker_intensities_zero_on_empty_mask():
    image = np.zeros((32, 32, 3), dtype=np.uint8)
    masks = np.zeros((32, 32), dtype=np.uint32)
    out = ihc._per_cell_marker_intensities(image, masks)
    assert set(out.keys()) == set(CELL_TYPES)
    for v in out.values():
        assert v == 0.0


def test_per_cell_marker_intensities_scales_into_region_profile_range():
    # Single cell covering top-left 8x8 with full red intensity
    image = np.zeros((32, 32, 3), dtype=np.uint8)
    image[:8, :8, 0] = 255  # all red in cell area
    masks = np.zeros((32, 32), dtype=np.uint32)
    masks[:8, :8] = 1

    out = ihc._per_cell_marker_intensities(image, masks)
    # PanCK = R / 51 -> 255/51 = 5.0 (max RegionProfile-acceptable value)
    assert out["PanCK"] == pytest.approx(5.0, abs=1e-6)
    assert out["CD3"] == 0.0
    assert out["CD8"] == 0.0
    # FoxP3 = mean of CD3 + CD8 = 0
    assert out["FoxP3"] == 0.0
    # Result must build a valid RegionProfile
    region = RegionProfile(
        CD3_density=out["CD3"], CD8_density=out["CD8"],
        FoxP3_density=out["FoxP3"], PanCK_density=out["PanCK"],
    )
    assert region.PanCK_density == 5.0


def test_run_ihc_arm_missing_data_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        ihc.run_ihc_arm(tmp_path / "data", tmp_path / "out")
