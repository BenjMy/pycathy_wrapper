"""Bounds of the updated parameters (A2), bounded ic prior (A3), sorted ens_valid (B6)."""
from types import SimpleNamespace

import numpy as np

from pyCATHY.DA import perturbate
from pyCATHY.DA.cathy_DA import DA, _nonphysical_values


def _fake_da(entry, ens_valid=(3, 1, 2), count=1):
    """Minimal stand-in for a DA object (only what the tested methods use)."""
    da = SimpleNamespace(
        dict_parm_pert={"ZROOT0": dict(entry)},
        count_DA_cycle=count,
        ens_valid=list(ens_valid),
        NENS=4,
    )
    da._clip_updated_parameter = lambda n, v: DA._clip_updated_parameter(da, n, v)
    da._check_after_analysis = lambda *a: {0, 1, 2, 3}
    return da


def test_negative_root_depth_is_nonphysical():
    mask, label = _nonphysical_values("ZROOT0", np.array([-0.2, 0.4, 20.0, np.nan]), lambda: 5.0)
    assert mask.tolist() == [True, False, True, False]  # NaN is not flagged
    assert "root depth" in label


def test_clip_after_analysis_keeps_members_and_nan():
    da = _fake_da({"clip_min": 0, "clip_max": 1.5})
    out = DA._clip_updated_parameter(da, "ZROOT0", [-5.6, 0.3, 17.5, np.nan])
    np.testing.assert_allclose(out[:3], [0, 0.3, 1.5])
    assert np.isnan(out[3])


def test_no_bounds_no_clip():
    da = _fake_da({"clip_min": None, "clip_max": None})
    out = DA._clip_updated_parameter(da, "ZROOT0", [-5.6, 17.5])
    np.testing.assert_allclose(out, [-5.6, 17.5])


def test_update_pert_parm_dict_clips_and_sorts():
    da = _fake_da({"clip_min": 0, "clip_max": 1.5})
    DA.update_pert_parm_dict(da, "update_nb1", ["St. var.", "ZROOT0"],
                             np.array([[-1.0, 0.4, 9.0, 0.5]]))
    vals = da.dict_parm_pert["ZROOT0"]["update_nb1"]
    assert vals.min() >= 0 and vals.max() <= 1.5
    assert da.ens_valid == sorted(da.ens_valid)


def test_ic_prior_respects_per_bounds():
    """A bounded ic prior has no positive pressure head (they crash CATHY with IPOND=0)."""
    parm = {
        "type_parm": "ic", "nominal": -1.0, "mean": -1.0, "sd": 1.75, "sampling_type": "normal",
        "ensemble_size": 500, "per_type": None, "clip_min": -5.0, "clip_max": -0.05,
    }
    out = perturbate.perturbate_parm({}, parm=parm, type_parm="ic", mean=-1.0, sd=1.75,
                                     sampling_type="normal", ensemble_size=500,
                                     per_type=None, savefig=False)
    samples = np.asarray(out["ic"]["sampling"])
    assert samples.max() <= -0.05 and samples.min() >= -5.0


def test_zone_index_ignores_nodata():
    """loci must index the zones without -9999 (it indexes mesh3d_nodes_valid)."""
    import numpy as np
    from pyCATHY.DA import localisation
    from pyCATHY.DA.cathy_DA import _valid_zone_values

    ys, xs = np.meshgrid(np.arange(3), np.arange(3), indexing="ij")
    pts = [[y, x, z] for z in (0, 1) for y, x in zip(ys.ravel(), xs.ravel())]
    grid3d = {"mesh3d_nodes": np.array(pts, float), "nnod": len(pts)}
    raster = np.array([[-9999, 1, 1], [-9999, 2, 2], [-9999, 2, 2]])
    groups = localisation.mesh_nodes_local_valid(raster, grid3d)
    sizes = {int(v): len(groups[i]) for i, v in enumerate(_valid_zone_values(raster))}
    assert sizes == {1: 4, 2: 8}  # 2 surface nodes x 2 depths, 4 x 2 depths


def test_feddes_update_is_not_lost_under_copy_on_write():
    """B9: the chained assignment `df[col].iloc[i] = v` does nothing with pandas >= 3."""
    import pandas as pd
    from pyCATHY.DA.cathy_DA import _set_feddes_value

    feddes = pd.DataFrame({"ZROOT": [0.4, 0.6, 0.8], "PCREF": [-4.0, -4.0, -4.0]}, index=[1, 2, 3])
    other = feddes.copy()
    _set_feddes_value(feddes, "ZROOT", 1, 9.9)
    assert feddes["ZROOT"].tolist() == [0.4, 9.9, 0.8]  # written
    assert feddes["PCREF"].equals(other["PCREF"])        # other columns untouched
    assert feddes.index.tolist() == [1, 2, 3]            # position, not label
