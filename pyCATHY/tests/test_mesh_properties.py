import numpy as np
import pytest
import pyvista as pv

from pyCATHY.mesh_properties import property_summary, resolve_location, set_mesh_property


def _mesh():
    return pv.ImageData(dimensions=(3, 3, 3)).cast_to_unstructured_grid()  # 27 pts, 8 cells


def test_scalar_broadcast_to_points():
    m = _mesh()
    assert set_mesh_property(m, "Ks", 1e-5) == "point"
    assert m.point_data["Ks"].shape == (m.n_points,) and np.all(m.point_data["Ks"] == 1e-5)


def test_auto_location_point_cell_and_errors():
    m = _mesh()
    assert set_mesh_property(m, "a", np.arange(m.n_points)) == "point"
    assert set_mesh_property(m, "b", np.arange(m.n_cells)) == "cell"
    assert "b" in m.cell_data and "a" in m.point_data
    with pytest.raises(ValueError):
        set_mesh_property(m, "c", np.arange(5))
    with pytest.raises(ValueError):
        set_mesh_property(m, "d", np.arange(m.n_points), location="cell")
    with pytest.raises(ValueError):
        resolve_location(m, [1.0], "nodes")


def test_rejects_nan_and_bad_name():
    m = _mesh()
    bad = np.ones(m.n_points); bad[3] = np.nan
    with pytest.raises(ValueError):
        set_mesh_property(m, "x", bad)
    with pytest.raises(ValueError):
        set_mesh_property(m, "", 1.0)


def test_field_legacy_and_summary(tmp_path):
    m = _mesh()
    assert set_mesh_property(m, "legacy", [1, 2, 3], location="field") == "field"
    set_mesh_property(m, "zroot", np.linspace(0.1, 0.9, m.n_points))
    s = property_summary(m)
    assert s["zroot"][0] == "point" and "legacy" not in s
    f = tmp_path / "m.vtk"; m.save(f, binary=False)
    back = pv.read(f)
    assert np.allclose(back.point_data["zroot"], m.point_data["zroot"])


# ── soil / root mapping ──────────────────────────────────────
import pandas as pd

from pyCATHY.mesh_properties import (add_root_properties, add_soil_properties,
                                     feddes_root_profile, interpolate_to_points,
                                     node_geometry, surface_to_nodes)


def _layered(nx=3, ny=2, nz=4, dz=0.5):
    """Terrain-following mesh: surface z = 10 + x*0.1, nz layers, dz thick."""
    xs, ys = np.arange(nx + 1) * 1.0, np.arange(ny + 1) * 1.0
    X, Y, L = np.meshgrid(xs, ys, np.arange(nz + 1), indexing="ij")
    Z = 10 + 0.1 * X - L * dz
    return pv.StructuredGrid(X, Y, Z).cast_to_unstructured_grid()


def test_node_geometry():
    m = _layered()
    col, depth, level = node_geometry(m)
    assert col.max() + 1 == 4 * 3 and depth.min() == 0 and np.isclose(depth.max(), 2.0)
    assert set(level) == {0, 1, 2, 3, 4} and np.allclose(depth, level * 0.5)  # slope-independent


def test_surface_to_nodes_origin():
    m = _layered()
    a = np.array([[1, 2, 3], [4, 5, 6]])                  # (Ny=2, Nx=3)
    lo = surface_to_nodes(m, a, 0.0, 0.0, 1.0, 1.0, "lower")
    up = surface_to_nodes(m, a, 0.0, 0.0, 1.0, 1.0, "upper")
    p = m.points
    i = np.where((p[:, 0] == 2) & (p[:, 1] == 0))[0][0]   # cell (iy=0, ix=2)
    assert lo[i] == 3 and up[i] == 6
    with pytest.raises(ValueError):
        surface_to_nodes(m, a, 0, 0, 1, 1, "middle")


def test_feddes_profile():
    b = feddes_root_profile([0.0, 0.5, 1.0, 2.0], [1.0] * 4, [1.0] * 4)
    assert np.isclose(b[0], 1.0) and 0 < b[1] < 1 and b[2] == 0 and b[3] == 0
    assert feddes_root_profile([0.1], [0.0], [1.0])[0] == 0.0


def test_add_root_and_soil_and_interp():
    m = _layered()
    veg = np.array([[1, 1, 2], [1, 1, 2]])
    FP = pd.DataFrame({"ZROOT": [0.4, 1.2], "PZ": [4.0, 4.0]}, index=[1, 2])
    add_root_properties(m, veg, FP, (0.0, 0.0, 1.0, 1.0))
    z = m.point_data["ZROOT"]
    assert set(np.unique(z)) == {0.4, 1.2}
    d = node_geometry(m)[1]
    assert np.all((m.point_data["root_zone"] == 1) == (d < z))
    assert np.all(m.point_data["root_beta"][d >= z] == 0)
    with pytest.raises(ValueError):
        add_root_properties(m, veg + 5, FP, (0.0, 0.0, 1.0, 1.0))

    idx = pd.MultiIndex.from_product([[1, 2], [1, 2, 3, 4]], names=["zone", "layer"])
    SPP = pd.DataFrame({"PERMX": np.arange(1, 9) * 1e-6}, index=idx)
    zone = np.array([[1, 1, 2], [1, 1, 2]])
    add_soil_properties(m, zone, SPP, (0.0, 0.0, 1.0, 1.0), nstr=4)
    layer = m.point_data["soil_layer"]
    assert layer.min() == 1 and layer.max() == 4          # bottom nodes -> last layer
    k = np.where((m.points[:, 0] == 3) & (m.points[:, 1] == 0))[0]
    assert np.allclose(sorted(m.point_data["PERMX"][k]), [5e-6, 6e-6, 7e-6, 8e-6, 8e-6])

    pts = m.points[:5] + 0.01
    lin = interpolate_to_points(m, ["ZROOT", "PERMX"], pts)
    assert np.all(np.isfinite(lin["ZROOT"])) and np.all(np.isfinite(lin["PERMX"]))
    far = interpolate_to_points(m, "veg_class", [[100.0, 100.0, 0.0]], "nearest")
    assert far["veg_class"][0] in (1.0, 2.0)
    with pytest.raises(ValueError):
        interpolate_to_points(m, "ZROOT", pts, "cubic")


# ── xarray (layer, y, x) view ────────────────────────────────
from pyCATHY.mesh_properties import mesh_to_xarray, xarray_to_mesh


def test_xarray_round_trip():
    m = _layered(nx=3, ny=2, nz=4)
    set_mesh_property(m, "ZROOT", np.linspace(0.1, 2.0, m.n_points))
    ds = mesh_to_xarray(m)
    assert ds["ZROOT"].dims == ("layer", "y", "x") and ds.sizes == {"layer": 5, "y": 3, "x": 4}
    assert float(ds["depth"].sel(layer=0).max()) == 0 and np.isclose(float(ds["depth"].max()), 2.0)
    assert np.allclose(ds["z"].sel(layer=0).values, 10 + 0.1 * ds["x"].values[None, :])   # terrain
    assert not np.isnan(ds["ZROOT"].values).any() and (ds["node_id"].values >= 0).all()
    # modify in xarray, write back: node-wise inverse
    ds["double"] = ds["ZROOT"] * 2
    written = xarray_to_mesh(ds, m)
    assert set(written) == {"ZROOT", "double"} and np.allclose(m.point_data["double"], 2 * m.point_data["ZROOT"])
    with pytest.raises(ValueError):
        xarray_to_mesh(ds.isel(x=slice(0, 2)), m)             # not this mesh
    # masked column -> NaN in the view, ids -1
    keep = ~((m.points[:, 0] == 0) & (m.points[:, 1] == 0))
    sub = m.extract_points(keep, adjacent_cells=False)
    ds2 = mesh_to_xarray(sub, names=[])
    assert (ds2["node_id"].values < 0).sum() == 5 and np.isnan(ds2["z"].values).sum() == 5
    assert sub.n_points == (ds2["node_id"].values >= 0).sum()


def test_feddes_table_accepts_update_soil_array():
    """soil_FP["FP"] is an (n_veg, 6) ndarray (PCANA PCREF PCWLT ZROOT PZ OMGC), not a DataFrame."""
    import numpy as np
    from pyCATHY.mesh_properties import _as_feddes_table
    FP = np.array([[0, 0, 0, 1.0, 1.0, 0], [0, 0, 0, 2.5, 3.0, 0]], dtype=float)
    t = _as_feddes_table(FP)
    assert list(t.index) == [1, 2] and list(t["ZROOT"]) == [1.0, 2.5] and list(t["PZ"]) == [1.0, 3.0]
    assert _as_feddes_table(t) is t
    try:
        _as_feddes_table(np.zeros((2, 3)))
        raise AssertionError("short table must raise")
    except ValueError:
        pass
