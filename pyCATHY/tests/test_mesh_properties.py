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
