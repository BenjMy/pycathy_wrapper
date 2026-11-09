"""surface_node_areas must reproduce the fort.111 written by the processor (real cases)."""
import os

import numpy as np
import pytest

from pyCATHY.irrigation import surface_node_areas

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "node_areas")


@pytest.mark.parametrize("case,dx,n_nodes,n_cells", [
    ("venezia", 0.1, 154, 130),       # 13 x 10 cells, regular
    ("salerno", 80.0, 1692, 1536),    # 74 x 72 grid, irregular catchment (0 outside)
])
def test_node_areas_match_fort111_of_the_reference_runs(case, dx, n_nodes, n_cells):
    dem = np.loadtxt(os.path.join(DATA, f"{case}_dem"), skiprows=6)
    ref = np.loadtxt(os.path.join(DATA, f"{case}_fort.111"), ndmin=2)[0]
    areas, nodes = surface_node_areas(dem, dx, dx)
    assert len(areas) == len(ref) == n_nodes
    np.testing.assert_allclose(areas, ref, rtol=1e-10)
    np.testing.assert_allclose(areas.sum(), n_cells * dx * dx, rtol=1e-12)


def test_user_nodata_value_is_also_inactive():
    dem = np.ones((3, 4))
    masked = dem.copy()
    masked[0, :] = -9999
    a_ref, _ = surface_node_areas(dem[1:], 1.0, 1.0)
    a_masked, _ = surface_node_areas(masked, 1.0, 1.0)
    np.testing.assert_allclose(a_masked, a_ref)


def test_node_order_is_north_to_south_west_to_east():
    _, nodes = surface_node_areas(np.ones((2, 3)), 1.0, 1.0)
    assert nodes[0].tolist() == [0, 0] and nodes[3].tolist() == [0, 3] and nodes[4].tolist() == [1, 0]
    assert len(nodes) == 12


def test_cell_raster_to_node_rates_conserves_volume_and_uniform_fields():
    from pyCATHY.irrigation import cell_raster_to_node_rates

    rng = np.random.default_rng(0)
    rates = rng.random((3, 5, 6)) * 1e-6
    dem = np.ones((5, 6))
    dem[0, :2] = 0                                   # inactive cells
    rates[:, 0, :2] = 0
    areas, _ = surface_node_areas(dem, 25.0, 25.0)
    node_rates = cell_raster_to_node_rates(rates, 25.0, 25.0, dem=dem)
    assert node_rates.shape == (3, len(areas))
    np.testing.assert_allclose((node_rates * areas).sum(1), rates.sum((1, 2)) * 625.0, rtol=1e-12)
    uniform = cell_raster_to_node_rates(np.full((5, 6), 2e-7), 25.0, 25.0)
    np.testing.assert_allclose(uniform, 2e-7, rtol=1e-12)
