"""
Helpers to prepare the surface boundary conditions of a CATHY project (irrigation, ...).

The surface nodes of the 3-d mesh are those of the DEM grid touched by at least one active
cell, numbered row by row from the northern row, west to east (see ``triangoli.f``).
"""
import numpy as np


def surface_node_areas(dem, dx, dy, nodata=(0, -9999)):
    """
    Area attached to every surface node of the mesh (lumped, as in the processor).

    Every active cell is split in two triangles by its NW-SE diagonal and each triangle gives
    one third of its area to each of its vertices. The result reproduces the ``fort.111``
    file written by the processor (checked on the Venezia and Salerno cases: relative error
    below 1e-13 for all nodes).

    Parameters
    ----------
    dem : 2d array (nrow, ncol)
        DEM, row 0 = northern row. Cells equal to one of ``nodata`` are inactive.
    dx, dy : float
        Cell size along x (columns) and y (rows).
    nodata : tuple, optional
        Values marking cells outside the catchment (0 for the processor, -9999 for users).

    Returns
    -------
    areas : 1d array
        Area of the surface nodes, in the order of the nodes of the mesh (``output/xyz``).
    nodes : 2d int array (n_nodes, 2)
        (row, column) index of each node on the (nrow+1, ncol+1) grid of nodes, row 0 = north.
    """
    dem = np.asarray(dem)
    nrow, ncol = dem.shape
    active = ~np.isin(dem, nodata)
    area_grid = np.zeros((nrow + 1, ncol + 1))
    touched = np.zeros((nrow + 1, ncol + 1), dtype=bool)
    third = dx * dy / 2.0 / 3.0
    for r, c in zip(*np.nonzero(active)):
        nw, ne, sw, se = (r, c), (r, c + 1), (r + 1, c), (r + 1, c + 1)
        for tri in ((sw, se, nw), (se, ne, nw)):
            for v in tri:
                area_grid[v] += third
        for v in (nw, ne, sw, se):
            touched[v] = True
    nodes = np.argwhere(touched)  # row-major: north to south, west to east
    return area_grid[touched], nodes


# --------------------------------------------------------------------------- #
# irrigation_map: input of the irrigation module of the processor (unit IINDEM12)
# --------------------------------------------------------------------------- #
def read_irrigation_map(path, nrow, ncol):
    """
    Read an ``other-dem/irrigation_map`` file (format of ``irr_read_one.f`` / ``irr_read_nxt.f``).

    The file holds ``IRR_FLAG``, ``IRR_SW`` and ``IRR_HS`` and then, repeated for every time
    where the rates change, a time value followed by a raster of irrigation rates [L/T]
    (matrix convention, row 0 = north) if ``IRR_HS = 0``, or one single rate if ``IRR_HS = 1``.
    Rates are piecewise constant in time.

    Returns
    -------
    dict with ``irr_flag``, ``irr_sw``, ``irr_hs``, ``times`` (n,) and ``rates``
    ((n, nrow, ncol), or (n,) when ``irr_hs == 1``).
    """
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f.read().splitlines()]
    lines = [ln for ln in lines if ln]
    irr_flag, irr_sw, irr_hs = (float(lines[k].split()[0]) for k in range(3))
    times, rates, i = [], [], 3
    while i < len(lines):
        times.append(float(lines[i].split()[0]))
        i += 1
        if int(irr_hs) == 0:
            rates.append(np.array([[float(x) for x in lines[i + r].split()[:ncol]] for r in range(nrow)]))
            i += nrow
        else:
            rates.append(float(lines[i].split()[0]))
            i += 1
    return {
        "irr_flag": int(irr_flag),
        "irr_sw": irr_sw,
        "irr_hs": int(irr_hs),
        "times": np.array(times),
        "rates": np.array(rates),
    }


def write_irrigation_map(path, times, rates, irr_sw, irr_flag=1):
    """
    Write an ``other-dem/irrigation_map`` file.

    Parameters
    ----------
    times : 1d array
        Times [s] from which each rate applies (increasing, the first one is 0).
    rates : array (n, nrow, ncol) or (n,)
        Irrigation rates [L/T] (m/s) per time, rasters with row 0 = north (spatially variable,
        ``IRR_HS = 0``) or one value per time (homogeneous, ``IRR_HS = 1``).
    irr_sw : float
        Water saturation threshold: irrigation is applied at a surface node only when its
        saturation is below ``irr_sw``.
    irr_flag : int
        1 activates the module, 0 skips all irrigation processing.
    """
    times = np.asarray(times, dtype=float)
    rates = np.asarray(rates, dtype=float)
    if np.any(np.diff(times) <= 0):
        raise ValueError("irrigation times must be strictly increasing")
    spatial = rates.ndim == 3
    if len(times) != len(rates):
        raise ValueError("one rate (map) per time is required")

    def fmt(x):
        return "0" if x == 0 else "%.8g" % x

    with open(path, "w") as f:
        f.write("%d         IRR_FLAG (=1 for irrigation module, =0 to skip all irrigation input and processing)\n" % irr_flag)
        f.write("%-9s IRR_SW   (threshold S_w value for application or not of available irrigation water)\n" % fmt(irr_sw))
        f.write("%d         IRR_HS   (=0 for spatially variable irrigation input data, =1 for homogeneous)\n" % (0 if spatial else 1))
        for t, r in zip(times, rates):
            f.write("%.1f\n" % t)
            if spatial:
                for row in r:
                    f.write(" ".join(fmt(v) for v in row) + "\n")
            else:
                f.write(fmt(r) + "\n")


# --------------------------------------------------------------------------- #
# irrigation_map generator from field schedules
# --------------------------------------------------------------------------- #
def build_irrigation_map(fields, shape, t_end, days_per_year=365):
    """
    Build the (times, rates) of an ``irrigation_map`` from per-field irrigation schedules.

    Every field irrigates during an *event*: ``duration_h`` hours at ``rate_mm_h`` starting at
    ``start_hour`` of every ``interval_days``-th day between ``start_doy`` and ``end_doy``, each
    year. The raster at any time is the sum of the rates of the fields being irrigated; the
    result is piecewise constant and compressed (a new time is only written when the raster
    changes), as expected by ``write_irrigation_map``.

    Parameters
    ----------
    fields : list of dict
        Keys: ``mask`` (bool array ``shape``, True = field cells), ``start_doy``, ``end_doy``,
        ``interval_days``, ``start_hour``, ``duration_h``, ``rate_mm_h``.
    shape : (nrow, ncol)
    t_end : float
        End of the simulation [s]; events starting after it are ignored.
    days_per_year : int

    Returns
    -------
    times : 1d array [s], starts at 0
    rates : array (n, nrow, ncol) in m/s
    """
    day = 86400.0
    events = []  # (t_on, t_off, field index)
    n_years = int(np.ceil(t_end / (days_per_year * day)))
    for k, fld in enumerate(fields):
        for year in range(n_years):
            for doy in range(int(fld["start_doy"]), int(fld["end_doy"]) + 1, int(fld["interval_days"])):
                t_on = (year * days_per_year + doy - 1) * day + fld["start_hour"] * 3600.0
                t_off = t_on + fld["duration_h"] * 3600.0
                if t_on < t_end:
                    events.append((t_on, min(t_off, t_end), k))
    breaks = sorted({0.0} | {e[0] for e in events} | {e[1] for e in events})
    rates_m_s = [np.asarray(f["mask"], dtype=bool) * f["rate_mm_h"] / 1000.0 / 3600.0 for f in fields]
    times, maps = [], []
    for t in breaks:
        raster = np.zeros(shape)
        for t_on, t_off, k in events:
            if t_on <= t < t_off:
                raster += rates_m_s[k]
        if not maps or not np.array_equal(raster, maps[-1]):
            times.append(t)
            maps.append(raster)
    return np.array(times), np.array(maps)


def cell_raster_to_node_rates(rates, dx, dy, nodata=(0, -9999), dem=None):
    """
    Volume conserving transfer of cell rates [L/T] to surface node rates (same lumping as
    ``surface_node_areas``: NW and SE nodes get 1/3 of the cell area, NE and SW nodes 1/6).

    Used to express a rate given per DEM cell (e.g. an irrigation raster) as the per-node
    values of a spatially variable ``atmbc`` (HSPATM=0): ``sum(node_rate * node_area)`` equals
    ``sum(cell_rate * dx * dy)``.

    Parameters
    ----------
    rates : array (nrow, ncol) or (nt, nrow, ncol)
    dem : array (nrow, ncol), optional
        DEM defining the active cells (default: all cells active).

    Returns
    -------
    node_rates : array (n_nodes,) or (nt, n_nodes), order of ``surface_node_areas``
    """
    rates = np.asarray(rates, dtype=float)
    single = rates.ndim == 2
    if single:
        rates = rates[None]
    nt, nrow, ncol = rates.shape
    active = np.ones((nrow, ncol), bool) if dem is None else ~np.isin(np.asarray(dem), nodata)
    area, nodes = surface_node_areas(np.ones((nrow, ncol)) if dem is None else dem, dx, dy, nodata=nodata)
    index = -np.ones((nrow + 1, ncol + 1), dtype=int)
    index[nodes[:, 0], nodes[:, 1]] = np.arange(len(nodes))
    weights = np.zeros((len(nodes), nrow * ncol))
    cell_area = dx * dy
    for r, c in zip(*np.nonzero(active)):
        k = r * ncol + c
        weights[index[r, c], k] += cell_area / 3.0          # NW
        weights[index[r + 1, c + 1], k] += cell_area / 3.0  # SE
        weights[index[r, c + 1], k] += cell_area / 6.0      # NE
        weights[index[r + 1, c], k] += cell_area / 6.0      # SW
    node_rates = rates.reshape(nt, nrow * ncol) @ weights.T / area
    return node_rates[0] if single else node_rates
