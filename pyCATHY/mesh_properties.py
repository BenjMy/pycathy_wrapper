"""
mesh_properties.py — attach physical properties to the CATHY pyvista mesh
=========================================================================

Single, checked entry point to write a property (soil, root, initial
condition, ...) onto the mesh that ends up in the ``.vtk`` file.

Locations
---------
``point`` : one value per mesh node (``mesh.n_points``) -> ``point_data``
``cell``  : one value per mesh cell  (``mesh.n_cells``)  -> ``cell_data``
``field`` : free-length array, stored as ``field_data`` (legacy behaviour;
            not drawn by ParaView/pyvista as a scalar)

Nothing here depends on CATHY, so it can be tested with a bare pyvista mesh.
"""
from __future__ import annotations

import numpy as np

LOCATIONS = ("auto", "point", "cell", "field")


def resolve_location(mesh, values, location: str = "auto") -> str:
    """Return 'point' / 'cell' / 'field' for ``values`` on ``mesh``."""
    if location not in LOCATIONS:
        raise ValueError(f"location must be one of {LOCATIONS}, got {location!r}")
    n = int(np.size(values))
    if location in ("point", "cell"):
        expected = mesh.n_points if location == "point" else mesh.n_cells
        if n != expected:
            raise ValueError(
                f"{n} values given for location='{location}', but the mesh has "
                f"{expected} {location}s."
            )
        return location
    if location == "field":
        return "field"
    # auto: nodes win when both counts are equal (CATHY properties are nodal)
    if n == mesh.n_points:
        return "point"
    if n == mesh.n_cells:
        return "cell"
    raise ValueError(
        f"cannot place {n} values: mesh has {mesh.n_points} points and "
        f"{mesh.n_cells} cells. Pass location='field' to store them unchecked."
    )


def set_mesh_property(mesh, name: str, values, location: str = "auto") -> str:
    """
    Write ``values`` as array ``name`` on ``mesh`` (in place).

    A scalar is broadcast to every node. Returns the location used.
    Values must be finite numbers (NaN/inf are rejected so that a bad
    property cannot silently reach the ``.vtk``).
    """
    if not name or not isinstance(name, str):
        raise ValueError("property name must be a non-empty string")
    if np.isscalar(values):
        arr = np.full(mesh.n_points, float(values))
        loc = "point"
    else:
        arr = np.asarray(values, dtype=float).ravel()
        loc = resolve_location(mesh, arr, location)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"property '{name}' contains NaN/inf values")
    if loc == "point":
        mesh.point_data[name] = arr
    elif loc == "cell":
        mesh.cell_data[name] = arr
    else:
        mesh.add_field_data(arr, name)
    return loc


def property_summary(mesh) -> dict:
    """``{name: (location, min, max)}`` for every scalar array on the mesh."""
    out = {}
    for loc, store in (("point", mesh.point_data), ("cell", mesh.cell_data)):
        for k in store.keys():
            a = np.asarray(store[k])
            if a.ndim == 1 and a.dtype.kind in "fiu" and a.size:
                out[k] = (loc, float(a.min()), float(a.max()))
    return out


# ─────────────────────────────────────────────────────────────
# Geometry helpers: columns, depth below surface, layer index
# ─────────────────────────────────────────────────────────────
def node_geometry(mesh, decimals: int = 3):
    """
    Per-node ``(column_id, depth, level)`` for a terrain-following layered mesh.

    * ``column_id`` : nodes sharing (x, y) (rounded to ``decimals``) form a column
    * ``depth``     : surface elevation of the column (max z) minus node z (m, >= 0)
    * ``level``     : 0 at the surface, increasing downwards
    """
    pts = np.asarray(mesh.points)
    xy = np.round(pts[:, :2], decimals)
    _, col = np.unique(xy, axis=0, return_inverse=True)
    col = np.asarray(col).ravel()
    ztop = np.full(col.max() + 1, -np.inf)
    np.maximum.at(ztop, col, pts[:, 2])
    depth = ztop[col] - pts[:, 2]
    level = np.zeros(len(pts), dtype=int)
    for c in range(col.max() + 1):
        idx = np.where(col == c)[0]
        level[idx[np.argsort(depth[idx])]] = np.arange(len(idx))
    return col, depth, level


def surface_to_nodes(mesh, surface_2d, x0, y0, dx, dy, origin: str = "lower"):
    """
    Look up a 2-D raster ``surface_2d`` (Ny, Nx) at every node's (x, y).

    ``origin='lower'``: row 0 is at ``y0`` (same convention as
    ``meshtools.build_mesh_dataset``); ``'upper'``: row 0 is the top row
    (``y0`` is then the lower edge of the last row, as for most GeoTIFFs).
    Indices are clipped to the raster, so nodes on the outer edge take the
    nearest cell.
    """
    if origin not in ("lower", "upper"):
        raise ValueError("origin must be 'lower' or 'upper'")
    a = np.asarray(surface_2d)
    ny, nx = a.shape
    pts = np.asarray(mesh.points)
    ix = np.clip(((pts[:, 0] - x0) / dx).astype(int), 0, nx - 1)
    iy = np.clip(((pts[:, 1] - y0) / dy).astype(int), 0, ny - 1)
    if origin == "upper":
        iy = ny - 1 - iy
    return a[iy, ix]


# ─────────────────────────────────────────────────────────────
# Root properties
# ─────────────────────────────────────────────────────────────
def feddes_root_profile(depth, zroot, pz):
    """
    CATHY root-density profile (see ``update_soil`` docs):
    ``beta = (1 - d/ZROOT) * exp(-PZ * d / ZROOT)`` for ``d < ZROOT`` else 0.
    ``depth``, ``zroot`` and ``pz`` are per-node arrays; ``zroot <= 0`` gives 0.
    """
    depth, zroot, pz = (np.asarray(a, dtype=float) for a in (depth, zroot, pz))
    ratio = np.divide(depth, zroot, out=np.full_like(depth, np.inf), where=zroot > 0)
    r = np.minimum(ratio, 1.0)          # beta is 0 at r=1, so this is exact
    beta = (1.0 - r) * np.exp(-pz * r)
    return np.where(ratio < 1.0, beta, 0.0)


FEDDES_COLUMNS = ["PCANA", "PCREF", "PCWLT", "ZROOT", "PZ", "OMGC"]


def _as_feddes_table(FP):
    """DataFrame indexed by vegetation class (1..n). Accepts the (n, 6) array that
    ``update_soil`` stores in ``soil_FP["FP"]`` (columns PCANA PCREF PCWLT ZROOT PZ OMGC)."""
    import pandas as pd
    if hasattr(FP, "index") and hasattr(FP, "columns"):
        return FP
    a = np.asarray(FP, dtype=float)
    if a.ndim != 2 or a.shape[1] < len(FEDDES_COLUMNS):
        raise ValueError(f"Feddes table must be (n_veg, >=6); got shape {a.shape}")
    return pd.DataFrame(a[:, :len(FEDDES_COLUMNS)], columns=FEDDES_COLUMNS,
                        index=np.arange(1, a.shape[0] + 1))


def add_root_properties(mesh, veg_class_2d, FP, raster, origin="lower",
                        zroot_col="ZROOT", pz_col="PZ"):
    """
    Write ``veg_class`` (1-based), ``ZROOT`` (m), ``root_zone`` (0/1) and
    ``root_beta`` (-) on the nodes.

    ``FP`` is a table indexed by vegetation class (1..nveg) with ``ZROOT``
    and ``PZ`` columns (the Feddes table). ``raster`` is ``(x0, y0, dx, dy)``.
    Classes absent from ``FP`` raise instead of silently getting 0.
    """
    FP = _as_feddes_table(FP)
    cls = surface_to_nodes(mesh, veg_class_2d, *raster, origin=origin).astype(int)
    missing = sorted(set(np.unique(cls)) - set(int(i) for i in FP.index))
    if missing:
        raise ValueError(f"vegetation classes {missing} are not in the Feddes table "
                         f"(has {list(FP.index)})")
    zr = np.asarray(FP[zroot_col].astype(float).reindex(cls).values)
    pz = np.asarray(FP[pz_col].astype(float).reindex(cls).values)
    _, depth, _ = node_geometry(mesh)
    set_mesh_property(mesh, "veg_class", cls.astype(float), "point")
    set_mesh_property(mesh, "ZROOT", zr, "point")
    set_mesh_property(mesh, "root_zone", (depth < zr).astype(float), "point")
    set_mesh_property(mesh, "root_beta", feddes_root_profile(depth, zr, pz), "point")


# ─────────────────────────────────────────────────────────────
# Soil properties
# ─────────────────────────────────────────────────────────────
def add_soil_properties(mesh, zone_2d, SPP, raster, nstr, origin="lower",
                        columns=None):
    """
    Write soil hydraulic properties (e.g. PERMX, POROS) on the nodes.

    ``SPP`` is the table indexed by ``(zone, layer)`` (both 1-based) that
    ``update_soil`` builds (``simu.soil_SPP["SPP_map"]``). A node at level
    ``k`` takes the properties of layer ``min(k, nstr-1) + 1``, i.e. the layer
    just below it (the bottom nodes take the last layer).
    """
    zone = surface_to_nodes(mesh, zone_2d, *raster, origin=origin).astype(int)
    _, _, level = node_geometry(mesh)
    layer = np.minimum(level, nstr - 1) + 1
    cols = list(columns) if columns is not None else [
        c for c in SPP.columns if c not in ("layer", "zone")]
    keys = list(zip(zone.tolist(), layer.tolist()))
    missing = sorted(set(keys) - set(SPP.index))
    if missing:
        raise ValueError(f"(zone, layer) pairs {missing[:5]} not in the SPP table")
    for c in cols:
        vals = SPP[c].astype(float).reindex(keys).values
        set_mesh_property(mesh, c, vals, "point")
    set_mesh_property(mesh, "soil_zone", zone.astype(float), "point")
    set_mesh_property(mesh, "soil_layer", layer.astype(float), "point")


# ─────────────────────────────────────────────────────────────
# Transfer to another mesh (e.g. ERT)
# ─────────────────────────────────────────────────────────────
def interpolate_to_points(mesh, names, points, method: str = "linear"):
    """
    Sample node properties ``names`` at ``points`` (n, 3), e.g. ERT cell centres.

    ``method='linear'`` uses scipy's Delaunay interpolation and falls back to
    nearest node outside the convex hull (no NaN); ``'nearest'`` is a pure
    nearest-node lookup (use it for class-like properties: veg_class,
    soil_zone, soil_layer).
    """
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

    if method not in ("linear", "nearest"):
        raise ValueError("method must be 'linear' or 'nearest'")
    src = np.asarray(mesh.points)
    pts = np.asarray(points, dtype=float)
    out = {}
    for n in ([names] if isinstance(names, str) else names):
        vals = np.asarray(mesh.point_data[n], dtype=float)
        near = NearestNDInterpolator(src, vals)(pts)
        if method == "nearest":
            out[n] = near
        else:
            lin = LinearNDInterpolator(src, vals)(pts)
            out[n] = np.where(np.isnan(lin), near, lin)
    return out


# ─────────────────────────────────────────────────────────────
# xarray view: (layer, y, x) <-> mesh nodes
# ─────────────────────────────────────────────────────────────
def mesh_to_xarray(mesh, names=None, decimals: int = 3):
    """
    Regular ``(layer, y, x)`` xarray Dataset view of the node properties.

    Nodes must lie on a horizontal grid (same x/y values in every layer,
    as in CATHY's terrain-following mesh); columns missing from the grid
    (masked DEM cells) become NaN. Coordinates ``x``, ``y``, ``layer`` (0 =
    surface); data variables: ``z(layer, y, x)`` (terrain-following node
    elevation), ``depth(layer, y, x)`` (below the local surface) and
    ``node_id(layer, y, x)`` (mesh point index, -1 where there is no node)
    plus every requested ``point_data`` array (default: all scalar ones).
    ``xarray_to_mesh`` is the exact inverse for the properties.
    """
    import xarray as xr

    pts = np.asarray(mesh.points)
    xs = np.unique(np.round(pts[:, 0], decimals))
    ys = np.unique(np.round(pts[:, 1], decimals))
    ix = np.searchsorted(xs, np.round(pts[:, 0], decimals))
    iy = np.searchsorted(ys, np.round(pts[:, 1], decimals))
    _, depth, level = node_geometry(mesh, decimals)
    nl = int(level.max()) + 1
    shape = (nl, len(ys), len(xs))
    node_id = np.full(shape, -1, dtype=int)
    node_id[level, iy, ix] = np.arange(len(pts))
    if (node_id >= 0).sum() != len(pts):
        raise ValueError("two nodes map to the same (layer, y, x): the mesh is not a "
                         "single-valued layered grid")

    def _fill(vals):
        out = np.full(shape, np.nan)
        out[level, iy, ix] = np.asarray(vals, dtype=float)
        return out

    dims = ("layer", "y", "x")
    data = {"z": (dims, _fill(pts[:, 2])), "depth": (dims, _fill(depth)),
            "node_id": (dims, node_id)}
    keys = list(names) if names is not None else [
        k for k in mesh.point_data.keys()
        if np.asarray(mesh.point_data[k]).ndim == 1
        and np.asarray(mesh.point_data[k]).dtype.kind in "fiu"
        and len(mesh.point_data[k]) == mesh.n_points]
    for k in keys:
        data[k] = (dims, _fill(mesh.point_data[k]))
    return xr.Dataset(data, coords={"layer": np.arange(nl), "y": ys, "x": xs})


def xarray_to_mesh(ds, mesh, names=None) -> list:
    """
    Write the ``(layer, y, x)`` variables of ``ds`` (made by ``mesh_to_xarray``,
    possibly modified or extended) back onto ``mesh`` as point data. Returns the
    names written. NaN values at existing nodes are rejected.
    """
    nid = np.asarray(ds["node_id"].values)
    if nid.max() + 1 != mesh.n_points or (nid >= 0).sum() != mesh.n_points:
        raise ValueError("ds does not describe this mesh (node_id / n_points mismatch)")
    skip = {"z", "depth", "node_id"}
    keys = list(names) if names is not None else [
        v for v in ds.data_vars if v not in skip and ds[v].dims == ("layer", "y", "x")]
    ok = nid >= 0
    order = nid[ok]
    for k in keys:
        arr = np.empty(mesh.n_points)
        arr[order] = np.asarray(ds[k].values, dtype=float)[ok]
        set_mesh_property(mesh, k, arr, "point")
    return keys
