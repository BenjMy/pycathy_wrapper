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
