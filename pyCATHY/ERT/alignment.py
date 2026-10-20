"""ERT (pyGIMLi) <-> CATHY mesh alignment: frame conversion, checks, auto-suggestion.

Frame model: p_cathy = p_ert * signs + shift   (signs = (sx, sy, sz), each +1 or -1;
sz = -1 handles a z-positive-down ERT convention). The SAME transform must be applied to
the electrodes, the pg forward mesh and any electrode block written to data files.

Symptom this guards against: ~50% negative geometric factors k, rhoa far outside the
model range, systematic_misfit every cycle -> garbage ERT forward -> EnKF diverges.
"""
import warnings
import numpy as np


class ERTAlignmentError(RuntimeError):
    pass


def bounds_from_nodes(nodes):
    """[xmin, xmax, ymin, ymax, zmin, zmax] (same order as pyvista .bounds)."""
    n = np.asarray(nodes, float)
    return np.array([n[:, 0].min(), n[:, 0].max(), n[:, 1].min(), n[:, 1].max(),
                     n[:, 2].min(), n[:, 2].max()])


def transform_points(points, signs=(1, 1, 1), shift=(0, 0, 0)):
    return np.asarray(points, float) * np.asarray(signs, float)[None, :] \
        + np.asarray(shift, float)[None, :]


def _overlap_fraction(pb, cb):
    lo, hi = np.maximum(pb[0::2], cb[0::2]), np.minimum(pb[1::2], cb[1::2])
    return float(np.prod(np.clip(hi - lo, 0, None)) / np.prod(cb[1::2] - cb[0::2]))


def _inside_fraction(pb, cb):
    """Fraction of the pg-mesh bounding box that lies inside the CATHY volume (1 = contained)."""
    lo, hi = np.maximum(pb[0::2], cb[0::2]), np.minimum(pb[1::2], cb[1::2])
    vol = float(np.prod(pb[1::2] - pb[0::2]))
    return float(np.prod(np.clip(hi - lo, 0, None)) / vol) if vol > 0 else 0.0


def suggest_frame(electrodes, pg_points, cathy_bounds, z_signs=(1, -1)):
    """Try the z conventions and return the best (signs, shift, overlap).

    shift centres the electrode footprint on the CATHY footprint and puts the highest
    electrodes on the CATHY surface. Best = largest overlap of the transformed pg mesh with
    the CATHY volume (ties keep z_sign=+1). SUGGESTION ONLY: it assumes the array is
    centred on the plot; a known site offset always wins.
    """
    e, cb = np.asarray(electrodes, float), np.asarray(cathy_bounds, float)
    best = None
    for sz in z_signs:
        signs = (1, 1, sz)
        e2 = transform_points(e, signs)
        eb = bounds_from_nodes(e2)
        shift = ((cb[0] + cb[1]) / 2 - (eb[0] + eb[1]) / 2,
                 (cb[2] + cb[3]) / 2 - (eb[2] + eb[3]) / 2,
                 cb[5] - eb[5])
        ov = _inside_fraction(bounds_from_nodes(transform_points(pg_points, signs, shift)), cb)
        if best is None or ov > best[2] + 1e-9:
            best = (signs, tuple(float(s) for s in shift), ov)
    return best


def check_ert_cathy_alignment(electrodes, cathy_bounds, pg_points=None,
                              tol_xy=0.05, tol_z=0.05, tol_node=0.05, min_inside=0.95,
                              on_fail='raise', label='ERT'):
    """Check electrodes (n,3) and the pg mesh nodes (m,3), BOTH ALREADY in the CATHY frame.

    1. electrodes inside the CATHY volume; highest electrodes at the CATHY surface;
    2. array extent / centre consistent with the CATHY footprint;
    3. pg mesh lies INSIDE the CATHY volume (it may be much smaller: the forward mesh is a
       sub-volume around the array, CATHY states are interpolated onto it);
    4. every electrode coincides with a pg mesh node (pyGIMLi electrodes ARE mesh nodes):
       convention-independent test that mesh and electrodes were moved together.
    on_fail: 'raise' | 'warn' | 'ignore'. Returns a dict (ok, problems, ...).
    """
    e = np.asarray(electrodes, float)
    cb = np.asarray(cathy_bounds, float)
    eb = bounds_from_nodes(e)
    problems = []

    inx = (e[:, 0] >= cb[0] - tol_xy) & (e[:, 0] <= cb[1] + tol_xy)
    iny = (e[:, 1] >= cb[2] - tol_xy) & (e[:, 1] <= cb[3] + tol_xy)
    inz = (e[:, 2] >= cb[4] - tol_z) & (e[:, 2] <= cb[5] + tol_z)
    frac_in = float(np.mean(inx & iny & inz))
    if frac_in < 1.0:
        problems.append(f"{(1-frac_in)*100:.0f}% of {len(e)} electrodes outside the CATHY volume "
                        f"(x:{(~inx).sum()}, y:{(~iny).sum()}, z:{(~inz).sum()} out)")
    dz_surf = cb[5] - eb[5]
    if abs(dz_surf) > tol_z:
        problems.append(f"highest electrode z={eb[5]:.3f} vs CATHY surface z={cb[5]:.3f} "
                        f"(offset {dz_surf:+.3f} m: z sign or datum?)")
    c_e = np.array([(eb[0] + eb[1]) / 2, (eb[2] + eb[3]) / 2])
    c_c = np.array([(cb[0] + cb[1]) / 2, (cb[2] + cb[3]) / 2])
    ext_e = np.array([eb[1] - eb[0], eb[3] - eb[2]])
    ext_c = np.array([cb[1] - cb[0], cb[3] - cb[2]])
    if np.any(ext_e > ext_c + 2 * tol_xy):
        problems.append(f"array extent {ext_e.round(2)} m > CATHY footprint {ext_c.round(2)} m (units?)")

    overlap = nearest_max = None
    if pg_points is not None:
        pts = np.asarray(pg_points, float)
        pb = bounds_from_nodes(pts)
        overlap = _overlap_fraction(pb, cb)            # share of the CATHY volume covered (info)
        inside = _inside_fraction(pb, cb)
        if inside < min_inside:
            problems.append(f"only {inside*100:.0f}% of the pg mesh lies inside the CATHY volume "
                            f"(pg {pb.round(2)} vs CATHY {cb.round(2)})")
        # distance of the electrodes to the sides / bottom of the pg mesh (top = free surface):
        # a forward mesh hardly larger than the array biases rhoa through its boundaries
        margin = float(min((e[:, 0] - pb[0]).min(), (pb[1] - e[:, 0]).min(),
                           (e[:, 1] - pb[2]).min(), (pb[3] - e[:, 1]).min(), (e[:, 2] - pb[4]).min()))
        span = float(max(eb[1] - eb[0], eb[3] - eb[2], eb[5] - eb[4]))
        if margin < 0.5 * span and on_fail != 'ignore':
            warnings.warn(f"[{label}] pg forward mesh extends only {margin:.2f} m beyond the "
                          f"electrodes (array span {span:.2f} m): boundary effects can bias the "
                          f"forward apparent resistivities of weak-signal quadrupoles "
                          f"(see nc_model.unreliable_quadrupoles: homogeneous test, rhoa/rho must be 1)")
        from scipy.spatial import cKDTree
        d, _ = cKDTree(pts).query(e)
        nearest_max = float(d.max())
        if nearest_max > tol_node:
            problems.append(f"{(d > tol_node).sum()}/{len(e)} electrodes are not on a pg mesh node "
                            f"(max distance {nearest_max:.3f} m): mesh and electrodes were not "
                            f"moved with the same transform")

    report = dict(ok=not problems, margin_to_mesh_boundary=margin if pg_points is not None else None, problems=problems, frac_electrodes_inside=frac_in,
                  pg_overlap=overlap, max_electrode_node_dist=nearest_max,
                  electrode_bounds=eb, cathy_bounds=cb)
    if problems:
        msg = f"[{label}] mesh/electrode misalignment:\n  - " + "\n  - ".join(problems)
        if on_fail == 'raise':
            raise ERTAlignmentError(msg)
        if on_fail == 'warn':
            warnings.warn(msg)
    return report


def check_k_sign(k, max_frac_neg=0.02, on_fail='raise', label='ERT'):
    """Healthy arrays have geometric factors of one sign; mixed signs = wrong geometry."""
    k = np.asarray(k, float)
    frac = float(np.mean(k < 0))
    frac = min(frac, 1 - frac)
    if frac > max_frac_neg:
        msg = f"[{label}] {frac*100:.0f}% of k have the minority sign (expected ~0)"
        if on_fail == 'raise':
            raise ERTAlignmentError(msg)
        warnings.warn(msg)
    return frac


def transform_tldat_electrodes(src, dst, signs=(1, 1, 1), shift=(0, 0, 0)):
    """Copy a pygimli .tldat/.shm file transforming the electrode block only
    (line0 = n, line1 = '# x y z', then n electrode lines)."""
    with open(src) as f:
        lines = f.read().split('\n')
    nel = int(lines[0].split()[0])
    for i in range(2, 2 + nel):
        tok = lines[i].split()
        xyz = transform_points([[float(t) for t in tok[:3]]], signs, shift)[0]
        lines[i] = ' '.join(f"{v:.6f}" for v in xyz) + (' ' + ' '.join(tok[3:]) if len(tok) > 3 else '')
    with open(dst, 'w') as f:
        f.write('\n'.join(lines))
