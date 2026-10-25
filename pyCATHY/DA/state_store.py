"""Compact storage of the ensemble states of a DA run (replaces the growing long ``df_DA``).

Why: ``df_DA`` has one row per (cycle, member, node) with 7 float64 columns (56 bytes) of which five
repeat the same value for every node of a block, and ``pd.concat`` copied it at every cycle and every
member. For 1936 nodes x 100 members x 362 cycles that is ~4 GB (NC); for a 41.6k-node mesh ~85 GB.

Here the states are kept as float32 arrays, one chunk per cycle:

    psi, sw, analysis : [n_stored_nodes, n_members]      (12 bytes per node and member)
    rejected          : [n_members] int8

Two switches reduce the volume further (``DAConfig.store_nodes`` / ``store_full``):

    store_full  = 'all'       every node at every cycle (default, = the old content)
                = 'analysis'  every node only at the cycles with an analysis (survey cycles);
                              the other cycles keep only ``store_nodes``
    store_nodes = None | node ids kept at the cycles that are not stored in full

``to_frame()`` rebuilds the legacy long table (same columns, compact dtypes, plus ``idnode``), so the
plotting code written for ``df_DA`` keeps working. With ``chunk_dir`` every cycle is also written to
``cycle_<c>.npz`` (O(1) per cycle, a crash keeps the finished cycles); ``StateStore.from_chunks`` reads
them back.
"""
import glob
import os
import re

import numpy as np
import pandas as pd

FIELDS = ("psi", "sw", "analysis")
LEGACY_COLUMNS = ["time", "Ensemble_nb", "psi_bef_update", "sw_bef_update_", "analysis", "OL", "rejected"]


class StateStore:
    def __init__(self, nens, dtype="float32", nodes=None, full="all", chunk_dir=None):
        if full not in ("all", "analysis"):
            raise ValueError(f"store_full must be 'all' or 'analysis', got {full!r}")
        if full == "analysis" and nodes is None:
            raise ValueError("store_full='analysis' needs store_nodes (the nodes kept at the other cycles)")
        self.nens = int(nens)
        self.dtype = np.dtype(dtype)
        self.nodes = None if nodes is None else np.unique(np.asarray(nodes, dtype=np.int64))
        self.full = full
        self.chunk_dir = chunk_dir
        self.attrs = {}                               # str / number metadata written in the NetCDF (da_type, ...)
        self.cycles, self.node_ids, self.is_analysis = [], [], []
        self.data = {f: [] for f in FIELDS}
        self.rejected = []
        if chunk_dir:
            os.makedirs(chunk_dir, exist_ok=True)

    # ------------------------------------------------------------------ writing
    def add(self, cycle, psi, sw, analysis, rejected=None, is_analysis=None):
        """Store one cycle. psi, sw, analysis: [n_nodes, n_members] (all nodes of the mesh)."""
        psi, sw, analysis = (np.asarray(a) for a in (psi, sw, analysis))
        if psi.shape != sw.shape or psi.shape != analysis.shape or psi.shape[1] != self.nens:
            raise ValueError(f"states must be [n_nodes, {self.nens}], got {psi.shape}, {sw.shape}, {analysis.shape}")
        if is_analysis is None:                       # an analysis changed the pressure heads
            is_analysis = bool(np.any(analysis != psi))
        keep_all = self.full == "all" or bool(is_analysis)
        if keep_all or self.nodes is None:
            ids = np.arange(psi.shape[0])
        else:
            ids = self.nodes
            if ids.max() >= psi.shape[0]:
                raise ValueError("store_nodes contains a node id outside the mesh")
        rej = np.zeros(self.nens, np.int8) if rejected is None else np.asarray(
            [int(r) for r in rejected], dtype=np.int8)
        chunk = {f: np.ascontiguousarray(a[ids], dtype=self.dtype) for f, a in zip(FIELDS, (psi, sw, analysis))}
        self.cycles.append(int(cycle))
        self.node_ids.append(ids)
        self.is_analysis.append(bool(is_analysis))
        for f in FIELDS:
            self.data[f].append(chunk[f])
        self.rejected.append(rej)
        if self.chunk_dir:
            np.savez(os.path.join(self.chunk_dir, f"cycle_{int(cycle):05d}.npz"), cycle=int(cycle), node_ids=ids,
                     is_analysis=bool(is_analysis), rejected=rej, nens=self.nens, **chunk)

    # ------------------------------------------------------------------ reading
    def __len__(self):
        return len(self.cycles)

    @property
    def nbytes(self):
        return sum(a.nbytes for f in FIELDS for a in self.data[f]) + sum(i.nbytes for i in self.node_ids)

    def arrays(self, field="sw", nodes=None):
        """(cycles, members, node_ids, array[n_cycle, n_member, n_node]) for the nodes stored at EVERY cycle
        (or the requested ``nodes``, which must be among them). Rejected members are NaN."""
        common = self.node_ids[0]
        for ids in self.node_ids[1:]:
            common = np.intersect1d(common, ids)
        nodes = common if nodes is None else np.unique(np.asarray(nodes, dtype=np.int64))
        missing = np.setdiff1d(nodes, common)
        if missing.size:
            raise ValueError(f"{missing.size} requested nodes are not stored at every cycle (e.g. {missing[:5]})")
        out = np.full((len(self), self.nens, nodes.size), np.nan, dtype=self.dtype)
        for k, (ids, a, rej) in enumerate(zip(self.node_ids, self.data[field], self.rejected)):
            out[k] = a[np.searchsorted(ids, nodes)].T
            out[k][rej > 0] = np.nan
        return np.asarray(self.cycles), np.arange(self.nens), nodes, out

    def to_frame(self):
        """Legacy long table (cycle, member, node), compact dtypes, plus the node id ``idnode``."""
        if not len(self):
            return pd.DataFrame(columns=LEGACY_COLUMNS)
        parts = []
        for c, ids, p, s, a, rej, isa in zip(self.cycles, self.node_ids, self.data["psi"], self.data["sw"],
                                              self.data["analysis"], self.rejected, self.is_analysis):
            n = ids.size
            parts.append(pd.DataFrame({
                "time": np.full(n * self.nens, c, dtype=np.int32),
                "Ensemble_nb": np.repeat(np.arange(self.nens, dtype=np.int16), n),
                "psi_bef_update": p.T.reshape(-1),
                "sw_bef_update_": s.T.reshape(-1),
                "analysis": a.T.reshape(-1),
                "OL": np.zeros(n * self.nens, dtype=np.int8),
                "rejected": np.repeat(rej, n),
                "idnode": np.tile(ids.astype(np.int32), self.nens),
            }))
        return pd.concat(parts, ignore_index=True)

    def to_xarray(self, nodes=None):
        """``xarray.Dataset`` with ``psi``, ``sw``, ``analysis`` [cycle, member, node] (float32), ``rejected`` [cycle,
        member] and ``is_analysis`` [cycle]; coordinates ``cycle`` (the atmbc cycle counter), ``member``, ``node`` (ids).
        Only the nodes stored at EVERY cycle (or ``nodes``) are included, as in ``arrays``; rejected members are NaN."""
        import xarray as xr

        cycles, members, node_ids, _ = None, None, None, None
        data = {}
        for f in FIELDS:
            cycles, members, node_ids, data[f] = self.arrays(f, nodes)
        coords = dict(cycle=cycles.astype(np.int32), member=members.astype(np.int16), node=node_ids.astype(np.int32))
        ds = xr.Dataset(
            {f: (("cycle", "member", "node"), data[f]) for f in FIELDS},
            coords=coords,
        )
        ds["rejected"] = (("cycle", "member"), np.stack(self.rejected).astype(np.int8))
        ds["is_analysis"] = (("cycle",), np.asarray(self.is_analysis, dtype=bool))
        ds.attrs.update(self.attrs)
        ds.attrs.update(nens=self.nens, store_full=self.full,
                        description="pyCATHY ensemble states; cycle c <-> atm_times[c]; psi/sw before the update, analysis after")
        return ds

    def save_netcdf(self, path, nodes=None, complevel=4):
        """Write ``to_xarray()`` as a compressed NetCDF (zlib; needs netCDF4 or h5netcdf, else uncompressed scipy)."""
        return save_dataset(self.to_xarray(nodes), path, complevel)

    @classmethod
    def from_chunks(cls, chunk_dir):
        files = sorted(glob.glob(os.path.join(chunk_dir, "cycle_*.npz")),
                       key=lambda f: int(re.findall(r"cycle_(\d+)", f)[-1]))
        if not files:
            raise FileNotFoundError(f"no cycle_*.npz in {chunk_dir}: this run wrote no states (crashed before the first "
                                    "cycle, or a DA type that does not store them)")
        with np.load(files[0]) as z:
            nens, dtype = int(z["nens"]), z["psi"].dtype
        st = cls(nens, dtype=dtype)
        for f in files:
            with np.load(f) as z:
                st.cycles.append(int(z["cycle"]))
                st.node_ids.append(z["node_ids"])
                st.is_analysis.append(bool(z["is_analysis"]))
                st.rejected.append(z["rejected"])
                for k in FIELDS:
                    st.data[k].append(z[k])
        return st


def save_dataset(ds, path, complevel=4):
    """NetCDF with zlib compression when an engine supports it (netCDF4, h5netcdf), else uncompressed (scipy)."""
    import importlib.util

    for engine in ("netcdf4", "h5netcdf"):
        if importlib.util.find_spec("netCDF4" if engine == "netcdf4" else "h5netcdf"):
            enc = {v: dict(zlib=True, complevel=complevel) for v in ds.data_vars if ds[v].dtype.kind in "fiu"}
            ds.to_netcdf(path, engine=engine, encoding=enc)
            return path
    ds.to_netcdf(path)                                  # scipy: NetCDF3, no compression
    return path


def load_states(path):
    """``xarray.Dataset`` of the states of a run from ``<project>_states.nc`` (a file) or ``<project>_states/`` (the
    per-cycle chunks, a directory). Same layout as ``StateStore.to_xarray``."""
    import xarray as xr

    if os.path.isdir(path):
        return StateStore.from_chunks(path).to_xarray()
    return xr.load_dataset(path)
