"""Single place that defines where a DA run writes its files (same layout for EnKF, PF and ES-MDA).

    <workdir>/<project>/
        <project>_states.nc         final ensemble states (xarray Dataset)
        <project>_states/           per-cycle chunks (crash safety only)
        <project>_params.nc         ES-MDA: parameter ensemble of every iteration + predictions + misfit
        <project>_issues.log        light log of the problems met during the run
        <project>_rejections.json   members rejected, with the reason
        <project>_ERT_checks.json   ERT quality checks
        <project>_run.json          manifest written by ``DA.write_run_summary``
        console_log.txt             full tee of the console (only with ``DAConfig.console_log='full'``)
"""
from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass(frozen=True)
class RunOutputs:
    workdir: str
    project_name: str

    @property
    def root(self) -> str:
        return os.path.join(self.workdir, self.project_name)

    def _file(self, suffix: str) -> str:
        return os.path.join(self.root, f"{self.project_name}{suffix}")

    states_nc = property(lambda self: self._file("_states.nc"))
    states_chunks = property(lambda self: self._file("_states"))
    params_nc = property(lambda self: self._file("_params.nc"))
    issues_log = property(lambda self: self._file("_issues.log"))
    rejections_json = property(lambda self: self._file("_rejections.json"))
    ert_checks_json = property(lambda self: self._file("_ERT_checks.json"))
    run_json = property(lambda self: self._file("_run.json"))
    console_log = property(lambda self: os.path.join(self.root, "console_log.txt"))

    def fresh_chunk_dir(self) -> str:
        """``<project>_states/`` ready for a new run: created, and emptied of the stale chunks of a previous run."""
        import glob

        d = self.states_chunks
        if os.path.isdir(d):
            for f in glob.glob(os.path.join(d, "cycle_*.npz")):
                os.remove(f)
        return d

    def registry(self) -> dict:
        """Name -> path (relative to the project folder) of the outputs that exist now."""
        names = dict(states=self.states_nc, states_chunks=self.states_chunks, params=self.params_nc, issues=self.issues_log,
                     rejections=self.rejections_json, ert_checks=self.ert_checks_json,
                     console_log=self.console_log, run=self.run_json)
        return {k: os.path.relpath(p, self.root) for k, p in names.items() if os.path.exists(p)}
