# Testing `version="withIrr"` in the pyCATHY wrapper

> **Living document.** It is both a test plan and a hand-off log. Any AI assistant (or human)
> picking this work up must read sections 0, 3 and 4 first, and update sections 3, 4 and 10 before
> finishing a session.
>
> Last updated: 2026-11-09 (session 11) | Wrapper repo: `BenjMy/pycathy_wrapper` | Fortran src repo: `CATHY-Org/CATHY_src_WithIrrigation`

---

## 0. Ground rules

1. **Only touch the Python wrapper** (`pyCATHY/`). Do **not** modify the Fortran sources (`code/*.f`, `*.H`,
   `Makefile`) nor the source repo's `dataset/`. If a fix seems to require it, **stop and ask the owner** (section 9).
2. **No direct file-to-file comparison** between versions. The wrapper must absorb differences in input-file
   structure (see section 6). Compare *parsed values and physical outputs* with tolerances, never raw text diffs.
3. Prefer expressing version differences as **data in `pyCATHY/version_config.py`** (keyed by version) rather than
   `if self.version == "withIrr"` branches scattered through `cathy_tools.py`. This keeps the later merge small.
4. **`cathy.fnames` is never modified** (owner decision, session 2). The solver reads DEM-related inputs from
   `other-dem/` and the wrapper reads/writes them there for `withIrr` (see `CATHY_DIR_LAYOUT` in `version_config.py`).
5. No hardcoded local paths (`/home/...`) in the wrapper. Source files come from GitHub (section 5, T-01).
6. Do not mark an item done unless it was actually run; record the command and result in section 10.

---

## 1. Goal

Show that `CATHY(..., version="withIrr")` can reproduce, through the wrapper API only:

- **(A)** the wrapper's own smoke test (`test.py` of the source repo),
- **(B)** the existing wrapper tests and documentation examples (as far as they are version-agnostic),
- **(C)** two real cases that so far were run **without** the wrapper:
  - Venezia: `Venezia_2011_asp_W_sceORIG_sceua_2l_v4_KGE_flat_test_colt2Rosetta222748` (`.zip` + `.py`)
  - Salerno: `Salerno_20250131_ROSETTA_dirnew_mars_wobaseflow_Kzexp_hap_t6` (`.zip` + `.py`)

---

## 2. Branch comparison: `withIrrigation` vs `main` (measured 2026-11-06)

| Item | Result |
|---|---|
| Merge base | `7db40ec4` |
| `withIrrigation` ahead of `main` | **1 commit** (`337faba8`, 2025-11-11, "add config version manager for irrigation code") |
| `main` ahead of `withIrrigation` | **61 commits** (latest 2026-10-26: DA refactors, ES-MDA, Archie, ERT, perturbation scenarios, ...) |
| `cathy_tools.py` change on `main` since base | +2212 / -960 lines |
| Trial merge (`git merge-tree`) | **3 conflicts**: `cathy_tools.py` (content), `config_manager.py` (add/add), `version_config.py` (add/add). `importers/cathy_inputs.py` auto-merges. |

What `main` already has that the branch lacks:

- Working `recompileSrc()` (`-O -mcmodel=large -fPIC`, compile/link error handling). The branch version had `bashCommand` commented out.
- `version_config.py` / `config_manager.py` *also* contain `withIrr` entries, plus a `CATHY_SOIL_FORMAT` table (`default`, `withIrr`, `SCF_variable`) and a soil reader adapted to it.
- A GitHub-clone pattern for alternative sources: `_clone_cathy_src()` with a pinned `CATHY_SRC_URL`, and the `SCF_variable` version cloning `CATHY-Org/CATHY_SCF_LAI` into a temp dir. **This is the precedent to follow for `withIrr`.**
- A large pytest suite (60+ files in `pyCATHY/tests/`) and a docs-example test (`test_docs_examples.py`).

What `main` does **not** have (must be ported for `withIrr`):

- A `withIrr` source-fetch block in `CATHY.__init__` (main has none).
- The `prepro_input_dirname = "prepro-dem"` convention (main hardcodes `"prepro"`).
- The `dem_parameters` reader fix for the new `ZRATIO(I) I=1 NSTR` header (section 5, T-03).
- The `other-dem/` -> `input/` relocation and matching `cathy.fnames` rewrite (T-04).

> Python note: `main`'s `test_docs_examples.py` imports `tomllib` (Python >= 3.11). The author's env
> (`micromamba` env `pyCATHY`) showed Python 3.10 in a traceback. Check this before running the suite.

---

## 3. Branch strategy: DECIDED (session 2)

Owner chose **Option A**: new branch **`feat/withIrr-on-main`** created from `origin/main` (tip `8cfc2afc` on 2026-11-06).
The old `withIrrigation` branch is kept untouched as an archive. All work below happens on `feat/withIrr-on-main`.

To keep the final PR small and conflict-free: version differences live in `pyCATHY/version_config.py`
(`CATHY_DIR_LAYOUT`, plus the pre-existing `CATHY_H_*`, `PREPRO_UPDATE_PREPO_INPUTS`, ...); `cathy_tools.py` only asks the
config. Default-version behaviour must stay byte-identical (checked by comparing the same tests on clean `main` and on the branch).

---

## 4. Status tracker (update me)

Legend: `[x]` done and verified | `[~]` partly done / unverified | `[ ]` not started | `[!]` blocked

### 4.1 Preparation

- [x] P-01 Clone wrapper, list branches
- [x] P-02 Compare `withIrrigation` vs `main` (section 2)
- [x] P-03 Branch strategy decided: Option A, branch `feat/withIrr-on-main` (section 3)
- [~] P-04 Python version: owner env is 3.10.14; code written in session 2 is 3.10-compatible and was run on 3.12. `main`'s `test_docs_examples.py` needs `tomllib` (>= 3.11): not run yet, owner "may use a more recent version".
- [x] P-05 Venezia and Salerno **are in the source repo since commit `3234b01` ("add examples", 2026-11-06)**: `examples/<case>/` + `examples/<case>.py` (earlier checks were done on `f172f61`, before that push). See section 11.
- [ ] P-06 Owner to apply `feat_withIrr_on_main.patch`, commit and push (assistant has no GitHub credentials). Patch applies cleanly to `origin/main` (checked with `git apply --check`).

### 4.2 Fixes for `withIrr` on `main`

| ID | Item | Status | Where |
|---|---|---|---|
| T-01 | Fetch `CATHY_src_WithIrrigation` from GitHub (no local paths) | `[x]` | `_fetch_withIrr_src()`; pinned commit `f172f61a66b5ab566170498358e1228549b63be9`; override with `CATHY(..., src_ref="main")` or env `PYCATHY_WITHIRR_REF`. Reuses main's `_clone_cathy_src`. Legacy `prepro`/`input` moves skipped for withIrr. |
| T-02 | `recompileSrc` undefined `bashCommand` | `[x]` | Not needed: `main` already has the fixed version. |
| T-03 | `read_dem_parameters` header `ZRATIO(I) I=1 NSTR` (duplicate `NSTR` overwrote the real value) | `[x]` | `importers/cathy_inputs.py` (keys lowercased, tokens collapsed). Unit-tested. |
| T-04 | `other-dem/` vs `input/` | `[x]` (no `cathy.fnames` edit) | Files stay in `other-dem/`; `CATHY._input_path()` resolves logical names -> real folder/file per version. |
| T-05 | `parm` written in the old layout (no `NCOUT` in the new source) | `[x]` | `parm_fmt` per version in `CATHY_DIR_LAYOUT`; stray `]` removed. |
| T-05b | `dem_parameters` written in the old line grouping | `[x]` | `dem_parameters_fmt` per version. |
| T-06 | Infinite retry in `run_processor` | `[x]` | Not needed: `main` already retries once then raises with stderr. |
| T-07 | Solver stderr hidden | `[x]` | Not needed: `main` already surfaces it for `IPRT1=3`. |
| T-08 | `update_prepo_inputs` hardcoded the terrain parameter list of hap.in (new source adds *Strahler Stream Order*) | `[x]` | Now read from `PREPRO_UPDATE_PREPO_INPUTS[version]`. |
| T-09 | `update_cathyH` rewrote CATHY.H in the old format (header cut at 109 lines, missing `INTEGER` declarations) -> 150 compile errors | `[x]` | `cathyH_mode: "patch"` for withIrr: the source's own `src/CATHY.H` is kept, only derived values are replaced (`ROWMAX, COLMAX, DEMRES, MAXSTR, MAXTRM, MAXZON, MAXVEG, NRMAX, MAXPRT, MAXVP, N1MAX` + user kwargs). `MAXTRM` from the wrapper's formula because the source's 65000 is too small for this domain. |
| T-10 | Hardcoded `"prepro"` (26 places) and `"input"` paths | `[x]` | `CATHY._prepro_path()` / `_input_path()`; `zone -> zone_map`, `root_map -> veg-type_map`, `atmbc -> atmbc-flow`, `ic -> ic-flow`, ... |
| T-11 | `_check_grid3d_vs_dem` looked up key `factor` (withIrr calls it `rescale_factor`) | `[x]` | one-line fallback |
| T-12 | `cathy_DA.py` still has a hardcoded `"prepro"` (and DA copies `input/`) | `[ ]` | DA with `withIrr` is out of scope for now; will break for `withIrr` until done. |
| T-13 | Names not mapped yet (old `mesh`, transport vtk) | `[x]` | `mesh -> xyz-coords` added to `input_file_names` (file locations themselves are read from `cathy.fnames`); transport vtk names `celt<3000+k>.vtk` via `cathy_utils.transport_vtk_filename` (only written with TRAFLAG=1; not exercised by a run). Tests in `test_withIrr_layout.py`. |
| T-14 | New `gen3d.f` writes the mesh body (tetra + nodes, 2-d grid) to Fortran units 1199/1198, which have no file name in `cathy.fnames` -> `fort.1199` / `fort.1198` in the project folder; `output/grid3d` holds only a header (`NNOD N NT`), `grid2d.exp` is empty. `read_grid3d` then returned an empty mesh. | `[x]` (wrapper workaround) | `CATHY._assemble_mesh_files()` rebuilds `output/grid3d` (header + `fort.1199`) and `output/grid2d.exp` after each `IPRT1=3` run; config key `mesh_body_files`. Owner decided (Q8): **no Fortran change**, this workaround is the solution. Relies on gfortran's default `fort.N` naming. |
| T-15 | Flow vtk files are `flow2000.vtk, flow2001.vtk, ...` (`flow<2000+k>`) instead of `100.vtk, 101.vtk, ...`; plotters and Archie hardcoded the old names | `[x]` | `cathy_utils.flow_vtk_filename()` / `flow_vtk_glob()` detect the naming from the folder content (no version plumbing); used by `show_vtk`, `show_vtk_TL`, `petro_Archie`. Default naming unchanged (tested). Transport files are `celt<nnnn>.vtk` (not handled, T-13). |
| T-16 | Without a mesh, `run_processor(IPRT1=2)` ran, failed to read `output/grid3d`, then retried with `IPRT1=3` and **left `IPRT1=3` in `parm`** (next runs silently mesh-only) | `[x]` | For versions with `mesh_body_files`, `run_processor` builds the mesh first when `output/grid3d` is missing/empty, then restores the requested parameters. Default version untouched. |
| T-17 | Unnamed units writing `fort.1133` (recharge.f: TIME, SSTORE1, SSTORE2, REC), `fort.1178` (detout.f: detailed surface output), `fort.1200` (cathy_main.f: TIME, QTRAN, RWU series) | `[x]` | `CATHY._collect_unit_files()` moves them to `output/recharge_storage`, `output/detout_1178`, `output/rwu_series` after each run (config `unnamed_unit_files`, withIrr only). Verified on a real `test.py` run: no `fort.*` left in the project folder. |
| T-24 | **Bounds of `CATHY.H` that depend on the BC inputs were shrunk by the wrapper** (`NQMAX=1` while the Venezia case has 154 Neumann nodes -> the old-source solver stops silently right after reading the inputs, only the initial state is written; its log says `NODI NEUMANN NONATM. = 154 MAX. CONSENTITO=1`). In patch mode `MAXVEG=1` would do the same for a case with several crop classes (`NVEG is too large`). | `[x]` | `_existing_numeric_parameters` + rule "never go below the case's own CATHY.H" for `NPMAX, NPMAX_TRA, NQMAX, NSFMAX, NNSFMX` (rewrite mode) and for the derived bounds (patch mode). Unit-tested. |
| T-25 | Old source writes no `output/grid3d` (0 bytes; nodes in `output/xyz`): the post-run mesh read failed | `[x]` | `run_processor(..., read_grid=False)` (additive option, default unchanged). |
| T-26 | `_check_grid3d_vs_dem` on a project loaded from disk had no DEM and expected the full grid (5475 nodes for Salerno instead of 1692) | `[x]` | reads the project's own `dem` (nodata convention translated). |
| T-27 | `update_parm` crashed with a single print time (`min()` of an empty diff) | `[x]` | guarded. |
| T-21 | `read_parm` creates a bogus key when leftover numbers are glued to a label (Venezia `parm`: `TMAX608400.086400.0`, Fortran ignores them) -> `KeyError: 'TMAX'` | `[x]` | known labels are matched, leftovers dropped; test with the real Venezia file (`tests/data/hand_edited/parm_glued_label`). |
| T-18 | `update_dem` / `update_lakes_map` wrote header `rows`/`cols` swapped (`rows=N, cols=M`); invisible with square DEMs, Fortran end-of-file on the 75x10 DEM. **Also a bug in the default version for non-square DEMs.** | `[x]` | `rows=hapin["M"]`, `cols=hapin["N"]`; unit-tested. |
| T-19 | `update_prepo_inputs(DEM)` passed `DEMRES=min(delta_x, delta_y)` to CATHY.H. `DEMRES` is the grid-generation coarsening factor (1 = one node per DEM cell), not the cell size; with spacing 2 the node arrays were undersized -> `NODE n IS NOT CONNECTED TO ANY ELEMENTS`. **Also affects default with spacing != 1.** | `[x]` | override removed (default stays 1). Dataset DEM now meshes: 836 nodes, 40500 tetra. |
| T-20 | ~~Guard rejecting DEMs with `-9999` for withIrr~~ **withdrawn (session 6): its premise was wrong.** `triangoli.f` meshes a cell iff `dem .ne. 0`; Salerno's `dem` has 0 outside the catchment and its `dtm_13.val` (preprocessor input) -9999. `-9999` is therefore *not* "unsupported", it just is not the processor's mask value. | `[x]` replaced by T-22 | see T-22 |
| T-22 | Nodata convention differs between sources: default `-9999` in `dem`; withIrr `0` in the processor `dem` (-9999 would be meshed as an elevation) | `[x]` | user-facing convention stays -9999; `update_dem` writes `0` for masked cells in `dem` (config `dem_nodata_in_dem_file`) and keeps -9999 in `dtm_13.val`; `read_inputs('dem')` translates back. Verified end to end: DEM with a -9999 band of 10 rows -> mesh of 737 nodes = exactly the 650 active cells. Doc example `plot_3d` (-9999 mask) should now pass for withIrr (not re-run). |
| T-23 | Node areas (`fort.111`) | `[x]` (Q14 default) | new module `pyCATHY/irrigation.py::surface_node_areas(dem, dx, dy)`: lumped areas (two triangles per cell, NW-SE diagonal), node order = `output/xyz`. Reproduces the shipped `fort.111` of **Venezia (154 nodes, rel. err 7e-14) and Salerno (1692 nodes, 4e-16)** from the DEM alone, no Fortran run. Tests `test_irrigation_node_areas.py` (fixtures: both DEMs and `fort.111`). |

### 4.3 Test levels

| Level | Scope | Status |
|---|---|---|
| L0 | `test.py` of the source repo | `[x]` **runs to completion (exit 0)** on 3.12: fetch, parm/hap.in/dem_parameters/CATHY.H updates, preprocessor, compile, mesh run (`IPRT1=3`), final run. *Caveat:* the solver itself reports `CONVERGENCE NOT ACHIEVED IN 10 ITERATIONS` with the settings that `test.py` passes (`DELTAT=5, DTMIN=1e-2, DTMAX=1e2, TRAFLAG=0`; verified they are what ends up in `parm`). With the source's own `parm` the solver iterates normally. This is a property of the test settings, not of the wrapper; owner to decide whether `test.py` should change. `run_processor(verbose=False)` hides the message. |
| L1 | Fortran-free unit tests | `[x]` `pyCATHY/tests/test_withIrr_layout.py`: 22 tests pass (readers, line groupings vs parsed files, path resolution, CATHY.H patching, mesh assembly, vtk naming, raster headers, nodata guard, DEMRES). Fixtures in `pyCATHY/tests/data/withIrr/`. |
| L2 | `main` pytest suite (default version, regression guard) | `[~]` Subset (17 files; re-run in session 6: `main` 170 passed / 3 failed, branch 197 passed (= 170 + 27 new) / the same 3 failed) run on `main` (exported with `git archive`, **not** `git worktree`: the repo checkout is ~3 GB and filled the sandbox disk) **and** on the branch: identical, 170 passed / 3 failed (pre-existing: 2 in `test_soil_write_once.py`, 1 in `test_read_mbeconv.py`). Not collectable on either: `test_mapping.py`, `test_refactor_equivalence.py`, `test_hydro_surfsurbsurface.py`, `test_archie_downstream.py`, `test_archie_scenarios.py` (missing data/imports). Not run: DA/ES-MDA files, `test_docs_examples.py`. |
| L3 | Documentation examples with `version="withIrr"` | `[~]` harness `run_example.py` (examples unchanged, `version` injected, headless pyvista). **Pass: `plot_1`, `plot_3a`, `plot_3c`, `plot_4_pyCATHY_outputs`.** Not wrapper bugs (example/dataset specific): `plot_2` passes `np.ones([N, M])` (cols x rows, only valid for square DEMs); `plot_3_meshing_from_weill`, `plot_3_meshing_from_Camporese`, `plot_3b`, `plot_5a` read `prepro/dem` directly in the script (`prepro-dem/dem` in withIrr); `plot_3d`, `plot_3d_withnonodes` use a `-9999` mask (T-20); `plot_3f` needs missing data `../../data/dtmplot1`; `plot_4b_pyCATHY_outputs` reads `output/grid3d` without any run (default dataset ships one, withIrr does not); `plot_4b_waterTable` indexes output time 3 but only 2 outputs exist because the example time stepping does not converge on the withIrr dataset (see L3 note on `plot_1`). Doc examples were NOT modified (wrapper-only rule): Q10. Earlier `plot_1` note: with the example time stepping the solver does not converge on the withIrr dataset; with the source's own stepping it does. |
| L4 | Venezia | `[~]` **Level A done for an 8-day window** (section 14); Level C (new source, translated inputs) not done |
| L5 | Salerno | `[~]` files available, reconnaissance done (section 11), blocked on questions Q11-Q15 |

---

## 5. Environment and how to run (reproducible)

```bash
# wrapper
git clone https://github.com/BenjMy/pycathy_wrapper.git && cd pycathy_wrapper
git checkout feat/withIrr-on-main   # after applying feat_withIrr_on_main.patch (see P-06)
pip install -e .

# solver build deps (Ubuntu)
sudo apt-get install -y gfortran liblapack-dev libblas-dev

# smoke test, in an empty working dir (it creates dataset_newCATHY/ in the cwd)
mkdir -p /tmp/run && cd /tmp/run
curl -O https://raw.githubusercontent.com/CATHY-Org/CATHY_src_WithIrrigation/main/test.py
python -W ignore test.py
```

Debugging tips learned in session 1:

- To see the real solver error, run it by hand: `cd <workdir>/dataset_newCATHY && ./cathy`.
- `test.py` can loop forever (T-06). Run it with `timeout` or in the background with a log file.
- Sandbox limits for AI sessions: network allowlist includes `github.com`, `raw.githubusercontent.com`, `pypi.org`,
  `archive.ubuntu.com`; single commands are killed after 300 s, so run long jobs with `nohup ... &` and poll the log.
  Avoid `pkill -f <pattern>` with a pattern that also matches your own command line (it killed the shell twice).

---

## 6. Input-file structure variations (wrapper must absorb these)

Rule: the wrapper's **in-memory representation** (dicts / DataFrames from `read_*`) is the common ground. Each
version differs only in how files are **laid out on disk**. Each row below needs a `version_config.py` entry and a test.

| File / folder | `default` (CATHY_src_org) | `withIrr` (CATHY_src_WithIrrigation) | Status |
|---|---|---|---|
| Source tree | `src/`, `prepro/` | `src/` (= repo `code/`), `prepro-dem/src` (= `code/prepro/src`) | T-01 done |
| DEM-related inputs | `input/dem_parameters`, `input/base_map`, ... | listed in `cathy.fnames` under `other-dem/` (`dem_parameters`, `posizione_serb`, `livelli_iniz_s`, `effraininp`, `base_map`, `veg-type_map`, `irrigation_map`) | T-04 done (files stay there) |
| Prepro folder name | `prepro` | `prepro-dem` | T-10 done (`_prepro_path`) |
| `dem_parameters` layout | one value per line, header `zratio(i),i=1,nstr` | `DELTA_X DELTA_Y` on one line; header `ZRATIO(I) I=1 NSTR` | T-03, T-05b done |
| `parm` line 1 | `IPRT1 NCOUT TRAFLAG` (3 values, e.g. dataset_oldCATHY) | `IPRT1 TRAFLAG` (2 values; `datin.f` reads `IPRT1,TRAFLAG`) | T-05 done |
| `parm` rest | see `dataset_oldCATHY/input/parm` | see `dataset/input/parm` | T-05 done: grouping validated against the parsed dataset file (`sum(parm_fmt) == len(read_parm)`) |
| `soil` | per-zone/layer rows, FP params | extra vegetation params (`PCANA PCREF PCWLT ZROOT PZ OMGC`), `main` has `CATHY_SOIL_FORMAT["withIrr"]` | verify on `main` (not exercised by `test.py`) |
| Input names | `atmbc`, `ic`, `mesh`, `root_map`, `transp*` | `atmbc-flow`, `ic-flow`, `atmbc-transport`, `ic-transport`, `dirbc-transport`, `veg-type_map`, `irrigation_map`, `mmm_input`, `xyz-coords`, ... | atmbc/ic/transp family mapped (T-10); `mesh` unverified (T-13) |
| `CATHY.H` | default header params | `CATHY_H_PARAMS["withIrr"]` exists | T-09: source file patched instead (those config entries are currently unused by the writer) |

The source repo ships **`dataset_oldCATHY/`** (old layout: `input/atmbc`, `prepro/`, ...) next to **`dataset/`** (new layout).
Same physical problem, two layouts -> ideal for the equivalence test in L1/L2.

---

## 7. Comparison methodology ("not a direct comparison")

1. **Inputs**: read both decks through the wrapper readers and compare the resulting dict/DataFrame values
   (names normalised, e.g. lowercase keys). A difference in on-disk layout is **not** a failure; a difference in parsed value is.
2. **Outputs**: compare physical quantities, never files byte-for-byte. Use the wrapper's own output readers
   (psi, sw, hgraph/mass balance, grid3d, discharge/outlet if present).
3. **Proposed tolerances** (to confirm with the owner, Q4): fields `rtol=1e-6, atol=1e-8`; mass balance totals `rtol=1e-4`;
   number of time steps / Newton iterations reported but **not** asserted. Compiler flags differ between sources
   (`-O` vs `-O3 -mcmodel=medium` in the source `Makefile`), so bit-identity is not expected.
4. **Any loosening of a tolerance must be logged in section 10 with the reason.**
5. New feature outputs (irrigation) have no `default` counterpart: compare to the original run's results only.

---

## 8. Procedure per level

### L0: `test.py`
1. Run section 5 commands. Expected end state: `run_processor(IPRT1=3)` and the last `run_processor(... DTMAX=1e2 ...)` complete, `output/grid3d` exists.
2. State at end of session 2: completes with exit 0 on Python 3.12 (see 4.3 for the solver convergence caveat).

### L1: Fortran-free unit tests (add to `pyCATHY/tests/`, names suggestion `test_withIrr_*.py`)
- `read_dem_parameters` on both layouts returns identical dicts (T-03).
- parm writer/reader round-trip for `default` and `withIrr` (T-05).
- `dataset_oldCATHY` vs `dataset` inputs parse to equal values where the physics is identical (section 7.1).
- File placement test: after `CATHY(version="withIrr")` the tree has `src/`, `input/`, `prepro-dem/src/`, and `cathy.fnames` has no `other-dem/` (T-01, T-04). Mock the git clone so CI needs no network.

### L2: `main` pytest suite
- Run `pytest pyCATHY/tests -x -q` on the chosen branch with default version. Record pass/fail count in section 10.
- Mark tests that are version-agnostic (pure DA/ERT/stat) vs wrapper-I/O tests; only the latter are worth parametrising on `version`.
- Particularly relevant: `test_mapping.py`, `test_soil_write_once.py`, `test_hydro_surfsurbsurface.py`, `test_read_mbeconv.py`, `test_package_compiles.py`, `test_module_split.py`, `test_refactor_equivalence.py`.

### L3: documentation examples (`doc/content/SSHydro/`)
Order by increasing difficulty and record each in section 10:
1. `plot_1_pyCATHY_weilletal.py` -> 2. `plot_2_pyCATHY_inputs.py` -> 3. `plot_3a_create_soil_withzones.py`, `plot_3b_soil3d_from_weill.py` ->
4. `plot_3c/3d_spatial_atmbc_*` -> 5. `plot_3_meshing_*`, `plot_3f_meshing_from_subCatchment.py` ->
6. `plot_4_pyCATHY_outputs.py`, `plot_4b_*`, `plot_5a_infiltration_withBC.py`. DA examples (`doc/content/DA/`) last.

For each: run once with `version` unset (baseline), once with `version="withIrr"`, compare per section 7. Examples that need data
the repo does not ship are `skipped (reason)`, not failed.

### L4 / L5: Venezia and Salerno (blocked on P-05)
Plan once files are available:
1. Unzip the original run to `reference/<case>/`; keep `.py` untouched as `reference/<case>/original.py`.
2. Inventory the original script: which inputs are written by hand (soil, atmbc, parm, dem_parameters, ...), which parameters differ from the defaults, which outputs were analysed.
3. Refactor to a wrapper script `examples/withIrr/<case>.py` using only wrapper calls (`CATHY(..., version="withIrr")`, `update_parm`, `update_soil`, `update_atmbc`, `update_dem_parameters`, ...). Anything that cannot be expressed through the API becomes a numbered issue in section 4.2, **not** a hand-edited input file.
4. Run, then compare outputs with the reference per section 7 (same output variables the original analysis used, e.g. the KGE/discharge series for Venezia).
5. Record the diff of *parsed inputs* between the wrapper-written deck and the reference deck: expected differences (layout) vs unexpected (values).

---

## 9. Questions

Answered by the owner:
- ~~Q1 branch strategy~~ -> Option A, `feat/withIrr-on-main`.
- ~~Q2 rewrite `cathy.fnames`~~ -> **no**; keep `other-dem/`.

Still open:
- **Q3.** Keep the source pinned to a commit (current: `f172f61`, single commit of the source repo) or follow `main`? (override exists)
- **Q4.** Tolerances for output comparison (section 7.3).
- **Q5.** Re-send Venezia and Salerno `.zip` + `.py` (P-05); which outputs are the "truth"?
- **Q6.** `test.py`: should the last `run_processor` keep `DTMIN=1e-2, DELTAT=5` given that the solver does not converge with them? (Fortran untouched.)
- **Q7.** OK to leave DA (`cathy_DA.py`, T-12) for later?
- ~~Q8 Fortran change for the mesh files~~ -> **no, Fortran is never touched** (owner, session 3); the wrapper workaround of T-14 stays.
- **Q9.** Doc examples: for `withIrr`, should the wrapper warn when the solver reports `CONVERGENCE NOT ACHIEVED` (currently silent with `verbose=False`)?
- **Q10.** Doc examples hardcoding `prepro/dem` / `np.ones([N, M])` (section 4.3 L3): change the doc examples to `simu.read_inputs('dem')` and `np.ones([M, N])` (docs only, not the wrapper, not Fortran)?
- Note (session 4): `TRAFLAG` does **not** change the flow vtk names: `VTK_FLOW` is always called with `2000+kprt`; `TRAFLAG=1` only adds transport files (`VTK_TRAN_CELT`, `3000+kprt`, not handled, T-13). Naming is therefore a property of the source version and is detected from the folder content (T-15).

---

---

## 11. Study cases in the repo (session 5): reconnaissance

Repo `CATHY-Org/CATHY_src_WithIrrigation`, commit `3234b01` ("add examples"), folder `examples/` (1.2 GB: full reference runs). Reference outputs are in each case's `output/` and `vtk/`. **The `.py` scripts are stubs** (only `run_preprocessor`, and for Salerno `update_veg_map` + the two `run_processor` calls); inputs are NOT prepared by the wrapper yet.

### 11.1 Venezia (`Venezia_2011_asp_W_sceORIG_sceua_2l_v4_KGE_flat_test_colt2Rosetta222748`)
- Soil column experiment (SWC probes, Venice 2011), calibrated with SCE-UA in MATLAB (`sceua_main.m`, `functn_DLC.m`, `matlab/*`). Mesh **13 x 10 cells, dx=0.1 m, 36 layers**, 1 zone, TMAX 30 d, 287 prints (15 min). Tiny: fast to run.
- **OLD layout** (`input/`, `prepro/`, no `other-dem/`) and **OLD source**: its `src/` is identical to `dataset_oldCATHY/src` (270 files, +`fbckr.f`). `parm` has `NCOUT` and a hand-edited label (T-21). The script nevertheless says `version='withIrr'`.
- Inputs prepared by R: `R/make_atmbc_fake_ET_Rad_15min_warmup_homog.R` (and variants: uniform, oneline...), `R/make_nansfneubc_v2_R_Rain_ET*.R`, `R/makeIC*.R`, `input/hydparams.csv` (Rosetta), `input/soil_loam*.txt`.
- Reference plot targets (PNG, MATLAB): `matlab/R_SWC_1fig.png` (measurements vs simulated SWC, 6 probes, 24 d), `R_SWC_1fig_atm_sim.png`, `R_SWC_1fig_atm_mes.png`, `R_SWC_atm_sim_mes_rev.png`; simulated SWC data: `R_SWC_long.csv`, `matlab/R_SWC.csv`; probe nodes: `DLC_index_*.asc`, `fort.111` (node areas).
- Status: running the case with its **own old `src/`** through the wrapper (default version): preprocessor OK; `run_processor(IPRT1=3)` exits silently (rc 0, no message) and writes no `output/grid3d` (not diagnosed; old-source behaviour to check). Not run with the new source.

### 11.2 Salerno (`Salerno_20250131_ROSETTA_dirnew_mars_wobaseflow_Kzexp_hap_t6`)
- Catchment with greenhouses, channel, SWMM stormwater, irrigation. Mesh **74 x 72 cells, dx=80 m, 13 layers**, 1536 active cells, **TMAX 52,786,700 s (~611 d)**, `DTMAX=1400 s`. Too long for this sandbox (300 s per command, small disk): plan = full input generation, run a **truncated window** (first N days) and compare to the reference at the same times.
- **NEW layout** (`other-dem/`, `prepro-dem/`) but its `src/` **differs from the repo's `code/`**: 18 files changed (`cathy_main.f` 775 changed lines, `datin_tra.f` 42, `inital.f` 32, `datin.f` 27, `concini_node.f` 20, `volfin.f` 15, `tpnodi.f` 13, `openio.f` 13, ...) plus 1 extra file `bcnxt_DLC.f`; also `cathy_main_Matteo.f` in the case folder. It will not reproduce exactly with `code/` (Q12).
- **`input/atmbc` is excluded by `.gitignore`** (file too large for GitHub?): the processor cannot start (`Cannot open ./input/atmbc`). It must be regenerated or provided (Q13).
- Inputs prepared by R: `R/make_atmbc_moverain_shortevent_impervious_baseflow_spatial_MAY_cathy2_base(_FAST).R` (1297 lines; spatial atmbc from `input/Shiny_Rain_BC*.csv`, `Shiny_ET0_BC*.csv`, crop coefficient maps `OPAC/map_kc_*.csv`, land use, channel routing of impervious runoff), `R/make_nansfneubc_v2_R_ET_Irr_warmup_allnodes.R` (**irrigation**: surface Neumann BC `nansfneubc` from `R/Table_NEUMANN.csv` columns `Irrgreen_m3s`, `ET0green_m3s`, hourly; irrigation set to 0 on DOY 151-211 = "solarization"; volumes converted per node with the node areas `fort.111`; nodes of the greenhouse land use from `other-dem/veg-type_map`), `R/make_dirichlet_v2.R` (+`R/node_dirichlet.txt`), `R/make_Ks.R`.
- The R scripts read solver outputs: `output/xyz` (mesh run, available) and `fort.111` (**node area per node, produced by a commented block in the Fortran**: ReadMe says uncomment, run, stop). The wrapper can compute node areas itself from the mesh/DEM instead (Q14).
- Reference plot targets (JPG, MATLAB; no PNG): `hgraph_compare.jpg`, `hgraph_compare_atm.jpg`, `hgraph_compare_atm_Q_split.jpg`, `recharge_compare_hgatmsf*.jpg`, `matlab/*Q_WL_curve*.jpg`, `WL_correct_baseflow*.jpg`, `Q_correct_baseflow*.jpg`, `hgirrigation_irrnotused.jpg`; input maps: `Kc_1..21.jpg`, `Land_cover.jpg`, `Zone.jpg`, `Channel.jpg`, `DTM.jpg`, `Inflow_area.jpg`, `R_ET.jpg`, `R_Rain_ET.jpg`.
- **CORRECTED (session 8)**: irrigation IS an input raster. Both irrigation sources (B and C) read `other-dem/irrigation_map` (unit IINDEM12), see 11.7. The Neumann route of `make_nansfneubc_v2_R_ET_Irr_warmup_allnodes.R` is an alternative approach (`hgirrigation_irrnotused.jpg`), not what Salerno's reference run uses.

### 11.3 Proposed plan (to be confirmed with the owner)
1. **Level A (execution)**: run each case's *original, complete inputs* through the wrapper and compare outputs with the reference (Venezia SWC series at probes; Salerno hydrograph/recharge over a truncated window). Needs: Venezia old-source mesh issue solved; Salerno `atmbc` available.
2. **Level B (input generation)**: regenerate each input with the wrapper (new module, e.g. `pyCATHY/irrigation.py` + atmbc/nansfneubc/dirichlet helpers; Python port of the R logic) and compare **parsed values** with the original inputs (tolerances), then outputs.
3. **Level C (cross-version, Venezia)**: build a new-layout Venezia project from the old inputs using the layout tables of `version_config.py` and run it with the new source `code/`; compare with the old-source reference run.
4. Validation plots: re-plot in Python the quantities of the MATLAB figures above and compare visually/numerically.

### 11.4 Open questions (owner)
- ~~Q11~~ **answered (session 7): consider both** - see 11.6.
- ~~Q11-orig~~ **Q11.** Venezia was run with the OLD source (`src` = `dataset_oldCATHY/src`). Which reproduction do you want: (a) its own old `src/` via the default version, (b) the new irrigation source `code/` (cross-version regression, Level C), or (c) both? I suggest (c).
- ~~Q12~~ **answered (session 7): consider both** - see 11.6.
- ~~Q12-orig~~ **Q12.** Salerno's `src/` differs from `code/` (18 files, `cathy_main.f` +775 lines, `bcnxt_DLC.f`). Which is the reference for reproduction: Salerno's own `src/`, or `code/` (are those differences already merged / obsolete)? Using the case's own `src/` needs no Fortran edit.
- ~~Q13~~ **answered (session 7): owner will share `atmbc` later** (heavy file on their PC).
- ~~Q13-orig~~ **Q13.** Salerno `input/atmbc` is gitignored. Can you provide it (upload / Zenodo / LFS) to validate a regenerated one value by value? Otherwise I regenerate it from the CSV tables and cannot check it against the original.
- ~~Q14~~ answered (session 10): computing node areas in Python is OK. Evidence on the sources: the `write(111,...)` node-area block is **commented** in the old source (A, `cathy_main.f` / `etran.f` comments) and **absent** in B and C (no `areanod`/`write(111` at all), so `fort.111` cannot be produced from B or C without editing the Fortran. The shipped `fort.111` files come from an edited build; `surface_node_areas` reproduces them to 1e-13.
- ~~Q15~~ answered: **yes, validate Salerno on a truncated window** (first N days: a wet event + the irrigation period; N to be chosen when the `atmbc` is available). Acceptance figures: `hgraph_compare` (discharge) and `recharge_compare_*` first; `Kc`/land-cover/zone maps are input checks.
- ~~Q16~~ answered: **yes, in `pyCATHY/` (`pyCATHY/irrigation.py`) and callable as wrapper methods**. Done for the irrigation map: `simu.update_irrigation_map(times=..., rates=..., irr_sw=..., fields=[...])` (session 10). `update_atmbc` already exists; a Neumann-BC method is not needed since the reference runs use the irrigation module.

### 11.5 Progress (session 6)
- **Done, validated against the repo's reference runs**: surface node areas from the DEM (T-23). This removes the dependency on the commented Fortran block that writes `fort.111` (Q14 resolved by default: computed in Python, matches exactly).
- **Done**: per-source nodata convention (T-22) - needed to reproduce Salerno's irregular catchment (mask = `dem == 0` for the processor).
- **Found**: `dataset_oldCATHY/src` is incomplete (missing `fbckr.f`, the linker fails: `undefined reference to fbckr_`). Venezia's `src/` has that file (it is the "+1 file"). Also `main` crashes on an existing project without `output/` (legacy `rmtree(tmp_src)` block).
- ~~Irrigation script reading notes~~ **(partly obsolete, see 11.7: the reference run uses `irrigation_map`, not this Neumann script)** (`make_nansfneubc_v2_R_ET_Irr_warmup_allnodes.R`): irrigation is applied on the nodes of **layer `Nlayer <- 3`** of `output/xyz` (comment: "1 for atmbc, 2... subsurface; plastic says that Nlayer=2 has Kz=0"), i.e. a *subsurface* Neumann BC, not the surface; node land use is taken from the 1-4 surrounding cells (edge/corner cases); `Irrgreen_m3s/3`; irrigation zeroed for DOY 151-211; consecutive zero-irrigation steps are removed before writing. A faithful port + validation against the shipped `input/nansfneubc` (present in the repo, only `atmbc` is gitignored) is the next step (needs `R/Table_NEUMANN.csv`).
- Still open: Q11-Q13, Q15, Q16 (defaults assumed in the meantime: both Venezia variants; truncated Salerno window; module in `pyCATHY/`).

### 11.6 Which source is which (session 7, measured on commit `3234b01`)
There are **three** Fortran sources involved, none identical to another:

| Source | Where | Irrigation | MMM scheme | Input names | Folders |
|---|---|---|---|---|---|
| **A. old** | Venezia `src/` (= `dataset_oldCATHY/src`, 270 files, +`fbckr.f`) | no (0 mentions of IRR) | no | `atmbc`, `ic`, `transp*`, `root_map` | `input/`, `prepro/` |
| **B. irrigation, Sep 2024** | Salerno `src/` (279 files) | yes (289 IRR mentions; header "Sep/24: Added an option for input of irrigation data") | no (11 mentions) | `atmbc`, `ic`, `transp*` (still old names), `veg-type_map` | `other-dem/`, `prepro-dem/` |
| **C. irrigation + MMM, Sep 2025** | repo `code/` (283 files; this is what `version="withIrr"` fetches; unchanged between `f172f61` and `3234b01`) | yes | yes (header "Sep/25: Michaelis-Menten-Monod (MMM) scheme", 279 mentions, `MMM_CBLOCK.H`, `datin_mmm.f`, `mmm_*.f`) | `atmbc-flow`, `ic-flow`, `transport-parm`, ... | `other-dem/`, `prepro-dem/` |

- Venezia src vs `code/`: 70 files differ (3412 changed lines; `cathy_main.f` 1334, `openio.f` 252, `datin.f` 198, `inital.f` 194, `atmone.f` 127, ...); `code/` adds `irr_*.f`, `interp_irr.f`, `moistab_irr.f`, `vtk_flow.f`, `vtk_tran_*.f`, `mmm_*.f`, `celldetq.f`, `pondout.f`, ...; Venezia has 6 files `code/` dropped/renamed (`Conf_Fac.f`, `init0.f`, `vtkris3d*.f`).
- Salerno src vs `code/`: 18 files differ (1000 changed lines, mostly the Sep/25 MMM additions plus refactors); `code/` adds `datin_mmm.f`, `lel_kel.f`, `mmm_*.f`; Salerno has the author's `bcnxt_DLC.f` (a call to it is commented out in its `cathy_main.f`) and keeps some Sep/24 code. So **B is an earlier snapshot of C plus the DLC variant**; C is the newer code. Whether B and C give the same flow results (MMM concerns solutes; Salerno runs with `TRAFLAG=0`) is an *empirical* question for Level A vs Level C.
- **Consequence**: `cathy.fnames` of the project is the source of truth for file names/folders. The wrapper now reads it (`_input_path` / `_prepro_path`, session 7) and falls back on `CATHY_DIR_LAYOUT` only when the project has no `cathy.fnames`. Tested on the three real `cathy.fnames` (fixtures in `tests/data/fnames/`). The `input_file_names` table of `withIrr` therefore only describes source C.
- **Decision (owner asked to consider both, session 7)**: for each case run (1) the case's **own** `src/` (reproduction of the authors' results) and (2) `code/` (C) (does the current code reproduce them). Venezia (2) needs the old->new layout translation (Level C). Using the case's own `src/` needs a wrapper option to take `src/` + `cathy.fnames` from a local folder instead of fetching (to implement: e.g. `CATHY(..., src_dir=...)`; today the scripts work by pointing `dirName` at the case folder).
- Salerno `atmbc`: the owner will share it later (heavy, on their PC). Until then Salerno runs cannot start; Level A for Salerno is blocked; the input generators (irrigation, atmbc) can be developed and validated against the other shipped inputs (`nansfneubc`, ...).

### 11.7 Irrigation = `other-dem/irrigation_map` raster (session 8, from `irr_read_one.f`/`irr_read_nxt.f` and Salerno's real file)
- **Decision (owner): use the most updated source = `code/` (C, Sep/25)**; the case's own sources (A, B) remain for reproduction checks (11.6).
- Both B and C list `'other-dem/irrigation_map' unit IINDEM12` and `'output/hgirrigation' unit IOUT33` in `cathy.fnames`. The `code/` dataset ships it **disabled** (`IRR_FLAG=0`, no rasters).
- **Format**: `IRR_FLAG` (1 = module on) / `IRR_SW` (apply irrigation at a surface node only if its water saturation < IRR_SW) / `IRR_HS` (0 = raster per time, 1 = one homogeneous rate per time); then for each time of change: `TIMEIN` followed by NROW lines of NCOL rates [L/T = m/s], matrix convention (row 0 = north), piecewise constant in time; after the last time the rate stays constant. Only for DEM-input mode with the subsurface module on; irrigation output `output/hgirrigation`.
- **Salerno's file** (15.8 MB): `IRR_FLAG=1, IRR_SW=0.35, IRR_HS=0`, **1224 time blocks** (irregular times: 0, 32400, 36000, 118800, ...; 611 blocks irrigated, the others zero rasters), max 384 irrigated cells per block (greenhouses), rates ~1e-7 m/s (0.3-0.4 mm/h). `irrigation_map_0` = same with `IRR_FLAG=0`, `IRR_SW=0.65`.
- Produced by the **same R script as the atmbc** (`make_atmbc_..._cathy2_base(_FAST).R`, lines ~990-1110: `matirr` map x `Irrgreen_m3s[ii]`, with `Irrgreen` from `Table_NEUMANN.csv`/Irr table in mm per step -> m/s, `jfactor_warmup`, solarization logic, zero blocks compressed).
- **Done**: `pyCATHY/irrigation.py::read_irrigation_map / write_irrigation_map` (round trip exact; fixtures = excerpts of the real file; `tests/test_irrigation_map.py`).
- **To do**: (1) a generator `irrigation rates table + greenhouse mask raster -> irrigation_map` reproducing Salerno's file (validation target: the shipped 1224 blocks) - needs `R/Table_NEUMANN.csv`/the irrigation table and the `matirr` mask (probably `other-dem/veg-type_map` land-use class) from the R script lines ~990-1110; (2) a wrapper method (e.g. `simu.update_irrigation_map(...)`) that writes it to `other-dem/` through `_input_path`; (3) check against `output/hgirrigation` of the reference run.

## 12. Synthetic irrigation district (session 9): `synthetic_irrigation_district/`

Files (outputs folder): `run_synthetic_district.py` (single multi-year run), `plot_synthetic_district.py` (inputs, soil moisture, irrigation diagnostics, maps), `district_lib.py` + `run_scenarios.py` + `plot_scenarios.py` (scenario matrix, module vs classic). Everything is prepared and run **only through the wrapper** with the most updated source (`version="withIrr"`, `code/`).

- **Domain**: 16 x 16 cells of 25 m (16 ha), quasi flat (0.01 % slope: a strictly flat DEM makes the wrapper put a 0 in cell [0,0], which is *inactive* for this source), 12 layers (0.06 m ... 6 m), homogeneous default soil, `IC = -1.5/-2 m`, no lateral/bottom BC (sealed).
- **Forcing**: daily synthetic Mediterranean series (wet winters, dry summers, ~700 mm rain/y, ET0 ~1100 mm/y; potential ET = 0.7 ET0) -> homogeneous `atmbc` (HSPATM=1, IETO=1).
- **Irrigation**: `other-dem/irrigation_map` written with `write_irrigation_map`, schedules built with `build_irrigation_map` (new): 4 blocks of 6 x 6 cells: A drip (3 mm/d), B sprinkler (20 mm / 3 d), C furrow (48 mm / 10 d), D rainfed; season 15 Apr - 15 Sep; `IRR_SW = 0.8`.
- **Crops**: one vegetation class per field (`veg-type_map`), Feddes parameters per class (`FP_map`: PCANA, PCREF, PCWLT, **ZROOT**, PZ, OMGC); root depth is a per-class input and is used by the root water uptake (`etran.f`).
- **Classic comparison**: no irrigation module; irrigation is added to the rain in a spatially variable `atmbc` (HSPATM=0, per surface node) using the new `cell_raster_to_node_rates` (volume conserving cell -> node transfer).
- **Solver behaviours measured on this source** (important for users):
  1. Solver time steps do **not** land on `atmbc` times (checked on all 24 records of a run: none hit except TMAX); each step uses the record valid at the step end (`ATMTIM(2) < TIME <= ATMTIM(3)`; same for the irrigation records, `irr_read_nxt.f`). Piecewise-constant forcing/irrigation is therefore smeared by up to one step (per event -2 % ... +25 % with 1 h steps; +9 % over a 60 d window). Owner view (session 9): acceptable, what matters is that runs reach the end; kept as a documented property, same for module and classic.
  2. `hgirrigation` / `hgatmsf` repeat a step line when a step is rejected and redone: cumulative columns over-count (22.7 vs 15.7 mm); the scenario post-processing keeps the last line per step.
  3. `hgflag` logs every air-dry surface node at every step: GBs in dry summers (359 MB in 100 dry days). Keep it sparse during the run (truncate, Linux keeps the file sparse); `district_lib` does it with a thread.
  4. `IPRT=4` writes velocity fields (0.9 MB per print for 3757 nodes); `IPRT=2` does not write the 3-d `sw`.
  5. Background jobs die at session boundaries: runs must be short and **resumable** (`run_scenarios.py` skips finished cases).
- **Sandbox**: 1 CPU core; a 60-day window runs in ~2.5 min incl. two compilations.

### 12.1 Scenario results (session 9): module vs classic, 60-day summer window (days 170-230)

13 cases = 6 scenarios x {irrigation module, classic (irrigation added to rain, spatial atmbc)} + threshold case (module only). Default (highly permeable, Ks 1.9e-4 m/s) soil, 16 x 16 cells, `IC = -1.5 m`. Files: `synthetic_irrigation_district/plots_scenarios/` (`fig5_rootzone_<scenario>.png`, `fig6_scenarios_summary.png`, `scenario_table.md`, `scenario_metrics.json`), raw small results in `scenario_results/*.npz`.

Scenarios: S0 base (blocks, root 0.6 m, morning/evening events), S1 root 0.3 m, S2 root 1.0 m, S3 irrigation at night (22:00), S4 strips layout, S5 checkerboard layout, S6 threshold `IRR_SW = 0.5` (module only; classic has no threshold).

| scenario | mode | rz Sw irr | stress d | deep gain mm | irrigation mm: module measured / classic as given (map nominal) | not used mm | OVL mm | ACT mm |
|---|---|---|---|---|---|---|---|---|
| S0_base | module | 0.421 | 38.1 | 40.8 | 132.7 (122.1) | 0.9 | 0.00 | 166.1 |
| S0_base | classic | 0.414 | 40.0 | 32.7 | 122.1 (122.1) | 0.0 | 0.00 | 157.7 |
| S1_root0.3 | module | 0.420 | 37.4 | 55.2 | 131.4 (122.1) | 0.8 | 0.00 | 163.0 |
| S1_root0.3 | classic | 0.414 | 38.1 | 48.3 | 122.1 (122.1) | 0.0 | 0.00 | 156.1 |
| S2_root1.0 | module | 0.417 | 43.9 | 28.2 | 133.2 (122.1) | 1.0 | 0.00 | 168.6 |
| S2_root1.0 | classic | 0.413 | 45.8 | 19.6 | 122.1 (122.1) | 0.0 | 0.00 | 161.2 |
| S3_night | module | 0.422 | 36.1 | 39.6 | 132.9 (122.1) | 0.9 | 0.00 | 165.4 |
| S3_night | classic | 0.413 | 38.7 | 31.7 | 122.1 (122.1) | 0.0 | 0.00 | 156.0 |
| S4_strips | module | 0.444 | 33.9 | 63.0 | 174.4 (158.6) | 0.0 | 0.00 | 200.4 |
| S4_strips | classic | 0.438 | 33.9 | 50.4 | 158.6 (158.6) | 0.0 | 0.00 | 189.5 |
| S5_checker | module | 0.443 | 38.7 | 63.9 | 126.6 (112.5) | 0.0 | 0.00 | 167.3 |
| S5_checker | classic | 0.430 | 56.1 | 43.1 | 112.5 (112.5) | 0.0 | 0.00 | 150.0 |
| S6_thresh0.5 | module | 0.329 | 60.0 | -6.5 | 48.9 (122.1) | 82.4 | 0.00 | 80.2 |

Reading (what the numbers support, no more):
- Module and classic differ little here: root-zone saturation within 0.003-0.013, stress days within 2 d except the checkerboard (39 vs 56 d), and no overland flow in any case. With this very permeable soil neither approach ponds, and `IRR_SW = 0.8` almost never blocks irrigation (0.0-1.0 mm of ~130 mm not used).
- The module delivers 9-12 % more irrigation than the map nominally contains (132.7 vs 122.1 mm; 174.4 vs 158.6; 126.6 vs 112.5) because steps are not aligned with events (section 12, solver behaviours). The classic irrigation is reported **as given in the atmbc, not measured**, so the classic columns may also be smeared by the same amount; the "water gained below the roots" is therefore higher for the module partly because it receives more water.
- Root depth (S1/S2) changes the water below the roots (55 / 41 / 28 mm for 0.3 / 0.6 / 1.0 m in the module) and the stress days (37 / 38 / 44), as expected for deeper roots drawing the profile down; same ranking in classic.
- Timing (S3 night vs S0): small effect with this soil (36 vs 38 stress days).
- S6: with `IRR_SW = 0.5` the module withholds 82 of 131 mm, applies 49 mm and the irrigated fields are stressed the whole window (60 d, root-zone S_w 0.33): the threshold is a strong control, `IRR_SW` must be chosen above the saturation reached after irrigation.
- Limits: 60 days not multi-year (the 3-year run did not survive a session pause); homogeneous atmbc (same ET on all fields); one soil; classic not measured; stress threshold S_w < 0.45 is a fixed indicator, not the Feddes function. To see ponding/runoff differences a less permeable soil would be needed (Q21).
- **Q21.** Rerun the scenarios with a less permeable soil (loam, Ks ~ 3e-6 m/s) so that ponding, runoff and the threshold effect show up?

---

## 13. Toolchain without system gfortran / LAPACK

Moved out of this plan (owner decision, session 10): see `NO_SYSTEM_DEPS.md`, to be tried in a separate branch (`feat/no-system-deps`). **Nothing of it is in `feat/withIrr-on-main`** once the increment `feat_withIrr_increment.patch` is applied (the first push, `bd7c0dcf`, still contained the unused `pyCATHY/fortran/blas_lite.f`; the increment removes it).

## 14. Reproducing the study cases through the wrapper (session 11): `reproduce_cases/`

Files: `case_tools.py` (repo clone, trimmed case staging, nodal comparison), `reproduce_venezia.py`, `reproduce_salerno.py`, `results/` (Venezia plot and numbers). Both scripts use the case's **own** `src/` (as owner decided: consider both sources; the new `code/` run is the next step) and original inputs, through the wrapper; the shipped reference outputs of the repository are the target.

**Venezia, level A** (own old source, default version, layout `input/` + `prepro/`, 8-day window starting from the case's first print at day 6.08): the wrapper-driven run reproduces the shipped `output/sw`:

| quantity | value |
|---|---|
| print times compared | 24 (days 0 and 6.09-7.92) |
| S_w range (reference) | 0.397 - 0.866 |
| max / mean abs difference | 7.5e-4 / 1.6e-5 |
| rmse | 5.6e-5 |

Needed wrapper fixes found on the way: T-21 (glued label in `parm`), T-24 (BC bounds of `CATHY.H`), T-25 (`read_grid=False`), plus passing the original print times truncated to the window (`TMAX` alone is reset to the last print time) and `MAXVEG`. Run time of 8 days: ~3 min (1 core); the full 30 days (~11 min) was not run (sessions are short).

**Salerno** (own irrigation source of Sep 2024, new layout), stages `stage / check / run / compare`:
- `check` (done): the original `parm` has the 50 entries of the wrapper table; the shipped `irrigation_map` has `IRR_FLAG=1, IRR_SW=0.35, IRR_HS=0`, 1224 blocks, nominal 22565 m3 in the first 3 days.
- `run` (mechanics tested with a zero-flux `--smoke` atmbc, 6 h, 130 steps): the wrapper compiled the case's own source, built exactly the shipped mesh (1692 surface nodes, 23688 nodes), ran it, collected the `fort.*` units into `output/`, no stray files. **Results of the smoke run are not comparable** (no real forcing).
- `compare` (written, not exercised): `hgatmsf` cumulative fluxes, `hgraph` outlet series and `hgirrigation` (steps de-duplicated) vs the reference over the window; needs the real `atmbc`.
- Blocker: `input/atmbc` is gitignored in the repository. Owner will share it; then: `python reproduce_salerno.py --days N --atmbc /path/to/atmbc`.

**Salerno `irrigation_map` regenerated exactly** (`reproduce_cases/validate_salerno_irrigation_map.py`, no `atmbc` needed): hourly `Irr_m3s` of `R/Table_NEUMANN.csv`, land use `other-dem/veg-type_map` (classes 3 and 4 irrigated), monthly Kc maps `OPAC/map_kc_YYMM.csv`, record selection of the R script (keep an irrigated hour and the first zero hour after it). Result vs the shipped 15.8 MB file: 1224 / 1224 blocks, identical times, identical irrigated cells, max rate difference 0.0 (4e-14 after writing and reading back through `simu.update_irrigation_map`). Code in `pyCATHY/irrigation.py` (`select_irrigation_records`, `irrigation_map_from_series`), tested; `reproduce_salerno.py --regen-irrigation` uses it.

**Venezia level C (WIP, not working yet)**: `reproduce_cases/reproduce_venezia_levelC.py` translates the original inputs to the new source (DEM/layers/zone/veg/hap.in/dem_parameters/parm through the wrapper API; soil, ic, atmbc, nansf*, sfbc, retctab, nudging copied by logical name through `cathy.fnames`) and runs `code/`. First run: the project is built and the solver runs, but `output/sw` only holds the initial state (the same symptom the old source showed before T-24 was fixed), so the comparison step fails. Not diagnosed yet: check `output/risul` of the new run for a bound message (`MAX. CONSENTITO`) and the CATHY.H bounds of the new project (patch mode) first.

Not done: Level C for Venezia (same inputs translated to the new layout and run with `code/`), regeneration of Salerno's `atmbc` (the `Shiny_*` CSVs, Kc maps, land use, channel and impervious-runoff logic of the 1297-line R script are in the repo; not ported), acceptance figures (`hgraph_compare`, `recharge_compare_*`) once the atmbc is available, a full 30-day Venezia run.

---

## 10. Session log (append only)

### Session 1: 2026-11-06
- Cloned wrapper; compared branches (section 2).
- Ran `test.py` from the source repo against the `withIrrigation` branch; peeled off blockers T-01 to T-04 in order, each verified by re-running; stopped at T-05.
- Delivered `withIrrigation_fixes.patch` (T-01 to T-04). Not pushed by the assistant (no credentials). **Unknown whether the owner applied it.**
- Findings T-05 to T-07 recorded, no code written for them.
- Venezia/Salerno files not received -> L4/L5 blocked.
- Environment used: Ubuntu 24, Python 3.12, gfortran + liblapack/libblas installed via apt; pyvista 0.49.2; GitPython 3.2.1.
- Not run: `main` pytest suite, doc examples.

### Session 2: 2026-11-06
- Created `feat/withIrr-on-main` from `origin/main`; ported T-01, T-03, T-04 (without `cathy.fnames` edit), solved T-05, found and solved T-08, T-09, T-10, T-11 (section 4.2). T-02, T-06, T-07 were already fixed on `main`.
- Debug path (useful pattern): run `test.py`, read the first failure, run `./cathy` by hand in the project folder, compare with the source dataset, fix in the wrapper only.
- Verified: `test.py` exits 0; `cathy.fnames` of the project is identical to the source's; `test_withIrr_layout.py` 12 passed; regression subset identical on `main` and branch.
- Deliverable: `feat_withIrr_on_main.patch` (7 files, applies cleanly to `origin/main`). Not committed/pushed by the assistant.
- Not done: doc examples (L3), `main` DA/ERT tests, Venezia/Salerno (files missing), T-12, T-13.
- Env: Ubuntu 24, Python 3.12, gfortran + LAPACK/BLAS; owner env Python 3.10.14.

### Session 3: 2026-11-06
- Venezia/Salerno: searched all branches of the wrapper and the source repo (single commit `f172f61`), web search: **not on GitHub**. Still need upload or push (P-05).
- L3 started: example 1 through the harness. Found and fixed T-14, T-15, T-16 (wrapper only); observed T-17.
- Verified: L0 `test.py` exit 0 again with complete `output/grid3d` (1.8 MB); `cathy.fnames` identical to the source; 16 unit tests pass; regression subset identical on `main` and branch (170 passed / 3 pre-existing failures).
- Deliverable: updated `feat_withIrr_on_main.patch` (10 files, applies to `origin/main`).
- Sandbox lesson: do not `git worktree add` this repo (3 GB); use `git archive origin/main pyCATHY | tar -x -C <dir>`.
- Not done: L3 examples 2+, DA, Venezia, Salerno.

### Session 4: 2026-11-07
- Found and fixed T-18, T-19, T-20 (T-18/T-19 also affect the default version for non-square DEMs / spacing != 1). Read all L3 logs, classified them (L3 row).
- Verified TRAFLAG vs vtk names in the Fortran: no dependency on TRAFLAG for flow files (see note in section 9).
- Regression: `main` 170 passed / 3 failed; branch 190 passed / same 3 failed. `test_withIrr_layout.py`: 20 passed.
- Deliverable: updated `feat_withIrr_on_main.patch`.
- Not done: Venezia/Salerno (files missing), DA (T-12), T-13, remaining L3 follow-ups (Q10).

### Session 5: 2026-11-07
- Venezia/Salerno **found** in the source repo (commit `3234b01`), reconnaissance in section 11 (layouts, source versions, R/MATLAB inputs, plots, blockers).
- T-21 `read_parm` glued label fixed; regression subset: branch 191 passed / same 3 pre-existing failures; patch regenerated (11 files).
- Salerno script run on a trimmed copy: preprocessor and `update_veg_map` OK, solver stops at `Cannot open ./input/atmbc` (gitignored). Venezia with old src: preprocessor OK, mesh run silent exit.
- Not done: reproduction of either case (waiting on Q11-Q16), irrigation/atmbc module.

### Session 6: 2026-11-08
- Withdrew the nodata guard (T-20) after finding Salerno's own files contradict it; implemented the correct per-source convention (T-22) and verified masked meshing end to end.
- Added `pyCATHY/irrigation.py` (`surface_node_areas`), validated on both cases (T-23).
- Regression: `main` 170 passed / 3 failed; branch 197 passed / same 3 failed. Patch regenerated (applies to `origin/main`).
- Not done: irrigation nansfneubc port, atmbc port, Level A/C runs, plots.

### Session 7: 2026-11-08
- Answered "own src vs `code/`" with measurements (11.6): three different sources (old / irrigation Sep24 / irrigation+MMM Sep25); owner: consider both; `atmbc` of Salerno to be shared later.
- Found Salerno's `cathy.fnames` keeps the old input names (`atmbc`, `ic`), unlike `code/` (`atmbc-flow`): file locations are now read from the project's `cathy.fnames` (fallback: version table). 30 layout/area tests pass.
- Not done: `src_dir` option, irrigation `nansfneubc` port, atmbc port, Level A/C runs.

### Session 8: 2026-11-08
- Owner: use the most updated source (`code/`). "Is it an irrigation map raster?": **yes**, verified in the Fortran and in Salerno's real file; corrected two wrong notes in 11.2/11.5 (Neumann route is an unused alternative).
- Added `read_irrigation_map` / `write_irrigation_map` + 5 tests (real-file excerpts). Regression subset: see patch notes; patch regenerated.
- Not done: generator reproducing the shipped irrigation_map, `update_irrigation_map` method, atmbc port, Level A/C runs.

### Session 9: 2026-11-09
- Synthetic irrigation district (section 12) built and run only through the wrapper; 13 scenario cases (module vs classic) run, results/plots in `synthetic_irrigation_district/`.
- New in `pyCATHY/irrigation.py`: `surface_node_areas` (validated vs `fort.111` of both cases), `read_/write_irrigation_map` (validated on Salerno excerpts), `build_irrigation_map` (field schedules), `cell_raster_to_node_rates` (volume conserving). New `pyCATHY/fortran/blas_lite.f` (not yet used by `recompileSrc`).
- Toolchain study (section 13): LAPACK-free and system-gfortran-free builds are byte-identical on a 6-day case; oversized `CATHY.H` is byte-identical too.
- Fixes found along the way: nodata convention per source (-9999 vs 0 in the processor dem), cathy.fnames-driven file resolution, glued labels in `read_parm`.
- Regression: `main` 170 passed / 3 failed; branch 210 passed / the same 3 failed. Patch: 24 files, applies to `origin/main`.
- Not done: S5/S6 follow-ups, plot_synthetic_district.py not run on a real multi-year output, `blas_lite`/managed gfortran integration in `recompileSrc`, Venezia/Salerno reproduction (waiting for the Salerno atmbc), irrigation Neumann/atmbc port.

### Session 10: 2026-11-09
- **Upstream checked** (owner warned of changes): source repo unchanged (`3234b01`); wrapper `main` unchanged (`8cfc2afc`); **new remote branch `feat/withIrr-on-main` = owner's push of the session-8/9 patch (commit `bd7c0dcf`, identical 24 files, +2626/-152)**. All further work is an *increment on top of `bd7c0dcf`* (`feat_withIrr_increment.patch`), no longer a patch on `main`.
- T-13, T-17 done (above); `update_irrigation_map` added (Q16); `blas_lite.f` removed from this branch and the toolchain study moved to `NO_SYSTEM_DEPS.md` (separate branch).
- Regression: `main` 170 passed / 3 failed; branch 216 passed / the same 3 failed. `test.py` exit 0, no `fort.*` left in the project folder.

### Session 11: 2026-11-09
- Venezia level A reproduced numerically (section 14); Salerno scripts + smoke test; T-24...T-27 fixed in the wrapper. Regression: `main` 170 passed / 3 failed; branch 218 passed / the same 3 failed. Increment patch = `feat_withIrr_increment.patch` (on top of `bd7c0dcf`, applies cleanly).
- Salerno `irrigation_map` regenerated exactly from the repo tables (section 14). Next: Level C for Venezia; Salerno `atmbc` (port or owner's file) -> real comparison on a window; then `no-system-deps` branch if wanted.

<!-- Next session: Level C Venezia, Salerno atmbc; P-06 (apply patch), P-05 (Venezia/Salerno files) -> L4/L5; L3 examples 2+ (use dataset-appropriate time stepping); T-12/T-13; Q8. -->
