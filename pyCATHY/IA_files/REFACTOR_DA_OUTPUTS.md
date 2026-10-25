# Refactor: homogeneous DA outputs + light issue log

Scope: `pyCATHY/DA` (EnKF, PF = sequential path; ES-MDA = smoother path).
Baseline before any change: `pytest test_compact_states test_console_report test_esmda_cathy` -> 90 passed.

## 0. Diagnosis (what the code does today)

| Item | Sequential (EnKF + PF) | ES-MDA |
|---|---|---|
| Entry point | `run_DA_sequential` -> `_run_DA_sequential_impl` | `run_DA_esmda` -> `esmda_cathy.run_esmda_on_cathy` |
| `StateStore` used | yes (`_save_cycle_states`, via `compact_states=True`) | **no** (never imported in `esmda_cathy.py`) |
| `<project>_states/` (chunks) | only if `DAConfig.save_results=True` (default **False**) | never written by the library |
| `<project>_states.nc` | only if `save_results=True` (`_finalize_state_backup` returns early otherwise) | never written by the library |
| Console tee `console_log.txt` | everything on stdout/stderr, append, no size cap | same |
| Diagnostics | `_rejections.json`, `_ERT_checks.json`, 46 `logger.warning` calls, PF warnings returned in a list | `logger.warning` only |

Root causes of the inconsistency:

1. Writing the final states is **coupled to `save_results`** (a crash-safety option), so a run with the default config
   leaves no states at all. This is the most likely reason the EnKF/PF run has no `states` folder/file.
2. ES-MDA bypasses `StateStore`. If your ES-MDA runs do show a `states` folder, it comes from the notebook/script
   (not from this library); I could not see it.
3. Output paths are rebuilt with `os.path.join(workdir, project, project + "_xxx")` in ~8 places: no single contract.
4. `console_log.txt` is a raw tee of stdout (rich tables, per-cycle prints, member processor output): heavy and
   unfilterable. Real problems are scattered between it, 3 JSON files and the logger.

## 1. Decision: one layout for every DA type

```
<project>/
  <project>_run.json          manifest (outputs registry + issues summary)   [all DA types]
  <project>_states.nc         final states, one xarray Dataset               [all DA types, ALWAYS]
  <project>_states/           per-cycle npz chunks = crash safety only       [only if save_results]
  <project>_issues.log        light issue log (this refactor)                [all DA types, ALWAYS]
  <project>_rejections.json   unchanged
  <project>_ERT_checks.json   unchanged
  console_log.txt             only if DAConfig.console_log == "full"         [opt-in, was always on]
```

Why this one:

- `states.nc` is the single file plotting/post-processing should read (`load_states` already handles both forms).
- The chunk folder is an implementation detail of crash recovery, so it must not be what distinguishes DA types.
- Writing states must not depend on `save_results`.

## 2. Steps

- [x] **S1** `run_outputs.py`: `RunOutputs` dataclass = the only place that builds output paths. Replace the ad-hoc
      `os.path.join` in `cathy_DA.py`. DONE for states, chunks, rejections, ERT checks, issues, console log.
- [x] **S2** Decouple the final `states.nc` from `save_results` (new `DAConfig.write_states=True`); chunks stay tied to
      `save_results`. Applies to EnKF and PF (same code path).
- [x] **S3** `issues.py`: light issue log. A `logging.Handler` on the `pyCATHY` logger (WARNING+) writing one line per
      *distinct* issue (digits normalised, repeats counted, first/last cycle kept, hard cap on lines), plus
      `IssueLog.record(kind, ...)` for structured events. Output `<project>_issues.log`, summary dict for the manifest.
- [x] **S4** Structured detectors routed to the issue log:
      member rejection, ensemble collapse (spread ratio after/before of states and parameters = "inbreeding"),
      PF `N_eff` degeneracy / resampling, obs perturbation ~ 0, parameters out of bounds, rank-deficient update
      (n_param >= n_members), ERT check failures, psi clipping.
- [x] **S5** `DAConfig.console_log` = `"off"` (default) | `"full"` (legacy tee). Wire both `run_DA_sequential` and
      `run_DA_esmda` to the issue log.
- [x] **S6** `DA.run_summary()` returns `{outputs, issues, n_cycles, rejected}` so the notebook merges it into its run JSON
      (the manifest of the notebook is written outside the library). DONE: `DA.run_summary()`; the notebook still has to merge it.
- [x] **S7** Tests (no CATHY executable needed) + docs + changelog entry. Tests and docstrings DONE (`test_run_outputs_issues.py`,
      17 tests) + changelog entry (`doc/changes.md`, S15); no dedicated sphinx page.
- [x] **S8** ES-MDA states = states of the POSTERIOR forward run (`esmda_cathy.esmda_states_store`), same `StateStore` /
      `<project>_states.nc`, `analysis` = NaN, attrs `da_type="esmda"`, `states_source="posterior_forward_run"`.
      CORRECTION of section 5: no extra run is needed, `run_esmda(final_forward=True)` (default, and the only value the
      CATHY adapter uses) already runs the model with the final parameters; the member folders hold that run.
      Unit-tested with a fake simulation; NOT validated on a real CATHY run (reads `psi`/`sw` with `read_outputs`
      exactly as `predict_member_blocks` does).
- [x] **S9** `<project>_params.nc` for ES-MDA: every iterate in physical space [iteration, member(, dim)], predictions
      [run, obs, member], misfit, alpha, valid. Sequential DA does not write it (its parameters are in `dict_parm_pert`).
- [x] **S10** ES-MDA calls `backup_results_DA` (dict_parm_pert, df_obs, ET_DA_xr, df_Archie, meta_DA). An empty `df_DA`
      is never written any more. `df_performance` is still NOT produced for ES-MDA.
- [x] **S11** `meta_DA` is now written; `load_backup` loads `params`, prints ONE line when there are no states
      ("neither _states.nc nor a _states/ folder" / "this run wrote no states") instead of a stack of errors.
- [x] **S11b** ES-MDA issue log: spread collapse between iterates, reduced chi-square above 2 or increasing, members
      dropped (`esmda_issues`). Every ES-MDA output step is guarded: a failure there cannot lose the run.
- [ ] **S12 (nb8, script not in the repo, spec only)** `nb8_plots_synthetic.py`: per-run try/except with a ONE-line reason
      (no traceback dump), a final table `plotted / skipped (reason)`, capabilities by `da_type`
      (state-based scores need states; parameter figures need `_params.nc`), no abort of the other runs.
- [x] **S13 (run_twin.sh, DONE in the patched script)** `plots` pre-check: list the runs of the case with status (finished / missing states /
      missing params) before calling nb8, and exit non-zero when something was skipped.
- [~] **S14 (run_twin.sh DONE, real scripts unknown)** Real-data case in the same entry point:
      - `real` : DA of the real data, as today (`REAL_NB5`), but with a default/auto-detected script instead of an
        error when the variable is unset, and the same outputs layout as the twin runs (S1-S6).
      - `real-plots` : calls the real-data plotting script (`nb5_DA_ERT_soil3D_nc_plots`, name to confirm) on the finished
        real runs; it must load with the same `load_backup` / `load_states` as nb8 and skip with a one-line reason.
      - `plots` also runs `real-plots` when real runs exist (`REAL=0` to disable); `everything` gets an opt-in
        `WITH_REAL=1` (real data are slow and have no truth: they must never be mixed into the truth-based scores).
      - Real runs have NO truth: figures that need it (node errors, scores) must be skipped by capability, not crash.
      - Real-data ES-MDA has the same missing-files problem as section 5 (S8-S10 apply).


## 3. Review of the case `SYN_NC_ic_enkf_N64_s1_unbiased`

Consistent: 23/23 surveys assimilated, 0 skipped, no ensemble warnings, atm times equal truth, obs error 5 % equals the
truth noise 5 %, `log_obs` true with a lognormal multiplicative noise, runtime 40 min for N=64.

To check / improve:

1. `list_parm2update = ["St. var."]` and `n_updated_parameters = 0`: this is a **state-only** EnKF (no parameter
   updated). `param_bounds_mode = "reflect"` and `vgp_corr` are then irrelevant; record them as `null` to avoid
   suggesting they were active. Spread collapse must be monitored on the **states** (S4), not only on parameters.
2. `inflation = 1.1` with N=64 and a state-only update: with no parameter spread to maintain, collapse of the state
   ensemble over 23 analyses is the main risk. The manifest has no spread/ESS diagnostic: add one (S4/S6).
3. `store_full = "all"` stores every node at every cycle (float32): fine for this mesh, but record
   the resulting file size and `n_cycles` in the manifest.
4. `da_type = "enkf_Evensen2009"` but `argv` says `--da-type enkf`; `scenario = "NC_ic"` vs `--scenario ic`. Keep the
   resolved names, but store the raw CLI values under `argv_parsed` for reproducibility.
5. `archie_truth.a = 0.3`: unusual (Archie `a` is typically 0.6 to 1.0). Verify it is the intended convention of the code.
6. Reproducibility: absolute paths with a user e-mail (`argv`, `prior_cache`, `project_dir`, `../` not normalised).
   Store paths relative to `--root`. Add the git hash of pyCATHY and of the notebook, and the package versions.
7. Two seeds exist (`seed = 1` ensemble, `truth_noise.seed = 1234`): keep, but name them
   `ensemble_seed` / `truth_noise.seed` to avoid ambiguity.
8. Missing outputs registry: nothing points to `states.nc`, `rejections.json`, `issues.log`.
9. `status = "finished"` has no counts. Add `n_cycles`, `n_rejected_members`, `n_issues`, `wall_time_s`.
10. **`perturb_obs` is not in the manifest.** `AnalysisConfig.perturb_obs` defaults to `False` and its own docstring says
    this gives a spread smaller than the posterior and "the ensemble tends to collapse over the cycles"; `True` is
    recommended for `enkf_Evensen2009*`. This run is `enkf_Evensen2009` and the value is not recorded: check it, and add
    `perturb_obs`, `sakov`, `localisation`, `threshold_rejected`, `psi_bounds` to the manifest.

## 4. Status after the refactor

Verified: 472 tests pass (455 before + 17 new); the 9 failures and 15 errors are identical before and after
(they need the CATHY sources/data downloaded by the tests, absent in this sandbox).

NOT verified: a real EnKF / PF / ES-MDA run (no CATHY executable here). The hooks in the sequential loop
(`_record_analysis_issues` after `run_analysis`, `run_DA_sequential` reading `config` from `args[6]`/`kwargs`) are covered
by unit tests on shells, not by a full run. Run one short case and read `<project>_issues.log` before trusting it.

Update: S8-S11 done (see section 7). Still open: S12-S14 (nb8 / run_twin.sh / real data scripts are not in the repo).

## 5. Issue found with `./run_twin.sh plots` (ES-MDA runs cannot be plotted)

Symptom (log of `run_twin.sh plots`, case NC): the 2 EnKF runs load (`Reading states: ..._states.nc`), the 2 ES-MDA runs
(`SYN_NC_ic_VGP_esmda...`, `SYN_NC_ic_esmda...`) find NONE of `meta_DA, dict_parm_pert, df_obs, df_performance,
ET_DA_xr, df_Archie, states` and nb8 stops on `AssertionError: no states / df_DA ([])`
(`nb8_plots_synthetic.py:146`, `load_df`). `run_twin.sh` only forwards to nb8: it is not the cause.

Root causes in the library:

1. `backup_results_DA` is called only from `_run_window` (sequential). `run_esmda_on_cathy` never calls it: no
   `dict_parm_pert`, `df_obs`, `ET_DA_xr`, `df_Archie` for ES-MDA.
2. ES-MDA never feeds the `StateStore` (S8) -> no `states.nc`, and the loader falls back to the chunk folder, hence the
   misleading `no cycle_*.npz in .../_states`.
3. `meta_DA` is requested by `load_backup` but never written (noise on all runs).
4. I cannot tell from the log whether those two ES-MDA runs finished (nb7 writes a status in its run JSON): nb8 should
   read it and report "not finished" separately from "finished without states".

### Recommendation for ES-MDA states

Do both parts, not one:

- **Always** store the parameter ensemble of every iteration (S9): it is what a smoother estimates, it is tiny, and it
  is exact.
- **Plus one extra forward run** with the final posterior parameters, stored as states (S8). Do NOT store "the states of
  the last forward run": in ES-MDA that run uses the parameters from BEFORE the last update, so it is not the posterior.
  Plotting it next to the EnKF analysis (`compare-da`) would silently compare different things.
  Cost: +1 run of N members on top of the Na iterations (about +25 % for Na = 4).
- Without the posterior run, nb8's state-based scores (node errors against the truth: "scores on 1936 of 1936 nodes")
  are impossible for ES-MDA, and that is the main point of `compare-da`.

## 6. Real data in `run_twin.sh` (current state, from the script as uploaded)

| Command | Real data today |
|---|---|
| `real` | runs `$REAL_NB5` (your nb5 / nb6_ESMDA script) only; fails if `REAL_NB5` is unset |
| `everything` | NC + HC synthetic only, no real case |
| `plots` | `nb8_plots_synthetic.py` only, no real-data plotting |

Open questions (I only have `run_twin.sh`, not the nb5 scripts):
1. Exact name/arguments of the real-data plotting script (`nb5_DA_ERT_soil3D_nc_plots`?).
2. Real project folder/naming (no `SYN_` prefix?) so the runs can be discovered like the twin ones.
3. Which DA types are run on real data (EnKF, PF, ES-MDA).

## 7. Status after ES-MDA homogenisation

Tests: 480 pass (455 before the whole refactor); same 9 failures + 15 errors as the untouched repo (environment).
Files of a run, now for EnKF, PF and ES-MDA: `_states.nc`, `_issues.log`, `_dict_parm_pert`, `_df_obs`, `_ET_DA_xr`,
`_df_Archie`, `_meta_DA`; ES-MDA also `_params.nc`.

Differences that remain, on purpose or not yet fixed:
- `df_performance`: written by the sequential loop only.
- `_params.nc`: ES-MDA only.
- the empty `<project>_states/` folder seen for ES-MDA in the `plots` log is not created by this library (only
  `StateStore(chunk_dir=...)` creates it, used by the sequential path): probably nb7. To check in nb7.
- the sequential hooks and the ES-MDA states reader have not run on a real CATHY simulation.

## 8. Second `run_twin.sh plots` log: what really failed

The `File not found` lines were NOT what stopped the plots; they are noise printed by `load_backup`, one block per file
and per run. The plots stopped on:

```
nb8_helpers.py:148  row['rmse_ol'] = float(sc_ol['rmse'][sel].mean())
IndexError: Boolean index has wrong length: 23 instead of 362
```

`sel` is a boolean mask of length 23 (the 23 survey times of the DA scores) applied to the open-loop scores, which have
362 rows (every output time). It happens on the first run with an open loop (`SYN_NC_ic_VGP_enkf`), so it is independent
of ES-MDA. The fix belongs in `nb8_helpers.summarize_run`: select the open-loop rows by TIME value
(e.g. `sc_ol[sc_ol[time_col].isin(sc.loc[sel, time_col])]`), not by position. I do not have `nb8_helpers.py`: not patched.
Better still, store the open loop at the same cycles as the DA runs (nb7 / `--no-da` path, not in the repo).

The `File not found` lines, run by run:

| Run | What the log shows | Cause | Status |
|---|---|---|---|
| all runs | `meta_DA` not found | `load_backup` asked for it, nothing wrote it | writer added (new runs only) |
| `SYN_NC_ic_VGP_esmda` | nothing at all, empty `_states/` | no results: not finished / crashed / never started (the log cannot tell) | `status` command names it; nb8 skipped it correctly |
| `SYN_NC_ic_esmda` | has `dict_parm_pert` + `states.nc`, lacks `df_obs`, `df_performance`, `ET_DA_xr`, `df_Archie` | ES-MDA never called `backup_results_DA` | fixed for new runs (S10); `df_performance` still missing |

Note: your `SYN_NC_ic_esmda` already has a `_states.nc`, which the clone I analysed does not write for ES-MDA. Either nb7
writes it, or your installed pyCATHY is not this version: check `python -c "import pyCATHY; print(pyCATHY.__file__)"`
and apply the patch to THAT tree.

Done in this round:
- `load_backup` prints ONE line per run (`Loaded <run>: ... | missing: ... | errors: ...`); `verbose=True` for the detail.
  It no longer tries the chunk folder when it is empty (that produced the duplicated `no cycle_*.npz` error).
- `pyCATHY.DA.run_status` (`python -m pyCATHY.DA.run_status <simulations> --glob "SYN_NC_*" [--strict]`): verdict per run
  (`ok`, `partial`, `no states`, `no results`, `failed`), states kind, issue counts, missing files; reads the run json status.
- `run_twin.sh`: `status` command; `plots` prints that table first and still runs nb8, then `real-plots` when
  `REAL_PLOTS` is set (`REAL=0` disables); `everything` runs `real` + real plots with `WITH_REAL=1`; nb8's exit code is kept.

## 9. Round S15: library manifest `<project>_run.json`

Gap found: `RunOutputs.run_json` was documented as written by `DA.write_run_summary`, which did not exist; `run_summary()`
lacked the fields asked in section 3 (items 3, 8, 9, 10 and the state-only note of item 1).

- [x] **S15** `DA.write_run_summary(status, extra)`, called by `_run_logging` at the end of EVERY run (EnKF, PF, ES-MDA),
      also on failure (`library_status="failed"`; the original exception is never hidden). Atomic write (tmp + replace).
      Content: `outputs` (now incl. `run`), `issues`, `n_issues`, `n_cycles`, `n_stored_cycles`, `rejected_members`,
      `n_rejected_members`, `states_nc_bytes`, `wall_time_s`, `versions`, and `options` = `da_type`, `n_members`, `seed`,
      `inflation`, `threshold_rejected`, `psi_bounds`, `localisation`, `use_log_transformed_obs`, `store_full`, `save_results`,
      `write_states`, `console_log`, `perturb_obs`, `sakov`, `n_updated_parameters`, `updated_parameters`, `state_only`.
      When `state_only`, `param_bounds_mode` and `vgp_corr` are recorded as `null` (section 3, item 1).
      When `DAConfig.analysis` is None, `perturb_obs`/`sakov` are recorded as "AnalysisConfig default (analysis=None)"
      (the effective default is `perturb_obs=False`: check item 10 for `enkf_Evensen2009`).
- Key is `library_status`, NOT `status`: the notebook manifest keeps its own `status`; `run_status` reads the script's first
  and falls back to `library_status`.
- Tests: 4 new (109 in the 4 DA test files). Full suite: 489 passed, 9 failed + 15 errors identical to the untouched repo.
- Not done (outside the library): notebook-side items of section 3 (relative paths, git hash of the notebook, `argv_parsed`,
  seed naming), S12 (nb8), S14 real scripts, `nb8_helpers.summarize_run` fix (section 8), `df_performance` for ES-MDA.
  Still unverified on a real CATHY run.

## 10. Round S16: ES-MDA chunks, `states_kind`, nb8 homogenised

Question asked: StateStore + `_states/` chunks + `_states.nc` for ES-MDA too, and one processing/plotting pipeline.
State before this round: StateStore and `_states.nc` already existed for ES-MDA (S8); chunks did not.

- [x] **S16a** ES-MDA writes `<project>_states/` chunks under the SAME rule as the sequential DA (`DAConfig.save_results`),
      through `RunOutputs.fresh_chunk_dir()` (also used by the sequential path; removes stale chunks).
- [x] **S16b** `StateStore` chunks now persist `attrs` (`attrs_json`): a chunks-only reload keeps `da_type`, `states_kind`,
      `states_source`. Old chunks without it still load.
- [x] **S16c** `states_kind` attribute in `states.nc` / chunks: `forecast_before_update` (EnKF, PF: `psi`/`sw` are the
      values BEFORE the update at each cycle) vs `smoothed_posterior` (ES-MDA: posterior forward run, conditioned on ALL the
      data incl. future ones). `description` of the ES-MDA file no longer claims "before the update".
- [x] **S16d** nb8_helpers: `summarize_run` selects the open loop by TIME value (fixes `IndexError: Boolean index has wrong
      length: 23 instead of 362`); `nearest_index`; `parameter_evolution` / `parm_evolution_ds` / `parameter_table(params_ds=)`
      read every ES-MDA iterate from `_params.nc` (sequential: `dict_parm_pert`).
- [x] **S16e** nb8_plots_synthetic (S12): OL drawn on its own time axis (states) / matched by time (profiles); ONE-line reasons
      (`--debug` for tracebacks); final table plotted / skipped per run; exit 1 only for real skips (`n/a` reasons, e.g. a
      state-only run has no parameter figure, do not count); `states_kind` column in summary.csv and a warning in the heatmap
      title when forecast and smoother states are mixed.
- Verified: 4 new tests (ES-MDA chunks, stale chunks, attrs round trip, default kind); nb8 smoke-tested end to end on a
  synthetic project (EnKF run + ES-MDA run + 362-time open loop, stub loader). NOT run on real CATHY outputs.
- Open: a real crash-safety for ES-MDA would be a per-ITERATION checkpoint of the parameters (today `_params.nc` is written
  once at the end; the state chunks are written after the posterior run, so they only protect the NetCDF write).
  Like-for-like EnKF vs ES-MDA state scores would need the sequential `analysis` field (psi only) or a smoother-style
  re-run of the EnKF posterior: not done.

## 11. Verdict on section 5 (ES-MDA runs cannot be plotted), checked against the code

| Root cause of section 5 | Status | Where |
|---|---|---|
| 1. ES-MDA never called `backup_results_DA` | fixed for NEW runs | S10 |
| 2. ES-MDA never fed the StateStore / no `states.nc` | fixed for NEW runs (+ chunks with `save_results`) | S8, S16a |
| 3. `meta_DA` never written | fixed for NEW runs | S11 |
| 4. nb8 must separate "not finished" from "finished without states" | fixed in this round: unfinished runs are named with their status in the final table, exit 1 | S16f |
| Recommendation "one extra forward run with the posterior" | NOT needed: `run_esmda(final_forward=True)` (default, the only call of the CATHY adapter) already runs the final parameters; verified in `esmda.run_esmda` | S8 |

Correction made in S16f: `states_kind` was `smoothed_posterior` unconditionally. It is now `last_iteration_forward_run`
(description: "NOT the posterior") when the final forward run is absent, so such states can never be plotted next to the EnKF as if
they were posterior.

Still true: existing ES-MDA runs made BEFORE these patches are not repaired (they lack the files); they have to be re-run.
`df_performance` is still not written for ES-MDA but it is not in `run_status.REQUIRED` and nb8 does not read it.
NOT verified on a real CATHY run (reads `psi`/`sw` with `read_outputs`, unit-tested with a fake simulation only).

