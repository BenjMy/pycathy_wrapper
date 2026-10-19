# Ensemble crash report: pyCATHY DA with ERT (`DA_ERT_part2_runDA.py`, scenario `ic_ZROOT_upd_ZROOT`)

Date of run: 16/10/2026, 13:19. Error: `ValueError: % number of rejected ensemble is too high: 84.375` (27 of 32 members rejected), raised in `_raise_if_too_many_rejected` at the start of window 2 (atmbc update 2).


> **Status (patch `pycathy_fixes.patch` on `43f7126`).** **Legend: ✅ = implemented in the patch and covered by a synthetic unit test; it has NOT been validated on the real `ic_ZROOT_upd_ZROOT` run.** Items without ✅ are open or only partly done. 255 synthetic tests pass (numpy only; `perturbate()` is also exercised on a stub simulation). The 27/32 rejection is therefore not yet confirmed fixed. A progress tracker is in section 15; new findings are in section 8 (B13 to B37).

## 1. TL;DR

The run did not crash because of one bad member. It crashed because the **first analysis step destroyed the ensemble**:

1. The prior is already poor. `ic ~ N(-1, 1.75)` gives about 28% of members a **positive pressure head**, and 6 of 32 members fail in window 1.
2. The EnKF analysis (`enkf_Evensen2009_Sakov`) **drops the HPH^T term**, so the gain is not a proper Kalman gain, and 2415 fully correlated ERT data points make it overshoot by orders of magnitude.
3. The update pushed `ZROOT0` from 0.4 ± 0.005 m to values between -5.6 and +17.5 m.
4. `_check_after_analysis` only checks the **upper** bound of ZROOT. Negative root depths pass, CATHY writes NaN in `mbeconv['CUM.']`, and 14 more members are rejected in window 2.
5. 6 + 7 + 14 = 27 rejected, so 84.4% > the 80% threshold.

## 2. Evidence from the log

### Member bookkeeping

| Stage | Members | Reason |
|---|---|---|
| Window 1, before analysis | 5, 15, 20, 22, 23, 25 rejected (6) | CATHY stopped at t = 0.02 to 0.07 s (time-step collapse); member 20 has a malformed `mbeconv` (shape (19,1) vs (19,19)). All six have `ic >= 0` (0.51, 0.98, 0.02, 1.00, 0.58, 0.58), meaning positive pressure head with `IPOND=0`. |
| After analysis 1 | 1, 7, 9, 11, 19, 24, 30 dropped (7) | `unfeasible root depth ZROOT0` > mesh thickness (17.46, 17.48, 6.87, 17.47, 2.95, 12.5, 13.8). Leaves 19/32 valid. |
| Window 2, before analysis | 2, 3, 4, 6, 8, 12, 14, 16, 18, 21, 26, 28, 31, 32 rejected (14) | `mbeconv['CUM.']` is NaN. **These are exactly the 14 members whose updated ZROOT0 is negative** (-0.20 to -5.56). Survivors 10, 13, 17, 27, 29 have ZROOT0 between 0.15 and 1.50. |
| Total | 6 + 7 + 14 = 27 | 27/32 = 84.375% |

### The ZROOT update is a monotonic function of ic

`update_nb1` for ZROOT0, ordered by prior `ic`:

- ic of about 1.6 to 2.05 gives ZROOT0 of about +17.47 (members 11, 1, 7)
- ic = -0.44 gives +6.87 (member 9)
- ic = -0.93 gives +2.95 (member 19)
- ic of about -1.3 to -1.5 gives +0.15 to +1.1
- ic of about -1.7 to -2.9 gives -0.2 to -3.0 (negative)
- ic = -5.0 gives -5.56 (member 6)

The prior ZROOT spread (0.005) is irrelevant. The update is just `c * innovation_i`, with `c` far too large. ZROOT has almost no physical influence on a 24 h ERT survey of a (probably) bare, homogeneous soil, so `cov(ZROOT, y)` is pure sampling noise. A correct gain would leave it almost untouched.

### Other signals in the log

- `rank1_ensemble: 100.0% of the predicted-data variance is in one mode`. This is expected, since the whole state is a function of one scalar (homogeneous `ic`).
- `k_too_large` (55 quadrupoles, idx 2360 to 2414) appears in **the observation files** (meshLi0..9.data) and in every forward run.
- `rhoa_implausible` (about 100 to 116 quadrupoles) and `signal_below_noise` (360 to 880 of 2415 data points, i.e. 15 to 36%).
- Forward rhoa range is e.g. 0.139 to 16.1 for a homogeneous model of 14.7 to 17.6 Ohm.m. A homogeneous model should give rhoa = rho everywhere, so a third of the data is noise-dominated.
- `Assimilating all: 2415 data`. All of them are assimilated, including the corrupted ones.
- `Valid ensemble members: 19/32 (59%)` after the first update.
- `threshold_rejected=80` hides the problem. With the default (10) the run would have stopped in window 1.
- Warning `Could not copy processor files into cathy_origin (FileExistsError ... prepro)`. The project folder is reused from a previous run (see B7).

## 3. Probable bugs / design flaws (ranked)

### A. Critical (explain the crash)

✅ **A1. "Sakov" analysis omits HPH^T** *(fixed: exact ensemble-space solve, `enkf._solve_obs_space`)* (`enkf.py::enkf_analysis`, `if Sakov:` branch).

```python
COV = data_cov.T
inv_data_pert = data_pert / np.diag(COV)[:, None]
```

This computes `K = P_xo R^-1` instead of `K = P_xo (S S^T/(N-1) + R)^-1`.

- It is only valid when HPH^T << R, e.g. a single observation.
- With N_obs = 2415 and strongly correlated obs, the effective gain is inflated by about `1 + N_obs·σ_s²/R`. With σ_s of about 20 Ohm.m and R = 25 this is of order 10^4 (back-of-envelope estimate, please verify).
- Even for the "true" Sakov scheme, the observation-space inverse must stay in ensemble space (e.g. Woodbury / ETKF form). Diagonal R alone is not it.
- Most likely a primary cause of the ZROOT blow-up. It also affects the states (psi), which can become positive or unphysical.

✅ **A2. No lower bound (and no projection) on updated parameters** *(fixed: lower bound in `_nonphysical_values`, clip after analysis)* (`cathy_DA.py::_nonphysical_values`, `_check_after_analysis`).

- ZROOT is only tested with `values > root_depth_limit`. Negative values pass.
- `clip_min=0, clip_max=1.5` exist in `dict_parm_pert['ZROOT0']`, but are applied only at sampling and never after the analysis.
- Rejecting members (instead of clipping or transforming) collapses the ensemble: 32 to 26 to 19 to 5.
- `transf_type = None` in `scenario.toml`. A `log` or `log-ratio` transform (A=0, B=1.5) would keep ZROOT positive by construction.

✅ **A3. Prior ensemble contains non-physical `ic`** *(fixed: `per_bounds` honoured for `ic`, truncated sampling; bounds set in `scenario.toml`)* (`scenario.toml`, `ic` sd = 1.75, mean = -1).

- About 28% of draws are positive pressure heads, with no `clip_max`. These cause time-step collapse (members 5, 15, 20, 22, 23, 25) or fully saturated states (members 1, 7, 11 with res0 = 14.7).
- ic of about -5 m (member 6) is also extreme.

**A4. Updating a parameter that is not observable by the data.**

- ZROOT in the ERT window has no (or negligible) sensitivity. With any gain error, spurious correlations produce large updates.

### B. Important (degrade or hide the problem)

- ✅ **B1. Observation error unit mismatch.** *(fixed, opt-in: `meta["data_err_unit"]="percent"`)* `dataErr = 5` is documented as "error in %", but `make_data_cov` uses `ones * data_err**2`, i.e. an absolute σ = 5 Ohm.m for every datum (5% of 100 Ohm.m but 33% of 15 Ohm.m). Likely wrong for rhoa.
- ✅ **B2. Corrupted quadrupoles are assimilated.** *(fixed: flagged and excluded from the filter only, see section 12)* `k_too_large`, `rhoa_implausible` and `signal_below_noise` data should be masked in both obs and prediction (obs and forward model use the same sequence). Absolute noise 1e-4 V swamps small R = rhoa/k.
- **B3. Rank-1 ensemble.** One homogeneous scalar `ic` means the analysis can only apply a uniform shift. 2415 observations give no new information beyond one number.
- **B4. Possible time off-by-one** *(structure established by code reading, see B34; whether the data generation uses the same convention is still to verify with the real data, T9)*
- ✅ **B5. Wide use of `threshold_rejected=80`** masks early failure. Should log which stage rejected whom. *(fixed: a rejection ledger `<project>_rejections.json` records cycle, stage, member and reason; the error message lists the members by stage; the script uses 30)*
- ✅ **B6. `self.ens_valid = list(set.intersection(...))`** *(fixed: sorted)* in `update_pert_parm_dict` is not guaranteed sorted. Wrap in `sorted()`.
- ✅ **B7. Stale project folder.** *(handled in the run script, which deletes the folder before each run (T0); the library still copies only if the folder is absent)* `DA_ERT/<project>/DA_Ensemble/cathy_origin` and the `cathy_N` folders already existed (FileExistsError). `_create_subfolders_ensemble` only copies if the folder is absent, so old `input/`, `output/` (psi, mbeconv) and `prepro` may be reused. A failed member may then read **old** mbeconv or psi files.
- ✅ **B8. `read_mbeconv` fragility.** Member 20 raised `ValueError: Shape of passed values is (19, 1), indices imply (19, 19)`. *(diagnosed and fixed, see B33)*
- ✅ **B9. `_update_feddes_member` uses chained assignment** *(fixed with `_set_feddes_value`; verified on pandas 3.0.2 that the old form silently does nothing)* `df[col].iloc[i] = value`. Under pandas Copy-on-Write (default in pandas 3.x), this silently does nothing. The 14/14 match above suggests it works in your environment, but confirm the pandas version and that `input/soil` of member 2 really contains ZROOT = -1.30 after the update.
- ✅ **B10. `_check_positive_pressure_heads` is called on the forecast** (`ensemble_psi_valid`, before the analysis), not on the analysed psi. The warning therefore says nothing about the update's effect on psi. *(fixed: the analysed psi is checked too; the warning names the members, the maximum and the share of nodes. Positive heads are reported, not clipped.)*
- **B11. Archie / fluid conductivity assumed** (`Unknown fluid conductivity`, 588 µS/cm; a=0.3, m=n=1.7). If the synthetic observations were generated with different parameters, there is a systematic bias in H(x). Not a crash cause.
- ✅ **B12. `simu_solution.update_ic(...)` is applied to the solution project**, not to `simuWithDA` (script, after the file copy). Harmless since the DA overwrites ic per member, but probably unintended.

## 4. Tests to run (in this order)

Each test has a clear pass/fail criterion. Run on a cleaned project folder (see T0).

✅ **T0. Clean run directory.** Delete `DA_ERT/meshLi_withDA_ic_ZROOT_upd_ZROOT/` before each run. Pass: no `FileExistsError` warning; rejections are reproducible.

**T1. Isolate ZROOT.** Run scenario `ic_ZROOT_NOupdate` (listUpdateParm = ["St. var."]) with everything else identical. Pass/fail: if the crash disappears, the cause is the parameter update (A1, A2, A4). If it persists, the state update (A1) or the prior (A3) is responsible.

✅ **T2. Fix the prior only.** Clip `ic` to e.g. [-5, -0.05] (or use sd of about 0.5 or a bounded transform). Pass: 0 rejected members in window 1.

✅ **T3. Unit test of the Sakov branch.** Synthetic linear problem: x (n=50), H a random matrix (m=2000, correlated rows), R = σ²I, N=30 members. Compare `enkf_analysis(..., Sakov=True/False)` against the analytic Kalman gain. Pass: analysis mean error comparable to the full gain. Currently expect a blow-up in the Sakov case that scales with m. Also compare ZROOT-like unobservable parameters: their update should be about 0.

**T4. Print the analysis diagnostics for cycle 0.** Log `max|K·d|` for the states and the ZROOT row, and the ratio `trace(HPH^T)/trace(R)`. Pass: ratio is of order 1 (not 1e3 to 1e5) and the ZROOT update is within a few sd.

**T5. Use the full covariance form** (`enkf_Evensen2009`, not `_Sakov`), or an ensemble-space (ETKF/Woodbury) implementation. Pass: ZROOT0 after update stays within about ±5 sd of the prior; no member crashes.

✅ **T6. Bound and transform parameters.** Set `transf_type = "log-ratio"` with bounds (0, 1.5) for ZROOT, or clip after analysis. Add a lower bound check in `_nonphysical_values`. Pass: no negative ZROOT, no NaN in `CUM.`.

✅ **T7. Mask bad quadrupoles.** Remove the 55 `k_too_large` quadrupoles and those flagged `signal_below_noise` and `rhoa_implausible` from `dict_obs` data, `data_cov`, and the forward prediction. Pass: forward rhoa for a homogeneous model is within a few % of the model resistivity, and N_obs drops from 2415 to roughly 1300 to 2000.

✅ **T8. Fix R units.** Use σ_i = 0.05·|d_i| (or the reciprocal error `rec_err`) instead of 5 Ohm.m. Optionally set `DAConfig(use_log_transformed_obs=True, log_transform_cov=True)`. Pass: R diag range is consistent with the data range.

**T9. Time alignment.** *(partly done: each analysis logs the observation time and the model window, and the configuration is checked before the run, see B34; still to compare with how the synthetic data were generated)* Print the time of each observation, the window start/end, and the state time used in `map_states2Observations`. Check against the `sw_times` of the synthetic run. Pass: obs time == state time for every cycle; the last observation is assimilated.

**T10. Reproduce the member failures.** Run members 5, 15, 20, 22, 23, 25 standalone with their `ic` values and read `output/mbeconv` and the CATHY log. Pass: you know whether positive ic is the only cause (try ic = -0.05).

**T11. File check after analysis.** For member 2, read `DA_Ensemble/cathy_2/input/soil` and confirm that the ZROOT written equals `update_nb1[1]` (checks B9). Also read `input/ic` for min/max psi after the update.

**T12. Increase ensemble information.** Add perturbations that affect ERT (Archie parameters, Ks, spatially varying or layered `ic`) so that the data can constrain more than one mode (B3).

Optional diagnostics to add: log `ens_valid` after every stage and which check removed which member (`before_run`, `mbeconv`, `after_analysis`), and make `threshold_rejected` default to something small again once the run is stable.

## 5. Suggested fixes (smallest change first)

1. ✅ Clean the project folder, clip `ic`, and keep `threshold_rejected` low (T0, T2).
2. Do not update ZROOT with ERT only (T1). If it is required, bound it with a log-ratio transform (A2).
3. ✅ Replace the Sakov branch by the correct gain or an ensemble-space form (A1, T3, T5).
4. ✅ Mask corrupted quadrupoles and fix the observation error units (B1, B2).
5. ✅ Add a lower bound to `_nonphysical_values` for ZROOT and an `ic > 0` check, and clip instead of reject where possible.
6. Sort `ens_valid`, log per-member rejection reasons, and verify the time alignment (B4, B6).

## 6. Files to give the next AI

`DA_ERT_part2_runDA.py`, `scenario.toml`, `cathy_DA.py`, `enkf.py`, `observations.py`, `ERT_checks.py`, `mapper.py` (not provided yet), and the **full log of a clean run**. Also useful: `meshLi_withDA_ic_ZROOT_upd_ZROOT_ERT_checks.json`, `DA_Ensemble/cathy_5/output/mbeconv`, and `cathy_20/output/mbeconv`.

Suggested prompt: "Using the report above, implement T0 to T8 as small patches. First write the unit test for `enkf_analysis` (T3) and show that the Sakov branch fails, then fix it. Do not touch the CATHY Fortran. Report which ensemble members are rejected at each stage after each patch."

## 7. Confidence

| Claim | Confidence |
|---|---|
| Negative ZROOT causes the 14 window-2 rejections (exact 14/14 match) | High |
| Positive `ic` causes the 6 window-1 rejections (6/6 have ic >= 0) | High |
| Sakov branch makes the gain far too large (code reading plus the size of the update) | High, needs T3 to quantify |
| R unit mismatch (B1), time off-by-one (B4), pandas CoW (B9), stale folders (B7) | Plausible, unverified |


---

## 8. New findings in round 1

### ✅ B13. Observations are never perturbed (spread collapse), measured
`prepare_observations` carries a comment saying the filter is an ETKF (Bishop et al., 2001) that does not need perturbed observations, and `perturbate_obs()` is an unimplemented stub returning `None`. But the update actually applied is the stochastic-EnKF form (`X + P_xo (S S^T/(N-1) + R)^-1 d`, with the same `d` for every member), not a square-root transform. For that form the analysis spread equals the posterior spread **minus `K R K^T`**.

Measured on a linear-Gaussian test (n=5, m=40, N=4000, exact Kalman reference):

| Case | trace of analysis covariance |
|---|---|
| exact posterior | 0.650 |
| current code (no perturbation) | 0.095 (deficit 0.555 = `tr(K R K^T)`) |
| `perturb_obs=True` | within 10% of exact (unit test) |

Consequence: the ensemble loses spread at every cycle, which also feeds back into the gain of the next cycle. This is independent of the crash, but it will matter as soon as the run survives window 2.
**Implemented, opt-in** (default unchanged): `AnalysisConfig.perturb_obs`, passed through `DAConfig(analysis=AnalysisConfig(perturb_obs=True))`. Recommended for `enkf_Evensen2009*`. Alternative: a real ETKF / square-root update. Decision needed from the author.

### ✅ B14. `enkf_analysis_inflation` and `enkf_analysis_inflation_multiparm` computed a wrong ensemble mean *(fixed)*
`observation.reshape(m, N)` was applied to the transposed `(N, m)` array, which scrambles memory order. With `alpha = 1` these filters differed from `enkf_analysis` by 0.24; after the fix by 3e-15 (regression test added). Anyone who used these two `DA_type`s before got wrong analyses. `enkf_analysis_inflation_with_parm_cov` had the same line and was fixed too (it is unreferenced).

### ✅ B15. `DA/performance.py` is dead code and broken *(deleted)*
No module imported it (only the generated `egg-info` listed it); its functions already existed in `cathy_DA.py`, and it used `os` without importing it. The file is removed.

### B16. Unreferenced and hard-coded experimental filters
`enkf_assimilation_with_parm_cov` and `enkf_analysis_inflation_with_parm_cov` are not reachable from `ANALYSES`, and the first contains "example values" (`sigma_Ks = 1e-5`, `sigma_ZROOT = 0.1`) that are computed but never enter the gain. Remove or finish.

### ✅ B17. Mutable default arguments *(fixed in `DA/` and `ERT/`; 49 remain in other modules)*
`perturbate_parm(mean=[], sd=[], ...)`, `perturbate_parm_by_layers(mean=[], sd=[], nominal=[])` and 7 in `petro_Archie.py` (56 `B006` in the repo). Not a known cause of a bug here, but a function that later mutates one of them would leak state between calls.

### ✅ B18. `perturbate(..., nlayers=N)` crashed *(fixed)*
`if 'nlayers' in kwargs: nzones = kwargs.pop('nlayers')` assigned to the wrong variable, so `nlayers` was undefined and the first per-layer loop raised. One-word fix.

### ✅ B19. Zone localisation shifted by one when the map contains nodata *(fixed, verified on a toy map)*
`_analyse_by_zones` enumerated `np.unique(mask)` **including** `-9999` but indexed `mesh3d_nodes_valid`, which is built **without** it. On a 3x3 map with nodata and zones 1 and 2, zone 1 was given the nodes of zone 2 (8 nodes instead of 4) and zone 2 raised `IndexError`. The same index also selects the zone's parameter in `list_parm2update[loci + 1]`, so parameters were shifted too. Fixed with `_valid_zone_values()`; test added. Applies to `veg_map` and `zones` localisation whenever the map has nodata; a map without nodata was never affected.

### ✅ B25. `per_bounds` was parsed in two incompatible ways, and errors were swallowed *(fixed)*
`check4bounds` read `{min, max}` inside a bare `except: pass`, while `process_perturbation_block` read only a two-element list and ignored dicts. Demonstrated: a `[0, 1.5]` list, a `{mn, mx}` typo, and a shared `{min, max}` used with `het_size > 1` (several layers, zones or Feddes rows) **all returned `(None, None)`, i.e. no bounds and no message**. In your single-vegetation run the `ZROOT` bound happened to work (`het_size == 1`); with several vegetation classes it would have been dropped silently. Now one parser (`parse_bounds`) accepts `None`/`"None"`, `{min, max}` (either key optional) and `[min, max]`, shared or per entry, keeps the defaults for an open side, and raises `ValueError` naming the parameter for anything malformed or with `min >= max`. A `per_bounds` list shorter than `per_name` logs a warning instead of failing.

### ✅ B26. Key naming depends on the parameter family, and `listUpdateParm` was not checked against it *(fixed)*
Built on a stub with the real `perturbate()`: a scalar global `Ks` is keyed **`Ks`**, a per-layer `Ks` is `Ks0`, `Ks1`, and a scalar `ZROOT` is **`ZROOT0`**. `listUpdateParm = ["Ks0"]` for the scalar `Ks` matched nothing and was dropped without a message: the parameter was perturbed but never updated. `perturbate()` now checks every update name against the keys it actually built (`check_update_parameters`) and fails with a suggestion (`'Ks0' ... did you mean 'Ks'?`). The documentation must give the naming rule per family (to do, 14.2).

### ✅ B27. `root_map` entries overlap when a vegetation class spans several zones, and the later write wins *(reported by `describe_ensemble`, behaviour unchanged)*
With classes numbered 1 to 3 over two zones where class 2 covers both, `Ks0` (class 1, zone 1) is completely overwritten in zone 1 by `Ks1` (class 2), because entries are written in order into the same rows of the `(zone, layer)` soil table. The perturbation of `Ks0` then has **no effect**, and updating `Ks0` changes nothing. `root_map` and `zone` controls on the same parameter family, or `layers` combined with `zone`, can overlap the same way. The table of `describe_ensemble` now lists these overlaps. Whether the intended behaviour is "last wins", "average" or an error is a design decision for the author.

### B28. The pushed code lacks several changes that `IA_files/REFACTORING.md` records as DONE
Checked on `43f7126` (the remote `HEAD`), with markers and counts: **R4g** (reseeding) and **R7** (`NotImplementedError` in `sensors_measures`) were missing, the invalid-escape fixes were missing in `perturbate.py` and `cathy_plots.py`, and the `print` to logger conversion is absent for `cathy_tools` (46 `print` + 49 `console`), `cathy_utils` (10), `meshtools` (38) and `cathy_plots` (12), although the log lists them as converted. R4h (dense covariance) **is** present. Re-applied here: R4g, R7, the escape fixes, and the two prints of `perturbate.py`. **Decision needed from the author:** push the earlier local session changes, or ask for the four conversions to be redone (the `print2logger.py` tool named in the log is not in the repo). `REFACTORING.md` now carries a "Reality check" section and a session-3 table.

### B29. Undefined names (`F821`) that cannot run: 50 found, 5 fixed
✅ Fixed: `UNIT_DIVISOR` (now `cathy_utils.get_UNIT_divisor()`), `r2_score` (own `_r2_score`, no scikit-learn), `df_tensio`, `df_TDR`. Open: `sensitivity.py` (37) is a half-ported Morris/Sobol script whose functions use variables of the original script (`sample`, `SPP`, `Si`, `simu`, `morris_problem`); adding imports would hide the symptom. `plotters/cathyDA_plots.plot_perf` (8) relies on script globals. Both need a decision: rewrite with explicit arguments, or remove.

### ✅ B30. A `layers` scenario could not be sampled through the standard call *(fixed; pre-existing, not a regression)*
`perturbate_parm(..., parm=dp, ...)`, as called in the run script, raised `TypeError: 'NoneType' object cannot be interpreted as an integer` for any parameter with `pert_control_name = "layers"`: the layered sampler needs `nlayers`, which only existed as an undocumented `**kwargs` entry, and the entries built by `perturbate()` did not carry it. Checked on the untouched `43f7126`: same failure. The earlier example tests (E2) only exercised the building step, so this was never seen. Now `perturbate()` stores `nlayers` in each layered entry, `perturbate_parm` falls back to it, and all five layered branches raise a clear `ValueError` when it is really missing. The sampling path through `perturbate_parm` for layers is now tested.

### ✅ B31. `per_bounds` did not reach the prior sample of the layered and zoned parameters *(fixed; measured)*
Measured with `Ks`, prior sd 1e-5 and bounds [9e-5, 1.1e-4]:

| control | bounds stored in the entry | prior sample range | bounds respected |
|---|---|---|---|
| `layers` (before) | yes | 6.6e-5 to 1.33e-4 | no (the Sobol sampler had no notion of bounds) |
| `zone` (before) | **no, `(None, None)`** | 6.0e-5 to 1.30e-4 | no, and the post-analysis clipping (A2) had nothing to clip to |
| `root_map` (before) | no | unbounded | no |
| all three (after) | yes | 9.0e-5 to 1.10e-4 | yes |

So for a `zone` or `root_map` parameter, `per_bounds` was dropped silently, for the prior and for the A2 clipping. Fixes: (1) `_apply_scenario_bounds` stores the scenario bounds on every entry after the build, keeping the defaults an entry already has for an open side (e.g. ZROOT >= 0); (2) the layered sampler maps its Sobol points through the truncated quantile function; (3) a max-only bound is recognised (the samplers tested only `clip_min`); (4) `sampling_dist_trunc` accepts an open side and `sd = 0` (it divided by zero); (5) per-entry bound lists are resolved inside the layer and zone loops, and an open side no longer inherits the previous entry's bound. Bounds apply to the sampled distribution, i.e. to the parameter itself when `per_type` is `None`. Tests fail without the fix.

### B32. The combined `zone` + `root_map` + `layers` layout described in the code is not implemented
A comment in the soil-parameter branches says a scenario containing several controls "should contain a list of lists of lists". But `pert_control_name` is one string per parameter and each branch handles one control, so no combined layout exists. Worked example E5 of the plan (`zone` + `layers`, nested lists) therefore has nothing to test and is removed from the plan until the author decides whether to implement it. Today two controls on one parameter family need two different parameters; the same `per_name` twice is rejected by the validator because the keys would collide.

### ✅ B33. Member 20's obscure `mbeconv` error was a run that stopped after its first time step *(diagnosed and fixed)*
`np.loadtxt` returns a 1-D array when the file has exactly one data row, and building the table from it fails with the shape message seen in the log. Reproduced on a synthetic file and on the untouched reader. So member 20 had not produced a corrupted file: CATHY had written one time step and stopped (the time-step collapse seen for the other five members). `read_mbeconv` now reads one row correctly (the member is then rejected with the real reason by the existing checks) and raises `no time step` for an empty file. The rejection ledger (B5) records the reason either way.

### ✅ B34. Observation to analysis pairing: how it works, and a check before the run *(detection implemented; behaviour unchanged)*
Read from the code (`_run_window`, `_get_data2assimilate`, `petro_Archie.get_sw_ens_i`):
1. the model is run over **[t, t_next]** for `t` in `all_atmbc_times[:-1]`;
2. an analysis follows **if `t` is in `ENS_times`** (exact float comparison);
3. the data of the k-th analysis are the **k-th entry of `dict_obs`, chosen by counter, not by time**;
4. the predicted observation is computed from `df_sw.iloc[-1]`, the **last time of the member output**, i.e. the state at `t_next` (assuming the output ends at the end of the window, which CATHY writes).

Consequences: (a) an observation stamped `t` is compared with the model state at `t + window length` (the confirmed part of B4); (b) the observation at the final time is never assimilated; (c) **any observation time that fails the exact test (for example `3600.0000001`, or a time rounded as the script comment about "rounding to the nearest hundred" suggests) is skipped, and every later observation is paired with the wrong analysis**, silently. Whether (a) is wrong depends on how the synthetic data were generated (`tA=tt` with `tt` from the solution's `sw` times is the state at `tt`), which needs the real data.
`time_alignment.check_time_alignment(ENS_times, all_atmbc_times)` is a pure function that reports these cases (error for a skipped observation or unsorted times, warning for duplicates and for the one-window offset, with the offset in seconds). It runs at the start of `run_DA_sequential` and **only logs** (`observation times: ...`); it never stops a run. Tests include a mirror of the pairing loop. Opt-in alignment modes (assimilate at the window end, select the data by time) were deliberately not implemented: they change the assimilation timing and cannot be validated without a real run.

### ✅ B35. `alpha_VG` builds and samples but cannot be written to the members *(found while writing the guides; rejected up front)*
`perturbate()` builds `alpha_VG` and the sampler runs, but `_update_soil_params_from_perturbations` looks it up in `_SOIL_PARM_COLUMNS`, has no entry, and raises `ValueError: Unknown parameter 'alpha_VG'` when the member files are written, long after the scenario was accepted. My layouts test measured only build and sampling, so it was blind to this; it now also checks that every soil parameter the validator accepts has a column. The validator rejects `alpha_VG` with the reason (a conversion to `VGPSATCELL = -1/alpha` would be needed, a semantic change not made here; perturb `SATCELL_VG`). `VGP`, which is commented out in `perturbate.py` and builds nothing, is rejected too.

### ✅ B36. The Feddes parameters follow different layouts, now measured, validated and documented
60 layouts were run through the real `perturbate()` and `perturbate_parm` on a stub (`tests/test_scenario_layouts.py`), and the validator agrees with every outcome (it reports an error exactly when the real code fails; mutation-tested). Findings: `ZROOT` and `PCREF` take a scalar when the Feddes table has one row and **a list with one entry per row** when it has several (a per-class `ZROOT` works with that layout); `PCWLT`, `OMGC` and `PZ` take **a scalar only** (a list fails inside the sampler); `per_type` is a single value for all of them and `pert_control_name` is ignored; `porosity` and `ic` do not support `root_map`; `WTPOSITION` supports no control; a list without a control fails for every soil parameter and for `ic`. Checked against the untouched `HEAD`: 26 of 30 Feddes combinations behave identically, and the 4 differences are improvements (a scalar `PCWLT`/`PZ` crashed because the old truncated sampler did arithmetic on a `None` bound). `perturbate_dist` also raised an `UnboundLocalError` for an unknown `per_type`; it now says what is allowed.

### ✅ B37. The ERT quadrupole mask broke the covariance-localisation matrix *(my B2 change; fixed)*
The matrix has one column per observation, but the mask removes observations before the filter runs: with ERT, one excluded quadrupole made the localised filter fail with `Localization matrix shape (10, 20) does not match P_xo[:states, :] (10, 18)` (reproduced). `_mask_localisation_matrix` now drops the excluded columns (a matrix already the width of the kept data is left alone; any other width logs a warning). Also: `gaspari_cohn` returned `-2.8e-16` at exactly `2 L`; the weights are now clamped to [0, 1].

### Open, found while reading the layer, zone and localisation code (details in section 14)
- ✅ **B20. `root_map` assumes vegetation classes numbered 1..K without gaps.** *(now detected by `validate_scenario`, behaviour itself unchanged)* The scenario is indexed by position in `np.unique(veg_raster)`, but `_assign_soil_value` writes to `veg_map == index + 1`. A map with classes 1 and 3 silently applies the second entry to class 2 (which may not exist). Same assumption for `zone` perturbation (zone number = `index + 1`).
- ✅ **B21. Zone and class counts included the nodata value** *(fixed in `perturbate()` at its 7 counting sites with `_n_classes`; the new test fails with the old counting)*
- ✅ **B22. Zone-parameter pairing was positional and implicit** *(fixed: pairing by name, `scenario.select_zone_parameters`: `ZROOT0` belongs to zone 0, whatever its position in `listUpdateParm`; a zone without a parameter is analysed for the states only (the positional rule raised `IndexError`); the old positional rule remains only for names without a numeric suffix)*. (Suspected, not tested: with `localisation="nodes"` every analysis updates *all* parameters and the results are `hstack`ed; the shape of `analysis_param` for more than one parameter should be verified. The code logs `loc_shortcut` once.)
- ✅ **B23. Stale docstring.** *(rewritten, with the verified mapping rules)* `perturbate()` says "Perturbation per zone is not yet implemented", but `zone` and `root_map` controls exist in the code.
- ✅ **B24. Localisation conventions are undocumented** *(documented in the docstrings of `gaspari_cohn` and `build_localization_matrix` and pinned by `test_localisation.py`: weight zero beyond 2 L, `L > 0` enforced, matrix shape (n_grid, n_obs), column order irrelevant if consistent, add `z` for depth localisation; the extrusion of surface observations to all depths remains a design limitation)*: the Gaspari-Cohn weight is zero beyond `2 L` (so `L` is half the support, not the cut-off radius); `grid3d['mesh3d_nodes']` is ordered `[y, x, z]`; surface observations are extruded to **all depths** of the column (`mesh_nodes_local_valid`), which is wrong for ERT, whose sensitivity changes with depth.

### Also noted (not fixed)
- ✅ `sampling_dist` called `np.random.seed(1)` on every call, so two normal-sampled parameters got **identical z-scores** (a rank-1 parameter ensemble); the van Genuchten sampler did the same. Reseeding removed, as the author's own log (R4g) had recorded as done: call `np.random.seed(n)` once beforehand to reproduce (test `test_perturbate_seed.py`). `sampling_dist_trunc` already used the global generator.
- Updated states are still unbounded (psi can become positive after the analysis): `_check_positive_pressure_heads` only warns (B10).
- `enkf_dual_analysis` uses `COV[6]` of step 1 as the mean error variance of step 2; with the old Sakov branch this was `R`, with the full branch it is `HPH^T + R`. The meaning changes with the flag, check before using it.

---

## 9. Nomenclature and naming

✅ **Done for the filters** (`enkf.py`): the local variables of the six filters now use one vocabulary (rows `A`/`dA`/`dD`/`S`/`B`/`MeasAvg`/`COV`/`dAS` and the `augm_*`, `*_pert`, `obs_avg`, `data_pert`, `inv_data_pert` family of the table below). The table is also in the module as a comment. 267 tokens renamed in 6 functions, no argument and no public keyword touched, no collision with an existing name. Proven three ways: (1) the AST of the new file equals the original AST under the rename map, (2) all 17 filter/option combinations return bit-identical arrays, (3) `tests/test_enkf_golden.py` freezes the outputs of the original code (a mutation test confirms it fails if the solver is broken). **Not done**: the `Sakov` keyword and `DA_type` names (public, need a deprecated alias), `ens_valid`/`obs_valid`, `threshold_rejected`, the error-name family, file and module renames.

Principle: names say *what a quantity is*, not which paper it came from. Keep the standard EnKF symbols in docstrings, not in identifiers.

| Now | Proposed | Why |
|---|---|---|
| `A`, `Amean`, `dA` | `augmented_state`, `augmented_mean`, `state_anomalies` | `A` is also used for matrices elsewhere |
| `dD`, `data_pert` (innovation), `data_pert` (perturbed data in `make_data_cov`) | `innovations`, `perturbed_obs` | one name currently means two things |
| `S`, `obs_pert` | `obs_anomalies` | |
| `MeasAvg`, `obs_avg` | `obs_mean` | same quantity, 2 names (this is how B14 hid) |
| `COV`, `B`, `inv_data_pert`, `dAS` | `obs_space_cov` (`S S^T/(N-1)+R`), `weighted_innovations`, `cross_cov` | |
| `Sakov` flag, `DA_type="..._Sakov"` | `ensemble_space_solve` (keep `Sakov` as a deprecated alias) | it no longer means "drop HPH^T"; the name was misleading |
| `ens_valid`, `obs_valid`, `obs_mesh_nodes_valid` | `valid_members`, `valid_obs`, `valid_nodes` (0-based, sorted, always `ndarray[int]`) | one convention, documented |
| `mbeconv`, `CUM.` | keep (CATHY file names), wrap in `read_mass_balance()` | |
| `dict_parm_pert[name][update_key]` | small `ParameterEnsemble` dataclass (values per cycle, bounds, transform) | removes stringly-typed keys like `"update_nb1"` |
| `threshold_rejected` (a percent, default 10) | `max_rejected_percent` | the unit is in the name |
| `ERT_noise`, `dataErr`, `fwdNoiseLevel`, `data_err` | `obs_error`, `obs_error_unit`, `forward_noise_percent` | four spellings of related things |

Caution for the module renames: `ERT_checks.py` to `ert_checks.py` is a **case-only rename**. On case-insensitive filesystems (macOS and Windows defaults) the two names are the same file, so a shim named `ERT_checks.py` next to `ert_checks.py` collides, and `git mv` of a case-only change needs an intermediate name. Do the rename in two commits (`ERT_checks.py` to `_ert_checks_tmp.py` to `ert_checks.py`) and make the shim a differently named module or a package `__getattr__`, then test on macOS or Windows CI before releasing.

Module and file names (ruff `N999`: 16 modules; `N802`: 99 camelCase functions): `ERT_checks.py` to `ert_checks.py`, `petro_Archie.py` to `petro_archie.py`, `cathy_DA.py` to `da.py`. Do it with a one-release shim (`from .ert_checks import *` plus a `DeprecationWarning`) so notebooks keep working. Scenario keys (`per_name`, `per_nom`, `listUpdateParm`) are user-facing: rename only with a loader that accepts both.

## 10. Modularity

`cathy_DA.py` is 3831 lines with one `DA` class of 86 methods; 8 methods exceed 80 lines (`_DA_analysis` 167, `_analyse_by_zones` 101, `_run_DA_sequential_impl` 93). Suggested split, in this order (each step is a pure move, tests green before and after):

1. ✅ `analysis.py` (config and registry together): `AnalysisConfig`, `DAConfig`, legacy-kwarg handling, `AnalysisResult`, the adapters and `ANALYSES`, `run_analysis`. Moved without change (378 lines).
2. ✅ (partly) `param_tables.py`: `_SOIL_PARM_COLUMNS`, `_assign_soil_value`, `_nonphysical_values`, `_set_feddes_value`, `_valid_zone_values`, `_check_positive_pressure_heads`, `_update_key` (153 lines). Still in the `DA` class: `_clip_updated_parameter`, the back-transforms.
3. ✅ (merged with step 1) the `ANALYSES` registry and adapters. The filters stay in `enkf.py` / `pf.py`.
4. ☐ `da/obs_selection.py`: `_ert_quality_mask`, `_get_data2assimilate`, `_obs_key_select` (these are `DA` methods that use `self`: they need their inputs made explicit first).
5. ☐ `da/ensemble_io.py`: reading psi/sw/mbeconv, writing member inputs. This is where B7, B8 and B9 live, so give it its own tests with temporary folders.
6. ☐ `DA` becomes an orchestrator of the window loop only.

Also moved, without change: ✅ `runtime.py` (output tee, working-directory handling, logging set-up, 157 lines) and ✅ `performance_tables.py` (RMSE table and parquet storage, 80 lines). `cathy_DA.py` went from 3993 to 3384 lines. **How the move was proven**: a script cut the definitions by AST line ranges; afterwards (1) the AST of all 35 moved definitions equals the original, (2) every moved name is the same object when imported from `cathy_DA` (so no caller changes), (3) no undefined name in the five files, (4) the whole test suite passes. The new modules depend only on the standard library, numpy, pandas and `enkf`/`pf`, never on `cathy_DA` (guarded by `test_module_split.py`), and the 7 imports that became unused in `cathy_DA` were removed.

Inside `enkf.py`: six functions repeat the same 6 steps (means, anomalies, innovation, solve, cross-covariance, update). Extract `_anomalies()` and `_update()` and keep the six public functions as thin wrappers. Before touching them, freeze their current outputs in characterization tests (the `alpha=1` equivalence test is the first one; B14 shows why).

Interfaces: filters return a positional list of 11 items consumed by index (`result[9]`, `result[6]`); return a `NamedTuple` (the adapter already builds an `AnalysisResult`, so extend it). Methods that are plain functions with `self` (`performance.py`) go away with B15.

## 11. Linter, formatting and performance

### Linter baseline (ruff 0.17, config now in `pyproject.toml`)
Repo-wide with the new config: **1573 findings**, 851 auto-fixable. Top: `W293`/`W291` whitespace (736, safe to autofix), `N802` function names (99), `E402` imports (63), `F841` unused variables (58), `B006` mutable defaults (56), `F401` unused imports (52), **`F821` undefined names (50)**. The `F821` hits are latent `NameError`s and should be triaged first (4 were in `performance.py`, fixed; list the rest with `ruff check pyCATHY --select F821`).
Plan: (1) run `ruff format` and `--fix` for whitespace in one commit that does nothing else (so `git blame` stays usable via `.git-blame-ignore-revs`); (2) fix `F821`, `F841`, `B006`, `E722`; (3) add `ruff check` and `pytest -m "not needs_cathy"` to CI; (4) tighten the per-file ignores as modules are refactored. I did **not** autofix in round 1: it would touch hundreds of lines in files nobody could run.

### Tests
The suite cannot run without the CATHY executable and downloaded data (`test_mapping`, `test_DA_Noble`, `test_hydro_*` fail at collection even on the untouched checkout). Mark those `needs_cathy` and make the numerical core testable without them: the 24 new tests need only numpy and run in about 2 s.

### Performance (estimates from code reading, **not profiled**)
- ✅ **R as a dense `m x m` matrix per time step** (`make_data_cov`): 2415 data is 46.7 MB per cycle and 10 cycles were kept in `stacked_data_cov` (**466.6 MB**, growing as `cycles x m^2`). *(done, opt-in: `make_data_cov(..., diagonal=True)` keeps the diagonal, **193 kB** for the same case, i.e. `m` times less; `DA._data_cov_of_cycle` rebuilds the dense matrix of one cycle at analysis time (46.7 MB transient), so every filter and diagnostic receives exactly the matrix it always did, bit for bit, tested. The run script now passes `diagonal=True`. Default stays dense for compatibility. A fully diagonal path through the filters was deliberately not attempted: a vector reaching `S S^T/(N-1) + R` would broadcast along rows instead of the diagonal, and `np.diag(vector)` would build a matrix.)*
- ✅ **Localisation matrix**: `build_localization_matrix` looped over grid points in Python; now one `cdist`, identical output (max difference 0), about 2x faster on 20000 x 300.
- **Analysis solve**: previously `O(m^3)` for the non-Sakov path; now `O(m N^2)` when `m > N`, direct `m x m` otherwise (picked automatically, after my first version tried an `N x N` matrix and ran out of memory at N=20000).
- **Diagnostics**: `analyze_covariances` and the COV plots build `m x m` matrices; keep them behind `logger.isEnabledFor(DEBUG)` (partly done).
- **Ensemble means** use `np.tile`; `mean(axis=1, keepdims=True)` with broadcasting avoids the copies.
- **Member I/O**: `copy_files_recursively` of the full solution folder per member, and VTK files per time step. Profile before optimising; likely dominated by disk, not numpy.
- **Parallelism**: `parallel=True` runs members; `ERT_checks` and the forward ERT run per member per cycle. Profile one cycle (`cProfile` plus wall time per stage) to see whether forward ERT dominates.

## 12. Policy for the ERT quadrupole mask: flag, never delete

Requirement: **all quadrupoles stay in the data, in the outputs and in the reports; only the filter skips the bad ones, and each skip is traceable.** Implemented as follows, with the limits that remain.

1. `ERT_checks.quadrupole_quality_mask(df)` is a pure function: it reads the table and returns a keep-mask plus a count per reason. It never writes to the table (unit test compares the frame before and after).
2. `_DA_analysis` applies the mask to the **filter inputs only** (data, R, predictions together). The full-length data are still returned, so the performance metrics, the diagnostics tables and `dict_obs` keep all quadrupoles.
3. Every exclusion is recorded per cycle and per sensor in `<project>_ERT_checks.json` as `{"stage": "mask", "cycle", "sensor", "n_total", "n_excluded", "reasons", "excluded_idx"}`, with the **complete** index list (not truncated at 20). The record has no `ok` key, so it is not counted as a failed check.
4. If the stacked layout cannot be matched to the data (several sensors, localisation), nothing is excluded and a warning is logged. It never guesses an alignment: a mis-aligned mask would drop the *wrong* quadrupoles silently, which is worse than dropping none.
5. `DAConfig(mask_bad_ert=False)` restores the old behaviour. The default is `True`, which changes existing ERT results; say so in the release notes.

Residual risks to check on the first real run:
- The mask uses only observed `k` and `rhoa` (and `valid`). Quadrupoles that are fine in the observations but bad in the forward model are not flagged (same sequence, so unlikely, but `check_ensemble_forward` should be cross-checked).
- `k_outlier` is permissive (8 spreads, at least 0.5 decade): it removes gross outliers, not marginal ones. `signal_below_noise` and `rhoa_implausible` need a noise model and `res0`; they are **off** by default.
- The set of excluded quadrupoles can differ between cycles, so the number of data varies. R, S and the innovation are cut together, so this is consistent, but any code that assumes a fixed `m` across cycles must be checked.
- Never mask from the innovation (that is circular and biases the analysis toward the prior).

**Check on the next run:** `excluded_idx` of cycle 0 should contain the 55 `k_too_large` indices (about 2360 to 2414 in the log), and `n_total` should be 2415.

## 13. Next-run checklist (when the laptop is available)

1. Apply `pycathy_fixes.patch` to a clean checkout (then, optionally, `whitespace_only.patch`, and check it with `python verify_ast_identical.py` from the repository root: it must report the same number of identical files and no `DIFFERENT` line), copy the updated `scenario.toml` and script, run the tests that need no CATHY: `pytest pyCATHY/tests/test_enkf_analysis.py pyCATHY/tests/test_DA_param_bounds.py pyCATHY/tests/test_ert_quality_mask.py pyCATHY/tests/test_obs_error_variance.py pyCATHY/tests/test_scenario_*.py pyCATHY/tests/test_perturbate_*.py pyCATHY/tests/test_undefined_names.py pyCATHY/tests/test_package_compiles.py pyCATHY/tests/test_rejection_ledger.py pyCATHY/tests/test_localisation.py pyCATHY/tests/test_bounds_sampling.py` (255 tests: all the files named `test_*.py` added by the patch, `test_scenario_layouts.py` and `test_docs_examples.py` included).
2. First run `ic_ZROOT_NOupdate` (T1, control), then `ic_ZROOT_upd_ZROOT`.
3. Read, in the log: the `observation times:` lines at the start (an ERROR there means some observation is never assimilated or paired with the wrong analysis), the `Perturbed parameters:` table (which rows each parameter writes, and any `perturbation layout:` overlap warning), one `Analysis n: observation stamped t=... ; model window [...]` line per analysis (B4), rejected members per stage (expect none in window 1 with bounded `ic`), `quadrupoles excluded` line, any `clipped to` line for ZROOT0, ZROOT0 range after update 1 (expect within a few sd of 0.4).
4. If it survives: repeat with `AnalysisConfig(perturb_obs=True)` and compare the ensemble spread of psi between cycles.
5. Send back `<project>_ERT_checks.json`, `<project>_rejections.json` (who was rejected, in which cycle and stage, and why) and the full log.

### Confidence update

| Claim | Confidence |
|---|---|
| Sakov branch gain error scales with number of data (A1) | **Verified** on synthetic data: update of an unobservable parameter 0.09, 0.67, 8.9 for 10, 200, 2000 data (prior sd 0.005) |
| Missing perturbation shrinks the spread by `K R K^T` (B13) | **Verified** on synthetic data |
| Scrambled mean in the two inflation filters (B14) | **Verified** (0.24 to 3e-15) |
| The fixes remove the 27/32 rejection on the real case | **Not tested** (no real run yet) |
| Performance items in section 11 | Estimates, not profiled |


---

## 14. To do: documentation of perturbation, update and localisation by layer and by zone

Goal: someone who has never read `perturbate.py` can configure a layered or zoned run, know which parameters are updated where, and not hit B19 to B24. Everything below was checked against the code at `43f7126` unless marked *(to verify)*.

### 14.1 Deliverables (in order)

✅ **Written** as three reStructuredText pages in `doc/content/DA_guides/` (`perturbation.rst` for D1 to D4, `localisation.rst` for D5, `troubleshooting.rst` for D6), linked from `doc/index.rst`. The pages build with Sphinx under `-W` (no warning); the full site was not built (its gallery runs the examples). They are **tested**: `tests/test_docs_examples.py` extracts every `toml` example and runs it through `validate_scenario`, `perturbate` and the sampling, checks that the catalogue names every parameter, that the conventions stated for localisation are those the code enforces, and that the messages quoted in the troubleshooting table exist in the code. The support matrix is the table of `tests/test_scenario_layouts.py`.

| # | Page | Content | Done when |
|---|---|---|---|
| ✅ D1 | **Concepts: where a perturbation applies** | the four controls (`None`, `layers`, `zone`, `root_map`), one diagram of the soil table indexed `(zone, layer)` | a reader can say which rows of the soil table a given entry changes |
| ✅ D2 | **Scenario file reference** (`scenario.toml`) | every key, type, units, default, and the list-of-lists layout per control | each key has a table row and an example; a validator rejects unknown keys |
| ✅ D3 | **Parameter catalogue** | for each perturbable parameter: scenario name, soil column(s), allowed controls, units, default bounds, transform | generated from `_SOIL_PARM_COLUMNS` and the `perturbate()` entries so it cannot drift |
| ✅ D4 | **Updating parameters** | `listUpdateParm` syntax, `"St. var."`, bounds and clipping after the analysis, back-transform, what is written to the member input files | worked example showing prior, analysed value and the file written |
| ✅ D5 | **Localisation guide** | domain localisation (`nodes`, `zones`, `veg_map`), covariance localisation matrix, the interaction with parameters | a decision table (which one for which problem) and one runnable example each |
| ✅ D6 | **Gotchas and validation** | B19 to B24 as "rules", and the checks that enforce them | every rule has a failing-input test and an error message that names the rule |

### 14.2 Facts D1 to D4 must state (verified)

- **Naming**: a perturbed entry is `<name><index>`, e.g. `Ks0`, `porosity1`, `ZROOT0`; `index` is **0-based**. The index is parsed from the name (`re.split("(\d+)", ...)`), so a parameter name must not itself end in a digit.
- **Mapping of the index to the soil table** (`_assign_soil_value`, table indexed `(zone, layer)`):
  - `layers`: rows of layer `index` in every zone;
  - `zone`: rows of zone number `index + 1` in every layer (**zones are 1-based, the index is 0-based**);
  - `root_map`: the zones whose vegetation class is `index + 1` (`veg_map == index + 1`), all layers;
  - `None` or no index: every row.
- **Feddes (vegetation) parameters take another path.** `ZROOT`, and the other Feddes parameters, are not written through the soil table: `_update_feddes_member` sets `FP_map[<name>].iloc[index]`, so `index` is the **row position in the Feddes table** (vegetation class order, 0-based), not `veg_map == index + 1`. Document both mappings side by side, because `ZROOT0` and `Ks0` look alike and do not follow the same rule. (The assignment is chained `.iloc[...] =`, see B9.)
- **Scenario layout**: scalar control uses one number per parameter; `layers`, `zone` and `root_map` use one list entry per layer, zone or class (`per_nom`, `per_mean`, `per_sigma`, `per_type`; `per_bounds` is shared by the entries of a parameter). `pert_control_name` is a list aligned with `per_name`. Combinations (`zone` + `root_map` + `layers`) need nested lists: document one full example, because the layout is easy to get wrong.
- **Soil columns written** (`_SOIL_PARM_COLUMNS`): `Ks` writes `PERMX`, `PERMY` and `PERMZ` together (isotropy is imposed); `porosity` to `POROS`; `n_VG`, `thetar_VG`, `SATCELL_VG` to the van Genuchten columns. Say it explicitly: users will otherwise assume `Ks` is vertical only.
- **Soil file written only once complete**: a member's soil file is rewritten when no value is missing, i.e. after the last zone or layer is filled (`_update_soil_params_from_perturbations`). A parameter updated for only some zones or layers leaves the rest at their previous values.
- **Order of operations per cycle**: analysis, back-transform, **clip to bounds** (new), validity check, write member inputs. Say that bounds apply to the sampling *and* to the analysed values, and that clipping piles values on the bound (log line `clipped to`).
- **Bounds and transforms**: `per_bounds = { min = , max = }` per parameter; `"None"` means unbounded; `transf_type` (`None`, `log`, `log-ratio`, `hyperbolic`; applied in `cathy_DA` before and after the analysis) acts on the update space; a bounded prior uses a truncated normal.
- **Reproducibility**: seeding behaviour (see "Also noted" in section 8) must be documented as it is today, then fixed (one generator from `DAConfig.seed`).

### 14.3 Facts D5 must state

| Mode | Unit of analysis | Observations used | Parameters updated | Typical use | Caveats |
|---|---|---|---|---|---|
| none | whole domain | all | all in `listUpdateParm` | few parameters, small domain | spurious long-range correlations with a small ensemble |
| `nodes` | one surface node (all its depths) | observations at that node | **all** parameters in every analysis *(to verify the shapes of `analysis_param` for several parameters)* | point sensors, dense maps | one analysis per node: cost; per-node parameter updates are not combined |
| `zones` | one zone of the zone map | observations inside the zone | state and `list_parm2update[zone_number]` (positional) | soil units with their own parameter | needs the nodata fix (B19) and matching list order (B22) |
| `veg_map` | one vegetation class | observations inside the class | same as `zones`, with the class parameter (e.g. `ZROOT<k>`, to verify that the Feddes row order matches the class order) | root zone by vegetation | classes numbered 1..K (B20) |
| covariance localisation matrix | state x obs, element-wise product with the covariance | all, damped by distance | per matrix | ERT, EM with many data | radius convention: weights vanish beyond `2 L` (Gaspari-Cohn); coordinates are `[y, x, z]` |

Also document: only one domain localisation at a time (a list with more than one entry raises); `"EM"` repeats the surface nodes per coil; a zone with no observation is skipped silently *(add a log line)*; nodes no zone analysed keep their forecast state (`_assemble_localised_analysis` warns); surface observations are extruded to all depths, so for ERT prefer a covariance localisation matrix that depends on depth, and say so (B24).

### 14.4 Code changes that the documentation depends on (do before writing D6)

1. ✅ *(done: `pyCATHY/DA/scenario.py`, called by `perturbate()`; skip with `validate_scenario=False`)* A **scenario validator** run before the ensemble is built: unknown keys, list lengths equal to the number of layers, zones or classes, `pert_control_name` aligned with `per_name`, zone and class numbering contiguous from 1 (B20, B21), `listUpdateParm` entries that exist in the perturbation dictionary. Fail with one message per problem.
2. ✅ **Zone-parameter pairing by name instead of position** *(done)*: `zones: {1: "Ks0", 2: "Ks1"}` in the scenario, keeping the positional form as a deprecated fallback (B22).
3. ✅ *(done: `scenario.describe_ensemble(list_pert, listUpdateParm, soil_index, zone_raster, veg_raster)`; it uses the real `_assign_soil_value`, so it cannot disagree with the run; wired into `DA.describe_perturbation()` and logged once at the start of a run)* A `describe()` method printing the resolved configuration for a run: for each parameter, the control, the rows of the soil table it touches, the bounds, the prior sd and whether it is updated. It doubles as the quickest way to catch a layout mistake.
4. Rewrite the `perturbate()` docstring (B23) and give `localisation.py` module and function docstrings with the conventions of B24.
5. Log, once per cycle, which zones or nodes were analysed, how many observations each used, and which parameters each updated.

### 14.5 Worked examples to write and keep under test

- ✅ E1: homogeneous soil, `Ks` and `porosity` scalar (the baseline).
- ✅ E2: two layers, `Ks0` and `Ks1` perturbed and updated (`layers`).
- ✅ E3: two zones from a zone map, `Ks` by zone, with nodata in the map (regression for B19).
- ✅ E4: vegetation map with three classes: a soil parameter by `root_map` (`Ks0` to `Ks2`) **and** the Feddes parameter `ZROOT0` to `ZROOT2` (row position in the Feddes table), both updated per class with `veg_map` localisation. The pair makes the difference between the two mappings visible.
- E5: *(removed: no combined layout exists in the code, see B32)*
- ✅ E6: covariance localisation matrix: `test_localisation.py` (13 tests: Gaspari-Cohn properties, the matrix against the former loop implementation, depth localisation, and the effect on the analysis: weights of one reproduce the unlocalised analysis, weights of zero leave the states untouched while the parameters are still updated). Still to add: a figure of the weights against distance for the documentation.
Each example is a small `scenario.toml` plus a test that builds the ensemble **without running CATHY** (using `describe()` and the soil table) and asserts which rows changed. Run in CI without the CATHY executable.

### 14.6 Documentation tooling

- Sphinx or MyST pages under `doc/`, with the parameter catalogue (D3) and the scenario key table (D2) **generated** from code at build time.
- Docstrings in NumPy style on every public function in `perturbate.py`, `localisation.py` and the `DA` methods named above; enforce with `ruff` rule set `D` (start with `D1` missing docstrings on public API, then `D2xx`) once the module split of section 10 is done.
- One page of "common errors" mapping each validator message to its fix.
- Doctests for the pure helpers (`_valid_zone_values`, `_assign_soil_value`, `gaspari_cohn`).

### 14.6b What the validator checks today (so the docs can cite it)

Errors: missing `per_name`/`per_nom`/`per_mean`/`per_sigma`; any list not aligned with `per_name`; duplicate names; **unknown parameter name (with a case-insensitive suggestion: `Zroot` suggests `ZROOT`)**; unknown control, sampling, perturbation type or transform; a control without per-entry lists; per-layer, per-zone and per-class lists whose length differs from the model; `per_sigma < 0`; malformed bounds, `min >= max`, mean outside bounds; zone or vegetation classes not numbered 1..K when `zone` or `root_map` is used; `listUpdateParm` entries that are not perturbed or out of range. Warnings: bounds narrower than 2 sd, `"St. var."` listed twice. All problems are reported in one `ScenarioError`. Not checked yet: the order of `listUpdateParm` against the zones (B22), the number of Feddes rows for a scalar `ZROOT`, the layout of nested (`zone` + `layers`) lists.

Also checked, after the build: every name of `listUpdateParm` must be a key `perturbate()` really produced (B26).

### 14.6c Worked examples already under test (`pyCATHY/tests/test_scenario_examples.py`)

E1 global scalar (`Ks`, all rows); E2 two layers; E3 two zones (one-based zone numbers against zero-based index); E4 vegetation classes with an overlap warning, plus the Feddes row mapping of `ZROOT1`; error cases (a write that would fail in the run, `root_map` without maps, an update name that was not built). E6 (localisation) is in `test_localisation.py`; E5 has no basis in the code (B32).

### 14.7 Order of work

1. ~~B20 to B22 validator~~ (done, see 14.6b); ~~pairing by name (item 2)~~ done. 2. ~~`describe()`~~ (done as `describe_ensemble`; wire it into `DA`). 3. Examples E1 to E4 done, E5 and E6 open. 4. D1, D2, D3 pages. 5. D4, D5, D6 pages. 6. Docstring pass and doctests. Steps 1 to 3 are small and testable without CATHY; steps 4 to 6 are mostly writing once the behaviour is pinned down.


---

## 15. Refactoring progress tracker

✅ = done in the patch and unit-tested (not yet validated on the real run). ☐ = open.

### Correctness
- ✅ A1 gain of the Sakov branch (`_solve_obs_space`)
- ✅ A2 lower bound and clipping of analysed parameters
- ✅ A3 bounded `ic` prior
- ✅ B1 observation error unit (opt-in `data_err_unit`)
- ✅ B2 ERT mask, flag never delete
- ✅ B6 sorted `ens_valid`
- ✅ B7 stale project folder (in the run script)
- ✅ B9 chained assignment in the Feddes update
- ✅ B12 `update_ic` on the right object
- ✅ B13 observation perturbation (opt-in `perturb_obs`; ☐ default to be decided by the author)
- ✅ B14 wrong ensemble mean in the inflation filters
- ✅ B15 dead `performance.py` deleted
- ✅ B17 mutable defaults in `DA/` and `ERT/` (☐ 49 remain: `cathy_plots` 28, `cathy_tools` 17, `petro` 2, `meshtools` 2)
- ✅ B18 `nlayers` kwarg
- ✅ B19 zone localisation off by one with nodata
- ✅ B20 non-contiguous classes rejected by the validator
- ✅ B23 `perturbate()` docstring
- ✅ B5 rejection ledger (`_rejections.json`, members by stage in the error)
- ✅ B10 positive psi checked after the analysis (reported, not clipped)
- ✅ B21 nodata is not a zone in `perturbate()`
- ✅ B30 `layers` sampling through the standard call
- ✅ B34 time-alignment checker (`time_alignment.py`): skipped observations, wrong pairing, one-window offset, logged before the first window
- ✅ B8/B33 `read_mbeconv` reads a one-step file (the cause of member 20's error)
- ✅ B22 zone-parameter pairing by name
- ✅ B31 `per_bounds` honoured by the prior of `layers`, `zone` and `root_map` (and stored for the A2 clipping)
- ✅ B24 localisation conventions documented and pinned by tests
- ✅ `perturb_obs` available in the localised filter too (it was silently missing there), with a test through every registered adapter
- ✅ `build_localization_matrix` vectorised, `gaspari_cohn` rejects `L <= 0`
- ✅ B25 `per_bounds` parsing
- ✅ B26 update names checked against the built keys
- ✅ B27 overlapping writes reported
- ✅ Reseeding in `sampling_dist` (R4g re-applied)
- ✅ R7 `NotImplementedError` readers, ✅ invalid escapes (0 left), ✅ 4 of the `F821` names
- ☐ A4 update of an unobservable parameter (decide: do not update ZROOT from ERT alone, T1)
- ☐ B3 rank-1 ensemble (T12), ☐ B4 whether the data use the same time convention (T9, needs the real data; the pairing rules are now documented and checked, see B34), ☐ B11 Archie assumptions
- ☐ B16 unreferenced filters, ☐ B32 decide whether to implement combined controls
- ☐ B28 push the earlier local changes or redo the logger conversion of `cathy_tools`, `cathy_utils`, `meshtools`, `cathy_plots`
- ☐ B29 `sensitivity.py` and `plot_perf` (45 `F821`)
- ☐ `enkf_dual_analysis` semantics of `COV[6]`; ☐ forward-model side of the ERT mask

### Nomenclature (section 9)
- ✅ local variables of the six filters in `enkf.py` (`A`, `dD`, `S`, `MeasAvg`, `COV`, `augm_*`, `*_pert` ...), proven AST- and bit-identical, golden test added
- ☐ public names: the `Sakov` keyword and `DA_type` strings (with a deprecated alias), `ens_valid`, `obs_valid`
- ☐ file and module renames (two-step, see the case-insensitivity caution)
- ☐ `threshold_rejected` to `max_rejected_percent`, unified error names (`dataErr`, `data_err` ...)

### Modularity (section 10)
- ✅ new self-contained, pure-function module `DA/scenario.py` (validator, `describe_ensemble`)
- ✅ `test_module_split.py` guards the re-exports and forbids cycles
- ✅ small tested helpers: `_valid_zone_values`, `_set_feddes_value`, `_clip_updated_parameter`, `_solve_obs_space`, `quadrupole_quality_mask`
- ✅ first four blocks of the `cathy_DA.py` split moved out (analysis and config, parameter tables, runtime, performance tables), AST-identical and re-exported; ☐ remaining: observation selection, ensemble I/O, the 8 methods over 80 lines
- ☐ shared helpers in `enkf.py` and a `NamedTuple` result instead of positional indexing

### Linter, tests, performance (section 11)
- ✅ `ruff` configured in `pyproject.toml`
- ✅ `test_package_compiles.py` (syntax, invalid escapes as errors, no mutable defaults in `DA/` and `ERT/`)
- ✅ analysis solve picks the cheaper side (Woodbury or direct) automatically
- ✅ whitespace autofix as a **separate patch** (`whitespace_only.patch`: 672 fixes in 19 files, AST identical in every file, checked by `verify_ast_identical.py`; the 50 findings left are trailing spaces inside string literals, which a fix would change)
- ☐ `ruff format` (reflows code; do it as its own commit and verify with the same script), ☐ `F841`, `E402`, `N802`
- ☐ mark `needs_cathy` tests and run the rest in CI
- ✅ R stored as its diagonal (opt-in `diagonal=True`, 466.6 MB to 193 kB for 2415 data x 10 cycles; the filters still receive the dense matrix of the current cycle)
- ☐ profile one cycle, ☐ diagonal-aware filters (removes the transient `m x m` of one cycle)

### Documentation (section 14)
- ✅ scenario validator, ✅ `describe_ensemble`, ✅ worked examples E1 to E4 as tests (plus the layered sampling path end to end), ✅ `perturbate()` docstring
- ✅ `describe_ensemble` wired into `DA` (`describe_perturbation`, logged at the start of a run)
- ✅ pairing zones to parameters by name
- ☐ per-cycle log of analysed zones
- ✅ E6 localisation tests; E5 removed (B32)
- ✅ localisation docstrings
- ✅ pages D1 to D6 (`doc/content/DA_guides/`), tested by `test_docs_examples.py`
- ✅ B35 `alpha_VG` / `VGP` rejected up front; ✅ B36 Feddes and soil layouts measured (60 cases), validated and documented; ✅ B37 mask x localisation matrix; ✅ Gaspari-Cohn weights clamped to [0, 1]
- ☐ the full Sphinx site build (needs the gallery), ☐ doctests, ☐ a figure of the Gaspari-Cohn weights
