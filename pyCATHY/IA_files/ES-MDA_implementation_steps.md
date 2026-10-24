# ES-MDA implementation in `pycathy_wrapper` — step-by-step plan

> **How to use this file.** Paste it (or one step at a time) into a new AI request, e.g.
> *"Follow `ES-MDA_implementation_steps.md`, do **Step 3 only**, then stop and show me the diff and tests."*
> Each step is self-contained, lists the files to touch, and has an acceptance check.
> Do not start a step before the previous one passes its check.

Repo: `https://github.com/BenjMy/pycathy_wrapper` (cloned locally as `pycathy-wrapper/`).

---

## 0. Background

**ES-MDA** (Emerick & Reynolds, 2013) = Ensemble Smoother with Multiple Data Assimilation.
Instead of one Kalman update with all data, the *same* data are assimilated `Na` times, each with an
inflated observation-error covariance `alpha_i * R`, and the **full forward model is re-run from t = 0
between iterations**.

Constraint on inflation coefficients (so the posterior equals a single Kalman update in the linear-Gaussian case):

```
sum_{i=1..Na} 1 / alpha_i = 1
```

Simple default: `alpha_i = Na` for all i. Geometrically decreasing alternatives are also common.

Update for member `j` at iteration `i`:

```
d_obs,j = d_obs + sqrt(alpha_i) * R^(1/2) * z_j ,   z_j ~ N(0, I)       # fresh draw every iteration
m_j^{i+1} = m_j^i + C_MD (C_DD + alpha_i R)^(-1) (d_obs,j - g(m_j^i))
```

where `C_MD` = cross-covariance (parameters, predicted data), `C_DD` = covariance of predicted data,
`g(.)` = forward model (CATHY) run over the **whole** assimilation window.

---

## 1. What already exists in the repo (read before coding)

Explored in `pyCATHY/DA/`:

| File | Relevant content |
|---|---|
| `analysis.py` | `AnalysisConfig`, `DAConfig`, `AnalysisResult`, registry `ANALYSES` (dict `DA_type -> adapter`), `run_analysis(...)` |
| `enkf.py` | `enkf_analysis(data, data_cov, param, ensemble, predict_obs, **kw)` (supports `Sakov`, `perturb_obs`, `alpha`), localized/inflation/dual variants, helper `_solve_obs_space`, `_sample_obs_noise` |
| `cathy_DA.py` | `prepare_DA`, `map_states2Observations`, `_read_state_ensemble`, `_DA_analysis`, `_update_input_ensemble`, `run_DA_smooth` (**stub: runs one window and ends with `pass`**), `run_DA_sequential` / `_run_window` (sequential EnKF loop) |
| `perturbate.py` | parameter ensemble perturbation |
| `localisation.py` | Gaspari–Cohn localisation matrices |
| `tests/` | `test_enkf.py`, `test_enkf_analysis.py`, `test_enkf_golden.py`, `test_module_split.py`, ... |
| `examples/DA/` | `plot_3_run_sequentialDA_SMC.py` etc. — templates for a new example |

**Key observation:** `run_DA_smooth` already has the "run all DA times stacked, then analyse once" skeleton
(ES with one iteration). ES-MDA = loop that skeleton `Na` times with inflated `R` and a parameter-only
(or parameter + initial-state) update followed by a model re-run. Re-use, do not rewrite.

---

## 2. Design decisions (confirm with the user before Step 1)

- [ ] **What is updated?** Default: *parameters only* (+ optionally initial conditions). State is **not** updated
      mid-run, since the model is restarted from t = 0 each iteration.
- [ ] **Number of iterations `Na`** and **alpha schedule** (`constant` | `geometric` | user list).
- [ ] **Observation operator** reuse: `map_states2Observations` over *all* assimilation times at once.
- [ ] **Perturbed observations**: fresh draw per iteration and per member (required by ES-MDA).
- [ ] **Parameter transforms / bounds**: reuse existing log10 / bounds handling (see tests `test_log10_cov.py`,
      `test_bounds_sampling.py`, `test_reflect_bounds.py`).
- [ ] **Localisation / inflation**: optional, off by default.
- [ ] **Failed members**: reuse the existing rejection logic (`threshold_rejected`, `ens_valid`).

---

## 3. Implementation steps

### Step 1 — Inflation coefficient schedule
- **File:** `pyCATHY/DA/esmda.py` (new)
- Implement `compute_alphas(n_iter, schedule="constant", user_alphas=None) -> np.ndarray`.
  - `constant`: `alpha_i = n_iter`.
  - `geometric`: decreasing ratio (e.g. factor 0.5–0.8), then **normalise** so `sum(1/alpha) == 1`.
  - `user`: validate and normalise/raise if the constraint is violated (tolerance 1e-8).
- **Check:** unit test that `np.isclose(np.sum(1/alphas), 1)` for every schedule and several `n_iter`.

### Step 2 — Pure-numpy ES-MDA update (no CATHY dependency)
- **File:** `pyCATHY/DA/esmda.py`
- Implement
  ```python
  def esmda_update(param, predicted, data, data_cov, alpha, rng=None, localisation=None):
      """param: (n_param, n_ens); predicted: (n_obs, n_ens); data: (n_obs,);
      data_cov: (n_obs,) variances or (n_obs, n_obs). Returns updated param (n_param, n_ens)."""
  ```
  1. Perturb data: `d_j = data + sqrt(alpha) * L z_j` (use/extend `enkf._sample_obs_noise` or `data_cov` Cholesky; support diagonal).
  2. Anomalies: `dM = param - mean`, `dD = predicted - mean`.
  3. `C_MD = dM dD^T / (N-1)`, `C_DD = dD dD^T / (N-1)`.
  4. Solve `(C_DD + alpha R) X = (d_j - predicted)` with `np.linalg.solve` / `lstsq` (**never** explicit inverse; consider ensemble-space/SVD form from `enkf._solve_obs_space` when `n_obs >> n_ens`).
  5. `param_new = param + C_MD X`.
- **Check:** (a) linear-Gaussian toy problem (`g(m) = G m`): after `Na` iterations with `sum(1/alpha)=1`,
  ensemble mean/cov ≈ analytical Kalman posterior (large N); (b) shapes; (c) reproducibility with a seeded `rng`.

### Step 3 — Configuration object
- **Files:** `pyCATHY/DA/analysis.py` (or `esmda.py`)
- Add `ESMDAConfig` dataclass: `n_iter`, `alpha_schedule`, `user_alphas`, `seed`, `update_initial_conditions: bool`,
  `inflation_params`, `localisation_matrix`, `cleanup_between_iterations: bool`.
- Follow the style of `AnalysisConfig` / `DAConfig` (`__post_init__` validation, `from_legacy_kwargs` only if needed).
- **Check:** invalid configs raise `ValueError`.

### Step 4 — Register as an analysis type
- **File:** `pyCATHY/DA/analysis.py`
- Add an adapter `_esmda(data, data_cov, param, state, prediction, config)` returning `AnalysisResult`
  (`analysis` = unchanged state, `analysis_param` shaped `(n_members, n_param)` like the other EnKF filters, `extra={"alpha": ...}`).
- Add `"esmda": _esmda` to `ANALYSES`.
- **Check:** `run_analysis("esmda", ...)` works on synthetic arrays; existing tests (`test_module_split.py`, `test_refactor_equivalence.py`) still pass.

### Step 5 — Smoother driver with iteration loop
- **File:** `pyCATHY/DA/cathy_DA.py` (new method next to `run_DA_smooth`)
- Add `run_DA_esmda(self, dict_obs, list_parm2update, dict_parm_pert, list_assimilated_obs, config: ESMDAConfig, parallel=True, **kwargs)`:
  ```
  alphas = compute_alphas(...)
  prepare_DA(...)                                   # once
  for i, alpha in enumerate(alphas):
      _update_input_ensemble(list_parm2update)      # write current params to each ensemble folder
      run forward model over ALL assimilation times  # reuse the code path of run_DA_smooth
      prediction = map_states2Observations(...)      # all times stacked
      validate prediction / reject failed members    # _validate_prediction, ens_valid
      param = esmda_update(param, prediction, data, R, alpha, rng)
      apply bounds / transforms
      store iteration results (params, misfit, alpha)
  final: one more forward run with the last params (to get final predictions/states)
  ```
- Reuse the stacked-time logic from `run_DA_smooth`; remove its trailing `pass` only if safe (keep backward compatibility).
- Make sure ensemble folders are re-initialised from the **same initial conditions** every iteration.
- **Check:** dry run on the smallest test case in `pyCATHY/tests/` with `Na = 2`, 4–6 members; no crash, params change between iterations.

### Step 6 — Data handling for stacked times
- **Files:** `cathy_DA.py` (`_get_data2assimilate`, `_obs_key_select`), `time_alignment.py`
- Build the **stacked** `data` vector and `data_cov` for all observation times and observation types
  (SWC, ERT, ET, ...) in a consistent order matching `prediction`.
- Support heterogeneous error variances (block diagonal `R`).
- **Check:** `data.shape[0] == prediction.shape[0]` and ordering test using `test_time_alignment.py` as a model.

### Step 7 — Diagnostics and logging
- **Files:** `esmda.py`, `pyCATHY/DA/performance.py`, `plotters/cathyDA_plots.py`
- Per iteration, log/store: `alpha_i`, RMSE and weighted misfit `mean_j ||R^-1/2 (d - g(m_j))||^2`, ensemble spread of each parameter.
- Add plot: misfit vs iteration; parameter ensemble (prior vs each iteration vs posterior).
- Reuse `console.py` reporting if suitable.
- **Check:** misfit decreases (or stabilises) on a synthetic twin experiment.

### Step 8 — Tests
- **Folder:** `pyCATHY/tests/`
- `test_esmda_alphas.py` (Step 1), `test_esmda_update.py` (Step 2: linear-Gaussian vs analytical posterior, seed reproducibility, diagonal vs full `R`),
  `test_esmda_driver.py` (Step 5: mock forward model — **do not require the CATHY executable in CI**; monkeypatch the run + mapping functions).
- Run the whole suite to verify no regression: `pytest pyCATHY/tests -x -q`.

### Step 9 — Example script and docs
- **Folder:** `examples/DA/` → `plot_6_run_ESMDA.py` (copy structure of `plot_3_run_sequentialDA_SMC.py`).
- Add a docstring section on ES-MDA to the docs in `doc/` and a short entry in `README.md`.
- Include the alpha-sum constraint and the "parameters only, model restarted each iteration" caveat.

### Step 10 — Review checklist
- [ ] `sum(1/alpha) == 1` enforced everywhere
- [ ] Fresh observation perturbations for every iteration **and** every member
- [ ] Linear solves, no explicit matrix inverses
- [ ] Bounds / log-transform respected after each update
- [ ] Failed ensemble members handled (excluded from covariance, optionally resampled)
- [ ] RNG seeded through config for reproducibility
- [ ] Backward compatibility of `run_DA_smooth`, `run_DA_sequential`, `ANALYSES`
- [ ] Type hints + numpy-style docstrings (match repo)
- [ ] All existing tests + new tests pass

---

## Progress log

| Step | Status | Notes |
|---|---|---|
| 1 Alpha schedule | done | `compute_alphas` in `pyCATHY/DA/esmda.py` (constant / geometric / user) |
| 2 Numpy update | done | `esmda_update`; reuses `enkf._solve_obs_space`; validated against the analytical Kalman posterior |
| 3 Config | done | `ESMDAConfig` (alphas computed and validated in `__post_init__`) |
| 4 Registry | done | `"esmda"` in `ANALYSES`; `AnalysisConfig.esmda_alpha` / `esmda_rng`; state returned unchanged |
| 5a Generic driver `run_esmda` | done | `esmda.py`: forward model as a callable; drops failed members, drops non-finite observations, `project` hook for bounds, callback, history |
| 5b CATHY adapter `DA.run_DA_esmda` | done, **experimental** | `esmda_cathy.py` + thin method in `cathy_DA.py`; `tensio`, `swc` and `ERT`; tested against a stub, **not** against the CATHY executable |
| 6 Stacked data | done | `stack_data` (per-sensor values, ERT quadrupole mask, log10), `observation_plan`, `predict_member_blocks`, `predict_ert_blocks` (state read at every observation time) |
| 7 Diagnostics | done | `pyCATHY/plotters/esmda_plots.py`: `plot_misfit`, `plot_parameter_evolution`, `plot_esmda_diagnostics`; tested headless |
| 8 Tests | done (stub level) | `test_esmda.py` (43), `test_esmda_cathy.py` (83, stub simulation + real `_get_data2assimilate` comparison), `test_esmda_plots.py` (8); a real CATHY run is still needed |
| 9 Example / docs | done | `examples/DA/plot_6_ESMDA_synthetic.py` (runs without CATHY, safe for the gallery build) and `doc/content/DA_guides/esmda.rst` |

Branch status: rebased on upstream `main` (`526056a7`, "refactor pf"); one conflict, in `AnalysisConfig`
(upstream added `pf_temper_neff_frac` / `pf_liu_west_a`, this work adds `esmda_alpha` / `esmda_rng`): both kept.
Against upstream `main` alone the ES-MDA branch adds no failure (10 failed -> 9 failed, 15 collection errors unchanged, all
from the sandbox environment); passing tests 284 -> 416. The one test that now passes,
`test_docs_examples.py::test_the_guides_are_linked_from_the_documentation_index`, was failing upstream because
`doc/content/DA_guides/index.rst` did not exist; the index added here fixes it.

Hardening from reading the real `DA` code (the stub had hidden it): `grid3d` is loaded lazily by `DA` (not by
`prepare_DA`), so the adapter now does the same dry run as `selec_atmbc_window` before the first `swc` prediction;
observation entries are validated up front (`data`, `mesh_nodes`); `KeyError` is no longer treated as a member failure.

Still open: validate `DA.run_DA_esmda` on a real CATHY run (ERT included); EM / ET observation support;
`update_initial_conditions`.

Also done: `doc/content/DA_guides/index.rst` added (it was referenced by `doc/index.rst` but never existed), so the guide is reachable.

Baseline note: 10 tests fail and 15 modules error at collection on the untouched repo
(missing CATHY binaries, network, optional deps); the sets are identical with the ES-MDA changes.

### ERT and time convention (added after reading nb4/nb5/nc_model.py)

- ERT: the sequential operator is reused through a new optional `time_ass` argument of
  `mapper._map_ERT_parallel` (row of the member `sw` output; default None = unchanged behaviour, the only edit
  outside the ES-MDA modules). Quadrupole mask and log10 of data / predictions / R follow the sequential analysis.
- `ESMDAConfig.obs_time_convention`: `"state_time"` (default: `run_DA_sequential` since `d37036f9`, 2026-10-19,
  analysis when the window ENDS at the stamp) or `"window_start"` (up to `083f0ac`, which `nb5` / `obs_window_start`
  were written for). `pyCATHY.DA.time_alignment` still documents the older pairing.
- `list_parm2update` is matched by exact key or root of numbered keys (the old substring rule selected `a_Archie`
  when asking for `ic`).
- `DA._get_data2assimilate` cannot stack sensors of different lengths in one cycle (ragged `np.squeeze`): the smoother
  extracts per sensor instead.
- `nb6_ESMDA_ERT_soil3D_nc.py` = `nb5` with only the sections listed in its header changed (`nb5_to_nb6.diff`).
  Not run: it needs your data, pygimli and a CATHY executable.

- `run_DA_esmda(..., dry_run=True)` prepares and checks everything and returns the plan and cost without simulating;
  `nb6` exposes it as `DRY_RUN`. Run it first on any new setup.

### Found on the first real run of nb6 (time check, console log)

- The pre-run check failed on the first survey: `nb5`'s diagnostic `simuDA.read_inputs('atmbc')` **replaces**
  `simuDA.atmbc` by the float32 read-back of the file (`prepare_DA` copies `simuDA.atmbc['time']`), and the check
  used a fixed 1e-6 s tolerance. Fixed: `nb6` restores `simuDA.atmbc` after the diagnostic; the check is now
  `esmda_cathy.check_state_times` (relative tolerance, tested with a float32 vector, explains the likely cause).
  The same `read_inputs('atmbc')` call exists in `nb5`: its effect on the sequential run is not verified here.
- Console log: `run_DA_esmda` already logs like `run_DA_sequential` (`<project>/console_log.txt`, appended, tested);
  nb6 also logs the part before the run (`LOG_CONSOLE`, `<simulations>/console_log_nb6_ESMDA_ERT.txt`).

### Second real run of nb6: nothing left to estimate

- With `nc_ERT_ic_VGP_pert` the scenario's `listUpdateParm` is `['St. var.']` (written for the sequential filters, where the
  state update corrects the initial condition), so `resolve_parameter_names` correctly refused. `nb6` now adds `'ic'`
  (`ESMDA_ESTIMATE_IC`): with a `_pert` scenario `ic` is the only estimated parameter, the VGP keep their prior sample.
- `DA._update_ic_member` writes a global `ic` from `ini_perturbation` at cycle 0 (the key the smoother rewrites), but a
  LAYERED `ic` from `ini_ic_withLayers`: the adapter raises `NotImplementedError` for it instead of silently ignoring updates.

### Dry run of nb6 on real data worked; recap and progress logging added

- The real dry run reported 1 parameter (`ic`), 100 members, 72 718 data over 21 windows (= the kept quadrupoles of the
  filtered surveys), 5 model runs, 500 member simulations, 10 500 ERT forwards: about 5x the cost of the EnKF run
  (21 x 100 = 2 100 ERT forwards and one pass of the period). `NENS = 32, ESMDA_N_ITER = 2` costs about one EnKF run.
- Before the first model run (and in the dry run) the strategy, the atmbc period and a table of the assimilation windows are
  logged at INFO (console + `console_log.txt`); one INFO line per iteration gives valid members, misfit and alpha.
  `nb6` has `SAVE_PREDICTIONS = False` (the saved result is ~2 kB instead of ~290 MB).

### What the adapter assumes (check on the first real CATHY run)

- Observation times (`dict_obs` keys) are on the same absolute time axis as the atmbc times and are
  **printing times**: the adapter writes `TIMPRTi = [0, obs times - t0]` and `TMAX` to every member's `input/parm`
  (`update_parm(..., IPRT1=2, ...)`, copied from `selec_atmbc_window`). A missing time raises `ESMDATimeError`.
- Member folders still hold the **full** atmbc / initial condition after `prepare_DA`, and CATHY overwrites its
  `output/psi` and `output/sw` on every run.
- The iterate is written by overwriting `dict_parm_pert[name]["ini_perturbation"]` and calling
  `apply_ensemble_updates(cycle_nb=0)`; the prior is restored afterwards and the posterior stored as `update_nb1`.
- `_check_before_analysis` (mass-balance check) uses `self.parm[...TIMPRT...][-1]` as the expected end time: verify it
  reflects the window set by `update_parm`.
- Not supported yet (raise `NotImplementedError`): EM / ET observations, `update_initial_conditions`,
  `cleanup_between_iterations`.

---

## 4. Request template for future AI sessions

```
Context: repo pycathy-wrapper (pyCATHY/DA). Plan: ES-MDA_implementation_steps.md.
Task: implement Step <N> only.
Rules: reuse existing helpers (enkf._sample_obs_noise, _solve_obs_space, run_DA_smooth skeleton),
       keep backward compatibility, add tests, do not run CATHY in tests.
Output: unified diff + the acceptance check result. Stop after this step.
```

## 5. References

- Emerick, A. A., & Reynolds, A. C. (2013). *Ensemble smoother with multiple data assimilation.* Computers & Geosciences, 55, 3–15.
- Evensen, G. (2009). *Data Assimilation: The Ensemble Kalman Filter* (2nd ed.). Springer.
