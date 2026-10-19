# REFACTORING.md — pyCATHY DA / ERT tools

Living plan for AI sessions (Claude). **Read this first, update the Status column and the Changelog at the end of every session.**
Scope so far: `cathy_DA.py`, `mapper.py`, `observations.py`, `petro_Archie.py`, `simulate_ERT.py`, `sensors_measures.py` (+ new `ERT_checks.py`). Not seen yet: `enkf.py`, `cathy_tools.py`, `localisation.py`, perturbation module.

## Rules
1. Never `pip install pyCATHY` (pip is not the latest version). Ask the user for missing files; do not guess; do not fetch GitHub unless the user re-allows it.
2. Minimal, backward-compatible edits (exact-match single-occurrence replace; assert count == 1). `py_compile` after every change.
3. Every fix gets an assertion-based test (`test_*.py`, synthetic data, no pygimli/pyCATHY needed — extract functions via `ast` when the module imports heavy deps).
4. Logging: `logger = logging.getLogger(__name__)`; no `print`/`rich` in library code. INFO = one line per step/cycle/member; DEBUG = shapes, stats, captured pyGIMLi output; WARNING = failed checks / recoverable problems.
5. Failures are tracked and surfaced: ERT check reports -> `<project>_ERT_checks.json`.
6. Original files are in `/mnt/user-data/uploads/` (read-only); work in `/home/claude/refactored/`.
7. Not testable here: anything that imports pygimli / pyCATHY (not installed). Say so in the changelog.

## Status legend
DONE · TODO · NEEDS-SOURCE (file not provided) · NEEDS-USER (decision/data)

## P0 — correctness bugs
| # | Item | Status |
|---|------|--------|
| R2 | `use_log_transformed_obs=True`: data/predictions were log10 but R stayed linear (data effectively ignored). New `_log10_cov` (R_log = J R J, J=diag(1/(ln10 (d+eps)))), option `DAConfig.log_transform_cov` (default True; False = old behaviour). ERT check receives the same R. Test: `_log10_cov` gives 0.1/ln10 for 10 % errors | DONE |
| R3 | `transform_parameters`: back-transform result was discarded in `_DA_analysis` (analysed params stored in transformed space, then transformed again next cycle). New `_back_transform_analysis_param` (per column). Fixed `log-ratio` back-transform (was `exp((p-A)/(B-p))`, now `(A+B e^y)/(1+e^y)`), `hyperbolic` (was `log` of a missing key; now arcsinh / sinh). Removed dead `np.vstack`. Test: `test_transforms.py` | DONE |
| R3b | `update_pert_parm_dict` read row `-1` for the first parameter when no "St. var." came first; now counts only non-"St. var." rows | DONE |
| R4 | Forward noise: mapper never passed `fwdNoiseLevel`/`fwdNoiseAbs` from the obs dict to `simulate_ERT` (always 5 %). Now passed; real value logged (old message had a literal `{??}`) | DONE |
| R4b | Unverified: is the pyGIMLi noise seed identical across members (identical noise realisations)? Check `ert.simulate(..., seed=)` in pygimli | NEEDS-USER (run with pygimli) |
| R4c | rhoa <= 0 (decided by user): WARNING with count/min (`simulate_ERT`, plus the ERT check report), values kept tiny (1e-3 Ohm.m); log transform clips to 0 before `+eps` in `_log_transform_obs`, `_log10_cov`, `ERT_checks.check_innovation` (no NaN / -inf), histogram uses the same logs. Tests: `test_log10_cov.py` | DONE |
| R4d | Duplicate members (cathy_21 == cathy_16) seen in the log; `check_ensemble_predictions` flags it, root cause unknown (parameter perturbation / copy of folders?) | NEEDS-SOURCE (perturbation module, `cathy_tools.py`) |
| R4e | Archie noise: was `normal(0, (pert*sigma)**2)` (variance used as scale: ~1e-6 relative noise for pert=0.05; DA branch applied the last node's draw to the whole row). Now relative std: `sigma += N(0,1) * pert * |sigma|`, vectorised, floored at 1 % of the noise-free sigma (`pert_sigma_Archie=0.05` = 5 % noise on conductivity). Changes results only when `pert_sigma_Archie != 0`. Tested in `test_refactor_equivalence.py` | DONE (suggested by Claude, accepted by default: revisit if the thesis eq. 4.4 intends otherwise) |
| R4f | `rec_err` is an absolute error in Ohm.m (user): the DA path (`make_data_cov`, `data_err**2`) already used the variance; the legacy `init_obs_cov_matrice` (stored in `dict_obs['data_cov']`, not used by the DA) used `rec_err` / `data_err` unsquared: now squared | DONE |
| R4g | `perturbate.sampling_dist` (`np.random.seed(1)` per call) and the van Genuchten perturbation (`np.random.seed(seed=1)`) gave all parameters the same standard-normal draws (rank-1 parameter ensemble, overrode `DAConfig.seed`). Reseeding removed (user: reproducibility not needed; a seed set once beforehand still reproduces). Test: `test_perturbate_seed.py`. Unused `seed=True` arg of `perturbate_parm` left for API compat. Other samplers (`sampling_dist_trunc`: truncnorm `.rvs`) use the global RNG already | DONE |
| R4h | `enkf.analyze_covariances` built the dense `(n_state x n_state)` matrix (≈14 GB for 41616 nodes) and `np.linalg.cond(COV)` (SVD) on every analysis, results unused; the 'collapse' warning tested the mean of perturbations (0 by construction, always fired). Now: means computed from the perturbations (no dense matrix), diagnostics only at DEBUG, collapse test on the spread. Outputs identical (`test_enkf.py`) | DONE |

## P1 — robustness / bugs
| # | Item | Status |
|---|------|--------|
| R5 | ERT checks: scheme, k, rhoa, fwd-vs-obs quadrupoles (off-by-one), obs, R, Archie inputs, ensemble, innovation; WARNING via logger; JSON report | DONE (`ERT_checks.py`, `test_ERT_checks.py`) |
| R6 | `simulate_ERT`: bare `except:` mesh fallbacks -> loop over `.msh/.bms/.vtk` with `except Exception` + DEBUG; stdout muting with `os.dup2` left stdout muted on exceptions -> `_captured_stdout` | DONE (tested: `test_refactor_equivalence.py`; pygimli path untested) |
| R6b | Other bare excepts: `mapper._add_2_ensemble_Hx` (-> `AttributeError`), `observations.make_data_cov` and `petro_Archie.Archie_rho_DA` (-> `Exception`, same behaviour) | DONE |
| R7 | `sensors_measures.read_tensiometers/read_TDR_prob` returned undefined variables (UnboundLocalError) -> explicit `NotImplementedError` | DONE (need sample files to implement) |
| R7b | `observations.py` invalid escape `"$\Omega$"` -> raw string | DONE |
| R8 | `petro_Archie`: `SW_2_ER0_DA` NameError when `DA_cnb`/`Ens_nbi` not in kwargs (now defaults: None / from folder name); duplicate `Ens_nbi` keyword TypeError if present in kwargs; dead `DA_cnb = time_ass`. Still TODO: resipy branch, `Archie_ERT`, discharge branch | DONE (partial) |
| R8b | `create_ERT_survey_pg`: `res0` now required, cast to float array, shape-checked with a clear message (NaN/<=0 checks already existed). Explicit keyword args for `verbose/fwdNoiseLevel/fwdNoiseAbs` still TODO | DONE (res0) / TODO (kwargs) |
| R8d | Bare `except:` in the whole package (19 more in cathy_tools 9, cathy_plots 4, cathy_utils, sensors_measures, cathy_outputs, meshtools, localisation, perturbate) -> `except Exception:` (same behaviour, but Ctrl-C / SystemExit no longer swallowed). Next: narrow each to the real exception type | DONE (mechanical) / TODO (narrow) |
| R8c | Mutable defaults (62 remaining in the package: cathy_plots 28, cathy_tools 17, petro_Archie 7, perturbate 5...; harmless unless mutated - audit with `ast`, convert only the mutated ones to `None`):  `fwd_ERT_survey dump`, `prepare_observations dict_obs` -> None. `Archie_rho_DA` list defaults left (never mutated, positional API) | DONE (partial) |

## P2 — logging (print -> logger)
| # | Item | Status |
|---|------|--------|
| R1 | cathy_DA (≈35 console/print), mapper, observations, petro_Archie, simulate_ERT; package logger gets timestamped stdout handler if none configured; `_log_once` | DONE |
| R1b | `self.console` kept only to call `mapper.swc_mapper(console=...)` (API compat); remove when callers are updated | TODO |
| R1c | `enkf.py` (51 prints) and `cathy_tools.py` (90 print/console calls) -> logger via AST rewriter (levels by content: errors/warnings/INFO steps/DEBUG dumps; rich markup/emoji stripped; `self.console` kept for `rich_display` only); `_configure_package_logging(verbose)` in CATHY init keeps the old default (verbose=True shows INFO on stdout) unless the app configured logging. Tests: `test_enkf.py`, `test_cathy_tools_logging.py`. `perturbate.py` (2), `localisation.py` (2), `pf.py` (13), `cathy_utils.py` (10) also converted (AST rewriter). `performance.py` has no prints; also converted: `meshtools` (38), `importers/cathy_outputs` (26), `cathy_inputs` (2), `plotters/cathy_plots` (12), `cathyDA_plots` (3), `petro` (4), `sensitivity` (4; the Sobol/Morris result table stays at INFO). All via `print2logger.py` (AST rewriter, see Tools). Not scanned: `config_manager`, `ERT/` other files, `tests/` | DONE (whole package except unscanned) |

## P3 — docs / syntax / hygiene
| # | Item | Status |
|---|------|--------|
| R9 | `cathy_DA._nonphysical_values`: substring name checks -> `_has_name_token` (token at start or after a non-alphanumeric: keeps zone/layer suffixes, no longer matches inside words, e.g. `a_Archie` used to match `pert_sigma_Archie`). Test `test_name_matching.py`. Obs-key checks (`"EM"`, `"ERT"`, `"swc"`, `"all"` in `obskey2map`/`sensor` names) left as is (keys are `ERT_*`, `EM_*`...: confirm naming). `int ** array` and boxplot labels: not located | DONE (params) / TODO (obs keys) |
| R9b | (`SW_2_ER0_DA` done) NumPy-style docstrings with types for public functions (`SW_2_ERT_ERa_DA`, `Archie_rho_DA`, `create_ERT_survey_pg`, mapper functions); type hints | TODO |
| R9c | Remove dead/commented code: `ER0_to_EM_ECa_DA` block removed from petro_Archie; others in petro_Archie/simulate_ERT remain | DONE (partial) |
| R9d | `CATHY.__init__` ignored its `verbose` argument (`self.verbose = False`): now `self.verbose = verbose` (only gates a DEBUG dump of parm keys). `make_console`/`self.console` kept for `rich_display` and `swc_mapper(console=)` | DONE (verbose) / TODO (console removal) |

## P0 — crash-report round (session 3, from `DA_ERT_ensemble_crash_report.md`)
Tests: pytest under `pyCATHY/tests/` (numpy only, no CATHY run; `python -m pytest pyCATHY/tests/test_enkf_analysis.py ...`). Not executed on the real `ic_ZROOT_upd_ZROOT` case.
| # | Item | Status |
|---|------|--------|
| C1 | **A1** `enkf` Sakov branch used `diag(R)` only (dropped `S S^T/(N-1)`): gain inflated with the number of data (unobservable-parameter update 0.09 / 0.67 / 8.9 for 10 / 200 / 2000 data). New `_solve_obs_space` (exact, ensemble space, picks Woodbury or direct solve by size) used at the 4 sites. Test: `test_enkf_analysis.py` | DONE |
| C2 | **A2** lower bound for ZROOT in `_nonphysical_values`; analysed parameters clipped to `clip_min/clip_max` (`_clip_updated_parameter`); **B6** `ens_valid` sorted | DONE |
| C3 | **A3** `ic` prior honours `per_bounds` (truncated sampling); bounds set in `scenario.toml` | DONE |
| C4 | **B1** `data_err_unit='percent'` (`obs_error_variance`), default unchanged | DONE |
| C5 | **B2** ERT quadrupole mask: flag, never delete (`quadrupole_quality_mask`, `DAConfig.mask_bad_ert`, record in `_ERT_checks.json`) | DONE |
| C6 | **B13** observations not perturbed: analysis spread short by `K R K^T` (0.095 vs exact 0.65). Opt-in `AnalysisConfig.perturb_obs` | DONE (opt-in; author decides default) |
| C7 | **B14** `enkf_analysis_inflation*`: `reshape` of the transposed prediction scrambled the ensemble mean (diff vs plain EnKF 0.24 -> 3e-15) | DONE |
| C8 | **B9** chained `.iloc[...] =` in `_update_feddes_member` silently does nothing with pandas >= 3 -> `_set_feddes_value` | DONE |
| C9 | **B18** `perturbate(nlayers=)` assigned `nzones`; **B19** zone localisation off by one when the map has nodata (`_valid_zone_values`) | DONE |
| C10 | **B25** `check4bounds` swallowed errors and `process_perturbation_block` read bounds differently: one parser (`parse_bounds`), malformed bounds raise | DONE |
| C11 | **B20-B22, B26** `pyCATHY/DA/scenario.py`: `validate_scenario` (before build), `check_update_parameters` (after build), `describe_ensemble` (what each parameter touches; reports overlapping writes, B27) | DONE (pairing zones<->parameters by name still TODO) |
| C12 | **B15** dead `DA/performance.py` deleted; **B17** mutable defaults in `perturbate.py` and `Archie_rho_DA` (outputs identical, checked against git HEAD); **B23** `perturbate()` docstring rewritten | DONE |
| C13 | `F821` undefined names: `UNIT_DIVISOR`, `r2_score` (own `_r2_score`), `df_tensio`, `df_TDR` fixed (50 -> 45) | DONE for these; **TODO/NEEDS-USER**: `sensitivity.py` (37, half-ported script: rewrite with explicit arguments or remove), `cathyDA_plots.plot_perf` (8, relies on script globals) |
| C14 | `ruff` configured in `pyproject.toml` (baseline 1573 findings); `test_package_compiles.py` compiles every module with warnings as errors and rejects mutable defaults in `DA/`, `ERT/` | DONE |
| C15 | **B5** rejection ledger: `_record_rejection`, `rejection_summary`, `<project>_rejections.json`; the "too many rejected" error lists the members by stage (message prefix unchanged). Fed by `_check_before_analysis` and `_check_after_analysis` (tested with the negative-ZROOT0 pattern of the crash report) | DONE |
| C16 | **B10** positive psi is also checked on the **analysed** ensemble (members, max, share of nodes); **B4** each analysis logs the observation time and the model window (alignment itself not verified) | DONE (B4: log only) |
| C17 | **B21** `_n_classes` (nodata is not a zone) at the 7 counting sites of `perturbate()` | DONE |
| C18 | **B30** `layers` control could not be sampled through the standard `perturbate_parm` call (`nlayers` None -> TypeError; same on untouched HEAD): entries now carry `nlayers`, 5 branches raise a clear error without it | DONE |
| C19 | `DA.describe_perturbation()` (rows written by each parameter, overlaps) logged once at the start of a run (diagnostics, never stops a run) | DONE |
| C20 | **B31** `per_bounds` did not reach the prior of `layers` (Sobol sampler had no bounds) and was dropped for `zone`/`root_map` (entry stored `(None, None)`, so the A2 clipping had nothing either): `_apply_scenario_bounds` post-pass, truncated quantile in the layered sampler, max-only bounds recognised, `sampling_dist_trunc` accepts open sides and `sd = 0`, per-entry bounds resolved in the loops | DONE |
| C21 | **B24** localisation conventions in docstrings (`gaspari_cohn` support is 2 L, `L > 0` enforced; matrix shape and coordinates), `build_localization_matrix` vectorised with `cdist` (identical output, ~2x faster); `perturb_obs` also in `enkf_analysis_localized_with_inflation` (was silently ignored); tests `test_localisation.py`, `test_bounds_sampling.py` | DONE |
| C22 | **B32** the combined `zone`+`root_map`+`layers` layout described in a code comment is not implemented (one control string per parameter) | TODO (decision: implement or delete the comment) |
| C23 | **B22** zone-parameter pairing by name (`scenario.select_zone_parameters`: `ZROOT0` is zone 0 whatever its position; a zone without a parameter updates the states only instead of `IndexError`; positional rule kept for names without a numeric suffix) | DONE |
| C24 | **B8** `read_mbeconv` failed on a file with one single time step (1-D `loadtxt` result -> "Shape of passed values is (19, 1)...": the cause of member 20's error in the crash log, a run stopped after step 1); empty file raises "no time step" | DONE |
| C25 | whitespace-only cleanup (W291/W293/W292, safe fixes) delivered as a separate patch, AST identical in every file (`verify_ast_identical.py`); 50 findings inside string literals left on purpose | DONE (separate patch) |
| C26 | **Split of `cathy_DA.py` (3993 -> 3384 lines), step 1**: `analysis.py` (config + filter registry), `param_tables.py`, `runtime.py`, `performance_tables.py`; pure move proven by AST equality of the 35 definitions, identity of every re-exported name, no undefined names, tests; `scenario.py` now imports `param_tables` | DONE (obs selection, ensemble I/O and the long `DA` methods TODO) |
| C27 | **Nomenclature, `enkf.py`**: local variables of the 6 filters renamed to one vocabulary (documented in the module): 267 tokens, no argument/keyword touched; AST equal under the rename map, 17 filter/option cases bit-identical, `tests/test_enkf_golden.py` freezes the original outputs (mutation-tested) | DONE (public names such as the `Sakov` keyword TODO) |
| C28 | **B4/B34** pairing of observations and model states read from the code: window `[t, t_next]`, analysis if `t in ENS_times` (exact), data of the k-th analysis chosen **by counter**, state read with `iloc[-1]` (= `t_next`). New pure `DA/time_alignment.py` (`check_time_alignment`) logged before the first window: skipped observations (final time, rounded floats) and the wrong pairing they cause, unsorted/duplicated times, the one-window offset. Log only; behaviour unchanged | DONE (detection; whether the data use the same convention needs the real data) |
| C29 | **Performance**: `make_data_cov(..., diagonal=True)` stores R as its diagonal (466.6 MB -> 193 kB for 2415 data x 10 cycles); `DA._data_cov_of_cycle` rebuilds the dense matrix of the current cycle, so filters/diagnostics get the identical matrix (bit-identical, tested). Default dense. Run script passes `diagonal=True` and its `'$\Ohm$'` escape fixed | DONE (opt-in; diagonal-aware filters TODO) |
| C30 | **Documentation D1-D6**: `doc/content/DA_guides/{index,perturbation,localisation,troubleshooting}.rst`, linked from `doc/index.rst`; builds with Sphinx `-W`; every `toml` example run through `perturbate` by `tests/test_docs_examples.py` | DONE (full site build needs the gallery) |
| C31 | **B35/B36** 60 scenario layouts measured against the real `perturbate`+sampler, validator checked to agree both ways (`tests/test_scenario_layouts.py`, mutation-tested); `alpha_VG` (cannot be written to the members) and `VGP` (commented out) rejected up front; Feddes rules (ZROOT/PCREF one entry per row, PCWLT/OMGC/PZ scalar only, per_type scalar); `perturbate_dist` clear error for unknown `per_type` | DONE |
| C32 | **B37** the ERT quadrupole mask broke the covariance-localisation matrix (columns = observations): `_mask_localisation_matrix`; `gaspari_cohn` weights clamped to [0, 1] | DONE |

### Reality check: the pushed commit `43f7126` vs this log (found in session 3)
Items marked DONE above in sessions 1-2 that are **absent from the code on GitHub** (checked with markers and counts):
- R4g (reseeding) and R7 (`NotImplementedError` in `sensors_measures`): **re-applied in session 3** (`test_perturbate_seed.py`).
- Invalid `\d` / `\psi` escapes in `perturbate.py` and `cathy_plots.py`: re-fixed (0 left in the package).
- R1c print/console -> logger: present for `enkf`, `pf`, `localisation`, `cathy_outputs`, `cathy_inputs`, `cathyDA_plots`, `petro`, `sensitivity`; **absent** for `cathy_tools` (46 `print` + 49 `console`), `cathy_utils` (10), `meshtools` (38), `cathy_plots` (12). `perturbate.py` had 2 prints (converted in session 3). The `print2logger.py` tool is not in the repo.
- R4h was checked and **is** present (the `cond` call is guarded by `isEnabledFor(DEBUG)`).
**NEEDS-USER:** push the session-2 changes that were made locally, or ask for the print -> logger conversion of the four modules to be redone.

## How to resume
1. Re-run: `python -W error::RuntimeWarning test_ERT_checks.py; python test_refactor_equivalence.py; python test_transforms.py; for f in *.py; do python -m py_compile $f; done`
2. Pick the first non-DONE row in priority order, state which in one line, implement, add a test, update the table + Changelog.
3. Run-dependent checks (R4b, R4c, R4e): look at `<project>_ERT_checks.json` and the INFO log of a real run before deciding.

## Changelog
- Session 1: ERT_checks module + integration (simulate_ERT, petro_Archie, mapper, observations, cathy_DA); logging migration; R2, R3, R3b, R4, R6, R6b, R7, R7b. Tests: standalone only; pygimli/pyCATHY integration NOT executed.
- Session 1 (cont.): R8b (res0), R8c (partial). Next: R8, R9, R9b, R9c, then items needing user input.
- Session 1 (cont. 2): R8 partial (petro_Archie kwargs/NameError).

- Session 1 (cont. 3): R9 reviewed, blocked on naming list.

- Session 2: fetched enkf.py and cathy_tools.py (GitHub, read-only clone, repo cathy_DA/mapper identical to uploads); R1c done for them, R4h fixed, R4g found (perturbate.py seeds), R9d noted.
- Session 2 (cont.): print->logger also in perturbate, localisation, pf, cathy_utils. Pre-existing SyntaxWarning (invalid `\d` escapes) in perturbate.py lines ~1382/1417/1514: TODO raw strings.
- Session 2 (cont. 2): R4g fixed (no reseeding in perturbate.py).
- Session 2 (cont. 3): print->logger for meshtools, cathy_outputs, cathy_inputs, cathy_plots, cathyDA_plots, petro, sensitivity. Invalid-escape SyntaxWarnings fixed with raw strings (perturbate, cathy_plots, cathyDA_plots, observations); `test_cathy_tools_logging.py` now fails on any print / rich markup / invalid escape in the converted modules.

## Tools (reusable)
- `print2logger.py FILE [--apply]`: AST rewriter `print(...)`/`self.console.print/rule(...)` -> `logger.<level>(...)` (level by message content: error/warning keywords, INFO steps, else DEBUG; rich markup/emoji removed; dry-run lists every change). Afterwards add `import logging` + `logger = logging.getLogger(__name__)` and review the levels. (Kept in the session scratchpad; ask to regenerate if needed.)
- `rawfix.py FILE...`: prefixes with `r` the string literals whose backslash sequences are all invalid escapes (same value, no SyntaxWarning).
- Session 2 (cont. 4): R9 params done, R9d verbose fixed.
- Session 2 (cont. 5): R8d bare excepts -> except Exception (19).

- Session 2 (cont. 6): R4c, R4e, R4f decided/implemented.
- Session 3: crash-report round (C1-C14 above): Sakov gain, bounds/clipping, perturb_obs (opt-in), ERT mask (flag only), scenario validator + describe_ensemble, B9/B14/B18/B19/B25-27 fixed, F821 partly. Found that the pushed code lacks several session-1/2 items (see Reality check); R4g, R7 and the escape fixes re-applied. 102+ synthetic tests; nothing run on the real case.
- Session 3 (cont.): C15-C19 (rejection ledger, analysed-psi check, assimilation-time log, nodata not a zone, layered sampling fix, describe_perturbation). 115 synthetic tests; nothing run on the real case.
- Session 3 (cont. 2): C20-C22 (bounds reach every prior, localisation tests and docstrings, perturb_obs in the localised filter). 140 synthetic tests; nothing run on the real case.
- Session 3 (cont. 3): C23-C25 (zone pairing by name, one-step mbeconv, whitespace patch with AST proof). 150 synthetic tests; nothing run on the real case.
- Session 3 (cont. 4): C26 (first four modules split out of cathy_DA.py, AST-proven). 163 synthetic tests; nothing run on the real case.
- Session 3 (cont. 5): C27 (enkf.py vocabulary, proven identical, golden test). 172 synthetic tests; nothing run on the real case.
- Session 3 (cont. 6): C28 (observation/state time pairing documented and checked before the run). 181 synthetic tests; nothing run on the real case.
- Session 3 (cont. 7): C29 (diagonal storage of R). 186 synthetic tests; nothing run on the real case.
- Session 3 (cont. 8): C30-C32 (guides, measured layouts, mask x localisation). 255 synthetic tests; nothing run on the real case.
