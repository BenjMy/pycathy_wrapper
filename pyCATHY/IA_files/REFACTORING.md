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
