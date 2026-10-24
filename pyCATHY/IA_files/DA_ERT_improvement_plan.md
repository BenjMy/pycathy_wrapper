# DA + ERT improvement plan: Agripolis 3D NC (`nb5_DA_ERT_soil3D_nc.py`)

Status: **PLAN ONLY, nothing fixed yet** (as requested). Written 2026-10-20.
Sandbox: `/home/claude/sandbox/pycathy_wrapper` (clone of `https://github.com/BenjMy/pycathy_wrapper`, `main`, HEAD `4cb4d072` "Update cathy_tools.py").
This file is the single source of truth for restarting the AI request: every item has a checkbox and a **Check** (how to prove it is done). Update the boxes and the *Progress ledger* at the end of each session.

Legend of evidence tags: **[VERIFIED]** read in the code of the clone · **[MEASURED]** run in the sandbox · **[SUSPECTED]** needs a test before fixing · **[IA]** comes from `pyCATHY/IA_files/` (earlier AI sessions).

---

## 0. How to restart (copy this prompt into a new session)

> Read `DA_ERT_improvement_plan.md`. Clone `https://github.com/BenjMy/pycathy_wrapper` into a sandbox first (rule 1). Take the first unchecked item of section 6 (order of work), state which one in one line, write the failing test first, fix, tick the box, add a line to the ledger (section 10). Files I must provide if missing: see section 1.

### Rules (merged from `IA_files/REFACTORING.md` + your instructions)
1. **First action of every session: clone the repo into a sandbox** (you re-allowed GitHub; REFACTORING.md rule 1 said not to fetch it unless re-allowed). Never `pip install pyCATHY` (not the latest).
2. Ask for missing files, do not guess.
3. Minimal, backward-compatible edits; `py_compile` after each; **one assertion-based test per fix** (`pyCATHY/tests/test_*.py`, synthetic data).
4. Library code logs with `logging` (no `print`). Failures end up in `<project>_ERT_checks.json` / `<project>_rejections.json`.
5. Not testable in the sandbox: anything importing `pygimli` or running CATHY. Mark such items `NEEDS-RUN` and give the exact command for your laptop.
6. Tick the checkbox only when the **Check** passes. Items that need your real run stay `[ ]` + `NEEDS-RUN`.

---

## 1. Inputs still needed (ask before fixing the related items)

- [ ] `ert_alignment.py` (imported by nb5: `suggest_frame`, `transform_points`, `transform_tldat_electrodes`, `check_ert_cathy_alignment`, `bounds_from_nodes`, `ERTAlignmentError`). Not in the upload nor in the repo. **Check:** file received.
- [ ] Log of the last nb5 run (scenario `nc_ERT_ic_VGP_pert`, all members rejected at cycle 5) + `<project>_ERT_checks.json`, `<project>_rejections.json`.
- [ ] `nb5_plots_suggestions.md` (cited in the header of `nb5_DA_ERT_soil3D_nc_plots.py`).
- [ ] Real acquisition times of the 6 ERT surveys (nb5 `ERT_TIMES_H` holds placeholders: `{0:-1, 1:3, 2:9, 3:21, 4:45, 5:69}`).
- [ ] Output of `describe_ensemble` (the table printed in nb5 section 5) for the scenario you run.
- Already in the repo, no need to upload: `ERT/simulate_ERT.py`, `ERT/ERT_checks.py`, `DA/*.py`, `IA_files/*`.
- Note: `nb5_DA_ERT_soil3D_nc_1.py` is byte-identical to `nb5_DA_ERT_soil3D_nc.py` (checked with `diff`): ignore the copy.

---

## 2. State of the repo (what the AI sessions already did, checked in the clone)

`IA_files/REFACTORING.md` + `DA_ERT_ensemble_crash_report.md` list fixes C1–C32 / B1–B37. I grepped the clone for the marker function of each: **the session-3 patch is merged in HEAD** (`_solve_obs_space`, `AnalysisConfig.perturb_obs`, `parse_bounds`, `_apply_scenario_bounds`, `validate_scenario`, `check_update_parameters`, `select_zone_parameters`, `_record_rejection`, `check_time_alignment`, `_n_classes`, `_valid_zone_values`, `_set_feddes_value`, `quadrupole_quality_mask`, `obs_error_variance`, `_mask_localisation_matrix`, `_clip_updated_parameter`).
Baseline **[MEASURED]**: `pytest` on `test_scenario_layouts, test_scenario_validation, test_bounds_sampling, test_enkf_analysis, test_time_alignment, test_scenario_examples` → **119 passed** (numpy/scipy/pandas + `shapely pyvista rich pyyaml xarray ipywidgets natsort geopandas rioxarray gitpython` needed only to import the package).

Still open from the IA files (merged into the lists below): B16, B27 (decision), B28 (prints→logger in `cathy_tools`, `cathy_utils`, `meshtools`, `cathy_plots`), B29/C13 (`sensitivity.py`, `plot_perf`), B32/C22, C6/B13 (`perturb_obs` default), C29 (diagonal-aware filters), R1b, R9b, R8b kwargs, R4b (pyGIMLi noise seed).

Repo hygiene found now:
- [x] **H1 [VERIFIED] DONE s2** `doc/index.rst` toctree lists `content/DA_guides/index.rst` but that file **does not exist** (only `localisation.rst`, `perturbation.rst`, `troubleshooting.rst`): broken doc build. REFACTORING C30 says it was created, so it was lost. **Check:** `sphinx-build -W` passes (done on the 5 guides); `DA_guides/index.rst` + `overview.rst` created; `test_the_guides_are_linked_from_the_documentation_index` now passes. Also added the Archie parameters to the catalogue of `perturbation.rst` (pre-existing failing test `test_the_catalogue_mentions_every_parameter`).
- [x] **H2 [VERIFIED] DONE s2 (partly)** `tests/test_perturbate_seed.py` opens `"perturbate.py"` with a relative path: "no tests ran" / collection error unless run from `pyCATHY/DA/`. **Check:** `test_perturbate_seed.py` passes from the repo root (done). **Still open:** 8 more test modules fail to import from the repo root for the same reason (cwd-dependent `import ERT_checks` / `open("enkf.py")`): `test_ERT_checks, test_archie_downstream, test_archie_scenarios, test_enkf, test_log10_cov, test_mapping, test_refactor_equivalence, test_transforms`; plus `test_read_mbeconv::test_file_without_a_step_says_so` fails at HEAD. They were not run as a suite by the earlier sessions, so their "255 tests pass" count is not reproducible from the root (new item H2b).
- [ ] **H3 [VERIFIED]** `pyCATHY/IA_files/` (incl. a 235 kB zip holding older copies of `cathy_tools.py`, `cathy_DA.py`, `perturbate.py`) sits inside the package. Stale copies can be mistaken for the source. Move to `doc/ai_notes/` or exclude in `MANIFEST.in`. **Check:** `python -c "import pyCATHY, os; print(os.path.exists(os.path.join(os.path.dirname(pyCATHY.__file__),'IA_files')))"` → False.

---

## 3. Bugs in parameter perturbation

Vocabulary: **uniform / global** = one scalar for the whole domain (`pert_control_name = "None"`; also the `sampling_type = "uniform"` distribution, both treated below); **layers** = one value per soil layer (keys `<name>0..<name>NSTR-1`); **zones / root_map** = one value per soil zone / vegetation class (keys `<name>0..`, zone `k` is zone number `k+1`).

### 3.1 Global (scalar) and `sampling_type`

- [ ] **P01 [VERIFIED] `sampling_type="uniform"` crashes in the scalar path.** `sampling_dist()` reads `kwargs["minmax_uni"]`, but `perturbate_parm` calls `sampling_dist(sampling_type, mean, sd, N)` with no kwargs in the `ic` and generic branches → `KeyError`. Only `VG_pert_rules` supports `minmax_uni`. **Fix:** take `minmax_uni` from `per_bounds` (uniform on `[min, max]`) or from `mean ± sqrt(3)·sd`; document which. **Check:** new test samples `uniform` for `ic`, `Ks`, `thetar_VG` and the range equals the bounds.
- [ ] **P02 [VERIFIED] `sampling_type` is silently ignored when bounds exist.** In the `ic`, `ZROOT`, `porosity` and generic branches, `_has_bounds(...)` routes to `sampling_dist_trunc` (truncated **normal**) whatever `sampling_type` says. Same for `Archie_pert_rules` (fixed truncated normal). **Fix:** one sampler `sample_prior(sampling_type, mean, sd, bounds, N)` (normal/lognormal/uniform, truncated by bounds) used by every branch; `validate_scenario` warns/raises when a type cannot be honoured. **Check:** test per type × with/without bounds; sample support inside bounds; empirical mean/sd within tolerance.
- [ ] **P03 [VERIFIED] `lognormal` mixes units.** `np.random.lognormal(np.log(mean), sigma=sd)` takes `sd` as the std of log(x) while the scenario `per_sigma` is physical; for `Ks` (sd 1e-5) the spread is nil, for a negative `mean` (`ic`) the result is NaN. **Fix:** convert (`sigma_log = sqrt(log(1+(sd/mean)^2))`, `mu_log = log(mean)-sigma_log²/2`), reject `mean<=0`. **Check:** test recovers mean and sd of the physical sample.
- [ ] **P04 [VERIFIED] Archie ignores `per_bounds`.** `Archie_pert_rules` clips a∈[0,2.5], m,n∈[1.3,3], rFluid≥0 with hard-coded numbers; nb5 keeps `per_bounds` equal by hand. Also dispatch by substring (`"a" in type_parm`; works only because `"Archie"` has an upper-case A). **Fix:** use the entry's bounds, exact names (`a_Archie`, …). **Check:** scenario with bounds `[0.1, 1.0]` for `a_Archie` gives samples inside.
- [ ] **P05 [VERIFIED] Family dispatch by substring.** `perturbate_parm`: `"ic" in type_parm`, `"VG" in type_parm`; `cathy_DA._update_parameter_member`: `name.casefold() in "ks"`, `in "ic"`, `in "porosity"`, `in "atmbc"` (reverse substring: `"s"`, `"k"`, `"i"`, `"or"` all match). A new parameter name can be routed to the wrong family. **Fix:** exact-name tables (`param_tables`), error for unknown names. Same family as REFACTORING R9. **Check:** test with names `"s"`, `"icx"`, `"Ks_x"` → clear `ValueError`.
- [ ] **P06 [IA B31, VERIFIED in generator comment] Bounds act on the *factor* for `multiplicative` + `layers` but on the *physical value* when clipping the analysed values.** `make_scenario_nc_VGP.py` avoids it with `per_type="None"`. **Fix:** make the validator reject (or convert) `multiplicative` + `per_bounds`. **Check:** validator test.
- [ ] **P07 [IA B3/T12, VERIFIED by design] Homogeneous scalar `ic` makes a rank-1 ensemble** (all 2415 data constrain one number). The VGP-by-layers scenarios were created to break it. **Check:** `describe_ensemble` / `rank1_ensemble` log shows < 100 % of predicted-data variance in the first mode for the scenario you run.

### 3.2 By layers (`pert_control_name = "layers"`)

- [ ] **P10 [VERIFIED] `sampling_type` is ignored (always normal / truncated normal).** `perturbate_parm_by_layers` has no `sampling_type` argument. → merged into P02.
- [ ] **P11 [MEASURED] Sobol sampler: no benefit, not reproducible, wasteful.** `qmc.Sobol(d=NENS).random(n=nlayers)`: points = layers, dimensions = members (stratification would matter across members, not across layers). 200 repetitions (N=32, L=15, mean 1.4, sd 0.07): ensemble sd per layer 0.0696 vs 0.0697 for plain iid draws; mean |corr| between members 0.217 vs 0.217; two calls differ and `np.random.seed` has no effect (scrambled Sobol uses its own entropy); `n=15` is not a power of 2 (scipy balance warning, silenced in my run). Each key (45 for the VGP) re-draws the whole `(nlayers, N)` matrix and keeps one row; `ini_<key>_withLayers` differs per key. nb5 works around this with a pickled prior. **Fix:** one draw per parameter family from the global RNG (`rng.standard_normal((nlayers,N))` through the truncated quantile), stored once, seeded by `DAConfig.seed`. **Check:** same seed ⇒ identical prior; prior cache in nb5 no longer needed.
- [ ] **P12 [VERIFIED, existing FIXME] Layers, materials and parameters are drawn independently.** alpha/n/thetar of one layer are uncorrelated, and the 6+3+6 layers of a material (A 0–5, B 6–8, C 9–14) have independent values although they are one material. Botto et al. 2018 route in the FIXME of `make_scenario_nc_VGP.py`. **Fix (staged):** (a) new control `material` (shared draw across a layer group; nb5 goes from 45 to 9 parameters); (b) `per_corr` + Cholesky on Johnson-transformed variables (`transf_type`); (c) check the analysis keeps the correlation. **Check:** sample correlation matrix ≈ requested `per_corr`; EnKF update of a synthetic linear problem keeps correlation.
- [ ] **P13 [VERIFIED] Layer numbering is 1-based when the table is built, 0-based once read back.** `init_soil_SPP_map_df` → layers 1..NSTR; `read_soil` doc: "layers indexing starts from 0"; `_assign_soil_value` uses the 0-based key. `nc_model.setup_model` uses label slices `1:6`, `7:10`, `10:` (inclusive, **layer 10 is in two slices**, the last wins), which equals the generator's 0-based `0–5 / 6–8 / 9–14` only by write order. **Fix:** one `MATERIALS` constant shared by `nc_model.py` and `make_scenario_nc_VGP.py`, non-overlapping; keep the nominal-vs-soil-table check of nb5 (it is good). **Check:** test builds both and asserts equality layer by layer.
- [x] **P14 [CONFIRMED + FIXED s2] Soil-file write condition.** `_soil_tables_to_fill` sets NaN only in the columns of the **first** parameter handled; `_update_soil_params_from_perturbations` writes the file when `not df.isna().any().any()`. Consequences to test: (a) a `listUpdateParm` that lists only some layers leaves NaN → the soil file is **never written, silently**; (b) once the first parameter is complete, every later layer of `n_VG`/`thetar_VG` rewrites the whole file (≈45 writes × NENS × cycle, `backup=True`); (c) any NaN in an unrelated column blocks the write. **Fix:** fill from the existing table (no NaN), write once per member after the loop. **Check:** `tests/test_soil_write_once.py` (3 tests; 2 fail on the original code: partial layers gave 0 writes, SATCELL+n_VG gave 5 writes/member instead of 1). Fix: Ks/porosity/VG keys share one table per member (`_EnsembleUpdate.soil`), no NaN fill, written once after the key loop (`_write_updated_soil`); the Feddes branch reads that shared table and the final write uses the updated Feddes table, so Feddes+soil keys no longer overwrite each other in key order. `NEEDS-RUN`: Feddes + VG in the same `listUpdateParm` (not covered by the stub test).
- [ ] **P15 [SUSPECTED, test first] Layered `ic` is re-applied every cycle.** `_update_ic_member`: the non-layered branch is guarded by `cycle_nb == 0`, the `layers` branch is not, so at `cycle_nb > 0` it may overwrite the analysed pressure head with the prior (depends on the order with `_update_states_from_analysis`). Not hit by nb5 today (`ic` control is `"None"`), but `IC_TRUE` of the truth model is layered. **Check:** stub test calling both in the real order over 2 cycles.
- [ ] **P16 [IA, VERIFIED] `alpha_VG` and `VGP` are rejected** (cannot be written to the members / code commented out): perturb `SATCELL_VG`. Keep; document (done in the .rst). Optional: implement `alpha_VG` → `VGPSATCELL = -1/alpha`.
- [ ] **P17 [design, NEEDS-RUN] 45 (or 48) parameters with NENS = 32 and ERT only.** Rank-deficient covariance ⇒ spurious correlations. Options in order: P12(a) 9 parameters; `AnalysisConfig.inflation_params` / `jitter_params` (check they are really applied in `enkf_Evensen2009_Sakov`; nb5 says inflation is not used by this filter); NENS ≥ 100; update states first (scenario `_pert`) then parameters. **Check:** compare scenarios `_pert` vs `_upd` on the same prior (RMS ERT, parameter spread).

### 3.3 By zones / `root_map`

- [ ] **P20 [IA B20, validator only] Zones and vegetation classes must be numbered 1..K without gaps** (index `k` ↔ zone/class `k+1`). Detected, behaviour unchanged. **Check:** nb5-like case with classes {1,3} → `validate_scenario` error.
- [ ] **P21 [IA B27, decision needed] `root_map`/`zone`/`layers` entries can overlap on the same soil rows; the last write wins silently** (reported by `describe_ensemble`). **Decision:** error, average, or last-wins documented. **Check:** `describe_ensemble` warning turns into an error if you choose "error".
- [ ] **P22 [VERIFIED, SUSPECTED effect] `root_map` assigns all layers of the zones of a class** (`df.loc[zone[index+1 == veg_map], cols]`), and `zone` may contain the nodata value (-9999) that B19 only fixed in the *analysis*. **Check:** toy map with nodata: `_assign_soil_value` does not raise `KeyError` and writes only valid zones.
- [ ] **P23 [IA B32/C22] The combined `zone`+`root_map`+`layers` layout described in a code comment is not implemented.** Delete the comment or implement (nested lists). Interacts with P12(a).
- [ ] **P24 [IA B24] Localisation extrudes surface observations to all depths** (wrong for ERT, whose sensitivity varies with depth). Not used by nb5 (no localisation) but needed if the parameter count stays high.

---

## 4. Items for `nb5_DA_ERT_soil3D_nc.py`

- [ ] **N01 [VERIFIED, DECISION s3] `DAConfig(psi_bounds=...)` does not exist at HEAD** and you decided **not** to add it (clipping implemented in s2, then fully removed in s3: code, test, docs). **Remaining fix is in nb5 only:** delete `PSI_BOUNDS` and the `psi_bounds=PSI_BOUNDS` argument of `DA_CONFIG` (and its comment block). The positive-ψ problem stays reported by the B10 warning; address its cause (prior, gain, ERT alignment, N10) instead. Open design question: reject members with positive analysed ψ? (not decided). **Check:** nb5 starts on the library at HEAD.
- [ ] **N02 `ert_alignment.py` missing** (section 1). After receipt: move it to `pyCATHY/ERT/` with tests (see V3) so your colleague does not need a side file. **Check:** `from pyCATHY.ERT.alignment import ...` works; synthetic test of `check_ert_cathy_alignment` (inside/outside, z-sign).
- [ ] **N03 `ERT_TIMES_H` are placeholders.** The whole time alignment (window-start stamps, snapping ≤ 30 min) depends on them. **Fix:** read the acquisition time from the `.tldat`/log, or write the real ones; print the table (nb5 already does) and keep `check_time_alignment` errors fatal. **Check:** your sign-off + `check_time_alignment` has no error.
- [ ] **N04 [SUSPECTED] Frame of the forward model.** Header says electrodes and mesh are now transformed together; section 3 says the forward runs in the **raw** pygimli frame and the transform is used only for the alignment check and `mesh_nodes_modif`. Both can be right but a colleague cannot tell. **Fix:** read `ERT/simulate_ERT.py` (in repo) and state in one comment which frame each object is in; add an end-to-end check that maps a known CATHY saturation field to the ERT cells (a layered ρ field gives the expected apparent resistivity). **Check:** `NEEDS-RUN` (pygimli).
- [ ] **N05 Prior cache is keyed only by scenario name and NENS** (`prior_ensemble_<SCENARIO>_NENS<N>.pkl`): edit the toml (same name) and a stale prior is silently reused. **Fix:** hash the scenario table + seed in the file name; becomes optional after P11. **Check:** editing a sigma changes the cache name.
- [ ] **N06 `ERT_STRICT = False`** lets the run continue after a failed forward test. Fine for exploration, wrong default for a colleague. **Fix:** default `True`; the quadrupole mask only protects the *observations*, not bad *predictions*. **Check:** run stops on the old misaligned geometry (reproduce with a wrong shift).
- [ ] **N07 [IA B28] `np.random.seed(SEED)` is set before `perturbate()`; the Sobol part ignores it** (P11). After P11, set the seed once through `DAConfig.seed` only.
- [ ] **N08 [design] Filter choice.** `enkf_Evensen2009_Sakov` + `perturb_obs=True` (good: C1 + C6), inflation 1.0 (unused by this filter). With P17 decide inflation/localisation; compare `USE_LOG_OBS` True/False and `ERT_DATA_ERR` 3/5/10 % (R2, B1). **Check:** `NEEDS-RUN` table of RMS/spread per option.
- [ ] **N09 [IA known_issues] "ERT assimilate in non-parallel model leads to wrong observation size"**: keep `parallel=True`, add this to the run checklist; test the non-parallel path or document it as unsupported. **Check:** `NEEDS-RUN`.
- [ ] **N10 Cycle-5 rejection of the last run** (all members rejected): hypothesis = ERT/CATHY mesh misalignment (already addressed by the alignment gate) plus positive ψ after analysis (N01). **Check:** re-run `nc_ERT_ic_VGP_pert` with N01+N03 done: `<project>_rejections.json` shows 0 members rejected for `after_analysis`. `NEEDS-RUN`.
- [ ] **N11 Section 7 (results) is commented out**; results are in the plots script. Add a final `print` of the paths of `df_DA`, `ET_DA_xr`, `_ERT_checks.json`, `_rejections.json` (what the colleague must open).
- [ ] **N12 [fine, keep] nominal-vs-soil-table check, `describe_ensemble` + `assert` on `ERROR`, `check_update_parameters`, mesh reader consistency asserts.** Do not remove when refactoring.

---

## 5. Items for `nb5_DA_ERT_soil3D_nc_plots.py` (+ `nb5_plot_helpers.py`)

I only read the structure of these two files (functions listed; `ts_figure`, `profiles_figure`, `spread_figure`, `scatter_figure`, `rms_figure`, `ert_hist`, `parm_figure`; helpers for SMC/TDR). A line-by-line review is the first task of the plot session.
- [ ] **G01 Duplicated configuration.** The header says "set `path2prj` / `PARAM_NAMES` / `SCENARIO` exactly as in nb5": two copies to keep in sync. **Fix:** one `nb5_config.py` (or read the scenario toml + project name) imported by both. **Check:** changing `SCENARIO` once changes both.
- [ ] **G02 `PARAM_NAMES` must come from the scenario** (45 VGP keys: plot per material/layer group, not 45 panels). **Check:** `parm_figure` for `_upd` shows 3 parameters × layers as one panel each.
- [ ] **G03 Figures that answer the validation questions** (section 8): prior vs posterior vs open loop ERT RMS per cycle; ensemble spread of ψ/θ at the TDR depths vs sensor; predicted vs observed apparent resistivity scatter (log axes) per survey; parameter shift in units of prior sd (detects the A1/A4 blow-up pattern); rejected members per cycle from `_rejections.json`.
- [ ] **G04 Failure behaviour.** `safe(fn)` wrappers exist; make each failing figure print the missing input (file/column) instead of a bare traceback. **Check:** run with a deleted `df_DA` backup → readable message.
- [ ] **G05 Needs `nb5_plots_suggestions.md`** (section 1).

---

## 6. Order of work (first unchecked item first)

1. Section 1 inputs → 2. N01 (script cannot start) → 3. H1–H3 → 4. P14, P15, P22 **tests only** (decide real bugs) → 5. P01–P06 (one sampler) → 6. P11 + N05 + N07 → 7. P13 → 8. N02, N03, N04, N06 → 9. P12(a) then P17 → 10. G01–G05 → 11. docs (`DA_overview_implemented.rst`, V-protocol) → 12. IA leftovers (B16, B27, B28, C13).

---

## 7. Merge table with `IA_files` (so nothing is done twice)

| IA id | What | In HEAD? | Left to do (this plan) |
|---|---|---|---|
| A1/C1 | Sakov gain omitted HPHᵀ | yes (`_solve_obs_space`) | `NEEDS-RUN` confirmation on nb5 |
| A2/C2, B6 | bounds + clipping after analysis, sorted `ens_valid` | yes | P06 (factor vs physical) |
| A3/C3, C20/B31 | `ic` and layered/zoned priors honour `per_bounds` | yes | P02 (type ignored) |
| B1/C4 | ERT error in % | yes (`data_err_unit`) | nb5 already uses it |
| B2/C5, B37/C32 | quadrupole mask (flag, never delete) + localisation | yes | N06 (predictions not covered) |
| B13/C6 | `perturb_obs` | opt-in | nb5 enables it; author decides default |
| B14, B17, B18, B19, B21, B22, B25–B27, B30 | various perturbation fixes | yes | P21 decision (B27) |
| B20 | gaps in zone numbering | validator only | P20 |
| B28 | print→logger missing in 4 modules | no | IA leftovers |
| B32/C22 | combined controls | not implemented | P23 / P12(a) |
| B34/C28 | obs/state time pairing | detection only | N03 |
| B35/B36/C31 | `alpha_VG`, Feddes layouts | rejected/validated | P16 |
| C15/B5 | rejection ledger | yes | N10 |
| C29 | diagonal R | opt-in (used by nb5) | diagonal-aware filters |
| C30 | doc guides | **partly lost** (no `index.rst`) | H1 |
| R4b | pyGIMLi noise seed per member | open | irrelevant while `fwdNoiseLevel = 0` (nb5) |

---

## 8. Validation protocol for a colleague (V0–V10)

Goal: reproduce and judge a run without knowing the history. Each step has a pass criterion; stop at the first failure.

- **V0 Environment.** Clone the repo, `pip install -e .`, install pygimli. `pytest pyCATHY/tests -q` → all pass (baseline today: 119 passed on the subset above).
- **V1 Scenario.** `python make_scenario_nc_VGP.py` then open `scenario_nc_ERT.toml`. Pass: every list in a layered entry has 15 values; nominals equal the soil table (nb5 prints a WARNING otherwise).
- **V2 Time table.** nb5 prints "ERT file / acquired / DA stamp". Pass: stamps strictly increasing, no duplicate window, `check_time_alignment` has no `[ERROR]`.
- **V3 Geometry.** Alignment gate passes (electrodes on mesh nodes, inside the CATHY volume, top = soil surface). Pass: no `ERTAlignmentError`; view once with `SHOW_MESH_ALIGNMENT=True`.
- **V4 Forward test.** Homogeneous 30 Ω·m medium: median rhoa/ρ ≈ 1, > 50 % quadrupoles within ±10 %. Pass: no `ERT forward test FAILED` warning.
- **V5 Ensemble description.** The `describe_ensemble` table: each updated parameter touches the intended rows, no overlap warning.
- **V6 Prior sanity.** Histograms (`SAVE_PRIOR_FIGS=True`): inside bounds, correct spread; open-loop run (`WITH_DA=False`) uses the same prior.
- **V7 Cycle 1.** `<project>_rejections.json`: 0 members rejected before the first analysis; analysis log shows `max|K d|` reasonable; analysed ψ ≤ 0.
- **V8 Whole run.** Rejected < `threshold_rejected`; parameter shifts within ~3 prior sd; ensemble spread does not collapse (< 10 % of prior after cycle 2 = warning).
- **V9 Skill.** ERT RMS (log10) DA < open loop at every survey; TDR/SMC comparison (plots script) improves at the sensor depths.
- **V10 Report.** Save: toml, log, `_ERT_checks.json`, `_rejections.json`, figures; fill the table of section 10.

---

## 9. Deliverables

- `DA_ERT_improvement_plan.md` (this file).
- `DA_overview_implemented.rst` for `doc/content/DA_guides/` (add to the toctree, see H1). Written from the code at HEAD `4cb4d072`; where an item is only planned it says so (e.g. `psi_bounds`).

---

## 10. Progress ledger (update every session)

| Date | Session | Items ticked | Tests run | Notes |
|---|---|---|---|---|
| 2026-10-20 | 1 | none (plan only) | baseline 119 passed | clone + read `IA_files`, nb5, make_scenario, nc_model, plots structure; measured Sobol sampler (P11) |
| 2026-10-20 | 2 | H1, H2 (partly), P14 (N01 added then removed, see s3) | full synthetic suite: no new failure vs HEAD; 3 previously failing tests now pass; +8 new tests | patch `session2_fixes.patch` (5 modified + 4 new files; s3: 6 files after removing psi_bounds); Sphinx -W build of the guides OK |
| 2026-10-20 | 3 | N01 reverted by decision | no new failure vs HEAD; P14 tests pass | `psi_bounds` removed (code, test, docs); patch regenerated (6 files); nb5 must drop `psi_bounds` |

**Not done in session 1:** line-by-line review of `nb5_DA_ERT_soil3D_nc_plots.py`/`nb5_plot_helpers.py`; reading `ERT/simulate_ERT.py`, `mapper.py`, `enkf.py` internals beyond what is cited; anything needing pygimli or CATHY.
