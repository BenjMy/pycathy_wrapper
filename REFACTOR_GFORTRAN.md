# REFACTOR_GFORTRAN.md — pyCATHY without a system gfortran / LAPACK / BLAS

Branch: `refactor/no-gfortran-install` (created from `main`; no existing branch covered this topic).
Status: **steps 1, 3, 4A, 5 and the Linux part of 6/7 implemented** (see "Implementation status"); macOS and Windows still to be validated in CI.

## Implementation status

| Step | Status | Notes |
|---|---|---|
| 1. LAPACK/BLAS without touching the sources | DONE (Linux) | `link_libraries()` adds the conda `lib` dir + rpath, then `-llapack -lblas`; Fortran sources untouched. |
| 2. Compiler from conda | DONE (docs/env) | `environment.yml`, `installing.rst`: `fortran-compiler libblas liblapack`. Windows/macOS conda packages not tested here. |
| 3. `pyCATHY/fortran_build.py` | DONE | compiler lookup, platform flags, parallel `.f` compile, F90 module ordering (replaces the "compile twice" hack), no `chdir`/`os.system`. |
| 4A. Content-hash build cache | DONE | user cache dir, `$PYCATHY_CACHE`; key = sources + `*.H` + compiler version + flags. |
| 4B. Parallel build | DONE | thread pool over `*.f`. |
| 5. Wrapper clean-up | DONE | `run_preprocessor` raises with the compiler output on failure; `recompileSrc(verbose=False, force=False)` (backward compatible), `CATHY.last_build` keeps the `BuildResult`; `processor_name` is `cathy.exe` on Windows. |
| 6. CI + docs | PARTIAL | `.github/workflows/fortran-build.yml` (Linux, macOS Intel/ARM, Windows) **written but never run**. |
| Distribution | DRAFT | `conda-recipe/meta.yaml` (untested; sha256, PyPI url and noarch/`fortran-compiler` question are TODO). `documentation.yml` still uses `apt` (docs build only; not migrated). |
| 7. Tests | DONE | `pyCATHY/tests/test_fortran_build.py` (19 tests, real builds skipped without a compiler). |

Verified in the authoring sandbox (Linux x86-64, system gfortran + libopenblas, no conda): preprocessor builds, `cathy` builds from the Weill et al. sources in ~4.5 s (parallel), second build is a cache hit (0.0 s), and `CATHY.run_preprocessor()` + `run_processor()` run end to end. `CATHY_SCF_LAI` sources also build (cathy 4.6 s, prepro 1.7 s). **Not verified**: conda gfortran packages, macOS, Windows, a full `SCF_variable` simulation, results identical to the old build bit by bit (same sources/flags/libraries, so expected).

## Goal

`conda|mamba create -n pycathy ... && pip install pyCATHY` works on Linux, macOS (Intel + Apple Silicon) and Windows, with:

- no `apt install gfortran libopenblas-dev` / `dnf install ...`
- no WSL on Windows
- the Fortran compiler coming from the conda/mamba environment
- no `-llapack -lblas` at link time

## Current state (what the code and docs do today)

| Where | What |
|---|---|
| `doc/content/installing.rst` | Requires system Fortran + BLAS + LAPACK. Ubuntu: `apt install gfortran libopenblas-dev` (22) or `blas-dev lapack-dev` (20); Fedora: `gcc-gfortran blas-devel lapack-devel`; Windows: "install WSL". |
| `.github/workflows/documentation.yml` | `apt-get install gfortran liblapack-dev libblas-dev`. `main.yml` (tests) installs no compiler. |
| `environment.yml` | Python 3.10, git, pip, numpy. No compiler. |
| `CATHY.run_preprocessor` (`cathy_tools.py`) | Compiles `prepro/src/*.f90` with hard-coded `gfortran -O -mcmodel=large -fPIC`, **twice per call**, via `os.system(... > /dev/null 2>&1)`; runs `./pycppp`. |
| `CATHY.recompileSrc` (`cathy_tools.py`) | `gfortran -O -mcmodel=large -fPIC -c` on every `*.f` (280 files), then link with `-llapack -lblas`. Relies on `os.chdir(.../src)` and `glob("*.o")`. Called from `run_processor(recompile=True)` (default). |

### What actually needs LAPACK/BLAS (checked on `pyCATHY/tests/weil_exemple/my_cathy_prj/src`)

- LAPACK: only `iperplane.f` → `dgetrf` + `dgetrs` on a **4×4** system.
- BLAS: `dcopy` (called in `dxpay.f`, `normres.f`, `precond.f`; **not defined locally**).
- `ddot`, `dnrm2`, `daxpy` are already defined in `solscal-extended.f` (no external lib needed).

=> The dependency is tiny, but the sources stay untouched: the libraries come from conda (see step 1).
To be confirmed with a real build (no gfortran in the authoring sandbox): link without `-llapack -lblas` and list the undefined symbols.

### Why it recompiles

`CATHY.H` holds `PARAMETER` dimensions (`NODMAX`, `NTRMAX`, `MAXSTR`, `NMAX`, `NTEMAX`, `NPMAX`, `MAXZON`, `MAXVEG`, ...), written by `update_cathyH`. Static array sizes are compile-time, so any change of mesh/DEM size, layers, zones or vegetation classes requires a new `cathy` executable. The preprocessor (`prepro/src/*.f90`) reads a parameter file at run time: it does **not** depend on the case and only needs compiling once.

### Portability blockers found

1. Compiler name hard-coded to `gfortran`. The conda `gfortran_*` packages install a prefixed binary (e.g. `x86_64-conda-linux-gnu-gfortran`) and expose it through `$FC` after activation, so a bare `gfortran` is often not on `PATH`.
2. `-mcmodel=large` is x86-64 only; expected to fail on Apple Silicon (arm64). To verify.
3. `os.system(... > /dev/null 2>&1)` and `./pycppp` do not work on Windows (no `/dev/null`, executable is `.exe`).
4. Global `os.chdir` + `glob` on the current directory (fragile, not thread/process safe for DA ensembles).
5. Compiler failures in the preprocessor are silently swallowed.

## Vegetation: findings (answer to question 1)

Checked in `pycathy_wrapper` (`cathy_tools.py`) and in the Fortran of both `CATHY_src_org` (test project sources) and `CATHY_SCF_LAI` (`src/datin.f`).

- Monthly workflow (Python script: LAI -> lookup table -> root depth + Feddes parameters -> `update_soil` / `update_veg_map`) writes:
  - `input/soil`: per vegetation class `PCANA PCREF PCWLT ZROOT PZ OMGC` (+ `SCF` block in `SCF_variable`), via `_write_SOIL_file`.
  - `input/root_map`: vegetation class raster, via `_write_root_map`.
  - `src/CATHY.H`: **only `MAXVEG`** (number of classes), via `update_veg_map` -> `update_cathyH`.
- Fortran `datin.f` reads those values at run time (`READ(IIN4,*) PCANA(I),PCREF(I),PCWLT(I),ZROOT(I),...` for `I=1,NVEG`; `READ(IIN4,*) (SCF(I),I=1,NVEG)` in SCF_LAI). `MAXVEG` only dimensions arrays (`NVEG.GT.MAXVEG` -> error). No `.f` file holds the vegetation values.
- Hence changing LAI / root depth / Feddes / SCF month by month **never requires recompiling**, provided `MAXVEG` does not change (`update_veg_map(..., maxveg=N)` already keeps it fixed).
- Why it recompiles today anyway: `run_processor(recompile=True)` is the default and `recompileSrc()` rebuilds all `*.f` unconditionally. Immediate workaround: `run_processor(recompile=False)` after the first month. Proper fix: content-hash cache (option A), which makes a rewritten but identical `CATHY.H` a cache hit.
- Note: `CATHY_SCF_LAI` ships `MAXVEG=15` in its `CATHY.H` (template value); wrapper overwrites it per project.
- Not seen yet: your LAI loop script. Please send it to confirm that nothing else (e.g. `MAXZON`, layers, DEM) changes between months.

## Plan

### 1. LAPACK/BLAS without touching the Fortran sources (decision: **no edits to any `.f` / `.f90`**)
- Default: get `libblas` + `liblapack` (conda-forge, OpenBLAS/MKL/Accelerate depending on platform) from the **conda environment** and link with `-llapack -lblas` as today, adding `-L$CONDA_PREFIX/lib` (`%CONDA_PREFIX%\\Library\\lib` on Windows) when needed. Same routines as now, so results are unchanged; nothing is installed system-wide.
- Optional later (only if conda linking is painful on Windows): a wrapper-side shim file (`pycathy/fortran/lapack_shim.f`: `dcopy`, `dgetrf`, `dgetrs`) compiled *next to* the CATHY sources, never inside them. Not needed for the first release.
- Applies to both `CATHY_src_org` and `CATHY_SCF_LAI`; no upstream fork/PR, no vendored patched sources.
- Validate: Weill et al. case gives the same output as the current build.

### 2. Get the compiler from conda/mamba
- Add the conda-forge `fortran-compiler` metapackage (platform-specific gfortran/flang) to `environment.yml` and to the docs; keep `make`/`ninja` optional.
- Windows: decide which toolchain is supported (conda-forge gfortran via MSYS2/m2w64 vs flang). **To test.**
- Rewrite `installing.rst`: one conda path for all OS; keep apt/dnf only as "advanced / system compiler".

### 3. New module `pyCATHY/fortran_build.py` (no `chdir`, no `os.system`)
- `find_fortran_compiler()`: `$PYCATHY_FC` → `$FC` → `<prefix>-gfortran` in `$CONDA_PREFIX/bin` → `gfortran` → `flang`; clear error message with the conda install command if none found.
- `compile_fortran(src_dir, sources, flags, exe)`: `subprocess.run([...], cwd=...)` with list arguments, `pathlib`, `.exe` suffix on Windows, optional parallel `-c` compilation (`ThreadPoolExecutor`), returns a result object (return code, failed files, stderr).
- Flags by platform: `-mcmodel=large` only where supported; drop `-fPIC` for executables unless required.
- Failures are tracked and surfaced in the pipeline report (not swallowed); `recompileSrc` and `run_preprocessor` become thin callers.

### 4. Recompilation strategy (needed: compiling once is not enough)

| Option | Idea | Cost / risk |
|---|---|---|
| **A. Content-hash build cache** (proposed first) | Key = sha256(all `*.f` + `CATHY.H` + compiler version + flags). Hit → copy cached `cathy`; miss → build and store. Changing something that does not alter the key never recompiles; DA ensembles with identical dimensions compile once. | Low. Needs a cache location (project dir vs user cache dir). |
| **B. Incremental objects + parallel build** | Cache one `.o` per source keyed by hash(file + headers it includes); relink only. Compile the rest with N workers. | Most files include `CATHY.H`, so a dimension change still rebuilds nearly everything; parallelism helps most. |
| **C. Vegetation parameters at run time** | **Already the case** (see "Vegetation: findings"): Feddes/root/SCF values are read from `input/soil` at run time. Only `MAXVEG` (number of classes) is compile-time. The monthly loop only needs `recompile=False` (or a cache hit). | No Fortran change needed. |
| D. Dimension buckets | Round dimensions up to size classes so many projects share one executable. | Larger static arrays, more memory, interacts with `-mcmodel=large`. |
| E. Prebuilt binaries per OS | Ship executables for reference cases. | Does not cover arbitrary DEMs. |
| F. `ccache` (conda) | Compiler-level cache. | Easy add-on, little gain with A+B. |

Proposal: **A + B now; C is already satisfied** (monthly LAI loop: no recompile as long as `MAXVEG` is constant). Preprocessor: compile once per pyCATHY version (cache) and drop the double compile.

### 5. Wrapper clean-up
- `run_preprocessor` / `run_processor`: absolute paths, `subprocess.run(cwd=...)`, executable name from `fortran_build`, no `./`.
- Keep `recompile=True/False` and `recompileSrc()` signatures (backward compatible); add `force_recompile=False`.

### 6. CI and docs
- GitHub Actions matrix: `ubuntu-latest`, `macos-13` (Intel), `macos-14` (arm64), `windows-latest`, using `conda-incubator/setup-miniconda` with conda-forge; **no** `apt install gfortran`.
- Smoke test per OS: build prepro + `cathy` on the Weill example and compare against a reference output.
- Update `doc/content/installing.rst`, `environment.yml`, `documentation.yml`; remove the WSL instructions once Windows is validated.

### 7. Tests (assertion-based)
- `find_fortran_compiler` with mocked `PATH`/env (conda-prefixed name, missing compiler).
- Build cache: hit/miss on a changed `CATHY.H`, unchanged source, changed flag.
- Link-flag construction per platform (conda prefix present / absent).
- Failure tracking: a deliberately broken `.f` is reported, not swallowed.

## Decisions (from the user)

| Topic | Decision |
|---|---|
| Fortran sources | **Do not modify.** Both `CATHY_src_org` and `CATHY_SCF_LAI` are targets. |
| Windows | **Native Windows is mandatory** (conda install on Windows is the expected route). |
| Apple Silicon | Nice to have, not mandatory. |
| DA ensembles | All members share the same `CATHY.H` (to confirm with the LAI loop script / a DA run). |
| Distribution | pip **and** a conda-forge recipe (`fortran-compiler`, `libblas`, `liblapack` as run dependencies). |
| Numerical tolerance | Not an issue: real LAPACK from conda is kept (step 1), results unchanged. |
| Build cache location (my choice, revisit if needed) | User cache dir (`platformdirs`), override with `PYCATHY_CACHE`; ensemble members and projects share binaries. |

## Open points / risks to test in CI (not decided by questions)

1. **Windows native toolchain**: conda-forge Fortran on win-64 plus linking against conda's LAPACK/BLAS; the Python side (`/dev/null`, `./pycppp`, `.exe`) is handled in step 3. Possible separate Windows recipe (different compiler package/flags) if the unified one fails.
2. **`-mcmodel=large`**: x86-64 only; probably unsupported on Apple Silicon and of limited use on Windows (2 GB image limit). Question left for you: do your largest meshes need more than 2 GB of static arrays? If not, the flag can be dropped on those platforms.
3. **Shared `CATHY.H` in ensembles**: verify with a DA run before relying on one compile per ensemble.
4. **Your LAI loop script**: still needed to confirm that nothing but vegetation inputs changes between months.
