# REFACTOR_GFORTRAN.md — pyCATHY without a system gfortran / LAPACK / BLAS

Branch: `refactor/no-gfortran-install` (created from `main`; no existing branch covered this topic).
Status: **plan only, no code changed yet.**

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

=> The dependency is tiny and can be removed from the Fortran source itself (see step 1).
To be confirmed with a real build (no gfortran in the authoring sandbox): link without `-llapack -lblas` and list the undefined symbols.

### Why it recompiles

`CATHY.H` holds `PARAMETER` dimensions (`NODMAX`, `NTRMAX`, `MAXSTR`, `NMAX`, `NTEMAX`, `NPMAX`, `MAXZON`, `MAXVEG`, ...), written by `update_cathyH`. Static array sizes are compile-time, so any change of mesh/DEM size, layers, zones or vegetation classes requires a new `cathy` executable. The preprocessor (`prepro/src/*.f90`) reads a parameter file at run time: it does **not** depend on the case and only needs compiling once.

### Portability blockers found

1. Compiler name hard-coded to `gfortran`. The conda `gfortran_*` packages install a prefixed binary (e.g. `x86_64-conda-linux-gnu-gfortran`) and expose it through `$FC` after activation, so a bare `gfortran` is often not on `PATH`.
2. `-mcmodel=large` is x86-64 only; expected to fail on Apple Silicon (arm64). To verify.
3. `os.system(... > /dev/null 2>&1)` and `./pycppp` do not work on Windows (no `/dev/null`, executable is `.exe`).
4. Global `os.chdir` + `glob` on the current directory (fragile, not thread/process safe for DA ensembles).
5. Compiler failures in the preprocessor are silently swallowed.

## Plan

### 1. Remove LAPACK/BLAS from the Fortran source
- Add `pycathy_blas_min.f` with a plain `dcopy`.
- In `iperplane.f` replace `dgetrf/dgetrs` by a small 4×4 Gaussian elimination with partial pivoting (`solve4x4`).
- Drop `-llapack -lblas` from the link command.
- Apply to the sources fetched by `_clone_cathy_src` (`CATHY-Org/CATHY_src_org`) and to `CATHY_SCF_LAI` (`version="SCF_variable"`); bump `CATHY_SRC_COMMIT`.
- Validate: same outputs as the LAPACK build on the Weill et al. test case (tolerance to define).

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
| **C. Vegetation (and other physical) parameters at run time** | Read the monthly vegetation parameters from an input file instead of compiling them in; only true dimensions stay in `CATHY.H`. | Fortran change in CATHY-Org; removes your most frequent recompile reason. Depends on where those parameters live today (question 1). |
| D. Dimension buckets | Round dimensions up to size classes so many projects share one executable. | Larger static arrays, more memory, interacts with `-mcmodel=large`. |
| E. Prebuilt binaries per OS | Ship executables for reference cases. | Does not cover arbitrary DEMs. |
| F. `ccache` (conda) | Compiler-level cache. | Easy add-on, little gain with A+B. |

Proposal: **A + B now, C for the vegetation parameters**. Preprocessor: compile once per pyCATHY version (cache) and drop the double compile.

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
- `solve4x4` vs `numpy.linalg.solve` on random well-conditioned and pivot-requiring matrices.
- Failure tracking: a deliberately broken `.f` is reported, not swallowed.

## Open questions (need your input)

1. **Monthly vegetation parameters**: where are they today — hard-coded in a `.f` (e.g. a block-data / table) of `CATHY_SCF_LAI`, or only in `CATHY.H`? Which file do you edit when you change them? This decides whether option C is feasible.
2. Which sources are the target: `CATHY_src_org` only, or also `CATHY_SCF_LAI` (`version="SCF_variable"`)? Do I patch them upstream (fork + PR in CATHY-Org) or vendor patched sources inside pyCATHY?
3. Windows: native support is a hard requirement, or is "works on Linux/macOS conda, Windows best effort" acceptable for the first release?
4. Largest problems you run (nodes, layers): do you rely on `-mcmodel=large` (static arrays > 2 GB)? Is Apple Silicon a target?
5. In DA ensembles, do all members share the same `CATHY.H` (same dimensions)? Today each member may recompile.
6. Build cache location: inside the project folder (portable, visible) or a user cache dir (shared across projects)?
7. Distribution: pip only (pip cannot install a compiler) or also a conda-forge recipe with `fortran-compiler` as run dependency?
8. Numerical tolerance accepted when replacing LAPACK `dgetrf/dgetrs` by the 4×4 solver (bit-for-bit is not guaranteed).
