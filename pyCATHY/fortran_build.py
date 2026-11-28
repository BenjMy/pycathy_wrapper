"""
Portable Fortran build helpers for pyCATHY.

Replaces the hard-coded ``gfortran`` + ``os.system`` + ``os.chdir`` calls:

- the compiler is looked up (``$PYCATHY_FC``, ``$FC``, conda environment, ``PATH``);
- LAPACK/BLAS are taken from the active conda environment when there is one;
- no shell, no ``/dev/null``, no ``chdir`` (works on Linux, macOS and Windows);
- the CATHY sources are compiled in parallel and the executable is cached by
  content hash, so an unchanged ``CATHY.H`` / source set / compiler / flags never
  triggers a recompilation (the monthly vegetation loop, DA ensemble members...).

The Fortran sources are never modified.
"""

from __future__ import annotations

import hashlib
import os
import platform
import re
import shutil
import stat
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

ENV_COMPILER = "PYCATHY_FC"
ENV_CACHE = "PYCATHY_CACHE"
ENV_FFLAGS = "PYCATHY_FFLAGS"

# Order matters (Fortran 90 modules): same list as the former hard-coded command.
PREPRO_SOURCES = (
    "mpar.f90 mbbio.f90 wbb_sr.f90 csort.f90 qsort.f90 depit.f90 cca.f90 "
    "smean.f90 dsf.f90 facet.f90 hg.f90 mrbb_sr.f90 bb2shp_sr.f90 shape.f90 "
    "dbase.f90 streamer.f90 cppp.f90"
).split()


ENV_BLAS = "PYCATHY_BLAS"  # "shim" (default) or "system"

# CATHY only needs three external routines: DCOPY (BLAS), DGETRF and DGETRS (LAPACK, always
# called on a 4x4 system by IPERPLANE). The wrapper compiles this reference-algorithm
# implementation next to the CATHY sources (they stay untouched), so no BLAS/LAPACK library
# has to be installed or linked on any platform. PYCATHY_BLAS=system links the real ones.
LAPACK_SHIM_NAME = "pycathy_blas_lapack_min.f"
LAPACK_SHIM_SOURCE = """\
C     Minimal DCOPY / DGETRF / DGETRS (reference algorithms, unblocked)
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
      INTEGER N,INCX,INCY
      DOUBLE PRECISION DX(*),DY(*)
      INTEGER I,IX,IY
      IF (N.LE.0) RETURN
      IX=1
      IY=1
      IF (INCX.LT.0) IX=(-N+1)*INCX+1
      IF (INCY.LT.0) IY=(-N+1)*INCY+1
      DO 10 I=1,N
         DY(IY)=DX(IX)
         IX=IX+INCX
         IY=IY+INCY
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE DGETRF(M,N,A,LDA,IPIV,INFO)
      INTEGER INFO,LDA,M,N
      INTEGER IPIV(*)
      DOUBLE PRECISION A(LDA,*)
      INTEGER I,J,JP,K
      DOUBLE PRECISION PIV,TMP
      INFO=0
      IF (M.LT.0 .OR. N.LT.0 .OR. LDA.LT.MAX(1,M)) THEN
         INFO=-1
         RETURN
      END IF
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
      DO 60 J=1,MIN(M,N)
         JP=J
         PIV=ABS(A(J,J))
         DO 10 I=J+1,M
            IF (ABS(A(I,J)).GT.PIV) THEN
               PIV=ABS(A(I,J))
               JP=I
            END IF
   10    CONTINUE
         IPIV(J)=JP
         IF (A(JP,J).NE.0.0D0) THEN
            IF (JP.NE.J) THEN
               DO 20 K=1,N
                  TMP=A(J,K)
                  A(J,K)=A(JP,K)
                  A(JP,K)=TMP
   20          CONTINUE
            END IF
            DO 30 I=J+1,M
               A(I,J)=A(I,J)/A(J,J)
   30       CONTINUE
         ELSE IF (INFO.EQ.0) THEN
            INFO=J
         END IF
         DO 50 K=J+1,N
            DO 40 I=J+1,M
               A(I,K)=A(I,K)-A(I,J)*A(J,K)
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      RETURN
      END
C
      SUBROUTINE DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      CHARACTER TRANS
      INTEGER INFO,LDA,LDB,N,NRHS
      INTEGER IPIV(*)
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
      INTEGER I,J,K,IP
      DOUBLE PRECISION TMP
      LOGICAL NOTRAN
      INFO=0
      NOTRAN=(TRANS.EQ.'N' .OR. TRANS.EQ.'n')
      IF (.NOT.NOTRAN .AND. TRANS.NE.'T' .AND. TRANS.NE.'t'
     &    .AND. TRANS.NE.'C' .AND. TRANS.NE.'c') THEN
         INFO=-1
         RETURN
      END IF
      IF (N.LT.0 .OR. NRHS.LT.0 .OR. LDA.LT.MAX(1,N)
     &    .OR. LDB.LT.MAX(1,N)) THEN
         INFO=-2
         RETURN
      END IF
      IF (N.EQ.0 .OR. NRHS.EQ.0) RETURN
      IF (NOTRAN) THEN
         DO 20 I=1,N
            IP=IPIV(I)
            IF (IP.NE.I) THEN
               DO 10 K=1,NRHS
                  TMP=B(I,K)
                  B(I,K)=B(IP,K)
                  B(IP,K)=TMP
   10          CONTINUE
            END IF
   20    CONTINUE
         DO 60 K=1,NRHS
            DO 40 J=1,N
               IF (B(J,K).NE.0.0D0) THEN
                  DO 30 I=J+1,N
                     B(I,K)=B(I,K)-B(J,K)*A(I,J)
   30             CONTINUE
               END IF
   40       CONTINUE
            DO 55 J=N,1,-1
               IF (B(J,K).NE.0.0D0) THEN
                  B(J,K)=B(J,K)/A(J,J)
                  DO 50 I=1,J-1
                     B(I,K)=B(I,K)-B(J,K)*A(I,J)
   50             CONTINUE
               END IF
   55       CONTINUE
   60    CONTINUE
      ELSE
         DO 120 K=1,NRHS
            DO 80 J=1,N
               TMP=B(J,K)
               DO 70 I=1,J-1
                  TMP=TMP-A(I,J)*B(I,K)
   70          CONTINUE
               B(J,K)=TMP/A(J,J)
   80       CONTINUE
            DO 100 J=N,1,-1
               TMP=B(J,K)
               DO 90 I=J+1,N
                  TMP=TMP-A(I,J)*B(I,K)
   90          CONTINUE
               B(J,K)=TMP
  100       CONTINUE
  120    CONTINUE
         DO 150 I=N,1,-1
            IP=IPIV(I)
            IF (IP.NE.I) THEN
               DO 140 K=1,NRHS
                  TMP=B(I,K)
                  B(I,K)=B(IP,K)
                  B(IP,K)=TMP
  140          CONTINUE
            END IF
  150    CONTINUE
      END IF
      RETURN
      END
"""


def write_lapack_shim(cache_dir=None) -> Path:
    """Write the shim source (if needed) in the cache directory and return its path."""
    folder = (Path(cache_dir) if cache_dir else user_cache_dir()) / "shim"
    path = folder / LAPACK_SHIM_NAME
    try:
        if not path.exists() or path.read_text() != LAPACK_SHIM_SOURCE:
            folder.mkdir(parents=True, exist_ok=True)
            path.write_text(LAPACK_SHIM_SOURCE)
    except OSError:
        pass
    return path


def blas_mode(env=None, blas=None) -> str:
    """``"shim"`` (default: no BLAS/LAPACK needed) or ``"system"`` (link -llapack -lblas)."""
    env = os.environ if env is None else env
    mode = (blas or env.get(ENV_BLAS) or "shim").lower()
    if mode not in ("shim", "system"):
        raise ValueError(f"{ENV_BLAS} must be 'shim' or 'system', got {mode!r}")
    return mode


class FortranCompilerNotFound(RuntimeError):
    """No usable Fortran compiler could be located."""


@dataclass
class BuildResult:
    ok: bool
    exe: Path | None = None
    from_cache: bool = False
    key: str = ""
    failed: list = field(default_factory=list)  # [(source, stderr)]
    link_error: str = ""

    def report(self) -> str:
        """Human readable failure summary ('' when ok)."""
        if self.ok:
            return ""
        lines = [f"{len(self.failed)} file(s) failed to compile"] if self.failed else []
        lines += [f"  - {name}: {err.strip()}" for name, err in self.failed]
        if self.link_error:
            lines.append("link error: " + self.link_error.strip())
        return "\n".join(lines)


# ----------------------------------------------------------------------------
# platform helpers
# ----------------------------------------------------------------------------
def exe_suffix(system: str | None = None) -> str:
    system = system or platform.system()
    return ".exe" if system == "Windows" else ""


def exe_filename(name: str, system: str | None = None) -> str:
    """``name`` with the platform suffix, without doubling it (``cathy.exe`` stays as is)."""
    suffix = exe_suffix(system)
    return name if not suffix or name.lower().endswith(suffix) else name + suffix


def _conda_prefix(env) -> Path | None:
    prefix = env.get("CONDA_PREFIX")
    return Path(prefix) if prefix else None


def _conda_bin_dirs(prefix: Path) -> list[Path]:
    return [
        prefix / "bin",
        prefix / "Library" / "bin",
        prefix / "Library" / "mingw-w64" / "bin",
    ]


def find_fortran_compiler(env=None, which=shutil.which) -> str:
    """
    Locate a Fortran compiler.

    Order: ``$PYCATHY_FC``, ``$FC``, gfortran of the active conda environment
    (conda installs it with a platform prefix, e.g. ``x86_64-conda-linux-gnu-gfortran``),
    then ``gfortran`` / ``flang-new`` / ``flang`` on ``PATH``.
    """
    env = os.environ if env is None else env

    for var in (ENV_COMPILER, "FC"):
        value = env.get(var)
        if value:
            found = which(value) or (value if Path(value).exists() else None)
            if found:
                return str(found)

    prefix = _conda_prefix(env)
    if prefix is not None:
        for bindir in _conda_bin_dirs(prefix):
            if not bindir.is_dir():
                continue
            # gfortran first (linux/osx), then flang (what conda-forge ships on win-64)
            for pattern in ("*-gfortran", "*-gfortran.exe", "gfortran", "gfortran.exe",
                            "flang-new.exe", "flang-new", "flang.exe", "flang"):
                matches = sorted(bindir.glob(pattern))
                if matches:
                    return str(matches[0])

    for name in ("gfortran", "flang-new", "flang"):
        found = which(name)
        if found:
            return str(found)

    raise FortranCompilerNotFound(
        "No Fortran compiler found. Install one in your conda/mamba environment:\n"
        "    conda install -c conda-forge fortran-compiler\n"
        f"or set ${ENV_COMPILER} to the compiler path."
    )


def compiler_family(compiler: str | None) -> str:
    """``"flang"`` or ``"gfortran"`` (default) from the compiler file name."""
    name = Path(compiler).name.lower() if compiler else ""
    return "flang" if "flang" in name else "gfortran"


def default_compile_flags(system=None, machine=None, env=None, compiler=None) -> list[str]:
    """
    ``-O`` everywhere; ``-fPIC`` outside Windows; ``-mcmodel=large`` only on Linux x86-64
    (it is x86-64 specific; not supported on Apple Silicon, of no use on Windows).
    flang (conda-forge on Windows) only gets ``-O``.
    ``$PYCATHY_FFLAGS`` replaces the defaults entirely.
    """
    env = os.environ if env is None else env
    if env.get(ENV_FFLAGS):
        return env[ENV_FFLAGS].split()
    system = system or platform.system()
    machine = (machine or platform.machine()).lower()
    flags = ["-O"]
    if compiler_family(compiler) == "flang":
        return flags
    if system == "Linux" and machine in ("x86_64", "amd64"):
        flags.append("-mcmodel=large")
    if system != "Windows":
        flags.append("-fPIC")
    return flags


def runtime_link_flags(system=None, compiler=None) -> list[str]:
    """
    MinGW gfortran on Windows: link the Fortran runtime statically (no runtime DLL to
    find/relocate) and reserve a 256 MB stack (Windows default is 1 MB, Linux 8 MB: CATHY
    otherwise ends with 0xC00000FD, stack overflow).
    """
    system = system or platform.system()
    if system == "Windows" and compiler_family(compiler) == "gfortran":
        return ["-static-libgfortran", "-static-libgcc", "-Wl,--stack,268435456"]
    return []


def link_libraries(system=None, env=None, compiler=None) -> list[str]:
    """
    ``-llapack -lblas``, with the conda library directory (and rpath) when available.

    On Windows the conda libraries (``liblapack``/``libblas`` ``.lib`` import libraries, or
    ``.dll`` only) are passed by path: ``-llapack`` is not found by the MinGW linker.
    """
    env = os.environ if env is None else env
    system = system or platform.system()
    flags: list[str] = runtime_link_flags(system, compiler)
    prefix = _conda_prefix(env)
    if prefix is not None:
        libdir = prefix / "Library" / "lib" if system == "Windows" else prefix / "lib"
        if system == "Windows":
            # conda ships either import libraries (Library/lib/*.lib) or only the DLLs
            # (Library/bin/*.dll, which the MinGW linker can link against directly).
            for folder, ext in ((libdir, ".lib"), (prefix / "Library" / "bin", ".dll")):
                libs = [folder / ("liblapack" + ext), folder / ("libblas" + ext)]
                if all(lib.exists() for lib in libs):
                    return flags + [str(lib) for lib in libs]
        if libdir.is_dir():
            flags.append(f"-L{libdir}")
            if system != "Windows":
                flags.append(f"-Wl,-rpath,{libdir}")
    return flags + ["-llapack", "-lblas"]


def user_cache_dir(system=None, env=None, home=None) -> Path:
    """Build cache directory (``$PYCATHY_CACHE`` overrides the platform default)."""
    env = os.environ if env is None else env
    if env.get(ENV_CACHE):
        return Path(env[ENV_CACHE])
    system = system or platform.system()
    home = Path(home) if home else Path.home()
    if system == "Windows":
        base = Path(env.get("LOCALAPPDATA") or home / "AppData" / "Local")
    elif system == "Darwin":
        base = home / "Library" / "Caches"
    else:
        base = Path(env.get("XDG_CACHE_HOME") or home / ".cache")
    return base / "pyCATHY" / "build"


# ----------------------------------------------------------------------------
# build
# ----------------------------------------------------------------------------
def compiler_version(compiler: str, run=subprocess.run) -> str:
    try:
        p = run([compiler, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, timeout=60)
        return (p.stdout or "").splitlines()[0] if p.stdout else ""
    except (OSError, subprocess.SubprocessError):
        return ""


def build_key(src_dir, sources, compiler_id, compile_flags, link_flags) -> str:
    """sha256 over sources, headers (``*.H``), compiler identity and flags."""
    src_dir = Path(src_dir)
    h = hashlib.sha256()
    files = sorted({Path(s).name for s in sources} | {p.name for p in src_dir.glob("*.H")})
    for name in files:
        h.update(name.encode() + b"\0")
        h.update((src_dir / name).read_bytes())
        h.update(b"\0")
    for part in (compiler_id, " ".join(compile_flags), " ".join(link_flags)):
        h.update(part.encode() + b"\0")
    return h.hexdigest()


def _run(cmd, cwd):
    return subprocess.run(cmd, cwd=str(cwd), stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True, errors="replace")


def _make_executable(path: Path) -> None:
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


_MODULE_RE = re.compile(r"^\s*module\s+(?!procedure\b)(\w+)\s*$", re.I | re.M)
_USE_RE = re.compile(r"^\s*use\s+(\w+)", re.I | re.M)


def order_by_module_dependencies(src_dir, sources) -> list[str]:
    """
    Sort Fortran 90 sources so that a module is compiled before the files that ``use`` it
    (the former code relied on compiling the whole list twice). Ties and cycles keep the
    input order.
    """
    src_dir = Path(src_dir)
    provides, uses = {}, {}
    for name in sources:
        text = (src_dir / name).read_text(errors="replace")
        for mod in _MODULE_RE.findall(text):
            provides[mod.lower()] = name
        uses[name] = {m.lower() for m in _USE_RE.findall(text)}

    ordered, done, visiting = [], set(), set()

    def visit(name):
        if name in done or name in visiting:
            return
        visiting.add(name)
        for mod in sorted(uses[name]):
            dep = provides.get(mod)
            if dep and dep != name:
                visit(dep)
        visiting.discard(name)
        done.add(name)
        ordered.append(name)

    for name in sources:
        visit(name)
    return ordered


def build_preprocessor(src_dir, exe_name="pycppp", sources=PREPRO_SOURCES, compiler=None,
                       flags=None, use_cache=True, cache_dir=None) -> BuildResult:
    """Compile the CATHY preprocessor (single command: Fortran 90 module order matters)."""
    src_dir = Path(src_dir)
    compiler = compiler or find_fortran_compiler()
    flags = list(flags) if flags is not None else default_compile_flags(compiler=compiler)
    exe = src_dir / exe_filename(exe_name)
    cache = Path(cache_dir) if cache_dir else user_cache_dir()
    sources = order_by_module_dependencies(src_dir, list(sources))

    key = build_key(src_dir, sources, "prepro|" + compiler_version(compiler), flags, [])
    cached = cache / key / exe.name
    if use_cache and cached.exists():
        shutil.copy2(cached, exe)
        return BuildResult(True, exe, True, key)

    p = _run([compiler, *flags, "-o", str(exe), *sources], src_dir)
    if p.returncode != 0 or not exe.exists():
        return BuildResult(False, None, False, key, link_error=p.stderr)
    if use_cache:
        _store(exe, cached)
    return BuildResult(True, exe, False, key)


def build_cathy(src_dir, exe_name="cathy", compiler=None, flags=None, link_flags=None,
                workers=None, use_cache=True, force=False, cache_dir=None,
                blas=None) -> BuildResult:
    """
    Compile every ``*.f`` of ``src_dir`` (in parallel) and link ``exe_name``.

    ``blas="shim"`` (default, or ``$PYCATHY_BLAS``): DCOPY/DGETRF/DGETRS come from a small
    wrapper-provided file compiled next to the sources, nothing is linked.
    ``blas="system"``: link ``-llapack -lblas`` (conda libraries when available).

    Cache hit (same sources, ``*.H``, compiler, flags) -> the stored executable is copied
    and nothing is compiled. ``force=True`` ignores the cache (and refreshes it).
    """
    src_dir = Path(src_dir)
    compiler = compiler or find_fortran_compiler()
    flags = list(flags) if flags is not None else default_compile_flags(compiler=compiler)
    mode = blas_mode(blas=blas)
    if link_flags is None:
        link_flags = link_libraries(compiler=compiler) if mode == "system" else runtime_link_flags(compiler=compiler)
    exe = src_dir / exe_filename(exe_name)
    sources = sorted(p.name for p in src_dir.glob("*.f"))
    if not sources:
        return BuildResult(False, link_error=f"no *.f source file found in {src_dir}")

    cache = Path(cache_dir) if cache_dir else user_cache_dir()
    compile_list = list(sources)
    shim_id = ""
    if mode == "shim":
        compile_list.append(str(write_lapack_shim(cache)))
        shim_id = "shim-" + hashlib.sha256(LAPACK_SHIM_SOURCE.encode()).hexdigest()[:12]
    key = build_key(src_dir, sources, f"cathy|{mode}|{shim_id}|" + compiler_version(compiler), flags, link_flags)
    cached = cache / key / exe.name
    if use_cache and not force and cached.exists():
        shutil.copy2(cached, exe)
        return BuildResult(True, exe, True, key)

    for old in src_dir.glob("*.o"):
        old.unlink()

    def compile_one(name):
        obj = Path(name).with_suffix(".o").name
        p = _run([compiler, *flags, "-c", name, "-o", obj], src_dir)
        return Path(name).name, p.returncode, p.stderr

    workers = workers or min(32, os.cpu_count() or 1)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        outcomes = list(pool.map(compile_one, compile_list))
    failed = [(name, err) for name, rc, err in outcomes if rc != 0]
    if failed:
        return BuildResult(False, None, False, key, failed=failed)

    objects = [Path(n).with_suffix(".o").name for n in compile_list]
    p = _run([compiler, *flags, *objects, *link_flags, "-o", str(exe)], src_dir)
    if p.returncode != 0 or not exe.exists():
        return BuildResult(False, None, False, key, link_error=p.stderr)
    if use_cache:
        _store(exe, cached)
    return BuildResult(True, exe, False, key)


def _store(exe: Path, cached: Path) -> None:
    """Best effort: a read-only or full cache must never break a build."""
    tmp = cached.with_name(f".{cached.name}.{os.getpid()}.{id(exe)}.tmp")
    try:
        cached.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(exe, tmp)
        _make_executable(tmp)
        os.replace(tmp, cached)  # atomic: concurrent jobs (DA ensembles) never see a partial file
    except OSError:
        try:
            tmp.unlink()
        except OSError:
            pass
