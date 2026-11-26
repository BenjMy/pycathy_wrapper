"""Tests of pyCATHY.fortran_build (compiler lookup, flags, cache key, real mini build)."""

import shutil
from pathlib import Path

import pytest

from pyCATHY import fortran_build as fb


def _fake_exe(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\n")
    return path


# ----------------------------------------------------------------------------
# compiler lookup
# ----------------------------------------------------------------------------
def test_pycathy_fc_has_priority(tmp_path):
    fc = _fake_exe(tmp_path / "my-fc")
    env = {fb.ENV_COMPILER: str(fc), "FC": "other"}
    assert fb.find_fortran_compiler(env=env, which=lambda n: None) == str(fc)


def test_fc_env_resolved_through_which():
    env = {"FC": "x86_64-conda-linux-gnu-gfortran"}
    found = fb.find_fortran_compiler(env=env, which=lambda n: "/opt/" + n)
    assert found == "/opt/x86_64-conda-linux-gnu-gfortran"


def test_conda_prefixed_gfortran_is_found(tmp_path):
    exe = _fake_exe(tmp_path / "bin" / "x86_64-conda-linux-gnu-gfortran")
    env = {"CONDA_PREFIX": str(tmp_path)}
    assert fb.find_fortran_compiler(env=env, which=lambda n: None) == str(exe)


def test_windows_conda_layout(tmp_path):
    exe = _fake_exe(tmp_path / "Library" / "bin" / "gfortran.exe")
    env = {"CONDA_PREFIX": str(tmp_path)}
    assert fb.find_fortran_compiler(env=env, which=lambda n: None) == str(exe)


def test_path_fallback_and_flang():
    which = lambda n: "/usr/bin/flang-new" if n == "flang-new" else None  # noqa: E731
    assert fb.find_fortran_compiler(env={}, which=which) == "/usr/bin/flang-new"


def test_missing_compiler_gives_conda_hint():
    with pytest.raises(fb.FortranCompilerNotFound) as exc:
        fb.find_fortran_compiler(env={}, which=lambda n: None)
    assert "conda install -c conda-forge fortran-compiler" in str(exc.value)


# ----------------------------------------------------------------------------
# flags / names / cache dir
# ----------------------------------------------------------------------------
def test_default_flags_per_platform():
    assert fb.default_compile_flags("Linux", "x86_64", env={}) == ["-O", "-mcmodel=large", "-fPIC"]
    assert fb.default_compile_flags("Linux", "aarch64", env={}) == ["-O", "-fPIC"]
    assert fb.default_compile_flags("Darwin", "arm64", env={}) == ["-O", "-fPIC"]
    assert fb.default_compile_flags("Windows", "AMD64", env={}) == ["-O"]


def test_fflags_env_replaces_defaults():
    assert fb.default_compile_flags("Linux", "x86_64", env={fb.ENV_FFLAGS: "-O3 -g"}) == ["-O3", "-g"]


def test_exe_names():
    assert fb.exe_suffix("Windows") == ".exe"
    assert fb.exe_suffix("Linux") == ""
    assert fb.exe_filename("cathy", "Windows") == "cathy.exe"
    assert fb.exe_filename("cathy.exe", "Windows") == "cathy.exe"  # never doubled
    assert fb.exe_filename("cathy", "Darwin") == "cathy"


def test_link_libraries_with_and_without_conda(tmp_path):
    assert fb.link_libraries("Linux", env={}) == ["-llapack", "-lblas"]
    (tmp_path / "lib").mkdir()
    flags = fb.link_libraries("Linux", env={"CONDA_PREFIX": str(tmp_path)})
    assert flags[:2] == [f"-L{tmp_path / 'lib'}", f"-Wl,-rpath,{tmp_path / 'lib'}"]
    assert flags[-2:] == ["-llapack", "-lblas"]
    (tmp_path / "Library" / "lib").mkdir(parents=True)
    win = fb.link_libraries("Windows", env={"CONDA_PREFIX": str(tmp_path)})
    assert win == [f"-L{tmp_path / 'Library' / 'lib'}", "-llapack", "-lblas"]


def test_user_cache_dir(tmp_path):
    assert fb.user_cache_dir("Linux", env={fb.ENV_CACHE: str(tmp_path)}) == tmp_path
    assert fb.user_cache_dir("Linux", env={}, home=tmp_path) == tmp_path / ".cache/pyCATHY/build"
    assert fb.user_cache_dir("Darwin", env={}, home=tmp_path) == tmp_path / "Library/Caches/pyCATHY/build"
    assert fb.user_cache_dir("Windows", env={"LOCALAPPDATA": str(tmp_path)}) == tmp_path / "pyCATHY/build"


# ----------------------------------------------------------------------------
# module ordering and cache key
# ----------------------------------------------------------------------------
def test_module_dependency_order(tmp_path):
    (tmp_path / "main.f90").write_text("program p\nuse b_mod\nend program p\n")
    (tmp_path / "b.f90").write_text("module b_mod\nuse a_mod\nend module b_mod\n")
    (tmp_path / "a.f90").write_text("module a_mod\nend module a_mod\n")
    order = fb.order_by_module_dependencies(tmp_path, ["main.f90", "b.f90", "a.f90"])
    assert order == ["a.f90", "b.f90", "main.f90"]


def test_module_order_keeps_input_order_without_dependencies(tmp_path):
    for n in ("x.f90", "y.f90", "z.f90"):
        (tmp_path / n).write_text("subroutine s\nend\n")
    names = ["z.f90", "x.f90", "y.f90"]
    assert fb.order_by_module_dependencies(tmp_path, names) == names


def test_module_order_cycle_does_not_hang(tmp_path):
    (tmp_path / "a.f90").write_text("module am\nuse bm\nend module am\n")
    (tmp_path / "b.f90").write_text("module bm\nuse am\nend module bm\n")
    assert sorted(fb.order_by_module_dependencies(tmp_path, ["a.f90", "b.f90"])) == ["a.f90", "b.f90"]


def test_build_key_sensitivity(tmp_path):
    (tmp_path / "a.f").write_text("      end\n")
    (tmp_path / "CATHY.H").write_text("      PARAMETER (MAXVEG=1)\n")
    base = fb.build_key(tmp_path, ["a.f"], "gfortran 13", ["-O"], ["-llapack"])
    assert base == fb.build_key(tmp_path, ["a.f"], "gfortran 13", ["-O"], ["-llapack"])
    assert base != fb.build_key(tmp_path, ["a.f"], "gfortran 14", ["-O"], ["-llapack"])
    assert base != fb.build_key(tmp_path, ["a.f"], "gfortran 13", ["-O2"], ["-llapack"])
    (tmp_path / "CATHY.H").write_text("      PARAMETER (MAXVEG=2)\n")
    assert base != fb.build_key(tmp_path, ["a.f"], "gfortran 13", ["-O"], ["-llapack"])


# ----------------------------------------------------------------------------
# real build (skipped when no compiler is installed)
# ----------------------------------------------------------------------------
def _compiler_or_skip():
    try:
        return fb.find_fortran_compiler()
    except fb.FortranCompilerNotFound:
        pytest.skip("no Fortran compiler available")


def _mini_project(root: Path) -> Path:
    src = root / "src"
    src.mkdir()
    (src / "CATHY.H").write_text("      INTEGER MAXVEG\n      PARAMETER (MAXVEG=1)\n")
    (src / "main.f").write_text(
        "      PROGRAM MAIN\n      INCLUDE 'CATHY.H'\n      CALL SUB(MAXVEG)\n      END\n"
    )
    (src / "sub.f").write_text(
        "      SUBROUTINE SUB(N)\n      INTEGER N\n      WRITE(*,*) N\n      END\n"
    )
    return src


def test_real_build_cache_and_rebuild(tmp_path):
    compiler = _compiler_or_skip()
    src = _mini_project(tmp_path)
    cache = tmp_path / "cache"
    kw = dict(compiler=compiler, flags=["-O"], link_flags=[], cache_dir=cache)

    first = fb.build_cathy(src, "prog", **kw)
    assert first.ok and not first.from_cache and first.exe.exists()

    first.exe.unlink()
    again = fb.build_cathy(src, "prog", **kw)
    assert again.ok and again.from_cache and again.exe.exists()
    assert again.key == first.key

    forced = fb.build_cathy(src, "prog", force=True, **kw)
    assert forced.ok and not forced.from_cache

    (src / "CATHY.H").write_text("      INTEGER MAXVEG\n      PARAMETER (MAXVEG=2)\n")
    changed = fb.build_cathy(src, "prog", **kw)
    assert changed.ok and not changed.from_cache and changed.key != first.key


def test_real_build_failure_is_reported_not_swallowed(tmp_path):
    compiler = _compiler_or_skip()
    src = _mini_project(tmp_path)
    (src / "broken.f").write_text("      THIS IS NOT FORTRAN (\n")
    res = fb.build_cathy(src, "prog", compiler=compiler, flags=["-O"], link_flags=[],
                         cache_dir=tmp_path / "cache")
    assert not res.ok and res.exe is None
    assert [name for name, _ in res.failed] == ["broken.f"]
    assert "broken.f" in res.report()


def test_real_build_without_sources(tmp_path):
    res = fb.build_cathy(tmp_path, "prog", compiler="gfortran", cache_dir=tmp_path / "c")
    assert not res.ok and "no *.f source" in res.report()


def test_real_preprocessor_module_order(tmp_path):
    compiler = _compiler_or_skip()
    (tmp_path / "main.f90").write_text("program p\nuse b_mod\ncall hello()\nend program p\n")
    (tmp_path / "b.f90").write_text(
        "module b_mod\ncontains\nsubroutine hello()\nprint *, 'ok'\nend subroutine\nend module b_mod\n"
    )
    res = fb.build_preprocessor(tmp_path, "pp", sources=["main.f90", "b.f90"], compiler=compiler,
                                flags=["-O"], cache_dir=tmp_path / "cache")
    assert res.ok, res.report()
    assert res.exe.exists()
    assert shutil.which(str(res.exe)) or res.exe.is_file()
