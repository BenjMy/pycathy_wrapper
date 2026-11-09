"""
Fortran-free tests of the version dependent layout of a CATHY project
(version="withIrr" vs the default one): readers, line groupings, file names, CATHY.H patching.

Fixtures in data/withIrr/ are copies of files of the CATHY_src_WithIrrigation dataset.
"""
import os
from types import SimpleNamespace

import pytest

from pyCATHY.cathy_tools import CATHY, _patch_fortran_parameters
from pyCATHY.config_manager import VersionConfigManager
from pyCATHY.importers import cathy_inputs as in_CT

HERE = os.path.dirname(os.path.abspath(__file__))
WITHIRR = os.path.join(HERE, "data", "withIrr")
DEFAULT = os.path.join(HERE, "data", "mapping", "input")


def layout(version):
    return VersionConfigManager("cathy", "dir_layout", version)


def stub(version):
    """Minimal object carrying what _input_path/_prepro_path need."""
    cfg = layout(version)
    ns = SimpleNamespace(
        workdir="/w/",
        project_name="prj",
        input_dirname="input",
        prepro_input_dirname=cfg.get("prepro_dirname", "prepro"),
        other_dem_dirname=cfg.get("other_dem_dirname"),
        _other_dem_files=list(cfg.get("other_dem_files", [])),
        _input_file_names=dict(cfg.get("input_file_names", {})),
        _prepro_file_names=dict(cfg.get("prepro_file_names", {})),
        _fnames_entries=CATHY._fnames_entries,
    )
    ns.__dict__["_resolve_from_fnames"] = lambda pd, cand: CATHY._resolve_from_fnames(ns, pd, cand)
    return ns


# ---------------------------------------------------------------- dem_parameters
def test_dem_parameters_withIrr_header_does_not_overwrite_nstr():
    dem = in_CT.read_dem_parameters(os.path.join(WITHIRR, "other-dem", "dem_parameters"))
    assert dem["nzone"] == 750
    assert dem["nstr"] == 9
    assert dem["n1"] == 25
    assert len(dem["zratio(i),i=1,nstr"]) == dem["nstr"]
    assert dem["delta_x"] == dem["delta_y"] == 2.0


def test_dem_parameters_default_layout_unchanged():
    dem = in_CT.read_dem_parameters(os.path.join(DEFAULT, "dem_parameters"))
    assert len(dem["zratio(i),i=1,nstr"]) == dem["nstr"]


@pytest.mark.parametrize("version,path", [
    ("1.0.0", os.path.join(DEFAULT, "dem_parameters")),
    ("withIrr", os.path.join(WITHIRR, "other-dem", "dem_parameters")),
])
def test_dem_parameters_fmt_covers_every_entry(version, path):
    dem = in_CT.read_dem_parameters(path)
    assert sum(layout(version).get("dem_parameters_fmt")) == len(dem)


# ------------------------------------------------------------------------- parm
@pytest.mark.parametrize("version,path", [
    ("1.0.0", os.path.join(DEFAULT, "parm")),
    ("withIrr", os.path.join(WITHIRR, "input", "parm")),
])
def test_parm_fmt_covers_every_entry(version, path):
    parm = in_CT.read_parm(path)
    assert sum(layout(version).get("parm_fmt")) == len(parm)


def test_parm_withIrr_has_no_ncout():
    parm = in_CT.read_parm(os.path.join(WITHIRR, "input", "parm"))
    assert "NCOUT" not in parm
    assert "TRAFLAG" in parm


# ------------------------------------------------------------------- file paths
def test_default_layout_paths():
    s = stub("1.0.0")
    assert CATHY._input_path(s, "atmbc") == "/w/prj/input/atmbc"
    assert CATHY._input_path(s, "root_map") == "/w/prj/input/root_map"
    assert CATHY._input_path(s, "dem_parameters") == "/w/prj/input/dem_parameters"
    assert CATHY._prepro_path(s, "zone") == "/w/prj/prepro/zone"
    assert CATHY._prepro_path(s) == "/w/prj/prepro"


def test_withIrr_layout_paths_follow_cathy_fnames_of_the_sources():
    s = stub("withIrr")
    assert CATHY._input_path(s, "atmbc") == "/w/prj/input/atmbc-flow"
    assert CATHY._input_path(s, "ic") == "/w/prj/input/ic-flow"
    assert CATHY._input_path(s, "parm") == "/w/prj/input/parm"
    assert CATHY._input_path(s, "soil") == "/w/prj/input/soil"
    # DEM related files stay where the solver reads them (cathy.fnames is never modified)
    assert CATHY._input_path(s, "dem_parameters") == "/w/prj/other-dem/dem_parameters"
    assert CATHY._input_path(s, "root_map") == "/w/prj/other-dem/veg-type_map"
    assert CATHY._prepro_path(s, "zone") == "/w/prj/prepro-dem/zone_map"
    assert CATHY._prepro_path(s, "hap.in") == "/w/prj/prepro-dem/hap.in"


# -------------------------------------------------------------------- CATHY.H
def test_patch_parameters_replaces_only_requested_values():
    src = open(os.path.join(WITHIRR, "CATHY.H")).readlines()
    new = _patch_fortran_parameters(src, {"ROWMAX": 30, "MAXTRM": 202500.0, "MAXPRT": 9})
    changed = [(a.strip(), b.strip()) for a, b in zip(src, new) if a != b]
    assert len(new) == len(src)
    assert len(changed) == 3
    text = "".join(new)
    assert "ROWMAX=30,COLMAX=10,DEMRES=1" in text
    assert "MAXTRM=202500)" in text
    assert "NRMAX=1,MAXPRT=9,MAXVP=40,MAXQOUT=1" in text
    # derived sizes keep their expressions, declarations and comments are untouched
    assert "NODMAX=(ROWMAX / DEMRES + 1)*(COLMAX / DEMRES + 1)" in text
    assert "cxcx  PARAMETER (NIAUXMAX=NFACEMAX + MAXTRM + 1)" in text


def test_patch_parameters_is_idempotent():
    src = open(os.path.join(WITHIRR, "CATHY.H")).readlines()
    once = _patch_fortran_parameters(src, {"MAXSTR": 12})
    assert _patch_fortran_parameters(once, {"MAXSTR": 12}) == once


def test_patch_parameters_rejects_lines_beyond_column_72():
    with pytest.raises(ValueError):
        _patch_fortran_parameters(
            ["      PARAMETER (AAAAAAAAAA=1,BBBBBBBBBB=2,CCCCCCCCCC=3,DDDDDDDDDD=4)\n"],
            {"AAAAAAAAAA": 123456789012345678901234567890},
        )


# ------------------------------------------------------------ mesh body files
def test_withIrr_declares_mesh_body_files_default_does_not():
    assert layout("withIrr").get("mesh_body_files") == {
        "grid3d": "fort.1199",
        "grid2d.exp": "fort.1198",
    }
    assert layout("1.0.0").get("mesh_body_files") == {}


def test_assemble_mesh_files_rebuilds_grid3d_from_header_and_body(tmp_path):
    (tmp_path / "output").mkdir()
    (tmp_path / "output" / "grid3d").write_text("      4      8      2\n")
    (tmp_path / "output" / "grid2d.exp").write_text("")
    (tmp_path / "fort.1199").write_text("1 2 3 4 1\n5 6 7 8 1\n")
    (tmp_path / "fort.1198").write_text("2 4 1 0\n")
    s = stub("withIrr")
    s.workdir, s.project_name, s.version = str(tmp_path), "", "withIrr"
    CATHY._assemble_mesh_files(s)
    assert (tmp_path / "output" / "grid3d").read_text() == "      4      8      2\n1 2 3 4 1\n5 6 7 8 1\n"
    assert (tmp_path / "output" / "grid2d.exp").read_text() == "2 4 1 0\n"
    assert not (tmp_path / "fort.1199").exists()


def test_assemble_mesh_files_is_a_noop_for_default_version(tmp_path):
    (tmp_path / "output").mkdir()
    (tmp_path / "output" / "grid3d").write_text("header\nbody\n")
    s = stub("1.0.0")
    s.workdir, s.project_name, s.version = str(tmp_path), "", "1.0.0"
    CATHY._assemble_mesh_files(s)
    assert (tmp_path / "output" / "grid3d").read_text() == "header\nbody\n"


# ------------------------------------------------------------------ vtk names
def test_flow_vtk_names_follow_the_files_present(tmp_path):
    from pyCATHY.cathy_utils import flow_vtk_filename, flow_vtk_glob

    d = str(tmp_path)
    # default sources
    assert flow_vtk_filename(d, 1) == "101.vtk"
    assert flow_vtk_filename(d, 12) == "112.vtk"
    assert flow_vtk_filename(d, 105) == "205.vtk"
    assert flow_vtk_glob(d) == ("1*.vtk", "100.vtk")
    # withIrr sources
    (tmp_path / "flow2000.vtk").write_text("")
    assert flow_vtk_filename(d, 1) == "flow2001.vtk"
    assert flow_vtk_glob(d) == ("flow*.vtk", "flow2000.vtk")


# ------------------------------------------------- rasters of a non-square DEM
@pytest.mark.parametrize("method,arg", [("update_dem", "DEM"), ("update_lakes_map", "lakes_map")])
def test_raster_header_rows_cols_match_a_non_square_dem(tmp_path, method, arg):
    import numpy as np

    s = stub("1.0.0")
    s.workdir, s.project_name = str(tmp_path), ""
    (tmp_path / "prepro").mkdir()
    s.hapin = {"N": 10, "M": 75, "xllcorner": 0.0, "yllcorner": 0.0}  # 75 rows, 10 columns
    s._prepro_path = lambda name="": CATHY._prepro_path(s, name)
    s.version = "1.0.0"
    s._dem_nodata_to_file = lambda DEM: CATHY._dem_nodata_to_file(s, DEM)
    getattr(CATHY, method)(s, **{arg: np.ones((75, 10))})
    name = "dem" if method == "update_dem" else "lakes_map"
    lines = (tmp_path / "prepro" / name).read_text().splitlines()
    assert lines[4].split()[1] == "75" and lines[5].split()[1] == "10"
    assert len(lines) - 6 == 75


# ------------------------------------------------------------------- DEMRES
def test_demres_is_not_the_cell_size():
    """DEMRES of CATHY.H is a coarsening factor (1), it must not follow the grid spacing."""
    import inspect

    src = inspect.getsource(CATHY.update_prepo_inputs)
    assert "DEMRES=min" not in src


# --------------------------------------------------- hand edited parm (Venezia)
def test_read_parm_tolerates_leftover_numbers_glued_to_a_label():
    # real file: '... DELTAT DTMIN  DTMAX  TMAX608400.086400.0' (Fortran ignores the leftovers)
    parm = in_CT.read_parm(os.path.join(HERE, "data", "hand_edited", "parm_glued_label"))
    assert parm["TMAX"] == 2592000.0
    assert parm["DTMAX"] == 900.0
    assert not any(k.startswith("TMAX") and k != "TMAX" for k in parm)
    assert len(parm["(TIMPRT(I),I=1,NPRT)"]) == parm["NPRT"]


# --------------------------------------------- nodata convention of the dem file
@pytest.mark.parametrize("version,expected", [("1.0.0", -9999.0), ("withIrr", 0.0)])
def test_update_dem_translates_nodata_for_the_processor(tmp_path, version, expected):
    import numpy as np

    s = stub(version)
    s.workdir, s.project_name, s.version = str(tmp_path), "", version
    (tmp_path / layout(version).get("prepro_dirname", "prepro")).mkdir()
    s.hapin = {"N": 3, "M": 2, "xllcorner": 0.0, "yllcorner": 0.0}
    s._prepro_path = lambda name="": CATHY._prepro_path(s, name)
    s._dem_nodata_to_file = lambda DEM: CATHY._dem_nodata_to_file(s, DEM)
    DEM = np.array([[1.0, -9999.0, 2.0], [3.0, 4.0, -9999.0]])
    CATHY.update_dem(s, DEM)
    d = tmp_path / layout(version).get("prepro_dirname", "prepro")
    dem = np.loadtxt(d / "dem", skiprows=6)
    dtm = np.loadtxt(d / "dtm_13.val")
    assert dem[0, 1] == expected and dem[1, 2] == expected   # processor file
    assert dtm[0, 1] == -9999.0 and dtm[1, 2] == -9999.0     # preprocessor file keeps -9999
    assert dem[0, 0] == 1.0 and dem[1, 1] == 4.0


# --------------------------------------------- locations read from cathy.fnames
FNAMES = os.path.join(HERE, "data", "fnames")


def project_with(tmp_path, fnames_file, version):
    import shutil

    shutil.copy(os.path.join(FNAMES, fnames_file), tmp_path / "cathy.fnames")
    s = stub(version)
    s.workdir, s.project_name = str(tmp_path), ""
    return s


def rel(s, path):
    return os.path.relpath(path, s.workdir)


def test_fnames_of_the_irrigation_sources_of_2024_keep_the_old_input_names(tmp_path):
    # Salerno's own src: new folders (other-dem, prepro-dem) but atmbc / ic / transp names
    s = project_with(tmp_path, "salerno_2024_irrigation.fnames", "withIrr")
    assert rel(s, CATHY._input_path(s, "atmbc")) == "input/atmbc"
    assert rel(s, CATHY._input_path(s, "ic")) == "input/ic"
    assert rel(s, CATHY._input_path(s, "dem_parameters")) == "other-dem/dem_parameters"
    assert rel(s, CATHY._input_path(s, "root_map")) == "other-dem/veg-type_map"
    assert rel(s, CATHY._prepro_path(s, "zone")) == "prepro-dem/zone_map"


def test_fnames_of_the_irrigation_sources_of_2025_use_flow_names(tmp_path):
    # code/ of CATHY_src_WithIrrigation
    s = project_with(tmp_path, "withIrr_code_2025.fnames", "withIrr")
    assert rel(s, CATHY._input_path(s, "atmbc")) == "input/atmbc-flow"
    assert rel(s, CATHY._input_path(s, "ic")) == "input/ic-flow"
    assert rel(s, CATHY._input_path(s, "dem_parameters")) == "other-dem/dem_parameters"
    assert rel(s, CATHY._input_path(s, "root_map")) == "other-dem/veg-type_map"
    assert rel(s, CATHY._prepro_path(s, "zone")) == "prepro-dem/zone_map"


def test_fnames_of_the_old_sources_venezia(tmp_path):
    s = project_with(tmp_path, "venezia_old.fnames", "1.0.0")
    assert rel(s, CATHY._input_path(s, "atmbc")) == "input/atmbc"
    assert rel(s, CATHY._input_path(s, "root_map")) == "input/root_map"
    assert rel(s, CATHY._input_path(s, "dem_parameters")) == "input/dem_parameters"
    assert rel(s, CATHY._prepro_path(s, "zone")) == "prepro/zone"


def test_without_fnames_the_table_of_the_version_is_used(tmp_path):
    s = stub("withIrr")
    s.workdir, s.project_name = str(tmp_path), ""
    assert rel(s, CATHY._input_path(s, "atmbc")) == "input/atmbc-flow"
