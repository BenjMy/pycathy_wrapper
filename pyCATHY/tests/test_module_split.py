"""The code moved out of cathy_DA.py stays importable from there, and the new modules stay
free of cycles (they must not import cathy_DA)."""
import ast
import importlib
import pathlib

import pytest

DA_DIR = pathlib.Path(__file__).resolve().parents[1] / "DA"

MOVED = {
    "analysis": ["AnalysisResult", "AnalysisConfig", "DAConfig", "ANALYSES", "run_analysis",
                 "_from_enkf", "_enkf", "_enkf_localized_with_inflation", "_particle_filter",
                 "_LEGACY_DA_KEYS", "_ENKF_DIAGNOSTICS"],
    "param_tables": ["_SOIL_PARM_COLUMNS", "_assign_soil_value", "_nonphysical_values",
                     "_has_name_token", "_set_feddes_value", "_valid_zone_values",
                     "_check_positive_pressure_heads", "_update_key"],
    "runtime": ["_Tee", "_tee_run_output", "_working_directory", "_restores_cwd", "_configure_logging"],
    "performance_tables": ["_compute_rmse", "_compute_diff_matrix", "_append_to_performance_df",
                           "_save_performance_to_parquet", "_load_performance_from_parquet"],
}


@pytest.mark.parametrize("module, names", MOVED.items())
def test_moved_names_are_the_same_objects_in_cathy_DA(module, names):
    new = importlib.import_module(f"pyCATHY.DA.{module}")
    old = importlib.import_module("pyCATHY.DA.cathy_DA")
    for name in names:
        assert getattr(old, name) is getattr(new, name), name


@pytest.mark.parametrize("module", MOVED)
def test_new_modules_do_not_import_cathy_DA(module):
    tree = ast.parse((DA_DIR / f"{module}.py").read_text())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
            imported |= {f"{node.module}.{a.name}" for a in node.names}
        elif isinstance(node, ast.Import):
            imported |= {a.name for a in node.names}
    assert not any(name.endswith("cathy_DA") for name in imported), imported


def test_scenario_module_does_not_need_the_whole_DA_class():
    """scenario.py reads the soil-table helpers from param_tables, not from cathy_DA."""
    text = (DA_DIR / "scenario.py").read_text()
    assert "from pyCATHY.DA.param_tables import" in text
    assert "from pyCATHY.DA.cathy_DA import _SOIL_PARM_COLUMNS" not in text
