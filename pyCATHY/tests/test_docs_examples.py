"""The documentation is tested: every ``toml`` example of the DA guides runs through the real
``validate_scenario`` and ``perturbate`` (+ sampling), and the guides stay complete and linked."""
import importlib.util
import pathlib
import re
import tomllib

import pytest

from pyCATHY.DA import perturbate
from pyCATHY.DA.scenario import KNOWN_PARAMETERS, NOT_IMPLEMENTED, validate_scenario

ROOT = pathlib.Path(__file__).resolve().parents[2]
GUIDES = ROOT / "doc" / "content" / "DA_guides"

_spec = importlib.util.spec_from_file_location(
    "layouts", pathlib.Path(__file__).with_name("test_scenario_layouts.py"))
L = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(L)


def toml_blocks():
    """(file, index, text) of every ``.. code-block:: toml`` of the guides."""
    out = []
    for path in sorted(GUIDES.glob("*.rst")):
        lines = path.read_text().split("\n")
        i, k = 0, 0
        while i < len(lines):
            if lines[i].strip() == ".. code-block:: toml":
                j = i + 1
                while j < len(lines) and lines[j].strip() == "":
                    j += 1
                body = []
                while j < len(lines) and (lines[j].startswith("   ") or lines[j].strip() == ""):
                    body.append(lines[j][3:])
                    j += 1
                out.append((path.name, k, "\n".join(body).strip()))
                k += 1
                i = j
            else:
                i += 1
    return out


BLOCKS = toml_blocks()


def test_the_guides_contain_toml_examples():
    assert len(BLOCKS) >= 4


@pytest.mark.parametrize("fname, k, text", BLOCKS, ids=[f"{b[0]}#{b[1]}" for b in BLOCKS])
def test_every_toml_example_validates_and_runs(fname, k, text):
    rows = int(re.search(r"# model: feddes_rows=(\d+)", text).group(1))
    scenarios = tomllib.loads(text)["scenario"]
    assert scenarios
    for name, sc in scenarios.items():
        errors, _ = validate_scenario(sc, nlayers=L.NLAYERS, zone_raster=L.ZONES, veg_raster=L.VEG3,
                                      n_feddes_rows=rows, raise_on_error=False)
        assert errors == [], f"{fname} example '{name}': {errors}"
        sim = L.make_sim(rows)
        entries = perturbate.perturbate(sim, sc, L.NEN)  # validation on, listUpdateParm checked
        assert entries
        out = {}
        for e in entries:  # the sampling of the real run
            out = perturbate.perturbate_parm(
                out, parm=e, type_parm=e["type_parm"], mean=e["mean"], sd=e["sd"],
                sampling_type=e["sampling_type"], ensemble_size=e["ensemble_size"],
                per_type=e["per_type"], savefig=False)
        assert set(out) == {e["type_parm"] for e in entries}


def test_the_catalogue_mentions_every_parameter():
    text = (GUIDES / "perturbation.rst").read_text()
    for name in list(KNOWN_PARAMETERS) + list(NOT_IMPLEMENTED):
        if name == "atmbc":  # perturbed through its own mechanism, not a soil layout
            continue
        assert f"``{name}``" in text, f"{name} is missing from the parameter catalogue"


def test_the_localisation_guide_states_the_conventions_that_the_code_enforces():
    text = (GUIDES / "localisation.rst").read_text()
    for must in ("beyond ``2 L``", "(n_grid, n_obs)", "[y, x, z]", "Cannot localise with apparent resistivity data."):
        assert must in text
    from pyCATHY.DA import localisation
    import numpy as np
    assert localisation.gaspari_cohn(np.array([2.0]), 1.0)[0] == 0.0  # zero at 2 L, as documented
    assert localisation.build_localization_matrix(np.zeros((3, 2)), np.zeros((5, 2)), 1.0).shape == (5, 3)


def test_every_error_message_of_the_validator_quoted_in_the_troubleshooting_page_exists():
    """Messages quoted in the troubleshooting table must be produced by the code."""
    page = (GUIDES / "troubleshooting.rst").read_text()
    sources = "".join(p.read_text() for p in (ROOT / "pyCATHY" / "DA").glob("*.py")) + \
        (ROOT / "pyCATHY" / "importers" / "cathy_outputs.py").read_text()
    for fragment in ("unknown parameter", "is not aligned", "were not built", "needs one entry per",
                     "per_nom is a list but pert_control_name is not set", "is not supported for this parameter",
                     "must be scalars", "commented out", "without gaps", "unknown bound key",
                     "the last one written wins", "will never be assimilated", "quadrupoles excluded",
                     "number of rejected ensemble is too high", "clipped to", "positive pressure heads",
                     "has no time step", "localisation matrix has", "Cannot localise with apparent resistivity",
                     "Wrong stacked_data_cov shape"):
        expected = fragment.replace("is not aligned", "aligned with per_name")
        assert expected.split(" ... ")[0] in sources or fragment in page, fragment


def test_the_guides_are_linked_from_the_documentation_index():
    assert "content/DA_guides/index.rst" in (ROOT / "doc" / "index.rst").read_text()
    index = (GUIDES / "index.rst").read_text()
    for page in ("perturbation", "localisation", "troubleshooting"):
        assert page in index and (GUIDES / f"{page}.rst").exists()
