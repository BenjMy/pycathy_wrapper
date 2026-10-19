"""per_bounds parsing: every accepted form works, malformed input raises (never silent)."""
import pytest

from pyCATHY.DA.perturbate import check4bounds, parse_bounds


def _sc(bounds, name="ZROOT"):
    return {"per_name": [name], "per_bounds": [bounds]}


@pytest.mark.parametrize(
    "spec, expected",
    [
        ({"min": 0, "max": 1.5}, (0, 1.5)),
        ([0, 1.5], (0, 1.5)),
        ({"min": 0}, (0, None)),
        ({"max": -0.05}, (None, -0.05)),
        ("None", (None, None)),
        (None, (None, None)),
    ],
)
def test_accepted_forms(spec, expected):
    assert parse_bounds(spec) == expected


@pytest.mark.parametrize(
    "spec",
    [{"mn": 0, "mx": 1.5}, {"min": 2, "max": 1}, {"min": 1, "max": 1}, [0, 1, 2], "abc", {"min": "a"}],
)
def test_malformed_raises(spec):
    with pytest.raises(ValueError):
        parse_bounds(spec)


def test_check4bounds_keeps_defaults_for_open_sides():
    assert check4bounds(_sc({"min": 0}), 0, -1, 9) == (0, 9)
    assert check4bounds(_sc("None"), 0, 0, 5) == (0, 5)
    assert check4bounds({"per_name": ["a"]}, 0, 0, 5) == (0, 5)  # per_bounds is optional


def test_shared_bounds_apply_to_every_layer_or_class():
    """A single {min, max} used to be dropped silently when het_size > 1."""
    sc = _sc({"min": 0, "max": 1.5})
    for nb in range(3):
        assert check4bounds(sc, 0, None, None, het_size=3, het_nb=nb) == (0, 1.5)
    sc = _sc([0, 1.5])  # a [min, max] pair is shared too, even when het_size == 2
    assert check4bounds(sc, 0, None, None, het_size=2, het_nb=1) == (0, 1.5)


def test_per_entry_bounds():
    sc = _sc([{"min": 0, "max": 1}, "None", {"min": 5, "max": 6}])
    assert check4bounds(sc, 0, None, None, het_size=3, het_nb=0) == (0, 1)
    assert check4bounds(sc, 0, 7, 8, het_size=3, het_nb=1) == (7, 8)
    assert check4bounds(sc, 0, None, None, het_size=3, het_nb=2) == (5, 6)
    with pytest.raises(ValueError):  # wrong number of entries
        check4bounds(sc, 0, None, None, het_size=2, het_nb=0)


def test_typo_in_scenario_is_an_error_not_silent():
    with pytest.raises(ValueError, match="ZROOT"):
        check4bounds(_sc({"mn": 0, "mx": 1.5}), 0, None, None)


def test_missing_entry_leaves_unbounded():
    sc = {"per_name": ["a", "b"], "per_bounds": [{"min": 0, "max": 1}]}
    assert check4bounds(sc, 1, None, None) == (None, None)
