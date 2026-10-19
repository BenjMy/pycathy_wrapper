"""Names that used to raise NameError."""
import numpy as np
import pytest

from pyCATHY import petro


def test_r2_matches_the_textbook_definition():
    y = np.array([1.0, 2.0, 3.0, 4.0])
    assert petro._r2_score(y, y) == pytest.approx(1.0)
    assert petro._r2_score(y, np.full(4, y.mean())) == pytest.approx(0.0)
    assert petro._r2_score(y, y + 1) == pytest.approx(1 - 4 / 5)


def test_print_stats_ensemble_ET_rejects_unknown_unit_instead_of_nameerror():
    xr = pytest.importorskip("xarray")
    from pyCATHY.importers.cathy_outputs import print_stats_ensemble_ET

    ds = xr.Dataset({"ACT. ETRA": (("time", "ensemble"), np.ones((2, 3)))},
                    coords={"time": np.array([0, 3600], dtype="timedelta64[s]").astype("timedelta64[ns]"),
                            "ensemble": [0, 1, 2]})
    with pytest.raises(ValueError, match="time_unit"):
        print_stats_ensemble_ET(ds, time_unit="weeks")
    print_stats_ensemble_ET(ds, time_unit="hours")  # runs
