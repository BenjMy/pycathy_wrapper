"""Reader/writer of the irrigation_map input (format of irr_read_one.f), real Salerno excerpts."""
import os

import numpy as np
import pytest

from pyCATHY.irrigation import read_irrigation_map, write_irrigation_map

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "irrigation")
NROW, NCOL = 72, 74


def test_read_header_and_blocks_of_the_real_file():
    m = read_irrigation_map(os.path.join(DATA, "salerno_irrigation_map_first6"), NROW, NCOL)
    assert (m["irr_flag"], m["irr_sw"], m["irr_hs"]) == (1, 0.35, 0)
    assert m["times"].tolist() == [0.0, 32400.0, 36000.0, 118800.0, 122400.0, 205200.0]
    assert m["rates"].shape == (6, NROW, NCOL)


def test_read_irrigated_block_has_greenhouse_rates():
    m = read_irrigation_map(os.path.join(DATA, "salerno_irrigation_map_first_irrigated_block"), NROW, NCOL)
    r = m["rates"][0]
    assert (r > 0).sum() > 0 and (r > 0).sum() <= 384
    assert 1e-8 < r[r > 0].min() <= r.max() < 1e-5          # m/s, order of 0.1-10 mm/h


def test_round_trip_is_exact(tmp_path):
    m = read_irrigation_map(os.path.join(DATA, "salerno_irrigation_map_first_irrigated_block"), NROW, NCOL)
    out = tmp_path / "irrigation_map"
    write_irrigation_map(out, m["times"], m["rates"], irr_sw=m["irr_sw"], irr_flag=m["irr_flag"])
    back = read_irrigation_map(out, NROW, NCOL)
    np.testing.assert_allclose(back["rates"], m["rates"], rtol=1e-7)
    assert back["times"].tolist() == m["times"].tolist() and back["irr_sw"] == 0.35


def test_homogeneous_rates_and_large_times_keep_their_digits(tmp_path):
    out = tmp_path / "irr"
    write_irrigation_map(out, [0.0, 52786700.0], [0.0, 2.5e-8], irr_sw=0.65)
    m = read_irrigation_map(out, 3, 4)
    assert m["irr_hs"] == 1 and m["times"][1] == 52786700.0 and m["rates"][1] == pytest.approx(2.5e-8)


def test_times_must_increase(tmp_path):
    with pytest.raises(ValueError):
        write_irrigation_map(tmp_path / "x", [0.0, 0.0], [0.0, 1.0], irr_sw=0.5)


# ------------------------------------------------------------ schedule generator
def test_build_irrigation_map_events_and_compression():
    from pyCATHY.irrigation import build_irrigation_map

    mask_a = np.zeros((4, 4), bool); mask_a[0:2, 0:2] = True
    mask_b = np.zeros((4, 4), bool); mask_b[2:4, 2:4] = True
    fields = [
        dict(mask=mask_a, start_doy=2, end_doy=3, interval_days=1, start_hour=6, duration_h=2, rate_mm_h=3.6),
        dict(mask=mask_b, start_doy=2, end_doy=2, interval_days=1, start_hour=7, duration_h=2, rate_mm_h=7.2),
    ]
    t, r = build_irrigation_map(fields, (4, 4), t_end=5 * 86400)
    day = 86400.0
    # A on 06-08 h of days 2 and 3; B on 07-09 h of day 2 (overlap 07-08 h)
    assert t.tolist() == [0.0, day + 6 * 3600, day + 7 * 3600, day + 8 * 3600, day + 9 * 3600,
                          2 * day + 6 * 3600, 2 * day + 8 * 3600]
    assert r[0].sum() == 0
    np.testing.assert_allclose(r[1][0, 0], 1e-6)                       # 3.6 mm/h = 1e-6 m/s
    np.testing.assert_allclose(r[2][0, 0], 1e-6); np.testing.assert_allclose(r[2][3, 3], 2e-6)
    np.testing.assert_allclose(r[3][3, 3], 2e-6); assert r[3][0, 0] == 0   # A stopped, B still on
    assert r[4].sum() == 0 and r[-1].sum() == 0


def test_build_irrigation_map_repeats_each_year_and_stops_at_t_end():
    from pyCATHY.irrigation import build_irrigation_map

    mask = np.ones((2, 2), bool)
    fld = dict(mask=mask, start_doy=10, end_doy=10, interval_days=1, start_hour=0, duration_h=1, rate_mm_h=1.0)
    t, r = build_irrigation_map([fld], (2, 2), t_end=2 * 365 * 86400)
    on = t[(r.reshape(len(r), -1).sum(1) > 0)]
    assert on.tolist() == [9 * 86400.0, (365 + 9) * 86400.0]
    t2, _ = build_irrigation_map([fld], (2, 2), t_end=365 * 86400)
    assert (t2 < 365 * 86400).all()


def test_generated_map_is_readable_by_the_writer_and_reader(tmp_path):
    from pyCATHY.irrigation import build_irrigation_map

    mask = np.zeros((3, 3), bool); mask[1, 1] = True
    t, r = build_irrigation_map(
        [dict(mask=mask, start_doy=1, end_doy=5, interval_days=2, start_hour=6, duration_h=3, rate_mm_h=2.0)],
        (3, 3), 10 * 86400)
    write_irrigation_map(tmp_path / "irr", t, r, irr_sw=0.8)
    back = read_irrigation_map(tmp_path / "irr", 3, 3)
    np.testing.assert_allclose(back["rates"], r, rtol=1e-7)
    assert back["times"].tolist() == t.tolist()
