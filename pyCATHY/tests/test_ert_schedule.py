import json
import pandas as pd
import pytest
import pyCATHY.ERT.schedule as S


def test_times_relative_to_irrigation(tmp_path):
    sp = tmp_path / 'sp.json'
    sp.write_text(json.dumps({'INFILTRATION_START': {'NC': '2025-05-12 10:00', 'HC': '2025-06-01 09:00'}}))
    t0 = S.load_infiltration_start(sp, 'no_comp')
    assert t0 == pd.Timestamp('2025-05-12 10:00')
    c = tmp_path / 'a.csv'
    c.write_text("file,datetime\nst5t00.tldat,12/05/2025 09:00\nst5t01.tldat,12/05/2025 13:30\n")
    acq = S.read_acquisition_times(c)
    assert S.times_h_since_irrigation(acq, t0) == {0: -1.0, 1: 3.5}


def test_survey_column_and_duplicates(tmp_path):
    c = tmp_path / 'b.csv'
    c.write_text("survey,datetime\n0,2025-05-12 09:00\n0,2025-05-12 10:00\n")
    with pytest.raises(ValueError):
        S.read_acquisition_times(c, dayfirst=False)


def test_unknown_label_is_an_error(tmp_path):
    with pytest.raises(KeyError):
        S.load_infiltration_start('x.json', 'light_comp')


def test_check_schedule():
    assert S.check_schedule({0: -1, 1: 3, 2: 9}, irr_start_hour=2, sim_end_h=72) == []
    assert S.check_schedule({0: -3, 1: 3}, 2)                 # t_cat = -1 h
    assert S.check_schedule({0: 3, 1: 2}, 2)                  # not increasing
    assert S.check_schedule({0: 1, 1: 80}, 2, sim_end_h=72)   # after the end
