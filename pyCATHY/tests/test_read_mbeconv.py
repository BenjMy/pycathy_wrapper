"""B8: a mbeconv file with one single time step used to fail with an obscure shape error."""
import numpy as np
import pytest

from pyCATHY.importers.cathy_outputs import read_mbeconv

HEADER = "MASS BALANCE AND CONVERGENCE BEHAVIOUR\nNSTEP DELTAT TIME ...\n"


def _write(tmp_path, nrows):
    rows = [" ".join(f"{(i + 1) * (j + 1) / 10:.4f}" for j in range(19)) for i in range(nrows)]
    path = tmp_path / "mbeconv"
    path.write_text(HEADER + "\n".join(rows) + ("\n" if rows else ""))
    return str(path)


@pytest.mark.parametrize("nrows", [1, 2, 5])
def test_any_number_of_steps_is_read(tmp_path, nrows):
    df = read_mbeconv(_write(tmp_path, nrows))
    assert len(df) == nrows
    assert list(df.columns)[:3] == ["NSTEP", "DELTAT", "TIME"]
    assert df["NSTEP"].iloc[0] == pytest.approx(0.1)


def test_file_without_a_step_says_so(tmp_path):
    with pytest.raises(ValueError, match="no time step"):
        read_mbeconv(_write(tmp_path, 0))


# Field widths of CATHY's mbeconv FORMAT 1240: I6,2E13.6,I6,E10.3,2E13.6,2E13.5,2E10.3,7E11.3,E10.3
_WIDTHS = [6, 13, 13, 6, 10, 13, 13, 13, 13, 10, 10] + [11] * 7 + [10]
_DECS = [None, 6, 6, None, 3, 6, 6, 5, 5, 3, 3] + [3] * 7 + [3]


def _fortran_row(values):
    out = ""
    for v, w, d in zip(values, _WIDTHS, _DECS):
        out += (f"{int(v)}" if d is None else f"{v:.{d}E}").rjust(w)
    return out


def _values(i):
    v = [float(i + 1), 1.0 + i, 10.0 * (i + 1), 3.0, 1.5] + [0.1 * (j + 1) for j in range(14)]
    return v


def _write_fortran(tmp_path, rows, tail=""):
    path = tmp_path / "mbeconv"
    path.write_text(HEADER + "\n".join(rows) + "\n" + tail)
    return str(path)


def test_negative_E10_3_field_fused_with_previous_column_is_recovered(tmp_path):
    # VIN (col 10) is E10.3: a negative value fills the field and fuses with CUM.DSTORE
    vals = _values(1)
    vals[9] = -2.345e-01
    line = _fortran_row(vals)
    assert len(line.split()) == 18  # the failure mode: 18 tokens instead of 19
    rows = [_fortran_row(_values(0)), line, _fortran_row(_values(2))]
    df = read_mbeconv(_write_fortran(tmp_path, rows))
    assert len(df) == 3
    assert df["VIN"].iloc[1] == pytest.approx(-2.345e-01)
    assert df["CUM.DSTORE"].iloc[1] == pytest.approx(vals[8], rel=1e-4)


def test_e_less_three_digit_exponent_is_recovered(tmp_path):
    vals = _values(0)
    vals[5] = 2.121996e-314  # written by gfortran as "2.121996-314"
    line = _fortran_row(vals).replace("2.121996E-314", "2.121996-314")
    assert "2.121996-314" in line
    df = read_mbeconv(_write_fortran(tmp_path, [line, _fortran_row(_values(1))]))
    assert len(df) == 2
    assert df["STORE1"].iloc[0] == pytest.approx(0.0, abs=1e-300)


def test_truncated_last_line_is_skipped_and_located(tmp_path, caplog):
    rows = [_fortran_row(_values(i)) for i in range(3)]
    truncated = _fortran_row(_values(3))[:60]
    with caplog.at_level("WARNING"):
        df = read_mbeconv(_write_fortran(tmp_path, rows, tail=truncated))
    assert len(df) == 3
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "1 malformed row(s) skipped" in msgs and "last line" in msgs


def test_real_three_line_header_gives_no_warning(tmp_path, caplog):
    """CATHY writes 3 '#' header lines to mbeconv (volume, names, units); readers skip 2: line 3 is not 'malformed'."""
    header = ("#INITIAL VOLUME OF WATER IN THE SUBSURFACE =   1.72937E+02\n"
              "#NSTEP       DELTAT         TIME  NLIN   AVG.LIN ...\n"
              "#                                 ITER.    ITER.      (>0 x inc)\n")
    path = tmp_path / "mbeconv"
    path.write_text(header + "\n".join(_fortran_row(_values(i)) for i in range(3)) + "\n")
    with caplog.at_level("WARNING"):
        df = read_mbeconv(str(path))
    assert len(df) == 3
    assert not [r for r in caplog.records if "malformed" in r.getMessage()]
