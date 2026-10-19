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
    with pytest.warns(UserWarning):  # numpy: empty input
        with pytest.raises(ValueError, match="no time step"):
            read_mbeconv(_write(tmp_path, 0))
