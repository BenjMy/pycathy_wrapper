import numpy as np
import pytest
from pyCATHY.ERT.ert_data import filter_tldat, read_electrodes

TL = """3
# x y z
0 0 0
1 0 0
2 0 0.5
4
# a b m n r err
1 2 3 3 10.0 1.0
1 2 3 3 0.0001 2.0
1 2 3 3 5.0 999.99
1 2 3 3 7.0 3.0
0
"""


@pytest.fixture
def f(tmp_path):
    p = tmp_path / 's.tldat'
    p.write_text(TL.replace('a b m n', 'a b m n').replace('# a b m n r err', '# a b m n r err'))
    return p


def test_electrodes(f):
    e = read_electrodes(f)
    assert e.shape == (3, 3) and e[2].tolist() == [2, 0, 0.5]


def test_filter_err_and_min_r(f, tmp_path):
    n, tot = filter_tldat(f, tmp_path / 'o.tldat', min_abs_r=1e-3)
    assert (n, tot) == (2, 4)
    assert read_electrodes(tmp_path / 'o.tldat').shape == (3, 3)


def test_drop_mask_and_errors(f, tmp_path):
    n, _ = filter_tldat(f, tmp_path / 'o.tldat', drop=[True, False, False])
    assert n == 2
    with pytest.raises(ValueError):
        filter_tldat(f, tmp_path / 'o.tldat', drop=[True])
    with pytest.raises(ValueError):
        filter_tldat(f, tmp_path / 'o.tldat', err_max=0.0)
