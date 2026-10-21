import pytest
from pyCATHY.DA.cathy_DA import DA


class Ctx:
    dict_parm_pert, update_key, cycle_nb = {}, "update_nb1", 1


@pytest.mark.parametrize("name", ["s", "icx", "Ks_x", "k", "i", "or", "atm", "WT"])
def test_unknown_or_partial_names_raise(name):
    d = DA.__new__(DA)
    with pytest.raises(ValueError):
        d._update_parameter_member(Ctx, name, [name, "", ""], 0)


@pytest.mark.parametrize("name,method", [("atmbc", "_update_atmbc_member"), ("WTPOSITION", "_update_wtposition_member"),
                                         ("ic", "_update_ic_member")])
def test_known_names_dispatch(name, method):
    d = DA.__new__(DA)
    called = []
    setattr(d, method, lambda *a, **k: called.append(name))
    d._update_parameter_member(Ctx, name, [name, "0", ""], 0)
    assert called == [name]


def test_ks_and_porosity_case_insensitive():
    d = DA.__new__(DA)
    hit = []
    d._update_soil_params_from_perturbations = lambda *a, **k: hit.append(1)
    d.update_archie_parameter = lambda *a, **k: None
    ctx = type("C", (), {"dict_parm_pert": {}, "update_key": "u", "cycle_nb": 1, "ks": None, "porosity": None})()
    for n in ("Ks", "ks", "Porosity"):
        d._update_parameter_member(ctx, n, [n, "0", ""], 0)
    assert len(hit) == 3
