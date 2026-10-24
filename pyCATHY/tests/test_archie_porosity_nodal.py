"""SW_2_ER0_DA must apply the porosity of EACH NODE (layers / zones) in the Archie law.

Input shapes reaching SW_2_ER0_DA:
  A  ERT ensemble (cathy_DA._ert_ensemble_inputs): object array, one entry per member = nodal list;
     entries of invalid members are floats (np.ones(NENS, dtype=object) * [mean])
  B  EM mapper (mapper._map_EM): the nodal list of ONE member
  C  legacy: one porosity per zone, e.g. [0.4]
Run:  python test_archie_porosity_nodal.py
"""
import numpy as np
import pandas as pd
from pyCATHY.ERT import petro_Archie as Archie

NN = 12                                    # nodes
layer_phi = np.repeat([0.415, 0.41, 0.35, 0.41], NN // 4)          # layered / compacted-zone porosity
sw = pd.Series(np.linspace(0.3, 0.9, NN))
P = dict(porosity=[0.40], rFluid_Archie=[17.0], a_Archie=[0.3], m_Archie=[1.7], n_Archie=[1.7],
         pert_sigma_Archie=[0.0])


def expected(phi, P=P):
    return P["rFluid_Archie"][0] * P["a_Archie"][0] * np.asarray(phi) ** (-P["m_Archie"][0]) * sw.values ** (-P["n_Archie"][0])


def run(poro, ens=0, P=P):
    rho, df_Archie, _ = Archie.SW_2_ER0_DA("prj", P, poro, {}, "x/cathy_%d" % (ens + 1), df_sw=sw, Ens_nbi=ens, DA_cnb=1)
    return rho, df_Archie


# A: ensemble, member 1 asked; member 0 invalid (float, as built by _ert_ensemble_inputs); member 2 different porosity
ens = np.ones(3, dtype=object) * [0.40]
ens[1], ens[2] = list(layer_phi), list(layer_phi * 0.9)
rho, df = run(ens, ens=1)
assert np.allclose(rho, expected(layer_phi)), "A: member 1 must use ITS nodal porosity"
assert np.allclose(df["porosity"], layer_phi), "A: df_Archie must report the nodal porosity used"
rho2, _ = run(ens, ens=2)
assert np.allclose(rho2, expected(layer_phi * 0.9)), "A: member 2"
# B: one member, nodal list / array
for poro in (list(layer_phi), np.asarray(layer_phi)):
    assert np.allclose(run(poro)[0], expected(layer_phi)), "B: nodal porosity of a single member"
# C: one porosity per zone (legacy) -> broadcast
assert np.allclose(run([0.40])[0], expected(np.full(NN, 0.40))), "C: scalar porosity"
# wrong size must fail loudly
try:
    run(list(layer_phi[:5])); raise SystemExit("size mismatch accepted")
except ValueError:
    pass
print("porosity is nodal in every input shape")
