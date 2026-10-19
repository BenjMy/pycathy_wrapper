import sys, re, tomllib, types, numpy as np
exec(open('test_archie_scenarios.py').read().split("scen = tomllib")[0])
from pyCATHY.DA.cathy_DA import DA
from pyCATHY.DA.param_tables import _ARCHIE_PARM_NAMES
from pyCATHY.ERT import petro_Archie as Archie

NENS = 64
sc = tomllib.load(open("scenario.toml","rb"))["scenario"]["ic_ZROOT_ArchiePert_upd"]
np.random.seed(42)
entries = perturbate.perturbate(Sim(), sc, NENS)
out = {}
for e in entries:
    out = perturbate.perturbate_parm(out, parm=e, type_parm=e["type_parm"], mean=e["mean"], sd=e["sd"],
            sampling_type=e["sampling_type"], ensemble_size=e["ensemble_size"], per_type=e["per_type"], savefig=False)

# same state as after set_Archie_parm(...) in the script
ns = types.SimpleNamespace(NENS=NENS, Archie_parms_mat_ens=[],
        Archie_parms=dict(porosity=[0.4], rFluid_Archie=[17.0068], a_Archie=[0.3],
                          m_Archie=[1.7], n_Archie=[1.7], pert_sigma_Archie=[0.0]))
for key in [k for k in out if "Archie" in k]:
    kr = re.split(r"(\d+)", key)
    for ens in range(NENS):
        DA.update_archie_parameter(ns, _ARCHIE_PARM_NAMES, kr, "ini_perturbation", ens, out, key)

rho = {}
for sw in (0.3, 0.6, 0.9):
    r = []
    for ens in range(NENS):
        P = Archie.get_Archie_ens_i(ns.Archie_parms, ens)
        r.append(Archie.Archie_rho_DA(rFluid_Archie=P["rFluid_Archie"], sat=[np.array([sw])], porosity=[0.4],
                 a_Archie=P["a_Archie"], m_Archie=P["m_Archie"], n_Archie=P["n_Archie"],
                 pert_sigma_Archie=[0.0])[0])
    rho[sw] = np.array(r)
print("member-wise Archie params -> distinct values:",
      {k: len(set(np.round(ns.Archie_parms[k],6))) for k in ("rFluid_Archie","m_Archie","n_Archie","a_Archie")})
print("a_Archie (unperturbed) first 3 members:", ns.Archie_parms["a_Archie"][:3])
print("\nrho spread over the 64 members (porosity 0.4, no noise):")
for sw, r in rho.items():
    print(f"  Sw={sw}: mean={r.mean():7.2f} Ohm.m  CV={r.std(ddof=1)/r.mean()*100:5.1f} %  range=[{r.min():.1f}, {r.max():.1f}]")
print("\n(for comparison: pert_sigma_Archie noise = 2 % / 5 % CV, per node & per run; ERT_noise in R = 5 %)")
