"""Class managing ERT data simulation and inversion + petrophysical transformation
"""


import contextlib
import io
import logging
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv  # if not installed : pip install with conda install pyvista (in a conda terminal)

from pyCATHY.ERT import ERT_checks

try:
    import pygimli as pg
    import pygimli.meshtools as mt
    from pygimli.physics import ert
except ImportError:
    pygimli = None

try:
    from resipy import Project  # geophysics tools
    #from resipy import R2  # geophysics tools
except ImportError:
    resipy = None

logger = logging.getLogger(__name__)


@contextlib.contextmanager
def _captured_stdout(verbose):
    """Keep pyGIMLi/BERT stdout out of the console: it is sent to the logger (DEBUG).

    Python-level prints are captured, and (when stdout is a real file descriptor,
    i.e. not in a notebook) the C-level output is muted. stdout is always restored,
    even if the block raises. ``verbose=True`` lets everything through.
    """
    if verbose:
        yield
        return
    buf = io.StringIO()
    saved_fd = devnull_fd = None
    try:
        sys.stdout.flush()
        saved_fd = os.dup(sys.stdout.fileno())
        devnull_fd = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull_fd, sys.stdout.fileno())
        fd = sys.stdout.fileno()
    except (AttributeError, OSError, ValueError):  # no real fd (notebook, pytest capture)
        fd = None
    try:
        with contextlib.redirect_stdout(buf):
            yield
    finally:
        if fd is not None:
            os.dup2(saved_fd, fd)
        for handle in (saved_fd, devnull_fd):
            if handle is not None:
                os.close(handle)
        text = buf.getvalue().strip()
        if text:
            logger.debug("pyGIMLi output:\n%s", text)

# report of the last `create_ERT_survey_pg` call of this process (read by petro_Archie)
_STATE = {"last_report": None}


def get_last_check_report():
    """ERTCheckReport of the last forward run done in this process (None if none)."""
    return _STATE["last_report"]


def create_ERT_survey_pg(pathERT, 
                         sequence, 
                         mesh,
                         **kwargs
                         ):

    verbose = False
    if "verbose" in kwargs:
        verbose = kwargs["verbose"]

    fwdNoiseLevel = 5
    if "fwdNoiseLevel" in kwargs:
        fwdNoiseLevel = kwargs["fwdNoiseLevel"]
    logger.debug("fwd ERT noise level: %s%%", fwdNoiseLevel)
    pre, ext = os.path.splitext(mesh)
    logger.debug("forward mesh: %s%s", pre, ext)

    mesh3d = None
    for reader, ext_try in (
        (lambda f: mt.readGmsh(f, verbose=verbose), ".msh"),
        (lambda f: pg.load(f, verbose=verbose), ".bms"),
        (lambda f: pg.load(f, verbose=verbose), ".vtk"),
    ):
        try:
            mesh3d = reader(pre + ext_try)
            logger.debug("forward mesh read from %s%s", pre, ext_try)
            break
        except Exception as err:  # missing file or unreadable format: try the next one
            logger.debug("cannot read %s%s: %s", pre, ext_try, err)
    if mesh3d is None:
        raise ValueError(
            "Cannot read {}: valid extensions are .msh, .bms, or .vtk".format(mesh)
        )
    
    sequence_file_extension = os.path.splitext(sequence)[1]

    if "shm" in sequence_file_extension:
        scheme = pg.load(sequence)
    elif "dat" in sequence_file_extension:
        scheme = pg.load(sequence)
    elif "txt" in sequence_file_extension:
        shm = pd.read_csv(sequence, delimiter=" ", header=None)
        shm = shm.to_numpy()
        shm_new = np.subtract(shm, 1)
        scheme = pg.DataContainerERT()
        if "dict_ERT" in kwargs:
            scheme.setSensorPositions(kwargs["dict_ERT"]["elecs"])
        for i, elec in enumerate("abmn"):
            scheme[elec] = shm_new[:, i]
    else:
        raise ValueError(
            "Sequence file format not recognized use .shm or .txt as a list of quadrupoles"
        )
    
    n_elec = scheme.sensorCount() if hasattr(scheme, "sensorCount") else None
    ERT_checks.check_scheme(
        scheme["a"], scheme["b"], scheme["m"], scheme["n"],
        n_elec=n_elec, name="fwd ERT scheme",
    ).warn()

    if "res0" not in kwargs:
        raise ValueError("create_ERT_survey_pg needs res0 (one resistivity per mesh cell)")
    res0 = np.asarray(kwargs["res0"], dtype=float)

    if res0.ndim != 1 or len(res0) != len(mesh3d.cells()):
        raise ValueError(
            f"wrong initial resistivity input: shape {res0.shape} for "
            f"{len(mesh3d.cells())} mesh cells"
        )
    # Check for negative values
    if np.any(res0 < 0):
        raise ValueError("res0 contains negative values. Resistivity must be positive.")
        
    
    # print('+'*20)
    # # print(sequence_file_extension)
    # # print(sequence)
    # print(scheme)
    # print(mesh3d)
    # print(fwdNoiseLevel)
    # print('+'*20)

    if np.isnan(res0).any():
        raise ValueError("⚠️ Error: res0 contains NaN values! Stopping execution.")
    # Check for negative values
    if np.any(res0 <= 0):
        raise ValueError("res0 contains negative values. Resistivity must be positive.")
        
        

    # optional absolute voltage noise of the simulation (pyGIMLi default if not given)
    sim_kwargs = {}
    if "fwdNoiseAbs" in kwargs:
        sim_kwargs["noiseAbs"] = kwargs["fwdNoiseAbs"]

    with _captured_stdout(verbose):
        het = ert.simulate(
            mesh3d,
            res=res0,
            scheme=scheme,
            calcOnly=False,
            verbose=False,
            noiseLevel=fwdNoiseLevel,
            **sim_kwargs,
        )
    # pg.show(mesh3d)
    # pg.show(het)
    # het.exportVTK('testmeshpg.vtk')
    # het.saveResult('testpg')
    
    logger.debug("fwd rhoa (raw) min/max: %.4g / %.4g", min(het['rhoa']), max(het['rhoa']))
    # checks on the raw forward data (before the rhoa <= 0 correction below)
    # 1e-4 V is the assumed pyGIMLi default of `noiseAbs` when not set explicitly
    report = ERT_checks.check_forward_data(
        het,
        res0=res0,
        noise_abs=kwargs.get("fwdNoiseAbs", 1e-4),
        name="fwd ERT data",
    )
    report.warn()
    _STATE["last_report"] = report
    n_fixed = int(np.sum(np.asarray(het['rhoa']) <= 0))
    if n_fixed:
        logger.warning(
            "fwd ERT: %d of %d apparent resistivities are <= 0 (min %.3g Ohm.m): "
            "set to 1e-3 Ohm.m so that they can be log-transformed; check the "
            "quadrupoles / noise level (see the ERT check report)",
            n_fixed,
            len(het['rhoa']),
            float(np.min(np.asarray(het['rhoa']))),
        )
    # het.loc[het['rhoa'] <= 0, 'rhoa'] = 1e-3
    het['rhoa'][het['rhoa'] <= 0] = 1e-3

    # print('*0+'*20)
    # print(res0)
    # print('*0+'*20)



    return het


def create_ERT_survey_Resipy(pathERT, elecsXYZ, sequence, mesh, **kwargs):
    """
    #https://hkex.gitlab.io/resipy/api.html
    """

    noise_level = 5
    if "noise_level" in kwargs:
        noise_level = kwargs["noise_level"]

    # os.chdir(pathERT)

    isExist = os.path.exists(pathERT)

    if not isExist:

        # Create a new directory because it does not exist
        os.makedirs(pathERT)

    ERT = Project(pathERT, typ="R3t")  # + 'ERT_fwdmodel'
    ERT.setTitle("Rhizo_synth")

    if type(elecsXYZ) is str:
        # elecsXYZ = np.loadtxt(elecsXYZ, delimiter='\t')
        elecsXYZ = np.genfromtxt(elecsXYZ, delimiter=",", skip_header=1)

    ERT.setElec(np.c_[elecsXYZ[:, 0], elecsXYZ[:, 2], elecsXYZ[:, 1], elecsXYZ[:, 3]])

    ERT.importMesh(mesh)

    if "res0" in kwargs:
        ERT.setRefModel(kwargs["res0"])

    ERT.importSequence(sequence)
    # -----------------------------------------------

    ERT.forward(noise=noise_level)

    return ERT


def fwd_ERT_survey(ERT, noise, show=False, dump=None):
    """
    Fwd ERT model

    Parameters
    ----------
    ERT : TYPE
        DESCRIPTION.
    noise : TYPE
        DESCRIPTION.
    show : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    ERT : TYPE
        DESCRIPTION.

    """

    # ----------#
    # fwd modelling
    # ---------------------------------------------------------#

    if dump is None:
        dump = []
    ERT.forward(
        noise=0.05, iplot=show, dump=dump
    )  # forward modelling with 5 % noise added to the output

    return ERT


def invert_ERT_survey(ERT, show=False):

    # ---------------------------------------------------------#
    # inversion
    # ---------------------------------------------------------#
    ERT.param["num_xy_poly"] = 0
    ERT.param["zmin"] = -np.inf
    ERT.param["zmax"] = np.inf
    ERT.param["data_type"] = 1  # using log of resistitivy
    ERT.err = False  # if we want to use the error from the error models fitted before
    ERT.param["a_wgt"] = 0.001
    ERT.param["b_wgt"] = 0.05

    ERT.invert()

    if show == True:
        ERT.showResults(index=0)  # show the initial model
        ERT.showResults(index=1, sens=False)  # show the inverted model


def Archie_ERT(ERT, rFluid, porosity, m, n):
    """
    Archie transformation
    (only possible on inverted values)

    Parameters
    ----------
    ERT : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    # (rho/rFluid/porosity^-m)^(-1/n)
    rFluid = 1
    porosity = 0.55
    m = 1
    n = 1

    str_formula = f"(x['Resistivity']* {rFluid} * {porosity} **(-m))**(-1/ {n})"
    # print(str_formula)

    # only possible on inverted mesh (meshResults)
    ERT.computeAttribute(str_formula, name="saturation")

    pass
