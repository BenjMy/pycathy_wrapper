"""Class managing ERT data simulation and inversion + petrophysical transformation
"""


import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv  # if not installed : pip install with conda install pyvista (in a conda terminal)

try:
    import pygimli as pg
    import pygimli.meshtools as mt
    from pygimli.physics import ert
except ImportError:
    pygimli = None

try:
    from resipy import R2  # geophysics tools
except ImportError:
    resipy = None


def create_ERT_survey_pg(pathERT, sequence, mesh, **kwargs):

    verbose = False
    if "verbose" in kwargs:
        verbose = kwargs["verbose"]

    noise_level = 5
    if "noise_level" in kwargs:
        noise_level = kwargs["noise_level"]

    pre, ext = os.path.splitext(mesh)

    try:
        mesh3d = mt.readGmsh(pre + ".msh", verbose=verbose)
    except:
        try:
            mesh3d = pg.load(pre + ".bms", verbose=verbose)
        except:
            raise ValueError(
                "Cannot read " + mesh + ": valid extensions are .msh and .bms"
            )

    if "shm" in sequence:
        scheme = pg.load(sequence)
    elif "txt" in sequence:
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

    res0 = 1
    if "res0" in kwargs:
        res0 = kwargs["res0"]

    if len(res0) != len(mesh3d.cells()):
        raise ValueError("wrong initial resistivity input")

    try:
        if not verbose:
            devnull = open("/dev/null", "w")
            oldstdout_fno = os.dup(sys.stdout.fileno())
            os.dup2(devnull.fileno(), 1)
    except:
        pass
    het = ert.simulate(
        mesh3d,
        res=res0,
        scheme=scheme,
        calcOnly=False,
        verbose=verbose,
        noiseLevel=noise_level,
    )

    try:
        if not verbose:
            os.dup2(oldstdout_fno, 1)
    except:
        pass

    # # Reset the standard output
    # sys.stdout = sys.__stdout__

    # # Discard the captured output
    # captured_output.close()

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

    ERT = R2(pathERT, typ="R3t")  # + 'ERT_fwdmodel'
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


def fwd_ERT_survey(ERT, noise, show=False, dump=[]):
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
