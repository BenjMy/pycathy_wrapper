#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""For the manager look at BERT https://gitlab.com/resistivity-net/bert."""

import numpy as np
import matplotlib.pyplot as plt
import os
import pyvista as pv
from pyCATHY import meshtools as mt 
from pyCATHY.rhizo_tools import CATHY_2_Resipy
from pyCATHY.ERT import simulate_ERT as simuERT

def SW_2_ERa(df_sw,workdir,project_name,
             porosity,
             path_CATHY,
             pathERT, meshERT, elecs, sequenceERT):
    
    
   
    ER_converted_ti = Archie_rho(rFluid=1, 
                                sat = df_sw[-1],
                                porosity=porosity, 
                                a=1.0, m=2.0, n=2.0)

    # read in CATHY mesh data
    # ------------------------------------------------------------------------
    # path_CATHY = '/home/ben/Documents/CATHY/pyCATHY/mary_rhizo_withDA/rhizo_ET_irr_PRD_withVp/vtk/'
    # path_CATHY = os.path.join(workdir, project_name , 'vtk/')
    mesh_CATHY = pv.read(path_CATHY + '100.vtk')
    
    # print(ER_converted_ti)

    # add attribute converted to CATHY mesh
    # ------------------------------------------------------------------------
    mesh_CATHY_new_attr = mt.add_attribute_2mesh(ER_converted_ti,
                            mesh_CATHY,
                            'ER_converted',
                            overwrite=True, 
                            path = path_CATHY)

    # print(mesh_CATHY)
    # print(mesh_CATHY_new_attr)

    # copy attribute to resipy mesh
    # ------------------------------------------------------------------------
    mesh_Resipy_new_attr, scalar_new = CATHY_2_Resipy(mesh_CATHY_new_attr,meshERT,scalar='ER_converted',
                   show=False, path= os.path.join(workdir, project_name, 'vtk/'))

    # fwd ERT data
    # ------------------------------------------------------------------------
    res0 = mesh_Resipy_new_attr.get_array(scalar_new)
    # mesh_Resipy_new_attr.cell_data[scalar_new] 
    # mesh_Resipy_new_attr.array_names

    # ERT.mesh
    # res = ERT.mesh.df['res0']
    ERT_predicted = simuERT.create_ERT_survey(os.path.join(pathERT, project_name,'predicted'), 
                                                elecs, 
                                                sequenceERT, 
                                                meshERT, 
                                                res0=res0)
    ERT_predicted = simuERT.fwd_ERT_survey(ERT_predicted, noise=10)
    df_ERT_predicted = ERT_predicted.surveys[0].df

    # save to dataframe and export file
    # ------------------------------------------------------------------------
    filename = 'ER_predicted.csv'
    
    isExist = os.path.exists(os.path.join(pathERT, project_name))
    if not isExist:
        os.makedirs(os.path.join(pathERT, project_name))

    df_ERT_predicted.to_csv(os.path.join(pathERT, project_name, filename))
    # ERT_predicted[]
    # simuERT.invert_ERT_survey(ERT)
    
    
    return df_ERT_predicted
  

def Archie_rho(rFluid, sat, porosity, a=1.0, m=2.0, n=2.0):
    '''
    Compute ER values at each mesh nodes

    Parameters
    ----------
    rFluid : int
        conductivity of the pore fluid.
    sat : np.array
        water saturation.
    porosity : int or np.array
        DESCRIPTION.
    a : float, optional
        tortuosity factor. The default is 1.0.
    m : float, optional
        cementation exponent. The default is 2.0.(usually in the range 1.3 -- 2.5 for sandstones)
    n : float, optional
        saturation exponent. The default is 2.0. 

    Returns
    -------
    TYPE
        Electrical Resistivity (Ohm.m).

    '''

    return rFluid * a * porosity**(-m) * sat**(-n)



def Archie_sat(rho, rFluid, porosity, a=1.0, m=2.0, sat=1.0, n=2.0):
    '''
    
    rho: resistivity
    𝑆𝑤 : water saturation
    𝜙: the porosity of the soil
    𝜎_{𝑤} is the conductivity of the pore fluid
    𝑎, 𝑚, and 𝑛 are empirically derived parameters
    𝑎 is the tortuosity factor
    𝑚 is the cementation exponent
    𝑛 is the saturation exponent
    Returns
    -------
    𝑆𝑤 : water saturation


    '''
    return  (rho/rFluid/porosity^-m)^(-1/n)


# import pygimli as pg
# import pygimli.meshtools as mt

# def Archie_sat(rho, rFluid, porosity, a=1.0, m=2.0, sat=1.0, n=2.0):
#     '''
#     rho: resistivity
#     𝑆𝑤 : water saturation
#     𝜙: the porosity of the soil
#     𝜎_{𝑤} is the conductivity of the pore fluid
#     𝑎, 𝑚, and 𝑛 are empirically derived parameters
#     𝑎 is the tortuosity factor
#     𝑚 is the cementation exponent
#     𝑛 is the saturation exponent
#     Returns
#     -------
#     𝑆𝑤 : water saturation


#     '''
#     return  (rho/rFluid/porosity^-m)^(-1/n)
    

# def Archie_rho(rFluid, porosity, a=1.0, m=2.0, sat=1.0, n=2.0):
#     '''
#     sat 𝑆𝑤 : water saturation
#     porosity phi 𝜙: the porosity of the soil
#     rFluid is the conductivity of the pore fluid
#     𝑎, 𝑚, and 𝑛 are empirically derived parameters
#     𝑎 is the tortuosity factor
#     𝑚 is the cementation exponent
#     𝑛 is the saturation exponent
#     Returns
#     -------
#     Resistivity

#     '''
#     return rFluid * a * porosity**(-m) * sat**(-n)



# [docs]def resistivityArchie(rFluid, porosity, a=1.0, m=2.0, sat=1.0, n=2.0,
#                       mesh=None, meshI=None, fill=None, show=False):
#     r"""
#     Resistivity of rock for the petrophysical model from Archies law.

#     Calculates resistivity of rock for the petrophysical model from Archie's law.
#     :cite:`Archie1942`

#     .. math::
#         \rho = a\rho_{\text{fl}}\phi^{-m} S^{-n}

#     * :math:`\rho` - the electrical resistivity of the fluid saturated rock in
#       :math:`\Omega\text{m}`
#     * :math:`\rho_{\text{fl}}` - rFluid: electrical resistivity of the fluid in
#       :math:`\Omega\text{m}`
#     * :math:`\phi` - porosity 0.0 --1.0
#     * :math:`S` - fluid saturation 0.0 --1.0 [sat]
#     * :math:`a` - Tortuosity factor. (common 1)
#     * :math:`m` - Cementation exponent of the rock (usually in the range 1.3 -- 2.5 for sandstones)
#     * :math:`n` - is the saturation exponent (usually close to 2)

#     If mesh is not None the resulting values are calculated for each cell of the mesh. All
#     parameter can be scalar, array of length mesh.cellCount() or callable(pg.cell). If rFluid is
#     non-steady n-step distribution than rFluid can be a matrix of size(n, mesh.cellCount()) If
#     meshI is not None the result is interpolated to meshI.cellCenters() and prolonged (if fill
#     ==1).

#     Notes ----- We experience some unstable nonlinear behavior. Until this is clarified all results
#     are rounded to the precision 1e-6.

#     Examples -------- >>> #

#     WRITEME
#     """
#     phi = porosity
#     if isinstance(porosity, list):
#         phi = np.array(porosity)

#     if mesh is None:
#         return rFluid * a * phi**(-m) * sat**(-n)

#     rB = None

#     if isinstance(rFluid, float):
#         rB = pg.Matrix(1, mesh.cellCount())
#         rB[0] = pg.solver.parseArgToArray(rFluid, mesh.cellCount(), mesh)

#     elif isinstance(rFluid, pg.Vector):
#         rB = pg.Matrix(1, len(rFluid))
#         rB[0] = pg.solver.parseArgToArray(rFluid, mesh.cellCount(), mesh)

#     elif hasattr(rFluid, 'ndim') and rFluid.ndim == 1:
#         rB = pg.Matrix(1, len(rFluid))
#         rB[0] = pg.solver.parseArgToArray(rFluid, mesh.cellCount(), mesh)

#     elif hasattr(rFluid, 'ndim') and rFluid.ndim == 2:
#         rB = pg.Matrix(len(rFluid), len(rFluid[0]))
#         for i, rFi in enumerate(rFluid):
#             rB[i] = rFi

#     phi = pg.solver.parseArgToArray(phi, mesh.cellCount(), mesh)
#     a = pg.solver.parseArgToArray(a, mesh.cellCount(), mesh)
#     m = pg.solver.parseArgToArray(m, mesh.cellCount(), mesh)
#     S = pg.solver.parseArgToArray(sat, mesh.cellCount(), mesh)
#     n = pg.solver.parseArgToArray(n, mesh.cellCount(), mesh)

#     if show:
#         pg.show(mesh, S, label='S')
#         pg.show(mesh, phi, label='p')
#         pg.wait()

#     r = pg.Matrix(len(rB), len(rB[0]))
#     for i, _ in enumerate(r):
#         r[i] = rB[i] * a * phi**(-m) * S**(-n)

#     r.round(1e-6)

#     if meshI is None:
#         if len(r) == 1:
#             return r[0].copy()
#         return r

#     rI = pg.Matrix(len(r), meshI.cellCount())
#     if meshI:
#         pg.interpolate(mesh, r, meshI.cellCenters(), rI)

#     if fill:
#         for i, ri_ in enumerate(rI):
#             # slope == True produce unstable behavior .. check!!!!!!
#             rI[i] = mt.fillEmptyToCellArray(meshI, ri_, slope=False)

#     rI.round(1e-6)

#     if len(rI) == 1:
#         # copy here because of missing refcounter TODO
#         return rI[0].array()
#     return rI



# [docs]def transFwdArchiePhi(rFluid=20, m=2):
#     r"""
#     Archies law transformation function for resistivity(porosity).

#     .. math::
#         \rho & = a\rho_{\text{fl}}\phi^{-m}\S_w^{-n} \\
#         \rho & = \rho_{\text{fl}}\phi^(-m) =
#         \left(\phi/\rho_{\text{fl}}^{-1/n}\right)^{-n}

#     See also :py:mod:`pygimli.physics.petro.resistivityArchie`

#     Returns ------- trans : :gimliapi:`GIMLI::RTransPower` Transformation function

#     Examples -------- >>> from pygimli.physics.petro import * >>> phi = 0.3 >>> tFAPhi =
#     transFwdArchiePhi(rFluid=20) >>> r1 = tFAPhi.trans(phi) >>> r2 = resistivityArchie(rFluid=20.0,
#     porosity=phi, ... a=1.0, m=2.0, sat=1.0, n=2.0) >>> print(r1-r2 < 1e-12) True >>> phi = [0.3]
#     >>> tFAPhi = transFwdArchiePhi(rFluid=20) >>> r1 = tFAPhi.trans(phi) >>> r2 =
#     resistivityArchie(rFluid=20.0, porosity=phi, ... a=1.0, m=2.0, sat=1.0, n=2.0) >>> print((r1-r2
#     < 1e-12)[0]) True
#     """
#     return pg.trans.TransPower(-m, rFluid**(1./m))



# [docs]def transInvArchiePhi(rFluid=20, m=2):  # phi(rho)
#     """
#     Inverse Archie transformation function porosity(resistivity).

#     # rFluid/rho = phi^m ==> phi = (rFluid/rho)^(1/m) = (rho/rFluid)^(-1/m) See ---
#     :py:mod:`pygimli.physics.petro.transFwdArchiePhi`
#     """
#     return pg.trans.TransPower(-1/m, rFluid)



# [docs]def transFwdArchieS(rFluid=20, phi=0.4, m=2, n=2):  # rho(S)
#     """Inverse Archie transformation function resistivity(saturation)."""
#     # rho = rFluid * phi^(-m) S^(-n)
#     return pg.trans.TransPower(-n, (rFluid*phi**(-m))**(1/n))



# [docs]def transInvArchieS(rFluid=20, phi=0.4, m=2, n=2):  # S(rho)
#     """Inverse Archie transformation function saturation(resistivity)."""
#     # rFluid/rho = phi^m S^n => S=(rFluid/rho/phi^m)^(1/n)
#     # S = (rho/rFluid/phi^-m)^(-1/n)
#     return pg.trans.TransPower(-1/n, rFluid*phi**(-m))



# def test_Archie():
#     """Test Archie."""
#     import unittest

#     dx = 0.01
#     phivec = np.arange(dx, 0.5, dx)
#     swvec = np.arange(dx, 1, dx)

#     phi0 = 0.4  # 40% porosity
#     rhow = 20  # 20 Ohmm tap water
#     tFAPhi = transFwdArchiePhi(rFluid=rhow)
#     tFAS = transFwdArchieS(rFluid=rhow, phi=phi0)
#     tIAPhi = transInvArchiePhi(rFluid=rhow)
#     tIAS = transInvArchieS(rFluid=rhow, phi=phi0)

#     ax = plt.subplots()[1]
#     # direct function
#     rA = resistivityArchie(rFluid=rhow, porosity=phivec)
#     rS = resistivityArchie(rFluid=rhow, porosity=phi0, sat=swvec)
#     ax.semilogy(phivec, rA, 'b-')
#     ax.semilogy(swvec, rS, 'r-')

#     # forward transformation
#     fA = tFAPhi.trans(phivec)
#     fS = tFAS.trans(swvec)
#     np.testing.assert_allclose(rA, fA, rtol=1e-12)
#     np.testing.assert_allclose(rS, fS, rtol=1e-12)

#     ax.semilogy(phivec, fA, 'bx', markersize=10)
#     ax.semilogy(swvec, fS, 'rx', markersize=10)

#     # inverse transformation
#     iA = tIAPhi.invTrans(phivec)
#     iS = tIAS.invTrans(swvec)
#     np.testing.assert_allclose(rA, iA, rtol=1e-12)
#     np.testing.assert_allclose(rS, iS, rtol=1e-12)

#     ax.semilogy(phivec, iA, 'bo')
#     ax.semilogy(swvec, iS, 'ro')

#     plt.show()


# if __name__ == "__main__":
#     pass