#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update and show Boundary conditions
============================================

This example shows how to use pyCATHY object to update a 3d BC properties from a DEM and run the hydrological model.

Questions:
    - what can of boundary condition are set at the outlet point by default?
    - Default side boundary conditions?

For a good reading see: Surface-subsurface flow modeling with path-based runoff 
routing, boundary condition-based coupling, and assimilation of multisource observation data
WATER RESOURCES RESEARCH, VOL. 46, W02512, doi:10.1029/2008WR007536, 2010

Exemples treated in the notebook:
    
    - 1st: 
    - 2nd example: the entire outlet face is a no-flow boundary and water is allowed to leave the system only by exfiltration



*Estimated time to run the notebook = 5min*

"""


# !! run preprocessor change the DEM shape !
# dtm_13 does not have the same shape anymore!

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pyCATHY.meshtools as mt
from pyCATHY import cathy_tools
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.plotters import cathy_plots as cplt

#%% Init CATHY model
# ------------------------
path2prj = "../SSHydro/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, prj_name="soil3d_from_weil", clear_src=True)

rootpath = os.path.join(simu.workdir + simu.project_name)
#%% Fetch and show initial DEM
# the dimension of the mesh is squared (20,20)


dem_mat, str_hd_dem = in_CT.read_dem(
    os.path.join(simu.workdir, simu.project_name, "prepro/dem"),
    os.path.join(simu.workdir, simu.project_name, "prepro/dtm_13.val"),
)

simu.show_input(prop="dem")

maxdepth = 10
# # linear z depth
# # -------------------------------------------------------------
zb = np.linspace(0, maxdepth, 10)
nstr = len(zb)
zr = list((np.ones(len(zb))) / (nstr))

sum(zr)

# zb = np.geomspace(1e-1, maxdepth, num=15)
# nstr=len(zb)
# zr = [abs(zb[0]/maxdepth)]
# zr.extend(list(abs(np.diff(zb)/maxdepth)))

# np.shape(dem_mat)
simu.update_prepo_inputs(
    DEM=dem_mat,
    nstr=nstr,
    zratio=zr,
    base=max(zb),
)
fig = plt.figure()
ax = plt.axes(projection="3d")
simu.show_input(prop="dem", ax=ax)

# simu.update_soil()
simu.create_mesh_vtk(verbose=True)



#%%

# .. note:
#     The boundary conditions are defined in the nansfdirbc (Dirichlet),
#     nansfneubc (Neumann), and sfbc (seepage face) files.

#     We have two types of boundary conditions (BC):
#     - Neumann BC (or specifed flux)
#     - Dirichlet BC (or pressure).


# .. note:
#     - Pioggia: condizioni di Neumann. Quando non ci può più essere
#     infiltrazione metto Dirichlet.
#     - Evaporazione: si indica un limite di pressione minimo ( Pmin ) al di
#     sotto del quale si ha uno switch da Neumann a Dirichlet
#     (in quanto al di sotto di questo valore non si ha più evapotraspirazione).

# .. note:
#     The boundary condition for any given surface node can switch between a
#     Dirichlet condition and a Neumann condition depending on the saturation
#     (or pressure) state of that node.

# .. note:
#     A Neumann (or specified flux) boundary condition corresponds to
#     atmosphere-controlled infiltration or exfiltration, with the flux equal
#     to the rainfall or potential evaporation rate given by the atmospheric input data.
#     When the surface node reaches a threshold level of saturation or moisture deficit,
#     the boundary condition is switched to a Dirichlet (specified head) condition,
#     and the infiltration or exfiltration process becomes soil limited [1].



#%%

# exemple provided by Laura B.
# ----------------------------

# C     Write dirbc
#       write(33,*) 0.0, 'time'
#       write(33,*) '0', a
#       do i=1,nnod3
#          if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
#      1       (y(i).eq.5))then
#          write(33,*) i
#          endif
#       enddo
#       do i=1,nnod3
#          if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
#      1       (y(i).eq.5))then
#          write(33,*) -z(i)-WTdepth
#          endif
#       enddo

#       write(33,*) 2e+20, 'time'
#       write(33,*) '0', a
#       do i=1,nnod3
#          if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
#      1       (y(i).eq.5))then
#          write(33,*) i
#          endif
#       enddo
#       do i=1,nnod3
#          if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
#      1       (y(i).eq.5))then
#          write(33,*) -z(i)-WTdepth
#          endif
#       enddo

# modicare il valore di NPMAX nel file 27 CATHY.H nel caso
# in cui si inseriscano dei NDIRC ed il valore di NP2MAX nel caso si inseriscano dei
# NDIR. I valori di NPMAX e NP2MAX corrispondono al numero massimo
# di nodi NDIRC e NDIR che si possono inserire.

#%%
# Dirichlet Boundary conditions (or specified pressure) at time t

# - To simulate the no-flow boundaries conditions for the bottom and
#   vertical sides of the domain it is necessary to set NDIR and NDIRC
#   equal to zero.
# - To simulate different boundary conditions, it is necessary to
#   indicate the number of selected nodes through NDIR or NDIRC,
#   then to specify the node ID’s that you want to consider and
#   eventually the value of pressure head or flux that you want to assign.
%matplotlib auto

try:
    del simu.mesh_bound_cond_df
except:
    pass
simu.update_nansfdirbc(no_flow=True)

#%%
# Neumann boundary conditions (or specifed flux) at time t
try:
    del simu.mesh_bound_cond_df
except:
    pass
simu.update_nansfneubc(no_flow=True)

#%%
simu.update_sfbc()
#%%
simu.show_bc()

#%% Run  hydrological model

simu.run_processor(IPRT1=2, verbose=True)
