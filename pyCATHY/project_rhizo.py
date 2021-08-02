# cathy project rhizo lab
from cathy_tools import CATHY
import plottools as cplt
#import cathy_plot as ctplot

import os
import numpy as np

path2prj ='/home/ben/Documents/CATHY/CathyGitbucket/Test_Ben/'
os.chdir(path2prj)


# %% Create the mesh

# we first need to create a 2d surface mesh with the correct attribute to inject into CATHY;
# Usually created using gmsh

#os.chdir('C:/Users/Benjamin/Documents/Simulation/simu-cathy/rhizoLab/create_mesh')

# import the gmsh mesh and create the mesh for CATHY (=grid) (add header_instr)
# header_instr give instruction about how to extrude the 2d mesh
#mesh_dict = ct.create_3dmesh_CATHY(gmsh_mesh='rhizo2d.msh',header_instr=None)

# show mesh
# ctplot.

# %% Prepare input files
os.chdir('/home/ben/Documents/CATHY/rhizoLab/create_mesh/')

#print(path2prj)
survey = CATHY(dirName=path2prj)

survey.create_inputs()

#mesh_dict = test.mshParse('rhizo2d.msh')
zratio=np.ones([1,10])*0.025
survey.create_3dmesh_CATHY('rhizo2d.msh',NZONE=1,NSTR=25,N1=30,NNOD=357,
                            NTRI=700,ZRATIO=zratio,Z1=0)


# test.processor_name= 'cathy'

# %% Visualise Data

# saturation evolution over time
# relationship between relative hydraulic conductivity and pressure head
cplt.showvtk(path2prj+'my_cathy_prj/vtk/100.vtk')
