#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Some tips to use pybert for closed geometries (e.g. rhizotron).
"""
import os 
#os.chdir('e:/Padova/Z_GitlabLocal/Exemple_Luca/3dinv_2dplot')
#os.chdir('E:/Padova/Experiments/8_Rhizotron_DAPHNAE_Lancaster/Z_Codes
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.mplviewer import drawSensors

import pybert as pb
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

MainPath='E:/Padova/Experiments/8_Rhizotron_DAPHNAE_Lancaster/'
os.chdir(MainPath)
## ----------------- Define rhizotron configuration------------------------- ##
split = 1

## ----------------- IMPORT Scheme --------------------------------- ##
if split==1:
    orig_scheme = pb.load("./1a_CreateSequences/SequenceERT_Rhizo_71ElecsSplit.shm")
else:
    orig_scheme = pb.load("./1a_CreateSequences/SequenceERT_Rhizo_71Elecs.shm")


## ----------------- IMPORT MESH --------------------------------- ##

if split==1:
#    os.chdir('./1b_CheckMesh/Splitted/BERT/')
    mesh3d= mt.readGmsh('./1_Mesh/BaseRhizo_Splitted_Vrte.msh', verbose=True)
#    mesh3d= mt.readGmsh('./1b_CheckMesh/Splitted/BERT/BaseRhizo_Splitted_withVoid.msh', verbose=True)
    #mesh3d = pg.load("mesh_Splitted.bms", verbose=True) # if you wanna use your tet mesh instead
    #\1b_CheckMesh\Splitted\BERT
    mesh3d.save('MeshR_S.bms')
else:
    mesh3d= mt.readGmsh('./1_Mesh/BaseRhizo.msh', verbose=True)
#    mesh3d= mt.readGmsh('./1b_CheckMesh/BERT/BaseRhizo.msh', verbose=True)
    # OR
#    mesh3d = pg.load("./1b_CheckMesh/BERT/mesh_BaseRhizo.bms", verbose=True) 

print('Mesh: Nodes:', mesh3d.nodeCount(),
      'Cells:', mesh3d.cellCount(),
      'Boundaries:', mesh3d.boundaryCount())

mesh3d.exportVTK('mesh3d_Markers')
#pg.show(mesh3d, notebook=True)
cmap = plt.cm.get_cmap("viridis", 4)
mesh3dVTK = pv.read('mesh3d_Markers.vtk')
p = pv.Plotter(notebook=False)
p.add_mesh(mesh=mesh3dVTK,scalars='_Marker',cmap=cmap)    # add a dataset to the scene
p.add_bounding_box()
_ = p.show_bounds(grid='front', location='outer', all_edges=True)
p.show()

sensors=[]
sensorsNp=[]
for node in mesh3d.nodes():
    if node.marker() == -99:
        sensors.append(node.pos())
        sensorsNp.append(node.pos().array())
#scheme.setSensorPositions(sensors)
#scheme.set("valid", np.ones(len(meas)))
#scheme.save('SequenceERT_Rhizo_72Elecs.shm')
sensorsNp = np.vstack(sensorsNp)
height = 0.50
Effheight = max(sensorsNp[:,1])
width = 0.45 
thickness = 0.03
geom = mt.createRectangle(start=[0,0], end=[width, Effheight], marker=1,
                          boundaryMarker=0)
pg.show(geom, markers=True)
fig, ax = plt.subplots()
pg.show(geom, markers=True, ax=ax, hold=True)
drawSensors(ax, sensorsNp, diam=0.002)
##ax.figure.savefig(directoryFig +'Geom'+'.png')

# %%
## ----------------- SIMULATE ERT --------------------------------- ##

# Set refernce nodes in corners (necessary for closed geometries)
lower_left_node = mesh3d.findNearestNode([mesh3d.xmin(), mesh3d.ymin(), 0.0])
mesh3d.node(lower_left_node).setMarker(pg.MARKER_NODE_CALIBRATION)
#mesh.node(node_number).setMarker(-999)
#mesh.node(node_number2).setMarker(-1000)
lower_right_node = mesh3d.findNearestNode([mesh3d.xmax(), mesh3d.ymin(), 0.0])
mesh3d.node(lower_right_node).setMarker(pg.MARKER_NODE_REFERENCEELECTRODE)    
                    
# Initialize the ERTManager
ert = pb.Resistivity()
ert.setMesh(mesh3d, omitBackground=True) 

# Create meaningless resistivity model 
rho = []
for cell in mesh3d.cells():   
    if cell.center().x() < height / 3:
        rho.append(150)
    elif cell.center().x() > (2 * height / 3):
        rho.append(150)
    else:
        rho.append(150)
rho = np.array(rho)

# if splitted mesh fix values of split plane
if split==1:
    ert.fop.region(3).setFixValue(1e9)
#    ert.fop.region(2).setStartValue(100)
    #ert.fop.region(3).setStartValue(1e9)
    ert.fop.region(3).setConstraintType(0)
    ert.fop.region(3).setConstraintsWeight(10.)
    ert.fop.regionManager().setInterRegionConstraint(2, 3, 0.0) # 10 or 0?

    rhoS = []
    for cell in mesh3d.cells():   
        if cell.marker()==2:
            if cell.center().x() < 0.22:
                rhoS.append(10)
            elif cell.center().x() > 0.23:
                rhoS.append(500)
        if cell.marker()==3:
           rhoS.append(1e9)       
    rhoS = np.array(rhoS)   

mesh3d.addExportData("True model (Ohm*m)", rhoS)
mesh3d.exportVTK("mesh3d_TrueModel.vtk")
    
data = ert.simulate(mesh=mesh3d, res=rhoS, scheme=orig_scheme,
                    noiseLevel=0.02, noiseAbs=0, verbose=True)
data.save('DataS.dat')