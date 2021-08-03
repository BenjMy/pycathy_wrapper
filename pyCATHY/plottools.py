#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:12:41 2021

@author: ben
"""

import pyvista as pv
import glob
import time 
import os

def showvtk(filename=None,unit=None,timeStep=0,notebook=False,path=None):
    """Short summary.

    Parameters
    ----------
    filename : type
        Description of parameter `filename`.
    unit : str
        ['pressure', 'TIME', 'saturation', 'permeability', 'velocity']
    timeStep : type
        Description of parameter `timeStep`.

    Returns
    -------
    type
        Description of returned object.

    """


    if path is None:
       path = os.getcwd()
       
    if filename is None:
        if unit == 'pressure':
            filename = '10' + str(timeStep) + '.vtk'
        if unit == 'saturation':
            filename = 'cele20' + str(timeStep) + '.vtk'
            
    mesh = pv.read(os.path.join(path,filename))
    
    if unit in list(mesh.array_names):
        print('plot '+ str(unit))
    else:
        print('physcial property not existing')


    # plotter = pv.PlotterITK()

    plotter = pv.Plotter(notebook=notebook)
    _ = plotter.add_mesh(mesh,show_edges=True,scalars=unit)
    legend_entries = []
    legend_entries.append(['Time='+ str(mesh['TIME']), 'w'])
    _ = plotter.add_legend(legend_entries)
    plotter.show_grid()
    #plotter.add_mesh_clip_box(mesh, color='white')
    cpos = plotter.show()

    return

def showvtkTL(filename=None,unit=None,timeStep='All',notebook=False,path=None):
    """Short summary.

    Parameters
    ----------
    filename : type
        Description of parameter `filename`.
    unit : type
        Description of parameter `unit`.
    timeStep : type
        Description of parameter `timeStep`.

    Returns
    -------
    type
        Description of returned object.

    """

    if path is None:
       path = os.getcwd()
        
        
    if filename is None:
        if unit == 'pressure':
            filename = '10*.vtk'
            filename0 = '100.vtk'
        if unit == 'saturation':
            filename = 'cele20*.vtk'
            filename0 = 'cele200.vtk'

    mesh = pv.read(os.path.join(path,filename0))

    if unit in list(mesh.array_names):
        print('plot '+ str(unit))
    else:
        print('physcial property not existing')
        
        
    plotter = pv.Plotter(notebook=notebook)
    mesh = pv.read(filename0)
    _ = plotter.add_mesh(mesh,show_edges=True)
    legend_entries = []
    legend_entries.append(['Time='+ str(mesh['TIME']), 'w'])
    plotter.show_grid()
    cpos = plotter.show(interactive_update=True, auto_close=False)


    # def TimeLapse():
    for files in glob.glob(filename):
        mesh = pv.read(files)
        _ = plotter.add_mesh(mesh,show_edges=True,scalars=unit)
        legend_entries = []
        legend_entries.append(['Time='+ str(mesh['TIME']), 'w'])
        plotter.show_grid()
        plotter.show_grid()
        plotter.update(force_redraw=True)
        time.sleep(1)

    cpos = plotter.show()

    return
