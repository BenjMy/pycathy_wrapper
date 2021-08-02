#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:12:41 2021

@author: ben
"""

import pyvista as pv
import glob
import time 

def showvtk(filename=None,unit=None,timeStep=0):
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

    mesh = pv.read(filename)

    if unit in list(mesh.array_names):
        print('plot '+ str(unit))
    else:
        print('physcial property not existing')

    if filename is None:
        if unit is 'pressure':
            filename = '10' + timeStep + '.vtk'


    plotter = pv.Plotter()
    _ = plotter.add_mesh(mesh,show_edges=True)
    legend_entries = []
    legend_entries.append(['Time='+ str(mesh['TIME']), 'w'])
    _ = plotter.add_legend(legend_entries)
    plotter.show_grid()
    #plotter.add_mesh_clip_box(mesh, color='white')
    cpos = plotter.show()

    return

def showvtkTL(filename,unit=None,timeStep='all'):
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


    plotter = pv.Plotter()
    mesh = pv.read('./rhizo_prj/vtk/100.vtk')
    _ = plotter.add_mesh(mesh,show_edges=True)
    legend_entries = []
    legend_entries.append(['Time='+ str(mesh['TIME']), 'w'])
    plotter.show_grid()
    cpos = plotter.show(interactive_update=True, auto_close=False)


    # def TimeLapse():
    for filename in glob.glob('./rhizo_prj/vtk/10*.vtk'):
        print(filename)
        mesh = pv.read(filename)
        _ = plotter.add_mesh(mesh,show_edges=True)
        legend_entries = []
        legend_entries.append(['Time='+ str(mesh['TIME']), 'w'])
        plotter.show_grid()
        plotter.show_grid()
        plotter.update(force_redraw=True)
        time.sleep(1)

    cpos = plotter.show()

    return
