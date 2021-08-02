#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:12:41 2021

@author: ben
"""

import pyvista as pv

def showvtk(filename,unit=None,timeStep=0):
    
    mesh = pv.read(filename)
    plotter = pv.Plotter()
    plotter.add_mesh(mesh,show_edges=True)
    plotter.show_grid()    
    #plotter.add_mesh_clip_box(mesh, color='white')    
    cpos = plotter.show()
    
    return 

def showvtkTL(filename,unit=None):


    plotter = pv.Plotter()
    mesh = pv.read('./rhizo_prj/vtk/100.vtk')
    plotter.add_mesh(mesh,show_edges=True)
    plotter.show_grid()
    cpos = plotter.show(interactive_update=True, auto_close=False)
    
    
    # def TimeLapse():
    for filename in glob.glob('./rhizo_prj/vtk/10*.vtk'):
        print(filename)
        mesh = pv.read(filename)
        plotter.add_mesh(mesh,show_edges=True)
        plotter.show_grid()
        plotter.update(force_redraw=True)
        time.sleep(1)
    
    cpos = plotter.show()
    
    return


