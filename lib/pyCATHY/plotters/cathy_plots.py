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
import matplotlib.pyplot as plt
import panel as pn
import numpy as np
import ipywidgets as widgets
import imageio
import natsort
from decimal import Decimal

# from pyvirtualdisplay import Display

from pyCATHY.importers import cathy_outputs as out_CT

def convert_time_units(t, x_units):

    xlabel = " (s)"
    if x_units == "days":
        xlabel = " (days)"
        t_new = [x / (24 * 60 * 60) for x in t]
        t_new = round(t_new[0], 2)
    if x_units == "hours":
        xlabel = " (h)"
        t_new = [x / (60 * 60) for x in t]
        t_new = round(t_new[0], 1)

    return xlabel, t_new


def showvtk(
    filename=None,
    unit=None,
    timeStep=0,
    notebook=False,
    path=None,
    savefig=False,
    **kwargs,
):
    """
    Short summary.

    Parameters ---------- filename : type Description of parameter `filename`. unit : str
    ['pressure', 'TIME', 'saturation', 'permeability', 'velocity'] timeStep : type Description of
    parameter `timeStep`.

    Returns ------- type Description of returned object.

    """

    if path is None:
        path = os.getcwd()

    if filename is None:
        if unit == "pressure":
            filename = "10" + str(timeStep) + ".vtk"
        if unit == "saturation":
            filename = "cele20" + str(timeStep) + ".vtk"

    mesh = pv.read(os.path.join(path, filename))

    if unit in list(mesh.array_names):
        print("plot " + str(unit))
    else:
        print("physcial property not existing")

    if notebook == True:
        # pn.extension('vtk')  # this needs to be at the top of each cell for some reason

        out = widgets.Output()

        def on_value_change(change):
            with out:

                PhysScalars = "pressure"
                time_step = timeStep
                # print(change['new'])
                if hasattr(change.owner, "options"):  # =='Phys. prop:':
                    PhysScalars = change.owner.options[change["new"] - 1]
                else:
                    time_step = change["new"]

                out.clear_output()
                mesh = pv.read("./my_cathy_prj/vtk/10" + str(time_step) + ".vtk")
                # display = Display(visible=0, size=(600, 400))
                # display.start()

                plotter = pv.Plotter(notebook=True)
                _ = plotter.add_mesh(mesh, scalars=PhysScalars[0])
                # plotter.show(True)

                legend_entries = []
                legend_entries.append(["Time=" + str(mesh["TIME"]), "w"])
                _ = plotter.add_legend(legend_entries)
                plotter.show_grid()
                # plotter.add_mesh_clip_box(mesh, color='white')
                cpos = plotter.show(True)

        slider = widgets.IntSlider(
            min=0, max=10, step=1, continuous_update=True, description="Time step #:",
        )
        # play = widgets.Play(min=1, interval=2000)
        choice = widgets.Dropdown(
            options=[("pressure", 1), ("saturation", 2)],
            value=1,
            description="Phys. prop:",
        )
        slider.observe(on_value_change, "value")
        choice.observe(on_value_change, "value")

        plotvtk = widgets.VBox([choice, slider, out])
        plotvtk

    else:

        plotter = pv.Plotter(notebook=True)
        _ = plotter.add_mesh(mesh, show_edges=True, scalars=unit)
        legend_entries = []
        legend_entries.append(["Time=" + str(mesh["TIME"]), "w"])
        _ = plotter.add_legend(legend_entries)
        plotter.show_grid()
        # plotter.add_mesh_clip_box(mesh, color='white')
        for key, value in kwargs.items():
            print(f"key: {key} | value: {value}")
            if key == "elecs":
                poly_elecs = pv.PolyData(value)
                poly_elecs["My Labels"] = [
                    f"Label {i}" for i in range(poly_elecs.n_points)
                ]
                plotter.add_point_labels(
                    poly_elecs, "My Labels", point_size=20, font_size=36
                )

        if savefig is True:
            # The supported formats are: ‘.svg’, ‘.eps’, ‘.ps’, ‘.pdf’, ‘.tex’
            plotter.save_graphic(
                filename + str(".pdf"), title="", raster=True, painter=True
            )

        cpos = plotter.show()

    return


def showvtkTL(
    filename=None,
    unit=None,
    timeStep="All",
    notebook=False,
    path=None,
    savefig=False,
    show=True,
    **kwargs,
):
    """
    Short summary.

    Parameters ---------- filename : type Description of parameter `filename`. unit : type
    Description of parameter `unit`. timeStep : type Description of parameter `timeStep`.

    Returns ------- type Description of returned object.

    """

    x_units = None
    xlabel = "s"
    for key, value in kwargs.items():
        if key == "x_units":
            x_units = value
            print(x_units)

    if path is None:
        path = os.getcwd()

    if filename is None:
        if unit == "pressure":
            filename = "10*.vtk"
            filename0 = "100.vtk"
        if unit == "saturation":
            filename = "cele20*.vtk"
            filename0 = "cele200.vtk"

    mesh = pv.read(os.path.join(path, filename0))

    if unit in list(mesh.array_names):
        print("plot " + str(unit))
    else:
        print("physcial property not existing")

    offscreen = False
    if show == False:
        offscreen = True

    plotter = pv.Plotter(notebook=False, off_screen=offscreen)
    plotter.add_mesh(mesh, show_edges=True)

    if savefig == True:
        plotter.open_gif(unit + ".gif")

    if x_units is not None:
        xlabel, t_lgd = convert_time_units(mesh["TIME"], x_units)
    legend_entry = "Time=" + str(t_lgd) + xlabel

    plotter.show_grid()
    # cpos = plotter.show(interactive_update=True, auto_close=False)
    cpos = plotter.show(auto_close=False)
    # plotter.add_legend(legend_entry)
    plotter.add_text(legend_entry, name="time-label")

    files = []
    for file in glob.glob(os.path.join(path, filename)):
        files.append(file)

    import natsort

    print(natsort.natsorted(files, reverse=False))

    for ff in natsort.natsorted(files, reverse=False):
        mesh = pv.read(ff)
        array_new = mesh.get_array(unit)

        if x_units is not None:
            xlabel, t_lgd = convert_time_units(mesh["TIME"], x_units)
        legend_entry = "Time=" + str(t_lgd) + xlabel
        plotter.update_scalars(array_new, render=True)
        plotter.add_text(legend_entry, name="time-label")

        plotter.render()
        plotter.write_frame()

    if savefig == True:
        # gif_original = filename + '.gif'
        # gif_speed_down = filename + 'new.gif'
        gif_original = unit + ".gif"
        gif_speed_down = unit + "_slow.gif"
        gif = imageio.mimread(gif_original)
        imageio.mimsave(gif_speed_down, gif, fps=0.8)

    plotter.close()

    return


def atmbc_inputs_plot(t_atmbc, v_atmbc, **kwargs):

    xlabel = "time (s)"
    for key, value in kwargs.items():
        if key == "x_units":
            x_units = value

    if x_units == "days":
        xlabel = "time (days)"
        t_atmbc = [x / (24 * 60 * 60) for x in t_atmbc]
    if x_units == "hours":
        xlabel = "time (h)"
        t_atmbc = [x / (60 * 60) for x in t_atmbc]

    print(t_atmbc)
    # https://matplotlib.org/stable/gallery/lines_bars_and_markers/stairs_demo.html#sphx-glr-gallery-lines-bars-and-markers-stairs-demo-py
    vdiff = v_atmbc[0] - v_atmbc[1]

    fig, ax = plt.subplots()
    ax.plot(t_atmbc, vdiff, "k*")
    ax.set(xlabel=xlabel, ylabel="Q (m/s)", title="atmbc inputs")
    ax.grid()

    plt.step(t_atmbc, v_atmbc[0], color="blue", where="post", label="Rain/Irr")
    plt.step(t_atmbc, v_atmbc[1], color="red", where="post", label="ET")
    plt.step(t_atmbc, vdiff, color="black", where="post", label="Diff")
    # ax.legend(['Rain/Irr','ET','diff'])
    ax.legend()

    pass


def indice_veg_plot(veg_map, **kwargs):
    '''
    View from top of the vegetation type (equivalent somehow to root map)

    Parameters
    ----------
    veg_map : np.array([])
        Indice of vegetation. The dimension of the vegetation map must match 
        the dimension of the DEM.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''

    fig, ax = plt.subplots()
    cf = ax.pcolormesh(veg_map, edgecolors="black")
    fig.colorbar(cf, ax=ax, label='indice of vegetation')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('view from top (before extruding)')
    plt.show()
    
    
    return fig, ax

def dem_plot_2d_top(df,parameter, **kwargs):
    '''
    View from top of the a given parameter

    Parameters
    ----------
    veg_map : np.array([])
        Indice of vegetation. The dimension of the vegetation map must match 
        the dimension of the DEM.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''

    fig, ax = plt.subplots()
    cf = ax.pcolormesh(veg_map, edgecolors="black")
    fig.colorbar(cf, ax=ax, label='indice of vegetation')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('view from top (before extruding)')
    plt.show()
    
    


def dem_plot(workdir, project_name, **kwargs):
    """
    DEM3D creates a 3D representation from a Grass DEM file

    """
    length = 1
    width = 1

    # Read the Header
    # str_hd_dem = {'north':0,'south':0,'east':0,'west':0,'rows':0,'cols':0}
    str_hd_dem = {}
    with open(
        os.path.join(workdir, project_name, "prepro/dem"), "r"
    ) as f:  # open the file for reading
        count = 0
        for line in f:  # iterate over each line
            if count < 6:
                str_hd, value_hd = line.split()  # split it by whitespace
                str_hd_dem[str_hd.replace(":", "")] = value_hd
            count += 1

    dem_file = open(os.path.join(workdir, project_name, "prepro/dem"), "r")
    dem_mat = np.loadtxt(dem_file, skiprows=6)
    dem_file.close()

    x = np.zeros(int(str_hd_dem["rows"]))
    y = np.zeros(int(str_hd_dem["cols"]))

    for a in range(int(str_hd_dem["rows"])):
        x[a] = float(str_hd_dem["west"]) + length * a

    for a in range(int(str_hd_dem["cols"])):
        y[a] = float(str_hd_dem["south"]) + width * a

    # x=x-width/2
    # y=y-length/2

    fig = plt.figure()
    ax = plt.axes(projection="3d")
    # Make data.
    X, Y = np.meshgrid(x, y)
    # Plot the surface.
    surf = ax.plot_surface(X, Y, dem_mat.T, cmap="viridis")
    # Add a color bar which maps values to colors.
    cbar = fig.colorbar(surf, shrink=0.25, orientation='horizontal',
                        label='Elevation (m)')

    ax.set(xlabel="Easting (m)", ylabel="Northing (m)", zlabel="Elevation (m)")
    plt.show()


def COCumflowvol(workdir, project_name):
    """
    Processes Cumflowvol.

    Returns ------- type Description of returned object.

    """
    
    CT_out.read_cumflowvol(os.path.join(workdir, project_name, "output", "cumflowvol"))
    


    fig, ax = plt.subplots()
    ax.plot(CUMFLOWVOL[:, 2], -CUMFLOWVOL[:, 7], "b-.")
    ax.plot(CUMFLOWVOL[:, 2] / 3600, CUMFLOWVOL[:, 7])
    ax.set_title("Cumulative flow volume")
    ax.set(xlabel="Time (s)", ylabel="Net Flow Volume (m^3)")
    ax.legend(["Total flow volume", "nansfdir flow volume"])

    return
