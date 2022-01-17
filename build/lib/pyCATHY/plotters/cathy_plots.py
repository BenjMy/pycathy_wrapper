#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:12:41 2021

@author: ben
"""

import glob
import time
import os


import pyvista as pv
# pv.global_theme.background = 'white'
pv.set_plot_theme('document')

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.dates import DateFormatter


import panel as pn
import ipywidgets as widgets
import imageio
import natsort
from decimal import Decimal

import pandas as pd
import numpy as np

import natsort

# from pyvirtualdisplay import Display

from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.cathy_utils import label_units, transform2_time_delta, convert_time_units



def show_hgsfdet(df_hgsfdeth=[], workdir=[], project_name=[], **kwargs):
    '''
    plot hgsfdet 
    
    

    Parameters
    ----------
    workdir : TYPE
        DESCRIPTION.
    project_name : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    fig, ax 

    '''
    
    # read hgraph file if df_hgraph not existing
    # ------------------------------------------------------------------------
    if len(df_hgsfdeth)==0:
        df_dtcoupling = out_CT.read_hgsfdet(filename='hgsfdeth')


    # create fig axis if not existing
    # ------------------------------------------------------------------------
    # if kwargs['ax'] == False:
    #     fig, ax = plt.subplots()
    
    fig, ax = plt.subplots(2,1)
    
    


    if 'delta_t' in kwargs:
        df_hgsfdeth['time'] = pd.to_timedelta(df_hgsfdeth['time'],unit='s') 
        
    df_hgsfdeth.pivot_table(values='NET SEEPFACE VOL',index='time').plot(ax=ax[0],
                                                                 ylabel='NET SEEPFACE VOL',
                                                                 xlabel='time (s)')

    df_hgsfdeth.pivot_table(values='NET SEEPFACE FLX',index='time').plot(ax=ax[1],
                                                                 ylabel='NET SEEPFACE FLX',
                                                                 xlabel='time (s)')
    
    return fig, ax 


def show_hgraph(df_hgraph=[], workdir=[], project_name=[],
                x='time', y='SW', **kwargs):
    '''
    plot hgraph 
    
    

    Parameters
    ----------
    workdir : TYPE
        DESCRIPTION.
    project_name : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    fig, ax 

    '''
    
    # read hgraph file if df_hgraph not existing
    # ------------------------------------------------------------------------
    if len(df_hgraph)==0:
        df_dtcoupling = out_CT.read_dtcoupling(filename='dtcoupling')


    # create fig axis if not existing
    # ------------------------------------------------------------------------
    # if kwargs['ax'] == False:
    #     fig, ax = plt.subplots()

    nstep = len(df_dtcoupling['Atmpot-d'])
    fig, ax = plt.subplots()
    
    jmax=max(df_dtcoupling['Atmact-vf'])
    
    timeatm = [df_dtcoupling['Deltat'][0]]
    for i in np.arange(1,nstep+1):
        # timeatm[i]=timeatm[i-1]+df_dtcoupling[i-1][1];
        timeatm.append(timeatm[i-1]+df_dtcoupling['Deltat'][i-1])
    
    plt.step(timeatm, df_dtcoupling['Atmact-vf'], color="blue", where="post", label="test")
    
    
    plt.step(timeatm, df_dtcoupling['Atmact-v'], color="blue", where="post", label="test")

    return fig, ax 


def show_vp(df_vp=[], workdir=[], project_name=[], 
            index='time', x='time', y='SW', **kwargs):
    '''
    plot vp 
    
    

    Parameters
    ----------
    workdir : TYPE
        DESCRIPTION.
    project_name : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    fig, ax 

    '''
    
    # read vp file if df_atmbc not existing
    # ------------------------------------------------------------------------
    if len(df_vp)==0:
        df_vp = out_CT.read_vp(filename='vp')


    # create fig axis if not existing
    # ------------------------------------------------------------------------
    # if kwargs['ax'] == False:
    #     fig, ax = plt.subplots()


    node_uni = df_vp['node'].unique()

    # plot panda df
    # ------------------------------------------------------------------------
    fig, ax = plt.subplots(1,len(node_uni))

    if len(node_uni)<2:
        ax = [ax]
            
    colors = cm.Reds(np.linspace(0.2, 0.8, len(df_vp['str_nb'].unique())))

    for i, n in enumerate(node_uni):
                  
        # select the node and the y attribute
        df_vp_pivot = pd.pivot_table(df_vp,values=y,index=['time','node','str_nb'])
        df_vp_pivot_nodei = df_vp_pivot.xs(n, level="node")
        df_vp_pivot_nodei.unstack('str_nb').plot(ax=ax[i], 
                                                 title='node: ' + str(n),
                                                 color = colors)
        
        # time_form = DateFormatter("%d-%h")
        # ax[i].xaxis.set_major_formatter(time_form)



        ax[i].legend([]);
        
        label = label_units(y)
        if i==0:
            ax[i].set_ylabel(label);

    return fig, ax 


def show_hgraph_2(df_hgraph=[], workdir=[], project_name=[],
                x='time', y='SW', **kwargs):
    '''
    plot hgraph 
    
    

    Parameters
    ----------
    workdir : TYPE
        DESCRIPTION.
    project_name : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    fig, ax 

    '''
    
    # read hgraph file if df_hgraph not existing
    # ------------------------------------------------------------------------
    if len(df_hgraph)==0:
        df_hgraph = out_CT.read_hgraph(filename='hgraph')


    # create fig axis if not existing
    # ------------------------------------------------------------------------
    # if kwargs['ax'] == False:
    #     fig, ax = plt.subplots()

    fig, ax = plt.subplots()
    
    if 'delta_t' in kwargs:
        df_hgraph['time'] = pd.to_timedelta( df_hgraph['time'],unit='s') 
        
    df_hgraph.pivot_table(values='streamflow',index='time').plot(ax=ax,
                                                                 ylabel='streamflow ($m^3/s$)',
                                                                 xlabel='time (s)')


    return fig, ax 


def show_atmbc(t_atmbc, v_atmbc, **kwargs):
    '''
    Plot atmbc=f(time)

    Parameters
    ----------
    t_atmbc : np.array
        time where atmbc change.
    v_atmbc : list of 1 or 2 arrays (when available)
        v_atmbc[0] is the array of Rain/Irrigation change over the time;
        v_atmbc[1] is the array of ET demand over the time;
    **kwargs 

    Returns
    -------
    fig, ax

    '''

    # NOT YET IMPLEMETED
    # read atmbc file if df_atmbc not existing
    # ------------------------------------------------------------------------
    # df_atmbc = []
    # if len(df_atmbc)==0:
    #     df_atmbc = out_CT.read_atmbc(filename='vp')
        
        
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

    # https://matplotlib.org/stable/gallery/lines_bars_and_markers/stairs_demo.html#sphx-glr-gallery-lines-bars-and-markers-stairs-demo-py
    
    
    if np.shape(t_atmbc) != np.shape(v_atmbc):
        v_atmbc_p = v_atmbc[0] # positif
        v_atmbc_n = v_atmbc[1] # negatif
        v_atmbc = v_atmbc[0] - v_atmbc[1]

    

    fig, ax = plt.subplots()
    ax.plot(t_atmbc, v_atmbc, "k*")
    ax.set(xlabel=xlabel, ylabel="Q (m/s)", title="atmbc inputs")
    ax.grid()

    if 'IETO' in kwargs:
        if kwargs['IETO'] != 0:
            plt.step(t_atmbc, v_atmbc, color="k", where="post")
        elif kwargs['IETO'] == 0: # case of linear interpolation between points
            ax.plot(t_atmbc, v_atmbc, "k--")

        

    if len(v_atmbc_p)>0:
        plt.step(t_atmbc, v_atmbc_p, color="blue", where="post", label="Rain/Irr")
        plt.step(t_atmbc, v_atmbc_n, color="red", where="post", label="ET")
        plt.step(t_atmbc, v_atmbc, color="black", where="post", label="Diff")
    ax.legend()

    return fig, ax 



def show_atmbc_3d(df_atmbc):
    '''
    Temporary (must exist throught show_vtk only)

    Parameters
    ----------
    df_atmbc : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    
    # df_atmbc.set_index('time',inplace=True)
    # df_atmbc.loc[0]

    # filename = './vtk/100.vtk'
    # mesh = pv.read(filename)
    # mesh.add_field_data(df_atmbc.loc[0]['value'], 'atmbc')

    # plotter = pv.Plotter(notebook=True)
    # _ = plotter.add_mesh(mesh, show_edges=True, scalars='atmbc')

    # # _ = plotter.add_legend(legend_entries)
    # plotter.show_grid()
    # cpos = plotter.show()
    
    show_vtk(new_field='atmbc')
    
    pass
                
                


def show_vtk(filename=None,unit='pressure',timeStep=0,notebook=False,path=None,
             savefig=False,**kwargs):
    '''
    Plot pyvista vtk file.


    Parameters
    ----------
    filename : TYPE, optional
        DESCRIPTION. The default is None.
    unit : str, optional
        ['pressure', 'TIME', 'saturation', 'permeability', 'velocity'] . The default is pressure.
    timeStep : TYPE, optional
        DESCRIPTION. The default is 0.
    notebook : TYPE, optional
        DESCRIPTION. The default is False.
    path : TYPE, optional
        DESCRIPTION. The default is None.
    savefig : TYPE, optional
        DESCRIPTION. The default is False.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''


    my_colormap = 'viridis'
    
    if path is None:
        path = os.getcwd()


    # Parse physical attribute + cmap from unit
    # -------------------------------------------------------------------------
    if filename is None:
        if unit == "pressure":
            filename = "10" + str(timeStep) + ".vtk"
            my_colormap = 'autumn'

        elif unit == "saturation":
            filename = "cele20" + str(timeStep) + ".vtk"
            my_colormap = 'Blues'

        elif 'ER' in unit:
            filename = "ER" + str(timeStep) + ".vtk"
            my_colormap = 'viridis'



    mesh = pv.read(os.path.join(path, filename))



    if unit in list(mesh.array_names):
        print("plot " + str(unit))
    else:
        print("physcial property not existing")



    # notebook = activate widgets 
    # -------------------------------------------------------------------------
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
                _ = plotter.add_mesh(mesh, scalars=PhysScalars[0], cmap=my_colormap)
                # plotter.show(True)

                if unit == "saturation":
                    plotter.update_scalar_bar_range([0,1])
                
                if 'clim' in kwargs:
                    plotter.update_scalar_bar_range([kwargs['clim'][0],1])

                    print(kwargs['clim'][0])
            
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


    # No notebook interaction
    # -------------------------------------------------------------------------
    else:

        plotter = pv.Plotter(notebook=True)
        _ = plotter.add_mesh(mesh, show_edges=True, scalars=unit, cmap=my_colormap)
        
        if unit == "saturation":
            plotter.update_scalar_bar_range([0,1])
        
        # options to colorbar
        # ---------------------------------------------------------------------
        if 'clim' in kwargs:
            plotter.update_scalar_bar_range([kwargs['clim'][0],kwargs['clim'][1]])

        # add time stamp as legend
        # ---------------------------------------------------------------------       
        legend_entries = []
        
        time_delta = transform2_time_delta(mesh["TIME"],'s')
        # print(mesh["TIME"])
        
        # legend_entries.append(["Time=" + str(mesh["TIME"]), "w"])
        legend_entries.append(["Time=" + str(time_delta[0]), "w"])

        
        _ = plotter.add_legend(legend_entries)
        plotter.show_grid()
        
        # plotter.add_mesh_clip_box(mesh, color='white')
        
        # add scatter points to the plot 
        # --------------------------------------------------------------------- 
        for key, value in kwargs.items():


            # add electrodes positions
            # ----------------------------------------------------------------- 
            if key == "elecs":
                poly_elecs = pv.PolyData(value)
                poly_elecs["My Labels"] = [
                    f"Label {i}" for i in range(poly_elecs.n_points)
                ]
                plotter.add_point_labels(
                    poly_elecs, "My Labels", point_size=20, font_size=36
                )

            # add tensiometers positions
            # ----------------------------------------------------------------- 
            
            # add TDR probe positions
            # ----------------------------------------------------------------- 
            
            
    # savefig 
    # --------------------------------------------------------------------- 
    if savefig is True:
        # The supported formats are: ‘.svg’, ‘.eps’, ‘.ps’, ‘.pdf’, ‘.tex’
        # print(os.path.join(path, 'vtk', filename + '.svg'))
        plotter.view_xz()
        plotter.save_graphic(os.path.join(path,filename + '.svg'),
                             title="", raster=True, painter=True)

        print('figure saved' + os.path.join(path,filename + '.svg'))
    cpos = plotter.show()

    return


def show_vtk_TL(
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
            my_colormap = 'autumn'
        elif unit == "saturation":
            my_colormap = 'Blues'
            filename = "cele20*.vtk"
            filename0 = "cele200.vtk"
        elif 'ER' in unit:
            filename = "ER" + str(timeStep) + ".vtk"
            my_colormap = 'viridis'
            
            
    mesh = pv.read(os.path.join(path, filename0))

    if unit in list(mesh.array_names):
        print("plot " + str(unit))
    else:
        print("physcial property not existing")

    offscreen = False
    if show == False:
        offscreen = True

    plotter = pv.Plotter(notebook=False, off_screen=offscreen)
    plotter.add_mesh(mesh, show_edges=True,cmap=my_colormap)

    if savefig == True:
        plotter.open_gif(os.path.join(path + unit + ".gif"))


    # options to colorbar
    # ---------------------------------------------------------------------
    if 'clim' in kwargs:
        plotter.update_scalar_bar_range([kwargs['clim'][0],kwargs['clim'][1]])

    legend_entry = "Time= " +str(mesh["TIME"])
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
        gif_original = os.path.join(path + unit + ".gif")
        gif_speed_down = os.path.join(path + unit + "_slow.gif")
        gif = imageio.mimread(gif_original)
        imageio.mimsave(gif_speed_down, gif, fps=0.8)
        print('gif saved' + os.path.join(path,gif_original))

    plotter.close()

    return



def indice_veg_plot(veg_map, **kwargs):
    '''
    View from top of the vegetation type (equivalent somehow to root map)

    Parameters
    ----------
    veg_map : np.array([])
        Indice of vegetation. The dimension of the vegetation map must match 
        the dimension of the DEM.
    **kwargs

    Returns
    -------
    fig, ax

    '''

    fig, ax = plt.subplots()
    cf = ax.pcolormesh(veg_map, edgecolors="black")
    fig.colorbar(cf, ax=ax, label='indice of vegetation')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('view from top (before extruding)')
    plt.show()
    
    
    return fig, ax

def dem_plot_2d_top(parameter, **kwargs):
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
    
    return fig, ax 
    
    


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
    '''
    plot COCumflowvol

    Parameters
    ----------
    workdir : TYPE
        DESCRIPTION.
    project_name : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''

    
    # read file if df_cumflowvol not existing
    # ------------------------------------------------------------------------
    df_cumflowvol = []
    if len(df_cumflowvol)==0:
        df_cumflowvol = out_CT.read_cumflowvol(os.path.join(workdir, project_name, 
                                                         "output", "cumflowvol"))
    

    # plot Net Flow Volume (m^3) = f(time)
    # ------------------------------------------------------------------------
    fig, ax = plt.subplots()
    ax.plot(df_cumflowvol[:, 2], -df_cumflowvol[:, 7], "b-.")
    ax.plot(df_cumflowvol[:, 2] / 3600, df_cumflowvol[:, 7])
    ax.set_title("Cumulative flow volume")
    ax.set(xlabel="Time (s)", ylabel="Net Flow Volume (m^3)")
    ax.legend(["Total flow volume", "nansfdir flow volume"])

    return fig, ax







# def mesh3d():


#     % MESH3D creates a 3D representation of the 3D grid
#     close all
#     clear all
#     load cape

#     fgrid = fopen('/Users/campo/Work/papers/transcathy-num-diff/cathy/bea/output/grid3d','r');
#     TETRA=[];
#     NODES=[];
#     NNOD=0;
#     N=0;
#     NT=0;

#     A = fscanf(fgrid,'%u %u %u',3);
#     NNOD = A(1);
#     N = A(2);
#     NT = A(3);
#     TETRA = fscanf(fgrid,'%u',[5,NT]); % Read data
#     TETRA = TETRA';
#     TETRA = TETRA(:,1:4);
#     NODES = fscanf(fgrid,'%g',[3,N]);
#     NODES = NODES';
#     % tetramesh(TETRA,NODES);%,'CData',NODES(:,3));
#     % % colormap(map)
#     % axis image;
#     % set(gca,'FontSize',14)
#     % xlabel('Easting (m)')
#     % ylabel('Northing (m)')
#     % zlabel('Elevation (m)')
#     clear S

#     S.Vertices = NODES;
#     S.Faces = TETRA;
#     S.FaceVertexCData = NODES(:,3);
#     S.FaceColor = 'interp';
#     % S.FaceAlpha = 0.5;
#     % S.LineStyle = ':';
#     S.LineWidth = 0.25;
#     % S.EdgeColor = 'red';
#     % S.LineWidth = 2;
#     figure
#     patch(S)
#     axis image
#     view(3)
#     set(gca,'FontSize',14);
#     xlabel('Easting (m)');
#     ylabel('Northing (m)');
#     zlabel('Elevation (m)');
#     % colormap(map);
#     colorbar;

# return 



# def dtcoupling(self):
#     """
#     processes the file dtcoupling and compares potential and actuale.

#     Returns ------- type Description of returned object.

#     """
#     dtcoupling_file = open(os.path.join(self.workdir ,'output' ,'dtcoupling'), 'r')
#     Lines = dtcoupling_file.readlines()
#     count = len(Lines)
#     dtcoupling_file.close()
    
#     # Using readline()
#     dtcoupling_file = open(os.path.join(workdir ,'output' ,'dtcoupling'), 'r')
#     nstep = count-31 # Number of timesteps
#     # DT = np.(file,'%g',[22,nstep]); % Read data
#     DT = np.loadtxt(dtcoupling_file,skiprows=22,max_rows=22+nstep)
#     dtcoupling_file.close()
    
#     DT[-1]
    
#     fig, axs = plt.subplots(3, 2)
    
#     axs[0,0].plot(DT[:,2],DT[:,8],'k:')
#     axs[0,0].set(xlabel='Time (h)', ylabel='Pot. & act. atm. fluxes (m^3/h)')
#     # axs[0,0].set_xlim([0,max(time)])
#     # axs[0,0].set_ylim([1000*(min(min(DT(:,11),DT(:,15)))),0])
#     axs[0,0].plot(DT[:,2],DT[:,13],'k-')
#     axs[0,0].legend(['Potential','Actual'])
    
#     axs[1,0].plot(DT[:,2],DT[:,16],'k:')
#     axs[1,0].plot(DT[:,2],DT[:,17],'k--')
#     axs[1,0].plot(DT[:,2],DT[:,18],'k-.')
#     axs[1,0].plot(DT[:,2],DT[:,19],'k-')
#     axs[1,0].set(xlabel='Time (h)', ylabel='Surface saturation fractions')
#     axs[1,0].legend(['Horton','Dunne','Ponded','Saturated'])
    
#     axs[2,0].plot(DT[:,2],DT[:,6],'k-')
#     axs[2,0].plot(DT[:,2],DT[:,7],'k-')
#     axs[2,0].set(xlabel='Time (h)', ylabel='Surface routing time steps')
#     axs[2,0].legend(['No backstep','Backstep'])
    
#     axs[0,1].plot(DT[:,2],DT[:,21]/DT[:,20],'k-')
#     # axs[1,0].plot(DT[:,2],DT[:,21],'k-')
#     axs[0,1].set(xlabel='Time (h)', ylabel='Surface/subsurface CPU')
    
#     axs[1,1].plot(DT[:,2],DT[:,1],'k-')
#     # axs[1,0].plot(DT[:,2],DT[:,21],'k-')
#     axs[1,1].set(xlabel='Time (h)', ylabel='Subsurface step size (h)')
    
#     axs[2,1].plot(DT[:,2],DT[:,4],'k-')
#     axs[2,1].plot(DT[:,2],DT[:,5],'k:')
#     axs[2,1].legend(['No backstep','Backstep'])
#     axs[2,1].set(xlabel='Time (h)', ylabel='Nonlinear iterations')



#     return


def show_DA_process_ens(EnsembleX,Data,DataCov,dD,dAS,B,Analysis, 
                       savefig=False,**kwargs):
    

    fig = plt.figure(figsize=(8, 6), dpi=300)
    ax1 = fig.add_subplot(2,5,1)
    cax = ax1.matshow(EnsembleX, aspect='auto') #,
          #cmap=cm.rainbow, norm=colors.LogNorm())
    ax1.set_title('Prior')
    ax1.set_ylabel('Parameters #')
    ax1.set_xlabel('Members #')
    cbar = fig.colorbar(cax, location='bottom')
    
    ax = fig.add_subplot(2,5,6)
    cax = ax.matshow(np.cov(EnsembleX), 
                      aspect='auto',cmap='gray')
    ax.set_title('cov(Prior)')
    ax.set_xlabel('Parameters #')
    ax.set_ylabel('Parameters #')
    cbar = fig.colorbar(cax, location='bottom')
    ax.set_yticks([])  
    
    
    ax = fig.add_subplot(2,5,2)
    cax = ax.matshow(np.tile(Data,(np.shape(EnsembleX)[1],1)).T, 
                      aspect='auto')
    ax.set_title('App. Res')
    ax.set_ylabel('Meas')
    ax.set_xlabel('Members #')
    cbar = fig.colorbar(cax, location='bottom')
    ax.set_yticks([])  
    
    # DataCov = np.diag(key_value[1]['data']['recipError'].to_numpy() )
    # DataCov.max()
    # DataCov.mean()
    # DataCov.min()
    ax = fig.add_subplot(2,5,7)
    cax = ax.matshow(DataCov, 
                      aspect='auto',cmap='gray')
                      # cmap=cm.rainbow, norm=colors.LogNorm())
                      # vmin=0, vmax=1e-29)
    ax.set_title('cov(meas)')
    ax.set_ylabel('Meas')
    ax.set_xlabel('Meas')
    cbar = fig.colorbar(cax, location='bottom')
    ax.set_yticks([])  
    
    
    ax = fig.add_subplot(2,5,3)
    cax = ax.matshow(dD.T, 
                     aspect='auto',
                     cmap='jet')
    ax.set_title('Meas - Sim')
    ax.set_ylabel('Meas')
    ax.set_xlabel('Members #')
    cbar = fig.colorbar(cax, location='bottom')
    ax.set_yticks([])  
    
    
    ax = fig.add_subplot(2,5,4,sharey=ax1)
    cax = ax.matshow(np.dot(dAS,B), 
                      aspect='auto',
                      cmap='jet')
    ax.set_title('Correction')
    ax1.set_ylabel('Parameters #')
    ax.set_xlabel('Members #')
    cbar = fig.colorbar(cax, location='bottom')
    ax.set_yticks([])  
    
    
    ax = fig.add_subplot(2,5,5,sharey=ax1)
    cax = ax.matshow(Analysis, 
                     aspect='auto')
    ax.set_title('Posterior')
    ax1.set_ylabel('Parameters #')
    ax.set_xlabel('Members #')
    cbar = fig.colorbar(cax, location='bottom')
    ax.set_yticks([])  
    
    
    savename = 'showDA_process_ens'
    if 'savename' in kwargs:
        savename = kwargs['savename']
        
    if savefig==True:
    
        plt.savefig(savename +'.png', dpi=300)

    
    return fig, ax
    


def DA_plot_time_dynamic(nodes_of_interest, DA, savefig=False, **kwargs):
        
    isOL = DA.loc[DA['OL']==True]
    isENS = DA.loc[DA['OL']==False]
    
    isENS_time_Ens = isENS.set_index(['time','Ensemble_nb'])
    isOL_time_Ens = isOL.set_index(['time','Ensemble_nb'])
    
    
    # isENS_time_Ens.take([-1])
    # simu_DA.grid3d
    nodes_of_interest = [0]
    print(nodes_of_interest)
    # take nodes of interests
    
    NENS = np.shape(DA)[1]
    # -----------------------------#
    if nodes_of_interest:
    
        isOL.insert(2, "idnode", np.tile(np.arange(int(len(isOL)/(max(isOL['time'])+1))),int(max(isOL['time']+1))), True)
        select_isOL =isOL[isOL["idnode"].isin(nodes_of_interest)]
        select_isOL = select_isOL.set_index(['time','idnode'])
        select_isOL = select_isOL.reset_index()
    
        # valuei_time_ens = DA_time_Ens.groupby(level=['time','Ensemble_nb'])['aft_update_'].take([0])
        # is_OL =  mean_time_ens['Ensemble_nb']==2
        
        
        # nodes_isOL_time_ens = isOL_time_Ens.groupby(level=['time','Ensemble_nb'])['aft_update_']
        # nodes_isOL_time_ens
        
        
        # nodes_isOL_time_ens = (nodes_isOL_time_ens.reset_index(level=1,drop=True)
        #                                             .reset_index(level=1,drop=False,name="idx")
        #                                             .reset_index(name="aft_update_")
        #                                             )
                               
                                
        isENS.insert(2, "idnode", np.tile(np.arange(int(len(isENS)/(max(isENS['time'])*NENS))),int(max(isENS['time']))*NENS), True)
        select_isENS =     isENS[isENS["idnode"].isin(nodes_of_interest)]
    
        select_isENS = select_isENS.set_index(['time','idnode'])
        # nodes_isENS_time_ens = nodes_isENS_time_ens.groupby(level=['time','idnode'])['aft_update_']
        select_isENS = select_isENS.reset_index()
    
                    
        # nodes_isENS_time_ens = isENS_time_Ens.groupby(level=['time','Ensemble_nb'])['aft_update_'].take(nodes_of_interest)
        # nodes_isENS_time_ens = (nodes_isENS_time_ens.reset_index(level=1,drop=True)
        #                                             .reset_index(level=1,drop=True)
        #                                             .reset_index(name="aft_update_")
        #                         )                                           
        # nodes_isENS_time_ens = nodes_isENS_time_ens.reset_index()
            
    
    
    else:
        # take the spatial average mean
        # -----------------------------#
        
        spatial_mean_isOL_time_ens = isOL_time_Ens.groupby(level=['time','Ensemble_nb'])['aft_update_'].mean()
        spatial_mean_isOL_time_ens = spatial_mean_isOL_time_ens.reset_index()
        
        spatial_mean_isENS_time_ens = isENS_time_Ens.groupby(level=['time','Ensemble_nb'])['aft_update_'].mean()
        spatial_mean_isENS_time_ens = spatial_mean_isENS_time_ens.reset_index()
        
        select_isOL = spatial_mean_isOL_time_ens
        select_isENS = spatial_mean_isENS_time_ens
    
    
    # take the ensemble mean
    # -----------------------------#
    ens_mean_isENS_time = (select_isENS.set_index(['time','Ensemble_nb','idnode'])
                            .groupby(level=['time','idnode'])['aft_update_']
                            .mean()
                            )
    ens_mean_isENS_time = ens_mean_isENS_time.reset_index(name='mean(ENS)')
    
    ens_min_isENS_time = (select_isENS.set_index(['time','Ensemble_nb','idnode'])
                            .groupby(level=['time','idnode'])['aft_update_']
                            .min()
                            )
    ens_min_isENS_time = ens_min_isENS_time.reset_index(name='min(ENS)')
    
    
    ens_max_isENS_time = (select_isENS.set_index(['time','Ensemble_nb','idnode'])
                            .groupby(level=['time','idnode'])['aft_update_']
                            .max()
                            )
    ens_max_isENS_time = ens_max_isENS_time.reset_index(name='max(ENS)')
    
    
    
    # ens_mean_isENS_time.rename({'aft_update_': 'mean(ENS)'},inplace=True)
    select_isOL.rename({'aft_update_': 'open_loop'}, axis=1, inplace=True)
    # ens_mean_isENS_time.rename({'aft_update_': 'mean(ENS)'}, axis=1, inplace=True)
    
    
    ax = select_isOL.pivot(index="time", 
                      # columns=["Ensemble_nb",'idnode'], 
                      columns=['idnode'], 
                      values=["open_loop"]).plot(style=['-'],color='red',label='open loop',
                                                    ylabel='pressure head $\psi$ (m)',) # water saturation (-)
    
    ens_mean_isENS_time.pivot(index="time", 
                        columns=['idnode'], 
                      # columns=['idnode'], 
                      values=["mean(ENS)"]).plot(ax=ax,style=['.-'],color='grey')
    
    ens_max_isENS_time.pivot(index="time", 
                        columns=['idnode'], 
                      # columns=['idnode'], 
                      values=["max(ENS)"]).plot(ax=ax,style=['.-'],color='magenta')
    
    ens_min_isENS_time.pivot(index="time", 
                        columns=['idnode'], 
                      # columns=['idnode'], 
                      values=["min(ENS)"]).plot(ax=ax,style=['.-'],color='blue',
                                                    xlabel= '(assimilation) time - (h)',
                                                    ylabel='pressure head $\psi$ (m)')

    savename = 'showDA_dynamic'
    if 'savename' in kwargs:
        savename = kwargs['savename']
        
    if savefig==True:
    
        plt.savefig(savename +'.png', dpi=300)

    
    return  ax