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

from matplotlib.colors import LogNorm

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

import matplotlib.font_manager
import matplotlib.style
import matplotlib as mpl
import matplotlib.dates as mdates

mpl.style.use('default')

mpl.rcParams['grid.color'] = 'k'
mpl.rcParams['grid.linestyle'] = ':'
mpl.rcParams['grid.linewidth'] = 0.25
# mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 0.75


# These parameters can also be put into the style or matplotlibrc.
# This is the dynamic approach of changing parameters.
nice_fonts = {
    "font.family": "serif",
    "font.serif" : "Times New Roman",
}
matplotlib.rcParams.update(nice_fonts)




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
        df_hgsfdeth = out_CT.read_hgsfdet(filename='hgsfdeth')


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

    if 'datetime' in kwargs:
        t_atmbc = kwargs['datetime']
        xlabel = "date"

        
    # https://matplotlib.org/stable/gallery/lines_bars_and_markers/stairs_demo.html#sphx-glr-gallery-lines-bars-and-markers-stairs-demo-py
    
    v_atmbc_n =[]
    if isinstance(v_atmbc, list):
        v_atmbc_p = v_atmbc[1] # positif
        v_atmbc_n = v_atmbc[0] # negatif
        v_atmbc = v_atmbc[0] - v_atmbc[1]

    

    



        

    # if np.shape(t_atmbc) != np.shape(v_atmbc):
    if len(v_atmbc_n)>0:
        fig, ax = plt.subplots(2,1, sharex=True)
        
        (ax1,ax2) = (ax[0], ax[1])
        
        color = 'tab:blue'
        ax1.set_xlabel('time (h)')
        ax1.set_ylabel('Rain/Irr', color=color)

        # ax1.step(t_atmbc, v_atmbc, color="green", where="post", label="net (diff)",marker='.')
        ax1.step(t_atmbc, v_atmbc_p, color="blue", where="post", label="Rain/Irr")
        
       
        color = 'tab:red'
        ax2.set_ylabel('ET', color=color)  # we already handled the x-label with ax1
        ax2.step(t_atmbc, -v_atmbc_n, color=color, where="post", label="ET")
        # ax2.bar(t_atmbc, -v_atmbc_n, color=color)
        
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show(block=False)
    else:
        # v_atmbc = 
        fig, ax = plt.subplots(figsize=(6,3))
        ax.plot(t_atmbc, v_atmbc, "k.")
        ax.set(xlabel=xlabel, ylabel="net Q (m/s)", title="atmbc inputs")
        ax.grid()
        plt.show(block=False)
        
        if 'IETO' in kwargs:
            if kwargs['IETO'] != 0:
                plt.step(t_atmbc, v_atmbc, color="k", where="post")
            elif kwargs['IETO'] == 0: # case of linear interpolation between points
                ax.plot(t_atmbc, v_atmbc, "k.")




        
    # ax.legend()
    
    # dateFormat = '%d-%H'
    if 'dateFormat' in kwargs:
        dateFormat = kwargs['dateFormat']
    
    if 'date_locator' in kwargs:
        interval = kwargs['date_locator'][1]
        if 'HourLocator' in kwargs['date_locator'][0]:
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=8))

        ax.xaxis.set_major_formatter(mdates.DateFormatter(dateFormat))
        plt.gcf().autofmt_xdate()
        
    
    # ax.set_yscale('symlog')
    plt.show(block=False)

    # if 'datetime' in kwargs:
    #     ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
    #     ax.xaxis.set_minor_locator(mdates.MonthLocator())
    #     ax.grid(True)

    return plt, ax 



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
        ['pressure', 'saturation', 'ER', 'permeability', 'velocity'] . The default is pressure.
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
        
        
        # for surface/subsurface hydrology
        # --------------------------------------------------------------------
        if unit == "pressure":
            my_colormap = 'autumn'

        elif unit == "saturation":
            my_colormap = 'Blues'



        if timeStep<10:    
            filename = "10" + str(timeStep) + ".vtk"
        elif timeStep<100:
            filename = "1" + str(timeStep) + ".vtk"
        elif timeStep<200:
            newnb = [int(x) for x in str(timeStep)]
            filename = "2" + str(newnb[1]) + str(newnb[2])  + ".vtk"
        elif timeStep<300:
            newnb = [int(x) for x in str(timeStep)]
            filename = "3" + str(newnb[1]) + str(newnb[2])  + ".vtk"
            
            
        
        # for transport
        # --------------------------------------------------------------------
        elif unit == "celerity":
            raise NotImplementedError('Transport vtk file output not yet implemented')

        
        
        # for ERT
        # --------------------------------------------------------------------
        elif 'ER' in unit:
            filename = "ER_converted" + str(timeStep) + ".vtk"
            my_colormap = 'viridis'
            unit = 'ER_converted' + str(timeStep)

        elif 'ER_int' in unit:
            filename = "ER_converted" + str(timeStep) + '_nearIntrp2_pg_msh.vtk'
            my_colormap = 'viridis'
            unit = 'ER_converted' + str(timeStep) + '_nearIntrp2_pg_msh'


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

        plotter = pv.Plotter(notebook=False)
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
        #plotter.show_grid()
        
        _ = plotter.show_bounds(minor_ticks=True,font_size=1)
        #plotter.add_mesh_clip_box(mesh, color='white')
        
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
        plotter.save_graphic(os.path.join(path,filename + unit + '.svg'),
                              title="", raster=True, painter=True)

        print('figure saved' + os.path.join(path,filename + '.svg'))
    cpos = plotter.show()

    return


def show_vtk_TL(filename=None,unit=None,timeStep="all", notebook=False, 
                path=None,savefig=False, show=True,**kwargs):
    '''
    Time lapse animation of selected time steps

    Parameters
    ----------
    filename : str, optional
        DESCRIPTION. The default is None.
    unit : str, optional
        pressure/saturation/ER. The default is None.
    timeStep : list, optional
        DESCRIPTION. The default is "all".
    notebook : Bool, optional
        pyvista notebook option. The default is False.
    path : str, optional
        path to the dir containing the vtk files. The default is None.
    savefig : Bool, optional
        save figure. The default is False.
    show : Bool, optional
        show plot. The default is True.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    

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
            filename = "1*.vtk"
            filename0 = "100.vtk"
            my_colormap = 'autumn'
        elif unit == "saturation":
            my_colormap = 'Blues'
            filename = "1*.vtk"
            filename0 = "100.vtk"
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

    plotter = pv.Plotter(notebook=notebook, off_screen=offscreen)
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


    # print(natsort.natsorted(files, reverse=False))

    for ff in natsort.natsorted(files, reverse=False):
        print(ff)
        mesh = pv.read(ff)
        array_new = mesh.get_array(unit)
        print(mesh["TIME"])

        if x_units is not None:
            xlabel, t_lgd = convert_time_units(mesh["TIME"], x_units)
            legend_entry = "Time=" + str(t_lgd) + xlabel
        plotter.update_scalars(array_new, render=True)
        plotter.add_text(legend_entry, name="time-label")

        plotter.render()
        plotter.write_frame()
        if unit == "saturation":
            plotter.update_scalar_bar_range([0,1])
            
            if 'clim' in kwargs:
                plotter.update_scalar_bar_range([kwargs['clim'][0],kwargs['clim'][1]])

            
            
    if savefig == True:
        # gif_original = filename + '.gif'
        # gif_speed_down = filename + 'new.gif'
        gif_original = os.path.join(path + unit + ".gif")
        gif_speed_down = os.path.join(path + unit + "_slow.gif")
        gif = imageio.get_reader(gif_original)
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

    # colors =  plt.cm.Vega20c( (4./3*np.arange(20*3/4)).astype(int) )
    
    cmap='tab10'
    if 'cmap' in kwargs: 
        cmap = kwargs['cmap']

    fig, ax = plt.subplots()
    cf = ax.pcolormesh(veg_map, edgecolors="black", cmap=cmap)
    fig.colorbar(cf, ax=ax, label='indice of vegetation')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('view from top (before extruding)')
    plt.show(block=False)
    plt.close()
    
    
    return fig, ax
  


def dem_plot_2d_top(parameter, label='', **kwargs):
    '''
    View from top of the a given parameter

    Parameters
    ----------
    parameter : np.array([]) or dict
        The dimension of the vegetation map must match 
        the dimension of the DEM.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''

    if label=='all':
        
        fig, axs = plt.subplots(int(len(parameter.keys())/2),
                                int(len(parameter.keys())/3),
                                sharex=True,sharey=True)
        # fig.set_title('view from top (before extruding)')

        for ax, p in zip(axs.reshape(-1),parameter.keys()): 
            cf = ax.pcolormesh(parameter[p]) # edgecolors="black"
            fig.colorbar(cf, ax=ax, label=p,fraction=0.046, pad=0.04, shrink=0.8)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_aspect('auto', 'box')

        plt.tight_layout()
        plt.show(block=False)
        
    else:
        fig, ax = plt.subplots()
        cf = ax.pcolormesh(parameter, edgecolors="black")
        fig.colorbar(cf, ax=ax, label=label)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('view from top (before extruding)')
        plt.show(block=False)

    return fig, ax 
    
    


def dem_plot(workdir, project_name, **kwargs):
    """
    DEM3D creates a 3D representation from a Grass DEM file

    """
    length = 1
    width = 1
    
    delta_x = 1 
    if 'delta_x' in kwargs:
        delta_x = kwargs['delta_x']
        
    
    delta_y = 1 
    if 'delta_y' in kwargs:
        delta_y = kwargs['delta_y']
        
        

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

    # dem_file = open(os.path.join(workdir, project_name, "prepro/dem"), "r")
    # dem_mat = np.loadtxt(dem_file, skiprows=6)
    # dem_file.close()

    dem_file = open(os.path.join(workdir, project_name, "prepro/dtm_13.val"), "r")
    dem_mat = np.loadtxt(dem_file, skiprows=0)
    dem_file.close()
    
    # x = np.zeros(int(str_hd_dem["rows"]))
    # y = np.zeros(int(str_hd_dem["cols"]))

    x = np.zeros(dem_mat.shape[0])
    y = np.zeros(dem_mat.shape[1])

    # for a in range(int(str_hd_dem["rows"])):
    #     x[a] = float(str_hd_dem["west"]) + length * a

    # for a in range(int(str_hd_dem["cols"])):
    #     y[a] = float(str_hd_dem["south"]) + width * a


    for a in range(dem_mat.shape[0]):
        x[a] = float(str_hd_dem["west"]) + delta_x * a

    for a in range(dem_mat.shape[1]):
        y[a] = float(str_hd_dem["south"]) + delta_y * a

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
    plt.show(block=False)
    plt.close()


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



def plot_mesh_bounds(mesh_bound_cond_df):


    m = np.array(['o','+'])
    
    mvalue = []
    alpha = []
    for bound_val in mesh_bound_cond_df['bound']:
        if bound_val == True:
            mvalue.append(1)
            alpha.append(1)
        else:
            mvalue.append(0)
            alpha.append(0.1)
    
            
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(mesh_bound_cond_df['x'], 
               mesh_bound_cond_df['y'], 
               mesh_bound_cond_df['z'], 
               c=mvalue)
                                            
    
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    
    plt.show(block=False)
          

    return fig, ax
   

# -----------------------------------------------------------------------------
# -----------------Plot results DATA ASSIMILATION------------------------------
# -----------------------------------------------------------------------------


def show_DA_process_ens(EnsembleX,Data,DataCov,dD,dAS,B,Analysis, 
                       savefig=False,**kwargs):
    

    label_sensor = 'raw data'
    if 'label_sensor' in kwargs:
        label_sensor = kwargs['label_sensor']
        
    fig = plt.figure(figsize=(12, 6), dpi=300)
    ax1 = fig.add_subplot(2,5,1)
    cax = ax1.matshow(EnsembleX, aspect='auto') #,
          #cmap=cm.rainbow, norm=colors.LogNorm())
    ax1.set_title('Prior')
    ax1.set_ylabel(r'$\psi$ params #')
    ax1.set_xlabel('Members #')
    cbar = fig.colorbar(cax, location='bottom')
    
    ax = fig.add_subplot(2,5,6)
    # cax = ax.matshow(np.cov(EnsembleX), 
    #                   aspect='auto',cmap='gray',
    #                   norm=LogNorm(vmin=np.matrix(EnsembleX).min(), 
    #                                vmax=np.matrix(EnsembleX).max())
    #                   )
    cax = ax.matshow(np.cov(EnsembleX), 
                      aspect='auto',cmap='gray',
                      norm=LogNorm()
                      )
    ax.set_title('cov(Prior)')
    ax.set_xlabel(r'$\psi$ params #')
    ax.set_ylabel(r'$\psi$ params #')
    # cbar = fig.colorbar(cax, location='bottom')
    cbar = fig.colorbar(cax, format="$%.1f$", location='bottom')

    ax.set_yticks([])  
    
    
    ax = fig.add_subplot(2,5,2)
    cax = ax.matshow(np.tile(Data,(np.shape(EnsembleX)[1],1)).T, 
                      aspect='auto')
    ax.set_title(label_sensor)
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
                      aspect='auto',cmap='gray_r')
                      # cmap=cm.rainbow, norm=colors.LogNorm())
                      # vmin=0, vmax=1e-29)
    ax.set_title('cov(meas)')
    ax.set_ylabel('Meas')
    ax.set_xlabel('Meas')
    cbar = fig.colorbar(cax, location='bottom')
    ax.set_yticks([])  
    
    
    ax = fig.add_subplot(2,5,3)
    cax = ax.matshow(dD, 
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
    ax.set_ylabel(r'$\psi$ params #')
    ax.set_xlabel('Members #')
    cbar = fig.colorbar(cax, location='bottom')
    ax.set_yticks([])  
    
    
    ax = fig.add_subplot(2,5,5,sharey=ax1)
    cax = ax.matshow(Analysis, 
                     aspect='auto')
    ax.set_title('Posterior')
    ax.set_ylabel(r'$\psi$ params #')
    ax.set_xlabel('Members #')
    cbar = fig.colorbar(cax, location='bottom')
    ax.set_yticks([])  
    
    plt.tight_layout()
    plt.show(block=False)
    savename = 'showDA_process_ens'
    if 'savename' in kwargs:
        savename = kwargs['savename']
        
    if savefig==True:
    
        fig.savefig(savename +'.png', dpi=300)

    plt.close()

    return fig, ax



def DA_RMS(df_performance,sensorName):
    
    
    header = ['time', 'ObsType','RMSE'+sensorName,'NMRMSE'+sensorName,'OL']   
    df_perf_plot = df_performance[header]   
    df_perf_plot['RMSE'+sensorName] = df_perf_plot['RMSE'+sensorName].astype('float64')
    df_perf_plot['NMRMSE'+sensorName] = df_perf_plot['NMRMSE'+sensorName].astype('float64')
    df_perf_plot.OL = df_perf_plot.OL.astype('str')
    
    # df_perf_plot['RMSE']
    df_perf_plot.dropna(inplace=True)
    
    # sensorName = ['swc','swc1','sw2']
    p0 = df_perf_plot.pivot(index='time',columns='OL', values='RMSE'+sensorName)
    # p0 = df_perf_plot.pivot(index='time',columns='OL', values='ObsType')
    p1 = df_perf_plot.pivot(index='time',columns='OL', values='NMRMSE'+sensorName)
    # p0.plot()
    
    fig, ax = plt.subplots(2,1)
    p0.plot(xlabel='time',ylabel='RMSE'+sensorName, ax=ax[0],style=['.-'])
    p1.plot(xlabel='time',ylabel='NMRMSE'+sensorName, ax=ax[1],style=['.-'])
    
    
    plt.savefig('performanceDA.png', dpi=300)
    
    return ax, plt


def DA_plot_parm_dynamic(parm = 'ks', 
                         dict_parm_pert={}, 
                         list_assimilation_times = [],
                         savefig=False, 
                         **kwargs):
    
    ensemble_size = len(dict_parm_pert[parm]['ini_perturbation'])
    # nb_times = len(df_DA['time'].unique())

    fig = plt.figure(figsize=(6, 3), dpi=150)
    ax = fig.add_subplot()
    ax.hist(dict_parm_pert[parm]['sampling'],
              ensemble_size, alpha=0.5, label='sampling')
    
    
    # bins=dict_parm_pert[parm]['ini_perturbation'].quantile([0,.05,0.1,0.15,0.20,0.25,0.3,0.35,0.40,0.45,0.5,0.55,0.6,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1]).to_list()

    ax.hist(dict_parm_pert[parm]['ini_perturbation'],
              ensemble_size, alpha=0.5, label='ini_perturbation')
    plt.legend(loc='upper right')
    plt.ylabel('Probability')
    # plt.axvline(x=dict_parm_pert[parm]['ZROOT_nominal'],linestyle='--', color='red')

    for nt in list_assimilation_times:
        try:
            ax.hist(dict_parm_pert[parm]['update_nb'+str(nt+1)],
                      ensemble_size, alpha=0.5, label='update nb' + str(nt+1))
        except:
            pass
        plt.legend(loc='upper right')
        plt.ylabel('Probability')
        

    if 'log' in kwargs:
        if kwargs['log']:
            plt.xscale('log')
        
    plt.show()
    
    return ax



def DA_plot_parm_dynamic_scatter(parm = 'ks', 
                                 dict_parm_pert={}, 
                                 list_assimilation_times = [],
                                 savefig=False,
                                 **kwargs):
    
    
    if 'ax' in kwargs:
        ax = kwargs['ax']
    else:
        fig = plt.figure(figsize=(6, 3), dpi=350)
        ax = fig.add_subplot()
        
    ensemble_size = len(dict_parm_pert[parm]['ini_perturbation'])



    mean_t = [np.mean(dict_parm_pert[parm]['ini_perturbation'])]
    mean_t_yaxis=np.mean(dict_parm_pert[parm]['ini_perturbation'])
    cloud_t = np.zeros([ensemble_size,len(list_assimilation_times)])
    cloud_t[:,0]= np.array(dict_parm_pert[parm]['ini_perturbation'])
    np.shape(cloud_t)
    
    try:
        for nt in list_assimilation_times[1:]:
            # print(np.shape(np.hstack(dict_parm_pert[parm]['update_nb'+str(int(nt))])))
            # mean_t.append(np.mean(dict_parm_pert[parm]['update_nb'+str(int(nt+1))]))
            print(int(nt))
            cloud_t[:,int(nt+1)]= np.hstack(dict_parm_pert[parm]['update_nb'+str(int(nt))])
    except:
        pass

    np.shape(cloud_t)
    

    # -------------------------------    
    

    
    dict_parm_new = {}
    name= []
    i = 0
    for k in dict_parm_pert[parm].keys():
        if 'upd' in k:
            dict_parm_new[parm+'_'+k]=np.hstack(dict_parm_pert[parm][k])
            name.append(str(i+1))
            i = i +1

            
    df = pd.DataFrame()      
    df = pd.DataFrame(data=dict_parm_new)
    df.index.name = 'Ensemble_nb'
    # df=df.reset_index()

    boxplot = df.boxplot()  
    
    
    boxplot.set_xticklabels(name, rotation=90)
    plt.ylabel(parm)
    plt.xlabel('assimilation time')
    
    if 'log' in kwargs:
        if kwargs['log']:
            plt.yscale('log')

    # -------------------------------    

    # fig = plt.figure(figsize=(6, 3), dpi=350)
    # ax = fig.add_subplot()
    
    # for nt in range(len(cloud_t)-1):
    #     plt.scatter(np.arange(0,len(cloud_t[0,:])),cloud_t[nt,:],label='ens_nb' + str(nt+1))
    
    
    # plt.ylabel(parm)
    # # plt.ylim([mean_t_yaxis-mean_t_yaxis/2,mean_t_yaxis+mean_t_yaxis/2])
    # plt.grid(which='minor')
    # plt.xlabel('assimilation time')
    # # plt.legend(loc='upper right')
    # plt.grid(visible=True, which='major', axis='both')

    # plt.show()
    # plt.savefig('update_parm.png', dpi=300)

    return ax
    
    

def DA_plot_time_dynamic(DA, 
                         state='psi', 
                         nodes_of_interest=[], 
                         savefig=False, 
                         **kwargs):
    
    if 'ax' in kwargs:
        ax = kwargs['ax']
    else:
        fig = plt.figure(figsize=(6, 3), dpi=350)
        ax = fig.add_subplot()
        
    isOL = DA.loc[DA['OL']==True]
    isENS = DA.loc[DA['OL']==False]
    
    isENS_time_Ens = isENS.set_index(['time','Ensemble_nb'])
    isOL_time_Ens = isOL.set_index(['time','Ensemble_nb'])
    
    if type(nodes_of_interest) != list:
        nodes_of_interest = [nodes_of_interest]
    # nodes_of_interest = [0]
    print(nodes_of_interest)
    # take nodes of interests
    
    NENS = int(max(DA['Ensemble_nb'].unique()))
    
    
    # key2plot = 'psi_bef_update'
    
    ylabel = 'pressure head (m)'
    key2plot = 'analysis'
    if 'sw' in state:
        key2plot = 'sw_bef_update_'
        ylabel = 'water saturation'


    
    # isENS_time_Ens.xs((1, 1), level=('time', 'Ensemble_nb'), axis=0)
    
    # -----------------------------#
    if len(nodes_of_interest)>0:
    
        if len(isOL)>0:
            isOL.insert(2, "idnode", np.tile(np.arange(int(len(isOL)/(max(isOL['time']+1)*NENS))),
                                             int(max(isOL['time']+1)*NENS)), True)
            select_isOL =isOL[isOL["idnode"].isin(nodes_of_interest)]
            select_isOL = select_isOL.set_index(['time','idnode'])
            select_isOL = select_isOL.reset_index()
            
            # mean, min and max of Open Loop
            # --------------------------------------------------------------------------
            ens_mean_isOL_time = (select_isOL.set_index(['time','Ensemble_nb','idnode'])
                                    .groupby(level=['time','idnode'])[key2plot]
                                    .mean()
                                    )
            ens_mean_isOL_time = ens_mean_isOL_time.reset_index(name='mean(ENS)_OL')
            
            ens_min_isOL_time = (select_isOL.set_index(['time','Ensemble_nb','idnode'])
                                    .groupby(level=['time','idnode'])[key2plot]
                                    .min()
                                    )
            ens_min_isOL_time = ens_min_isOL_time.reset_index(name='min(ENS)_OL')
            
            
            ens_max_isOL_time = (select_isOL.set_index(['time','Ensemble_nb','idnode'])
                                    .groupby(level=['time','idnode'])[key2plot]
                                    .max()
                                    )
            ens_max_isOL_time = ens_max_isOL_time.reset_index(name='max(ENS)_OL')

                                    
        
        if len(isENS)>0:

            if len(isOL)>0:

                isENS.insert(2, "idnode", np.tile(np.arange(int(len(isENS)/(max(isENS['time'])*(NENS)))),
                                          int(max(isENS['time']))*(NENS)), True)
            else:
            
                isENS.insert(2, "idnode", np.tile(np.arange(int(len(isENS)/(max(isENS['time'])*(NENS+1)))),
                                          int(max(isENS['time']))*(NENS+1)), True)
            
            # isENS.insert(2, "idnode", np.tile(np.arange(int(len(isENS)/(max(isENS['time'])*(NENS+1)))),
            #                                   int(max(isENS['time']))*(NENS)), True)

                
            
            # try:
                # isENS.insert(2, "idnode", np.tile(np.arange(int(len(isENS)/(max(isENS['time'])*(NENS+1)))),
                #                                   int(max(isENS['time']))*(NENS)), True)
            # except:
            #     pass
            
            # try:
            # isENS.insert(2, "idnode", np.tile(np.arange(int(len(isENS)/(max(isENS['time'])*(NENS+1)))),
            #                               int(max(isENS['time']))*(NENS+1)), True)
            # except:
            #     pass
                  
               
                        
            select_isENS =     isENS[isENS["idnode"].isin(nodes_of_interest)]
            select_isENS = select_isENS.set_index(['time','idnode'])
            select_isENS = select_isENS.reset_index()
                  
            # mean, min and max of Open Loop
            # --------------------------------------------------------------------------
            ens_mean_isENS_time = (select_isENS.set_index(['time','Ensemble_nb','idnode'])
                                    .groupby(level=['time','idnode'])[key2plot]
                                    .mean()
                                    )
            ens_mean_isENS_time = ens_mean_isENS_time.reset_index(name='mean(ENS)')
            
            ens_min_isENS_time = (select_isENS.set_index(['time','Ensemble_nb','idnode'])
                                    .groupby(level=['time','idnode'])[key2plot]
                                    .min()
                                    )
            ens_min_isENS_time = ens_min_isENS_time.reset_index(name='min(ENS)')
            
            
            ens_max_isENS_time = (select_isENS.set_index(['time','Ensemble_nb','idnode'])
                                    .groupby(level=['time','idnode'])[key2plot]
                                    .max()
                                    )
            ens_max_isENS_time = ens_max_isENS_time.reset_index(name='max(ENS)')
            
     
    
    else:
        # take the spatial average mean
        # -----------------------------#
        
        spatial_mean_isOL_time_ens = isOL_time_Ens.groupby(level=['time','Ensemble_nb'])['analysis'].mean()
        spatial_mean_isOL_time_ens = spatial_mean_isOL_time_ens.reset_index()
        
        spatial_mean_isENS_time_ens = isENS_time_Ens.groupby(level=['time','Ensemble_nb'])['analysis'].mean()
        spatial_mean_isENS_time_ens = spatial_mean_isENS_time_ens.reset_index()
        
        select_isOL = spatial_mean_isOL_time_ens
        select_isENS = spatial_mean_isENS_time_ens

    


    
    # Plot
    # --------------------------------------------------------------------------
    
    if len(isENS)>0:
        ens_mean_isENS_time.pivot(index="time", 
                            columns=['idnode'], 
                          # columns=['idnode'], 
                          values=["mean(ENS)"]).plot(ax=ax, style=['.-'],color='blue')
        
        ens_max_isENS_time.pivot(index="time", 
                            columns=['idnode'], 
                          # columns=['idnode'], 
                          values=["max(ENS)"]).plot(ax=ax,style=['.-'],color='blue')
               
        ens_min_isENS_time.pivot(index="time", 
                            columns=['idnode'], 
                          # columns=['idnode'], 
                          values=["min(ENS)"]).plot(ax=ax,style=['.-'],color='blue',
                                                        xlabel= '(assimilation) time - (h)',
                                                        ylabel='pressure head $\psi$ (m)')
                                                    
        lgd= ax.fill_between(ens_min_isENS_time['time'], 
                        ens_min_isENS_time['min(ENS)'], 
                        ens_max_isENS_time['max(ENS)'], 
                        alpha=0.2,
                        color='blue',
                        label='minmax DA')

                                                
    if len(isOL)>0:
        ens_mean_isOL_time.pivot(index="time", 
                          # columns=["Ensemble_nb",'idnode'], 
                          columns=['idnode'], 
                          values=["mean(ENS)_OL"]).plot(ax=ax, style=['-'],color='grey',label=False,
                                                       ylabel='pressure head $\psi$ (m)',) # water saturation (-)
        ens_min_isOL_time.pivot(index="time", 
                            columns=['idnode'], 
                          # columns=['idnode'], 
                          values=["min(ENS)_OL"]).plot(ax=ax,style=['--'],
                                                       color='grey',
                                                       label=False)
        
        ens_max_isOL_time.pivot(index="time", 
                            columns=['idnode'], 
                          # columns=['idnode'], 
                          values=["max(ENS)_OL"]).plot(ax=ax,style=['--'],
                                                       color='grey',
                                                       label=False)
    
        ax.fill_between(ens_mean_isOL_time['time'], 
                        ens_min_isOL_time['min(ENS)_OL'], 
                        ens_max_isOL_time['max(ENS)_OL'], 
                        alpha=0.2,
                        color='grey',
                        label='minmax OL')
        
    ax.set_xlabel('time (h)')
    # ax.set_ylabel('pressure head (m)')
    ax.set_ylabel(ylabel)

    savename = 'showDA_dynamic'
    if 'savename' in kwargs:
        savename = kwargs['savename']
        
    if savefig==True:
    
        plt.savefig(savename +'.png', dpi=300)

    
    return  ax, plt


def DA_plot_Archie(df_Archie,savefig=False,**kwargs):
    
    plt.scatter(df_Archie['sw'],df_Archie['ER_converted'])
    plt.xlabel('saturation')
    plt.ylabel('ER_converted')
    plt.show()


    if 'porosity' in kwargs:
        plt.scatter(df_Archie['sw']*kwargs['porosity'],
                    df_Archie['ER_converted'])
        plt.xlabel('swc')
        plt.ylabel('ER_converted')
    
    
    if savefig==True:
        plt.savefig('Archie.png', dpi=300)
        
    
    