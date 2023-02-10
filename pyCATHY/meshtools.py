#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Meshing tools
"""

import os

import numpy as np
import pyvista as pv

from pyCATHY.plotters import cathy_plots as cplt

try:
    import pygimli as pg
    import pygimli.meshtools as mt
except ImportError:
    pygimli = None

import matplotlib.pylab as plt


def CATHY_2_Simpeg(
    mesh_CATHY, ERT_meta_dict, scalar="saturation", show=False, **kwargs
):
    pass


def CATHY_2_pg(mesh_CATHY, ERT_meta_dict, scalar="saturation", show=False, **kwargs):
    """
    Interpolate CATHY mesh attribute to pygimli mesh.
    Add a new [`scalar`] attribute to the pygimli mesh (create a new mesh)

    .. Note:
        Need to flip axis because convention for CATHY and pygimli are different

    Parameters
    ----------
    mesh_CATHY : pvmesh
        CATHY mesh to transform to pygimli.
    ERT_meta_dict :dict
        Dictionnary containing ERT metadata (mesh, format, ..).
    scalar : str, optional
        scalar attribute to interpolate. The default is 'saturation'.
    show : bool, optional
        show the result of the interpolation using pyvista. The default is False.
    **kwargs : TYPE
        path : path of the mesh to overwrite

    Returns
    -------
    mesh_new_attr : TYPE
        DESCRIPTION.
    scalar_new : TYPE
        DESCRIPTION.

    """

    if type(ERT_meta_dict["forward_mesh_vtk_file"]) is str:
        mesh_OUT = pv.read(ERT_meta_dict["forward_mesh_vtk_file"])

    # flip y and z axis as CATHY and pg have different convention for axis
    # ------------------------------------------------------------------------
    in_nodes_mod = np.array(mesh_CATHY.points)
    in_nodes_mod_pg = np.array(mesh_OUT.points)

    idx = np.array([0, 2, 1])  # reorder xyz to xzy
    in_nodes_mod_m = in_nodes_mod[:, idx]
    # in_nodes_mod_m = in_nodes_mod[:, idx]
    in_nodes_mod_m[:, 2] = -np.flipud(in_nodes_mod_m[:, 2])
    in_nodes_mod_m[:, 1] = -np.flipud(in_nodes_mod_m[:, 1])

    # for i, axesi in enumerate(['x','y','z']):
    #     print('check {} consistency between the 2 meshes'.format(axesi))
    #     print(max(in_nodes_mod_m[:,i]),max(in_nodes_mod_pg[:,i]))
    #     print(min(in_nodes_mod_m[:,i]),min(in_nodes_mod_pg[:,i]))

    if "mesh_nodes_modif" in ERT_meta_dict.keys():
        in_nodes_mod_m = ERT_meta_dict["mesh_nodes_modif"]

    path = os.getcwd()
    if "path" in kwargs:
        path = kwargs["path"]

    # mesh_OUT, meshOUT = trace_mesh(mesh_CATHY,mesh_OUT,
    #                     scalar=scalar,
    #                     threshold=1e-1,
    #                     in_nodes_mod=in_nodes_mod_m
    #                     )

    data_OUT, warm_0 = trace_mesh(
        mesh_CATHY, mesh_OUT, scalar=scalar, threshold=1e-1, in_nodes_mod=in_nodes_mod_m
    )

    if len(warm_0) > 0:
        print(warm_0)

    scalar_new = scalar + "_nearIntrp2_pg_msh"
    if "time" in kwargs:
        time = kwargs["time"]
        mesh_new_attr, name_new_attr = add_attribute_2mesh(
            data_OUT, mesh_OUT, scalar_new, overwrite=True, time=time, path=path
        )
    else:
        mesh_new_attr, name_new_attr = add_attribute_2mesh(
            data_OUT, mesh_OUT, scalar_new, overwrite=True, path=path
        )

    if show:

        p = pv.Plotter(window_size=[1024 * 3, 768 * 2], notebook=False)
        p.add_mesh(mesh_new_attr, scalars=scalar_new)
        _ = p.add_bounding_box(line_width=5, color="black")
        cpos = p.show(True)

        p = pv.Plotter(window_size=[1024 * 3, 768 * 2], notebook=True)
        p.add_mesh(mesh_CATHY, scalars=scalar)
        _ = p.add_bounding_box(line_width=5, color="black")
        cpos = p.show(True)

    return mesh_new_attr, scalar_new


def CATHY_2_Resipy(mesh_CATHY, mesh_Resipy, scalar="saturation", show=False, **kwargs):

    # flip y and z axis as CATHY and Resipy have different convention for axis
    # ------------------------------------------------------------------------
    in_nodes_mod = np.array(mesh_CATHY.points)
    in_nodes_mod[:, 2] = -np.flipud(in_nodes_mod[:, 2])
    in_nodes_mod[:, 1] = -np.flipud(in_nodes_mod[:, 1])

    # check with a plot positon of the nodes for both meshes
    # ------------------------------------------------------------------------
    # p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
    # p.add_mesh(mesh_CATHY)
    # _ = p.add_points(np.array(mesh_CATHY.points), render_points_as_spheres=True,
    #                         color='red', point_size=20)
    # _ = p.add_points(in_nodes_mod, render_points_as_spheres=True,
    #                         color='blue', point_size=20)
    # _ = p.show_bounds(grid='front', all_edges=True, font_size=50)
    # cpos = p.show(True)

    path = os.getcwd()
    if "path" in kwargs:
        path = kwargs["path"]

    data_OUT = trace_mesh(
        mesh_CATHY,
        mesh_Resipy,
        scalar=scalar,
        threshold=1e-1,
        in_nodes_mod=in_nodes_mod,
    )

    # print(len(data_OUT))
    scalar_new = scalar + "_nearIntrp2Resipymsh"
    # print(mesh_Resipy)
    # get_array(mesh, name, preference='cell'

    if "time" in kwargs:
        time = kwargs["time"]
        mesh_new_attr, name_new_attr = add_attribute_2mesh(
            data_OUT, mesh_Resipy, scalar_new, overwrite=True, time=time, path=path
        )
    else:
        mesh_new_attr, name_new_attr = add_attribute_2mesh(
            data_OUT, mesh_Resipy, scalar_new, overwrite=True, path=path
        )

    if show == True:

        p = pv.Plotter(window_size=[1024 * 3, 768 * 2], notebook=True)
        p.add_mesh(mesh_new_attr, scalars=scalar_new)
        _ = p.add_bounding_box(line_width=5, color="black")
        cpos = p.show(True)

        p = pv.Plotter(window_size=[1024 * 3, 768 * 2], notebook=True)
        p.add_mesh(mesh_CATHY, scalars=scalar)
        _ = p.add_bounding_box(line_width=5, color="black")
        cpos = p.show(True)

    # if type(meshERT) is str:
    #     meshERTpv = pv.read(meshERT)

    # if savefig == True:

    #     plotter = pv.Plotter(notebook=True)
    #     _ = plotter.add_mesh(mesh_new_attr,show_edges=True)
    #     plotter.view_xz(negative=False)
    #     plotter.show_grid()
    #     plotter.save_graphic(path_CATHY + 'ERT' + str(DA_cnb) + str('.ps'),
    #                           title='ERT'+ str(DA_cnb),
    #                           raster=True,
    #                           painter=True)

    return mesh_new_attr, scalar_new


def trace_mesh(meshIN, meshOUT, scalar, threshold=1e-1, **kwargs):
    """
    Trace meshIN on meshOUT using nearest neigbour interpolation

    Parameters
    ----------
    meshIN : TYPE
        DESCRIPTION.
    meshOUT : TYPE
        DESCRIPTION.
    threshold : TYPE, optional
        DESCRIPTION. The default is 1e-1.

    Returns
    -------
    out_data : TYPE
        DESCRIPTION.
    """

    in_nodes_mod = np.array(meshIN.points)
    if "in_nodes_mod" in kwargs:
        in_nodes_mod = kwargs["in_nodes_mod"]
    meshIN.set_active_scalars(scalar)
    meshIN.points = in_nodes_mod

    # set_interpolation_radius()

    rd = max(np.diff(meshIN.points[:, 0])) / 1
    meshOUT_interp = meshOUT.interpolate(meshIN, radius=rd, pass_point_data=True)

    # plot_2d_interpolation_quality(meshIN,scalar,meshOUT,meshOUT_interp)

    result = meshOUT_interp.point_data_to_cell_data()
    out_data = result[scalar]

    warm_0 = ""
    if len(np.where(out_data == 0)) > 0:
        warm_0 = "interpolation created 0 values - replacing them by min value of input CATHY predicted ER mesh"

    out_data = np.where(out_data == 0, 1e-3, out_data)

    return out_data, warm_0


def set_interpolation_radius():

    # rd= min([abs(min(np.diff(meshIN.points[:,0]))),
    #      abs(min(np.diff(meshIN.points[:,1]))),
    #      abs(min(np.diff(meshIN.points[:,2])))
    #      ]
    #     )

    pass


def plot_2d_interpolation_quality(meshIN, scalar, meshOUT, result):

    # fig = plt.figure()
    # ax1 = plt.subplot(131)
    # # print(max(meshIN[scalar]))
    # # print(min(meshIN[scalar]))

    # # meshIN.points[:,0].min()
    # # meshOUT.points[:,0].min()
    # # meshOUT.points[:,0].max()
    # # meshOUT.points[:,1].min()
    # # meshOUT.points[:,1].max()

    # cm = plt.cm.get_cmap('RdYlBu')
    # # result = meshOUT.interpolate(meshIN, radius=rd, pass_point_data=True)
    # sc = ax1.scatter(meshIN.points[:,0],meshIN.points[:,1],c=meshIN[scalar],label='meshIN[scalar]',
    #             s=35, cmap=cm)
    # plt.colorbar(sc)
    cm = plt.cm.get_cmap("RdYlBu")
    # # fig = plt.figure()

    fig, axs = plt.subplots(2)

    # ax2 = plt.subplot(121)
    sc = axs[0].scatter(
        meshIN.points[:, 0],
        meshIN.points[:, 1],
        c=meshIN[scalar],
        label="meshIN[scalar]",
        s=35,
        cmap=cm,
    )
    # ,vmin=min(meshIN[scalar]), vmax=max(meshIN[scalar])
    axs[0].set_xlim([min(meshIN.points[:, 0]), max(meshIN.points[:, 0])])
    axs[0].set_ylim([min(meshIN.points[:, 1]), max(meshIN.points[:, 1])])

    # fig = plt.figure()
    # ax3 = plt.subplot(122)
    sc = axs[1].scatter(
        meshOUT.points[:, 0],
        meshOUT.points[:, 1],
        c=result[scalar],
        label="result[scalar]",
        s=35,
        cmap=cm,
    )

    axs[1].set_xlim([min(meshIN.points[:, 0]), max(meshIN.points[:, 0])])
    axs[1].set_ylim([min(meshIN.points[:, 1]), max(meshIN.points[:, 1])])

    import matplotlib.ticker as ticker

    def fmt(x, pos):
        a, b = "{:.2e}".format(x).split("e")
        b = int(b)
        return r"${} \times 10^{{{}}}$".format(a, b)

    plt.colorbar(sc, format=ticker.FuncFormatter(fmt))
    # plt.tight_layout()

    def uniquify(path):
        filename, extension = os.path.splitext(path)
        counter = 1

        while os.path.exists(path):
            path = filename + str(counter) + extension
            counter += 1

        return path

    savedir = os.getcwd()
    savename_test = os.path.join(savedir, "interpolation_q.png")
    savename = uniquify(savename_test)
    # print(savename)

    plt.savefig(savename)
    # plt.show()

    pass


# def find_nearest_cellcenter(node_coord,meshIN_nodes_coords,threshold=1e-1,
#                        **kwargs):
#     '''
#     Find nearest mesh node between two meshes

#     Parameters
#     ----------
#     node_coord : np.array

#     meshIN_nodes_coords : np.array

#     threshold : float
#         if distance > threshold --> closest = nan

#     Returns
#     -------
#     closest_idx : list
#         Node indice in the mesh_node.
#     closest : list
#         Node coordinate in the mesh_node.

#     '''

#     closest_idx = []
#     closest = []
#     # for i, nc in enumerate(cell_coords):
#         # euclidean distance

#     d = ( (meshIN_nodes_coords[:,0] - node_coord[0]) ** 2 +
#            (meshIN_nodes_coords[:,1]  - node_coord[1]) ** 2 +
#            (meshIN_nodes_coords[:,2]  - node_coord[2]) ** 2
#            ) ** 0.5

#     closest_idx.append(np.argmin(d))
#     closest.append(np.vstack(meshIN_nodes_coords[closest_idx,:]))

#     if d[np.argmin(d)]>threshold:
#         closest = 'nan'

#     return closest_idx, closest


def find_nearest_node(node_coord, meshIN_nodes_coords, threshold=1e-1, **kwargs):
    """
    Find nearest mesh node between two meshes

    Parameters
    ----------
    node_coord : np.array

    meshIN_nodes_coords : np.array

    threshold : float
        if distance > threshold --> closest = nan

    Returns
    -------
    closest_idx : list
        Node indice in the mesh_node.
    closest : list
        Node coordinate in the mesh_node.

    """

    closest_idx = []
    closest = []
    # for i, nc in enumerate(cell_coords):
    # euclidean distance

    d = (
        (meshIN_nodes_coords[:, 0] - node_coord[0]) ** 2
        + (meshIN_nodes_coords[:, 1] - node_coord[1]) ** 2
        + (meshIN_nodes_coords[:, 2] - node_coord[2]) ** 2
    ) ** 0.5

    closest_idx.append(np.argmin(d))
    closest.append(np.vstack(meshIN_nodes_coords[closest_idx, :]))

    if d[np.argmin(d)] > threshold:
        closest = "nan"

    return closest_idx, closest


def add_attribute_2mesh(
    data, mesh, name="ER_pred", overwrite=True, saveMesh=True, **kwargs
):
    """
    add a new mesh attribute to a vtk file

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    mesh : TYPE
        DESCRIPTION.
    name : TYPE
        DESCRIPTION.
    overwrite : TYPE, optional
        DESCRIPTION. The default is True.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    # for k in kwargs:
    #     print(k)
    # print(mesh)

    if type(mesh) is str:
        mesh = pv.read(mesh)

    try:
        mesh.point_data[name] = data
    except:
        mesh.cell_data[name] = data

    meshname = name + ".vtk"

    if saveMesh:
        path = os.getcwd()
        if "path" in kwargs:
            path = kwargs["path"]

        if "time" in kwargs:
            time = kwargs["time"]
            meshname = name + str(time) + ".vtk"
            mesh.save(path + meshname, binary=False)
        else:
            mesh.save(path + meshname, binary=False)

        # if overwrite==True:
        #     mesh.save(full_path)

    return mesh, name


def trace_mesh_pg(meshIN, meshOUT, method="spline", **kwargs):
    """
    Interpolate CATHY mesh (structured) into pygimli mesh (structured) using pygimli meshtools
    # Specify interpolation method 'linear, 'spline', 'harmonic'

    """
    meshIN = pg.load(meshIN)
    meshOUT = pg.load(meshOUT)
    out_data = pg.interpolate(meshIN["ER_converted_CATHYmsh*"], meshOUT, method=method)

    return out_data


#%%


def map_layers_2_DEM(layers, DEM, zone, dem_parameters):

    zone3d_layers_top, zone3d_layer_bot = get_zone3d_layer_depths(zone, dem_parameters)

    dem_mat3d_layers_top = [DEM - zz for zz in zone3d_layers_top]
    dem_mat3d_layers_bot = [DEM - zz for zz in zone3d_layer_bot]

    zone3d_topflag = []
    for li in range(dem_parameters["nstr"]):
        zone3d_topflag_li = np.ones(np.shape(dem_mat3d_layers_top[0]))

        bool_top_lli = []
        for ll in layers.keys():
            layers_adj_top = DEM - layers[ll][0]
            layers_adj_bot = DEM - layers[ll][1]
            # differences between top of the layer i of the mesh and top of the desired layer
            diff_top = dem_mat3d_layers_top[li] - layers_adj_top
            cond1 = diff_top <= 0
            diff_bot = dem_mat3d_layers_bot[li] - layers_adj_bot
            cond2 = diff_bot >= 0
            # if li<dem_parameters["nstr"]-1:
            #     cond2 = abs(diff)<= abs(dem_mat3d_layers_top[li+1] - dem_mat3d_layers_top[li])
            # else:
            #     cond2 = np.ones(np.shape(diff),dtype=bool)
            bool_top_lli = cond1 & cond2
            zone3d_topflag_li[bool_top_lli] = ll
        zone3d_topflag.append(zone3d_topflag_li)
    zone3d_topflag = np.array(zone3d_topflag)

    return zone3d_topflag


def get_zone3d_layer_depths(zone_raster, dem_parameters):

    zone3d_top = []
    zone3d_bot = []
    zone3d_raster = np.array([zone_raster] * dem_parameters["nstr"])

    for li in range(dem_parameters["nstr"]):
        top, bot = get_layer_depth(dem_parameters, li)
        zone3d_top.append(zone3d_raster[li] * top)
        zone3d_bot.append(zone3d_raster[li] * bot)

    return np.array(zone3d_top), np.array(zone3d_bot)


def get_layer_depths(dem_parameters):

    layers_top = []
    layers_bottom = []
    for li in range(dem_parameters["nstr"]):
        layeri_top, layeri_bottom = get_layer_depth(dem_parameters, li)
        layers_top.append(layeri_top)
        layers_bottom.append(layeri_bottom)
    return layers_top, layers_bottom


def get_layer_depth(dem_parameters, li):

    dempar = dem_parameters["zratio(i),i=1,nstr"].split("\t")
    dempar_ratio = [float(d) for d in dempar]

    # layeri_top = np.round(dempar_ratio[li]*simu.dem_parameters["base"]/simu.dem_parameters['nstr'],2)
    # layeri_top = dempar_ratio[li]*simu.dem_parameters["base"]/simu.dem_parameters['nstr']
    if li == 0:
        layeri_top = 0
    else:
        # layeri_top = np.cumsum(dempar_ratio[0:li])[-1]*(simu.dem_parameters["base"]/(simu.dem_parameters['nstr']-1))
        layeri_top = np.cumsum(dempar_ratio[0:li])[-1] * (dem_parameters["base"])

        # layeri_top = np.cumsum(dempar_ratio[0:li])[-1]

    if (li + 1) < len(dempar_ratio):
        # layeri_bottom = np.round(dempar_ratio[li+1]*simu.dem_parameters["base"]/simu.dem_parameters['nstr'],2)
        # layeri_bottom = layeri_top + dempar_ratio[li+1]*simu.dem_parameters["base"]/(simu.dem_parameters['nstr']-1)
        # layeri_bottom = np.cumsum(dempar_ratio[0:li+1])[-1]*(simu.dem_parameters["base"]/(simu.dem_parameters['nstr']-1))
        layeri_bottom = np.cumsum(dempar_ratio[0 : li + 1])[-1] * (
            dem_parameters["base"]
        )
    else:
        layeri_bottom = dem_parameters["base"]

    # print(layeri_top,layeri_bottom)
    return layeri_top, layeri_bottom


def zone3d(zone, dem_parameters):

    # define zone in 3dimension - duplicate number of layer (nstr) times
    # ---------------------------------------------------------------
    zones3d = [zone] * dem_parameters["nstr"]

    fig, axs = plt.subplots(
        int(dem_parameters["nstr"] / 2) + 1,
        2,
        sharex=True,
        sharey=(True),
        constrained_layout=False,
    )
    plt.tight_layout()
    axs = axs.ravel()
    for li in range(dem_parameters["nstr"]):
        layeri_top, layeri_bottom = get_layer_depth(dem_parameters, li)
        layer_str = (
            "[" + str(f"{layeri_top:.2E}") + "-" + str(f"{layeri_bottom:.2E}") + "]"
        )

        pmesh = cplt.show_raster(
            zones3d[li], prop=layer_str, cmap="jet", ax=axs[li], vmin=0, vmax=1
        )

    plt.colorbar(pmesh, ax=axs[:], location="right", shrink=0.6, cmap="jet")
    plt.tight_layout()

    return zones3d


def create_layers_inzones3d(simu, zones3d, layers_names, layers_depths=[[0, 1e99]]):

    # Loop over layers and zones and change flag if depth conditions is not respected
    # ----------------------------------------------------------------------------
    zones3d_layered = np.ones(np.shape(zones3d))
    for li in range(simu.dem_parameters["nstr"]):

        layeri_top, layeri_bottom = get_layer_depth(simu, li)
        # layer_str = '[' + str(layeri_top) + '-' + str(layeri_bottom) + ']'

        for zi in range(len(layers_names)):
            print("layers %i analysis" % zi)
            # print(layers_depths[zi][0])
            # print(layeri_top)
            # print(layers_depths[zi][1])
            # print(layeri_bottom)
            zi = 0

            # print(layers_depths[zi][0]<=layeri_top)
            # if depth of zone i is sup to mesh layers depth
            # ---------------------------------------------------------------------
            # if (layers_depths[zi][0]>=layeri_top) & (layers_depths[zi][1]<layeri_bottom):
            if (layers_depths[zi][0] <= layeri_top) & (
                layers_depths[zi][1] >= layeri_bottom
            ):
                # if (depths_ordered[zi][1]>=layeri_bottom):
                print("conds ok --> zi:" + str(zi))
                print(layers_depths[zi])
                print(layeri_top, layeri_bottom)

                zones3d_layered[li][zones3d[li] == zi + 1] = zi + 2
                print(
                    "replacing "
                    + str(np.sum(zones3d[li] == zi + 1))
                    + " values by"
                    + str(zi + 1)
                )

                if 10.5 * layers_depths[zi][1] < layeri_bottom:
                    raise ValueError(
                        "Required layer is finer than mesh layers - refine mesh"
                    )
            else:
                print("conds not ok")
                print(layers_depths[zi])
                print(layeri_top, layeri_bottom)

    return zones3d_layered


def plot_zones3d_layered(simu, zones3d_layered):

    fig, axs = plt.subplots(
        int(simu.dem_parameters["nstr"] / 2) + 1,
        2,
        sharex=False,
        sharey=(False),
        constrained_layout=False,
    )
    axs = axs.flat
    for li in range(simu.dem_parameters["nstr"]):

        layeri_top, layeri_bottom = get_layer_depth(simu, li)
        layer_str = "[" + str(layeri_top) + "-" + str(layeri_bottom) + "]"

        pmesh = cplt.show_raster(
            zones3d_layered[li],
            prop=layer_str,  # , cmap='jet',
            ax=axs[li],
            vmin=0,
            vmax=2,
        )

    plt.colorbar(pmesh, ax=axs[:], location="right", shrink=0.6, cmap="jet")
    plt.tight_layout()


def het_soil_layers_mapping_generic(
    simu, propertie_names, SPP, layers_names, layers_depths
):

    # extend to 3d the zone raster file
    # -----------------------------------
    zones3d = zone3d(simu)

    # insert layers flag into the 3d the zone raster file
    # -----------------------------------
    zones3d_layered = create_layers_inzones3d(
        simu, zones3d, layers_names, layers_depths
    )

    # plot 3d zones files layered
    # ------------------------------------------
    plot_zones3d_layered(simu, zones3d_layered)

    # np.shape(zones3d_layered)
    # np.shape(zones3d_axis_swap)

    # Loop over soil properties names
    # -------------------------------------------------
    index_raster = np.arange(0, simu.hapin["N"] * simu.hapin["M"])
    zones3d_axis_swap = np.swapaxes(zones3d_layered, 0, 2)
    prop_df = []  # properties dataframe
    layers_id = ["L" + str(i + 1) for i in range(simu.dem_parameters["nstr"])]
    layers_id = np.flipud(layers_id)

    for i, p in enumerate(propertie_names):

        p_df = np.zeros(np.shape(zones3d_axis_swap))

        # Loop over soil layers and assign value of soil properties
        # -------------------------------------------------
        for k, lname in enumerate(layers_names):
            p_df[zones3d_axis_swap == k + 1] = SPP[k][i]

        df_tmp = pd.DataFrame(
            np.vstack(p_df),
            columns=layers_id,
            index=index_raster,
        )

        prop_df.append(df_tmp)

    SoilPhysProp_df_het_layers_p = pd.concat(prop_df, axis=1, keys=propertie_names)
    SoilPhysProp_df_het_layers_p.index.name = "id raster"
    SoilPhysProp_df_het_layers_p.columns.names = ["soilp", "layerid"]

    SPP_map_dict = {}
    for p in propertie_names:
        SPP_map_dict[p] = []
        for li in range(simu.dem_parameters["nstr"]):
            v = SoilPhysProp_df_het_layers_p.iloc[0][p].loc["L" + str(li + 1)]
            SPP_map_dict[p].append(v)

    return SoilPhysProp_df_het_layers_p, SPP_map_dict


def _subplot_cellsMarkerpts(mesh_pv_attributes, xyz_layers0, xyzlayers1):

    pl = pv.Plotter(shape=(1, 2))
    # pl.add_mesh(mesh_pv_attributes,
    #             show_edges=True,
    #             )
    pl.show_grid()

    actor = pl.add_points(
        xyz_layers0[:, :-1],
        point_size=10,
        scalars=xyz_layers0[:, -1],
    )
    pl.subplot(0, 1)

    pl.show_grid()

    actor = pl.add_points(
        xyzlayers1[:, :-1],
        point_size=10,
        scalars=xyzlayers1[:, -1],
    )

    pl.show()


def _plot_cellsMarkerpts(mesh_pv_attributes, xyz_layers, workdir, project_name):

    # pl = pv.Plotter(off_screen=(True))
    pl = pv.Plotter()
    pl.add_mesh(mesh_pv_attributes, show_edges=True, opacity=0.4)
    pl.show_grid()

    actor = pl.add_points(
        xyz_layers[:, :-1],
        point_size=10,
        scalars=xyz_layers[:, -1],
    )
    pl.show(screenshot=os.path.join(workdir, project_name, "layersMarkers.png"))


def _build_xyz_marker_mat_squareDEM(
    zones3d_layered,
    mesh_pv_attributes,
    hapin,
    dem_parameters,
    workdir,
    project_name,
    to_nodes,
):

    # top, bot = get_layer_depth(dem_parameters,li)

    # get layers properties
    # ----------------------------------------------------------------
    dempar = dem_parameters["zratio(i),i=1,nstr"].split("\t")
    dempar_ratio = [float(d) for d in dempar]

    layeri_top = []
    layeri_top = [
        np.cumsum(dempar_ratio[0 : li + 1])[-1] * (dem_parameters["base"])
        for li in range(dem_parameters["nstr"])
    ]
    layeri_top.insert(0, 0)
    layeri_center = [
        (layeri_top[lti + 1] + layeri_top[lti]) / 2
        for lti in range(len(layeri_top) - 1)
    ]

    # get dem coordinates
    # ----------------------------------------------------------------
    y, x, dem = cplt.get_dem_coords(
        workdir=workdir, project_name=project_name, hapin=hapin
    )

    # build dem coordinates of the size of the DEM
    # ----------------------------------------------------------------
    xn = [np.ones(len(y)) * xu for xu in np.unique(x)]
    xn = np.hstack(xn)
    yn = np.tile(y, len(x))

    xyz = np.c_[
        xn - dem_parameters["delta_x"] / 2,
        yn - dem_parameters["delta_y"] / 2,
        np.hstack(dem.T),
    ]

    xyz_celldown = np.c_[
        xn - dem_parameters["delta_x"] / 4,
        yn - dem_parameters["delta_y"] / 4,
        np.hstack(dem.T),
    ]

    xyz_cellup = np.c_[
        xn + dem_parameters["delta_x"] / 4,
        yn + dem_parameters["delta_y"] / 4,
        np.hstack(dem.T),
    ]

    xyz_cells = np.r_[xyz_cellup, xyz_celldown]
    # len(xyz_cellup)

    xyz_nodes = np.c_[xn, yn, np.hstack(dem.T)]

    meshpts = mesh_pv_attributes.points
    cellpts = mesh_pv_attributes.cell_centers().points

    xyz_layers_nodes = []
    for i, li in enumerate(layeri_top[:-1]):
        marker_zone = np.ravel(zones3d_layered[i])
        xyz_layers_nodes.append(
            np.c_[xyz_nodes[:, 0:2], xyz_nodes[:, 2] - li, marker_zone]
        )

    xyz_layers = []
    for i, li in enumerate(layeri_center[:]):
        marker_zone = np.ravel(zones3d_layered[i])
        xyz_layers.append(np.c_[xyz[:, 0:2], xyz[:, 2] - li, marker_zone])

    xyz_layers_cells_top = []
    for i, li in enumerate(layeri_top[:-1]):
        marker_zone = np.hstack([np.ravel(zones3d_layered[i])] * 2)
        xyz_layers_cells_top.append(
            np.c_[xyz_cells[:, 0:2], xyz_cells[:, 2] - li, marker_zone]
        )

    xyz_layers_cells_cent = []
    for i, li in enumerate(layeri_center[:]):
        marker_zone = np.hstack([np.ravel(zones3d_layered[i])] * 2)
        xyz_layers_cells_cent.append(
            np.c_[xyz_cells[:, 0:2], xyz_cells[:, 2] - li, marker_zone]
        )

    xyz_layers = np.vstack(xyz_layers)
    xyz_layers_nodes = np.vstack(xyz_layers_nodes)
    xyz_layers_cells = np.vstack(np.r_[xyz_layers_cells_top, xyz_layers_cells_cent])

    _plot_cellsMarkerpts(mesh_pv_attributes, xyz_layers, workdir, project_name)

    _find_nearest_point2DEM(
        to_nodes,
        mesh_pv_attributes,
        workdir,
        project_name,
        xyz_layers_cells,
        xyz_layers_nodes,
    )


def _find_nearest_point2DEM(
    to_nodes,
    mesh_pv_attributes,
    workdir,
    project_name,
    xyz_layers_cells=[],
    xyz_layers_nodes=[],
):

    if to_nodes:
        # loop over mesh cell centers and find nearest point to dem
        # ----------------------------------------------------------------
        node_markers = []
        for nmesh in mesh_pv_attributes.points:
            # euclidean distance
            d = (
                (xyz_layers_nodes[:, 0] - nmesh[0]) ** 2
                + (xyz_layers_nodes[:, 1] - nmesh[1]) ** 2
                + (abs(xyz_layers_nodes[:, 2]) - abs(nmesh[2])) ** 2
            ) ** 0.5
            node_markers.append(xyz_layers_nodes[np.argmin(d), 3])
            # dbackup_nodes.append(min(d))
        # add data to the mesh
        # ----------------------------------------------------------------
        mesh_pv_attributes["node_markers"] = node_markers

        # simu.mesh_pv_attributes.save('mesh_with_markers.vtk',
        #                              binary=False)

    else:
        # loop over mesh cell centers and find nearest point to dem
        # ----------------------------------------------------------------
        cell_markers = []
        for nmesh in mesh_pv_attributes.cell_centers().points:
            # euclidean distance
            d = (
                (xyz_layers_cells[:, 0] - nmesh[0]) ** 2
                + (xyz_layers_cells[:, 1] - nmesh[1]) ** 2
                + (abs(xyz_layers_cells[:, 2]) - abs(nmesh[2])) ** 2
            ) ** 0.5
            cell_markers.append(xyz_layers_cells[np.argmin(d), 3])
            # dbackup_cell.append(min(d))

        # add data to the mesh
        # ----------------------------------------------------------------
        mesh_pv_attributes["cell_markers"] = cell_markers

    mesh_pv_attributes.save(
        os.path.join(
            workdir,
            project_name,
            "vtk/",
            project_name + ".vtk",
        ),
        binary=False,
    )

    # plt.plot(dbackup_cell)
    # plt.plot(dbackup_nodes)


def add_markers2mesh(
    zones3d_layered,
    mesh_pv_attributes,
    dem_parameters,
    workdir,
    project_name,
    hapin,
    to_nodes=False,
):

    x, y, dem = cplt.get_dem_coords(
        workdir=workdir,
        project_name=project_name,
        hapin=hapin,
        tranposeDEM=True,
    )

    if len(x) == len(y):

        _build_xyz_marker_mat_squareDEM(
            zones3d_layered,
            mesh_pv_attributes,
            hapin,
            dem_parameters,
            workdir,
            project_name,
            to_nodes,
        )

    else:

        # Create a regular mesh from the DEM x and y coordinates and elevation
        # ------------------------------------------------------------------
        xgrid, ygrid = np.meshgrid(x, y)
        XY = np.meshgrid(x, y)
        dem_flip = np.flipud(np.fliplr(dem.T))
        grid_coords_dem = np.array(
            [
                np.ravel(xgrid),
                np.ravel(ygrid),
                np.ravel(dem_flip),
                np.ones(len(np.ravel(dem_flip))),
            ]
        ).T

        # Get layer top and bottom
        # ------------------------------------------------------------------
        zone3d_top = []
        zone3d_bot = []
        for li in range(dem_parameters["nstr"]):
            top, bot = get_layer_depth(dem_parameters, li)
            zone3d_top.append(zones3d_layered[li] * top)
            zone3d_bot.append(zones3d_layered[li] * bot)

        zone3d_top = np.array(zone3d_top)
        zone3d_bot = np.array(zone3d_bot)

        dem_mat3d_layers = [
            dem_flip - zz for zz in zone3d_top
        ]  # -(zone3d_top-zone3d_bot)/2
        # dem_mat3d_layers = [dem_flip - zz for zz in (zone3d_top-zone3d_bot)/2] #-(zone3d_top-zone3d_bot)/2

        # Reduce all to 1D
        # ------------------------------------------------------------------
        dem_mat_stk = np.ravel(dem_mat3d_layers)

        grid_coords_stk_rep = np.vstack(
            np.array([grid_coords_dem] * dem_parameters["nstr"])
        )

        zones3d_col_stk = np.ravel(zones3d_layered)
        # xyz_layers  = np.c_[grid_coords_stk_rep[:,:-1],zones3d_col_stk]
        xyz_layers = np.c_[grid_coords_stk_rep[:, :-2], dem_mat_stk, zones3d_col_stk]

        # Plot to check position of points VS mesh
        # ------------------------------------------------------------------
        _plot_cellsMarkerpts(mesh_pv_attributes, xyz_layers, workdir, project_name)

        # Assign marker to mesh and overwrite it
        # ------------------------------------------------------------------
        _find_nearest_point2DEM(
            to_nodes,
            mesh_pv_attributes,
            workdir,
            project_name,
            xyz_layers,
            xyz_layers,
        )
