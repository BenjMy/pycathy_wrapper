"""
Input plots
===========

This example shows how to use pyCATHY object to plot inputs of the hydrological model.

*Estimated time to run the notebook = 5min*

"""

# map_prop_veg ?
# map_prop2zone



#%% Import packages

from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt
import numpy as np

#%% add CATHY object

path2prj ='weil_exemple_inputs_plot' # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj)
# simu.run_preprocessor()

#%%
simu.show_input(prop='dem')

#%%

# show time atmbc
# simu.show_input(prop='atmbc')

# In progress --> show spatial atmbc

# simu.update_dem_parameters()
# simu.update_prepo_inputs()
#%%
# In progress --> This will automatically create a new vtk mesh containing the 
simu.show_input(prop='root_map')


#%% update root map: add properties for zone 2
# Add a new zone

simu.update_prepo_inputs()
simu.update_veg_map() # calling without args to get the default values
simu.update_soil()

#%%
veg_map = simu.veg_map
veg_map[2:6,5:14] = 2
simu.update_veg_map(veg_map)
simu.show_input(prop='root_map')

#%% Update Feddes parameters
# Feddes is a dictionnary with 6 entries, and for each a list

FP_map_1zone = simu.soil_FP['FP_map'] # read existing mapping
FP_map_2zones = {}
for k in FP_map_1zone:
    if k == 'ZROOT':
        ZROOT_zone2 = FP_map_1zone['ZROOT'][0]/2
        FP_map_2zones[k] = [FP_map_1zone[k][0],ZROOT_zone2]
    else:
        FP_map_2zones[k] = [FP_map_1zone[k][0],FP_map_1zone[k][0]]
    
# simu.show_input(prop='soil', yprop='ZROOT', layer_nb=12)

#%%
simu.update_soil(FP_map=FP_map_2zones,
                 show=True)

simu.update_zone()
simu.show_input(prop='soil', yprop='PERMX', layer_nb=4)

# Here we can imaging to get a more complexe vegetation map from remote sensing data instead

#%%

simu.update_prepo_inputs()

#%%
# This will automatically create a new vtk mesh containing the zone flags
# error --> number of tretra in grid3d < n of tretra in the mesh (mission one element)
simu.update_zone()

#%%

simu.show_input(prop='soil', yprop='PERMX', layer_nb=1)
simu.show_input(prop='soil', yprop='POROS', layer_nb=2)

#%%
# Show layer number 10

simu.show_input(prop='soil', yprop='VGNCELL', layer_nb=10)

#%%
simu.update_soil()
df_soil, _ = simu.read_inputs('soil')
df = simu.read_inputs('soil')

#%% Add a new zone
zones = simu.zone
simu.update_prepo_inputs()
zones[5:14,5:14] = 2
simu.update_zone(zones)
simu.show_input(prop='zone')

#%% update soil: add properties for zone 2
# we just need to build a dictionnary as: {property: [value_zone1, value_zone2]}

# what if dimension of the heteregeneity is 3d?

SPP_map_1zone = simu.soil_SPP['SPP_map'] # read existing mapping
SPP_map_2zones = {}
for k in SPP_map_1zone:
    if k == 'PERMX':
        PERMX_zone2 = SPP_map_1zone['PERMX'][0]/2
        SPP_map_2zones[k] = [SPP_map_1zone[k][0],PERMX_zone2]
    else:
        SPP_map_2zones[k] = [SPP_map_1zone[k][0],SPP_map_1zone[k][0]]
    

simu.update_soil(SPP_map=SPP_map_2zones)

#%%
simu.show_input(prop='soil', yprop='PERMX', layer_nb=2)

simu.show_input(prop='soil', yprop='PERMX', layer_nb=12)

#%%
# Create grid
# simu.run_processor(IPRT1=3) # only generate the grid

#%% 
# In progress
# simu.show_bc(layer_nb=1, time=0)
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots(1,3)


def add_markers2mesh(simu,zones3d_layered,to_nodes=False):
    
    # simu.mesh_pv_attributes.clear_point_data()
    # simu.mesh_pv_attributes.clear_cell_data()
    # simu.mesh_pv_attributes.point_data.keys()
    # simu.mesh_pv_attributes.cell_data.keys()
    # simu.mesh_pv_attributes
    
    
    # get layers properties
    # ----------------------------------------------------------------
    dempar = simu.dem_parameters['zratio(i),i=1,nstr'].split('\t')
    dempar_ratio = [float(d) for d in dempar]
    
    layeri_top = []
    layeri_top = [np.cumsum(dempar_ratio[0:li+1])[-1]*(simu.dem_parameters["base"]) 
                  for li in range(simu.dem_parameters["nstr"])]
    layeri_top.insert(0, 0)
    
    layeri_center = [(layeri_top[lti+1]+layeri_top[lti])/2 for lti in range(len(layeri_top)-1)]
    
    # get dem coordinates 
    # ----------------------------------------------------------------
    y, x, dem = cplt.get_dem_coords(workdir=simu.workdir,
                        project_name=simu.project_name,
                        hapin=simu.hapin
                        )
    
    
    # build dem coordinates of the size of the DEM
    # ----------------------------------------------------------------
    xn = [np.ones(len(y))*xu for xu in np.unique(x)]
    xn = np.hstack(xn)
    yn = np.tile(np.flipud(y),len(x))
    
    xyz = np.c_[xn+simu.dem_parameters['delta_x']/2,
                yn+simu.dem_parameters['delta_y']/2,
                np.hstack(simu.DEM.T)
                ]
    xyz_layers = []
    for i, li in enumerate(layeri_center):
        marker_zone = np.ravel(zones3d_layered[i])
        xyz_layers.append(np.c_[xyz[:,0:2], xyz[:,2] -li, marker_zone])
    
    xyz_layers = np.vstack(xyz_layers)
    
    if to_nodes:
    # loop over mesh cell centers and find nearest point to dem
    # ----------------------------------------------------------------
        node_markers = []
        for nmesh in simu.mesh_pv_attributes.points:
            # euclidean distance
            d=((xyz_layers[:, 0] - nmesh[0]) ** 2 +
                (xyz_layers[:, 1] - nmesh[1]) ** 2 +
                (abs(xyz_layers[:, 2]) - abs(nmesh[2])) ** 2
                ) ** 0.5
            node_markers.append(xyz_layers[np.argmin(d),3])
    
        # add data to the mesh
        # ----------------------------------------------------------------
        simu.mesh_pv_attributes['node_markers'] = node_markers
        
        # simu.mesh_pv_attributes.save('mesh_with_markers.vtk',
        #                              binary=False)
    
    else:
        # loop over mesh cell centers and find nearest point to dem
        # ----------------------------------------------------------------
        cell_markers = []
        for nmesh in simu.mesh_pv_attributes.cell_centers().points:
            # euclidean distance
            d=((xyz_layers[:, 0] - nmesh[0]) ** 2 +
                (xyz_layers[:, 1] - nmesh[1]) ** 2 +
                (abs(xyz_layers[:, 2]) - abs(nmesh[2])) ** 2
                ) ** 0.5
            cell_markers.append(xyz_layers[np.argmin(d),3])
    
        # add data to the mesh
        # ----------------------------------------------------------------
        simu.mesh_pv_attributes['cell_markers'] = cell_markers
        
        simu.mesh_pv_attributes.save('mesh_with_markers.vtk',
                                     binary=False)
        

# #%% 
# simu.create_mesh_vtk()
# # simu.grid3d = in_CT.read_grid3d(simu.project_name)

            
# #%%
# # extend to 3d the zone raster file
# # -----------------------------------
# zones3d = zone3d(simu)

# # insert layers flag into the 3d the zone raster file
# # -----------------------------------
# zones3d_layered = create_layers_inzones3d(simu,
#                                         zones3d,
#                                         zone_names,
#                                         layers_depths
#                                         )
    
# #%%
# # add zones markers to mesh
# # ------------------------------------------

# add_markers2mesh(simu,zones3d_layered, to_nodes=False)
# add_markers2mesh(simu,zones3d_layered, to_nodes=True)
              
# ic_cells, id_cells2nodes = simu.map_dem_prop_2mesh('ic',scenario['ic'], to_nodes=False)
# ic_nodes = simu.map_dem_prop_2mesh('ic',scenario['ic'], to_nodes=True)






