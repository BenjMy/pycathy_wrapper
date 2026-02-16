"""
ERA5 ETp spatially variable as atmbc read_inputs
=================================================

"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pyCATHY.meshtools as mt
from pyCATHY import cathy_tools
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.plotters import cathy_plots as cplt
import pooch

import pyCATHY.petro as petro
# import utils

import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)

# import utils_Bousval
import pyvista as pv
from pathlib import Path
import xarray as xr
import geopandas as gpd

rootPath = Path(os.getcwd())

#%% Init CATHY model
# ------------------------
# path2prj = rootPath / "../data/solution_ET/"  # add your local path here
# simu = cathy_tools.CATHY(dirName=path2prj,
#                          prj_name="ERA5_ETp_spatially_from_weill"
#                          )


path2prj = "../SSHydro/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj,
			prj_name="ERA5_ETp_spatially_from_weill"
			)

# figpath = path2prj
# dimveg = 2

#%%



rootPathData= Path('/home/z0272571a@campus.csic.es/Nextcloud/BenCSIC/Codes/Tech4agro_org/test_Majadas_centum_dataset/')

print('loading datasets')

# ETa_Majadas_Buffer20k
# ETa_Majadas_H2Bassin

Majadas_ETa_dataset = rootPathData/'ETa_Majadas_H2Bassin.netcdf'
Majadas_ETp_dataset = rootPathData/'ETp_Majadas_H2Bassin.netcdf'
Majadas_CLC_dataset = rootPathData/'CLCover_Majadas.netcdf'
Majadas_rain_dataset = rootPathData/'RAIN_Majadas_H2Bassin.netcdf'

Majadas_CLC_gdf = gpd.read_file(rootPathData/'BassinH2_Majadas_corrected.shp')

ETa_ds = xr.load_dataset(Majadas_ETa_dataset,engine='netcdf4')
ETa_ds = ETa_ds.rename({"__xarray_dataarray_variable__": "ETa"})  # Rename the main variable to 'ETa'

ETp_ds = xr.load_dataset(Majadas_ETp_dataset)
ETp_ds = ETp_ds.rename({"__xarray_dataarray_variable__": "ETp"})  # Rename the main variable to 'ETa'

Rain_ds = xr.load_dataset(Majadas_rain_dataset)
Rain_ds = Rain_ds.rename({"__xarray_dataarray_variable__": "Rain"})  # Rename the main variable to 'ETa'

Rain_ds = Rain_ds.sel(time=slice(ETp_ds.time[0], None))

ETa_aligned, ETp_aligned, Rain_aligned = xr.align(ETa_ds,
                                                  ETp_ds,
                                                  Rain_ds,
                                                  join='outer',
                                                  fill_value=0
                                                  )




# Create the net variable directly
grid_xr_atmbc = xr.Dataset({
    'net_atmbc': Rain_aligned['Rain'] - ETp_aligned['ETp'],
    'Rain': Rain_aligned['Rain'],
    'ETp': ETp_aligned['ETp'],
    'ETa': ETa_aligned['ETa']
})
grid_xr_atmbc = grid_xr_atmbc.isel(band=0)

# grid_xr_atmbc = xr.merge([ETa_ds, ETp_ds, Rain_ds])
# grid_xr_atmbc['net_atmbc'] = grid_xr_atmbc['Rain'] - grid_xr_atmbc['ETp']
CLC = xr.load_dataset(Majadas_CLC_dataset)  # Load the CLC dataset


t0 = ETp_ds.time[0]
elapsed_sec = (ETp_ds.time - t0).astype('timedelta64[s]').values
ETp_ds.isel(band=0)['ETp'].to_numpy()

ETp_ds.isel(band=0,time=2)['ETp'].plot.imshow()
# ETp_ds.isel(band=0,time=20)['ETp'].plot.imshow()

grid_xr_atmbc.isel(time=2)['ETp'].plot.imshow()
grid_xr_atmbc.isel(time=2)['net_atmbc'].plot.imshow()
grid_xr_atmbc.isel(time=2)['Rain'].plot.imshow()



#%% spatially variable atmospheric boundary condition inputs

grid_xr_atmbc

# DEM, dem_header = simu.read_inputs('dem')
# DEM_new = np.ones(np.shape(DEM)) #np.mean(DEM)*
dim_y, dim_x = len(grid_xr_atmbc.y), len(grid_xr_atmbc.x)

mask = np.isnan(grid_xr_atmbc['ETa'].isel(time=0))
# # Create DEM array
DEM_new = np.ones((dim_y, dim_x), dtype=np.float32)  # initialize with 1 (or mean(DEM))
DEM_new[mask] = -9999  # replace NaNs with -9999

# Optional: slightly modify the last node to avoid preprocessing errors
# Apply 1e-3 to the first non-masked value
non_masked_indices = np.where(~mask)
if len(non_masked_indices[0]) > 0:
    first_y, first_x = non_masked_indices[0][0], non_masked_indices[1][0]
    DEM_new[first_y, first_x] -= 1e-3



DEM_new[-1,-1] = DEM_new[-1,-1] - 1e-3
simu.update_prepo_inputs(DEM_new,
                         xllcorner=int(grid_xr_atmbc.x.min().values),
                         yllcorner=int(grid_xr_atmbc.y.min().values),
                         delta_x=1,
                         delta_y=1,
                         ivert=1,
                         )

simu.update_dem_parameters()


fig, ax = plt.subplots(1)
img = ax.imshow(DEM_new,vmin=0,vmax=500)
plt.colorbar(img)
simu.show_input(prop="dem",vmin=0,vmax=500)

#%%

simu.update_veg_map()
#%%
simu.run_preprocessor(verbose=True)
simu.run_processor(IPRT1=3, verbose=True)

simu.create_mesh_vtk(verbose=True)

#%%
# os.getcwd()
# pathtoscenrario = Path('../../../EOMAJI/prepro/')
# grid_xr = xr.open_dataset(pathtoscenrario/'grid_xr_baseline_AquaCrop_sc0_weather_reference.netcdf')


#%%
grid_xr_sliced = grid_xr_atmbc.isel(time=slice(0, 10))
grid_xr_sliced.isel(time=2)['ETp'].plot.imshow()

t_atmbc = list(grid_xr_sliced['time'].values.astype('timedelta64[s]').astype(int))
# net_atmbc1d = grid_xr_sliced.net_atmbc.mean(dim=['x','y'])

# net_atmbc1d = grid_xr_sliced.net_atmbc
# Slice the first 400 time steps and modify 'net_atmbc'
# Get the midpoint of the x dimension
# mid_x = grid_xr_sliced.sizes['x'] // 2

# Create a mask to select half of the x range
# mask = grid_xr_sliced['x'] < grid_xr_sliced['x'].isel(x=mid_x)

# Divide 'net_atmbc' values by 2 for the selected range
# grid_xr_sliced['net_atmbc'] = grid_xr_sliced['net_atmbc'].where(~mask, grid_xr_sliced['net_atmbc'] / 2)

# Verify or plot
# grid_xr_sliced = grid_xr.isel(time=slice(0, 10))
# grid_xr_sliced.net_atmbc.plot.imshow(col="time", col_wrap=4)


#%% PAD RAIN and ETp
# ETp_meshnodes = Majadas_utils.xarraytoDEM_pad(grid_xr_sliced['net_atmbc'])

# grid_xr_sliced['net_atmbc'].isel(time=1).plot.imshow()
grid_xr_atmbc['ETp'].isel(time=29).plot.imshow()
grid_xr_atmbc['ETa'].isel(time=67).plot.imshow()

net_atmbc_meshnodes = mt.xarraytoDEM_pad(grid_xr_sliced['net_atmbc'])
t_atmbc = list(grid_xr_sliced['time'].values.astype('timedelta64[s]').astype(int))
mask = ~np.isnan(net_atmbc_meshnodes).any(axis=1)  # keep nodes with no NaNs
# net_atmbc_meshnodes_noNan = net_atmbc_meshnodes[:, mask]


# Convert xarray DataArray to NumPy array
net_atmbc_np = net_atmbc_meshnodes.values  # shape: (time, nodes)

# Create mask for nodes (columns) without any NaNs across all times
mask = ~np.isnan(net_atmbc_np).any(axis=0)  # axis=0 = time dimension

# Apply mask to remove NaN columns
net_atmbc_meshnodes_noNan = net_atmbc_np[:, mask]

# Time array
t_atmbc = grid_xr_sliced['time'].values.astype('timedelta64[s]').astype(int)

print(f"Original shape: {net_atmbc_np.shape}")
print(f"Shape after removing NaN nodes: {net_atmbc_meshnodes_noNan.shape}")



plt.imshow(net_atmbc_meshnodes_3d[0])
# np.unique(grid_xr_sliced.ETp_daily.isel(time=0))
# Count NaNs per column
import numpy as np

# Count NaNs per row
nan_counts_per_row = np.isnan(net_atmbc_3d_reshaped).sum(axis=1)

# Optional: print stats
for i, n_nan in enumerate(nan_counts_per_row):
    print(f"Time step {i}: {n_nan} NaNs out of {net_atmbc_3d_reshaped.shape[1]} columns")


simu.update_atmbc(
                  HSPATM=0, # spatially variable atmbc
                  IETO=1,
                  time=t_atmbc,
                  netValue=ETp_3d_reshaped,
                )
simu.grid3d['nnod']
np.shape(ETp_3d_reshaped)


#%%
fig, ax= plt.subplots()
simu.show_input('atmbc',ax=ax)
fig.savefig(os.path.join(simu.workdir + "/ERA5_ZROOT_spatially_from_weill/", 'atmbc.png'),
            dpi=300)



#%%
indice_veg = np.ones((np.shape(simu.DEM)), dtype=int)
raster_zone = np.ones((np.shape(simu.DEM)), dtype=int)
indice_zones = utils.root_updownhill_zones(raster_zone)

simu.update_veg_map(indice_veg)
simu.update_zone(indice_zones)

fig, ax = plt.subplots()
simu.show_input('root_map',ax=ax)
fig.savefig(os.path.join(simu.workdir,
                         simu.project_name,
                         simu.project_name+'root_map.png'
                         )
            , dpi=300)


#%%

df_SPP_map = simu.init_soil_SPP_map_df(nzones=2,nstr=15)
SPP_map = simu.set_SOIL_defaults(SPP_map_default=True)
df_FP_map = simu.init_soil_FP_map_df(nveg=dimveg)
FP_map = simu.set_SOIL_defaults(FP_map_default=True)

#%%

simu.update_soil(
                PMIN=-1e35,
                FP_map=FP_map,
                SPP_map=SPP_map,
                )

# simu.update_ic(INDP=0,
#                 IPOND=0,
#                 pressure_head_ini=-5 #.015
#                 )
simu.update_ic(INDP=4,
                IPOND=0,
                WTPOSITION=2 #.015
                )

simu.update_nansfdirbc(no_flow=True)
simu.update_nansfneubc(no_flow=True)
simu.update_sfbc(no_flow=True)

#%%
simu.update_parm(
                          TIMPRTi=list(t_atmbc),
                      )
# simu.workdir
# simu.project_name

#%%
simu.run_processor(
                    DTMIN=1e-2,
                    DELTAT=1,
                    DTMAX=1e3,
                    IPRT1=2,
                    TRAFLAG=0,
                    VTKF=2,
                    verbose=True
                   )

#%%

dtcoupling = simu.read_outputs('dtcoupling')

fig, ax = plt.subplots()
dtcoupling.plot(y='Atmpot-vf', ax=ax, color='k')
dtcoupling.plot(y='Atmact-vf', ax=ax, color='k', linestyle='--')
# ax.set_ylim(-1e-9,-5e-9)
ax.set_xlabel('Time (s)')
ax.set_ylabel('ET (m)')
plt.tight_layout()
fig.savefig(os.path.join(simu.workdir,
                         simu.project_name,
                         simu.project_name+'dtcoupling.png'
                         )
            , dpi=300)

# ETact_ref = dtcoupling

#%%

fig, axs = plt.subplots(4,4,
                        figsize=(8,5),
                        sharey=True,
                        sharex=True,
                        )
axs = axs.ravel()
for i in range(4*4):
    cmap = simu.show('spatialET',
                ax=axs[i],
                ti=i+2,
                scatter=True,
                vmin=5e-8,
                vmax=1e-7,
                colorbar=False
            )
    axs[i].set_ylabel('')
    axs[i].set_xlabel('')
    axs[i].set_title(f'#{i}')

cbar = fig.colorbar(cmap, ax=axs,
                    orientation='vertical',
                    fraction=0.02,
                    pad=0.04
                    )
cbar.set_label('ET (m/s)')

fig.savefig(os.path.join(simu.workdir,
                         simu.project_name,
                         simu.project_name+'spatialET.png'
                         )
            , dpi=300,
            bbox_inches='tight')

#%%

cplt.show_vtk(
    unit="pressure",
    timeStep=30,
    notebook=False,
    path=simu.workdir + "/ERA5_ZROOT_spatially_from_weill/vtk/",
    savefig=True,
)


#%% Read outputs
sw, sw_times = simu.read_outputs('sw')
df_psi = simu.read_outputs('psi')
# df_psi.index
#%% Choose 2 points uphill and downhill

depths = [0.05,0.15,0.25,0.75]

nodeIds,closest_positions = utils.get_mesh_node(simu,
                                        node_pos = [5,2.5,2],
                                        depths = depths
                                         )

nodeIds2,closest_positions2 = utils.get_mesh_node(simu,
                                        node_pos = [5,7.5,2] ,
                                        depths = depths
                                         )

#%%

pl = pv.Plotter(notebook=True)
mesh = pv.read(os.path.join(simu.workdir,
                                simu.project_name,
                                'vtk/121.vtk',
                                )
        )
pl.add_mesh(mesh,
            opacity=0.1
            )
pl.add_points(closest_positions,
              color='red'
              )
pl.add_points(closest_positions2,
              color='red'
              )
pl.show_grid()
pl.show()


#%%


# Generate a colormap
colors = plt.cm.tab10(np.linspace(0, 1, len(nodeIds)))

fig, ax = plt.subplots()

# Plot the first set of data with '.' markers
for i, nn in enumerate(nodeIds):
    ax.plot(sw.index / 86400, sw.iloc[:, nn], label=depths[i], marker='.', color=colors[i])

# Plot the second set of data with 'v' markers
for i, nn in enumerate(nodeIds2):
    ax.plot(sw.index / 86400, sw.iloc[:, nn], label=depths[i], marker='v', color=colors[i])

ax.set_xlabel('Time (Days)')
ax.set_ylabel('Saturation (-)')

# Create a unique legend
handles, labels = ax.get_legend_handles_labels()
unique_labels = {label: handle for handle, label in zip(handles, labels)}
ax.legend(unique_labels.values(), unique_labels.keys())

plt.title('ZROOT Zone 1 (.) / Zone 2 (v)')
fig.savefig(os.path.join(simu.workdir,
                         simu.project_name,
                         simu.project_name+'sw.png'
                         )
            , dpi=300)


#%%

# Generate a colormap
colors = plt.cm.tab10(np.linspace(0, 1, len(nodeIds)))

fig, ax = plt.subplots()

# Plot the first set of data with '.' markers
for i, nn in enumerate(nodeIds):
    ax.plot(df_psi.index / 86400, df_psi.iloc[:, nn], label=depths[i], marker='.', color=colors[i])

# Plot the second set of data with 'v' markers
for i, nn in enumerate(nodeIds2):
    ax.plot(df_psi.index / 86400, df_psi.iloc[:, nn], label=depths[i], marker='v', color=colors[i])

ax.set_xlabel('Time (Days)')
ax.set_ylabel('Pressure head (m)')

# Create a unique legend
handles, labels = ax.get_legend_handles_labels()
unique_labels = {label: handle for handle, label in zip(handles, labels)}
ax.legend(unique_labels.values(), unique_labels.keys())

plt.show()
plt.title('ZROOT Zone 1 (.) / Zone 2 (v)')

fig.savefig(os.path.join(simu.workdir,
                         simu.project_name,
                         simu.project_name+'psi.png'
                         )
            , dpi=300)
