"""
Nb of outlets >1
==================

The notebook illustrate how to work interactively: execute single cell, see partial results at different processing steps (preprocessing, processing, output)... You can share it to work collaboratively on it by sharing the link and execute it from another PC without any installation required.


*Estimated time to run the notebook = 5min*

"""

# In[2]:

from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt
import pyvista as pv

# In[13]:

path2prj = "."  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
			prj_name="outlets"
			)

#%%
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Define the dimensions of the DEM
dem_size = (15, 15)

# Create a synthetic DEM with a general slope
x = np.linspace(0, 1, dem_size[0])
y = np.linspace(0, 1, dem_size[1])
X, Y = np.meshgrid(x, y)
slope = np.exp(Y)**2  # Elevation decreases from top to bottom

# Define the outlet location and create a catchment
outlet_row, outlet_col = 14, 7  # Position of the single outlet (bottom row)

# Create the DEM with a single outlet
dem = slope
# dem[outlet_row, outlet_col] = 0.2  # Lower elevation at the outlet

# Create a catchment area around the outlet (higher elevation in the catchment)
catchment_radius = 4
for i in range(dem_size[0]):
    for j in range(dem_size[1]):
        if np.sqrt((i - outlet_row)**2 + (j - outlet_col)**2) < catchment_radius:
            dem[i, j] = min(dem[i, j], dem[outlet_row, outlet_col] + 0.2)  # Elevation increases towards the outlet

# Create an xarray DataArray for the DEM
dem_xr = xr.DataArray(dem, dims=['y', 'x'], coords={'y': np.arange(dem_size[0]), 'x': np.arange(dem_size[1])})

# Plot the synthetic DEM
fig, ax = plt.subplots(figsize=(8, 6))
dem_xr.plot.imshow(ax=ax, cmap='terrain')
ax.set_title('Synthetic DEM with Single Outlet')
plt.show()

#%%

simu.update_ic(INDP=0,
               pressure_head_ini=-50
               )
#%%

simu.update_prepo_inputs(
                        DEM=dem_xr.values,
                        delta_x = 1,
                        delta_y = 1,
                        base=6,
                        )

#%%


simu.run_preprocessor(verbose=False)

simu.create_mesh_vtk()

#%%
# simu.run_processor(IPRT1=3,verbose=True)
times = [0,1800,7200,7200*5]

simu.read_inputs('atmbc')
simu.update_atmbc(
                HSPATM=1,
                IETO=1,
                time=times,
                netValue = [0,4e-7,0,0]
                )

# # simu.atmbc
# simu.parm

simu.update_parm(
                IPRT=4,
                VTKF=4,
                )

#%%


# simu.update_mesh_boundary_cond()

simu.create_mesh_bounds_df(
                            'nansfdirbc', 
                            simu.grid3d["mesh3d_nodes"], 
                            times=times, 
                            )

print(simu.mesh_bound_cond_df.head())
print(simu.mesh_bound_cond_df.head())

# np.sum(simu.mesh_bound_cond_df['ymin_bound'])

simu.update_nansfneubc(
                        time = times,
                        ymin_bound=1.2e-8
                        )


#%% Show Boundary Conitions at time 0
# simu.show_bc(time=0)
# meshbc = simu.mesh_bound_cond_df
# cplt.plot_mesh_bounds('nansfneubc',meshbc, time=0,
#                       vmin=0)

#%%
# # simu.grid3d
# # len(simu.grid3d["mesh_tetra"])
simu.run_processor(IPRT1=2, 
                    DTMIN=1e-2,
                    DTMAX=1e2,
                    DELTAT=5,
                    TRAFLAG=0,
                    verbose=False
                    )

#%%

pl = pv.Plotter(notebook=False)
cplt.show_vtk(unit="saturation", 
              timeStep=1, 
              path=simu.workdir + "/outlets/vtk/",
              ax=pl,
              )
pl.show()


#%%

cplt.show_vtk_TL(
                unit="saturation",
                notebook=False,
                path=simu.workdir + "/outlets/vtk/",
                show=True,
                x_units='days',
                clim = [0,1],
                savefig=True,
            )

#%%
simu.show(prop="hgraph")



