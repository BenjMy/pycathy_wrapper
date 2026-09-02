"""
Output plots part 1
===================

Weill, S., et al. « Coupling Water Flow and Solute Transport into a Physically-Based Surface–Subsurface Hydrological Model ». 
Advances in Water Resources, vol. 34, no 1, janvier 2011, p. 128‑36. DOI.org (Crossref), 
https://doi.org/10.1016/j.advwatres.2010.10.001.

This example shows how to use pyCATHY object to plot the most common ouputs of the hydrological model.

*Estimated time to run the notebook = 5min*

"""


#%% Import packages
# Here we need to import `cathy_tools` class that control the CATHY core files preprocessing and processing
# We also import `cathy_plots` to render the results

from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt
import pyvista as pv
import os 
import matplotlib.pyplot as plt

#%% run processor
# if you add True to verbose, the processor log will be printed in the window shell
path2prj = "weil_exemple_outputs_plot1"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj)
simu.run_preprocessor()
simu.run_processor(IPRT1=2, 
                    DTMIN=1e-2,
                    DTMAX=1e2,
                    DELTAT=5,
                    TRAFLAG=0,
                    verbose=False
                    )



#%% read saturation file

df_sw, _ = simu.read_outputs('sw')
df_sw.head()

#%% Search for nearest mesh points

node, node_pos = simu.find_nearest_node([5,5,-1])
node2, node_pos2 = simu.find_nearest_node([5,5,1])
print(node_pos[0])

#%% plot node position on the mesh

pl = pv.Plotter(notebook=False)
cplt.show_vtk(unit="pressure", 
              timeStep=1, 
              path=os.path.join(simu.workdir,
                                simu.project_name,
                                'vtk'
                                ),
              style='wireframe',
              opacity=0.1,
              ax=pl,
              )
pl.add_points(node_pos[0],
              color='red'
              )
pl.add_points(node_pos2[0],
              color='red'
              )
pl.show()

#%% plot the saturation with time at a given mesh point

fig, ax = plt.subplots()
df_sw[node].plot(ax=ax)
df_sw[node2].plot(ax=ax)
ax.set_xlabel('time (s)')
ax.set_ylabel('saturation (-)')


#%% read pressure head file

df_psi = simu.read_outputs('psi')
# df_psi.head()
fig, ax = plt.subplots()
ax.plot(df_psi.index, df_psi.iloc[:,node[0]])
ax.plot(df_psi.index, df_psi.iloc[:,node2[0]])
ax.set_xlabel('time (s)')
ax.set_ylabel('pressure head (m)')
