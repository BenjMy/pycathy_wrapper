"""
Output plots part 2
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



#%% run processor
# if you add True to verbose, the processor log will be printed in the window shell



path2prj = "../SSHydro/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
			prj_name="weil_exemple_outputs_plot"
			)

simu.run_preprocessor()
simu.run_processor(IPRT1=2, 
                    DTMIN=1e-2,
                    DTMAX=1e2,
                    DELTAT=5,
                   TRAFLAG=0,
                   verbose=False
                   )


#%% plot NET SEEPFACE VOL and NET SEEPFACE FLX over the time t
simu.show(prop="hgsfdet")

#%% plot Atmact-vf = f (time)
simu.show(prop="dtcoupling", yprop="Atmpot-d")

#%% Another interesting graph looking at the **streamflow = f(time)**
simu.show(prop="hgraph")

#%% Plot the "Total flow volume" and the "nansfdir flow volume" = f(time)
simu.show(prop="cumflowvol")

#%% 3d visualiation of the pressure head for the time step 1
# To select another time step change the value in the function argument
cplt.show_vtk(
    unit="pressure",
    timeStep=1,
    notebook=False,
    path=simu.workdir + "/weil_exemple_outputs_plot/vtk/",
)

#%%  3d visualiation of the water saturation for the time step 1
# cplt.show_vtk(
#     unit="saturation",
#     timeStep=1,
#     notebook=False,
#     path=simu.workdir + "/my_cathy_prj/vtk/",
# )


