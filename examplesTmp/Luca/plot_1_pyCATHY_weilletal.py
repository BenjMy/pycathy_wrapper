"""
Weill et al example
==================

Weill, S., et al. « Coupling Water Flow and Solute Transport into a Physically-Based Surface–Subsurface Hydrological Model ». 
Advances in Water Resources, vol. 34, no 1, janvier 2011, p. 128‑36. DOI.org (Crossref), 
https://doi.org/10.1016/j.advwatres.2010.10.001.

The CATHY gitbucket repository provides the Weill et al. dataset example to test the installation. On top of that, we provide a computational notebook code to reproduce the results using the **pyCATHY wrapper** (https://github.com/BenjMy/pycathy_wrapper). 

The notebook illustrate how to work interactively: execute single cell, see partial results at different processing steps (preprocessing, processing, output)... You can share it to work collaboratively on it by sharing the link and execute it from another PC without any installation required.


*Estimated time to run the notebook = 5min*

"""

# In[2]:

from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt
import pyvista as pv

# In[13]:

path2prj = "../SSHydro/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
			prj_name="weill_exemple"
			)


#%%


simu.run_preprocessor(verbose=False)
# simu.run_processor(IPRT1=3,verbose=True)

# simu.read_inputs('atmbc')
simu.update_parm(TIMPRTi=[1800,7200])

# simu.atmbc
simu.parm


# simu.grid3d
# len(simu.grid3d["mesh_tetra"])
simu.run_processor(IPRT1=2, 
                    DTMIN=1e-2,
                    DTMAX=1e2,
                    DELTAT=5,
                   TRAFLAG=0,
                   verbose=False
                   )

#%%

pl = pv.Plotter(notebook=False)
cplt.show_vtk(unit="pressure", 
              timeStep=1, 
              path=simu.workdir + "/weill_exemple/vtk/",
              ax=pl,
              )
pl.show()


#%%
pl = pv.Plotter(notebook=True)
cplt.show_vtk(unit="pressure", 
              timeStep=1, 
              path=simu.workdir + "/weill_exemple/vtk/",
              ax=pl,
              )
pl.show()



