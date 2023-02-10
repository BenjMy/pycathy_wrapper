"""
Weil et al example
==================

The CATHY gitbucket repository provides the Weill et al. dataset example to test the installation. On top of that, we provide a computational notebook code to reproduce the results using the **pyCATHY wrapper** (https://github.com/BenjMy/pycathy_wrapper). 

The notebook illustrate how to work interactively: execute single cell, see partial results at different processing steps (preprocessing, processing, output)... You can share it to work collaboratively on it by sharing the link and execute it from another PC without any installation required.


*Estimated time to run the notebook = 5min*

"""

# In[2]:

from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt

# In[13]:


path2prj = "weil_exemple"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj)

simu.run_processor(IPRT1=3, verbose=True)
simu.grid3d
len(simu.grid3d["mesh_tetra"])
simu.run_processor(IPRT1=2, verbose=True)

cplt.show_vtk(unit="pressure", timeStep=1, notebook=False, path="./my_cathy_prj/vtk/")
