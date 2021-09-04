#!/usr/bin/env python
# coding: utf-8

# # Weill et al notebook CATHY exemple
# 
# The CATHY gitbucket repository provides the Weill et al. dataset example to test the installation. On top of that, we provide a computational notebook code to reproduce the results using the **pyCATHY wrapper** (https://github.com/BenjMy/pycathy_wrapper). 
# 
# The notebook illustrate how to work interactively: execute single cell, see partial results at different processing steps (preprocessing, processing, output)... You can share it to work collaboratively on it by sharing the link and execute it from another PC without any installation required.
# 

# We first need to load the pyCATHY package. The pyCATHY package is calling fortran codes via python subroutines and ease the process

# In[1]:


get_ipython().run_cell_magic('capture', '', '!pip install git+https://github.com/BenjMy/pycathy_wrapper.git')


# In[2]:


# import cathy package and dependencies
%%capture
get_ipython().system('pip install numpy GitPython pyvistaqt itkwidgets ')
from pyCATHY import cathy_tools
from pyCATHY import plot_tools as cplt


# We initiate a CATHY object; if the CATHY src files are not included within the 'path2prj' folder, they are automatically fetched from the gitbucket repository (providing the notebook is initiated with an internet connection). 
# 
# 
# 
# > Files are temporarily located by default in the 'my_cathy_prj' folder. We recommand to make a local backup copy of this folder or save it to the drive.
# 
# 
# 
# > Calling `cathy_tools.CATHY` returns an python-object `simu` where all the parameters are accessible through a python dict
# 
# 
# 

# In[ ]:


path2prj ='weil_exemple' # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj,
             clear_src=False, clear_outputs=True)


# # Step 1: preprocessing
# 
# ---
# We run the **preprocessor** `simu.run_preprocessor` with the default files inputs (weill et al example). 
# This step create all the prepro files required for the processor later on.
# 
# 
# > If you need to change preprocessor parameters, pass the variable as an argument and change its value to update the corresponding file. For instance I want to ***change the resolution x to 0.1m*** of the DEM: ``` # simu.run_preprocessor(delta_x=0.1) ```
# 
# 
# 
# 
# 
# 

# In[ ]:


simu.run_preprocessor(delta_x=0.1, verbose=False)


# # Step 2: processing
# ---
# 
# 
# We run the **processor** with the default files inputs (weill et al example). This can last up to 3minutes

# In[ ]:


simu.run_processor(verbose=False)


# # Visualisation
# ---
# Plot the results (for now using pyvista for 3d visualisation, need to update some dependencies before plotting)

# In[6]:


get_ipython().run_cell_magic('capture', '', "!apt-get install -qq xvfb\n!pip install pyvista panel -q\nimport os\nos.system('/usr/bin/Xvfb :99 -screen 0 1024x768x24 &')\nos.environ['DISPLAY'] = ':99'\nimport panel as pn\npn.extension('vtk')\nimport ipywidgets as widgets\nimport pyvista as pv")


# In[ ]:


cplt.showvtk(unit='pressure',timeStep=1,notebook=True, path='./my_cathy_prj/vtk/')

