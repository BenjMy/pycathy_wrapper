#!/usr/bin/env python
# coding: utf-8

# # Weill et al notebook CATHY exemple
# 
# The CATHY gitbucket repository provides the Weill et al. dataset example to test the installation. On top of that, we provide a computational notebook code to reproduce the results using the **pyCATHY wrapper** (https://github.com/BenjMy/pycathy_wrapper). 
# 
# The notebook illustrate how to work interactively: execute single cell, see partial results at different processing steps (preprocessing, processing, output)... You can share it to work collaboratively on it using the link and execute it from another PC without any installation required.
# 

# ## Step 0: package installation (not required if already installed)
# 
# ---
# 
# We first need to load the pyCATHY package. The pyCATHY package is calling fortran codes via python subroutines and ease the process

# In[1]:


get_ipython().run_cell_magic('capture', '', '!pip install git+https://github.com/BenjMy/pycathy_wrapper.git\n!pip install numpy GitPython pyvistaqt itkwidgets \n')


# In[2]:


get_ipython().run_cell_magic('capture', '', 'from pyCATHY import cathy_tools\nfrom pyCATHY.plotters import cathy_plots as cplt\n')


# In[15]:


#help(cathy_tools.CATHY)


# **We initiate a CATHY object**; 
# if the CATHY src files are not included within the 'path2prj' folder, they are automatically fetched from the gitbucket repository (providing the notebook is initiated with an internet connection). 
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

# In[3]:


path2prj ='weil_exemple' # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj,
                         clear_src=False, clear_outputs=True)


# # Step 1: preprocessing
# 
# ---
# We run the **preprocessor** using the command `simu.run_preprocessor` with the default files inputs (weill et al example). 
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

# In[4]:


help(simu.run_preprocessor)


# In[5]:


get_ipython().run_cell_magic('capture', '', 'simu.run_preprocessor(delta_x=0.1, verbose=False)\n')


# # Step 2: processing
# ---
# 
# 
# We run the **processor** with the default files inputs (weill et al example). This can last up to 3minutes

# In[6]:


help(simu.run_processor)


# In[17]:


simu.run_processor(IPRT1=2,verbose=True)


# # Visualisation
# ---
# Plot the results (for now using pyvista for 3d visualisation, need to update some dependencies before plotting)

# In[18]:


get_ipython().run_cell_magic('capture', '', "#!apt-get install -qq xvfb\n#!pip install pyvista panel -q\nimport os\n#os.system('/usr/bin/Xvfb :99 -screen 0 1024x768x24 &')\n#os.environ['DISPLAY'] = ':99'\n#import panel as pn\n#pn.extension('vtk')\n")


# In[9]:


help(cplt.show_vtk)
#print(os.getcwd())


# In[20]:


cplt.show_vtk(unit='pressure',
              timeStep=1,
              notebook=True, 
              path='./my_cathy_prj/vtk/'
             )

