"""
Test surface-subsurface hydrology modules
"""


# Weil et al example
# ==================
from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt

path2prj = "weil_exemple"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj)
simu.run_processor(IPRT1=3, verbose=True)
simu.run_processor(IPRT1=2, verbose=True)

cplt.show_vtk(unit="pressure", timeStep=1, notebook=False, path="./my_cathy_prj/vtk/")
