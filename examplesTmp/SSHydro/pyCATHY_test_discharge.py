"""

"""

# In[2]:

from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt
import pyvista as pv

# In[13]:
prjname= 'Venezia_2011_asp_W_sceORIG_sceua_2l_v4_KGE_flat_test_colt2Rosetta222748'
# path2prj = "../examplesTmp/SSHydro/"  # add your local path here
path2prj = "."  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
			prj_name="Venezia_2011_asp_W_sceORIG_sceua_2l_v4_KGE_flat_test_colt2Rosetta222748"
			)


#%%


simu.run_preprocessor(verbose=False)
# simu.run_processor(IPRT1=3,verbose=True)

# # simu.read_inputs('atmbc')
# simu.update_parm(TIMPRTi=[1800,7200],
#                  VTKF=2
#                  )

# simu.atmbc
# simu.parm
# simu.update_parm(VTKF=2)

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
cplt.show_vtk(unit="saturation", 
              timeStep=10, 
              path=simu.workdir +  f"/{prjname}/vtk/",
              ax=pl,
              )
pl.show()


#%%
pl = pv.Plotter(notebook=True)
cplt.show_vtk(unit="pressure", 
              timeStep=1, 
              path=simu.workdir + f"/{prjname}/vtk/",
              ax=pl,
              )
pl.show()

#%%

cplt.show_vtk_TL(
                unit="pressure",
                notebook=False,
                path=simu.workdir +  f"/{prjname}/vtk/",
                show=False,
                x_units='days',
                # clim = [0.55,0.70],
                savefig=True,
            )

#%%

# simu.show(prop="hgsfdet")

# simu.show(prop="dtcoupling", yprop="Atmpot-d")
simu.show(prop="hgraph")



# simu.show(prop="cumflowvol")
