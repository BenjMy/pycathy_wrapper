# -*- coding: utf-8 -*-
"""
@author: bmary
"""
from __future__ import print_function
import numpy as np

from pyCATHY import cathy_tools as CT
from pyCATHY.plotters import cathy_plots as pltCT
from pyCATHY import meshtools as mt

import matplotlib.pyplot as plt


import os
import pyvista as pv

def atmbc_PRD(workdir,project_name,dict_PRD,show=False,**kwargs):
    '''
    

    Parameters
    ----------
    workdir : TYPE
        DESCRIPTION.
    project_name : TYPE
        DESCRIPTION.
    dict_PRD : TYPE
        DESCRIPTION.
    show : TYPE, optional
        DESCRIPTION. The default is False.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    t_irr : TYPE
        DESCRIPTION.
    t_atmbc : TYPE
        DESCRIPTION.
    v_atmbc : TYPE
        DESCRIPTION.

    '''

    
    grid = dict_PRD['grid']
    time_drying0 = dict_PRD['time_drying0']
    irr_days = dict_PRD['irr_days']
    flux = dict_PRD['flux']
    lg_PRD = dict_PRD['lg_PRD']
    lg_PRD = dict_PRD['lg_PRD']

    
    x_min = max(grid["nodes_idxyz"][:, 1]) / 2
    x_max = max(grid["nodes_idxyz"][:, 1])
    y_min = min(grid["nodes_idxyz"][:, 2])
    y_max = max(grid["nodes_idxyz"][:, 2])

    nodes_xyz = grid["nodes_idxyz"]

    xmesh = grid["nodes_idxyz"][:, 1]
    ymesh = grid["nodes_idxyz"][:, 2]
    zmesh = grid["nodes_idxyz"][:, 3]

    HSPATM = 0
    IETO = 1

    # flux = 1.11111E-06
    totarea = 1  # 0.0

    sec_h = 3600
    days2sec = 24 * 3600



    time_drying0 = (time_drying0 * sec_h)  # days2sec # let the system dry out during 30days
    # irr_days = 1 # irrigate on left side during 10days
    # lg_PRD = 2 # 2 hours length of irrigation during the morning
    no_irr = 0
    t_irr = time_drying0  # initiate t_irr
    t_atmbc = [0]
    v_atmbc = [0]

    with open(os.path.join(workdir, project_name, "input/atmbc"), "w+") as atmbcfile:

        # write no irrigation during the first delta_drying0 days
        atmbcfile.write(
            str(HSPATM) + "\t" + str(IETO) + "\t" + "HSPATM" + "\t" + "IETO" + "\n"
        )
        atmbcfile.write(str(0) + "\t" + "TIME" + "\n")

        count_nmax = 0
        for i in range(int(grid["nnod3"])):
            if round(zmesh[i], 3) == round(max(zmesh), 3):
                count_nmax += 1
        print(count_nmax)

        # loop over the nodes of the mesh
        count = 0
        for i in range(int(grid["nnod3"])):
            if round(zmesh[i], 3) == round(max(zmesh), 3):
                # if i==max(range(int(grid['nnod3']))):
                if count == count_nmax:
                    # atmbcfile.write('0'.format('%f'))
                    atmbcfile.write(str(no_irr))
                else:
                    # atmbcfile.write('0'.format('%f')+ "\n")
                    atmbcfile.write(str(no_irr) + "\n")
            # atmbcfile.write('0'.format('%f')+ "\n")
            # atmbcfile.write(str(no_irr))

        # write irrigation once a day during 2h lasting 10 irr_days
        for k in range(irr_days):
            atmbcfile.write(str(t_irr) + "\t" + "TIME" + "\n")

            # loop over the nodes of the mesh
            count = 0
            for i in range(int(grid["nnod3"])):
                if round(zmesh[i], 3) == round(max(zmesh), 3):
                    count += 1
                    # loop over right or left side (odd or not)
                    if (k % 2) == 0:
                        if ((xmesh[i] > x_min) and (xmesh[i] < x_max)) and (
                            (ymesh[i] > y_min) and (ymesh[i] < y_max)
                        ):
                            # if count==count_nmax:
                            #     atmbcfile.write(str(flux/totarea))
                            # else:
                            atmbcfile.write(str(flux / totarea) + "\n")
                        else:
                            # if count==count_nmax:
                            #     atmbcfile.write(str(no_irr))
                            # else:
                            atmbcfile.write(str(no_irr) + "\n")

                    else:
                        if ((xmesh[i] > x_min) and (xmesh[i] < x_max)) and (
                            (ymesh[i] > y_min) and (ymesh[i] < y_max)
                        ):
                            # if count==count_nmax:
                            #     atmbcfile.write(str(no_irr))
                            # else:
                            atmbcfile.write(str(no_irr) + "\n")
                        else:
                            # if count==count_nmax:
                            #     atmbcfile.write(str(flux/totarea))
                            # else:
                            atmbcfile.write(str(flux / totarea) + "\n")

            t_atmbc.append(t_irr)
            t_irr += lg_PRD * sec_h
            v_atmbc.append(flux / totarea)

            # break irrigation during the rest of the day i.e. 24h -2h = 22h
            # atmbcfile.write("\n")
            atmbcfile.write(str(t_irr) + "\t" + "TIME" + "\n")
            t_atmbc.append(t_irr)
            # break during the resting hours
            count = 0
            for i in range(int(grid["nnod3"])):
                if round(zmesh[i], 3) == round(max(zmesh), 3):
                    count += 1
                    # if count==count_nmax:
                    #     atmbcfile.write(str(no_irr))
                    # else:
                    atmbcfile.write(str(no_irr) + "\n")

                    # atmbcfile.write(str(no_irr) +  "\n")

            v_atmbc.append(no_irr)
            t_irr += (24 - lg_PRD) * sec_h
            # t_atmbc.append(t_irr)
            # atmbcfile.write(str(no_irr) +  "\n")

        atmbcfile.close()

    if show == True:
        x_units = "sec"
        for key, value in kwargs.items():
            if key == "x_units":
                x_units = value
          
        pltCT.show_atmbc(
            t_atmbc, [v_atmbc, np.zeros(len(v_atmbc))], x_units=x_units)

    return t_irr, t_atmbc, v_atmbc


# class RHIZO(object):
#     '''Main RHIZO object.'''
# def __init__(self,dirName,prjName='my_cathy_prj'):

# rhizo.create_infitration(drippersPos=[],RWU=False)
#     # closestnode
#     # simu.create_parm()
#     # simu.create_soil()
#     # simu.create_ic()
#     # simu.create_atmbc()
#     # simu.create_nansfdirbc()

# rhizo.create_DA(drippersPos=[],RWU=False)
#     # closestnode
#     # simu.create_parm()
def update_rhizo_inputs(simu_DA, **kwargs):
    '''
    

    Returns
    -------
    None.

    '''
    
    # from pyCATHY.rhizo_tools import atmbc_PRD #create_infitration # create_DA

    pmin = -1.0E+20
    dem_synth=np.ones([9,10])*1e-3
    dem_synth[-1,-1]=0.99e-3
    zb = np.arange(0,0.5+0.1,0.05)
    nstr=len(zb)-1
    zr=list((np.ones(len(zb)))/(nstr))
    
    simu_DA.update_prepo_inputs(delta_x=0.05,delta_y=0.003, DEM=dem_synth,
                              N=np.shape(dem_synth)[1],
                              M=np.shape(dem_synth)[0],
                              N_celle=np.shape(dem_synth)[0]*np.shape(dem_synth)[1],
                              nstr=nstr, zratio=zr,  base=max(zb),dr=0.05,show=True)
    

    simu_DA.run_preprocessor(verbose=True,
                             KeepOutlet=False) #,show=True)
    
    simu_DA.update_parm()
    simu_DA.update_veg_map()
    simu_DA.update_soil()
    
   
    
    simu_DA.run_processor(IPRT1=3,
                          TRAFLAG=0,
                          verbose=False)

    
    
    
    if 'tobs' in kwargs:
        # time_of_interest=kwargs['tobs']
        time_of_interest= list(np.arange(0,2*1200,1200))
        # simu_DA.update_parm(TIMPRTi=time_of_interest,
        #                     TMAX=time_of_interest[-1])
        
    else:
        time_of_interest = list(np.arange(0,12*3600,3600))
    
    time_of_interest = list(np.arange(0,12*3600,3600))
    
    simu_DA.update_ic(INDP=3,WTPOSITION=0)
    flux=1.11111e-06
    
    
    # time_of_interest = list(np.arange(0,2000,1000))
    # np.arange(0,1,2)
    # fake time of interest to quickly create the mesh
    # time_of_interest = [0,0.1]
    # simu_DA.update_parm(TIMPRTi=time_of_interest,
    #                     TMAX=time_of_interest[-1])
            
    simu_DA.update_atmbc(HSPATM=1,IETO=1,TIME=time_of_interest,
                      VALUE=[np.zeros(len(time_of_interest)),np.ones(len(time_of_interest))*flux],
                      show=True,x_units='hours',diff=True) # just read new file and update parm and cathyH


    # input('press')
    # simu_DA.workdir
    # - here we just want to create the mesh


    

    """### 3- Boundary conditions"""
    
    simu_DA.update_nansfdirbc(noflow=True,TIME=time_of_interest) # Dirichlet uniform by default
    simu_DA.update_nansfneubc(TIME=time_of_interest) # Neumann (rain)
    simu_DA.update_sfbc(TIME=time_of_interest)
    
    
    """### 4- Soil and roots inputs"""
    
    PERMX = PERMY = PERMZ =[1.88E-04]
    ELSTOR = [1.00E-05]
    POROS = [0.55]
    VGNCELL = [1.46]
    VGRMCCELL = [0.15]
    VGPSATCELL =  [0.03125]
    
    SoilPhysProp = {'PERMX':PERMX,'PERMY':PERMY,'PERMZ':PERMZ,
                    'ELSTOR':ELSTOR,'POROS':POROS,
                    'VGNCELL':VGNCELL,'VGRMCCELL':VGRMCCELL,'VGPSATCELL':VGPSATCELL}
    
    """We could also used the predefined Van-Genuchten function (located in cathy_utils) to generate a set of soil physical properties"""
    
    #SoilPhysProp = utils.VanG()
    
    """The trick to defining a heterogeneous density of roots into the 3d medium from a single plant, 
    we defined different vegetation areas which are independent with varying root depth and Feddes parameters."""
    
    PCANA=[0.0]        
    PCREF=[-4.0]
    PCWLT=[-150]
    ZROOT=[0.2]
    PZ = [1.0]
    OMGC =[1.0] 
    
    
    FeddesParam = {'PCANA':PCANA,'PCREF':PCREF,'PCWLT':PCWLT,
                   'ZROOT':ZROOT,'PZ':PZ,
                   'OMGC':OMGC}
    
    veg_map = np.ones([int(simu_DA.hapin['M']),int(simu_DA.hapin['N'])])
    simu_DA.update_veg_map(root_depth=veg_map, show=True)
    
    """The choice of PMIN conditionne the switching condition"""
    
    simu_DA.update_soil(  PMINsimu_DA=pmin,
                          SPP=SoilPhysProp,
                          FP=FeddesParam,
                          verbose=True
                          
                        )


    
    # if 'tobs' in kwargs:
    #     print('tobs')
    #     simu_DA.update_parm(TIMPRTi=kwargs['tobs'],
    #                         TMAX=kwargs['TMAX'])
    #     time_of_interest=kwargs['tobs']
    # else:
    #     time_of_interest = list(np.arange(0,12*3600,3600))


    # if 'TMAX' in kwargs:
    #     simu_DA.update_parm()


    return simu_DA


#%% meshes interpolation

def CATHY_2_Simpeg(mesh_CATHY,mesh_Simpeg,scalar='saturation',show=False,**kwargs):
    pass

def CATHY_2_pg(mesh_CATHY,mesh_pg,scalar='saturation',show=False,**kwargs):
    
    
    
    # flip y and z axis as CATHY and pg have different convention for axis
    # ------------------------------------------------------------------------
    in_nodes_mod = np.array(mesh_CATHY.points)
    # in_nodes_mod_pg = np.array(mesh_pg.points)
    
    idx = np.array([0, 2, 1])
    in_nodes_mod_m = in_nodes_mod[:, idx]
        
    in_nodes_mod_m = in_nodes_mod[:, idx]
    in_nodes_mod_m[:,2] = -np.flipud(in_nodes_mod_m[:,2])
    in_nodes_mod_m[:,1] = -np.flipud(in_nodes_mod_m[:,1])


    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # ax.plot(in_nodes_mod_m[:,0],in_nodes_mod_m[:,1],in_nodes_mod_m[:,2])
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    
    
    # fig = plt.figure()
    # ax = fig.add_subplot(projection='3d')
    # ax.plot(in_nodes_mod_pg[:,0],in_nodes_mod_pg[:,1],in_nodes_mod_pg[:,2])
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')    
    
    
    
    path = os.getcwd()
    if 'path' in kwargs:
        path = kwargs['path']

    data_OUT = mt.trace_mesh(mesh_CATHY,mesh_pg,
                            scalar=scalar,
                            threshold=1e-1,
                            in_nodes_mod=in_nodes_mod_m)
    
    # print(data_OUT)
    scalar_new = scalar + '_nearIntrp2_pg_msh'
    # print(mesh_Resipy)
    # get_array(mesh, name, preference='cell'
              

    if 'time' in kwargs:
        time = kwargs['time']
        mesh_new_attr, name_new_attr = mt.add_attribute_2mesh(data_OUT, 
                                                                mesh_pg, 
                                                                scalar_new, 
                                                                overwrite=True,
                                                                time=time,
                                                                path=path)
    else:
        mesh_new_attr, name_new_attr = mt.add_attribute_2mesh(data_OUT, 
                                                                mesh_pg, 
                                                                scalar_new, 
                                                                overwrite=True,
                                                                path=path)
    
    if show == True:

        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
        p.add_mesh(mesh_new_attr,scalars=scalar_new)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
        
        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
        p.add_mesh(mesh_CATHY, scalars=scalar)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
        
    # if type(meshERT) is str:
    #     meshERTpv = pv.read(meshERT)
    
    # if savefig == True:
    
    #     plotter = pv.Plotter(notebook=True)
    #     _ = plotter.add_mesh(mesh_new_attr,show_edges=True)
    #     plotter.view_xz(negative=False)
    #     plotter.show_grid()
    #     plotter.save_graphic(path_CATHY + 'ERT' + str(DA_cnb) + str('.ps'), 
    #                           title='ERT'+ str(DA_cnb), 
    #                           raster=True, 
    #                           painter=True)
            
    
    return mesh_new_attr, scalar_new


def CATHY_2_Resipy(mesh_CATHY,mesh_Resipy,scalar='saturation',show=False,**kwargs):
 

    # flip y and z axis as CATHY and Resipy have different convention for axis
    # ------------------------------------------------------------------------
    in_nodes_mod = np.array(mesh_CATHY.points)
    in_nodes_mod[:,2] = -np.flipud(in_nodes_mod[:,2])
    in_nodes_mod[:,1] = -np.flipud(in_nodes_mod[:,1])


    # check with a plot positon of the nodes for both meshes
    # ------------------------------------------------------------------------
    # p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
    # p.add_mesh(mesh_CATHY)
    # _ = p.add_points(np.array(mesh_CATHY.points), render_points_as_spheres=True,
    #                         color='red', point_size=20)
    # _ = p.add_points(in_nodes_mod, render_points_as_spheres=True,
    #                         color='blue', point_size=20)
    # _ = p.show_bounds(grid='front', all_edges=True, font_size=50)
    # cpos = p.show(True)
    
    path = os.getcwd()
    if 'path' in kwargs:
        path = kwargs['path']

    data_OUT = mt.trace_mesh(mesh_CATHY,mesh_Resipy,
                            scalar=scalar,
                            threshold=1e-1,
                            in_nodes_mod=in_nodes_mod)
    
    # print(len(data_OUT))
    scalar_new = scalar + '_nearIntrp2Resipymsh'
    # print(mesh_Resipy)
    # get_array(mesh, name, preference='cell'
              

    if 'time' in kwargs:
        time = kwargs['time']
        mesh_new_attr, name_new_attr = mt.add_attribute_2mesh(data_OUT, 
                                                                mesh_Resipy, 
                                                                scalar_new, 
                                                                overwrite=True,
                                                                time=time,
                                                                path=path)
    else:
        mesh_new_attr, name_new_attr = mt.add_attribute_2mesh(data_OUT, 
                                                                mesh_Resipy, 
                                                                scalar_new, 
                                                                overwrite=True,
                                                                path=path)
    
    if show == True:

        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
        p.add_mesh(mesh_new_attr,scalars=scalar_new)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
        
        p = pv.Plotter(window_size=[1024*3, 768*2], notebook=True)
        p.add_mesh(mesh_CATHY, scalars=scalar)
        _ = p.add_bounding_box(line_width=5, color='black')
        cpos = p.show(True)
        
    # if type(meshERT) is str:
    #     meshERTpv = pv.read(meshERT)
    
    # if savefig == True:
    
    #     plotter = pv.Plotter(notebook=True)
    #     _ = plotter.add_mesh(mesh_new_attr,show_edges=True)
    #     plotter.view_xz(negative=False)
    #     plotter.show_grid()
    #     plotter.save_graphic(path_CATHY + 'ERT' + str(DA_cnb) + str('.ps'), 
    #                           title='ERT'+ str(DA_cnb), 
    #                           raster=True, 
    #                           painter=True)
            
    
    return mesh_new_attr, scalar_new





# -------------------------------------------------------------------#
#%% Infitration DATA

def create_infitration(self, dirfiles):
    self.set_drippers(dirfiles)

    pass

def set_drippers(self, dirfiles, drip_pos="drippers.txt"):

    print(os.getcwd())
    if isinstance(drip_pos, str):
        self.drippers = np.loadtxt(
            os.path.join(self.project_name, dirfiles, drip_pos),
            skiprows=1,
            delimiter=",",
        )
    else:
        self.drippers = drip_pos

    # check drippers position against DEM
    self.hapin
    mesh_x_max = float(self.hapin["xllcorner"]) + float(
        self.hapin["delta_x"]
    ) * float(self.hapin["N"])
    mesh_y_max = float(self.hapin["yllcorner"]) + float(
        self.hapin["delta_y"]
    ) * float(self.hapin["M"])

    mesh_x = []
    for xx in range(int(self.hapin["N"])):
        mesh_x.append(
            float(self.hapin["xllcorner"]) + float(self.hapin["delta_x"]) * xx
        )

    mesh_y = []
    for yy in range(int(self.hapin["M"])):
        mesh_y.append(
            float(self.hapin["yllcorner"]) + float(self.hapin["delta_y"]) * yy
        )

    print(mesh_x)
    print(mesh_y)

    print(mesh_x_max)
    print(max(self.drippers[:, 0]))

    if mesh_x_max < max(self.drippers[:, 0]):
        print(
            "Error: max mesh_x="
            + str(mesh_x_max)
            + "; max dripper x pos="
            + str(max(self.drippers[:, 0]))
        )

    if mesh_y_max < max(self.drippers[:, 1]):
        print(
            "Error: max mesh_x="
            + str(mesh_y_max)
            + "; max dripper y pos="
            + str(max(self.drippers[:, 1]))
        )

    # for dd in self.drippers:
    #     dd==

    #  (A==B).all()

    self.drippers_nodes = []

    pass
