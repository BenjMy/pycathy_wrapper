# -*- coding: utf-8 -*-
"""
@author: bmary
"""
from __future__ import print_function
import numpy as np
import cathy_tools as CT

import os

def atmbc_PRD(workdir,project_name,
              grid,
              flux=1.11111E-06,time_drying0=2,irr_days=1,lg_PRD=2):

    x_min = max(grid['nodes_idxyz'][:,1])/2
    x_max = max(grid['nodes_idxyz'][:,1])
    y_min = min(grid['nodes_idxyz'][:,2])
    y_max = max(grid['nodes_idxyz'][:,2])
    
    nodes_xyz =  grid['nodes_idxyz']
    
    xmesh = grid['nodes_idxyz'][:,1]
    ymesh = grid['nodes_idxyz'][:,2]
    zmesh = grid['nodes_idxyz'][:,3]   
    
    
    HSPATM=0
    IETO=1
    
    # flux = 1.11111E-06
    totarea=1 #0.0
    
    sec_h=3600
    days2sec=24*3600
    
    time_drying0 = time_drying0*days2sec # let the system dry out during 30days
    # irr_days = 1 # irrigate on left side during 10days
    # lg_PRD = 2 # 2 hours length of irrigation during the morning
    no_irr = 0
    t_irr =  time_drying0 # initiate t_irr
    t_atmbc = [0]
    v_atmbc = [0]
    
    with open(os.path.join(workdir , project_name, 'input/atmbc'), 'w+') as atmbcfile:
        
        # write no irrigation during the first delta_drying0 days
        atmbcfile.write(str(HSPATM) + "\t" + str(IETO) + "\t"
                        + 'HSPATM' + "\t" + 'IETO' + "\n")
        atmbcfile.write(str(no_irr) + "\t" + 'TIME' + "\n")
        for i in range(int(grid['nnod3'])):
             if round(zmesh[i], 3) == round(max(zmesh), 3):
                atmbcfile.write('0'.format('%f')+ "\n")
    
        # write irrigation once a day during 2h lasting 10 irr_days
        for k in range(irr_days):
            atmbcfile.write(str(t_irr) + 
                             "\t" + 'TIME' + "\n") 
            for i in range(int(grid['nnod3'])):
                if round(zmesh[i], 3) == round(max(zmesh), 3):
                    if (((xmesh[i]>x_min) and (xmesh[i]<x_max)) and ((ymesh[i]>y_min) and (ymesh[i]<y_max))):
                        atmbcfile.write(str(flux/totarea) + "\n")
                    else:
                        atmbcfile.write(str(no_irr)+ "\n")
            t_irr += lg_PRD*sec_h
            t_atmbc.append(t_irr)
            v_atmbc.append(flux/totarea)

            # break irrigation during the rest of the day i.e. 24h -2h = 22h
            atmbcfile.write("\n")
            atmbcfile.write(str(t_irr) + 
                             "\t" + 'TIME' + "\n")
            for i in range(int(grid['nnod3'])):
                 if round(zmesh[i], 3) == round(max(zmesh), 3):
                    atmbcfile.write(str(no_irr) +  "\n")
                                    
            t_irr += (24-lg_PRD)*sec_h
            t_atmbc.append(t_irr)
            v_atmbc.append(no_irr)
            
        atmbcfile.close()



    # # C Conto nodi per atmbc e li scrivo nel rispettivo file
    # for i in range(grid['nnod3']):
    #     if zmesh[i] == 0:
    #         if (((xmin[i]>x_min) and (xmin[i]<x_max)) and ((ymin[i]>y_min) and (ymin[i]<y_max))):
    #             counter=counter+1
    #             if (grid['nodes_idxyz'][:,0][i] == nn2d(i)) then
    #               write(25,*) nn(i),nn2d(i),areanod2d(i)
    #               totarea=totarea+areanod2d(i)
    #             end if
    #         end if
    #     end if
    #   end do
      # write(*,*) counter
      # write(*,*) totarea
  
  
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
    
    