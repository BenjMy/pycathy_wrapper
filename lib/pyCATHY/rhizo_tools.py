# -*- coding: utf-8 -*-
"""
@author: bmary
"""
from __future__ import print_function
import numpy as np

from pyCATHY import cathy_tools as CT
from pyCATHY.plotters import cathy_plots as pltCT


import os


def atmbc_PRD(
    workdir,
    project_name,
    grid,
    x_min=[],
    x_max=[],
    y_min=[],
    y_max=[],
    flux=1.11111e-06,
    time_drying0=2,
    irr_days=1,
    lg_PRD=2,
    show=False,
    **kwargs
):

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

    time_drying0 = (
        time_drying0 * sec_h
    )  # days2sec # let the system dry out during 30days
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
            print(count)

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

        pltCT.atmbc_inputs_plot(
            t_atmbc, [v_atmbc, np.zeros(len(v_atmbc))], x_units=x_units
        )

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
