# -*- coding: utf-8 -*-
"""
@author: bmary
"""
from __future__ import print_function
import sys
import numpy as np
import os
import warnings
import subprocess
import glob
from os import listdir
from os.path import isfile, join
import shutil

from collections import OrderedDict

import git
from git import Repo

import pandas as pd
# import resipy
# from resipy import Project

# import meshtools as mt
# import pyCATHY

import time
import rich.console
from rich.progress import track

# import cathy_plots as pltCT
# from pltCT import *

from pyCATHY.plotters import cathy_plots as plt_CT
# from pyCATHY.DA.cathy_DA import DA
from pyCATHY.DA import enkf


from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import sensors_measures as in_meas

from pyCATHY.ERT import petro_Archie as Archie


from functools import partial
import multiprocessing

import matplotlib.pyplot as plt

def make_console(verbose):
    """
    Start up the :class:`rich.console.Console` instance we'll use.

    Parameters
    ----------
    verbose : bool
        Whether or not to print status messages to stderr.
    """
    return rich.console.Console(stderr=True, quiet=not verbose)



# multiprocessor functions outside main CATHY object
# -----------------------------------------------------------------------------
def subprocess_run_multi(pathexe_list):
    '''
    

    Parameters
    ----------
    pathexe_list : TYPE
        DESCRIPTION.

    Returns
    -------
    p : TYPE
        DESCRIPTION.

    '''
    # https://stackoverflow.com/questions/44144584/typeerror-cant-pickle-thread-lock-objects
    print(f"x= {pathexe_list}, PID = {os.getpid()}")
    # self.console.print(":athletic_shoe: [b]nudging type: [/b]" + str(self.DAFLAG))

    os.chdir(pathexe_list)
    callexe = "./" + 'cathy'
    p = subprocess.run([callexe], text=True, capture_output=True)

    return p


def ERTmapping_run_multi(workdir,
                         dict_obs, 
                         count_DA_cycle, 
                         porosity, 
                         # path_CATHY,
                         forward_mesh_vtk_file,
                         output_dirname,
                         project_name):
    '''
    

    Parameters
    ----------
    list : TYPE
        DESCRIPTION.

    Returns
    -------
    p : TYPE
        DESCRIPTION.

    '''

    # not yet implemented\
        
    # case 4: ERT
    # --------------------------------------------------------------------
    # if key[0] in 'ERT':
    pathERT = os.path.split(dict_obs[0]['filename'])[0]
    print('*** ERT MAPPING // run ****')
    print(f"x= {workdir}, PID = {os.getpid()}")

    # print(count_DA_cycle)
    forward_mesh_vtk_file = glob.glob(
    pathERT + "/**/*forward_model*", recursive=True)
    # electrodes = glob.glob(pathERT + "/**/*electrodes.dat", recursive = True)
    electrodes = glob.glob(
    pathERT + "/**/*elecsXYZ.csv", recursive=True)
    # seq = glob.glob(pathERT + "/**/*protocol.dat", recursive = True)
    seq = glob.glob(pathERT + "/**/*ERT72.csv", recursive=True)
    
    # df_sw = CATHY.read_outputs(filename='sw',
    #                   path=os.path.join(os.getcwd(),
    #                                     output_dirname))
    
    df_sw = out_CT.read_sw(os.path.join(workdir,'output/sw'))
    # # df_sw[-1]
    
    # path_CATHY = os.path.join(os.getcwd(), 'vtk/')
    
    # porosity = CATHY.soil['SPP'][:, 4][0]
    Hx_ERT, df_Archie = Archie.SW_2_ERa(df_sw, 
                                        workdir, 
                                        project_name,
                                        porosity,
                                        # path_CATHY,
                                        pathERT,
                                        meshERT=forward_mesh_vtk_file[0],
                                        elecs=electrodes[0],
                                        sequenceERT=seq[0],
                                        DA_cnb=count_DA_cycle)
    
    
    # print(Hx_ERT)
    # self._add_2_ensemble_Hx(Hx, Hx_ERT['resist'])

    # p = []

    return Hx_ERT

    
# -----------------------------------------------------------------------------

class CATHY():  # IS IT GOOD PRACTICE TO PASS DA CLASS HERE ? I think we sould better pass main CATHY into children classes
    """ Main CATHY object

      When instantiated it creates the tree project directories with 'prjName' as root folder.
      The src files are fetched from the online repository if not existing
      (note that it is possible to call a specific version).

      """

    def __init__(self, dirName=None, prjName="my_cathy_prj", notebook=False,
                 version="1.0.0",
                 verbose=True,
                 **kwargs):
        """
        Create CATHY object.


        ..note:
            All variables in CAPITAL LETTERS use the semantic from the CATHY legacy fortran codes,
            while the others are created exclusively for the wrapper.

        ..note:
            All file variable are self dictionnary objects; example: soil file is contained in self.soil
        """

        # create project dir
        # ---------------------------------------------------------------------
        self.console = make_console(verbose)
        self.console.print(":checkered_flag: [b]Initiate CATHY object[/b]")

        # flag for notebook execution
        # ---------------------------------------------------------------------
        self.notebook = notebook  # flag if the script is run in a notebook

        # create working and project dir
        # ---------------------------------------------------------------------
        if dirName is None:
           self.workdir = os.path.join(os.getcwd())
        else:
            self.workdir = os.path.join(os.getcwd(), dirName)

        # change working directory
        if not os.path.exists(os.path.join(self.workdir)):
            os.makedirs(self.workdir, exist_ok=True)
        os.chdir(self.workdir)

        self.project_name = prjName

        if not os.path.exists(os.path.join(self.workdir, self.project_name)):
            os.makedirs(os.path.join(
                self.workdir, self.project_name), exist_ok=True)

        # convention for directory tree
        # ---------------------------------------------------------------------
        self.processor_name = "cathy"
        self.input_dirname = "input"
        self.output_dirname = "output"

        # dict related to CATHY inputs files
        # ---------------------------------------------------------------------
        self.parm = {}  # dict of parm input parameters
        self.soil = {}  # dict of soil input parameters
        self.ic = {}  # dict of ic input parameters
        self.cathyH = {}  # dict of cathyH C header parameters
        # self.nudging = {'NUDN': 0}   # Temporary

        # dict related to Data assimilation
        # ---------------------------------------------------------------------
        self.DAFLAG = False # Flag to trigger data assimilation process
        self.dict_obs = OrderedDict() # dictionnary containing all the observation data
        self.count_DA_cycle = None  # counter of the Assimilation Cycle
        
        # dict related to the mesh
        # ---------------------------------------------------------------------
        self.grid3d = {}

        for key, value in kwargs.items():

        # clear src files if required
        # ---------------------------------------------------------------------
            if key == "clear_src":  # clear src files
                if value == True:
                    if os.path.exists(
                        os.path.join(self.workdir, self.project_name, "src")
                    ):
                        print("clear src files")
                        shutil.rmtree(
                            os.path.join(
                                self.workdir, self.project_name, "src")
                        )
                        shutil.rmtree(
                            os.path.join(
                                self.workdir, self.project_name, "input")
                        )
                        shutil.rmtree(
                            os.path.join(
                                self.workdir, self.project_name, "tmp_src")
                        )
                        shutil.rmtree(
                            os.path.join(
                                self.workdir, self.project_name, "prepro")
                        )

        # fetch src files if not existing from Gitbucket repo
        # ---------------------------------------------------------------------
        if not os.path.exists(os.path.join(self.project_name, "src")):
            self.console.print(":worried_face: [b]src files not found[/b]")
            self.console.print("working directory is:" + str((self.workdir)))

            if version == "1.0.0":
                try:
                    Repo.clone_from(
                        "https://bitbucket.org/cathy1_0/cathy.git",
                        os.path.join(
                            self.workdir, self.project_name, "tmp_src"),
                        branch="master",
                    )
                    self.console.print(
                        ":inbox_tray: [b]Fetch cathy src files[/b]")
                    shutil.move(
                        os.path.join(
                            self.workdir, self.project_name, "tmp_src/src"),
                        os.path.join(self.workdir, self.project_name, "src"),
                    )

                    self.console.print(
                        ":inbox_tray: [b]Fetch cathy prepro src files[/b]")
                    shutil.move(
                        os.path.join(
                            self.workdir,
                            self.project_name,
                            "tmp_src/runs/weilletal/prepro",
                        ),
                        os.path.join(
                            self.workdir, self.project_name, "prepro"),
                    )

                    self.console.print(
                        ":inbox_tray: [b]Fetch cathy inputfiles[/b]")
                    shutil.move(
                        os.path.join(
                            self.workdir,
                            self.project_name,
                            "tmp_src/runs/weilletal/input",
                        ),
                        os.path.join(self.workdir, self.project_name, "input"),
                    )

                    pathsrc = os.path.join(
                        os.getcwd(), self.project_name, "tmp_src/runs/weilletal/"
                    )

                    onlyfiles = [
                        f for f in listdir(pathsrc) if isfile(join(pathsrc, f))
                    ]

                    for (
                        file
                    ) in (
                        onlyfiles
                    ):  # You could shorten this to one line, but it runs on a bit.
                        shutil.move(
                            os.path.join(pathsrc, file),
                            os.path.join(self.project_name, file),
                        )
                except:
                    print("no internet connection to fetch the files")
                    sys.exit()
                    pass

            if version == "G. Manoli":
                print("fetch cathy G. Manoli src files")
                # /home/ben/Documents/CATHY/CathyGitbucket/Test_Gabriele/1_Gabriele_Piante_NON_modificato/CATHY_RWU_ABL_1D/
                path_manoli = "/home/ben/Documents/CATHY/CathyGitbucket/Test_Gabriele/1_Gabriele_Piante_NON_modificato/CATHY_RWU_ABL_1D/"
                shutil.copytree(
                    path_manoli, os.path.join(
                        self.workdir, self.project_name, "src")
                )

        # clear output files if required
        # ---------------------------------------------------------------------
        for key, value in kwargs.items():
            if key == "clear_outputs":
                if value == True:
                    if not os.path.exists(
                        os.path.join(self.workdir, self.project_name, "output")
                    ):
                        self.create_output(output_dirname="output")
                    else:
                        shutil.rmtree(
                            os.path.join(
                                self.workdir, self.project_name, "output")
                        )
                        shutil.rmtree(
                            os.path.join(
                                self.workdir, self.project_name, "vtk")
                        )
                        self.create_output(output_dirname="output")

        pass


# --------------------------------------------------------------------------- #
# replace console cmd by a subprocess call
# --------------------------------------------------------------------------- #


    def run_preprocessor(self, KeepOutlet=True, verbose=False, **kwargs):
        """
        Run cppp.exe

        1. updates cathy parameters prepo

        2. recompile the source files prepo

        3. run the preprocessor using bash cmd

        Returns -------

        Running the executable file has allowed to generate a complete set of
        files describing physiographic features of the drainage system, as shown in Table 2.

        """

        for loopi in range(2):  # run it twice (to avoid the first error)
            os.chdir(os.path.join(self.workdir, self.project_name, "prepro/src/"))

            # clean all files compiled
            for file in glob.glob("*.o"):
                os.remove(file)

            # if self.notebook==False:
            bashCommand = "gfortran -O -o pycppp mpar.f90 mbbio.f90 wbb_sr.f90 csort.f90 qsort.f90 depit.f90 cca.f90 smean.f90 dsf.f90 facet.f90 hg.f90 mrbb_sr.f90 bb2shp_sr.f90 shape.f90 dbase.f90 streamer.f90 cppp.f90"

            try:
                # bashCommand = 'gcc -O -o pycppp mpar.f90 mbbio.f90 wbb_sr.f90 csort.f90 qsort.f90 depit.f90 cca.f90 smean.f90 dsf.f90 facet.f90 hg.f90 mrbb_sr.f90 bb2shp_sr.f90 shape.f90 dbase.f90 streamer.f90 cppp.f90'

                p = os.system(bashCommand)
                # run it twice (to avoid the first error)
                p = os.system(bashCommand)

                # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                # process.communicate()
                # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                # process.communicate()

                # if verbose:
                #     output, error = process.communicate()

                # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

                # if verbose:
                #     output, error = process.communicate()
                self.console.print(":cooking: [b]gfortran compilation[/b]")

            except:
                print("bash cmd not recognized")

            os.chdir(self.workdir)

        # try:
        # move to the directory where the source FORTRAN files are contained (cathy_main.f)
        shutil.move(
            os.path.join(self.workdir, self.project_name, "prepro/src/pycppp"),
            os.path.join(self.workdir, self.project_name, "prepro/pycppp"),
        )

        self.console.print(":athletic_shoe: [b]Run preprocessor[/b]")
        os.chdir(os.path.join(self.workdir, self.project_name, "prepro"))

        bashcmd = "./pycppp"
        my_data = "2\n0\n1\n"  # user input data
        p = subprocess.run([bashcmd], text=True,
                           input=my_data, capture_output=True)

        # try:
        #     p = subprocess.run([bashcmd], text=True, input=my_data, capture_output=True)
        # except:
        #     shutil.move(os.path.join(self.workdir,self.project_name, 'prepro/src/cppp'),
        #     os.path.join(self.workdir,self.project_name, 'prepro/cppp'))
        #     p = subprocess.run([bashcmd], text=True, input=my_data, capture_output=True)

        if verbose:
            print(p.stdout)
            print(p.stderr)

        if "catchment with more than one outlet cell!" in p.stdout:
            print("catchment with more than one outlet cell!")
            sys.exit()

        if KeepOutlet == False:
            print("remove outlet")
            self.DEM[0] = self.DEM[1]
            with open(
                os.path.join(self.workdir, self.project_name,
                             "prepro/dtm_13.val"), "w+"
            ) as f:
                np.savetxt(f, self.DEM, fmt="%1.4e")  # use exponential

        os.chdir(self.workdir)

        # -------------------------------------------- #
        # run processor only to build the 3d mesh
        # self.run_processor(verbose=True,IPRT1=2,TRAFLAG=0)
        # self.read_grid3d()
        # -------------------------------------------- #

        return

    def recompileSrc(self, verbose=False):
        """
        Other option is self.run_processor(runProcess=False)

        """

        self.console.print(":hammer_and_wrench: [b]Recompile src files[/b]")

        # clean all files previously compiled
        for file in glob.glob("*.o"):
            os.remove(file)

        # list all the fortran files to compile and compile
        for file in glob.glob("*.f"):
            bashCommand = "gfortran -c " + str(file)
            # gfortran -c *.f
            # gfortran *.o -L\MinGW\lib -llapack -lblas -o cathy

            process = subprocess.Popen(
                bashCommand.split(), stdout=subprocess.PIPE)

            # if verbose==True:
            #     output, error = process.communicate()
            output, error = process.communicate()

        # list all the fortran compiled files to compile and run
        files = ""
        for file in glob.glob("*.o"):
            files += " " + str(file)

        bashCommand = (
            "gfortran" + files + " -llapack -lblas -o " + self.processor_name
        )
        self.console.print(":cooking: [b]gfortran compilation[/b]")

        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        # if verbose==True:
        #     output, error = process.communicate()

        # os.chdir(os.path.join(self.workdir))
        # try to move created executable to the project folder
        try:
            shutil.move(
                os.path.join(self.workdir, self.project_name, "src", self.processor_name
                ),
                os.path.join(self.workdir, self.project_name,
                             self.processor_name)
            )
        except:
            # print("cannot find the new processsor:" + str(self.processor_name))
            self.console.print(
                ":pensive_face: [b]Cannot find the new processsor[/b]")

        pass

    def run_processor(self, recompile=True, runProcess=True, verbose=False, **kwargs):
        """
        Run cathy.exe


        1. updates cathy parameters and cathyH based on **kwargs

        2. recompile the source files (set False for notebook) and create the executable

        3. run the processor using bash cmd if runProcess is True

        Returns ------- type Description of returned object.

        """

        # set parallel flag for mutirun (used for DA)
        # --------------------------------------------------------------------
        parallel = False
        # for k, value in kwargs.items():
        #     if k == 'parallel':
        #         parallel = value
        if 'parallel' in kwargs:
                parallel = kwargs['parallel']

        if 'DAFLAG' in kwargs:
            self.DAFLAG = kwargs['DAFLAG']

        # update parm and cathyH
        # --------------------------------------------------------------------
        if len(kwargs.items()) > 0:  # if new arguments update parm and CATHYH
            self.update_parm(**kwargs, verbose=verbose)
            self.update_cathyH(**kwargs, verbose=verbose)

        # recompile
        # --------------------------------------------------------------------
        os.chdir(os.path.join(self.workdir, self.project_name, "src"))

        if recompile == True:

            self.recompileSrc()
            # self.console.print(":hammer_and_wrench: [b]Recompile src files[/b]")

        # back to root dir
        os.chdir(os.path.join(self.workdir, self.project_name))

        # some checkings
        # --------------------------------------------------------------------
        try:
            os.path.exists(os.path.join("inputs"))
        except OSError:
            print("input folder missing")
            sys.exit()

        # check if output folder
        if not os.path.exists("output"):
            os.makedirs("output", exist_ok=True)

        # check if vtk folder
        if not os.path.exists("vtk"):
            os.makedirs("vtk", exist_ok=True)

        # run the processor
        # --------------------------------------------------------------------
        if runProcess == True:

            t0 = time.time()  # executation time estimate

            self.console.print(":athletic_shoe: [b]Run processor[/b]")
            callexe = "./" + self.processor_name

            # case of Data Assimilation DA
            # ----------------------------------------------------------------
            if self.DAFLAG:

                # define Data Assimilation parameters
                # -------------------------------------------------------------

                # type of assimillation
                DA_type = 'Sakov'
                if 'DA_type' in kwargs:
                    DA_type = kwargs['DA_type']

                # perturbated parameters dict
                dict_parm_pert = []
                if 'dict_parm_pert' in kwargs:
                    dict_parm_pert = kwargs['dict_parm_pert']

                # observations dict
                dict_obs = []
                if 'dict_obs' in kwargs:
                    dict_obs = kwargs['dict_obs']

                # what needs to be updated ? only model states or also model parameters?
                list_update_parm = ['St. var.']
                if 'list_update_parm' in kwargs:
                    list_update_parm = kwargs['list_update_parm']

                # what needs to be assimilated ? default is ALL
                list_assimilated_obs = []
                if 'list_assimilated_obs' in kwargs:
                    list_assimilated_obs = kwargs['list_assimilated_obs']

                self._run_DA(callexe, parallel, verbose,
                             DA_type,
                             dict_obs,
                             list_update_parm,
                             dict_parm_pert,
                             list_assimilated_obs
                             )

            # case of simple simulation
            # ----------------------------------------------------------------
            else:

                p = subprocess.run([callexe], text=True, capture_output=True)

                if verbose == True:
                    print(p.stdout)
                    print(p.stderr)

                os.chdir(os.path.join(self.workdir))

            t1 = time.time()
            self.total = t1 - t0

                # try:
                #    subprocess.call([callexe_path, '-ARG'], shell=True)
                #    subprocess.Popen(["cathyEnv/bin/python"])
                # except:
                #    pass
                # process = subprocess.Popen(callexe.split(), stdout=subprocess.PIPE)
                # if verbose:
                #     output, error = process.communicate()
                # print('timer')

        return

    def _run_DA(self, callexe, parallel, verbose,
                        DA_type,
                        dict_obs,
                        list_update_parm,
                        dict_parm_pert,
                        list_assimilated_obs
                        ):
        '''
        

        Parameters
        ----------
        callexe : TYPE
            DESCRIPTION.
        parallel : TYPE
            DESCRIPTION.
        verbose : TYPE
            DESCRIPTION.
        DA_type : TYPE
            DESCRIPTION.
        dict_obs : TYPE
            DESCRIPTION.
        list_update_parm : TYPE
            DESCRIPTION.
        dict_parm_pert : TYPE
            DESCRIPTION.
        list_assimilated_obs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        # check that prepare DA has been properly done
        # to write

        # Initiate
        # -------------------------------------------------------------------
        self.dict_obs = dict_obs
        self.dict_parm_pert = dict_parm_pert
        update_key = 'ini_perturbation'
        
        # Infer ensemble size from perturbated parameter dictionnary
        # -------------------------------------------------------------------
        for name in self.dict_parm_pert:
            NENS = len(self.dict_parm_pert[name]['ini_perturbation'])

        # Infer ensemble update times from observation dictionnary
        # -------------------------------------------------------------------
        ENKFT = []
        for ti in self.dict_obs:
            ENKFT.append(float(ti))
            
        if self.parm['TMAX']> max(ENKFT):
            ENKFT.append(self.parm['TMAX'])
        else:
            exit()
        
        # print(ENKFT)
        # sys.exit()
        # start DA cycle counter
        # -------------------------------------------------------------------
        self.count_DA_cycle = 0
        # the counter is incremented during the update analysis

        # initiate DA
        # -------------------------------------------------------------------
        self._DA_init(NENS=NENS, # ensemble size
                      ENKFT=ENKFT, # assimilation times
                      parm_pert=dict_parm_pert)

        # Loop over assimilation cycle
        # -------------------------------------------------------------------
        for DA_cnb in range(len(ENKFT)):
            self.console.print('count DA cycle nb: ' +
                               str(self.count_DA_cycle))
            
            # multi run from the independant folders composing the ensemble
            # ----------------------------------------------------------------
            if parallel == True:
                pathexe_list = []
                for ens_i in range(self.NENS):
                    path_exe = os.path.join(self.workdir,
                                            self.project_name,
                                            'DA_Ensemble/cathy_' + str(ens_i+1))
                    pathexe_list.append(path_exe)

                with multiprocessing.Pool(processes=4) as pool:
                    result = pool.map(subprocess_run_multi, pathexe_list)
                    if verbose == True:
                        self.console.print(result)
                        # self.console.print(result.stderr)

            # process each ensemble folder one by one
            # ----------------------------------------------------------------
            else:

                self.console.print(
                    ":athletic_shoe: [b]nudging type: [/b]" + str(self.DAFLAG))

                for ens_i in track(range(self.NENS), description="Processing..."):

                    self.console.print(":keycap_number_sign: [b]ensemble nb:[/b]" + str(ens_i+1) +
                                        '/' + str(self.NENS))
                    os.chdir(os.path.join(self.workdir,
                                          self.project_name,
                                          'DA_Ensemble/cathy_' + str(ens_i+1)))
                    p = subprocess.run([callexe], text=True,
                                       capture_output=True)

                    if verbose == True:
                        self.console.print(p.stdout)
                        self.console.print(p.stderr)

            os.chdir(os.path.join(self.workdir))

            # check scenario
            # ----------------------------------------------------------------
            # and then create and applied the filter
            self._check_before_analysis(update_key)

            # analysis step
            # ----------------------------------------------------------------
            # the counter is incremented here 
            [EnsembleX, 
            self.Analysis,  
            AnalysisParam, 
            Data, 
            DataCov] = self._DA_analysis(DA_type,
                                        list_update_parm,
                                        list_assimilated_obs='all')
            

            self.count_DA_cycle = self.count_DA_cycle + 1

            # update step
            # ----------------------------------------------------------------

            if self.count_DA_cycle > 0:
                update_key = 'update_nb' + str(self.count_DA_cycle)
    
            for pp in enumerate(list_update_parm[:]):
                if 'St. var.' in pp[1]:
                    pass
                self.dict_parm_pert[pp[1]][update_key] = AnalysisParam

            
        
            # check analysis quality
            # ----------------------------------------------------------------
            # df_RMSE = self._RMSE()
            self._check_after_analysis(update_key)
    
            # create dataframe _DA_var_pert_df holding the results of the DA update
            # ---------------------------------------------------------------------
            self._DA_df(EnsembleX, self.Analysis)
    
            # update input files ensemble (state variable = psi)
            # ---------------------------------------------------------------------
            # self._update_input_ensemble(NENS,ENKFT,var_per)  #var_per[0]=='ic'
            # self._update_input_ensemble(self.NENS,self.ENKFT,self.dict_parm_pert,
            #                             update_parm_list=['St. var.'])
    
            # overwrite input files ensemble (perturbated variables)
            # ---------------------------------------------------------------------
             
            if self.count_DA_cycle<len(ENKFT):
                self._update_input_ensemble(self.NENS, ENKFT, self.dict_parm_pert,
                                    update_parm_list=list_update_parm)
    
            print('end of update')
            print(str(self.count_DA_cycle) + '/' + str(range(len(ENKFT))))
    
            # return self.DA_var_pert_df
            # return Analysis, AnalysisParam
            

            

    def create_output(self, output_dirname="output"):
        '''
        Create output directories

        Parameters
        ----------
        output_dirname : str, optional
            Name of the dir. The default is "output".

        Returns
        -------
        None.

        '''

        # create project_name/output folders
        if not os.path.exists(os.path.join(self.project_name, output_dirname)):
            os.mkdir(os.path.join(self.workdir,
                     self.project_name, output_dirname))

        if not os.path.exists(os.path.join(self.project_name, "vtk")):
            os.mkdir(os.path.join(self.workdir, self.project_name, "vtk"))

        # compute the processor (make clean make) using the Makefile
        pass

    def display_time_run(self):

        return self.total


# --------------------------------------------------------------------------- #
# update input files
# --------------------------------------------------------------------------- #


    def update_cathyH(self, verbose=False, **kwargs):
        """
        Update CathyH file input based on **kwargs parameters.

        Parameters ---------- **kwargs : all variable parameters included into cathyH. For instance
        ROWMAX # maximum NROW, with NROW = number of rows in the DEM.

        Returns ------- type Write a modified cathyH file based on new parameters parsed

        """

        indent = ''
        for kk, value in kwargs.items():
            if kk == 'indent':
                indent = value

        # if verbose == True:
        self.console.print(indent + ":black_nib: [b]Update cathyH files[/b]")

        CATHYH_file = open(
            os.path.join(self.workdir, self.project_name,
                         "src", "CATHY.H"), "r"
        )
        Lines0_109 = CATHYH_file.readlines()
        CATHYH_file.close()

        # check if hapin and parm exist
        # create them if not existing
        if hasattr(self, "hapin") is False:
            self.update_prepo_inputs()
        if "NR" not in self.parm:
            self.update_parm()

        DEMRES = 1

        if len(self.cathyH) == 0:

            self.cathyH = {
                "ROWMAX": self.hapin[
                    "M"
                ],  # maximum NROW, with NROW = number of rows in the DEM
                "COLMAX": self.hapin[
                    "N"
                ],  # maximum NCOL, with NCOL = number of columns in the DEM
                # 'COARSE': ,
                "MAXCEL": int(self.hapin["N"]) * int(self.hapin["M"]),
                "MAXRES": 1,
                "DEMRES": DEMRES,
                "NODMAX": (int(self.hapin["N"]) / DEMRES + 1)
                * (int(self.hapin["M"]) / DEMRES + 1),
                "NTRMAX": 2
                * (int(self.hapin["N"]) * int(self.hapin["M"]))
                / (DEMRES * DEMRES),
                "NP2MAX": 1,
                "MAXSTR": self.dem_parameters["nstr"],
                "NPMAX": 1,
                "NQMAX": 1,
                "NSFMAX": 1,
                "NNSFMX": 1,
                "MAXNUDN": 1,  # self.nudging["NUDN"],
                "MAXNUDT": 1,
                "MAXNUDC": 1,
                # 'MAXNENS': ,
                "MAXZON": self.dem_parameters[
                    "nzone"
                ],  # maximum NZONE, with NZONE = number of material types in the porous medium
                "MAXTRM": 52111,
                "MAXIT": 30,
                "NRMAX": self.parm["NR"],
                # maximum NPRT (ref. parm file), with NPRT = number of time values for detailed output
                "MAXPRT": self.parm["NPRT"],
                # maximum NUMVP (ref. parm input file), NUMVP = number of surface nodes for vertical profile output
                "MAXVP": int(self.parm["NUMVP"]),
                "N1MAX": self.dem_parameters[
                    "n1"
                ],  # maximum N1 (it is good to have N1 ≤ 20), N1 = maximum number of element connections to a node
                "MAXBOT": 1,
                "INTBOT": 1,
                "NTAB": 100,
                "MAXQOUT": 1,
                # Is related to data assimilation (DA) (to do not considered if you do not have DA)
                "MAXENNOUT": 52560,
                "MAXVEG": 1,
            }
        # self.run_processor(verbose=True,IPRT1=3,TRAFLAG=0)
        # self.read_grid3d()
        # # print('nnod' + str(self.nnod))
        # # print('NODMAX' + str(self.cathyH['NODMAX']))
        # print('N1MAX' + str(self.cathyH['N1MAX']))
        # print('MAXVP' + str(self.cathyH['MAXVP']))
        # sys.exit()

        # create dictionnary from kwargs
        for kk, value in kwargs.items():
            # if verbose == True:
                # self.console.print("modified:" + {keykwargs} + " | value:" + {value})
                # self.console.print(f"modified: {keykwargs} | value: {value}")
                # print(f"modified: {keykwargs} | value: {value}")
            if kk in self.cathyH.keys():
                self.cathyH[kk] = value
                self.console.print(f"modified: {kk} | value: {value}")

        # ---------------------------------------------------------------------
        # write new cathy H
        with open(
            os.path.join(self.workdir, self.project_name,
                         "src", "CATHY.H"), "w+"
        ) as CATHYH_file:
            for i, l in enumerate(Lines0_109):
                if i < 109:
                    CATHYH_file.write(l)

            CATHYH_file.write(
                "      PARAMETER (ROWMAX={},COLMAX={},DEMRES={})\n".format(
                    self.cathyH["ROWMAX"], self.cathyH["COLMAX"], self.cathyH["DEMRES"]
                )
            )
            CATHYH_file.write(
                "      PARAMETER (MAXCEL=ROWMAX*COLMAX,MAXRES=1)\n")
            CATHYH_file.write(
                "      PARAMETER (NODMAX=(ROWMAX/DEMRES+1)*(COLMAX/DEMRES+1))\n"
            )
            CATHYH_file.write(
                "      PARAMETER (NTRMAX=2*MAXCEL/(DEMRES*DEMRES))\n".format()
            )
            CATHYH_file.write(
                "      PARAMETER (NP2MAX=1,MAXSTR={})\n".format(
                    self.cathyH["MAXSTR"])
            )
            CATHYH_file.write("      PARAMETER (NFACEMAX=74000)\n".format())
            CATHYH_file.write(
                "      PARAMETER (NMAX=NODMAX*(MAXSTR + 1),NTEMAX=3*NTRMAX*MAXSTR)\n".format()
            )
            CATHYH_file.write(
                "      PARAMETER (NPMAX={},NPMAX_TRA=1,NQMAX={},NSFMAX={})\n".format(
                    self.cathyH["NPMAX"], self.cathyH["NQMAX"], self.cathyH["NSFMAX"]
                )
            )
            CATHYH_file.write(
                "      PARAMETER (NNSFMX={},MAXDIR=NODMAX+NPMAX+NSFMAX*NNSFMX)\n".format(
                    self.cathyH["NNSFMX"]
                )
            )
            CATHYH_file.write(
                "      PARAMETER (MAXNUDN={},MAXNUDT={},MAXNUDC={})\n".format(
                    self.cathyH["MAXNUDN"],
                    self.cathyH["MAXNUDT"],
                    self.cathyH["MAXNUDC"],
                )
            )
            CATHYH_file.write(
                "      PARAMETER (MAXZON={},MAXTRM={},MAXIT={},MAXVEG={})\n".format(
                    self.cathyH["MAXZON"],
                    self.cathyH["MAXTRM"],
                    self.cathyH["MAXIT"],
                    self.cathyH["MAXVEG"],
                )
            )
            CATHYH_file.write(
                "      PARAMETER (NRMAX={},MAXPRT={},MAXVP={})\n".format(
                    self.cathyH["NRMAX"], self.cathyH["MAXPRT"], self.cathyH["MAXVP"]
                )
            )
            CATHYH_file.write(
                "      PARAMETER (N1MAX={},NTPMAX=N1MAX*NMAX)\n".format(
                    self.cathyH["N1MAX"]
                )
            )
            CATHYH_file.write(
                "      PARAMETER (MAXBOT={},INTBOT={},MAXQOUT={})\n".format(
                    self.cathyH["MAXBOT"], self.cathyH["INTBOT"], self.cathyH["MAXQOUT"]
                )
            )
            CATHYH_file.write(
                "cxcx  PARAMETER (NIAUXMAX=NFACEMAX + MAXTRM + 1)\n".format()
            )
            CATHYH_file.write(
                "cxcx  PARAMETER (NRAUXMAX=5*NFACEMAX + MAXTRM,NQMAX_TRA=NODMAX)\n".format()
            )
            CATHYH_file.write(
                "      PARAMETER (NIAUXMAX=NMAX + MAXTRM + 1)\n".format())
            CATHYH_file.write(
                "      PARAMETER (NRAUXMAX=5*NMAX + MAXTRM,NQMAX_TRA=NODMAX)\n".format()
            )
            CATHYH_file.write("      PARAMETER (MAXVTKPRT=9)\n".format())
            CATHYH_file.write("      PARAMETER (MAXFCONTONODE=100,MAXLKP=3)\n")

        CATHYH_file.close()
        pass

    def update_prepo_inputs(self, DEM=None, verbose=False, show=False, **kwargs):
        """
        Update default prepro inputs i.e. hap.in and dtm_13.val files based on kwargs

        Parameters ---------- DEM : type Description of parameter `DEM`. verbose : type Description
        of parameter `verbose`. **kwargs : type Description of parameter `**kwargs`.

        Returns ------- type Description of returned object.

        """

        # print("update hap.in")
        self.console.print(":black_nib: [b]Update hap.in file[/b]")

        structural_parameter = [
            "delta_x",
            "delta_y",
            "N",
            "M",
            "N_celle",
            "xllcorner",
            "yllcorner",
        ]

        terrain_parameter = [
            "pt",
            "imethod",
            "lambda",
            "CC_threshold",
            "ndcf",
            "nchc",
            "A_threshold",
            "ASk_threshold",
            "kas",
            "DN_threshold",
            "local_slope_t",
            "p_outflow_vo",
            "bcc",
            "cqm",
            "cqg",
        ]

        # 'RIVULET NETWORK PARAMETERS (HYDRAULIC GEOMETRY OF THE SINGLE RIVULET)'
        rivulet_parameter = [
            "dr",
            "As_rf",
            "(Qsf_rf,w_rf)",
            "(Wsf_rf,b1_rf,b2_rf)",
            "(kSsf_rf,y1_rf,y2_rf)",
            "Qsi_rf",
        ]

        channel_parameter = [
            "As_cf",
            "(Qsf_cf,w_cf)",
            "(Wsf_cf,b1_cf,b2_cf)",
            "(kSsf_cf,y1_cf,y2_cf)",
            "Qsi_cf",
        ]

        # idhapin = np.ones(22) + [1,1,2,3,3,1] #+ [1,2,3,3,1]
        hapin = (
            structural_parameter
            + terrain_parameter
            + rivulet_parameter
            + channel_parameter
        )

        hap_file = open(
            os.path.join(self.workdir, self.project_name, "prepro/hap.in"), "r"
        )
        Lines = hap_file.readlines()
        hap_file.close()

        # read current values in hap.in

        tmp_param_value = []
        tmp_lnb = []  # save line nb where values are existing
        count = 0
        for line in Lines:
            x = line.strip().split("=")
            if len(x) > 1:  # skip lines with no =
                xs = " ".join(x[1].split())
                xsr = xs.replace(" ", ",")
                l = xsr.split(",")
                if isinstance(l, list):
                    l = "/t".join(l)
                tmp_param_value.append(l)
                tmp_lnb.append(count)
                count += 1

        self.hapin = {}

        for i in range(len(tmp_param_value)):
            self.hapin[hapin[i]] = tmp_param_value[i]

        # print('--'*30)
        # print(self.hapin)

        if self.hapin["dr"] != self.hapin["delta_x"]:
            print("adapt rivulet param to dem resolution")
            self.hapin["dr"] = self.hapin["delta_x"]

        # if (mod(delta_x,dr).gt.epsilon(delta_x/dr)) then:

        # check = self.hapin['delta_x'] % self.hapin['dr']
        # print(check)

        # check2 = self.hapin['delta_x']/self.hapin['dr']
        # print(check2)

        # # print('--'*30)
        # # print(self.hapin)

        # # iterate over the lines
        # # check from kwargs if values are changed
        # # save new line with new value if needed

        L = []
        tmp_param_value_new = []
        for key, value in kwargs.items():
            self.hapin[key] = value

        # # print('--'*30)
        # print(self.hapin)

        count = 0
        for i, line in enumerate(Lines):
            xnew = line.strip().split("=")
            if len(xnew) > 1:  # skip lines with no =
                for key, value in self.hapin.items():
                    if count < len(hapin):  # only consider structural parameters
                        # print(count)
                        if key == hapin[tmp_lnb[count]]:
                            # print(f'key: {key} | value: {value}')
                            if value != tmp_param_value[count]:
                                if isinstance(value, list):
                                    value_str = "/t".join(value)
                                    xnew[1] = value_str
                                else:
                                    xnew[1] = value
                                line = xnew[0] + "=              " + \
                                    str(xnew[1]) + "\n"
                                # print(line)
                        tmp_param_value_new.append(value)

                count += 1  # count line nb
            L.append(line)

        # # write the new hap.in file

        hap_file = open(
            os.path.join(self.workdir, self.project_name,
                         "prepro/hap.in"), "w+"
        )
        hap_file.writelines(L)
        hap_file.close()

        # print(self.hapin)

        # self.hapin = {}

        # for i in range(len(structural_parameter)):
        #     key = structural_parameter[i]
        #     self.hapin[key] = tmp_param_value_new[i]

        # %% dtm_13.val
        # If we start with a DEM file ("dtm_13.val") for an already delineated
        # catchment (i.e., a "catchment" DEM file instead of a "full" DEM file), then
        # only the last step in the above procedure (i.e., run "cppp" just once) is
        # needed (make sure that "Boundary channel construction" is set to 0 in
        # "hap.in").

        if DEM is not None:
            self.DEM = DEM

            # check consistency with hapin

            # check presence of the outlet
            if len(np.unique(DEM)) == 1:
                print("Error: outlet not defined")
                DEM_withOutlet = DEM
                DEM_withOutlet[0, 0] = 0
                self.DEM = DEM_withOutlet

            print("update dtm_13.val")
            with open(os.path.join(self.project_name, "prepro/dtm_13.val"), "w+") as f:
                # use exponential notation
                np.savetxt(f, self.DEM, fmt="%1.4e")

        self.update_dem_parameters(**kwargs)
        # self.update_transp(**kwargs)

        if show == True:
            # print("pltCT.dem_plot(self.workdir, self.project_name)")
            plt_CT.dem_plot(self.workdir, self.project_name)

        pass

    def update_dem_parameters(self, **kwargs):
        """
        Short summary.

        # of material types in the porous
        Parameters ---------- **kwargs : type NZONE : type
        medium. NSTR : type The number of vertical layers. N1 : type The maximum number of element
        connections to a node. NNOD : type Description of parameter `NNOD`. NTRI : type Description
        of parameter `NTRI`. ZRATIO : type The thickness of vertical layers or the fraction of
        total grid height that each layer is to occupy (ZRATIO (1) is for the surface‐most layer.
        ZRATIO values must sum to 1.). Z1 : type Description of parameter `Z1`. IVERT : type =0
        each layer will be parallel to the surface, including the base of the 3‐d grid. `ZRATIO` is
        applied to each vertical cross section. =1 base of the 3‐d grid will be flat, and `ZRATIO`
        is applied to each vertical cross section =2 base of the 3‐d grid will be flat, as will the
        NSTR‐1 horizontal cross sections above it. `ZRATIO` is applied only to the vertical cross
        section having the lowest elevation. =3 for each cell of the dem a single depth value is
        read in file input IIN60 (basement). `ZRATIO` is applied to each vertical cross section. =4
        the first NSTR‐1 layers from the surface will be parallel to the surface and the base of
        the 3‐d grid will be flat. `ZRATIO` is applied only to the vertical cross section having
        the lowest elevation. ISP : type =0 for flat surface layer (only one Z value is read in,
        and is replicated to all surface nodes); otherwise surface layer is not flat (Z values read
        in for each surface node); (for ISP=0, IVERT=0, 1, and 2 yield the same 3‐d mesh, given the
        same values of BASE and ZRATIO). BASE : type Value which defines the thickness or base of
        the 3‐d mesh. For `IVERT`=0, BASE is subtracted from each surface elevation value, so that
        each vertical cross section will be of thickness BASE, and the base of the 3‐d mesh will be
        parallel to the surface. For IVERT=1 or 2, BASE is subtracted from the lowest surface
        elevation value, say ZMIN, so that each vertical cross section will be of thickness (Z ‐
        ZMIN) + BASE, where Z is the surface elevation for that cross section. The base of the 3‐d
        mesh will thus be flat.

        Returns ------- type Description of returned object.

        """

        self.console.print(":black_nib: [b]update dem_parameters file [/b]")

        # set default parameters
        if hasattr(self, "dem_parameters") == False:
            self.console.print(
                ":pensive_face: [b]cannot find existing dem paramters[/b]")

            ltmp = [
                0.002,
                0.004,
                0.006,
                0.008,
                0.01,
                0.01,
                0.02,
                0.02,
                0.05,
                0.05,
                0.1,
                0.1,
                0.2,
                0.2,
                0.22,
            ]

            self.dem_parameters = {
                "delta_x": self.hapin["delta_x"],
                "delta_y": self.hapin["delta_y"],
                "factor": 1.0e0,
                "dostep": 1,
                "nzone": 1,
                "nstr": 15,
                "n1": 25,
                "ivert": 0,
                "isp": 1,
                "base": 3.0,
                "zratio(i),i=1,nstr": ltmp,
            }

        # create dictionnary from kwargs
        for keykwargs, value in kwargs.items():
            print(keykwargs, value)
            if keykwargs == "zratio":
                key = "zratio(i),i=1,nstr"
                if sum(value) != 1:
                    print("sum is not equal to 1 -->" + str(sum(value)))
            else:
                key = keykwargs

            try:
                self.dem_parameters[key]
                self.dem_parameters[key] = value
            except:
                pass

        for key, value in self.dem_parameters.items():
            if isinstance(value, list):
                strlst = "\t".join(str(e) for e in value)
                print(strlst)
                self.dem_parameters[key] = strlst

        # write file
        header_fmt = [1, 1, 1, 1, 3, 3, 1]
        with open(
            os.path.join(
                self.workdir, self.project_name, self.input_dirname, "dem_parameters"
            ),
            "w+",
        ) as dem_parametersfile:

            counth = 0
            for h in header_fmt:
                if h == 3:
                    dem_parametersfile.write(
                        str(list(self.dem_parameters.values())[counth])
                        + "\t"
                        + str(list(self.dem_parameters.values())[counth + 1])
                        + "\t"
                        + str(list(self.dem_parameters.values())[counth + 2])
                        + "\t"
                        + "\n"
                    )
                    counth += 3
                if h == 1:
                    dem_parametersfile.write(
                        str(list(self.dem_parameters.values())
                            [counth]) + "\t" + "\n"
                    )
                    counth += 1

            counth = 0
            for h in header_fmt:
                if h == 3:
                    dem_parametersfile.write(
                        str(list(self.dem_parameters.keys())[counth])
                        + "\t"
                        + str(list(self.dem_parameters.keys())[counth + 1])
                        + "\t"
                        + str(list(self.dem_parameters.keys())[counth + 2])
                        + "\t"
                        + "\n"
                    )
                    counth += 3
                if h == 1:
                    dem_parametersfile.write(
                        str(list(self.dem_parameters.keys())
                            [counth]) + "\t" + "\n"
                    )
                    counth += 1

        dem_parametersfile.close()

        pass

    # %% Ouput/INPUT FILES

    def update_zone(self, zone_xyz=[]):
        """
        Short summary.

        Parameters ---------- input_dirname : type Description of parameter `input_dirname`.

        Returns ------- type Description of returned object.

        """
        with open(
            os.path.join(self.workdir, self.project_name, "prepro/zone"), "w+"
        ) as zonefile:

            zonefile.write("north:     0" + "\n")
            zonefile.write("south:     0" + "\n")
            zonefile.write("east:     0" + "\n")
            zonefile.write("west:     0" + "\n")
            zonefile.write("rows:     " + str(self.hapin["M"]) + "\n")
            zonefile.write("cols:     " + str(self.hapin["N"]) + "\n")
            if len(zone_xyz) == 0:
                zone_xyz = np.c_[np.ones([self.hapin["M"], self.hapin["N"]])]
                np.savetxt(zonefile, zone_xyz, fmt="%i")
            else:
                # if np.shape(zone_xyz)== :
                np.savetxt(zonefile, zone_xyz, fmt="%i")

        zonefile.close()

        # update number of zone in the dem parameter file
        self.update_dem_parameters(nzone=len(np.unique(zone_xyz)))
        self.update_parm()
        self.update_cathyH(MAXZON=len(np.unique(zone_xyz)))

    # %% INPUT FILES

    def update_parm(self, verbose=False, **kwargs):

        self.console.print(":black_nib: [b]update parm file [/b]")

        # set default parameters
        # --------------------------------------------------------------------
        self.parm = {
            "IPRT1": 2,  # Flag for output of input and coordinate data
            "NCOUT": 0,
            "TRAFLAG": 1,
            "ISIMGR": 2,  # Flag for type of simulation and type of surface grid
            "PONDH_MIN": 0.00,  # Minimum ponding head
            "VELREC": 0,
            "KSLOPE": 0,
            "TOLKSL": 0.01,
            "PKRL": -3.0,
            "PKRR": -1.0,
            "PSEL": -3.0,
            "PSER": -1.0,
            "PDSE1L": -3.0,
            "PDSE1R": -2.5,
            "PDSE2L": -1.5,
            "PDSE2R": -1.0,
            "ISFONE": 0,
            "ISFCVG": 0,
            "DUPUIT": 0,
            "TETAF": 1.0,
            "LUMP": 1,
            "IOPT": 1,
            "NLRELX": 0,
            "OMEGA": 0.8,
            "L2NORM": 0,
            "TOLUNS": 1.0e-4,
            "TOLSWI": 1.0e30,
            "ERNLMX": 1.0e30,
            "ITUNS": 10,
            "ITUNS1": 5,
            "ITUNS2": 7,
            "ISOLV": 2,
            "ITMXCG": 500,
            "TOLCG": 1.0e-10,
            "DELTAT": 0.01,
            "DTMIN": 0.00001,  # Minimum FLOW3D time step size allowed
            "DTMAX": 10.0,  # Maximum FLOW3D time step size allowed
            # Time at end of simulation (TMAX is set to 0.0 for steady state problem)
            "TMAX": 3600.0,
            "DTMAGA": 0.0,
            "DTMAGM": 1.1,
            "DTREDS": 0.0,
            "DTREDM": 0.5,
            "IPRT": 4,
            "VTKF": 1,
            "NPRT": 3,
            "(TIMPRT(I),I=1,NPRT)": [1800.0, 3600.0, 7200.0],
            "NUMVP": 1,
            "(NODVP(I),I=1,NUMVP)": [441],
            "NR": 0,
            "NUM_QOUT": 0,
            "(ID_QOUT(I),I=1,NUM_QOUT)": [441],
        }

        # add DAFLAG ??
        # --------------------------------------------------------------------
        # Flag for the choice of the data assimilation scheme:
        # = 0 nudging (if NUDN=0, no data assimilation)
        # = 1 EnKF with Evensen's algorithm (Ocean Dynamics, 2004)
        # = 2 EnKF with parameters update
        # = 3 Particle filtering (SIR algorithm)
        # = 4 Particle filtering (SIR algorithm) with parameters update

        # if hasattr(self.parm, 'DAFLAG') is False:
        #     # temporary create DAFLAG here
        #     self.parm = {'DAFLAG': 0}  # dict of parm input parameters

        # create self.parm dictionnary from kwargs
        # --------------------------------------------------------------------
        for kk, value in kwargs.items():
            if kk in self.parm.keys():
                if verbose == True:
                    print(f"key: {kk} | value: {value}")

            # times of interest TIMPRTi
            # ----------------------------------------------------------------
            if kk == "TIMPRTi":
                key = "(TIMPRT(I),I=1,NPRT)"

                self.parm[key] = value

                # check if consistency between times of interest and
                # number of times of interest
                # ------------------------------------------------------------
                if len(value) != self.parm["NPRT"]:
                    self.parm["NPRT"] = len(value)                    
                    
                # if len(value)>1:
                #     value =  ' '.join(map(str, value))


            # points of interest NODVP
            # ----------------------------------------------------------------
            elif kk == "NODVP":
                key = "(NODVP(I),I=1,NUMVP)"
                self.parm[key] = value

                # check if consistency between node of interest and
                # number of nodes of interest
                # ------------------------------------------------------------
                if len(value) != self.parm["NUMVP"]:
                    self.parm["NUMVP"] = len(value)

            # other type of kwargs
            # ----------------------------------------------------------------
            else:
                if kk in self.parm.keys():
                    
                    # times of interest TIMPRTi
                    # ----------------------------------------------------------------
                    if kk == "TMAX":
                        # check if consistency between times of interest and
                        # TMAX
                        # ------------------------------------------------------------
                        if value < max(self.parm['(TIMPRT(I),I=1,NPRT)']):
                            self.parm["TMAX"] = max(self.parm['(TIMPRT(I),I=1,NPRT)']) 
                            
                            
                # print(self.parm)
                    # print(kk)
                # print(hasattr(self.parm,kk))
                # if hasattr(self.parm,kk):
                    self.parm[kk] = value
            # print('tmax=' + str(self.parm["TMAX"]))
        # transform array args to list
        # --------------------------------------------------------------------
        for kk, value in self.parm.items():
            if 'numpy.array' in str(type(value)):
                value = list(value)
                self.parm[kk] = value

        # update CATHYH
        # --------------------------------------------------------------------
        self.update_cathyH(MAXPRT=self.parm["NPRT"],
                           MAXVP=self.parm["NUMVP"],
                           indent='           :arrow_right_hook:')

        # write parm file
        # --------------------------------------------------------------------
        if 'filename' in kwargs:
            file2write = kwargs['filename']
        else:
            file2write = os.path.join(
                self.workdir, self.project_name, self.input_dirname, "parm")

        # backup file during DA scheme cycle
        # --------------------------------------------------------------------
        backup = False
        if 'backup' in kwargs:
            backup = kwargs['backup']
            
        if backup == True:
           dst_dir = file2write + str(self.count_DA_cycle)
           shutil.copy(file2write,dst_dir) 
           
           
        self._write_parm_file(file2write)

        pass

    def _write_parm_file(self, file2write):
        '''
        Overwrite existing parm file

        Returns
        -------
        New overwritten file.

        '''

        # header_fmt_parm = [3, 3, 2, 4, 4, 3, 3, 2, 4, 3, 3, 4, 4, 3, 2, 1, 2]
        header_fmt_parm = [3, 3, 2, 4, 4, 3, 3, 2, 4, 3, 3, 4, 4, 4, 2, 1, 2]
        counth = 0

        # with open(
        #     os.path.join(self.workdir, self.project_name, self.input_dirname, "parm"),
        #     "w+",
        # ) as parmfile:
        with open(file2write, "w+") as parmfile:

            # Write line by line according to header format
            # ----------------------------------------------------------------
            for i, h in enumerate(header_fmt_parm):
                left = right = []

                # left = values
                # ------------------------------------------------------------
                left = str(list(self.parm.values())[counth: counth + h])
                left = left.strip("[]").replace(",", "")
                left = left.strip("[]").replace("[", "")

                # right = keys
                # ------------------------------------------------------------
                right = str(list(self.parm.keys())[counth: counth + h])
                right = right.strip("[]").replace(",", "")
                right = right.replace("'", "")

                # add left + right
                # ------------------------------------------------------------
                line = left + "\t" + right + "\n"
                counth += h
                parmfile.write(str(line))

        parmfile.close()

        pass

    def update_ic(self, INDP=2, IPOND=0, WTPOSITION=[], verbose=False, **kwargs):
        '''
        The initial conditions file contains the pressure heads distribution for the study area (INDP)
        For example, to simulate a uniform water table depth or 0.5 m or 1.0 m from the ground surface,
        INDP=3 and WTHEIGHT=4.5 are selected


        Parameters
        ----------
        INDP : int, optional
            Flag for pressure head initial conditions (all nodes). The default is 2.

            - =0 for input of uniform initial conditions (one value read in)
            - =1 for input of non-uniform IC's (one value read in for each node)
            - =2 for calculation of fully saturated vertical hydrostatic equilibrium IC's
              (calculated in subroutine ICVHE). In the case of IPOND>0, the fully saturated
              hydrostatic IC is calculated (in subroutine ICVHEPOND) starting from the ponding head
              values at the surface nodes, rather than surface pressure heads of 0.
            - =3 for calculation of partially saturated vertical hydrostatic equilibrium IC's
              (calculated in subroutine ICVHWT) with the water table height (relative to the base
              of the 3‐d grid) given by parameter WTHEIGHT
            - =4 for calculation of partially saturated vertical hydrostatic equilibrium IC's
              (calculated in subroutine ICVDWT) with the water table depth (relative to the surface
              of the 3‐d grid) given by parameter WTHEIGHT

        IPOND : int, optional
            Flag for ponding head initial conditions (surface nodes). The default is 0.

            - =0 no input of ponding head initial conditions; otherwise (IPOND = 1 or 2) ponding
              head initial conditions are read into PONDNOD, and, where PONDNOD > 0, these values
              are used to update the surface node values in PTIMEP read in according to the
              previous INDP flag
            - =1 uniform ponding head initial conditions (one value read in)
            - =2 non-uniform ponding head initial conditions (one value read in for each node)

        WTPOSITION : float, optional
            For the case INDP=3, specifies the initial water table height relative to
            the base of the 3‐d grid. The default is [].

        Returns
        -------
        new ic file written/overwritten.

        '''

        #  check value of WTPOSITION
        # --------------------------------------------------------------------
        # For the case INDP=3, specifies the initial water table
        # height relative to the base of the 3‐d grid
        # if WTPOSITION>0:
        #     print('WTPOSITION must be negative - water table height
                # relative to the base of the 3‐d grid')
        #     sys.exit()

        # set default parameters
        # --------------------------------------------------------------------
        self.ic = {
            "INDP": INDP,
            "IPOND": IPOND,
            # For the case INDP=3, specifies the initial water table
            # height relative to the base of the 3‐d grid
            "WTPOSITION": WTPOSITION
                }

        # write ic file
        # --------------------------------------------------------------------

        if 'filename' in kwargs:
            filename = kwargs['filename']
        else:
            filename = os.path.join(
                self.workdir, self.project_name, self.input_dirname, "ic")
            
        backup = False
        if 'backup' in kwargs:
            backup = kwargs['backup']
            
        if backup == True:
           dst_dir = filename + str(self.count_DA_cycle-1)
           shutil.copy(filename,dst_dir) 
           
            
        if INDP == 0:
            with open(filename, "w+") as icfile:
                icfile.write(str(INDP) + "\t" + str(IPOND) +
                             "\t" + "INDP" + "\t" + "IPOND" + "\n")
                icfile.write(str(WTPOSITION) + "\t" + "WTPOSITION" + "\n")

            icfile.close()
        elif INDP == 1:

            pressure_head_ini = kwargs['pressure_head_ini']
            # print(pressure_head_ini)

            with open(filename, "w+") as icfile:
                icfile.write(str(INDP) + "\t" + str(IPOND) +
                             "\t" + "INDP" + "\t" + "IPOND" + "\n")
                # icfile.write(pressure_head_ini)
                np.savetxt(icfile, pressure_head_ini, fmt="%1.3e")

            icfile.close()
            

           

        pass

    def update_atmbc(self, HSPATM=0, IETO=0, TIME=None, VALUE=[None, None],
                     show=False, verbose=False, **kwargs):
        '''
        Atmospheric forcing term (atmbc - IIN6)


        ..note


        Parameters
        ----------
        HSPATM : int, optional

            - =0 for spatially variable atmospheric boundary condition inputs;
            blank or =9999 if unit IIN6 input is to be ignored; otherwise atmospheric BC's are
            homogeneous in space.
        IETO : TYPE, optional
            - =0 for linear interpolation of the atmospheric boundary condition inputs between different
            - otherwise the inputs are assigned as a piecewise constant function (ietograph).
            The default is 0.
        TIME : TYPE, optional
            DESCRIPTION. The default is None.
        VALUE : TYPE, optional
            DESCRIPTION. The default is [None, None].
        show : TYPE, optional
            DESCRIPTION. The default is False.
        verbose : TYPE, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        if 'diff' in kwargs:
            v_atmbc = VALUE[0] - VALUE[1]
        else:
            v_atmbc = VALUE

        # set default parameters
        # --------------------------------------------------------------------

        self.atmbc = {"HSPATM": HSPATM, "IETO": IETO,
            "TIME": TIME, "VALUE": VALUE}

        # overwrite existing input atmbc file
        # --------------------------------------------------------------------

        if 'filename' in kwargs:
            filename = kwargs['filename']
        else:
            filename = os.path.join(
                self.workdir, self.project_name, self.input_dirname, "atmbc")

        # print('????????????????????????????/')
        # print(filename)
        # print(TIME)
        # print(v_atmbc)
        backup = False
        if 'backup' in kwargs:
            backup = kwargs['backup']
            
        if backup == True:
           dst_dir = filename + str(self.count_DA_cycle-1)
           shutil.copy(filename,dst_dir) 
           
                   
        with open(filename, "w+") as atmbcfile:
            atmbcfile.write(str(HSPATM) + "\t" + str(IETO) +
                            "\t" + "HSPATM" + "\t" + "IETO" + "\n"
                            )

            # if v_atmbc is a scalar meaning that atmbc are homoegenous
            # -----------------------------------------------------------------
            for t, v in zip(TIME, v_atmbc):

                if verbose == True:
                    print(t, v)
                atmbcfile.write(str(t) + "\t" + "TIME" + "\n")
                # atmbcfile.close()

                if isinstance(v, float):
                    atmbcfile.write(str(v) + "\t" + "VALUE" + "\n")
                else:
                    print(t, v)
                    # atmbcfile.write(str(v) + "\t" + "VALUE" + "\n")  
                    np.savetxt(atmbcfile, v, fmt="%1.3e")

                    
        atmbcfile.close()

        self.update_parm(NPRT=len(TIME))

        # don't need to update if sequential DA as the cathy.exe is already created
        # ---------------------------------------------------------------------
        if 'omit_cathyH' not in kwargs:
            self.update_cathyH(MAXPRT=len(TIME))

        if show == True:
            # if HSPATM !=0:
            #     print('impossible to plot for non homogeneous atmbc')
            #     # sys.exit()
            # else:
            x_units = "sec"
            for key, value in kwargs.items():
                if key == "x_units":
                    x_units = value
            plt_CT.show_atmbc(TIME, VALUE, x_units=x_units)

        pass

    def update_nansfdirbc(
        self,
        TIME=[],
        NDIR=0,
        NDIRC=0,
        NQ3=None,
        noflow=True,
        bound_xyz=None,
        pressureHead=[],
    ):
        """
        Boundary conditions (nansfdirbc - IIN8, nansfneubc - IIN9, sfbc - IIN7) The boundary
        conditions are defined in the nansfdirbc (Dirichlet), nansfneubc (Neumann), and sfbc
        (seepage face) files. To simulate the no-flow boundaries conditions for the bottom and
        vertical sides of the domain it is necessary to set NDIR and NDIRC equal to zero. To
        simulate different boundary conditions, it is necessary to indicate the number of selected
        nodes through NDIR or NDIRC, then to specify the node ID’s that you want to consider and
        eventually the value of pressure head or flux that you want to assign.

         !!! update_nansfdirbc use the grid3d

        Parameters ---------- NDIR : int Number of non-atmospheric, non‐seepage face Dirichlet
        nodes in 2-d mesh. The BC's assigned to these surface nodes are replicated vertically
        (compare NDIRC) NDIRC : int Number of 'fixed' non-atmospheric, non-seepage face Dirichlet
        nodes in 3‐d mesh ('fixed' in the sense that these BC's are not replicated to other nodes ‐
        compare NDIR) NQ3 : int Number of non-atmospheric, non‐seepage face Neumann nodes in 3‐d
        mesh.

        Returns ------- type Description of returned object.

        """

        # read existing input nansfdirbc file
        # --------------------------------------------------------------------

        try:
            self.grid3d = in_CT.read_grid3d(self.project_name)
        except OSError:
            print("grid3d missing - need to run the processor with IPRT1=2 first")
            # self.run_processor(verbose=True,IPRT1=2,TRAFLAG=0)
            # self.read_grid3d()
            # # sys.exit()

        with open(
            os.path.join(
                self.workdir, self.project_name, self.input_dirname, "nansfdirbc"
            ),
            "w+",
        ) as nansfdirbcfile:

            if noflow:
                if len(TIME) == 0:
                    TIME = self.atmbc["TIME"]
                for tt in TIME:
                    nansfdirbcfile.write(str(tt) + "\t" + "TIME" + "\n")
                    nansfdirbcfile.write(
                        str(NDIR)
                        + "\t"
                        + str(NDIRC)
                        + "\t"
                        + "NDIR"
                        + "\t"
                        + "NDIRC"
                        + "\n"
                    )
            # elif bound_xyz is not None:
            #     for tt in TIME:
            #         #Se avessi delle variazioni dovrei indicare il nodo ed il valore di pressione
            #         nansfdirbcfile.write(str(tt) + "\t" + 'TIME' + "\n")
            #         for i in range(self.nnod3):
            #             if self.xmesh[i] == xb_left or self.xmesh[i] == xb_right or
            #                self.ymesh[i] == yb_left or self.ymesh[i] == yb_right:
            #                    nansfdirbcfile.write(str(i) + "\n")

            #         for i in range(self.nnod3):
            #             if self.xmesh[i] == xb_left or self.xmesh[i] == xb_right or
            #                self.ymesh[i] == yb_left or self.ymesh[i] == yb_right:
            #                    nansfdirbcfile.write(str(self.zmesh[i]-self.ic['WTPOSITION']) + "\n")

            # exemple provided by Laura B.

            # C     Write dirbc
            #       write(33,*) 0.0, 'TIME'
            #       write(33,*) '0', a
            #       do i=1,nnod3
            #          if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
            #      1       (y(i).eq.5))then
            #          write(33,*) i
            #          endif
            #       enddo
            #       do i=1,nnod3
            #          if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
            #      1       (y(i).eq.5))then
            #          write(33,*) -z(i)-WTdepth
            #          endif
            #       enddo

            #       write(33,*) 2e+20, 'TIME'
            #       write(33,*) '0', a
            #       do i=1,nnod3
            #          if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
            #      1       (y(i).eq.5))then
            #          write(33,*) i
            #          endif
            #       enddo
            #       do i=1,nnod3
            #          if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
            #      1       (y(i).eq.5))then
            #          write(33,*) -z(i)-WTdepth
            #          endif
            #       enddo

        # modicare il valore di NPMAX nel file 27 CATHY.H nel caso
        # in cui si inseriscano dei NDIRC ed il valore di NP2MAX nel caso si inseriscano dei
        # NDIR. I valori di NPMAX e NP2MAX corrispondono al numero massimo
        # di nodi NDIRC e NDIR che si possono inserire.
        nansfdirbcfile.close()

        self.update_parm()
        self.update_cathyH()

        # set default parameters
        # self.nansfdirbc =	{
        #               "NDIR": HSPATM,
        #               "NDIRC": IETO,
        #              }

        pass

    def update_nansfneubc(self, TIME=[], NQ=0, ZERO=0):
        '''


        Parameters
        ----------
        TIME : TYPE, optional
            DESCRIPTION. The default is [].
        NQ : TYPE, optional
            DESCRIPTION. The default is 0.
        ZERO : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        None.

        '''
        # read existing input nansfneubc file
        # --------------------------------------------------------------------

        with open(
            os.path.join(
                self.workdir, self.project_name, self.input_dirname, "nansfneubc"
            ),
            "w+",
        ) as nansfneubcfile:

            if len(TIME) == 0:
                TIME = self.atmbc["TIME"]
            for tt in TIME:
                nansfneubcfile.write(str(tt) + "\t" + "TIME" + "\n")
                nansfneubcfile.write(
                    str(ZERO) + "\t" + str(NQ) + "\t" +
                        "ZERO" + "\t" + "NQ" + "\n"
                )

        nansfneubcfile.close()

        # set default parameters
        # --------------------------------------------------------------------

        # set default parameters
        # self.nansfneubc =	{
        #               "NDIR": HSPATM,
        #               "NDIRC": IETO,
        #              }

        pass

    def update_sfbc(self, TIME=[]):

        with open(
            os.path.join(self.workdir, self.project_name,
                         self.input_dirname, "sfbc"),
            "w+",
        ) as sfbcfile:

            if len(TIME) == 0:
                TIME = self.atmbc["TIME"]
            for tt in TIME:
                sfbcfile.write(str(tt) + "\n")
                sfbcfile.write("0" + "\n")

        sfbcfile.close()

        # set default parameters
        # self.nansfneubc =	{
        #               "NDIR": HSPATM,
        #               "NDIRC": IETO,
        #              }

        pass

    def update_soil(self, IVGHU=[], FP=[], SPP=[], verbose=False, **kwargs):
        '''
        Soil parameters (soil - IIN4). The porous media properties.

        ..note:
            The first thing that must be decides is the type of relationship to describe the hydraulic
            characteristics of the unsaturated soil (i.e. retention curves). This can be done through
            the choice of the parameter **IVGHU** amongst the several options.


        Parameters
        ----------
        IVGHU : int, optional

            = -1 table look up for moisture curves

            = 0 for van Genuchten moisture curves

            = 1 for extended van Genuchten moisture curves

            = 2 for moisture curves from Huyakorn et al
                (WRR 20(8) 1984, WRR 22(13) 1986)
                with Kr=Se**n conductivity relationship

            = 3 for moisture curves from Huyakorn et al
                (WRR 20(8) 1984, WRR 22(13) 1986)
                with conductivity relationship from Table 3 of 1984 paper (log_10 Kr(Se) curve)

            = 4 for Brooks‐Corey moisture curves.

            The default is [].

        FP : list, optional

            Feddes Parameters. The default is [].

            [PCANA PCREF PCWLT ZROOT PZ OMGC]
            - 'PCANA': float anaerobiosis point
            - 'PCREF': float field capacity
            - 'PCWLT': float wilting point
            - 'ZROOT': float root depth
            - 'PZ': float ??
            - 'OMGC': float ??

            .. note:
                For details, see http://dx.doi.org/10.1002/2015WR017139

        SPP : list, optional
            Soil Physical Properties. The default is [].
            - 'PERMX' (NSTR, NZONE): saturated hydraulic conductivity - xx
            - 'PERMY' (NSTR, NZONE): saturated hydraulic conductivity - yy
            - 'PERMZ' (NSTR, NZONE): saturated hydraulic conductivity - zz
            - 'ELSTOR' (NSTR, NZONE): specific storage
            - 'POROS'  (NSTR, NZONE): porosity (moisture content at saturation)

            retention curves parameters VGN, VGRMC, and VGPSAT
            - 'VGNCELL' (NSTR, NZONE)
            - 'VGRMCCELL' (NSTR, NZONE)
            - 'VGPSATCELL' (NSTR, NZONE):

        verbose : TYPE, optional
            DESCRIPTION. The default is False.
        **kwargs

        Returns
        -------
        None.

        '''

        # set default parameters if "IVGHU": 0 and if SPP is not existing yet
        # --------------------------------------------------------------------
        if len(SPP) == 0:
            SPP = self.set_SOIL_defaults()
        # print(SPP)

        # read function arguments kwargs and udpate soil and parm files
        # --------------------------------------------------------------------
        for keykwargs, value in kwargs.items():
            if verbose == True:
                print(f"new: {keykwargs} | value: {value}")
            self.soil[keykwargs] = value
            self.parm[keykwargs] = value

        # loop over Feddes parameters FP mandatory fct arg
        # --------------------------------------------------------------------
        for fp in FP:  # loop over fedded parameters
            if verbose == True:
                print(fp, FP[fp])
            self.soil[fp] = FP[fp]

        # create prepro inputs if not existing (containing info about the DEM)
        # --------------------------------------------------------------------
        if hasattr(self, "dem_parameters") is False:
            self.update_prepo_inputs()

        # Soil Physical Properties strat by strat
        # --------------------------------------------------------------------
        SoilPhysProp = self._prepare_SPP_tb(SPP)
        self.soil['SPP'] = SoilPhysProp

        # Vegetation properties (PCANA,PCREF,PCWLT,ZROOT,PZ,OMGC)
        # --------------------------------------------------------------------
        # Read only if IVGHU=0
        FeddesParam = self._prepare_SOIL_vegetation_tb(FP)
        # print(FeddesParam)
        self.soil['FP'] = FeddesParam

        # write soil file
        # --------------------------------------------------------------------
        self._write_SOIL_file(SoilPhysProp, FeddesParam)

        pass

    def set_SOIL_defaults(self):

        self.soil = {
            "PMIN": -5.0,
            "IPEAT": 0,
            "SCF": 1.0,
            "CBETA0": 0.4,
            "CANG": 0.225,
            # Feddes parameters default values
            "PCANA": [0.0],
            "PCREF": [-4.0],
            "PCWLT": [-150],
            "ZROOT": [1.0],
            "PZ": [1.0],
            "OMGC": [1.0],
            "IVGHU": 0,  # The first thing that must be decided!
            "HUALFA": 0.02,
            "HUBETA": 2,
            "HUGAMA": 2,
            "HUPSIA": 0,
            "HUSWR": 0.333,
            "HUN": 1,
            "HUA": -5,
            "HUB": 1,
            "BCBETA": 1.2,
            "BCRMC": 0,
            "BCPSAT": -0.345,
        }

        # set Soil Physical Properties defaults parameters
        # --------------------------------------------------------------------

        PERMX = PERMY = PERMZ = 1.88e-04
        ELSTOR = 1.00e-05
        POROS = 0.55
        VGNCELL = 1.46
        VGRMCCELL = 0.15
        VGPSATCELL = 0.03125

        SPP = {
            "PERMX": PERMX,
            "PERMY": PERMY,
            "PERMZ": PERMZ,
            "ELSTOR": ELSTOR,
            "POROS": POROS,
            "VGNCELL": VGNCELL,
            "VGRMCCELL": VGRMCCELL,
            "VGPSATCELL": VGPSATCELL,
        }

        # set Soil Feddes Properties defaults parameters
        # --------------------------------------------------------------------
        # if len(FP)==0:
        #     FP = {'PCANA':self.soil['PCANA'],
        #           'PCREF':self.soil['PCREF'],'PCWLT':self.soil['PCWLT'],
        #           'ZROOT':self.soil['ZROOT'],'PZ':self.soil['PZ'],
        #           'OMGC':self.soil['OMGC']}

        return SPP

    def _prepare_SPP_tb(self, SPP):
        '''
        _prepare_SOIL_Physical_Properties_tb

        Parameters
        ----------
        SPP : TYPE
            DESCRIPTION.

        Returns
        -------
        np.array describing the SoilPhysProp with rows corresponding to the layer.

        '''
        # check number of zones
        if self.dem_parameters["nzone"] > 1:

            # initiate a table to fill with dimension of number of zone * number of strates
            SoilPhysProp = np.ones(
                [self.dem_parameters["nzone"] * self.dem_parameters["nstr"], 8]
                )
            k = 0

            # loop over strates
            # -----------------------------------------------------------------
            for istr in range(self.dem_parameters["nstr"]):
                izoneSoil = np.zeros([self.dem_parameters["nzone"], 8])

                #  loop over zones within a strate
                # --------------------------------------------------------------
                for izone in range(self.dem_parameters["nzone"]):
                    izoneSoil_tmp = []
                    for spp in SPP:
                        izoneSoil_tmp.append(SPP[spp][izone])
                    izoneSoil_tmp = np.hstack(izoneSoil_tmp)
                    izoneSoil[izone, :] = izoneSoil_tmp
                    ki = k
                    ke = k + self.dem_parameters["nzone"]
                    SoilPhysProp[ki:ke, :] = izoneSoil
                    # SoilPhysProp[self.dem_parameters['nzone']*istr*izone:self.dem_parameters['nzone']*istr*izone+self.dem_parameters['nzone'],:]=izoneSoil
                k += self.dem_parameters["nzone"]

        # case if there is only one zone in the mesh
        else:
            izoneSoil = []
            for spp in SPP:
                izoneSoil.append(SPP[spp])
            izoneSoil = np.hstack(izoneSoil)
            SoilPhysProp = np.tile(izoneSoil, (self.dem_parameters["nstr"], 1))

        return SoilPhysProp

    def _prepare_SOIL_vegetation_tb(self, FP):
        '''
        _prepare_SOIL_vegetation_tb

        Parameters
        ----------
        FP : dict
            dict containing Feddes parameters.
            - 'PCANA': anaerobiosis point
            - 'PCREF': field capacity
            - 'PCWLT': wilting point
            - 'ZROOT': float root depth
            - 'PZ': ??
            - 'OMGC': float ??
            For details, see http://dx.doi.org/10.1002/2015WR017139
        Returns
        -------
        FeddesParam: numpy array
            table or array describing Feddes parameters for a given DEM




        '''
        # Vegetation properties (PCANA,PCREF,PCWLT,ZROOT,PZ,OMGC)
        # --------------------------------------------------------------------
        # Read only if IVGHU=0

        # check number of vegetation
        if self.cathyH["MAXVEG"] > 1:
            print(self.cathyH["MAXVEG"])
            FeddesParam = np.zeros([self.cathyH["MAXVEG"], 6])
            for iveg in range(
                self.cathyH["MAXVEG"]
            ):  # loop over veg zones within a strate
                izoneVeg_tmp = []
                for sfp in FP:
                    izoneVeg_tmp.append(FP[sfp][iveg])
                    # if iveg==1:
                    #     izoneVeg_tmp.append('\t PCANA PCREF PCWLT ZROOT PZ OMGC')

                izoneVeg_tmp = np.hstack(izoneVeg_tmp)
                FeddesParam[iveg, :] = izoneVeg_tmp

        # case where unique vegetation type
        else:
            FeddesParam = np.c_[
                self.soil["PCANA"],
                self.soil["PCREF"],
                self.soil["PCWLT"],
                self.soil["ZROOT"],
                self.soil["PZ"],
                self.soil["OMGC"],
            ]

        return FeddesParam

    def _write_SOIL_file(self, SoilPhysProp, FeddesParam, **kwargs):
        '''
        _write_SOIL_file

        Parameters
        ----------
        SoilPhysProp : TYPE
            DESCRIPTION.
        FeddesParam : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        # number of side header for each row
        header_fmt_soil = [1, 2, 2, 6, 1, 5, 1, 2, 3]

        # open soil file
        # --------------------------------------------------------------------
        if self.DAFLAG:
            soil_filepath = os.path.join(
                os.getcwd(), self.input_dirname, "soil")
        else:
            soil_filepath = os.path.join(
                self.workdir, self.project_name, self.input_dirname, "soil")

        with open(os.path.join(soil_filepath), "w+") as soilfile:

            counth = 0  # count header index

            # Write line by line according to header format
            # ----------------------------------------------------------------
            for i, h in enumerate(header_fmt_soil):

                # left = values
                # ------------------------------------------------------------
                left = right = []
                left = str(list(self.soil.values())[counth: counth + h])
                left = left.strip("[]").replace(",", "")

                # right = keys
                # ------------------------------------------------------------
                right = str(list(self.soil.keys())[counth: counth + h])
                right = right.strip("[]").replace(",", "")
                right = right.replace("'", "")

                # Line 4: write Feddes Parameters Table
                # ------------------------------------------------------------
                if i == 3:  # Feddes parameters
                    np.savetxt(soilfile, FeddesParam, fmt="%1.3e")
                    counth += h
                else:
                    line = left + "\t" + right + "\n"
                    counth += h
                    soilfile.write(str(line))

            # End of soil file: write Soil Properties Table
            # ------------------------------------------------------------
            np.savetxt(soilfile, SoilPhysProp, fmt="%1.3e")
            soilfile.write(
                "PERMX PERMY  PERMZ  ELSTOR POROS,VGNCELL,VGRMCCELL,VGPSATCELL"
                + "\n"
            )

        soilfile.close()
        # print('new SOIL file written')

    def update_veg_map(self, indice_veg=1, show=False, **kwargs):
        '''
        Contains the raster map describing which type of vegetation every cell belongs to.


        Parameters
        ----------
        indice_veg : TYPE, optional
            DESCRIPTION. The default is 1.
        show : TYPE, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        indice_veg : TYPE
            DESCRIPTION.

        '''

        if isinstance(indice_veg, int):
            indice_veg = float(indice_veg)

        if hasattr(self, "hapin") is False:
            self.update_prepo_inputs()

        print("update root map")
        with open(
            os.path.join(
                self.workdir, self.project_name, self.input_dirname, "root_map"
            ),
            "w+",
        ) as rootmapfile:

            rootmapfile.write("north:     0" + "\n")
            rootmapfile.write("south:     0" + "\n")
            rootmapfile.write("east:     0" + "\n")
            rootmapfile.write("west:     0" + "\n")
            rootmapfile.write("rows:     " + str(self.hapin["M"]) + "\n")
            rootmapfile.write("cols:     " + str(self.hapin["N"]) + "\n")

            if isinstance(indice_veg, float):
                # if  root_depth>self.dem_parameters['base']:
                #     print('max root mesh > max mesh depth')
                #     sys.exit()
                indice_veg = (
                    np.c_[
                        np.ones([int(self.hapin["M"]), int(self.hapin["N"])])]
                    * indice_veg
                )
                np.savetxt(rootmapfile, indice_veg, fmt="%1.2e")
            else:
                # if np.shape(zone_xyz)== :
                # if  max(max(root_depth))>self.dem_parameters['base']:
                # print('max root mesh > max mesh depth')
                # sys.exit()
                np.savetxt(rootmapfile, indice_veg, fmt="%1.2e")

        rootmapfile.close()

        self.update_cathyH(MAXVEG=len(np.unique(indice_veg)))

        if show is not None:
            plt_CT.indice_veg_plot(indice_veg)

        return indice_veg

    # -------------------------------------------------------------------#
    # %% DATA ASSIMILATION FCTS

    def _DA_init(self, NENS=[], ENKFT=[], parm_pert=[]):
        """
        THIS SHOULD BE MOVED TO DA CLASS

        Initial preparation for DA


        .. note::

            1. update cathyH file given NENS, ENKFT + flag self.parm['DAFLAG']==True

            2. create the processor cathy_origin.exe

            3. create directory tree for the ensemble

            4. create panda dataframe for each perturbated variable

            5. update input files

        Parameters
        ----------
        NENS : int
            # Number of realizations in the ensemble kalman filter (EnKF)
        ENKFT : pandas df
            # Observation times for ensemble kalman filter (EnKF).

        Returns
        -------
        None.

        Note: for now only implemented for EnkF DA
        """

        # to do create a readme file inside DA folder to explain how the dir tree works
        # self.DAFLAG=True

        # update nudging file
        # self.update_nudging(NUDN=NUDT)
        # self.nudging = {'NUDN': NUDN}   # THIS IS TEMPORARY
        self.NENS = NENS  # THIS IS TEMPORARY

        # update cathyH DA parameters
        # ---------------------------------------------------------------------
        # NUDN  = number of observation points for nudging or EnKF (NUDN=0 for no nudging)
        # NUDN = len(ENKFT)  # this is true only when DA is made inside CATHY
        # self.update_cathyH(MAXNUDN=NUDN,ENKFT=ENKFT, verbose=True)

        # when DA is made outside CATHY we run step by step

        # TIMPRTi # time of interest for outputs
        # TIMPRTi=times_of_interest,
        # NODVP=nodes_of_interest,

        # THIS IS TEMPORARY, assimilation time should be infer from observation data directly
        # self.ENKFT = ENKFT

        # ?? is this necessary ??
        # self.update_cathyH(MAXNUDN=1,ENKFT=ENKFT, verbose=True)

        # run processor to create the cathy_origin.exe to paste in every folder
        # ---------------------------------------------------------------------
        # runProcess=True means open loop simulation 
    
        self.run_processor(runProcess=False)
        # OR self.recompileSrc(runProcess=False, NUDN=NUDN)

        # create sub directories for each ensemble
        # ---------------------------------------------------------------------
        self._create_subfolders_ensemble(NENS)

        # create initial dataframe DA_results_df for results
        # ---------------------------------------------------------------------
        # X, N_col, M_rows = self._read_state_ensemble()      
        
        self._DA_df()

        # update input files ensemble
        # ---------------------------------------------------------------------
        self._update_input_ensemble(NENS, ENKFT, parm_pert)

        # return self.DA_var_pert_df
        pass

    def _check_before_analysis(self,update_key):
        '''
        Filter is applied only on selected ensemble

        Parameters
        ----------
        update_key : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        
        

        
    def _check_after_analysis(self,update_key):
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        CHECK which scenarios are OK and which one to discard

        Returns
        -------
        None.

        '''
        print('check scenario')

        # ckeck if new value of Zroot is feasible
        if (self.dict_parm_pert['ZROOT'][update_key] < min(self.grid3d['nodes_idxyz'][:, -1])).any:
            print('impossible value of zroot')
            self.dict_parm_pert['ZROOT'][update_key] = self.dict_parm_pert['ZROOT']['ini_perturbation']-0.01*self.count_DA_cycle

            # print(self.dict_parm_pert)

        # condition 1: nrow(mbeconv)==0
        # condition 2: mbeconv[h,3]==(deltaT)

        pass

    def _DA_analysis(self, typ='Sakov',
                              list_update_parm=['St. var.'],
                              list_assimilated_obs=[]):
        """
        THIS SHOULD BE MOVED TO DA CLASS

        Analysis ensemble using ENKF


        1. prepare matrices

        2. apply filter

        3. checkings

        3. create panda dataframe for each perturbated variable

        4. update input files

        Parameters
        ----------
        NENS : int
            # Number of realizations in the ensemble kalman filter (EnKF)
        NUDN : int
            #  NUDN  = number of observation points for nudging or EnKF (NUDN=0 for no nudging)
        ENKFT : np.array([])
            # Observation time for ensemble kalman filter (EnKF).

        Returns
        -------
        None.

        Note: for now only implemented for EnkF DA
        """

        update_key = 'ini_perturbation'
        if self.count_DA_cycle > 0:
            update_key = 'update_nb' + str(self.count_DA_cycle)

        # map states to observation = apply H operator to state variable
        # ---------------------------------------------------------------------
        Hx = self._map_states2Observations(list_assimilated_obs)
        Observation = Hx  # (EnSize)

        # prepare ENKF
        # ---------------------------------------------------------------------
        # When updating the states only, the elements of X are the
        # pressure heads at each node of the finite element grid, while
        # the state augmentation technique is used when also updat-
        # ing the parameters

        EnsembleX, EnsembleX_augm, A, A_mean = self._prepare_ENKF(list_update_parm)

        print(list_update_parm)
        for pp in enumerate(list_update_parm[:]):
            if 'St. var.' in pp[1]:
                pass
            
            # (EnSize), self.dict_parm_pert
            Param = self.dict_parm_pert[pp[1]][update_key]
            print(pp)
        # apply ENKF
        # ---------------------------------------------------------------------
        # compute data covariance
        # ---------------------------------------------------------------------
        # for oo in list_assimilated_obs:
        #     print(oo)
        
        # dict_obs_cp = self.dict_obs.copy()
        tuple_list_obs = list(self.dict_obs.items())
        key_value = tuple_list_obs[self.count_DA_cycle]


        Data = key_value[1]['data']['resist'].to_numpy()  # (MeasSize)
        DataCov = key_value[1]['data']['recipError'].to_numpy()  # (MeasSize)

        [A, Amean, dA, 
         dD, MeasAvg, S, 
         COV, B, dAS, 
         Analysis, 
         AnalysisParam] = enkf.enkf_analysis(Data, 
                                            DataCov, 
                                            Param, 
                                            EnsembleX, 
                                            Observation)


        # Analysis, AnalysisParam = enkf.enkf_analysis(Data, 
        #                                              DataCov, 
        #                                              Param, 
        #                                              EnsembleX, 
        #                                              Observation)

        # return EnsembleX, Analysis, AnalysisParam, Data, DataCov
        return [A, Amean, dA, 
         dD, MeasAvg, S, 
         COV, B, dAS, 
         Analysis, 
         AnalysisParam]






    def _map_states2Observations(self, list_assimilated_obs, parallel=False,
                                 verbose = False):
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        Apply H mapping operator to convert model to predicted value

        Returns
        -------
        Hx : np.array
            Ensemble of the simulated observations.
        Hx_mean : TYPE
            Ensemble mean of the simulated observations.

        '''
        
        print(self.count_DA_cycle)
        self.count_DA_cycle = 1
        if list_assimilated_obs == 'all':
            items_dict = list(self.dict_obs.items())
            obs2map = [items_dict[self.count_DA_cycle][1]['data_type']]

        Hx = []
        # Hx_mean = []
        for obs_key in obs2map:
            # print(obs_key)

            if parallel == True:

                pathERT = os.path.split(self.dict_obs[0]['filename'])[0]
                forward_mesh_vtk_file = glob.glob(
                    pathERT + "/**/*forward_model*", recursive=True)

                if 'ERT' in obs_key:
                    for ens_nb in range(self.NENS):
                        path_fwd_ERT_list = []
                        path_fwd_ERT = os.path.join(self.workdir, self.project_name,
                                                    './DA_Ensemble/cathy_' + str(ens_nb+1))
                        path_fwd_ERT_list.append(path_fwd_ERT)
                        # path_CATHY = os.path.join(path_fwd_ERT, 'vtk/')

                with multiprocessing.Pool(processes=4) as pool:
                    
                    # Parallel run of a function with multiple arguments
                    # set constant values to all arguments which are not changed 
                    # during parallel processing
 
                    ERTmapping_args =partial(ERTmapping_run_multi, 
                                             dict_obs = self.dict_obs,
                                             count_DA_cycle = self.count_DA_cycle,
                                             # path_CATHY = path_CATHY,
                                             porosity = self.soil['SPP'][:, 4][0], 
                                             forward_mesh_vtk_file = forward_mesh_vtk_file[0],
                                             output_dirname =  self.output_dirname,
                                             project_name = self.project_name 
                                             ) 
                    Hx_ERT = pool.map(ERTmapping_args, 
                                      path_fwd_ERT_list, 
                                      )
                    print(Hx_ERT)
                    self._add_2_ensemble_Hx(Hx, Hx_ERT['resist'])
                    print(Hx)


                    
                    if verbose == True:
                        self.console.print(Hx)
                        # self.console.print(result.stderr)
            
        
        
            else:
                # loop over ensemble files
                # -----------------------------------------------------------------
                for ens_nb in range(self.NENS):
                # for ens_nb in track(range(NENS), description="updating ensemble file..."):
    
                    # change directory according to ensmble file nb
                    # os.chdir(os.path.join(self.workdir, self.project_name,
                    #                       './DA_Ensemble/cathy_' + str(ens_nb+1)))
                    path_fwd_ERT = os.path.join(self.workdir, self.project_name,
                                          './DA_Ensemble/cathy_' + str(ens_nb+1))
                    # need to read the vp file
                    # --------------------------------------------------------------------
                    # df_vp = self.read_outputs(filename='vp')
                        # PH: pressure head
                        # SW: Soil Water
                        # CKRW: Relative hydraulic conductivity output at all nodes
                    # groupby df_vp['node'] to elect a given node information
    
                    df_sw = self.read_outputs(filename='sw',
                                              path=os.path.join(path_fwd_ERT,
                                                                self.output_dirname))
    
                    # node = self.dict_obs[key]['position']
                    # time = self.dict_obs[key]['time']
                    if 'ERT' in obs_key:
                        # case 4: ERT
                        # --------------------------------------------------------------------
                        # if key[0] in 'ERT':
                        pathERT = os.path.split(self.dict_obs[0]['filename'])[0]
                        # print('********')
                        # print(self.count_DA_cycle)
                        # print(self.count_DA_cycle)
                        print(pathERT)
                        forward_mesh_vtk_file = glob.glob(
                            pathERT + "/**/*forward_model*", recursive=True)
                        # electrodes = glob.glob(pathERT + "/**/*electrodes.dat", recursive = True)
                        electrodes = glob.glob(
                            pathERT + "/**/*elecsXYZ.csv", recursive=True)
                        # seq = glob.glob(pathERT + "/**/*protocol.dat", recursive = True)
                        seq = glob.glob(pathERT + "/**/*ERT72.csv", recursive=True)
    
                        # df_sw = self.read_outputs(filename='sw',
                        #                           path=os.path.join(os.getcwd(),
                        #                                             self.output_dirname))
                        # df_sw[-1]
    
                        path_CATHY = os.path.join(os.getcwd(), 'vtk/')
    
                        porosity = self.soil['SPP'][:, 4][0]
                        Hx_ERT, df_Archie = Archie.SW_2_ERa(df_sw, 
                                                 path_fwd_ERT, 
                                                 self.project_name,
                                                 porosity,
                                                 # path_CATHY,
                                                 pathERT,
                                                 meshERT=forward_mesh_vtk_file[0],
                                                 elecs=electrodes[0],
                                                 sequenceERT=seq[0],
                                                 DA_cnb=self.count_DA_cycle,
                                                 Ens_nb=ens_nb)
    
                        # print(Hx_ERT)
                        self._add_2_ensemble_Hx(Hx, Hx_ERT['resist'])
                        
                        
                    if 'PH' in obs_key:
                        # case 1: pressure head assimilation (Hx_PH)
                        # --------------------------------------------------------------------
                        df_vp_PH = df_vp.table_pivot(
                            index=[time, node], value='PH')
                        # if key[0] in 'PH':
                        # filter the node if some of the instruments were not working
                        df_vp_PH_filt = df_vp_PH.iloc('node1', 'node2')
    
                        Hx_PH = df_vp_PH_filt
                          # Hx.vstack(Hx_PH)
    
                    if 'SW' in obs_key:
                        # case 2: sw assimilation (Hx_SW)
                        # --------------------------------------------------------------------
                        # if key[0] in 'SW':
    
                        df_vp_SW = df_vp.table_pivot(
                            index=[time, node], value='SW')
                        # filter the node if some of the instruments were not working
                        df_vp_SW_filt = df_vp_SW.iloc('node1', 'node2')
    
                        Hx_SW = df_vp_SW_filt * porosity
                        # Hx.vstack(Hx_SW)
    
                        # note: the value of the porosity can be unique or not depending on the soil physical properties defined
    
                    if 'discharge' in obs_key:
                        # case 3: discharge
                        # need to read the hgsfdet file (Hx_Q)
                        # --------------------------------------------------------------------
                        # if key[0] in 'discharge':
    
                        # derivation of the dircharge, Q from file 'hgsfdet'
                        Hx_Q = []
                        # Hx.vstack(Hx_Q)



        Hx = np.vstack(Hx).T

        return Hx  # , Hx_mean

    def _add_2_ensemble_Hx(self, Hx, Hx_2add):

        Hx.append(Hx_2add)

        pass


    def _read_state_ensemble(self):
        
        
        # read grid to infer dimension of the ensemble X
        # --------------------------------------------------------------------
        # self.grid3d = {}
        if len(self.grid3d) == 0:
            self.grid3d = in_CT.read_grid3d(
                os.path.join(self.workdir, self.project_name))

        M_rows = np.shape(self.grid3d['nodes_idxyz'])[0]
        N_col = self.NENS
        
        # Ensemble matrix X of M rows and N  + fill with psi values
        # --------------------------------------------------------------------
        X = np.zeros([M_rows, N_col])

        # state_var = 'pressure' # pressure head
        # read psi cathy output file to infer pressure head valies
        for j in range(self.NENS):
            df_psi = out_CT.read_psi(os.path.join(self.workdir, self.project_name,
                                               "DA_Ensemble/cathy_" + str(j+1),
                                               'output/psi'))
            # X[:,j] = df_psi['pressure']
            X[:, j] = df_psi[-1, :]

        # check if there is still zeros
        if np.count_nonzero(X) != M_rows*N_col:
            print('unconsistent filled X, missing values')
            
        return X, N_col, M_rows 

    def _prepare_ENKF(self, dict_obs, update_var_list=['St. var.']):
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        Matrice operations before ENKF

        Parameters
        ----------
        dict_obs : TYPE
            DESCRIPTION.
        update_var_list : TYPE
            DESCRIPTION.

        Returns
        -------
        A : np.array
            matrix of ensemble anomalies.
        D : np.array
            Difference between the measurements.
        X : np.array
            Ensemble matrix of M rows and N columns. N is the number of
            realizations and M is the state dimension, i.e., the number of nodes in the finite element grid
        X_augm : np.array
            Augmented by the number of parameters that are subject to up-
            date.
        x_mean : np.array
            mean of the ensemble.

        '''

        A = []
        D = []
        x_mean = []
        A_mean = []
        
    
        X, N_col, M_rows = self._read_state_ensemble()

        # Init Augmented Ensemble matrix
        # --------------------------------------------------------------------
        # combine the Ensemble and Param arrays.
        X_augm = X
        if len(update_var_list) > 1:
            for i in range(update_var_list-1):
                X_augm = X_augm.append(np.zeros([M_rows, 1]))

        # mean of the ensemble x_mean
        # --------------------------------------------------------------------
        # # Calculate mean
        # X_mean = (1./float(N_col))*np.tile(X.sum(1), (N_col,1)).transpose()

        ones = np.ones(N_col)
        XI = X_augm.dot(ones)
        # this is the mean value estimated over the whole number of scenarios
        X_augm_mean = (1/N_col)*(XI)

        # Matrix of ensemble anomalies A (ensemble perturbation from mean)
        # --------------------------------------------------------------------
        # these are the anomalies (difference)
        A = X_augm-(X_augm_mean.dot(np.ones(np.shape(X)[0])))

        # # Data perturbation from ensemble measurements
        # # --------------------------------------------------------------------
        # # D should be (MeasSize)x(EnSize)
        # D = Data - Observation
        # Ensemble mean of the simulated observations Hx
        # --------------------------------------------------------------------
        # Matrix of ensemble anomalies A
        # --------------------------------------------------------------------
        # Difference between the measurements D
        # --------------------------------------------------------------------

        return X, X_augm, A, A_mean

    def _create_subfolders_ensemble(self, NENS):
        '''
        THIS SHOULD BE MOVED TO DA CLASS


        Parameters
        ----------
        NENS : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        if not os.path.exists(os.path.join(self.workdir, self.project_name,
                                           "DA_Ensemble/cathy_origin")):
            os.makedirs(os.path.join(self.workdir, self.project_name,
                                     "DA_Ensemble/cathy_origin"))

            # copy input, output and vtk dir
            # ----------------------------------------------------------------
            for dir2copy in enumerate(['input', 'output', 'vtk']):
                shutil.copytree(os.path.join(self.workdir, self.project_name, dir2copy[1]),
                                os.path.join(self.workdir, self.project_name,
                                             "DA_Ensemble/cathy_origin", dir2copy[1])
                )

        try:
            # copy exe into cathy_origin folder
            # ----------------------------------------------------------------
            shutil.move(
                os.path.join(self.workdir, self.project_name,
                             self.processor_name),
                os.path.join(self.workdir, self.project_name,
                             "DA_Ensemble/cathy_origin", self.processor_name)
                        )

            # copy cathy.fnames into cathy_origin folder
            # ----------------------------------------------------------------
            shutil.copy(
                os.path.join(self.workdir, self.project_name,
                             'cathy.fnames'),
                os.path.join(self.workdir, self.project_name,
                             "DA_Ensemble/cathy_origin/cathy.fnames")
                        )

            # copy prepro into cathy_origin folder
            # ----------------------------------------------------------------
            shutil.copytree(os.path.join(self.workdir,  self.project_name,
                                         'prepro'),
                            os.path.join(self.workdir, self.project_name,
                                     "DA_Ensemble/cathy_origin/prepro"))

        except:
            self.console.print(":worried_face: [b]processor exe not found[/b]")

        path_origin = os.path.join(self.workdir, self.project_name,
                                   "DA_Ensemble/cathy_origin")

        # copy origin folder to each ensemble subfolders
        for i in range(NENS):
            path_nudn_i = os.path.join(self.workdir, self.project_name,
                                        "DA_Ensemble/cathy_" + str(i+1))

            if not os.path.exists(path_nudn_i):
                shutil.copytree(
                    path_origin, path_nudn_i
                )

    def _DA_df(self, X=None, Analysis=None):
        '''
        THIS SHOULD BE MOVED TO DA CLASS


        Parameters
        ----------
        NENS : TYPE
            DESCRIPTION.
        ENKFT : TYPE
            DESCRIPTION.
        var_per : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        DA_df_ti_ni = pd.DataFrame()
        DA_df_ti = pd.DataFrame()
        
        print(self.count_DA_cycle)
        if self.count_DA_cycle==0:
            self.DA_df = pd.DataFrame()

        #    for n in range(self.NENS):
        #         cols_root = ['time', 
        #                      'Ensemble_nb' +str(n),
        #                      'bef_update_' +str(self.count_DA_cycle), 
        #                      'aft_update_' +str(self.count_DA_cycle)]
                
                
        #         data_df_root = [self.count_DA_cycle*np.ones(len(X)), 
        #                         n*np.ones(len(X)), 
        #                         X[:,n],               
        #                         X[:,n]]
                
        #         DA_df_ti_ni = pd.DataFrame(np.transpose(data_df_root),
        #                               columns=cols_root)
                
        #         DA_df_ti= pd.concat([DA_df_ti, DA_df_ti_ni], axis=1)  
                
        else:

            for n in range(self.NENS):
                # cols_root = ['time', 
                #              'Ensemble_nb' +str(n),
                #              'bef_update_' +str(self.count_DA_cycle), 
                #              'aft_update_' +str(self.count_DA_cycle)]
                
                cols_root = ['time', 
                             'Ensemble_nb',
                             'bef_update_', 
                             'aft_update_']
                
                
                data_df_root = [self.count_DA_cycle*np.ones(len(X)), 
                                n*np.ones(len(X)), 
                                X[:,n],
                                Analysis[:,n]]
          
    
            
                DA_df_ti_ni = pd.DataFrame(np.transpose(data_df_root),
                                      columns=cols_root)
                
                DA_df_ti= pd.concat([DA_df_ti, DA_df_ti_ni], axis=0)   

               
               
            # print(self.DA_df)
            # print(DA_df_ti)
            # self.AA = self.DA_df
            # self.BB = DA_df_ti
            self.DA_df= pd.concat([self.DA_df, DA_df_ti], axis=0, ignore_index=True)
            # print(self.DA_df)

        # self.DA_var_pert_df['time'] = pd.to_timedelta(self.DA_var_pert_df['time'],unit='s')

        pass




    def _update_input_ensemble(self, NENS, ENKFT,
                               parm_pert=[],
                               update_parm_list=['St. var.'],
                               **kwargs):
        '''
        THIS SHOULD BE MOVED TO DA CLASS


        Select the time window of the hietograph.
        Write new variable perturbated/updated value into corresponding input file

        Parameters
        ----------
        NENS : int
            size of the ensemble.
        parm_pert : dict
            dict of all perturbated variables holding values and metadata.

        Returns
        -------
        New file written/overwritten into the input dir.

        '''

        self.console.print('_update_input_ensemble')

        # resample atmbc file for the given DA window
        # ---------------------------------------------------------------------
        self.selec_atmbc_window(NENS, ENKFT)

        # ---------------------------------------------------------------------
        self.update_ENS_files(parm_pert, update_parm_list=update_parm_list,
                              cycle_nb=self.count_DA_cycle)

        pass

    def selec_atmbc_window(self, NENS, ENKFT):
        '''
        THIS SHOULD BE MOVED TO DA CLASS


        Select the time window of the hietograph
        == time between two assimilation observation

        Parameters
        ----------
        NENS : TYPE
            DESCRIPTION.
        ENKFT : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        if len(self.grid3d) == 0:
            self.run_processor(IPRT1=3, verbose=False, DAFLAG=0)
            self.grid3d = in_CT.read_grid3d(
                os.path.join(self.workdir, self.project_name))

        # resample atmbc file for the given DA window
        # ----------------------------------------------------------------------
        # cathy and data assimilation are performed for a single time step
        # definition of the part of hyetograph to be applied in this step
        # self.update_cathyH() ?

        # read full simulation atmbc and filter time window
        # ----------------------------------------------------------------------
        print(os.path.join(self.workdir, self.project_name,
                              'input', 'atmbc'))
        df_atmbc, HSPATM, IETO = in_CT.read_atmbc(os.path.join(self.workdir, self.project_name,
                              'input', 'atmbc'), grid=self.grid3d)
        
        if self.count_DA_cycle is not None:
            try:
                time_window_atmbc = [
                    ENKFT[self.count_DA_cycle], ENKFT[self.count_DA_cycle+1]]
            except:
                pass
        else:
                time_window_atmbc = [ENKFT[0], ENKFT[1]]

            
        # try:
        #     # time_window_atmbc_df = [ENKFT.iloc[self.count_DA_cycl],ENKFT.iloc[self.count_DA_cycl+1]]
        #     time_window_atmbc = [
        #         ENKFT[self.count_DA_cycl], ENKFT[self.count_DA_cycl+1]]

        # # case where the simulation is not yet running (initial preparation)
        # except:
        #     print('try failed use 0 and 1')
        #     # time_window_atmbc_df = [ENKFT.iloc[0], ENKFT.iloc[1]]
        #     time_window_atmbc = [ENKFT[0], ENKFT[1]]

        # print(time_window_atmbc)
        # df_atmbc = df_atmbc.set_index('time')

        df_atmbc_window = df_atmbc[(df_atmbc.index >= time_window_atmbc[0]) &
                      (df_atmbc.index <= time_window_atmbc[1])]

        # df_atmbc_window['time'] = pd.to_timedelta(df_atmbc_window['time'],unit='s')

        # print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
        # print(self.count_DA_cycle)
        # print(ENKFT)
        # print(time_window_atmbc)
        # print(df_atmbc_window['time'])
        # print(df_atmbc_window['value'].to_numpy())
        # print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')

        # print(df_atmbc)
        # # need to resample if nb of atmbc points is too low with respect to data assimilation frequency
        # # ---------------------------------------------------------------------
        # resampling_mean = ENKFT[1]
        # df_atmbc_int = df_atmbc.set_index('time').resample(resampling_mean).mean().interpolate('linear')
        # df_atmbc_int_filt = df_atmbc_int[(df_atmbc_int.index > time_window_atmbc[0]) &
        #              (df_atmbc_int.index < time_window_atmbc[1])]

        # print(df_atmbc_int)
        # print(df_atmbc_int_filt)

        for ens_nb in range(NENS):

            # change directory according to ensmble file nb
            os.chdir(os.path.join(self.workdir, self.project_name,
                                  './DA_Ensemble/cathy_' + str(ens_nb+1)))

            # (this is always done since Kalman Filter is Sequential)
            # update atmbc file
            # --------------------------------------------------------------
            # convert to numpy and convert to second / np.timedelta64(1, 's')
            # self.update_atmbc(HSPATM=0,IETO=0,TIME=df_atmbc_int_filt.index/np.timedelta64(1, 's'),
            #                   VALUE=df_atmbc_int_filt['value'].to_numpy(),
            #                   show=False,
            #                   filename=os.path.join(os.getcwd(),'input/atmbc')) # specify filename path
            # self.update_atmbc(HSPATM=0,IETO=0,TIME=df_atmbc_int_filt.index/np.timedelta64(1, 's'),
            #                   VALUE=df_atmbc_int_filt['value'].to_numpy(),
            #                   show=False,
            #                   filename=os.path.join(os.getcwd(),'input/atmbc')) # specify filename path

            # time_window_atmbc_abs= [number - max(time_window_atmbc) for number in time_window_atmbc]

            # df_atmbc_window_abs = 
            
            # self.update_atmbc(HSPATM=1, IETO=1, TIME=time_window_atmbc_abs,
            #                   VALUE= np.tile(df_atmbc_window['value'].to_numpy(),
            #                                  (len(time_window_atmbc),1)),
            #                   show=False,
            #                   filename=os.path.join(os.getcwd(), 'input/atmbc'),
            #                   backup=True)  # specify filename path
            # "(TIMPRT(I),I=1,NPRT)"
            
            diff_time = time_window_atmbc[1] - time_window_atmbc[0]
            
            self.update_parm(TIMPRTi=[0,diff_time],TMAX=diff_time,
                             filename=os.path.join(os.getcwd(), 'input/parm'),
                             backup=True)

    def update_ENS_files(self, parm_pert, update_parm_list, **kwargs):
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        Overwrite ensemble files after analysis step

        Returns
        -------
        Overwritten files.

        '''

        update_key = 'ini_perturbation'
        if 'cycle_nb' in kwargs:
            if kwargs['cycle_nb'] > 0:
                update_key = 'update_nb' + str(kwargs['cycle_nb'])



        # loop over ensemble files to update ic --> 
        # this is always done since psi are state variable to update
        # ------------------------------------------------------------------
        for ens_nb in range(self.NENS):

            # change directory according to ensmble file nb
            os.chdir(os.path.join(self.workdir, self.project_name,
                                  './DA_Ensemble/cathy_' + str(ens_nb+1)))



            # state variable update use ANALYSIS update or not
            # --------------------------------------------------------------
            if 'St. var.' in update_parm_list:
                if kwargs['cycle_nb'] > 0:
                        self.update_ic(INDP=1, IPOND=0,
                                       pressure_head_ini=self.Analysis,
                                       filename=os.path.join(os.getcwd(), 'input/ic'),
                                       backup=True)
            else:
                
                print('no state variable update - use last iteration as initial conditions')
                df_psi = self.read_outputs(filename='psi',
                                          path=os.path.join(os.getcwd(),
                                                            self.output_dirname))
                



                if kwargs['cycle_nb'] > 0:
                        self.update_ic(INDP=1, IPOND=0,
                                       pressure_head_ini=df_psi[-1, :],
                                       filename=os.path.join(os.getcwd(), 'input/ic'),
                                       backup=True)
               

        # loop over dict of perturbated variable
        # ----------------------------------------------------------------------
        for key in update_parm_list:  # loop over perturbated variables dict

            # print(key[0])
            print(key)

            # loop over ensemble files
            # ------------------------------------------------------------------
            for ens_nb in range(self.NENS):

                # change directory according to ensmble file nb
                os.chdir(os.path.join(self.workdir, self.project_name,
                                      './DA_Ensemble/cathy_' + str(ens_nb+1)))

            # check parameter to update and update
            # ------------------------------------------------------------------

                # atmbc update
                # --------------------------------------------------------------
                if key in 'hietograph':
                    print('hietograph perturbated not yet implemented')

                    self.update_atmbc(verbose=True)

                # ic update
                # --------------------------------------------------------------
                if key in 'ic':
                    # check if key contain self.ic args # NOT YET IMPLEMENTED
                    # self.ic['INDP'] = 0
                    # self.ic['INDP'] = 0atmbc
                    self.update_ic(INDP=0, IPOND=0,
                                   WTPOSITION=parm_pert[key[0]
                                       ][update_key][ens_nb],
                                   verbose=True,
                                   filename=os.path.join(os.getcwd(), 'input/ic'))  # specify filename path

                # kss update
                # --------------------------------------------------------------
                if key in 'kss':
                    print('kss perturbated not yet implemented')

                    self.update_soil(verbose=True,
                                     filename=os.path.join(os.getcwd(), 'input/soil'))  # specify filename path

                # FeddesParam update
                # --------------------------------------------------------------
                elif key in ['PCANA', 'PCREF', 'PCWLT', 'ZROOT', 'PZ', 'OMGC']:

                    # parm_pert[key[0]][update_key][ens_nb] =
                    # print()
                    # print(self.count_DA_cycle)
                    # print(update_key)
                    # print(parm_pert[key[0]])

                    FeddesParam = {'PCANA': self.soil["PCANA"],
                                   'PCREF': self.soil["PCREF"],
                                   'PCWLT': self.soil["PCWLT"],
                                   'ZROOT': [parm_pert[key][update_key][ens_nb]],
                                   'PZ': self.soil["PZ"],
                                   'OMGC': self.soil["OMGC"]}
                    # print(FeddesParam)
                    self.update_soil(FP=FeddesParam,
                                     verbose=False,
                                     filename=os.path.join(os.getcwd(), 'input/soil'))  # specify filename path

        pass

    def _RMSE():
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        (Normalized) root mean square errors (NRMSEs)

        Returns
        -------
        None.

        '''

        # N0 = number of observation at a given time

        # RMSE_WC
        # RMSE_PH
        # RMSE_Q

        # RMSE[t]

        pass

    def read_observations(self, filename, data_type, data_err, show=False, **kwargs):
        '''
        read measures (real observations)

        Infer assimilation time for ENKF

        Parameters
        ----------
        filename : str
            filename of the observation dataset.
        data_type : str
            key tring to identify what type of measure is to read.
        data_err : float
            % of error for a given measurement dataset. The default is 0.05.
        show : Bool
            plot measure graph. The default is False.
        kwargs**:
            position : list
                xyz position of the observation

        Returns
        -------
        dict_obs : dict
            dict merging all observations + metadatas.
            First key level is data assimilation time
            Loop over data assimilation times means:
            for ti in self.dict_obs:
                print(ti)

        '''

        # specify assimilation time if not contained in the file
        # ---------------------------------------------------------------------
        tA = []
        if 'tA' in kwargs:
            tA = kwargs['tA']

        # discharge type read
        # ---------------------------------------------------------------------
        if data_type == 'discharge':
            df = in_meas.read_discharge(filename)
            units = '$m^{3}/s$'

        # tensiometer type read
        # ---------------------------------------------------------------------
        elif data_type == 'tensiometers':

            df = in_meas.read_discharge(filename)

            # add time stamp if not contained in the file

        elif data_type == 'ERT':

            df = in_meas.read_ERT(filename)
            units = '$\Omega$'

            # add time stamp if not contained in the file

        # no file specified (raise error)
        # ---------------------------------------------------------------------
        else:

            print('no file specified')

        # add node stamp if not contained in the file
        # ---------------------------------------------------------------------
        # nodes_of_interest, closest = self.find_nearest_node([[0,10,1],[10,0,1]])

        # store measure data (df) + metadata into a dict
        # ---------------------------------------------------------------------

        # if hasattr(self,dict_obs) is False:

        dict_obs_2add = {'filename': filename,
                        'data_type': data_type,
                        'units': units,  # units
                        'data': df,
                        'data_err': data_err,
                        'mesh_nodes': [],
                        'assimilation_times': tA
                        }
        dict_obs_2add = OrderedDict(dict_obs_2add)

        # concatenate dict
        # ---------------------------------------------------------------------

        dict_obs_new_time = OrderedDict()
        dict_obs_new_time[tA] = OrderedDict({})

        # print(dict_obs_new_time)
        for key in dict_obs_2add.keys():
            # print(key)
            dict_obs_new_time[tA][key] = dict_obs_2add[key]

        self.dict_obs = self._add_to_obs_dict(dict_obs_new_time)

        # %% Might be interesting to add a second level of key based on the method

        # dict_test = {}
        # dict_obs = {}
        # dict_obs_2add = {'filename': 'gdg2',
        #                 'data_type': 'gdg',
        #                 'units': 'fdgg', # units
        #                 }

        # for l in ['0','1']:
        #     dict_obs[l]= {}
        #     for m in ['ERT','Tensio']:
        #         dict_obs[l][m] = {}
        #         for key in dict_obs_2add.keys():
        #             dict_obs[l][m][key] = dict_obs_2add[key]

        #         dict_test = dict_test | dict_obs
        # # dict_test.get('keyA_1')

        # for dt in dict_test:
        #     print(dt)

        return self.dict_obs

    def _add_to_obs_dict(self, dict_obs_2add):

        self.dict_obs = self.dict_obs | dict_obs_2add
        # self.dict_obs.update(dict_obs_2add)

        return self.dict_obs


    def prepare_observations(self, dict_obs, Bishop=False):
        '''
        THIS SHOULD BE MOVED TO DA CLASS


        prepare observations before DA


        1. Measurement perturbation
        2. Compute matrice covariance and normalise data if necessary

        Parameters
        ----------
        dict_obs : dict
            dict containing measure data + metadata.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''


        self.init_obs_cov_matrice(dict_obs)

        # The EnKF algorithm implemented here is actually an en-
        # semble transform Kalman filter (Bishop et al., 2001) that
        # does not require the perturbation of observations. On the
        # other hand, the measurement error covariance matrix, R,
        # must be assumed to be known a priori.

        # (Bishop et al., 2001) --> not require the perturbation of observations
        # taking advantage of the high time resolution of the collected data.
        # ---------------------------------------------------------------------



        # unit conversion (?)
        # ------------------------------------------------------------------


        # When assimilating multiple variables, proper normaliza-
        # tion of the measurement error covariance matrices, anoma-
        # lies of the simulated data, and innovation vectors were per-
        # formed, using values of 0.6 m, 0.58, and 4.17 × 10 −5 m 3 s −1
        # for pressure head, water content and subsurface outflow,
        # respectively. The normalization ensures that in multivari-
        # ate assimilation scenarios the covariance matrices in the
        # Kalman gain are not ill-conditioned (Evensen, 2003; Cam-
        # porese et al., 2009b).


        # normalisation
        # ------------------------------------------------------------------


        if Bishop == True:


            # define measurement error covariance matrix, R
            # ------------------------------------------------------------------

            print('bishop')


        else:
        # pertubate measurements
        # ----------------------------------------------------------------------
            self.perturbate_obs()
            self.console.print('need to perturbate the measurements first')


        return self.dict_obs


    def perturbate_obs(self, dict_obs):

            if (self.DataCov.ndim == 0):
                DataPerturbation=np.sqrt(self.DataCov)*rn.randn(1, self.EnSize)
            elif (self.DataCov.ndim == 2):
                # Compute SVD of Data Covariance to generate noise
                U, s, V=np.linalg.svd(self.DataCov, full_matrices=False)
                DataPerturbation=np.dot(np.dot(U, np.diag(np.sqrt(s))),
                                          rn.randn(self.Data.shape[1], self.EnSize))
            else:
                print('Data Covariance should be matrix or scalar', '\n')

            DataArray=np.tile(
                self.Data[i, :], (self.EnSize, 1)).transpose() + DataPerturbation





    def init_obs_cov_matrice(self, dict_obs):
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        compute measurement cov matrice and add it to dict_obs

        Returns
        -------
        None.

        '''

        # Attribute to define the data-error covariance matrix.
        # In this example the data noise is assumed to be a scalar
        # that is constant at each observation. We use an array so
        # that the 'shape' of the data covariance is understood by
        # numpy in a linear algebra sense.


        for key in dict_obs:
            print(key)

        print(len(dict_obs))

        self.R=np.zeros([len(dict_obs), len(dict_obs)])

        # R:  observation error covariance matrix, R
        self.dict_obs['obs_cov_mat']=np.array(
            [[np.power(self.dict_obs['data_err'], 2.0)]])


        pass


    def read_inputs(self, filename):
        '''


        Parameters
        ----------
        filename : str
            name of the input CATHY file to read.

        Returns
        -------
        A dataframe or a dict describing file data/entries.

        '''

        if filename == 'atmbc':
            df=in_CT.read_atmbc(os.path.join(
                self.workdir, self.project_name, 'input', filename))
            return df


        else:
            print('no file specified')


        pass


    def read_outputs(self, filename, **kwargs):
        '''


        Parameters
        ----------
        filename : str
            name of the output file to read.

        Returns
        -------
        A dataframe or a dict describing file data/entries.

        '''

        path=os.path.join(self.workdir, self.project_name, 'output', filename)
        if 'path' in kwargs:
            path=kwargs['path'] + '/' + filename

        if filename == 'vp':
            df=out_CT.read_vp(path)
            return df
        if filename == 'hgraph':
            df=out_CT.read_hgraph(path)
            return df
        if filename == 'dtcoupling':
            df=out_CT.read_dtcoupling(path)
            return df
        if filename == 'hgsfdet':
            df=out_CT.read_hgsfdet(path)
            return df
        if filename == 'psi':
            df=out_CT.read_psi(path)
            return df
        if filename == 'sw':
            df=out_CT.read_sw(path)
            return df
        else:
            print('no file specified')


        pass




    # -------------------------------------------------------------------#
    # %% utils

    def find_nearest_node(self, node_coords):
        '''
        Find nearest mesh node

        Parameters
        ----------
        node_coords : list
            List of coordinates [[x,y,z],[x2,y2,z2]].

        Returns
        -------
        closest_idx : list
            Node indice in the grid3d.
        closest : list
            Node coordinate in the grid3d.

        '''

        if len(node_coords) < 1:
            node_coords=[node_coords]

        self.grid3d=in_CT.read_grid3d(self.project_name)

        closest_idx=[]
        closest=[]
        for i, nc in enumerate(node_coords):
            # euclidean distance
            d=((self.grid3d['nodes_idxyz'][:, 1] - nc[0]) ** 2 +
                   (self.grid3d['nodes_idxyz'][:, 2] - nc[1]) ** 2 +
                   (self.grid3d['nodes_idxyz'][:, 3] - nc[2]) ** 2
                   ) ** 0.5

            closest_idx.append(np.argmin(d))
            closest.append(self.grid3d['nodes_idxyz'][closest_idx[i], 1:])

            threshold=1e-1
            if d[np.argmin(d)] > threshold:
                self.console.print('no node close to the required points')
                print(d[np.argmin(d)])
                print(self.grid3d['nodes_idxyz'][closest_idx[i], 1:])
                print(nc)

        return closest_idx, closest


    def rich_create(self, **kwargs):

        # from rich.console import Console
        # from rich.table import Table

        # table = Table(title="Star Wars Movies")

        # table.add_column("Released", justify="right", style="cyan", no_wrap=True)
        # table.add_column("Title", style="magenta")
        # table.add_column("Box Office", justify="right", style="green")

        # table.add_row("Dec 20, 2019", "Star Wars: The Rise of Skywalker", "$952,110,690")
        # table.add_row("May 25, 2018", "Solo: A Star Wars Story", "$393,151,347")
        # table.add_row("Dec 15, 2017", "Star Wars Ep. V111: The Last Jedi", "$1,332,539,889")
        # table.add_row("Dec 16, 2016", "Rogue One: A Star Wars Story", "$1,332,439,889")

        dict_name='test'

        dict_rich={
            'definition': None,
            'perturbation_dict': None,
            'damping_value': None,
            'assimilated': False,
            'assimilation_times': None
        }

        self.console.print(dict_rich)



    def rich_display(self, title="Star Wars Movies", **kwargs):
        """
        Describe the variable state and fate during the simulation with a rich table

        Returns
        -------
        None.

        """

        # Assimilated Yes/No

        # from rich.console import Console
        # from rich.table import Table

        # table = Table(title=title)

        # console = Console(record=True)
        self.console.print(eval('self.' + str(title)))

        # console.save_text('test.txt')
        # console.save_html('test.html')



    def create_3dmesh_CATHY(
        self,
        gmsh_mesh=[],
        NZONE=[],
        NSTR=[],
        N1=[],
        NNOD=None,
        NTRI=None,
        ZRATIO=[],
        Z1=[],
        IVERT=0,
        ISP=0,
        BASE=4,
        debug=False,
    ):
        """
        Create 3d mesh (grid file) from gmsh file.

        Parameters ---------- gmsh_mesh : type Description of parameter `gmsh_mesh`. NZONE : type #
        of material types in the porous medium. NSTR : type The number of vertical layers. N1 :
        type The maximum number of element connections to a node. NNOD : type Description of
        parameter `NNOD`. NTRI : type Description of parameter `NTRI`. ZRATIO : type The thickness
        of vertical layers or the fraction of total grid height that each layer is to occupy
        (ZRATIO (1) is for the surface‐most layer. ZRATIO values must sum to 1.). Z1 : type
        Description of parameter `Z1`. IVERT : type =0 each layer will be parallel to the surface,
        including the base of the 3‐d grid. `ZRATIO` is applied to each vertical cross section. =1
        base of the 3‐d grid will be flat, and `ZRATIO` is applied to each vertical cross section
        =2 base of the 3‐d grid will be flat, as will the NSTR‐1 horizontal cross sections above
        it. `ZRATIO` is applied only to the vertical cross section having the lowest elevation. =3
        for each cell of the dem a single depth value is read in file input IIN60 (basement).
        `ZRATIO` is applied to each vertical cross section. =4 the first NSTR‐1 layers from the
        surface will be parallel to the surface and the base of the 3‐d grid will be flat. `ZRATIO`
        is applied only to the vertical cross section having the lowest elevation. ISP : type =0
        for flat surface layer (only one Z value is read in, and is replicated to all surface
        nodes); otherwise surface layer is not flat (Z values read in for each surface node); (for
        ISP=0, IVERT=0, 1, and 2 yield the same 3‐d mesh, given the same values of BASE and
        ZRATIO). BASE : type Value which defines the thickness or base of the 3‐d mesh. For
        `IVERT`=0, BASE is subtracted from each surface elevation value, so that each vertical
        cross section will be of thickness BASE, and the base of the 3‐d mesh will be parallel to
        the surface. For IVERT=1 or 2, BASE is subtracted from the lowest surface elevation value,
        say ZMIN, so that each vertical cross section will be of thickness (Z ‐ ZMIN) + BASE, where
        Z is the surface elevation for that cross section. The base of the 3‐d mesh will thus be
        flat. debug : type Description of parameter `debug`.

        Returns ------- type Description of returned object.

        """
        k=Project(
            typ="R2"
        )  # create a Project object in a working directory (can also set using k.setwd())
        k.importMesh(gmsh_mesh)
        k.showMesh()
        mesh_dict, con_matrix, phys_entity=mt.mshParse(gmsh_mesh)

        # insert all the elements of the superficial mesh i.e.
        # the headers, the 4th columns desribing the triangles and the nodes coordinates

        if NNOD is None:
            print("NNOD")

        print(ZRATIO)

        with open("grid", "w+") as gridfile:
            gridfile.write(
                "\t"
                + str(NZONE)
                + "\t"
                + str(NSTR)
                + "\t"
                + str(N1)
                + "\t"
                + "NZONE"
                + "\t"
                + "NSTR"
                + "\t"
                + "N1"
                + "\n"
            )
            gridfile.write(
                "\t"
                + str(NNOD)
                + "\t"
                + str(NTRI)
                + "\t"
                + "NNOD"
                + "\t"
                + "NTRI"
                + "\n"
            )
            gridfile.write(
                "\t"
                + str(IVERT)
                + "\t"
                + str(ISP)
                + "\t"
                + str(BASE)
                + "\t"
                + "IVERT"
                + "\t"
                + "ISP"
                + "\t"
                + "BASE (m) "
                + "\n"
            )
            # gridfile.write("\t" + '\t'.join(ZRATIO[1:]) +
            #                + 'ZRATIO' + "\n")
            np.savetxt("grid", list(ZRATIO))
            gridfile.write("ZRATIO" + "\n")

            gridfile.write("\t" + str(Z1) + "\t" + "Z(1) m" + "\n")

            # np.savetxt('grid', np.transpose([mesh_dict['node_x'], mesh_dict['node_y']]),fmt='%1.2f')
            np.savetxt(
                "grid",
                np.transpose([mesh_dict["node_x"], mesh_dict["node_y"]]),
                fmt="%1.2f",
            )
            gridfile.close()

        # for k, v in mesh_dict.items():  # TypeError: 'list' object is not callable
        #     print(k)
        # print(mesh_dict['elm_id'])
        # gridfile.write("durin's day\n")
        # copy mesh attribute to 'grid' files
    pass
