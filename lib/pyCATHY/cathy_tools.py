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
import pickle 
import re

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
from pyCATHY.DA import cathy_DA, enkf

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
    # p.close()


    return p



class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
        

from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout
            
            
    
# -----------------------------------------------------------------------------

class CATHY():  # IS IT GOOD PRACTICE TO PASS DA CLASS HERE ? I think we sould better pass main CATHY into children classes
    """ Main CATHY object

      When instantiated it creates the tree project directories with 'prj_name' as root folder.
      The src files are fetched from the online repository if not existing
      (note that it is possible to call a specific version).

      """

    def __init__(self, dirName=None, prj_name="my_cathy_prj", notebook=False,
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

        self.project_name = prj_name

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
        self.stacked_data_cov = [] # merged data covariance matrice of all the observation data
        self.count_DA_cycle = None  # counter of the Assimilation Cycle
        # self.df_performance = pd.DataFrame() # pandas dataframe with performence evaluation metrics

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
            
        # p.close()
        if "DEM resolution is not a multiple of the rivulet spacing!" in p.stdout:
            print("DEM resolution is not a multiple of the rivulet spacing!")
            sys.exit()
        
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
        '''
        Run cathy.exe


        1. updates cathy parameters and cathyH based on **kwargs

        2. recompile the source files (set False for notebook) and create the executable

        3. run the processor using bash cmd if runProcess is True

        Parameters
        ----------
        recompile : TYPE, optional
            DESCRIPTION. The default is True.
        runProcess : TYPE, optional
            DESCRIPTION. The default is True.
        verbose : TYPE, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        # set parallel flag for mutirun (used for DA ensemble)
        # --------------------------------------------------------------------
        parallel = False
        if 'parallel' in kwargs:
                parallel = kwargs['parallel']

        if 'DAFLAG' in kwargs:
            self.DAFLAG = kwargs['DAFLAG']

        # update parm and cathyH
        # --------------------------------------------------------------------
        # if len(kwargs.items()) > 0:  # if new arguments update parm and CATHYH
        # self.update_parm(**kwargs, verbose=verbose)
        
        # VERY VERY IMPORTANT NEVER COMMENT !
        self.update_parm(**kwargs, verbose=verbose)
        self.update_cathyH(**kwargs)


        if recompile == True:

            # recompile
            # --------------------------------------------------------------------
            os.chdir(os.path.join(self.workdir, self.project_name, "src"))
            
            # print(os.getcwd())

            self.recompileSrc()
            # self.console.print(":hammer_and_wrench: [b]Recompile src files[/b]")

        if 'path_CATHY_folder' in kwargs:
            os.chdir(kwargs['path_CATHY_folder'])
        else:
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

                # define Data Assimilation default parameters
                # -------------------------------------------------------------

                # type of assimillation
                DA_type = 'enkf_Evensen2009_Sakov' # Default DA type
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
                list_assimilated_obs = 'all'
                if 'list_assimilated_obs' in kwargs:
                    list_assimilated_obs = kwargs['list_assimilated_obs']

                open_loop_run = True
                if 'open_loop_run' in kwargs:
                    open_loop_run = kwargs['open_loop_run']
                    
                    
                    
                self._run_DA(callexe, parallel, 
                             DA_type,
                             dict_obs,
                             list_update_parm,
                             dict_parm_pert,
                             list_assimilated_obs,
                             verbose=True,
                             open_loop_run=open_loop_run
                             )

            # case of simple simulation
            # ----------------------------------------------------------------
            else:
                p = subprocess.run([callexe], text=True, capture_output=True)

                if verbose == True:
                    print(p.stdout)
                    print(p.stderr)
                
                # p.close()

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

    
    def _parse_ERT_metadata(self,key_value):
        
        # Load ERT metadata information form obs dict
        # -------------------------------------------
        forward_mesh_vtk_file = key_value[1]['forward_mesh_vtk_file'] 
        pathERT = os.path.split(self.dict_obs[0]['filename'])[0]
        seq = key_value[1]['sequenceERT'] 
        electrodes = key_value[1]['elecs'] 
        porosity = self.soil_SPP['SPP'][:, 4][0]
        
        return forward_mesh_vtk_file, pathERT, seq, electrodes, porosity
        
        
    def _map_states2Observations_ERT_parallel(self, 
                                              path_fwd_CATHY_list,
                                              list_assimilated_obs='all', 
                                              verbose = False,
                                              **kwargs):
        

        savefig = True
        if 'savefig' in kwargs:
           savefig = kwargs['savefig'] 
        
        ENS_times = []
        if 'ENS_times' in kwargs:
           ENS_times = kwargs['ENS_times'] 
        
        DA_cnb = []
        if 'DA_cnb' in kwargs:
           DA_cnb = kwargs['DA_cnb'] 
           
           
        # search key value to identify time and method
        tuple_list_obs = list(self.dict_obs.items())
        key_value = tuple_list_obs[self.count_DA_cycle]
        
        
        # Load ERT metadata information form obs dict
        # -------------------------------------------
        (forward_mesh_vtk_file, 
         pathERT, seq, electrodes, 
         porosity) = self._parse_ERT_metadata(key_value)
        
        
        

        if len(ENS_times)>0: # case of the open Loop
        
        
            prediction = []
            for t in range(len(ENS_times)-1):
    
                
                # freeze fixed arguments of Archie.SW_2_ERa       
                # -----------------------------------------------------------------
                ERTmapping_args = partial(Archie.SW_2_ERa, 
                                          self.project_name,
                                          porosity, 
                                          pathERT,
                                          forward_mesh_vtk_file,
                                          electrodes,
                                          seq,
                                          data_format= self.dict_obs[0]['data_format'],
                                          time_ass = t,
                                          savefig=True
                                          ) 
                
            
                with suppress_stdout():

                    # // run using ensemble subfolders path as a list    
                    # -----------------------------------------------------------------
                    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
                                                  
                        results_mapping_time_i = pool.map(ERTmapping_args, 
                                                    path_fwd_CATHY_list
                                                      )
                        
                        print(f"x= {path_fwd_CATHY_list}, PID = {os.getpid()}")
            
                        # print(len(results_mapping_time))
                        # print(results_mapping_time[0][0].shape)
                    
                    # Hx_ERT_time_i =  results_mapping_time_i[1][0]
                    
                    Hx_ERT_ens = []
                    for ens_i in range(len(path_fwd_CATHY_list)):
                        
                        Hx_ERT_time_i =  results_mapping_time_i[ens_i][0]
            
                        if 'pygimli' in self.dict_obs[0]['data_format']:
                            Hx_ERT_ens = self._add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_time_i['rhoa'])
                        else:
                            Hx_ERT_ens = self._add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_time_i['resist'])
                           
                    prediction.append(np.vstack(Hx_ERT_ens).T)  # (times * data size * EnSize)
                           
                                           
        else:
            
             Hx_ERT_ens = []

             # freeze fixed arguments of Archie.SW_2_ERa       
             # -----------------------------------------------------------------
             ERTmapping_args = partial(Archie.SW_2_ERa, 
                                       self.project_name,
                                       porosity, 
                                       pathERT,
                                       forward_mesh_vtk_file,
                                       electrodes,
                                       seq,
                                       data_format= self.dict_obs[0]['data_format'],
                                       DA_cnb = DA_cnb,
                                       savefig=True
                                       ) 
             
            
             # // run using ensemble subfolders path as a list    
             # -----------------------------------------------------------------
             with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
                                           
                results_mapping = pool.map(ERTmapping_args, 
                                             path_fwd_CATHY_list
                                               )
                 
                print(f"x= {path_fwd_CATHY_list}, PID = {os.getpid()}")
                 
                len(results_mapping)
             
                Hx_ERT_ens = []
                for ens_i in range(len(path_fwd_CATHY_list)):
                    
                    Hx_ERT_i =  results_mapping[ens_i][0]
        
                    if 'pygimli' in self.dict_obs[0]['data_format']:
                        Hx_ERT_ens = self._add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_i['rhoa'])
                    else:
                        Hx_ERT_ens = self._add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_i['resist'])
                 
                prediction = np.vstack(Hx_ERT_ens).T  # (EnSize)
            
        return prediction
                
                
    def _DA_openLoop(self,ENS_times,list_assimilated_obs,
                     parallel=True,
                     verbose=False):
        '''
        Run open Loop (no update) hydro simulation for an ensemble of realisation
        Evaluate the performance by comparison with measured data after mapping

        Parameters
        ----------
        ENS_times : TYPE
            DESCRIPTION.
        list_assimilated_obs : TYPE
            DESCRIPTION.
        run : TYPE, optional
            DESCRIPTION. The default is True.
        parallel : TYPE, optional
            DESCRIPTION. The default is True.
        verbose : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        '''
        
             

        # multi run CATHY hydrological model from the independant folders 
        # composing the ensemble
        # ----------------------------------------------------------------
        if parallel == True:
            pathexe_list = []
            for ens_i in range(self.NENS):
                path_exe = os.path.join(self.workdir,
                                        self.project_name,
                                        'DA_Ensemble/cathy_' + str(ens_i+1))
                pathexe_list.append(path_exe)

        
    
            with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
                result = pool.map(subprocess_run_multi, pathexe_list)
                
                if verbose == True:
                    self.console.print(result)
                    # self.console.print(result.stderr)
        
        # Loop over ensemble realisations one by one
        # ------------------------------------------------------------
        else:
            for ens_i in range(self.NENS):
                os.chdir(os.path.join(self.workdir,
                                      self.project_name, 
                                      'DA_Ensemble/cathy_' + str(ens_i+1)))

                self._run_hydro_DA_openLoop(time_of_interest=ENS_times,
                                  nodes_of_interest=[],
                                  simu_time_max=max(ENS_times),
                                  ens_nb=ens_i+1                                 
                                  )
                

        # save into the DA_df dataframe
        # -----------------------------------------------------------------
        for ens_i in range(self.NENS):

            os.chdir(os.path.join(self.workdir,
                                  self.project_name, 
                                  'DA_Ensemble/cathy_' + str(ens_i+1)))
            
            df_psi = self.read_outputs(filename='psi',
                                       path=os.path.join(os.getcwd(),
                                                         'output'))
            
            shift = len(df_psi)-len(ENS_times)
            for t in range(len(ENS_times)):
                self._DA_df(state=df_psi[t+shift,:], 
                            t_ass=t, 
                            openLoop=True,
                            ens_nb=ens_i+1)
                

        # map states 2 observation for performance evaluation
        # -----------------------------------------------------------------
        parallel_mapping = False
        if parallel_mapping==True: # // for each ensemble steps (NOT FOR EACH time steps)
            
            # list path to run in parallel //
            path_fwd_CATHY_list = []
            for ens_i in range(self.NENS):
                
                path_fwd_CATHY_list.append(os.path.join(self.workdir,
                                            self.project_name, 
                                            'DA_Ensemble/cathy_' + str(ens_i+1)))
                
                           
            
            ERT_tmp = False
            if ERT_tmp == True:
                # _map_states2Observations_ERT_par is looping over times
                # -------------------------------------------------------------
                prediction = self._map_states2Observations_ERT_parallel(path_fwd_CATHY_list, 
                                                                         ENS_times=ENS_times,
                                                                         savefig=True)

            # prediction is a 3d matrice of dimension (data_size * ens_size * assimilation_times_size) 
    
            # print(np.shape(predictions))
            # input('press to continue')
    
            # performance assessement for each observation steps
            # ----------------------------------------------------------
            for ens_i in range(self.NENS):
                
                for t in range(len(ENS_times)-1):
                    
                    key_value =  list(self.dict_obs.items())[t]     
                    if 'pygimli' in self.dict_obs[0]['data_format']:
                        data = key_value[1]['data']['rhoa'].to_numpy()  # (MeasSize)
                    else:
                        data = key_value[1]['data']['resist'].to_numpy()  # (MeasSize)
                    

                    self._performance_assessement(list_assimilated_obs, 
                                                  data, 
                                                  prediction[t], 
                                                  t_obs=t,
                                                  openLoop=True)
                    
                    
            
            
            
        else: # Don't run parallel mapping
        
            Hx_ens = [] # predicted observations
            for ens_i in range(self.NENS):

                # performance assessement for each observation steps
                # -------------------------------------------------------------
                for t in range(len(ENS_times)-1):
                    
                    path_fwd_CATHY = os.path.join(self.workdir,
                                          self.project_name, 
                                          'DA_Ensemble/cathy_' + str(ens_i+1))      
                    df_sw = self.read_outputs(filename='sw',
                                              path=os.path.join(path_fwd_CATHY,
                                                                self.output_dirname)
                                                     )                    
                    
                    # shift between vtkfiles and real assimilation times
                    # CATHY produce 100.vtk and 101.vtk indentical
                    # ---------------------------------------------------------
                    shift = len(df_sw)-len(ENS_times)

                    # Apply operator H on simulated data x --> Hx
                    # ---------------------------------------------------------
                    # Hx must be sensor size
                    Hx = self._map_states2Observations([None, df_sw[t+shift,:]],
                                                       list_assimilated_obs,
                                                       path_fwd_CATHY=path_fwd_CATHY,
                                                       ens_nb=ens_i+1,
                                                       savefig=False) 
                    
                    # Hx_ens must (EnSize*ENS_times) * sensors nb 
                    Hx_ens = self._add_2_ensemble_Hx(Hx_ens, Hx)
                    
            prediction = np.vstack(Hx_ens) # (EnSize)
            
            # nb of data x ensemble size x nb of obs times
            prediction_times = np.reshape(Hx_ens,[self.NENS,len(Hx),len(ENS_times)-1])           
            
            
            # Performance evaluation
            # ----------------------
            for t in range(len(ENS_times)-1):
                 
                # Loop trought observation dictionnary for a given assimilation time (t)
                # -----------------------------------------------------------------
                data = []    # data dict 2 map
                if list_assimilated_obs == 'all':
                    items_dict = list(self.dict_obs.items())
                    # There can be multiple sensors measuring at the same time
                    # -----------------------------------------------------------------
                    for sensors in items_dict[t][1].keys():                
                        data.append(items_dict[t][1][sensors]['data'])
                
                # key_value =  list(self.dict_obs.items())[t]     
                # if 'pygimli' in self.dict_obs[0]['data_format']:
                #     data = key_value[1]['data']['rhoa'].to_numpy()  # (MeasSize)
                # else:
                #     data = key_value[1]['data']['resist'].to_numpy()  # (MeasSize)
                        

                self._performance_assessement(list_assimilated_obs, 
                                              data, 
                                              prediction_times[:,:,t], 
                                              t_obs=t,
                                              openLoop=True)
                
        # ------------------------------------------------------
        # END of Open Loop simulation and performance evaluation
        # ------------------------------------------------------     
        pass
    
    
    
    
    def _run_DA(self, callexe, parallel,
                        DA_type,
                        dict_obs,
                        list_update_parm,
                        dict_parm_pert,
                        list_assimilated_obs,
                        verbose,
                        **kwargs
                        ):
        '''
        
        Run Data Assimilation 
        
        .. note::

            1. DA init (create subfolders)

            2a. run CATHY hydrological model (open loop)
            2b. run CATHY hydrological model recursively using DA times

            3. check before analysis

            4. analysis

            5. check after analysis

            6. update input files

        Parameters
        ----------
        callexe : str
            processor exe filename
        parallel : bool
            True for multiple ensemble folder simulations.
        DA_type : str
            type of Data Assimilation.
        list_assimilated_obs : TYPE
            list of names of observations to assimilate.
        dict_obs : dict
            observations dictionnary.
        list_update_parm : list
            list of names of parameters to update.
        dict_parm_pert : dict
            parameters perturbated dictionnary.


        '''


        # Initiate
        # -------------------------------------------------------------------
        update_key = 'ini_perturbation'
        self.dict_obs = dict_obs # self.dict_obs is already assigned (in read observatio! change it to self.obs
        self.dict_parm_pert = dict_parm_pert
        self.df_DA = pd.DataFrame()

        
        # Infer ensemble size NENS from perturbated parameter dictionnary
        # -------------------------------------------------------------------
        for name in self.dict_parm_pert:
            NENS = len(self.dict_parm_pert[name]['ini_perturbation'])

        # Infer ensemble update times ENS_times from observation dictionnary
        # -------------------------------------------------------------------
        ENS_times = []
        for ti in self.dict_obs:
            ENS_times.append(float(ti))
            
        if self.parm['TMAX']> max(ENS_times):
            ENS_times.append(self.parm['TMAX'])
        # else:
        #     pass
        
        # start DA cycle counter
        # -------------------------------------------------------------------
        self.count_DA_cycle = 0
        # (the counter is incremented during the update analysis)

        # initiate DA
        # -------------------------------------------------------------------
        self._DA_init(NENS=NENS, # ensemble size
                      ENS_times=ENS_times, # assimilation times
                      parm_pert=dict_parm_pert)
        
        # input("Please press the Enter key to proceed")

        # open Loop for comparison
        # -------------------------------------------------------------------
        # openLoopDA = True
        
        # update the perturbated parameters 
        # --------------------------------------------------------------------
        self.update_ENS_files(dict_parm_pert, 
                              update_parm_list='all',
                              cycle_nb=self.count_DA_cycle)
        
        
        if kwargs['open_loop_run']:
            self._DA_openLoop(ENS_times,list_assimilated_obs)
       
                
        # end of Open loop - start DA
        # -------------------------------------------------------------------

        # update input files ensemble again (time-windowed)
        # ---------------------------------------------------------------------
        self._update_input_ensemble(NENS, 
                                    ENS_times, 
                                    dict_parm_pert, 
                                    update_parm_list='all') 
        
        
        # -----------------------------------
        # Run hydrological model sequentially
        # -----------------------------------
        
        # self.sequential_DA(update_key)       
        # def sequential_DA(self,update_key):
            
            
        # Loop over assimilation cycle
        # -------------------------------------------------------------------
        for DA_cnb in range(len(ENS_times)-1):


            # multi run CATHY hydrological model from the independant folders composing the ensemble
            # ----------------------------------------------------------------
            if parallel == True:
                pathexe_list = []
                for ens_i in range(self.NENS):
                    path_exe = os.path.join(self.workdir,
                                            self.project_name,
                                            'DA_Ensemble/cathy_' + str(ens_i+1))
                    pathexe_list.append(path_exe)

            
        
                with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
                    result = pool.map(subprocess_run_multi, pathexe_list)
                    
                    if verbose == True:
                        self.console.print(result)
                        # self.console.print(result.stderr)
                        
                        

            # process each ensemble folder one by one
            # ----------------------------------------------------------------
            else:
                # Loop over ensemble realisations
                # ------------------------------------------------------------
                for ens_i in track(range(self.NENS), description="Run hydrological fwd model..."):

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
                        
                    # p.close()

            os.chdir(os.path.join(self.workdir))
            
            

            # check scenario (result of the hydro simulation)
            # ----------------------------------------------------------------
            self._check_before_analysis(update_key)
            # and then create and applied the filter
            
            

            # map states to observation = apply H operator to state variable
            # ---------------------------------------------------------------------
            
            
            # self._map_states2Observations()
            parallel_mapping = False
            if parallel_mapping == True:
                print('not yet tested')
                
                path_fwd_CATHY_list = []
                for ens_i in range(self.NENS):
                    
                    path_fwd_CATHY_list.append(os.path.join(self.workdir,
                                                self.project_name, 
                                                'DA_Ensemble/cathy_' + str(ens_i+1)))
                
                # _map_states2Observations_ERT_par
                # -------------------------------------------------------------
                prediction = self._map_states2Observations_ERT_parallel(path_fwd_CATHY_list,
                                                                        savefig=True,
                                                                        DA_cnb=self.count_DA_cycle)
                
                
            else:
                
                Hx_ens = []
                for ens_nb in range(self.NENS): # loop over ensemble files
        
                    path_fwd_CATHY = os.path.join(self.workdir, self.project_name,
                                          './DA_Ensemble/cathy_' + str(ens_nb+1))
        
                    df_sw = self.read_outputs(filename='sw',
                                              path=os.path.join(path_fwd_CATHY,
                                                                self.output_dirname)
                                             )
                    backup = True
                    if backup == True:
                        dst_dir = os.path.join(path_fwd_CATHY,'output/sw') + str(self.count_DA_cycle)
                        shutil.copy(os.path.join(path_fwd_CATHY,'output/sw'),dst_dir) 
        
                        dst_dir = os.path.join(path_fwd_CATHY,'output/psi') + str(self.count_DA_cycle)
                        shutil.copy(os.path.join(path_fwd_CATHY,'output/psi'),dst_dir) 
                    
                    

                    # Apply operator H on measured data x --> Hx
                    # ---------------------------------------------------------         
                    Hx = self._map_states2Observations([None, df_sw[-1]], # take the last iteration of the state matrice
                                                       list_assimilated_obs,
                                                       path_fwd_CATHY=path_fwd_CATHY,
                                                       ens_nb=ens_nb)
                    
                    Hx_ens.append(Hx)

                prediction = np.vstack(Hx_ens).T 
                prediction = prediction.reshape(self.NENS,len(Hx)) # (EnSize*Data size)

        
        
            # analysis step
            # ----------------------------------------------------------------
            
            
            (
             ensembleX, 
             data, 
             analysis, 
             analysis_param) = self._DA_analysis(prediction,
                                                 DA_type,
                                                 list_update_parm,
                                                 list_assimilated_obs='all')
                
                
                
            # plt_CT.show_RMSE()
            self._check_after_analysis(update_key,list_update_parm)
            
            
            # check analysis quality
            # ----------------------------------------------------------------
            self._performance_assessement(list_assimilated_obs, 
                                            data, 
                                            prediction,
                                            t_obs=self.count_DA_cycle)
            
            
            # the counter is incremented here                           
            # ----------------------------------------------------------------
            self.count_DA_cycle = self.count_DA_cycle + 1
            
            
            # update step
            # ----------------------------------------------------------------
            if self.count_DA_cycle > 0:
                update_key = 'update_nb' + str(self.count_DA_cycle)
    
            
            # update parameter dict with new ones
            # ----------------------------------------------------------------
            for pp in enumerate(list_update_parm[:]):
                if 'St. var.' in pp[1]:
                    pass
                else:
                    self.dict_parm_pert[pp[1]][update_key] = analysis_param
                    # self.dict_parm_pert = self.dict_parm_pert | var_per_2add # only for python 3.9

            
    
            # create dataframe _DA_var_pert_df holding the results of the DA update
            # ---------------------------------------------------------------------
            self._DA_df(state=ensembleX, 
                        state_analysis=analysis)
    
    
            # overwrite input files ensemble (perturbated variables)
            # ---------------------------------------------------------------------
            
            if self.count_DA_cycle<len(ENS_times)-1: # -1 cause ENS_Times include TMAX
            
                self._update_input_ensemble(self.NENS, 
                                            ENS_times, 
                                            self.dict_parm_pert,
                                            update_parm_list=list_update_parm, 
                                            analysis=analysis)

            else:
                print('------ end of DA ------')
                pass

    
            print('------ end of update ------')
            
            
        # export summary results of DA
        # ----------------------------------------------------------------
        self.backup_results()
        
    
   
    def backup_results(self):
        '''
        

        Returns
        -------
        None.

        '''
        
        
        with open(os.path.join(self.workdir, self.project_name + '.pkl'), 'wb') as f:
            pickle.dump(self.dict_parm_pert,f)
            pickle.dump(self.df_DA,f)
            pickle.dump(self.df_performance,f)
        
        f.close()
        
        

    def _run_hydro_DA_openLoop(self,time_of_interest,nodes_of_interest,simu_time_max, ens_nb):
        '''
        Run openLoop and save result into DA dataframe (ensemble nb = 999)

        Parameters
        ----------
        time_of_interest : TYPE
            DESCRIPTION.
        nodes_of_interest : TYPE
            DESCRIPTION.
        simu_time_max : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        
        cwd = os.getcwd()

        self.console.print(
            ":unlock: [b]open loop call[/b]")
        self.run_processor(recompile=True,
                           TIMPRTi=time_of_interest, 
                           NODVP=nodes_of_interest, 
                           TMAX=simu_time_max,
                           DAFLAG=0,
                           path_CATHY_folder= cwd
                           ) #TMAX=simu_time_max
        
        df_psi = self.read_outputs(filename='psi',
                                   path=os.path.join(cwd, 'output'))
            
        for t in range(len(time_of_interest)):
            self._DA_df(state=df_psi[t,:], 
                        t_ass=t, 
                        openLoop=True,
                        ens_nb=ens_nb)
            
        

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
            # sys.exit()

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
        # Loop over hap.in file lines
        # --------------------------------------------------------------------
        for i, line in enumerate(Lines):
            xnew = line.strip().split("=")
            if len(xnew) > 1:  # skip lines with no =
            
                # loop over hapin dictionnary
                # -----------------------------------
                for key, value in self.hapin.items():
                    if count < len(hapin): 
                    
                        if key == hapin[tmp_lnb[count]]:
                            if value != tmp_param_value[count]:
                                print(key)
                                print(value)
                                if isinstance(value, list):
                                    value_str = "/t".join(value)
                                    xnew[1] = value_str
                                else:
                                    xnew[1] = value
                                line = xnew[0] + "=              " + \
                                    str(xnew[1]) + "\n"
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

            # print("update dem")
            # with open(os.path.join(self.project_name, "prepro/dem"), "w+") as f:
            #     # use exponential notation
            #     np.savetxt(f, self.DEM, fmt="%1.4e")
                
        self.update_dem_parameters(**kwargs)

        if show == True:
            plt_CT.dem_plot(self.workdir, self.project_name,
                            delta_x=self.hapin['delta_x'],
                            delta_y=self.hapin['delta_y'])

        pass

    def update_dem_parameters(self, **kwargs):
        '''
        

        Parameters
        ----------
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''        

        self.console.print(":black_nib: [b]update dem_parameters file [/b]")

        # set default parameters
        # --------------------------------------------------------------------
        if hasattr(self, "dem_parameters") == False:
            self.console.print(
                ":pensive_face: [b]cannot find existing dem paramters[/b]")

            # layers default
            ltmp = [0.002,0.004, 0.006,0.008,0.01,0.01,0.02,0.02,0.05,0.05,0.1,0.1,0.2,0.2,0.22]

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
            if keykwargs == "zratio":
                key = "zratio(i),i=1,nstr"
                if sum(value) != 1:
                    warnings.warn("sum is not equal to 1 -->" + str(sum(value)))
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
        '''
        Update zone file

        Parameters
        ----------
        zone_xyz : TYPE, optional
            DESCRIPTION. The default is [].

        Returns
        -------
        
        .. note:
            - create global zone_xyz variable

        .. note:
            - Updated dem_parameters file (number of zones).
            - Updated parm file.
            - Updated CATHYH file.

        '''

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

        self.zone_xyz = zone_xyz
        
        # update number of zone in the dem parameter file
        self.update_dem_parameters(nzone=len(np.unique(zone_xyz)))
        self.update_parm()
        self.update_cathyH(MAXZON=len(np.unique(zone_xyz)))

    # %% INPUT FILES

    def update_parm(self, verbose=False, **kwargs):
        '''
        Update parm file

        Parameters
        ----------
        verbose : TYPE, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        
        .. note:
            - create global parm variable

        .. note:
            - Updated CATHYH file (NPRT and NUMVP).

        '''

        self.console.print(":black_nib: [b]update parm file [/b]")

        if len(self.parm) == 0:
            
            # set default parameters
            # --------------------------------------------------------------------
            self.parm = {
                "IPRT1": 3,  # Flag for output of input and coordinate data
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
                "(NODVP(I),I=1,NUMVP)": [0],
                "NR": 0,
                "NUM_QOUT": 0,
                "(ID_QOUT(I),I=1,NUM_QOUT)": [441],
            }

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
                    self.parm[kk] = value
                    
                    # check if consistency between times of interest and TMAX
                    # ------------------------------------------------------------
                    if '(TIMPRT(I),I=1,NPRT)' in kk:
                        if self.parm["TMAX"] < max(self.parm['(TIMPRT(I),I=1,NPRT)']):
                            self.parm["TMAX"] = max(self.parm['(TIMPRT(I),I=1,NPRT)']) 
                    
                # check if consistency between DELTAT, DTMIN, and DTMAX
                # ------------------------------------------------------------
                if self.parm["DELTAT"] < self.parm["DTMIN"]:
                   warnings.warn('adjusting DTMIN == DELTAT')
                   self.parm["DTMIN"] = self.parm["DELTAT"]
               
                if self.parm["DELTAT"] > self.parm["DTMAX"]:
                   warnings.warn('adjusting DTMAX == DELTAT')
                   self.parm["DTMAX"] = self.parm["DELTAT"]       
                        
                    

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

    def update_ic(self, INDP=2, IPOND=0, WTPOSITION=0, verbose=False, **kwargs):
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
            the base of the 3‐d grid. The default is 0.

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
        
           
        with open(filename, "w+") as icfile:

            if INDP == 0:
                icfile.write(str(INDP) + "\t" + str(IPOND) +
                             "\t" + "INDP \n")
                if 'pressure_head_ini' in kwargs:
                    icfile.write(str(kwargs['pressure_head_ini']))
                    # self.update_mesh_vtk(prop='ic', prop_value=pressure_head_ini, savevtk=True)
                else:
                    raise ValueError('Missing initial pressure value')

            elif INDP == 1:
                pressure_head_ini = kwargs['pressure_head_ini']    
                icfile.write(str(INDP) + "\t" + str(IPOND) +
                             "\t" + "INDP" + "\t" + "IPOND" + "\n")
                np.savetxt(icfile, pressure_head_ini, fmt="%1.3e")
                self.update_mesh_vtk(prop='ic', prop_value=pressure_head_ini, savevtk=True)


            elif INDP in [2,3]:
                icfile.write(str(INDP) + "\t" + str(IPOND) +
                             "\t" + "INDP" + "\t" + "IPOND" + "\n")
                icfile.write(str(WTPOSITION) + "\t" + "WTPOSITION" + "\n")               

        icfile.close()
        
        
               

        pass

    def update_atmbc(self, HSPATM=0, IETO=0, TIME=None, VALUE=[None, None],
                     show=False, verbose=False, **kwargs):
        '''
        Atmospheric forcing term (atmbc - IIN6)


        ..note:
        
        
                1 1                HSPATM,IETO
                0.0000000e+00      TIME
                5.5e-06              VALUE
                12.000000e+03      TIME
                0.00                 VALUE
                18.000000e+03      TIME
                0.00                 VALUE
                
                The values are those of a 200-min rainfall event at a uniform 
                intensity of 3.3·10-4 m/min, followed by 100 min of drainage.

        ..note:
            
                In case of simultaneous precipitation and evaporation, we impose at 
                the surface the net flux, i.e., precipitation minus evaporation.





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
        
        ..note:
            
                - Update parm file (NPRT).
                - Update CATHYH file (MAXPRT).

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
                atmbcfile.write("{:.3e}".format(t) + "\t" + "TIME" + "\n")
                # atmbcfile.close()

                if isinstance(v, float) | isinstance(v, int):
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
            # x_units = "sec"
            # for key, value in kwargs.items():
            #     if key == "x_units":
            #         x_units = value
            plt, ax = plt_CT.show_atmbc(TIME, VALUE, 
                                        **kwargs)

            return plt, ax
        
        pass
    

              

    def init_boundary_conditions(self):
        '''
        .. note:
            The boundary conditions are defined in the nansfdirbc (Dirichlet), 
            nansfneubc (Neumann), and sfbc (seepage face) files. 

            We have two types of boundary conditions (BC): 
            - Neumann BC (or specifed flux)
            - Dirichlet BC (or pressure).
            
            
        .. note:
            - Pioggia: condizioni di Neumann. Quando non ci può più essere
            infiltrazione metto Dirichlet. 
            - Evaporazione: si indica un limite di pressione minimo ( Pmin ) al di 
            sotto del quale si ha uno switch da Neumann a Dirichlet
            (in quanto al di sotto di questo valore non si ha più evapotraspirazione).

        .. note:
            The boundary condition for any given surface node can switch between a
            Dirichlet condition and a Neumann condition depending on the saturation
            (or pressure) state of that node. 
            
        .. note:
            A Neumann (or specified flux) boundary condition corresponds to 
            atmosphere-controlled infiltration or exfiltration, with the flux equal 
            to the rainfall or potential evaporation rate given by the atmospheric input data. 
            When the surface node reaches a threshold level of saturation or moisture de cit, 
            the boundary condition is switched to a Dirichlet (specified head) condition, 
            and the infiltration or exfiltration process becomes soil limited [1].

        Returns
        -------
        None.

        '''
        # self.update_nansfdirbc()
        # self.update_nansfneubc()
        # self.update_sfbc()
        
        self.search_mesh_bounds(self.grid3d['nodes_idxyz'])
        plt_CT.plot_mesh_bounds(self.mesh_bound_cond_df)
        
        self.create_mesh_vtk()
        
        pass
        
        


    def update_nansfdirbc(self,TIME=[],NDIR=0,NDIRC=0,NQ3=None,noflow=True,
                            bound_xyz=None,
                            pressure_head=[],
                        ):
        '''
        Dirichlet Boundary conditions (or specified pressure) at TIME t
       
        - To simulate the no-flow boundaries conditions for the bottom and 
          vertical sides of the domain it is necessary to set NDIR and NDIRC 
          equal to zero. 
        - To simulate different boundary conditions, it is necessary to 
          indicate the number of selected nodes through NDIR or NDIRC, 
          then to specify the node ID’s that you want to consider and
          eventually the value of pressure head or flux that you want to assign.


        ..note :
            update_nansfdirbc use the grid3d to refer to mesh nodes
            
        Parameters
        ----------
        TIME : np.array([]), optional
            Absolute simulation time in sec. 
            The default is [].
        NDIR : int, optional
            Number of non-atmospheric, non‐seepage face Dirichlet
            nodes in 2-d mesh. The BC's assigned to these surface nodes are replicated vertically. 
            The default is 0.
        NDIRC : int, optional
            Number of 'fixed' non-atmospheric, non-seepage face Dirichlet
            nodes in 3‐d mesh ('fixed' in the sense that these BC's are not replicated to other nodes ‐
            compare NDIR). 
            The default is 0.
        NQ3 : int, optional
            Number of non-atmospheric, non‐seepage face Neumann nodes in 3‐d
            mesh. 
            The default is None.
        noflow : Bool, optional
            To simulate the no-flow boundaries conditions for the bottom and 
            # vertical sides of the domain. The default is True.
        bound_xyz : TYPE, optional
            DESCRIPTION. The default is None.
        pressure_head : TYPE, optional
            Specify a value of node pressure head to impose as Dirichlet boundary condition. 
            The default is [].
         : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        # check that the mesh exist
        # --------------------------------------------------------------------
        try:
            self.grid3d = in_CT.read_grid3d(self.project_name)
        except OSError:
            print("grid3d missing - need to run the processor with IPRT1=2 first")

        self.init_boundary_conditions()

        # read existing input nansfdirbc file
        # --------------------------------------------------------------------
        with open(os.path.join(self.workdir, 
                               self.project_name, 
                               self.input_dirname,
                               "nansfdirbc"
                               ),"w+") as nansfdirbcfile:


            # To simulate the no-flow boundaries conditions for the bottom and 
            # vertical sides of the domain --> NDIR and NDIRC equal to zero
            #-------------------------------------------------------------
            if noflow: # Dirichlet
                if len(TIME) == 0:
                    TIME = self.atmbc["TIME"]
                for tt in TIME:
                    nansfdirbcfile.write(str(tt) + "\t" + "TIME" + "\n")
                    nansfdirbcfile.write(str(NDIR) + "\t"+ str(NDIRC) + "\t" + "NDIR" + "\t" + "NDIRC" + "\n")
                    
                    
                self.update_parm()
                self.update_cathyH()
                
                bool_Dirichlet = []
                bool_Dirichlet_val = []                
                for id_node in self.mesh_bound_cond_df['id_node']:
                    if self.mesh_bound_cond_df['bound'].loc[int(id_node-1)] == True:
                        bool_Dirichlet.append('Dirichlet')
                        bool_Dirichlet_val.append(1)
                    else:
                        bool_Dirichlet.append('Neumann')
                        bool_Dirichlet_val.append(0)

                self.update_mesh_bounds(bound_type='bound_type',
                                        bound_bool=bool_Dirichlet)
                self.update_mesh_vtk('bound_type',bool_Dirichlet_val)
                
                # len(bool_Dirichlet)
                # len(self.mesh_bound_cond_df['id_node'])
                
                
            # elif bound_xyz is not None:
            #     for tt in TIME:
            #         #Se avessi delle variazioni dovrei indicare il nodo ed il valore di pressione
            #         nansfdirbcfile.write(str(tt) + "\t" + 'TIME' + "\n")
            #         for i in range(self.nnod3):
            #             if self.xmesh[i] == xb_left or self.xmesh[i] == xb_right or
            #                 self.ymesh[i] == yb_left or self.ymesh[i] == yb_right:
            #                     nansfdirbcfile.write(str(i) + "\n")

            #         for i in range(self.nnod3):
            #             if self.xmesh[i] == xb_left or self.xmesh[i] == xb_right or
            #                 self.ymesh[i] == yb_left or self.ymesh[i] == yb_right:
            #                     nansfdirbcfile.write(str(self.zmesh[i]-self.ic['WTPOSITION']) + "\n")

            # NPMAX = len(NDIRC)
            # NP2MAX = len(NDIR)
            #     self.update_cathyH(NPMAX, NP2MAX)


                # exemple provided by Laura B.
                # ----------------------------
    
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
        
        

        pass

    def update_nansfneubc(self, TIME=[], NQ=0, ZERO=0):
        '''
        Neumann boundary conditions (or specifed flux) at TIME t


        Parameters
        ----------
        TIME : np.array([]), optional
            Absolute simulation time in sec. 
            The default is [].
        NQ : TYPE, optional
            number of non-atmospheric, non seepage face Neumann nodes in
            3-D mesh. The default is 0.
        ZERO : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        None.

        '''
        # read existing input nansfneubc file
        # --------------------------------------------------------------------

        with open(os.path.join(self.workdir, 
                               self.project_name, 
                               self.input_dirname, 
                               "nansfneubc"),"w+") as nansfneubcfile:

            if len(TIME) == 0:
                TIME = self.atmbc["TIME"]
            for tt in TIME:
                nansfneubcfile.write(str(tt) + "\t" + "TIME" + "\n")
                nansfneubcfile.write(
                    str(ZERO) + "\t" + str(NQ) + "\t" +
                        "ZERO" + "\t" + "NQ" + "\n"
                )

        nansfneubcfile.close()

        pass

    def update_sfbc(self, TIME=[]):
        '''
        Seepage face boundary conditions at TIME t

        Parameters
        ----------
        TIME : np.array([]), optional
            Absolute simulation time in sec. 
            The default is [].

        Returns
        -------
        None.

        '''

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

        pass

    
    def search_mesh_bounds(self,grid3d):
        
        self.mesh_bound_cond_df = pd.DataFrame(grid3d)
        self.mesh_bound_cond_df = self.mesh_bound_cond_df.rename(
                                    columns={0: "id_node",
                                             1: "x", 
                                             2: "y", 
                                             3: "z"})
    
        minx = self.mesh_bound_cond_df['x'].min()
        maxx = self.mesh_bound_cond_df['x'].max()
        miny = self.mesh_bound_cond_df['y'].min()
        maxy = self.mesh_bound_cond_df['y'].max()
        minz = self.mesh_bound_cond_df['z'].min()
        maxz = self.mesh_bound_cond_df['z'].max()
        
        test = self.mesh_bound_cond_df['z'].to_numpy()
        
        # surface_nodes: surface node
        
        bound_idx_min =  self.mesh_bound_cond_df['id_node'][self.mesh_bound_cond_df['x']==minx]
        bound_idy_min =  self.mesh_bound_cond_df['id_node'][self.mesh_bound_cond_df['y']==miny]
        bound_idz_min =  self.mesh_bound_cond_df['id_node'][self.mesh_bound_cond_df['z']==minz]
        
        bound_idx_max =  self.mesh_bound_cond_df['id_node'][self.mesh_bound_cond_df['x']==maxx]
        bound_idy_max =  self.mesh_bound_cond_df['id_node'][self.mesh_bound_cond_df['y']==maxy]
        bound_idz_max =  self.mesh_bound_cond_df['id_node'][self.mesh_bound_cond_df['z']==maxz]
        
        
        bounds_id = np.unique(np.concatenate([bound_idx_min,bound_idy_min,bound_idz_min,
                        bound_idx_max,bound_idy_max,bound_idz_max]))
        
        bound_bool = []
        for idb in self.mesh_bound_cond_df['id_node']:
            if idb in bounds_id:
                bound_bool.append(True)
            else:
                bound_bool.append(False)
            
        self.mesh_bound_cond_df['bound'] = bound_bool
        
        pass
        

    def update_mesh_bounds(self, bound_type='bound_type', bound_bool=[]):
        '''
        update_mesh_bounds

        Parameters
        ----------
        bound_type : str
            Neumann or Dirichlet.
        bound_bool : TYPE
            Boolean for bound cond.
        '''
        self.mesh_bound_cond_df[bound_type] = bound_bool
        
        pass
            

    def create_mesh_vtk(self,topo=None):
        '''
        Create custum mesh

        Returns
        -------
        None.

        '''
        
        
        import numpy as np
        import pyvista as pv 
        
        if topo==None:
        
            si, sj, sk  = ( self.hapin["M"]*self.dem_parameters['delta_x'], 
                            self.hapin["N"]*self.dem_parameters['delta_y'], 
                            self.dem_parameters['base']
                            )
            
            ni, nj, nk = int(self.hapin["M"]), int(self.hapin["N"]), self.dem_parameters['nstr']
            
            xcorn = np.arange(0, (ni+1)*si, si)
            xcorn = np.repeat(xcorn, 2)
            xcorn = xcorn[1:-1]
            xcorn = np.tile(xcorn, 4*nj*nk)
            
            ycorn = np.arange(0, (nj+1)*sj, sj)
            ycorn = np.repeat(ycorn, 2)
            ycorn = ycorn[1:-1]
            ycorn = np.tile(ycorn, (2*ni, 2*nk))
            ycorn = np.transpose(ycorn)
            ycorn = ycorn.flatten()
            
            zcorn = np.arange(0, (nk+1)*sk, sk)
            zcorn = np.repeat(zcorn, 2)
            zcorn = zcorn[1:-1]
            zcorn = np.repeat(zcorn, (4*ni*nj))
            
            corners = np.stack((xcorn, ycorn, zcorn))
            corners = corners.transpose()
            
            dims = np.asarray((ni, nj, nk))+1
            self.mesh_pv_attributes = pv.ExplicitStructuredGrid(dims, corners)
            # self.mesh_pv = self.mesh_pv.compute_connectivity()
            # self.mesh_pv.plot(show_edges=True)
        
        else:
            # https://docs.pyvista.org/examples/00-load/create-structured-surface.html
            raise NotImplementedError('creation of mesh with topo not yet implemented')
    
        pass

        

    def update_mesh_vtk(self, prop='', prop_value=[], savevtk=True):
        '''
        https://docs.pyvista.org/api/core/_autosummary/pyvista.ExplicitStructuredGrid.add_field_data.html#pyvista.ExplicitStructuredGrid.add_field_data

        Parameters
        ----------
        prop : TYPE, optional
            DESCRIPTION. The default is ''.
        prop_value : TYPE, optional
            DESCRIPTION. The default is [].

        Returns
        -------
        None.

        '''
        
        self.mesh_pv_attributes.add_field_data(prop_value, prop)
        self.mesh_pv_attributes.save(os.path.join(self.workdir, 
                               self.project_name, 
                               'vtk/',
                               self.project_name +
                               '.vtk'))
        
        pass
            
            
            
                    
    def update_soil(self, IVGHU=[], FP=[], SPP=[], soil_het_dim=1,
                    verbose=False, show=False, **kwargs):
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
                
        heteregeneity_dim: int
            - dimension of the heteregeneity
            

        verbose : TYPE, optional
            DESCRIPTION. The default is False.
        **kwargs

        Returns
        -------
        - update parm file
        - update CATHY.H file
        - update mesh vtk file

        '''
                
        # set default parameters if SPP and/or FP args are not existing yet
        # --------------------------------------------------------------------
        
        if len(self.soil)==0:
            self.set_SOIL_defaults()
        if len(SPP) == 0:
            SPP = self.set_SOIL_defaults(SPP_default=True)
        if len(FP) == 0:
            FP = self.set_SOIL_defaults(FP_default=True)
            

        # check size of soil properties 
        # --------------------------------------------------------------------
            
        if isinstance(SPP['PERMX'], float):
            if self.dem_parameters["nzone"]!=1:
                raise ValueError
        else:
            if len(SPP['PERMX'])!=self.dem_parameters["nzone"]:
                raise ValueError("Wrong number of zones: PERMX size is" + str(len(SPP['PERMX'])) 
                                 + 'while nzone is' + str(self.dem_parameters["nzone"]))
            
            

            
        # read function arguments kwargs and udpate soil and parm files
        # --------------------------------------------------------------------
        for keykwargs, value in kwargs.items():
            if verbose == True:
                print(f"new: {keykwargs} | value: {value}")
            self.soil[keykwargs] = value
            self.parm[keykwargs] = value

        # loop over Feddes parameters FP mandatory fct arg
        # --------------------------------------------------------------------
        for fp in FP:  # loop over fedded parameterssoil_het_dim
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
        self.soil_SPP = {}
        self.soil_SPP['SPP'] = SoilPhysProp # matrice with respect to zones
        self.soil_SPP['SPP_map'] = SPP # mapping with respect to zones

        # Vegetation properties (PCANA,PCREF,PCWLT,ZROOT,PZ,OMGC)
        # --------------------------------------------------------------------
        # Read only if IVGHU=0
        
        FeddesParam = self._prepare_SOIL_vegetation_tb(FP)
        self.soil_FP = {}
        self.soil_FP['FP'] = FeddesParam
        self.soil_FP['FP_map'] = FP # mapping with respect to zones

        # write soil file
        # --------------------------------------------------------------------
        self._write_SOIL_file(SoilPhysProp, FeddesParam)
        
        if show:
            raise NotImplementedError('show for soil not yet implemented')
            
        #     SPP_zones = self.map_prop2zone(SPP,prop=show)
        #     plt, ax = plt_CT.dem_plot_2d_top(SPP_zones,
        #                                      label=show)
            
        #     FP_zones = self.map_prop2zone(FP,prop='ZROOT')
        #     plt, ax = plt_CT.dem_plot_2d_top(SPP_zones,
        #                                      label=show)

        #     return plt, ax        
        
        
        # map SPP to the mesh
        # --------------------------------------------------------------------
        self.map_prop2mesh(SPP)


        
        
            
        pass
    
    

        
        
        

    def set_SOIL_defaults(self,
                          FP_default=False, 
                          SPP_default=False):

        self.soil = {
            "PMIN": -5.0,
            "IPEAT": 0,
            "SCF": 1.0, # here we assume that all soil is covered by the vegetation
            "CBETA0": 0.4,
            "CANG": 0.225,
            # Feddes parameters default values
            "PCANA": [0.0],
            "PCREF": [-4.0],
            "PCWLT": [-150],
            "ZROOT": [0.1],
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


        if FP_default:
            
            FP = {
                    # Feddes parameters default values
                    "PCANA": [0.0],
                    "PCREF": [-4.0],
                    "PCWLT": [-150],
                    "ZROOT": [0.1],
                    "PZ": [1.0],
                    "OMGC": [1.0],
                }
                
            # self.soil.update(FP)
            
            return FP

            
        # # set Soil Physical Properties defaults parameters
        # # --------------------------------------------------------------------
            
        if SPP_default:
            
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
            
            return SPP


        pass 

    def _prepare_SPP_tb(self, SPP, **kwargs):
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
        
        # heteregeneity_dim = 1
        # if 'soil_het_dim' in kwargs:
        #     heteregeneity_dim = kwargs['soil_het_dim']
            
            
        # check number of zones
        if self.dem_parameters["nzone"] > 1:
            
            
            if len(SPP['PERMX'])<=1:
                for i, spp in enumerate(SPP):
                    SPP[spp] = SPP[spp]*np.ones(self.dem_parameters["nzone"])
            else:
                pass
                    

                
            
        # check size of the heteregeneity of SPP
        # ----------------------------------------
        # 1d --> uniform
        # 2d --> lateral variations due to zones defined in surface
        # 3d --> lateral + vertical variations due to zones and strates
        
        # if heteregeneity_dim > 1:
                       
            # if (heteregeneity_dim==2|heteregeneity_dim==3):
                # raise NotImplementedError
            SoilPhysProp = []


            # loop over strates
            # -----------------------------------------------------------------
            for istr in range(self.dem_parameters["nstr"]):
                    
                #  loop over zones (defined in the zone file)
                # --------------------------------------------------------------
                for izone in range(self.dem_parameters["nzone"]):
                    izoneSoil = np.zeros([self.dem_parameters["nstr"], 8])
                   
                    for i, spp in enumerate(SPP):
                        izoneSoil[:,i] = np.ones(self.dem_parameters["nstr"])*SPP[spp][izone]
                        
                SoilPhysProp.append(izoneSoil)
            SoilPhysProp = np.vstack(SoilPhysProp)


            # SoilPhysProp_zone = np.ones([self.dem_parameters["nstr"],8])
            # SoilPhysProp_zone = np.ones([8,self.dem_parameters["nstr"]])
            # SoilPhysProp = []
            # #  loop over zones (defined in the zone file)
            # # --------------------------------------------------------------
            # for izone in range(self.dem_parameters["nzone"]):
            #     # loop over strates
            #     # -----------------------------------------------------------------
            #     for istr in range(self.dem_parameters["nstr"]):
                    
            #         for i, spp in enumerate(SPP):
            #             SoilPhysProp_zone[i,:] = np.ones(self.dem_parameters["nstr"])*SPP[spp][izone]
            
            #     SoilPhysProp.append(SoilPhysProp_zone)
            
            # SoilPhysProp = np.reshape(SoilPhysProp,[8*self.dem_parameters["nstr"],self.dem_parameters["nzone"]])
                    
            # np.shape(SoilPhysProp)
                #     izoneSoil_tmp = []
                #     for spp in SPP:
                #         izoneSoil_tmp.append(SPP[spp][izone])
                #     izoneSoil_tmp = np.hstack(izoneSoil_tmp)
                #     izoneSoil[izone, :] = izoneSoil_tmp
                #     ki = k
                #     ke = k + self.dem_parameters["nzone"]
                #     SoilPhysProp[ki:ke, :] = izoneSoil
                #     # SoilPhysProp[self.dem_parameters['nzone']*istr*izone:self.dem_parameters['nzone']*istr*izone+self.dem_parameters['nzone'],:]=izoneSoil
                # k += self.dem_parameters["nzone"]
                
                
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
        
        # Check if root_map file exist and is updated
        # -------------------------------------------
        if hasattr(self,'veg_map') is False:
            if len(FP['PCANA'])==1:
                self.update_veg_map()
            else:
                raise ValueError('Found multiple values of Feddes zones' +
                                 'but vegetation map is not defined')
       
        
        # Check vegetation heterogeneity dimension
        # ----------------------------------------
        if self.cathyH["MAXVEG"] != len(FP['PCANA']):
            raise ValueError("Wrong number of vegetations: PCANA size is " + str(len(FP['PCANA'])) 
                             + 'while MAXVEG is ' + str(self.cathyH["MAXVEG"]))

        # check number of vegetation
        # --------------------------------------------------------------------
        if self.cathyH["MAXVEG"] > 1:
            print(self.cathyH["MAXVEG"])
            FeddesParam = np.zeros([self.cathyH["MAXVEG"], 6])
            for iveg in range(self.cathyH["MAXVEG"]):  # loop over veg zones within a strate
                izoneVeg_tmp = []
                for sfp in FP:
                    izoneVeg_tmp.append(FP[sfp][iveg])

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

        # self.update_zone(indice_veg)
        self.update_cathyH(MAXVEG=len(np.unique(indice_veg)))
        # ,                        MAXZON=len(np.unique(indice_veg))
        
        self.veg_map = indice_veg


        if show is not None:
            plt, ax = plt_CT.indice_veg_plot(indice_veg)
            
            return indice_veg, plt, ax
        

        return indice_veg

    # -------------------------------------------------------------------#
    # %% DATA ASSIMILATION FCTS

    def _DA_init(self, NENS=[], ENS_times=[], parm_pert=[]):
        """
        THIS SHOULD BE MOVED TO DA CLASS

        Initial preparation for DA


        .. note::

            1. update cathyH file given NENS, ENS_times + flag self.parm['DAFLAG']==True

            2. create the processor cathy_origin.exe

            3. create directory tree for the ensemble

            4. create panda dataframe for each perturbated variable

            5. update input files

        Parameters
        ----------
        NENS : int
            # Number of realizations in the ensemble
        ENS_times : int
            # Observation times

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
        # NUDN = len(ENS_Times)  # this is true only when DA is made inside CATHY
        # self.update_cathyH(MAXNUDN=NUDN,ENS_Times=ENS_Times, verbose=True)

        # when DA is made outside CATHY we run step by step

        # TIMPRTi # time of interest for outputs
        # TIMPRTi=times_of_interest,
        # NODVP=nodes_of_interest,

        # THIS IS TEMPORARY, assimilation time should be infer from observation data directly
        # self.ENS_Times = ENS_Times

        # ?? is this necessary ??
        # self.update_cathyH(MAXNUDN=1,ENS_Times=ENS_Times, verbose=True)

        # run processor to create the cathy_origin.exe to paste in every folder
        # ---------------------------------------------------------------------
    
        # self.run_processor(runProcess=False, recompile=False)
        # OR self.recompileSrc(runProcess=False, NUDN=NUDN)

        # create sub directories for each ensemble
        # ---------------------------------------------------------------------
        self._create_subfolders_ensemble(NENS)

        # create initial dataframe DA_results_df for results
        # ---------------------------------------------------------------------
        # X, N_col, M_rows = self._read_state_ensemble()      
        # self._DA_df()


        # update all for the initial state to form the ensemble

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
        print(' ---- check scenarii before analysis ----')

        

        
    def _check_after_analysis(self,update_key,list_update_parm):
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        CHECK which scenarios are OK and which one to discard

        Returns
        -------
        None.

        '''
        print(' ---- check scenarii post update ----')

        for pp in enumerate(list_update_parm[:]):
            if 'St. var.' in pp[1]:
                pass
            else:
                # self.dict_parm_pert[pp[1]][update_key] = AnalysisParam

                # ckeck if new value of Zroot is feasible
                if not (self.dict_parm_pert[pp[1]][update_key] > abs(min(self.grid3d['nodes_idxyz'][:, -1]))).any:
                    print('impossible value of zroot')
                    print(self.dict_parm_pert[pp[1]][update_key])
                    print(min(self.grid3d['nodes_idxyz'][:, -1]))
                    sys.exit()
                    # self.dict_parm_pert[pp[1]][update_key] = self.dict_parm_pert[pp[1]]['ini_perturbation']-0.01*self.count_DA_cycle

            # print(self.dict_parm_pert)

        # condition 1: nrow(mbeconv)==0
        # condition 2: mbeconv[h,3]==(deltaT)

        pass

    def _DA_analysis(self, prediction, 
                     typ='enkf_Evensen2009_Sakov',
                     list_update_parm=['St. var.'],
                     list_assimilated_obs=[]):
        """
        THIS SHOULD BE MOVED TO DA CLASS

        Analysis ensemble using DA


        1. map state variable 2 Observations

        2. apply filter

        3. checkings

        Parameters
        ----------
        typ : str
            type of Data Assimilation
        NUDN : int
            #  NUDN  = number of observation points for nudging or EnKF (NUDN=0 for no nudging)
        ENS_times : np.array([])
            # Observation time for ensemble kalman filter (EnKF).

        Returns
        -------
        None.

        Note: for now only implemented for EnkF DA
        """

        update_key = 'ini_perturbation'
        if self.count_DA_cycle > 0:
            update_key = 'update_nb' + str(self.count_DA_cycle)


        # prepare parameters
        # ---------------------------------------------------------------------
        # When updating the states only, the elements of X are the
        # pressure heads at each node of the finite element grid, while
        # the state augmentation technique is used when also updat-
        # ing the parameters

        print(list_update_parm)
        param = []
        for pp in enumerate(list_update_parm[:]):
            if 'St. var.' in pp[1]:
                pass 
            else:               
            # (EnSize), self.dict_parm_pert
                param = self.dict_parm_pert[pp[1]][update_key]
                
                if self.dict_parm_pert[pp[1]]['transf_type'] == 'Log':
                    print('log transformation')
                    param = np.log(self.dict_parm_pert[pp[1]][update_key])
 
     
        # prepare states
        # ---------------------------------------------------------------------
        ensembleX, ens_size, sim_size = self._read_state_ensemble()

        # self.ensembleX_test = ensembleX
        # print(ensembleX)
        # print(self.count_DA_cycle)
        # print(max(np.cov(self.ensembleX_test)[0]))
        
        # plt.plot(ensembleX[:,0])
        # plt.show()

        # input("Please once you checked the ensemble (before analysis)")

        # ---------------------------------------------------------------------
        # prepare data
        # ---------------------------------------------------------------------

        def get_data(dict_obs,list_assimilated_obs, assimilation_time):
            # Loop trought observation dictionnary for a given assimilation time (t)
            # -----------------------------------------------------------------
            data = []    # data dict 2 map
            if list_assimilated_obs == 'all':
                items_dict = list(self.dict_obs.items())
                # There can be multiple sensors measuring at the same time
                # -----------------------------------------------------------------
                for sensors in items_dict[assimilation_time][1].keys():                
                    data.append(items_dict[assimilation_time][1][sensors]['data'])
                            
            return data
            
        # Loop trought observation dictionnary for a given assimilation time (t)
        # -----------------------------------------------------------------
        data = []    # data dict 2 map
        if list_assimilated_obs == 'all':
            items_dict = list(self.dict_obs.items())
            # There can be multiple sensors measuring at the same time
            # -----------------------------------------------------------------
            for sensors in items_dict[self.count_DA_cycle][1].keys():                
                data.append(items_dict[self.count_DA_cycle][1][sensors]['data'])
          
        
        # check size of data_cov
        # ---------------------------------------------------------------------
        if len(self.stacked_data_cov) != len(data):
            raise ValueError('need to compute data covariance')
            
        # run Analysis
        # ---------------------------------------------------------------------
                
        result_analysis = cathy_DA.run_analysis(typ,
                                                data,
                                                self.stacked_data_cov,
                                                param,
                                                ensembleX,
                                                prediction)
          
        # plot ensemble covariance matrices and changes (only for ENKF)
        # ---------------------------------------------------------------------
        if len(result_analysis)>2:
            [A, Amean, dA, 
             dD, MeasAvg, S, 
             COV, B, dAS, 
             analysis, 
             analysis_param]  = result_analysis
        
            plt_CT.show_DA_process_ens(ensembleX,data,self.stacked_data_cov,dD,dAS,B,analysis, 
                                       savefig=True, 
                                       savename= os.path.join(self.workdir,self.project_name,'DA_Matrices_t' 
                                       +str(self.count_DA_cycle)))
        # particule filter
        # ----------------
        else:
            [analysis, 
             analysis_param]  = result_analysis
            

        # Back-transformation of the parameters
        # ---------------------------------------------------------------------
        for pp in enumerate(list_update_parm[:]):
            if 'St. var.' in pp[1]:
                pass 
            else:                               
                if self.dict_parm_pert[pp[1]]['transf_type'] == 'Log':
                    print('back log transformation')
                    analysis_param = np.exp(analysis_param)
                    # self.dict_parm_pert[pp[1]][update_key]=analysis_param
                    print(analysis_param)

                    
                    
        # return EnsembleX, Analysis, AnalysisParam, Data, data_cov
        # ---------------------------------------------------------------------
        return(ensembleX, data,  analysis, 
                 analysis_param)



                            
                            

    def _map_states2Observations(self, state=[None, None],
                                 list_assimilated_obs='all', 
                                 parallel=False,
                                 verbose = False,
                                 **kwargs):
        '''
        Parameters
        ----------
        state : np.array([])
            [pressure heads, saturation] for each mesh nodes. The default is [None, None].
        list_assimilated_obs : TYPE, optional
            DESCRIPTION. The default is 'all'.
        parallel : Bool, optional
            DESCRIPTION. The default is False.
        verbose : Bool, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.
            
            
        THIS SHOULD BE MOVED TO DA CLASS

        Apply H mapping operator: convert model to predicted value
        
        
        .. note:
        
            - For ERT data, H is Archie law, SWC --> ER0 --> ERapp
            - For SWC data, no mapping needed (only calibration)
            - For tensiometer data, no mapping needed 
            - For discharge, no mapping needed

        Returns
        -------
        Hx : np.array
            Ensemble of the simulated observations.
        '''
        

        
        # infer soil parameters properties
        # ---------------------------------
        porosity = self.soil_SPP['SPP'][:, 4][0]


        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------
        obskey2map = [] # data type 2 map
        obs2map = []    # data dict 2 map
        if list_assimilated_obs == 'all':
            items_dict = list(self.dict_obs.items())
            # There can be multiple sensors measuring at the same time
            # -----------------------------------------------------------------
            for sensors in items_dict[self.count_DA_cycle][1].keys():                
                obskey2map.append(sensors)
                obs2map.append(items_dict[self.count_DA_cycle][1][sensors])


        Hx= [] # Ensemble of the simulated observations
        # Loop over observations to map
        # ---------------------------------------------------------------------
        for i, obs_key in enumerate(obskey2map):
            print(obs_key)
            
            if 'PH' in obs_key:
                # case 1: pressure head assimilation (Hx_PH)
                # -------------------------------------------------------------
                # df_vp_PH = df_vp.table_pivot(
                #     index=[time, node], value='PH')
                # # if key[0] in 'PH':
                # # filter the node if some of the instruments were not working
                # df_vp_PH_filt = df_vp_PH.iloc('node1', 'node2')

                 Hx_PH = state[0][obs2map[i]['mesh_nodes']]
                 Hx.append(Hx_PH)
                 print('Not yet tested')


            if 'sw' in obs_key:
                # case 2: sw assimilation (Hx_SW)
                # --------------------------------------------------------------------
                Hx_SW = state[1][obs2map[i]['mesh_nodes']] * porosity
                Hx.append(Hx_SW)
                # note: the value of the porosity can be unique or not depending on the soil physical properties defined

            if 'discharge' in obs_key:
                # case 3: discharge
                # need to read the hgsfdet file (Hx_Q)
                # --------------------------------------------------------------------
                # if key[0] in 'discharge':

                # derivation of the dircharge, Q from file 'hgsfdet'
                # Hx_Q = []
                # Hx.vstack(Hx_Q)
                 print('Not yet implemented')

        

            if 'ERT' in obs_key:
                
                path_fwd_CATHY = ''
                if 'path_fwd_CATHY' in kwargs:
                   path_fwd_CATHY = kwargs['path_fwd_CATHY'] 

                ens_nb = []
                if 'ens_nb' in kwargs:
                   ens_nb = kwargs['ens_nb'] 
                   
                   
                savefig = True
                if 'savefig' in kwargs:
                   savefig = kwargs['savefig'] 
                
                
                # search key value to identify time and method
                tuple_list_obs = list(self.dict_obs.items())
                key_value = tuple_list_obs[self.count_DA_cycle]
                
                
                # Load ERT metadata information form obs dict
                # -------------------------------------------
                forward_mesh_vtk_file = key_value[1]['forward_mesh_vtk_file'] 
                pathERT = os.path.split(self.dict_obs[0]['filename'])[0]
                seq = key_value[1]['sequenceERT'] 
                electrodes = key_value[1]['elecs'] 

                            # project_name,
                            # porosity,
                            # pathERT, meshERT, elecs, sequenceERT,
                            # path_fwd_CATHY,
                            
                            
                Hx_ERT, df_Archie = Archie.SW_2_ERa(#state[1], 
                                                      path_fwd_CATHY, 
                                                      self.project_name,
                                                      porosity,
                                                      pathERT,
                                                      forward_mesh_vtk_file,
                                                      electrodes,
                                                      seq,
                                                      savefig=savefig,
                                                      DA_cnb=self.count_DA_cycle,
                                                      Ens_nb=ens_nb,
                                                      data_format= self.dict_obs[0]['data_format'],
                                                      df_sw = state[1]
                                                    )
                
                if 'pygimli' in self.dict_obs[0]['data_format']:
                    Hx.append(Hx_ERT['rhoa'])
                else:
                    Hx.append(Hx_ERT['resist'])


        return np.vstack(Hx) 
    
    
    def _add_2_ensemble_Hx(self, Hx, Hx_2add):

        try:
            Hx.append(Hx_2add)
        except:
            Hx = list(Hx)
            Hx.append(Hx_2add)


        return Hx


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
            
            print(j)
            print(os.path.join(self.workdir, self.project_name,
                                               "DA_Ensemble/cathy_" + str(j+1),
                                               'output/psi'))
            df_psi = out_CT.read_psi(os.path.join(self.workdir, self.project_name,
                                               "DA_Ensemble/cathy_" + str(j+1),
                                               'output/psi'))
            X[:, j] = df_psi[-1, :]
            
            # input('press again here')

        # check if there is still zeros
        if np.count_nonzero(X) != M_rows*N_col:
            print('unconsistent filled X, missing values')
            sys.exit()
            
        return X, N_col, M_rows 



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

        # if not os.path.exists(os.path.join(self.workdir, self.project_name,
        #                                    "DA_Ensemble/cathy_origin")):
        if os.path.exists(os.path.join(self.workdir, self.project_name,
                                 "DA_Ensemble/cathy_origin")):
            shutil.rmtree(os.path.join(self.workdir, self.project_name,
                                     "DA_Ensemble/cathy_origin"))           # Removes all the subdirectories!
        os.makedirs(os.path.join(self.workdir, self.project_name,
                                 "DA_Ensemble/cathy_origin"),exist_ok=True)

        # copy input, output and vtk dir
        # ----------------------------------------------------------------
        for dir2copy in enumerate(['input', 'output', 'vtk']):
            
            if os.path.exists(os.path.join(self.workdir, self.project_name,
                                           "DA_Ensemble/cathy_origin", dir2copy[1])
                              ):
                print('delete' + dir2copy)
                shutil.rmtree(os.path.join(self.workdir, self.project_name,
                                               "DA_Ensemble/cathy_origin", dir2copy[1]))     
                
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

            if os.path.exists(path_nudn_i):
                shutil.rmtree(path_nudn_i)  
                
            shutil.copytree(
                    path_origin, path_nudn_i
                )

    def _DA_df(self, state=None, state_analysis=None, **kwargs):
        '''
        THIS SHOULD BE MOVED TO DA CLASS
        
        Build a dataframe container to keep track of the changes before and after update for each ensemble.
        This dataframe is the support for DA results visualisation.


        Parameters
        ----------
        sw : np.array([])
            State variable mesh nodes values (before update).
        sw_analysis : np.array([])
            State variable mesh nodes values after DA analysis.

        Returns
        -------
        None.

        '''
        
        df_DA_ti_ni = pd.DataFrame() # nested df for a given ensemble/given time
        df_DA_ti = pd.DataFrame() # nested df for a given time
        
        cols_root = ['time', 
                     'Ensemble_nb',
                     'bef_update_', 
                     'aft_update_',
                     'OL'] # Open loop boolean


        # Open loop kwargs
        # --------------------------
        t_ass = []
        if 't_ass' in kwargs:
            t_ass = kwargs['t_ass']

        NENS = []
        if 'ens_nb' in kwargs:
            ens_nb = kwargs['ens_nb']
        # --------------------------


        if 'openLoop' in kwargs:

            data_df_root = [t_ass*np.ones(len(state)), 
                            ens_nb*np.ones(len(state)), 
                            state,
                            state,
                            True*np.ones(len(state))]
            df_DA_ti = pd.DataFrame(np.transpose(data_df_root),
                                  columns=cols_root)
                
        else:
 
            for n in range(self.NENS):               
                
                data_df_root = [self.count_DA_cycle*np.ones(len(state)), 
                                n*np.ones(len(state)), 
                                state[:,n],
                                state_analysis[:,n],
                                False*np.ones(len(state))]
         
                df_DA_ti_ni = pd.DataFrame(np.transpose(data_df_root),
                                      columns=cols_root)
        
                df_DA_ti= pd.concat([df_DA_ti, df_DA_ti_ni], axis=0)   

                       
        self.df_DA= pd.concat([self.df_DA, df_DA_ti], axis=0, ignore_index=True)


        pass




    def _update_input_ensemble(self, NENS, ENS_times,
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
        self.selec_atmbc_window(NENS, ENS_times)

        # ---------------------------------------------------------------------
        
        analysis = []
        if 'analysis' in kwargs:
            analysis = kwargs['analysis']
            
        self.update_ENS_files(parm_pert, 
                              update_parm_list=update_parm_list,
                              cycle_nb=self.count_DA_cycle,
                              analysis=analysis)

        pass

    def selec_atmbc_window(self, NENS, ENS_times):
        '''
        THIS SHOULD BE MOVED TO DA CLASS


        Select the time window of the hietograph
        == time between two assimilation observation

        Parameters
        ----------
        NENS : TYPE
            DESCRIPTION.
        ENS_times : TYPE
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
        df_atmbc, HSPATM, IETO = in_CT.read_atmbc(os.path.join(self.workdir, self.project_name,
                              'input', 'atmbc'), grid=self.grid3d)
        
        if self.count_DA_cycle is not None:
            try:
                time_window_atmbc = [
                    ENS_times[self.count_DA_cycle], ENS_times[self.count_DA_cycle+1]]
            except:
                pass
        else:
                time_window_atmbc = [ENS_times[0], ENS_times[1]]



        df_atmbc_window = df_atmbc[(df_atmbc.index >= time_window_atmbc[0]) &
                      (df_atmbc.index <= time_window_atmbc[1])]


        for ens_nb in range(NENS):

            # change directory according to ensmble file nb
            os.chdir(os.path.join(self.workdir, self.project_name,
                                  './DA_Ensemble/cathy_' + str(ens_nb+1)))
            
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

        if 'analysis' in kwargs:
            analysis = kwargs['analysis']


        if update_parm_list=='all':
            update_parm_list = ['St. var.']
            for pp in parm_pert:
                update_parm_list.append(pp)

                print(update_parm_list)
                
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
                                       pressure_head_ini=analysis[:,ens_nb],
                                       filename=os.path.join(os.getcwd(), 'input/ic'),
                                       backup=True)
            else:
                raise ValueError('no state variable update - use last iteration as initial conditions ?')
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

                # ic update (i.e. water table position update)
                # --------------------------------------------------------------
                if key in 'ic':                  
                    self.update_ic(INDP=3, IPOND=0,
                                   WTPOSITION=parm_pert[key
                                       ][update_key][ens_nb],
                                   verbose=True,
                                   filename=os.path.join(os.getcwd(), 'input/ic'))  # specify filename path

                # kss update
                # --------------------------------------------------------------
                if key in 'ks':
                    # print('Ks perturbated not yet implemented')
                    
                    # SPP = self.set_SOIL_defaults(SPP_default=True)
                    # check soil heterogeneity dimension 
                    # -----------------------------------
                    
                    # def check_soil_heteregeneity_dim(self):
                    #     dim_soil_het = []
                    #     return dim_soil_het
                    # dim_soil_het = self.check_soil_heteregeneity_dim()
                    
                    # check soil heterogeneity dimension 
                    # -----------------------------------
                    dim_soil_het = int(len(parm_pert[key][update_key])/self.NENS)
                        
                    SPP =  self.soil_SPP['SPP_map'] 
                    SPP.update(PERMX=parm_pert[key][update_key][ens_nb],
                               PERMY=parm_pert[key][update_key][ens_nb],
                               PERMZ=parm_pert[key][update_key][ens_nb])
                    
                    
                    # FeddesParam = {'PCANA': self.soil["PCANA"],
                    #                 'PCREF': self.soil["PCREF"],
                    #                 'PCWLT': self.soil["PCWLT"],
                    #                 'ZROOT': [parm_pert[key][update_key][ens_nb]],
                    #                 'PZ': self.soil["PZ"],
                    #                 'OMGC': self.soil["OMGC"]}

                    self.update_soil(SPP=SPP,
                                     FP=self.soil_FP['FP_map'],
                                     verbose=True,
                                     filename=os.path.join(os.getcwd(), 'input/soil')
                                     )  # specify filename path

                # FeddesParam update
                # --------------------------------------------------------------
                elif key in ['PCANA', 'PCREF', 'PCWLT', 'ZROOT', 'PZ', 'OMGC']:                   
                    
                    FeddesParam = {'PCANA': self.soil["PCANA"],
                                   'PCREF': self.soil["PCREF"],
                                   'PCWLT': self.soil["PCWLT"],
                                   'ZROOT': [parm_pert[key][update_key][ens_nb]],
                                   'PZ': self.soil["PZ"],
                                   'OMGC': self.soil["OMGC"]}
                    

                    self.update_soil(FP=FeddesParam,
                                     verbose=False,
                                     filename=os.path.join(os.getcwd(), 'input/soil'))  # specify filename path
                    

        pass

    def _performance_assessement(self,list_assimilated_obs, 
                                 data, 
                                 prediction, 
                                 t_obs,
                                 **kwargs):
        '''
        THIS SHOULD BE MOVED TO DA CLASS
        (Normalized) root mean square errors (NRMSEs) 
        RMSE is compute separately for each observation assimilated 


        Parameters
        ----------
        list_assimilated_obs : TYPE
            DESCRIPTION.
        Data : TYPE
            Refers to the measured data.
        Observation : TYPE
            Refers to the simulated data.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''


        OL_bool = False
        if 'openLoop' in kwargs:          
            OL_bool = True
            
            
        # NAs default     
        # -------------------------------------------------------------
        RMSE = np.nan # root mean square error (RMSE)
        NMRMSE = np.nan # time-averaged normalized root mean square error
        
       
        
        # # compute metrics for each observation variable
        # # ------------------------------------------
        
        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------
        obs2eval_key = [] # data type 2 eval
        obs2eval = []    # data dict 2 eval
        if list_assimilated_obs == 'all':
            items_dict = list(self.dict_obs.items())
            # There can be multiple sensors measuring at the same time
            # -----------------------------------------------------------------
            for sensors in items_dict[self.count_DA_cycle][1].keys():                
                obs2eval_key.append(sensors)
                obs2eval.append(items_dict[self.count_DA_cycle][1][sensors])


        # if list_assimilated_obs == 'all':
        #     items_dict = list(self.dict_obs.items())
        #     obs2map = [items_dict[self.count_DA_cycle][1]['data_type']]
            
        for s, name_sensor in enumerate(obs2eval_key): # case where there is other observations than ERT
            print(name_sensor)
  
            # number of observation at a given time
            # for ERT number of observation = is the number of grid cells
            n_obs = np.shape(prediction)[1]
    
            data_Obs_diff_mat = np.zeros(np.shape(prediction))
            for i in range(len(data)):
                for j in range(np.shape(prediction)[0]): # Loop over ensemble rows
                    data_Obs_diff_mat[j,i] = abs(data[i]-prediction[j,i])
    
            # average differences over the number of ensemble
            data_Obs_diff_avg = np.sum(data_Obs_diff_mat,axis=0)*(1/n_obs)
            
            # average differences over the observations
            RMSE_avg = np.sum(data_Obs_diff_avg,axis=0)*(1/len(data))
            RMSE_sensor = data_Obs_diff_avg[s]
    
    
    
            try:
                self.df_performance
            except:
                # initiate if simulation just started             
                # -------------------------------------------------------------
                self.df_performance = pd.DataFrame()   
                self.RMSE = []
            
    
            self.RMSE.append(RMSE) 
    
    
    
            # root names for the collumns name                
            # -------------------------------------------------------------
            cols_root = ['time', 
                         'ObsType',
                         'RMSE',
                         'RMSE_avg',
                         'OL']        
            
            # root data                
            # -------------------------------------------------------------
            data_df_root = [[t_obs], 
                            [name_sensor],
                            [RMSE_sensor],
                            [RMSE_avg],
                            [OL_bool]]
            
            
            df_performance_ti = pd.DataFrame(np.array(data_df_root).T,
                                  columns=cols_root)
            
    
            # concatenate with main RMSE dataframe
            # -------------------------------------------------------------
            self.df_performance= pd.concat([self.df_performance, df_performance_ti], 
                                       axis=0, 
                                       ignore_index=True)
            
            
            # add NMRSE when all assimilation times have been integrated
            # ----------------------------------------------
            NMRMSE = (1/(t_obs+1))*np.sum(self.RMSE)
            # self.df_performance['NMRMSE'] = NMRMSE
            

        return self.df_performance
    
    

    def read_observations(self, filename, data_type, data_err, show=False, **kwargs):
        '''
        read measures (real observations) from file 
        and prepare for DA analysis (covariance matrice and perturbation)
        Uses `pyCATHY.importers`, `sensors_measures` to read files with common standards
        Need to be call for each time/ each observation

        Parameters
        ----------
        filename : str or dataframe
            filename (or data) of the observation dataset.
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

        
        OrderedDict([(0.0, 
                      {'swc': {'filename': None, 
                               'data_type': 'swc', 
                               'units': '$m^{3}/m^{3}$', 
                               'data': 0.3391464625, 
                               'data_err': 0.01, 
                               'mesh_nodes': [159], 
                               'assimilation_times': 0.0, 
                               'data_cov': [], 
                               'dataPert': [], 
                               'sensor_name': 'swc'}, 
                       'swc1': {'filename': None, 
                                'data_type': 'swc', 
                                'units': '$m^{3}/m^{3}$', 
                                'data': 0.3391464625, 
                                'data_err': 0.01, 
                                'mesh_nodes': [159], 
                                'assimilation_times': 0.0, 
                                'data_cov': [], 
                                'dataPert': [], 
                                'sensor_name': 'swc'}, 
                       }
                     3600,
                      {'swc': {'filename': None, 
                               'data_type': 'swc', 
                               'units': '$m^{3}/m^{3}$', 
                               'data': 0.3391464625, 
                               'data_err': 0.01, 
                               'mesh_nodes': [159], 
                               'assimilation_times': 0.0, 
                               'data_cov': [], 
                               'dataPert': [], 
                               'sensor_name': 'swc'}, 
                       'swc1': {'filename': None, 
                                'data_type': 'swc', 
                                'units': '$m^{3}/m^{3}$', 
                                'data': 0.3391464625, 
                                'data_err': 0.01, 
                                'mesh_nodes': [159], 
                                'assimilation_times': 0.0, 
                                'data_cov': [], 
                                'dataPert': [], 
                                'sensor_name': 'swc'}, 
                       .
                       .
                       .
                     )])
        
        '''

        dict_obs_2add = {}
        # specify mesh node position for point measurements
        # ---------------------------------------------------------------------
        mesh_nodes = []
        if 'mesh_nodes' in kwargs:
            mesh_nodes = kwargs['mesh_nodes']
            
            
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


        # discharge type read
        # ---------------------------------------------------------------------
        if data_type == 'swc':
            if isinstance(filename, str):
                df = in_meas.read_swc(filename)
            else:
                df = filename
                filename = None
            units = '$m^{3}/m^{3}$'
            obs_cov_type = None 
            # point sensors --> covariance between the sensors is formed later

        # tensiometer type read
        # ---------------------------------------------------------------------
        elif data_type == 'tensiometer':
            if isinstance(filename, str):
                df = in_meas.read_tensiometers(filename)
            else:
                df = filename
                filename = None

            units = '$kPa$'
            obs_cov_type = None 

        # ERT type read
        # ---------------------------------------------------------------------
        elif data_type == 'ERT':

            obs_cov_type = 'reciprocal_err'
            if 'obs_cov_type' in kwargs:
                obs_cov_type = kwargs['obs_cov_type']
           

            data_format = 'resipy'
            if 'data_format' in kwargs:
                data_format = kwargs['data_format']
                dict_obs_2add.update(
                     data_format = data_format
                     )

            elecs = []
            if 'elecs' in kwargs:
                elecs = kwargs['elecs']
            
            dict_obs_2add.update(
                     elecs = elecs
                     )


            df = in_meas.read_ERT(filename)
            units = '$\Omega$'


        # no file specified (raise error)
        # ---------------------------------------------------------------------
        else:
            print('no file specified')

            
        dict_obs_2add.update(filename = filename,
                            data_type =  data_type,
                            units= units,  # units
                            data=  df,
                            data_err=  data_err,
                            mesh_nodes =  mesh_nodes,
                            assimilation_times=  tA
                            )
        
        
        # add optionnal data metadata to the dictionnary
        # ---------------------------------------------------------------------
        if 'meta' in kwargs:
            meta = kwargs['meta']
            
            dict_obs_2add.update(meta)
            print(dict_obs_2add)
            
        
            
        # data covariance and perturbation (if required)
        # ---------------------------------------------------------------------
        
        data_cov, dataPert = self.prepare_observations(dict_obs_2add,
                                                      perturbate = False,
                                                      obs_cov_type = obs_cov_type
                                                      )
        dict_obs_2add.update(
                             data_cov = data_cov,
                             dataPert = dataPert
                             )
                        
        
        dict_obs_2add = OrderedDict(dict_obs_2add)
        
        # check if assimilation already existing andincrement sensor name if so
        # ---------------------------------------------------------------------
        sensor_name =  data_type
        
        if tA in self.dict_obs.keys():
            print('already existing assimilation time')
            for k in self.dict_obs[tA]:
                k
            if data_type in k:
                match = re.match(r"([a-z]+)([0-9]+)", k, re.I)
                if match:
                    items = match.groups()
                    it = int(items[1]) +1
                else:
                    it = 1
                sensor_name = sensor_name + str(it)
                self.dict_obs[tA][sensor_name] = {} 
                dict_obs_2add.update(sensor_name = sensor_name)
        else:
            self.dict_obs[tA] = {} 
            self.dict_obs[tA][sensor_name] = {} 
            dict_obs_2add.update(sensor_name = sensor_name)


        self.dict_obs = self._add_to_obs_dict(tA,sensor_name,dict_obs_2add)
                
        return self.dict_obs

    def _add_to_obs_dict(self, tA, sensor_name, dict_obs_2add):
        '''
        Merge a new observation into an existing dict

        Parameters
        ----------
        dict_obs_2add : dict
            dict of the new observation.

        Returns
        -------
        dict
            updated dict with all observations.

        '''

        for key in dict_obs_2add.keys():
            self.dict_obs[tA][sensor_name][key] = dict_obs_2add[key]
      
        # self.dict_obs = self.dict_obs | dict_obs_2add
        # self.dict_obs.update(dict_obs_2add)

        return self.dict_obs


    def prepare_observations(self, 
                             dict_obs = [], 
                             perturbate = False,
                             obs_cov_type = None,
                             Bishop=False, **kwargs):
        '''
        THIS SHOULD BE MOVED TO DA CLASS


        prepare observations before DA


        1. Measurement perturbation (if required)
        2. Compute matrice covariance 
        3. normalise data (if necessary)

        Parameters
        ----------
        dict_obs : dict
            dict containing measure data + metadata.

        Returns
        -------
        TYPE
            data_cov, DataPert.

        '''


        # compute individual observation data covariance 
        # ----------------------------------------------------
        data_cov = []
        if obs_cov_type is not None:
            data_cov = self.init_obs_cov_matrice(dict_obs, obs_cov_type)
        
        
        # compute individual observation perturbation covariance 
        # ----------------------------------------------------
            
        data_pert = []
        if perturbate == True:
            data_pert = self.perturbate_obs(dict_obs)
            

        # stacked data covariance (multiple observations at the same time)
        # ----------------------------------------------------
        if 'stacked_data_cov' in kwargs:
            self.stacked_data_cov = kwargs['stacked_data_cov']
        else:

            self.stacked_data_cov.append(data_cov)
    





        # (Bishop et al., 2001) --> not require the perturbation of observations
        # taking advantage of the high time resolution of the collected data.
        # ---------------------------------------------------------------------
        # The EnKF algorithm implemented here is actually an en-
        # semble transform Kalman filter (Bishop et al., 2001) that
        # does not require the perturbation of observations. On the
        # other hand, the measurement error covariance matrix, R,
        # must be assumed to be known a priori.


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



        return data_cov, data_pert


    def perturbate_obs(self, dict_obs):
        '''
        Not yet implemented
        '''

            # if (self.data_cov.ndim == 0):
            #     DataPerturbation=np.sqrt(self.data_cov)*rn.randn(1, self.EnSize)
            # elif (self.data_cov.ndim == 2):
            #     # Compute SVD of Data Covariance to generate noise
            #     U, s, V=np.linalg.svd(self.data_cov, full_matrices=False)
            #     DataPerturbation=np.dot(np.dot(U, np.diag(np.sqrt(s))),
            #                               rn.randn(self.Data.shape[1], self.EnSize))
            # else:
            #     print('Data Covariance should be matrix or scalar', '\n')

            # DataArray=np.tile(
            #     self.Data[i, :], (self.EnSize, 1)).transpose() + DataPerturbation


        pass


    def init_obs_cov_matrice(self,dict_obs, obs_cov_type):
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

        if obs_cov_type == 'reciprocal_err':
            # print(dict_obs)
            data_cov = np.diag(1/abs(dict_obs['data']['recipError'].to_numpy()))
        elif obs_cov_type == 'data_err':
            
            try:
                err_arr = np.tile(dict_obs['data_err'],len(dict_obs['data']))
                data_cov = np.diag(err_arr)
            except:
                err_arr = np.tile(dict_obs['data_err'],len(dict_obs['data']['a']))

            
            


        return data_cov 





    # ------------------------------------------------------------------------
    # %% Mapping of properties to zones/mesh
    # ------------------------------------------------------------------------
    def map_prop2mesh(self,dict_props):
        '''
        Add a map of a given physical property to the mesh

        Parameters
        ----------
        dict_props : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        
        if hasattr(self,'mesh_pv_attributes') == False:
            self.create_mesh_vtk()
        
        nnodes = int(self.mesh_pv_attributes.points.size/3)
        for dp in dict_props.keys():
            print(dp)
            if type(dict_props[dp]) == float: # homogeneous properties
                self.update_mesh_vtk(prop=dp, prop_value=np.ones(nnodes)*dict_props[dp])
            elif type(dict_props[dp]) == list: # homogeneous properties
                if len(dict_props[dp]) == 1:
                    self.update_mesh_vtk(prop=dp, prop_value=np.ones(nnodes)*dict_props[dp])


        
        pass

    def map_prop2zone(self,dict_props, prop):
               
        if hasattr(self,'zone_xyz')==False:
            warnings.warn('no known existing zones.')
            pass
        
        else:
            prop_zones = np.zeros(np.shape(self.zone_xyz))
            for z in range(len(np.unique(self.zone_xyz))):
                prop_zones[self.zone_xyz==z+1]=dict_props[prop][z]
            
            return prop_zones
    
    # ------------------------------------------------------------------------
    #%% Plot call
    # ------------------------------------------------------------------------
    def plot(self, prop, **kwargs):
        '''
        Call and parse to cathy.plotter from the main CATHY class

        Parameters
        ----------
        prop : str
            property to plot.
        **kwargs : kwargs

        Returns
        -------
        None.

        '''

        raise NotImplementedError('plotter from main class not yet implemented')
        # prop_zones = self.map2zone(prop)
        
        plt, ax = plt_CT.dem_plot_2d_top(SPP_zones,
                                         label=show)
            
        if prop == 'ZROOT':
            df=out_CT.dem_plot_2d_top(path)
            
            
    #%% Read outputs/inputs 
    # ------------------------------------------------------------------------
    def read_outputs(self, filename, **kwargs):
        '''
        Read CATHY format output file


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


    def read_inputs(self, filename):
        '''
        Read CATHY format input file


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

    # -------------------------------------------------------------------#
    # %% utils
    # -------------------------------------------------------------------#

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
        

        if np.array(node_coords).ndim <= 1:
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

            min(self.grid3d['mesh3d_nodes'][:, 2])
            min(self.grid3d['mesh3d_nodes'][:, 1])
            closest_idx.append(np.argmin(d))
            closest.append(self.grid3d['nodes_idxyz'][closest_idx[i], 1:])

            threshold=1e-1
            if d[np.argmin(d)] > threshold:
                self.console.print('no node close to the required points')
                print(d[np.argmin(d)])
                print(self.grid3d['nodes_idxyz'][closest_idx[i], 1:])
                print(nc)

        return closest_idx, closest


    def rich_display(self, title="Star Wars Movies", **kwargs):
        """
        Describe the variable state and fate during the simulation with a rich table

        Returns
        -------
        None.

        """
        self.console.print(eval('self.' + str(title)))






