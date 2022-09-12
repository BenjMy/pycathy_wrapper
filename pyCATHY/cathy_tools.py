"""Main class controlling the wrapper
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
from git import Repo
import pyvista as pv 
import pandas as pd
import time
import rich.console
from rich.progress import track
from rich import print

from pyCATHY.plotters import cathy_plots as plt_CT
# from pyCATHY.DA.cathy_DA import DA
from pyCATHY.DA import cathy_DA

from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import sensors_measures as in_meas

from pyCATHY.ERT import petro_Archie as Archie
from pyCATHY import cathy_utils as utils_CT


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
    Run multiple exe files in //
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
        
            
    
# -----------------------------------------------------------------------------

class CATHY():  # IS IT GOOD PRACTICE TO PASS DA CLASS HERE ? I think we sould better pass main CATHY into children classes
    """ 

    Main CATHY object

    When instantiated it creates the tree project directories with 'prj_name' as root folder.
    The src files are fetched from the online repository if not existing (note that it is possible to call a specific version).

    Parameters
    ----------
    None

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
        # (self.dict_obs is assimiliatiom times size)
        self.dict_obs = OrderedDict() # dictionnary containing all the observation data
        # (self.stacked_data_cov is data_size * assimiliatiom times size)
        self.stacked_data_cov = [] # merged data covariance matrice of all the observation data and all times 
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
                p = os.system(bashCommand)
                # run it twice (to avoid the first error)
                p = os.system(bashCommand)
                self.console.print(":cooking: [b]gfortran compilation[/b]")

            except:
                print("bash cmd not recognized")

            os.chdir(self.workdir)
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
        if verbose:
            print(p.stdout)
            print(p.stderr)
        if "DEM resolution is not a multiple of the rivulet spacing!" in p.stdout:
            print("DEM resolution is not a multiple of the rivulet spacing!")
            sys.exit()
        if "catchment with more than one outlet cell!" in p.stdout:
            print("catchment with more than one outlet cell!")
            sys.exit()
        if KeepOutlet == False:
            print("remove outlet")
            idoutlet=np.where(self.DEM==min(np.unique(self.DEM)))
            self.DEM[idoutlet[0],idoutlet[1]] = max(np.unique(self.DEM))

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
        try:
            shutil.move(
                os.path.join(self.workdir, self.project_name, "src", self.processor_name
                ),
                os.path.join(self.workdir, self.project_name,
                             self.processor_name)
            )
        except:
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
        recompile : Bool, optional
            recomplile CATHY src files. The default is True.
        runProcess : Bool, optional
            exectute compile CATHY exe. The default is True.
        verbose : Bool, optional
            Output CATHY log. The default is False.
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
        # VERY VERY IMPORTANT NEVER COMMENT !
        self.update_parm(**kwargs, verbose=verbose)
        self.update_cathyH(**kwargs)

        if recompile == True:
            # recompile
            # --------------------------------------------------------------------
            os.chdir(os.path.join(self.workdir, self.project_name, "src"))
            self.recompileSrc()

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
                
                verbose = False
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
                    
                    
                threshold_rejected = 10
                if 'threshold_rejected' in kwargs:
                    threshold_rejected = kwargs['threshold_rejected']
                    
                
                self.damping = 1                  
                if 'damping' in kwargs:
                    self.damping = kwargs['damping']
                    
                    
                self._run_DA(callexe, parallel, 
                             DA_type,
                             dict_obs,
                             list_update_parm,
                             dict_parm_pert,
                             list_assimilated_obs,
                             open_loop_run,
                             threshold_rejected,
                             verbose
                             )

            # case of simple simulation
            # ----------------------------------------------------------------
            else:
                p = subprocess.run([callexe], text=True, capture_output=True)

                if verbose == True:
                    print(p.stdout)
                    print(p.stderr)

                os.chdir(os.path.join(self.workdir))
                
                self.grid3d = in_CT.read_grid3d(self.project_name)


            # computation time
            # ----------------------------------------------------------------
            t1 = time.time()
            self.total_computation_time = t1 - t0
            

        return


#%%  Observation Processing

    def _parse_ERT_metadata(self,key_value):
        '''
        Extract ERT metadata information form obs dict
        '''
        ERT_meta_dict = {}
        ERT_meta_dict['forward_mesh_vtk_file'] = key_value[1]['ERT']['forward_mesh_vtk_file'] 
        ERT_meta_dict['pathERT'] = os.path.split(key_value[1]['ERT']['filename'])[0]
        ERT_meta_dict['seq'] = key_value[1]['ERT']['sequenceERT'] 
        ERT_meta_dict['electrodes'] = key_value[1]['ERT']['elecs'] 
        ERT_meta_dict['noise_level'] = key_value[1]['ERT']['data_err']        
        ERT_meta_dict['porosity'] = self.Archie_parms['porosity']
        ERT_meta_dict['data_format'] = key_value[1]['ERT']['data_format'] 
        return ERT_meta_dict
        

    def _map_ERT(self,state,
                path_fwd_CATHY,
                ens_nb,
                **kwargs):
        ''' 
        Mapping of state variable to observation (predicted) 
        ERT using pedophysical transformation H
        '''
    
        savefig = True
        if 'savefig' in kwargs:
           savefig = kwargs['savefig'] 
        
        # search key value to identify time and method
        # --------------------------------------------
        tuple_list_obs = list(self.dict_obs.items())
        key_value = tuple_list_obs[self.count_DA_cycle]  
        
        # Load ERT metadata information form obs dict
        # -------------------------------------------
        ERT_meta_dict = self._parse_ERT_metadata(key_value)
        
        Hx_ERT, df_Archie = Archie.SW_2_ERa_DA(self.project_name,
                                            self.Archie_parms, 
                                            ERT_meta_dict['porosity'], 
                                            ERT_meta_dict['pathERT'],
                                            ERT_meta_dict['forward_mesh_vtk_file'],
                                            ERT_meta_dict['electrodes'],
                                            ERT_meta_dict['seq'],
                                            path_fwd_CATHY,
                                            df_sw = state[1], # kwargs
                                            data_format= ERT_meta_dict['data_format'] , # kwargs
                                            DA_cnb = self.count_DA_cycle, # kwargs
                                            Ens_nb=ens_nb, # kwargs
                                            savefig=savefig, # kwargs
                                            noise_level = ERT_meta_dict['noise_level'],# kwargs
                                            dict_ERT = key_value[1]['ERT']#  kwargs
                                            )

        df_Archie['OL'] = np.ones(len(df_Archie['time']))*False
        self._add_2_ensemble_Archie(df_Archie)
        
        
        return Hx_ERT
                        
    def _map_ERT_parallel_DA(self,
                             ENS_times,
                             ERT_meta_dict,
                             key_time,
                             path_fwd_CATHY_list,
                             DA_cnb,
                             ):
         ''' 
         Parallel mapping of ERT data using pedophysical transformation H
         '''                         
         Hx_ERT_ens = []      
                        
         # freeze fixed arguments of Archie.SW_2_ERa_DA       
         # -----------------------------------------------------------------
         ERTmapping_args = partial(Archie.SW_2_ERa_DA, 
                                   self.project_name,
                                   self.Archie_parms,
                                   ERT_meta_dict['porosity'],
                                   ERT_meta_dict['pathERT'],
                                   ERT_meta_dict['forward_mesh_vtk_file'],
                                   ERT_meta_dict['electrodes'],
                                   ERT_meta_dict['seq'],
                                   data_format= ERT_meta_dict['data_format'],
                                   DA_cnb = DA_cnb,
                                   savefig=True,
                                   noise_level = ERT_meta_dict['noise_level'],# kwargs
                                   dict_ERT = key_time[1]['ERT']#  kwargs
                                   ) 
        
         # // run using ensemble subfolders path as a list    
         # -----------------------------------------------------------------
         with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            results_mapping = pool.map(ERTmapping_args, 
                                         path_fwd_CATHY_list
                                           )
            print(f"x= {path_fwd_CATHY_list}, PID = {os.getpid()}")
         
         for ens_i in range(len(path_fwd_CATHY_list)):
             df_Archie = results_mapping[ens_i][1]       
             df_Archie['OL'] = np.zeros(len(df_Archie))
             self._add_2_ensemble_Archie(df_Archie)
             Hx_ERT_ens_i =  results_mapping[ens_i][0]

             if 'pygimli' in self.dict_obs[key_time[0]]['ERT']['data_format']:
                Hx_ERT_ens = self._add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_ens_i['rhoa'])
             else:
                Hx_ERT_ens = self._add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_ens_i['resist'])
                
         return Hx_ERT_ens
                
                    
    def _map_ERT_parallel_OL(self,
                             ENS_times,
                             ERT_meta_dict,
                             key_time,
                             path_fwd_CATHY_list,
                             ):
        ''' 
        Parallel mapping of ERT data using pedophysical transformation H
        case of the open Loop = nested loop with ensemble time
        ''' 

        Hx_ERT_ens = []
        for t in range(len(ENS_times)):
            print('t_openLoop mapping:' + str(t))

            # freeze fixed arguments of Archie.SW_2_ERa       
            # -----------------------------------------------------------------
            ERTmapping_args = partial(Archie.SW_2_ERa_DA, 
                                      self.project_name,
                                      self.Archie_parms, 
                                      ERT_meta_dict['porosity'], 
                                      ERT_meta_dict['pathERT'],
                                      ERT_meta_dict['forward_mesh_vtk_file'],
                                      ERT_meta_dict['electrodes'],
                                      ERT_meta_dict['seq'],
                                      data_format= ERT_meta_dict['data_format'],
                                      time_ass = t,
                                      savefig=True,
                                      noise_level = ERT_meta_dict['noise_level'] ,
                                      dict_ERT = key_time[1]['ERT']
                                      ) 
            # 
            # -----------------------------------------------------------------
            with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
                results_mapping_time_i = pool.map(ERTmapping_args, 
                                                  path_fwd_CATHY_list
                                                  )
                print(f"x= {path_fwd_CATHY_list}, PID = {os.getpid()}")

            for ens_i in range(len(path_fwd_CATHY_list)):
                 Hx_ERT_time_i =  results_mapping_time_i[ens_i][0]     
                 
                 # print(np.shape(results_mapping_time_i))
                 df_Archie = results_mapping_time_i[ens_i][1]     
                 
                 # print(df_Archie)
                 df_Archie['OL'] = np.ones(len(df_Archie))
                 self._add_2_ensemble_Archie(df_Archie)

    
                 if 'pygimli' in ERT_meta_dict['data_format']:
                    Hx_ERT_ens = self._add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_time_i['rhoa'])
                 else:
                    Hx_ERT_ens = self._add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_time_i['resist'])
                     
        # prediction_ERT = np.reshape(Hx_ERT_ens,[self.NENS,
        #                                         len(Hx_ERT_ens[0]),
        #                                         len(ENS_times)])  # (EnSize * data size * times)
        return Hx_ERT_ens
    
    
    def _map_ERT_parallel(self, 
                            path_fwd_CATHY_list,
                            list_assimilated_obs='all', 
                            default_state = 'psi',
                            verbose = False,
                            **kwargs):
        ''' Mapping of state variable to observation (predicted) ERT using pedophysical transformation H,
        // run using ensemble subfolders path as a list
        '''

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
        key_time = tuple_list_obs[self.count_DA_cycle]
        # Load ERT metadata information form obs dict
        # -------------------------------------------       
        ERT_meta_dict = self._parse_ERT_metadata(key_time)
        
        if len(ENS_times)>0: # case of the open Loop = nested loop with ensemble time
            Hx_ERT_ens = self._map_ERT_parallel_OL(
                                                    ENS_times,
                                                    ERT_meta_dict,
                                                    key_time,
                                                    path_fwd_CATHY_list,
                                                )                                      
        else:
            Hx_ERT_ens = self._map_ERT_parallel_DA(
                                                    ENS_times,
                                                    ERT_meta_dict,
                                                    key_time,
                                                    path_fwd_CATHY_list,
                                                    DA_cnb,
                                                )                          
                                                
        prediction_ERT = np.vstack(Hx_ERT_ens).T  # (EnSize)
        
        # self.dict_obs
            
        return prediction_ERT

    def _evaluate_perf_OL(self,
                         parallel,
                         list_assimilated_obs,
                         path_fwd_CATHY_list,
                         ENS_times
                         ):
        
        prediction_OL = []
        if parallel:
            if 'ERT' in list_assimilated_obs:
                prediction_OL_ERT = self._map_ERT_parallel(path_fwd_CATHY_list,
                                                    savefig = True,
                                                    DA_cnb = self.count_DA_cycle,
                                                    ENS_times=ENS_times,
                                                    )   
                # prediction_OL_ERT is meas_size * ens_size * ENS_times size
                np.shape(prediction_OL_ERT)
            else:
                for t in range(len(ENS_times)):
                # for t in track(range(len(ENS_times)), description="OL Mapping observations to predicted obs..."):
                    prediction_stacked = self.map_states2Observations(
                                                        list_assimilated_obs,
                                                        ENS_times=ENS_times,
                                                        savefig=False,
                                                        parallel=parallel,
                                                        ) 
                    prediction_OL.append(prediction_stacked)
                    
        else:    
            for t in range(len(ENS_times)):
            # for t in track(range(len(ENS_times)), description="OL Mapping observations to predicted obs..."):
                prediction_stacked = self.map_states2Observations(
                                                    list_assimilated_obs,
                                                    ENS_times=ENS_times,
                                                    savefig=False,
                                                    parallel=parallel,
                                                    ) 
                # prediction_stacked is meas_size * ens_size
                prediction_OL.append(prediction_stacked)
                # prediction_OL is meas_size * ens_size * ENS_times size


        if parallel:
            if 'ERT' in list_assimilated_obs:
                prediction_OL.append(prediction_OL_ERT)
                
                if np.shape(prediction_OL)[0]>1:
                    prediction_OL =  np.hstack(prediction_OL)
            
        prediction_OL = np.reshape(prediction_OL,[
                                                  np.shape(prediction_OL)[1],
                                                  self.NENS,
                                                  len(ENS_times),
                                                  ]
                                  )
        for t in range(len(ENS_times)):
            # print(str(t) + 't perf OL')
            data_t, _  = self._get_data2assimilate(
                                                    list_assimilated_obs,
                                                    time_ass=t,
                                                )                    
            self._performance_assessement(
                                            list_assimilated_obs, 
                                            data_t, 
                                            prediction_OL[:,:,t],
                                            t_obs=t,
                                            openLoop=True,
                                        )
        return prediction_OL
        
                
    def _DA_openLoop(self,
                     ENS_times,
                     list_assimilated_obs,
                     parallel,
                     verbose=False):
        '''
        Run open Loop (no update/assimilation) hydro simulation for an ensemble of realisation
        Evaluate the performance by comparison with measured data after mapping

        Parameters
        ----------
        ENS_times : list
            List of the time steps (in sec) where observations are assimilated.
        list_assimilated_obs : list
            List of assimilation observations.
        parallel : Bool, optional
            Parallelize. The default is True.
        verbose : Bool, optional
            Display prints. The default is False.

        Returns
        -------
        Prediction_OL.
        Performance dataframe.

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
        
        # Loop over ensemble realisations one by one
        # ------------------------------------------------------------
        else:
            for ens_i in range(self.NENS):
                os.chdir(os.path.join(self.workdir,
                                      self.project_name, 
                                      'DA_Ensemble/cathy_' + str(ens_i+1)))
                len(list(self.atmbc['time']))
                self._run_hydro_DA_openLoop(time_of_interest=list(self.atmbc['time']),
                                              nodes_of_interest=[],
                                              simu_time_max=max(list(self.atmbc['time'])),
                                              ens_nb=ens_i+1                                 
                                              )
                
        # save into the DA_df dataframe
        # -----------------------------------------------------------------
        path_fwd_CATHY_list = []
        for ens_i in range(self.NENS):
            path_fwd_CATHY_list.append(os.path.join(self.workdir,
                                    self.project_name, 
                                    'DA_Ensemble/cathy_' + str(ens_i+1)))

            os.chdir(os.path.join(self.workdir,
                                  self.project_name, 
                                  'DA_Ensemble/cathy_' + str(ens_i+1)))
            df_psi = self.read_outputs(filename='psi',
                                        path=os.path.join(os.getcwd(),
                                                          'output'))
            df_sw = self.read_outputs(filename='sw',
                                        path=os.path.join(os.getcwd(),
                                                          'output'))
            
            
            shift = len(df_psi)-self.parm["NPRT"]
            if shift<0 | shift>2:
                print('Error for the ensemble nb:' + str(ens_i))
                raise ValueError('Error on the simulation:' 
                                  'nb of times contained in the outputs files is too small;' 
                                  'Check ')                
            for t in range(np.shape(df_psi)[0]-2):
                self._DA_df(state=[df_psi[t+shift,:], df_sw[t+shift,:]],
                            t_ass=t, 
                            openLoop=True,
                            ens_nb=ens_i+1)
             
        
        prediction_OL = self._evaluate_perf_OL(                    
                                                parallel,
                                                list_assimilated_obs,
                                                path_fwd_CATHY_list,
                                                ENS_times
                                                )
                   
        # ------------------------------------------------------
        # END of Open Loop simulation and performance evaluation
        # ------------------------------------------------------     
        pass
    

    def _mark_invalid_ensemble(self,ens_valid,
                                  prediction,
                                  ensemble_psi,
                                  ensemble_sw,
                                  analysis,
                                  analysis_param):
        ''' mark invalid ensemble - invalid ensemble are filled with NaN values '''          

        ensemble_psi_valid = np.empty(ensemble_psi.shape)
        ensemble_psi_valid[:] = np.NaN
        ensemble_psi_valid[:,ens_valid]  = ensemble_psi[:,ens_valid]            

        ensemble_sw_valid = np.empty(ensemble_psi.shape)
        ensemble_sw_valid[:] = np.NaN
        ensemble_sw_valid[:,ens_valid]  = ensemble_sw[:,ens_valid]        
        
        analysis_valid = np.empty(ensemble_psi.shape)
        analysis_valid[:] = np.NaN
        analysis_valid[:,ens_valid]  = analysis         

        prediction_valid = np.empty([prediction.shape[0],
                                     ensemble_psi.shape[1]])
        prediction_valid[:] = np.NaN
        prediction_valid[:,ens_valid]  = prediction
        
        analysis_param_valid = []
        if len(analysis_param[0])>0:
            analysis_param_valid = np.empty([analysis_param.shape[1],
                                            ensemble_psi.shape[1]])
            analysis_param_valid[:] = np.NaN
            analysis_param_valid[:,ens_valid]  = analysis_param.T     
        
        return (prediction_valid, 
                ensemble_psi_valid, 
                ensemble_sw_valid, 
                analysis_valid, 
                analysis_param_valid
                )
    
    
    
    def set_Archie_parm(self,
                        porosity=[],
                        rFluid_Archie=[1.0],
                        a_Archie=[1.0],
                        m_Archie=[2.0],
                        n_Archie=[2.0], 
                        pert_sigma_Archie=[0]
                        ):
        '''
        Dict of Archie parameters. Each type of soil is describe within a list
        Note that if pert_sigma is not None a normal noise is 
        added during the translation of Saturation Water to ER

        Parameters
        ----------
        rFluid : TYPE, optional
            Resistivity of the pore fluid. The default is [1.0].
        a : TYPE, optional
            Tortuosity factor. The default is [1.0].
        m : TYPE, optional
            Cementation exponent. The default is [2.0]. 
            (usually in the range 1.3 -- 2.5 for sandstones)
        n : TYPE, optional
            Saturation exponent. The default is [2.0].
        pert_sigma_Archie : TYPE, optional
            Gaussian noise to add. The default is None.
        
        ..note:
            Field procedure to obtain tce he covariance structure of the model 
            estimates is described in Tso et al () - 10.1029/2019WR024964
            "Fit a straight line for log 10 (S) and log 10 (ρ S ) using the least-squares criterion. 
            The fitting routine returns the covariance structure of the model estimates, which can be used to de-
            termine the 68% confidence interval (1 standard deviation) of the model estimates.""

        '''
        if len(porosity)==0:
            porosity = self.soil_SPP['SPP'][:, 4][0]   
        if not isinstance(porosity, list):
            porosity = [porosity]
        self.Archie_parms = {
                             'porosity':porosity, 
                             'rFluid_Archie':rFluid_Archie, 
                             'a_Archie':a_Archie, 
                             'm_Archie':m_Archie, 
                             'n_Archie':n_Archie, 
                             'pert_sigma_Archie':pert_sigma_Archie
                             }
        
        pass


    def _mapping_petro_init(self):
        ''' Initiate Archie and VGP'''
        
        print('Initiate Archie and VGP')
        porosity = self.soil_SPP['SPP'][:, 4][0]
        if not isinstance(porosity, list):
            porosity = [porosity]
        if hasattr(self,'Archie_parms') == False:
            print('Archie parameters not defined set defaults')
            self.Archie_parms = {'porosity': porosity,
                                 'rFluid_Archie':[1.0],
                                 'a_Archie':[1.0],
                                 'm_Archie':[2.0],
                                 'n_Archie':[2.0],
                                 'pert_sigma_Archie':[0]
                                 }
        if hasattr(self,'VGN_parms') == False:
            print('VGN parameters not defined set defaults: empty dict {} - not implemented yet')
            self.VGN_parms = {}
            
        pass
         
    

    def run_ensemble_hydrological_model(self,parallel,verbose,callexe):
        ''' multi run CATHY hydrological model from the independant folders composing the ensemble '''
        
        # ----------------------------------------------------------------
        if parallel == True:
            pathexe_list = []
            # for ens_i in range(self.NENS):
                
            for ens_i in self.ens_valid:
                path_exe = os.path.join(self.workdir,
                                        self.project_name,
                                        'DA_Ensemble/cathy_' + str(ens_i+1))
                pathexe_list.append(path_exe)       
            with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
                result = pool.map(subprocess_run_multi, pathexe_list)

                                        
                if verbose == True:
                    self.console.print(result)
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
                    
        pass
                    
    
    def update_pert_parm_dict(self,update_key,list_update_parm,analysis_param_valid):
        ''' update dict of perturbated parameters i.e. add a collumn with new params'''
        
        index_update = 0
        if self.count_DA_cycle > 0:
           for pp in list_update_parm:
               if 'St. var.' in pp:
                   index_update = index_update + 1
                   pass
               else:
                   update_key = 'update_nb' + str(self.count_DA_cycle)
                   self.dict_parm_pert[pp][update_key] = analysis_param_valid[index_update-1]
                   index_update = index_update + 1
   
           # check after analysis
           # ----------------------------------------------------------------
           id_valid_after = self._check_after_analysis(
                                                       update_key,
                                                       list_update_parm
                                                       )
           intersection_set = set.intersection(
                                               set(id_valid_after), 
                                               set(self.ens_valid)
                                               )
           self.ens_valid = list(intersection_set)
           
        pass
                            
                            
    def _run_DA(self, callexe, parallel,
                        DA_type,
                        dict_obs,
                        list_update_parm,
                        dict_parm_pert,
                        list_assimilated_obs,
                        open_loop_run,
                        threshold_rejected,
                        verbose,
                        **kwargs
                        ):
                                    
                             
        '''
        
        Run Data Assimilation 
        
        .. note::

            Steps are:
            1. DA init (create subfolders)
            2a. run CATHY hydrological model (open loop)
            2b. run CATHY hydrological model recursively using DA times
                <- Loop-->
                    3. check before analysis
                    4. analysis
                    5. check after analysis
                    6. update input files
                <- Loop-->
 
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
        
        
        # check if dict_obs is ordered in time
        # ------------------------------------       
        if (all(i < j for i, j in zip(list(dict_obs.keys()), list(dict_obs.keys())[1:]))) is False:
            raise ValueError('Observation List is not sorted.')
       
        # dict_obs.keys()
        # self.dict_obs.keys()
        if hasattr(self,'dict_obs') is False:
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
            
        # data_measure_df = self.dictObs_2pd()     


        # start DA cycle counter
        # -------------------------------------------------------------------
        self.count_DA_cycle = 0
        self.count_atmbc_cycle = 0
        # (the counter is incremented during the update analysis)
        
        # initiate DA
        # -------------------------------------------------------------------
        list_update_parm = self._DA_init(
                                          NENS=NENS, # ensemble size
                                          ENS_times=ENS_times, # assimilation times
                                          parm_pert=dict_parm_pert,
                                          update_parm_list = list_update_parm
                                          )
                            
        # initiate mapping petro
        # -------------------------------------------------------------------
        self._mapping_petro_init()          
        
        # update the perturbated parameters 
        # --------------------------------------------------------------------
        self.update_ENS_files(dict_parm_pert, 
                              update_parm_list='all', #list_update_parm
                              cycle_nb=self.count_DA_cycle)
        
        all_atmbc_times = self.atmbc['time']        
        # -------------------------------------------------------------------
        if open_loop_run:
            self._DA_openLoop(
                              ENS_times,
                              list_assimilated_obs,
                              parallel)     
        # end of Open loop - start DA
        
        # -------------------------------------------------------------------
        # update input files ensemble again (time-windowed)
        # ---------------------------------------------------------------------
        self._update_input_ensemble(NENS, 
                                    list(self.atmbc['time']), 
                                    dict_parm_pert, 
                                    update_parm_list='all')  #list_update_parm

        # -----------------------------------
        # Run hydrological model sequentially = Loop over atmbc times (including assimilation observation times)
        # -----------------------------------
        # self.sequential_DA()
        
        # TO CHANGE HERE
        for t_atmbc in all_atmbc_times: #atmbc times include assimilation observation times
            print(t_atmbc)
            # t_atmbc = self.atmbc['time'][-2]
            
            self.run_ensemble_hydrological_model(parallel,verbose,callexe)  
            os.chdir(os.path.join(self.workdir))
            
            # self.plot_ini_state_cov()
            # def plot_ini_state_cov():
            #     ensemble_psi, ensemble_sw, ens_size, sim_size = self._read_state_ensemble()
                
            # check scenario (result of the hydro simulation)
            # ----------------------------------------------------------------       
            rejected_ens = self._check_before_analysis(
                                                       self.ens_valid,
                                                       threshold_rejected,
                                                       )
            id_valid= ~np.array(rejected_ens)          
            self.ens_valid = list(np.arange(0,self.NENS)[id_valid])
            
            
            # define subloop here 
            # if 'optimize_inflation' in DA_type:
            # threshold_convergence = 10
            # while self.df_performance < threshold_convergence:
                
            # self.atmbc['time']
        
            if t_atmbc in ENS_times:
                
                print('t=' + str(t_atmbc))

                # map states to observation = apply H operator to state variable
                # ----------------------------------------------------------------                               
                prediction = self.map_states2Observations( 
                                                            list_assimilated_obs,
                                                            default_state = 'psi',
                                                            parallel=parallel,
                                                           )
                # print(len(prediction))
                # print(np.shape(prediction))
                # ValueError: operands could not be broadcast together with shapes (612,) (1836,) 
                
                # data2test , obs2evaltest = self._get_data2assimilate(list_assimilated_obs)
                # # len(data2test)
                
                # x = np.linspace(prediction.min(),
                #                 prediction.max(),
                #                 100)
                # y = x
                # fig, ax = plt.subplots(2,1)
                # for ii in range(np.shape(prediction)[1]):
                #     ax[0].scatter(prediction[:,ii],data2test)
                #     ax[0].plot(x, y, '-r', label='y=x')
                # ax[0].set_title('nb of ens' + str(np.shape(prediction)[1]))
    
                # ax[1].plot(data2test)
                # plt.savefig('test' + str(DA_cnb))
    
                # plt.plot(data2test)
                # plt.plot(prediction[:,ii])
                # ax.scatter(prediction[:,ii],data2test)
                # ax.set_aspect('equal')
    
                # fig.tight_layout()
    
                # DA analysis
                # ----------------------------------------------------------------       
                (ensemble_psi,
                 ensemble_sw,
                 data, 
                 analysis, 
                 analysis_param) = self._DA_analysis(prediction,
                                                     DA_type,
                                                     list_update_parm,
                                                     list_assimilated_obs,
                                                     ens_valid=self.ens_valid)
                                                     
                # DA mark_invalid_ensemble
                # ----------------------------------------------------------------   
                (prediction_valid, 
                 ensemble_psi_valid, 
                 ensemble_sw_valid,
                 analysis_valid, 
                 analysis_param_valid) = self._mark_invalid_ensemble(self.ens_valid,
                                                                    prediction,
                                                                    ensemble_psi,
                                                                    ensemble_sw,
                                                                    analysis,
                                                                    analysis_param)
                                                                    
                                                                    
                                                                    
                # check analysis quality
                # ----------------------------------------------------------------
                
                self.console.print(':face_with_monocle: [b]check analysis performance[/b]')
                self._performance_assessement(list_assimilated_obs, 
                                              data, 
                                              prediction_valid,
                                              t_obs=self.count_DA_cycle)
                
                
                # the counter is incremented here                           
                # ----------------------------------------------------------------
                self.count_DA_cycle = self.count_DA_cycle + 1 
                
                
                # # update parameter dictionnary
                # ----------------------------------------------------------------
                def check_ensemble(ensemble_psi_valid,ensemble_sw_valid):
                    if np.any(ensemble_psi_valid>0):
                        # raise ValueError('positive pressure heads observed')
                        print('!!!!!!positive pressure heads observed!!!!!')
                        psi_2replace = np.where(ensemble_psi_valid>=0)
                        ensemble_psi_valid_new = ensemble_psi_valid
                        ensemble_psi_valid_new[psi_2replace]=-1e-3

                check_ensemble(ensemble_psi_valid,ensemble_sw_valid)
    
                self.update_pert_parm_dict(update_key,list_update_parm,analysis_param_valid)
                
                             
            else:
                self.console.print(':confused: No observation for this time - run hydrological model only')
                print('!!!!!!!!! shoetcutttt here ensemble are anot validated!!!!!!!!!! S')
                ensemble_psi_valid, ensemble_sw_valid, ens_size, sim_size = self._read_state_ensemble()
                # analysis_valid = np.empty(ensemble_psi_valid.shape)
                # analysis_valid[:] = np.NaN                
                analysis_valid = ensemble_psi_valid
                
            self.count_atmbc_cycle = self.count_atmbc_cycle + 1 

            print('------ end of time step (s) -------' + str(int(t_atmbc)) + '/' + str(int(all_atmbc_times[-1])) + '------')
            print('------ end of atmbc update --------' + str(self.count_atmbc_cycle) + '/' + str(len(all_atmbc_times)-1) + '------')

            # create dataframe _DA_var_pert_df holding the results of the DA update
            # ---------------------------------------------------------------------
            self._DA_df(state=[ensemble_psi_valid, ensemble_sw_valid],
                        state_analysis=analysis_valid,
                        rejected_ens=rejected_ens)
                
            # export summary results of DA
            # ----------------------------------------------------------------
            
            meta_DA = {'listAssimilatedObs': list_assimilated_obs,
                       'listUpdatedparm':list_update_parm,
                       # '':,
                       # '':,
                       # '':,
                       }   
                
            self.backup_results_DA(meta_DA)
            
            # test = self.Archie
            
            self.backup_simu()
    
            # overwrite input files ensemble (perturbated variables)
            # ---------------------------------------------------------------------

            if self.count_atmbc_cycle<len(all_atmbc_times)-1: # -1 cause all_atmbc_times include TMAX
                # len(all_atmbc_times)
                self._update_input_ensemble(self.NENS, 
                                            all_atmbc_times, 
                                            self.dict_parm_pert,
                                            update_parm_list=list_update_parm, 
                                            analysis=analysis_valid
                                            )
            else:
                print('------ end of DA ------')
                pass   
            # print('------ end of update -------' + str(self.count_atmbc_cycle) + '/' + str(len(all_atmbc_times)-1) + '------')
            print('------ end of DA update -------' + str(self.count_DA_cycle) + '/' + str(len(ENS_times)) + '------')
            print('% of valid ensemble is: ' + str((len(self.ens_valid)*100)/(self.NENS)))
            # print(self.ens_valid)
            
            plt.close('all')



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
                           path_CATHY_folder= cwd,
                           # filename = os.path.join(cwd, 'output')
                           ) #TMAX=simu_time_max
        
        # df_psi = self.read_outputs(filename='psi',
        #                            path=os.path.join(cwd, 'output'))
        
        # df_sw = self.read_outputs(filename='sw',
        #                            path=os.path.join(cwd, 'output'))
                
        # for t in range(len(time_of_interest)):
        #     self._DA_df(state=[df_psi[t,:],df_sw[t,:]], 
        #                 t_ass=t, 
        #                 openLoop=True,
        #                 ens_nb=ens_nb)
            
        

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

        if verbose :
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
                "MAXVTKPRT": self.parm["NPRT"],
                
                # Is related to data assimilation (DA) (to do not considered if you do not have DA)
                "MAXENNOUT": 52560,
                "MAXVEG": 1,
            }

        # create dictionnary from kwargs
        for kk, value in kwargs.items():
            if kk in self.cathyH.keys():
                self.cathyH[kk] = value
                if verbose:
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
            CATHYH_file.write(
                "      PARAMETER (MAXVTKPRT={})\n".format(self.cathyH["MAXVTKPRT"])
                )
                
                
            # CATHYH_file.write("      PARAMETER (MAXVTKPRT=9)\n".format())
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
            if re.search("/t", tmp_param_value[i]) is None:
                if tmp_param_value[i].isdigit():
                    tmp_param_value[i] = int(tmp_param_value[i])
                else:
                    tmp_param_value[i] = float(tmp_param_value[i])
            else:
                list_tmp = tmp_param_value[i].split('/t')
                tmp_param_value[i] = [float(x) for x in list_tmp]
        
            self.hapin[hapin[i]] = tmp_param_value[i]


        if self.hapin["dr"] != self.hapin["delta_x"]:
            print("adapt rivulet param to dem resolution")
            self.hapin["dr"] = self.hapin["delta_x"]
            # sys.exit()

        L = []
        tmp_param_value_new = []
        for key, value in kwargs.items():
            self.hapin[key] = value

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
                                # print(key)
                                # print(value)
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
        Update DEM parameters

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

            # Write values first
            # -----------------------
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
                    print(str(list(self.dem_parameters.values())[counth]))
                    dem_parametersfile.write(
                        str(list(self.dem_parameters.values())
                            [counth]) + "\t" + "\n"
                    )
                    counth += 1

           # Write keys 
           # -----------------------
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
        
        Notes
        -----

        .. note::
            - create global zone_xyz variable

        .. note::
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
        
        if verbose:
            self.console.print(":black_nib: [b]update parm file [/b]")

        if len(self.parm) == 0:
            
            # set default parameters
            # --------------------------------------------------------------------
            self.parm = {
                "IPRT1": 3,  # Flag for output of input and coordinate data
                "NCOUT": 0,
                "TRAFLAG": 0, # Flag for transport
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
                
                # What is included into the output vtk file, CATHY default is pressure heads 1; 
                # to add saturation VTKF>=2
                "VTKF": 2, 
                "NPRT": 3,
                
                # Output times to record
                "(TIMPRT(I),I=1,NPRT)": [1800.0, 3600.0, 7200.0],
                "NUMVP": 1,
                "(NODVP(I),I=1,NUMVP)": [1], # should be positive (first node is 1)
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
                    
                    
                    
                    
        # check if consistency between times of interest and
        # number of times of interest
        # ------------------------------------------------------------
        if len(self.parm['(TIMPRT(I),I=1,NPRT)']) != self.parm["NPRT"]:
            warnings.warn('adjusting NPRT with respect to time of interests requested')
            self.parm["NPRT"] = len(self.parm['(TIMPRT(I),I=1,NPRT)'])     
                
        # check if consistency between DELTAT, DTMIN, and DTMAX
        # ------------------------------------------------------------        
        if (self.parm["DELTAT"] > min(np.diff(self.parm['(TIMPRT(I),I=1,NPRT)'])) or
            self.parm["DELTAT"] > 1e-2*min(np.diff(self.parm['(TIMPRT(I),I=1,NPRT)']))
            ):
           warnings.warn('adjusting DELTAT with respect to time of interests requested')
           self.parm["DELTAT"] = min(np.diff(self.parm['(TIMPRT(I),I=1,NPRT)']))/1e2
           
           
        if self.parm["DELTAT"] < self.parm["DTMIN"]:
           warnings.warn('adjusting 2*DTMIN == DELTAT: ' + str(self.parm["DTMIN"]))
           self.parm["DTMIN"] = self.parm["DELTAT"]/2
       
        if self.parm["DELTAT"] >= self.parm["DTMAX"]:
           warnings.warn('adjusting DTMAX == 2*DELTAT')
           self.parm["DTMAX"] = 2*self.parm["DELTAT"]   
           

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
            if self.count_DA_cycle is not None:
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
              are used to update the surface node values in PtimeP read in according to the
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
            if self.count_DA_cycle is not None:
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

    def update_atmbc(self, HSPATM=0, IETO=0, time=None, VALUE=[None, None],
                     show=False, verbose=False, **kwargs):
        '''
        Atmospheric forcing term (atmbc - IIN6)


        ..note:
        
        
                1 1                HSPATM,IETO
                0.0000000e+00      time
                5.5e-06              VALUE
                12.000000e+03      time
                0.00                 VALUE
                18.000000e+03      time
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
        time : TYPE, optional
            DESCRIPTION. The default is None.
        VALUE : TYPE, optional
            DESCRIPTION. The default is [Precipitation, EvapoTranspiration].
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
            v_atmbc = VALUE[1] - abs(VALUE[0])
        else:
            v_atmbc = VALUE

        # set default parameters
        # --------------------------------------------------------------------

        self.atmbc = {"HSPATM": HSPATM, "IETO": IETO,
            "time": time, "VALUE": VALUE}

        # len(self.atmbc['time'])
        # overwrite existing input atmbc file
        # --------------------------------------------------------------------

        if 'filename' in kwargs:
            filename = kwargs['filename']
        else:
            filename = os.path.join(
                self.workdir, self.project_name, self.input_dirname, "atmbc")


        backup = True
        if 'backup' in kwargs:
            backup = kwargs['backup']
            
        if backup == True:
            if self.count_DA_cycle is not None:
                dst_dir = filename + str(self.count_DA_cycle-1)
                shutil.copy(filename,dst_dir) 
               
                   
        with open(filename, "w+") as atmbcfile:
            atmbcfile.write(str(HSPATM) + "\t" + str(IETO) +
                            "\t" + "HSPATM" + "\t" + "IETO" + "\n"
                            )

            # if v_atmbc is a scalar meaning that atmbc are homoegenous
            # -----------------------------------------------------------------
            for t, v in zip(time, v_atmbc):

                if verbose == True:
                    print(t, v)
                atmbcfile.write("{:.3e}".format(t) + "\t" + "time" + "\n")
                # atmbcfile.close()

                if isinstance(v, float) | isinstance(v, int):
                    atmbcfile.write("{:.3e}".format(v) + "\t" + "VALUE" + "\n")
                    # atmbcfile.write(str(v) + "\t" + "VALUE" + "\n")
                else:
                    # atmbcfile.write(str(v) + "\t" + "VALUE" + "\n")  
                    np.savetxt(atmbcfile, v, fmt="%1.3e")

                    
        atmbcfile.close()

        self.update_parm(TIMPRTi=list(time),NPRT=len(time),TMAX=max(time))

        # don't need to update if sequential DA as the cathy.exe is already created
        # ---------------------------------------------------------------------
        if 'omit_cathyH' not in kwargs:
            self.update_cathyH(MAXPRT=len(time))

        if show == True:
            # if HSPATM !=0:
            #     print('impossible to plot for non homogeneous atmbc')
            #     # sys.exit()
            # else:
            # x_units = "sec"
            # for key, value in kwargs.items():
            #     if key == "x_units":
            #         x_units = value
            plt, ax = plt_CT.show_atmbc(time, VALUE, 
                                        **kwargs)

            return plt, ax
        
        pass
    

              

    def init_boundary_conditions(self, times):
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
            When the surface node reaches a threshold level of saturation or moisture deficit, 
            the boundary condition is switched to a Dirichlet (specified head) condition, 
            and the infiltration or exfiltration process becomes soil limited [1].

        Returns
        -------
        None.

        '''
        # self.update_nansfdirbc()
        # self.update_nansfneubc()
        # self.update_sfbc()
        
        try:
            self.console.print('init boundary condition dataframe')
            self.create_mesh_bounds_df(self.grid3d['nodes_idxyz'],times)
            plt_CT.plot_mesh_bounds(self.mesh_bound_cond_df)
            self.create_mesh_vtk()
        except:
            raise ValueError('cannot init boundary conditions dataframe')
            pass
        
        pass
        
    def check_for_inconsistent_BC(self):
        pass
    
    def set_BC_laterals(self, time, BC_type='', val=0):
        ''' 
        Set all sides expect surface one        
        '''
        
        nnodes = len(self.mesh_bound_cond_df[self.mesh_bound_cond_df['time (s)']==0]['id_node'])
        for tt in time:
            BC_bool_name = []
            BC_bool_val = []
        
            for id_node in range(nnodes):
                if self.mesh_bound_cond_df['bound'].loc[int(id_node)] == True:
                    BC_bool_name.append(BC_type)
                    BC_bool_val.append(0)
                else:
                    BC_bool_name.append(None)
                    BC_bool_val.append(val)

            self.update_mesh_boundary_cond(
                                            time = tt,
                                            BC_bool_name=BC_bool_name,
                                            BC_bool_val=BC_bool_val
                                           )
                
        self.check_for_inconsistent_BC()

        #self.update_mesh_vtk(BC_bool_name,BC_bool_val)
        
        pass
        
        
        

    def update_nansfdirbc(self,time=[],NDIR=0,NDIRC=0,NQ3=None,no_flow=False,
                            bound_xyz=None,
                            pressure_head=[],
                        ):
        '''
        Dirichlet Boundary conditions (or specified pressure) at time t
       
        - To simulate the no-flow boundaries conditions for the bottom and 
          vertical sides of the domain it is necessary to set NDIR and NDIRC 
          equal to zero. 
        - To simulate different boundary conditions, it is necessary to 
          indicate the number of selected nodes through NDIR or NDIRC, 
          then to specify the node ID’s that you want to consider and
          eventually the value of pressure head or flux that you want to assign.


        ..note::
            update_nansfdirbc use the grid3d to refer to mesh nodes
            
        Parameters
        ----------
        time : np.array([]), optional
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

        # check that mesh_bound_cond_df exist
        # --------------------------------------------------------------------
        if hasattr(self,'mesh_bound_cond_df') is False:
            self.init_boundary_conditions(time)
        
        # apply BC
        # --------------------------------------------------------------------
        if no_flow: # Dirichlet  == 0    
            print('shortcut set_BC_laterals mesh dataframe')          
            # self.set_BC_laterals(time=time, BC_type='Dirichlet', val=0)
        else:
            raise ValueError('Non homogeneous Dirichlet Boundary conditions Not yet implemented')
            

        self.update_parm()
        self.update_cathyH()
        
        
        with open(os.path.join(self.workdir, 
                                self.project_name, 
                                self.input_dirname,
                                "nansfdirbc"
                                ),"w+") as nansfdirbcfile:

            # self.mesh_bound_cond_df
            

            # To simulate the no-flow boundaries conditions for the bottom and 
            # vertical sides of the domain --> NDIR and NDIRC equal to zero
            #-------------------------------------------------------------
            if no_flow: # Dirichlet
                if len(time) == 0:
                    time = self.atmbc["time"]
                for tt in time:
                    # nansfdirbcfile.write(str(tt) + "\t" + "time" + "\n")
                    nansfdirbcfile.write("{:.3e}".format(tt) + "\t" + "time" + "\n")

                    nansfdirbcfile.write(str(NDIR) + "\t"+ str(NDIRC) + "\t" + "NDIR" + "\t" + "NDIRC" + "\n")
            else:
                raise ValueError('Non homogeneous Dirichlet Boundary conditions Not yet implemented')
                    
        nansfdirbcfile.close()
                
        # read existing input nansfdirbc file
        # --------------------------------------------------------------------
        # with open(os.path.join(self.workdir, 
        #                        self.project_name, 
        #                        self.input_dirname,
        #                        "nansfdirbc"
        #                        ),"w+") as nansfdirbcfile:


        #     # To simulate the no-flow boundaries conditions for the bottom and 
        #     # vertical sides of the domain --> NDIR and NDIRC equal to zero
        #     #-------------------------------------------------------------
        #     if noflow: # Dirichlet
        #         if len(time) == 0:
        #             time = self.atmbc["time"]
        #         for tt in time:
        #             # nansfdirbcfile.write(str(tt) + "\t" + "time" + "\n")
        #             nansfdirbcfile.write("{:.3e}".format(tt) + "\t" + "time" + "\n")

        #             nansfdirbcfile.write(str(NDIR) + "\t"+ str(NDIRC) + "\t" + "NDIR" + "\t" + "NDIRC" + "\n")
        #     else:
        #         raise ValueError('Non homogeneous Dirichlet Boundary conditions Not yet implemented')
                    
        # nansfdirbcfile.close()
        
        self.update_parm()
        self.update_cathyH()

                

                # len(bool_Dirichlet)
                # len(self.mesh_bound_cond_df['id_node'])
                
                
            # elif bound_xyz is not None:
            #     for tt in time:
            #         #Se avessi delle variazioni dovrei indicare il nodo ed il valore di pressione
            #         nansfdirbcfile.write(str(tt) + "\t" + 'time' + "\n")
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
                #       write(33,*) 0.0, 'time'
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
    
                #       write(33,*) 2e+20, 'time'
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

            
        
        
        
        
        

        pass

    def update_nansfneubc(self, time=[], NQ=0, ZERO=0, 
                          fixed_pressure=False, 
                          no_flow=False):
        '''
        Neumann boundary conditions (or specifed flux) at time t


        Parameters
        ----------
        time : np.array([]), optional
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
        
        # check that the mesh exist
        # --------------------------------------------------------------------
        try:
            self.grid3d = in_CT.read_grid3d(self.project_name)
        except OSError:
            print("grid3d missing - need to run the processor with IPRT1=2 first")

        # check that mesh_bound_cond_df exist
        # --------------------------------------------------------------------
        if hasattr(self,'mesh_bound_cond_df') is False:
            self.init_boundary_conditions(time)
        
        # apply BC
        # --------------------------------------------------------------------
        if no_flow: # Neumanm               
            # self.set_BC_laterals(time=time, BC_type='Neumanm', val=0)
            print('shortcut set_BC_laterals mesh dataframe')          

        else:
            raise ValueError('Non homogeneous Neumanm Boundary conditions Not yet implemented')
            
            
        # read existing input nansfneubc file
        # --------------------------------------------------------------------

        with open(os.path.join(self.workdir, 
                               self.project_name, 
                               self.input_dirname, 
                               "nansfneubc"),"w+") as nansfneubcfile:

            if len(time) == 0:
                time = self.atmbc["time"]
            for tt in time:
                # nansfneubcfile.write(str(tt) + "\t" + "time" + "\n")
                nansfneubcfile.write("{:.3e}".format(tt) + "\t" + "time" + "\n")

                nansfneubcfile.write(
                    str(ZERO) + "\t" + str(NQ) + "\t" +
                        "ZERO" + "\t" + "NQ" + "\n"
                )

        nansfneubcfile.close()

        pass

    def update_sfbc(self, time=[], sfbc = False, no_flow=False):
        '''
        Seepage face boundary conditions at time t

        Parameters
        ----------
        time : np.array([]), optional
            Absolute simulation time in sec. 
            The default is [].

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

        # check that mesh_bound_cond_df exist
        # --------------------------------------------------------------------
        if hasattr(self,'mesh_bound_cond_df') is False:
            self.init_boundary_conditions(time)
        
        # apply BC
        # --------------------------------------------------------------------
        if no_flow:              
            # self.set_BC_laterals(time=time, BC_type='sfbc', val=0)
            print('shortcut set_BC_laterals mesh dataframe')          
            
        else:
            raise ValueError('Non homogeneous Neumanm Boundary conditions Not yet implemented')
            
            
            
        with open(
            os.path.join(self.workdir, self.project_name,
                         self.input_dirname, "sfbc"),
            "w+",
        ) as sfbcfile:

            if len(time) == 0:
                time = self.atmbc["time"]
            for tt in time:
                sfbcfile.write("{:.3e}".format(tt) + "\n")
                sfbcfile.write("0" + "\n")

        sfbcfile.close()

        pass

    
    def create_mesh_bounds_df(self,grid3d,times):
        '''
        Create a dataframe with flag for different boundary condtions assigned to each nodes

        Parameters
        ----------
        grid3d : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

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
        
        print('SKip time dependence init boundary condition dataframe - consequences (?)')
        # self.mesh_bound_cond_df = pd.concat([self.mesh_bound_cond_df]*len(times), ignore_index=True)

        # times_reshaped = [val for val in times for _ in np.arange(len(grid3d))]
        # self.mesh_bound_cond_df['time (s)'] = times_reshaped
        # empty_vec = np.empty(len(times_reshaped))
        # empty_vec[:] = np.nan
        # self.mesh_bound_cond_df['BC_type'] = empty_vec #np.nan(len(times_reshaped))
        # self.mesh_bound_cond_df['BC_type_val'] = empty_vec
        
        pass
        

    def update_mesh_boundary_cond(self, 
                                  time,
                                  BC_bool_name=[],
                                  BC_bool_val=[],):
        '''
        update_mesh_bounds

        Parameters
        ----------
        bound_type : str
            Neumann or Dirichlet.
        bound_bool : TYPE
            Boolean for bound cond.
        '''
        self.console.print(":sponge: [b]update boundary condition dataframe[/b]")
        self.mesh_bound_cond_df.loc[self.mesh_bound_cond_df['time (s)']==time, 'BC_type'] = BC_bool_name
        self.mesh_bound_cond_df.loc[self.mesh_bound_cond_df['time (s)']==time, 'BC_type_val'] = BC_bool_val
        pass
            

    def create_mesh_vtk(self,topo=None):
        '''
        Create custum mesh

        Returns
        -------
        None.

        '''
        
        

        
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

        

    def update_mesh_vtk(self, prop='', prop_value=[], savevtk=True, **kwargs):
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
        
        meshname = self.project_name
        if 'meshname' in kwargs:
            meshname = kwargs['meshname']
        
        self.mesh_pv_attributes.add_field_data(prop_value, prop)
        self.mesh_pv_attributes.save(os.path.join(self.workdir, 
                               self.project_name, 
                               'vtk/',
                               meshname +
                               '.vtk'), binary=False)
        
        pass
            
            
            
                    
    def update_soil(self, IVGHU=[], FP=[], SPP=[], soil_het_dim=1,
                    verbose=False, show=False, **kwargs):
        '''
        Soil parameters (soil - IIN4). The porous media properties.

        ..note::
            The first thing that must be decides is the type of relationship to describe the hydraulic
            characteristics of the unsaturated soil (i.e. retention curves). This can be done through
            the choice of the parameter **IVGHU** amongst the several options.


        Parameters
        ----------
        IVGHU : int, optional
            = -1 table look up for moisture curves
            = 0 for van Genuchten moisture curves
            = 1 for extended van Genuchten moisture curves
            = 2 for moisture curves from Huyakorn et al (WRR 20(8) 1984, WRR 22(13) 1986) with Kr=Se**n conductivity relationship
            = 3 for moisture curves from Huyakorn et al (WRR 20(8) 1984, WRR 22(13) 1986)with conductivity relationship from Table 3 of 1984 paper (log_10 Kr(Se) curve)
            = 4 for Brooks‐Corey moisture curves.

            The default is [].

        FP : list, optional
            Feddes Parameters. The default is [].
            [PCANA PCREF PCWLT ZROOT PZ OMGC]
            - 'PCANA': float, anaerobiosis point --> h2 (≈ -0.5m)
            - 'PCREF': float, field capacity --> h3 (≈ 4-m)
            - 'PCWLT': float, wilting point --> h4 (≈ -150m)
            - 'ZROOT': float, root depth
            - 'PZ': float, pz is an empirical parameter
            - 'OMGC': float, 0<OMGC<1 Compensatory mechanisms for root water uptake

            .. note:
                Calculate actual transpiration Ta (m d−1). 
                A multiplicative reduction factor is defined by four pressure heads (0 > h1 > h2 > h3 > h4) 
                delimiting five phases of uptake
                For details, see http://dx.doi.org/10.1002/2015WR017139
                BETA = (1-DEPTH/ZROOT(VEG_TYPE(I)))*DEXP(-1.0d0*PZ(VEG_TYPE(I))*DEPTH/ZROOT(VEG_TYPE(I)))
                BTRANI(K) = MAX(ZERO,BETA*DZ*GX)
                BTRAN(I)  = BTRAN(I) + BETA*DZ
                OMG(I)    = OMG(I) + GX*BETA*DZ
                QTRANIE(K)=ETP(I)*BTRANI(K)/BTRAN(I)/MAX(OMG(I),OMGC(VEG_TYPE(I)))

        SPP : list, optional
            Soil Physical Properties. The default is [].
            - 'PERMX' (NSTR, NZONE): saturated hydraulic conductivity - xx
            - 'PERMY' (NSTR, NZONE): saturated hydraulic conductivity - yy
            - 'PERMZ' (NSTR, NZONE): saturated hydraulic conductivity - zz
            - 'ELSTOR' (NSTR, NZONE): specific storage
            - 'POROS'  (NSTR, NZONE): porosity (moisture content at saturation) = \thetaS

            retention curves parameters VGN, VGRMC, and VGPSAT
            - 'VGNCELL' (NSTR, NZONE): van Genuchten curve exponent  = n
            - 'VGRMCCELL' (NSTR, NZONE): residual moisture content = \thetaR
            - 'VGPSATCELL' (NSTR, NZONE): van Genuchten curve exponent --> 
                                          VGPSAT == -1/alpha (with alpha expressed in [L-1]);
                
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
        # if hasattr(self, 'soil_SPP') == 0:
            SPP = self.set_SOIL_defaults(SPP_default=True)
        if len(FP) == 0:
        # if hasattr(self, 'soil_FP') == 0:
            FP = self.set_SOIL_defaults(FP_default=True)
            

        # check size of soil properties 
        # --------------------------------------------------------------------
        if isinstance(SPP['PERMX'], float):
            if self.dem_parameters["nzone"]!=1:
                raise ValueError
        else:
            if len(SPP['PERMX'])!=self.dem_parameters["nzone"]:
                raise ValueError("Wrong number of zones: PERMX size is " + str(len(SPP['PERMX'])) 
                                 + 'while nzone is ' + str(self.dem_parameters["nzone"]))
                    
        # read function arguments kwargs and udpate soil and parm files
        # --------------------------------------------------------------------
        for keykwargs, value in kwargs.items():
            self.soil[keykwargs] = value
            self.parm[keykwargs] = value

        # loop over Feddes parameters FP mandatory fct arg
        # --------------------------------------------------------------------
        for fp in FP:  # loop over fedded parameterssoil_het_dim
            self.soil[fp] = FP[fp]

        # check consistency between parameters
        # --------------------------------------------------------------------
        if SPP['VGRMCCELL']>=SPP['POROS']:
            raise ValueError("residual water content is" + str(SPP['VGRMCCELL']) 
                             + '> porosity ' + str(SPP['POROS']))
        
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
        
        if show:
            update_map_veg = self.map_prop_veg(FP)
            fig, ax = plt_CT.dem_plot_2d_top(update_map_veg,
                                              label='all')
            fig.savefig(os.path.join(self.workdir,self.project_name,'map_veg.png'), dpi=400)
            

        # write soil file
        # --------------------------------------------------------------------
        self._write_SOIL_file(SoilPhysProp, FeddesParam, **kwargs)
        
        # if show:
        #     plt.savefig(os.path.join(self.workdir,self.project_name,'map_veg.png'), dpi=400)
        #     plt.close()  
            
        #     raise NotImplementedError('show for soil not yet implemented')
            
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
            

            SoilPhysProp = []


            # loop over strates
            # -----------------------------------------------------------------
            for istr in range(self.dem_parameters["nstr"]):
                izoneSoil = np.zeros([self.dem_parameters["nzone"], 8])
    
                #  loop over zones (defined in the zone file)
                # --------------------------------------------------------------
                for izone in range(self.dem_parameters["nzone"]):                  
                    for i, spp in enumerate(SPP):
                        izoneSoil[izone,i] = SPP[spp][izone]
                        
                SoilPhysProp.append(izoneSoil)
            SoilPhysProp = np.vstack(SoilPhysProp)
                
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
            if len(FP[list(FP)[0]])==1:
                self.update_veg_map()
            else:
                raise ValueError('Found multiple values of Feddes zones' +
                                 'but vegetation map is not defined')
       
        
        # Check vegetation heterogeneity dimension
        # ----------------------------------------
        if self.cathyH["MAXVEG"] != len(FP[list(FP)[0]]):
            raise ValueError("Wrong number of vegetations: PCANA size is " + str(len(FP[list(FP)[0]])) 
                             + 'while MAXVEG is ' + str(self.cathyH["MAXVEG"]))

        # check number of vegetation
        # --------------------------------------------------------------------
        if self.cathyH["MAXVEG"] > 1:
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
        '''

        # backup file during DA scheme cycle
        # --------------------------------------------------------------------
        backup = True
        if 'backup' in kwargs:
            backup = kwargs['backup']
            

            
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
            
        if 'path' in kwargs:
            soil_filepath =  os.path.join(kwargs['path'], "soil")

        if backup:
            if self.count_DA_cycle is not None:
               dst_dir = soil_filepath + str(self.count_DA_cycle)
               shutil.copy(soil_filepath,dst_dir) 
           

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

        #print("update root map")
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
                np.savetxt(rootmapfile, indice_veg, fmt="%1.2e")

        rootmapfile.close()

        # self.update_zone(indice_veg)
        self.update_cathyH(MAXVEG=len(np.unique(indice_veg)))
        # ,                        MAXZON=len(np.unique(indice_veg))
        
        self.veg_map = indice_veg


        if show is not None:
            plt, ax = plt_CT.indice_veg_plot(indice_veg,**kwargs)
            
            return indice_veg, plt, ax
        

        return indice_veg

    # -------------------------------------------------------------------#
    # %% DATA ASSIMILATION FCTS

    def _DA_init(self, NENS=[], ENS_times=[], parm_pert=[],update_parm_list='all'):
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

        self.NENS = NENS  # THIS IS TEMPORARY
        # updated_parm_list - update_parm_list


        if hasattr(self,'ens_valid') is False:
            self.ens_valid = list(np.arange(0,self.NENS))
                
        # create sub directories for each ensemble
        # ---------------------------------------------------------------------
        self._create_subfolders_ensemble(NENS)
        
        
        # update list of updated parameters based on problem heterogeneity
        # ---------------------------------------------------------------------
        
        updated_parm_list = ['St. var.']
        for d in self.dict_parm_pert.keys():
            for l in update_parm_list:
                if l in d:
                    updated_parm_list.append(d)
        
        return updated_parm_list


    def _check_before_analysis(self, ens_valid=[],
                               threshold_rejected=10):
        '''
        Filter is applied only on selected ensemble

        Parameters
        ----------
        update_key : TYPE
            DESCRIPTION.

        Returns
        -------
        rejected_ens

        '''      
        
        self.console.print(":white_check_mark: [b]check scenarii before analysis[/b]")
            
        ###CHECK which scenarios are OK and which are to be rejected and then create and applied the filter 
        ###only on the selected scenarios which are to be considered.
                   
        rejected_ens_new = []
        for n in range(self.NENS):
            
            if n in ens_valid:
                # print('ensemble_nb ' + str(n) + ' before update analysis')
                try:
                    df_mbeconv = out_CT.read_mbeconv(os.path.join(self.workdir, self.project_name,
                                                       "DA_Ensemble/cathy_" + str(n+1),
                                                       'output/mbeconv'))
                    if df_mbeconv['CUM.'].isnull().values.any():
                        rejected_ens_new.append(True)
                        self.console.print(":x: [b]df_mbeconv['CUM.'][/b]" + str(df_mbeconv['CUM.'].isnull().values.any())+ ', ens_nb:' + str(n+1))
        
                    elif (np.round(df_mbeconv['TIME'].iloc[-1]))<np.round(self.parm['(TIMPRT(I),I=1,NPRT)'][-1]) - self.parm['DELTAT']:
                        rejected_ens_new.append(True)
                        self.console.print(":x: [b]df_mbeconv['time'][/b]" + str(df_mbeconv['TIME'].iloc[-1]) + ', ens_nb:' + str(n+1))
        
                    elif len(df_mbeconv) == 0:
                        rejected_ens_new.append(True)
                        self.console.print(":x: [b]len(df_mbeconv)['TIME'][/b]" + str(len(df_mbeconv))+ ', ens_nb:' + str(n+1))
        
                    else:
                        rejected_ens_new.append(False)                  
                    
                except:
                    rejected_ens_new.append(True)
                    print('cannot read mbeconv')
                    pass
                    # UnboundLocalError: local variable 'df_mbeconv' referenced before assignment

            else:
                rejected_ens_new.append(True)

        #check whether the number of discarded scenarios is major than the 10% of the total number [N]; 
        #if the answer is yes than the execution stops
        # --------------------------------------------------------------------
        if sum(rejected_ens_new)>=(threshold_rejected/100)*self.NENS:
            print('new parameters (if updated)')
            print(self.dict_parm_pert)

            raise ValueError('% number of rejected ensemble is too high:' + str((sum(rejected_ens_new)*100)/(self.NENS)))
        
        return  rejected_ens_new

        

        
    def _check_after_analysis(self,update_key,list_update_parm):
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        CHECK which scenarios parameters are OK and which one to discard

        Returns
        -------
        None.

        '''
        self.console.print(":white_check_mark: [b]check scenarii post update[/b]")

        id_valid = [list(np.arange(0,self.NENS))]
        id_nonvalid = []
    
    
        def test_negative_values(update_key,pp):
            id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]<0)[0]))
            id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]>0)[0]))
            
            for i in id_nonvalid:
                self.console.print(":x: [b]negative" + str(pp) + ":[/b]" + 
                                   str(self.dict_parm_pert[pp[1]][update_key][i])+ 
                                   ', ens_nb:' + str(i))
            return id_nonvalid, id_valid
                          
        def test_range_values(update_key,pp,min_r,max_r):
            id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]<min_r)[0]))
            id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]>max_r)[0]))
            id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]>min_r)[0]))
            id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]<max_r)[0]))
            for i in id_nonvalid:
                self.console.print(":x: [b]out of range" + str(pp) + ":[/b]" + 
                                   str(self.dict_parm_pert[pp[1]][update_key][i])+ 
                                   ', ens_nb:' + str(i))
            return id_nonvalid, id_valid
        
            
        
        for pp in enumerate(list_update_parm[:]):
            if 'St. var.' in pp[1]:
                pass
            else:
                if 'rFluid'.casefold() in pp[1].casefold():
                    id_nonvalid, id_valid = test_negative_values(update_key,pp)
                elif 'porosity'.casefold() in pp[1].casefold():
                    id_nonvalid, id_valid = test_negative_values(update_key,pp)
                elif 'a_Archie'.casefold() in pp[1].casefold():
                    id_nonvalid, id_valid = test_range_values(update_key,pp, 0, 3)
                elif 'Ks'.casefold() in pp[1].casefold():
                                                                  
                    id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]<0)[0]))
                    id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]>0)[0]))
                    for i in id_nonvalid:
                        self.console.print(":x: [b]negative Ks:[/b]" + 
                                           str(self.dict_parm_pert[pp[1]][update_key][i])+ 
                                           ', ens_nb:' + str(i))
                elif 'PCREF'.casefold() in pp[1].casefold():
                    id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]>0)[0]))
                    id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]<0)[0]))
                    for i in id_nonvalid:
                        self.console.print('Positive values encountered in PCREF:' + 
                              str(self.dict_parm_pert[pp[1]][update_key][i])+ 
                                           ', ens_nb:' + str(i))
                # ckeck if new value of Zroot is feasible
                # ------------------------------------------------------------
                elif 'Zroot'.casefold() in pp[1].casefold():
                    id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]> 
                                                abs(min(self.grid3d['nodes_idxyz'][:, -1])))[0]))
                    id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]< 
                                             abs(min(self.grid3d['nodes_idxyz'][:, -1])))[0])
                                    )
                    for i in id_nonvalid:
                        self.console.print(":x: [b]unfeasible root depth:[/b]" + 
                                           str(self.dict_parm_pert[pp[1]][update_key][i])+ ', ens_nb:' + str(i))
    
        id_nonvalid_flat = [item for sublist in id_nonvalid for item in sublist]
        id_valid = set.difference(set(list(np.arange(0,self.NENS))), set(list(np.unique(id_nonvalid_flat))))
        return id_valid


    def spatialize_parameters(self,list_update_parm,parm):
        '''
        Extend the number of parameters according to the number of zones
        (useful for the heteregeneous case)
        '''
        i=0
        parm_extended = []
        for l in list_update_parm:
            print(l)
            if 'St. var.' in l:
                pass
            else:
                nb_of_zones = 1
                if 'surf_zones_param' in self.dict_parm_pert[l]:
                    nb_of_zones =  self.dict_parm_pert[l]['surf_zones_param']
                    print(nb_of_zones)
                    parm_extended.append(np.tile(parm[i,:],(nb_of_zones,1)))
                i = i + 1
        parm_extended = np.vstack(parm_extended)
        return parm_extended
        

    def transform_parameters(self,list_update_parm,
                             update_key=None,
                             back=False,
                             param=None):
        '''
        Parameter (back) transformation before analysis (log, bounded log, ...)

        ..note:: parameters are log transform so that they become Gaussian distributed 
                 to be updated by the Ensemble Kalman Filter
        '''
        param_new = []
        for pp in enumerate(list_update_parm[:]):
            param_trans = []

            if 'St. var.' in pp[1]:
                pass 
            else:           
                if update_key:
                    param = self.dict_parm_pert[pp[1]][update_key]
                if self.dict_parm_pert[pp[1]]['transf_type'] == 'log':
                    print('log transformation')
                    if back:
                        param_trans = np.exp(param)
                    else:
                        param_trans = np.log(param)
                elif self.dict_parm_pert[pp[1]]['transf_type'] == 'log-ratio':
                   print('bounded log transformation')
                   range_tr = self.dict_parm_pert[pp[1]]['transf_bounds']
                   A = range_tr['A'] # min range
                   B = range_tr['B'] # max range
                   if back:
                       param_trans = (
                                       np.exp(  (param - A)/
                                                (B - param)
                                              )
                                    )
                   else:
                       param_trans = (
                                       np.log(  (param - A)/
                                                (B - param)
                                              )
                                    )                  
                elif self.dict_parm_pert[pp[1]]['transf_type'] == 'hyperbolic':
                   print('hyperbolic transformation')
                   # Y = sinh−1(U )
                   param_trans = np.log(self.dict_parm_pert[pp[1]][update_key])
                   
                else:
                    param_trans = param
                param_new.append(param_trans)
            if len(param_new)>0:
                np.vstack(param_new)
                
        return np.array(param_new)
                
    
                                

    def _get_data2assimilate(self,list_assimilated_obs,time_ass=None,match=False):
        '''
        Loop over observation dictionnary and select corresponding data for a given assimilation time

        Parameters
        ----------
        list_assimilated_obs : list
            list of observation to assimilate.

        Returns
        -------
        data : list
            list of data values.
        '''
        
        # find the method to map for this time step
        # --------------------------------------------------------
        obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)
        
        def extract_data(sensor, time_ass, data):
            data2add = []
            if 'ERT' in sensor:
                if 'pygimli' in items_dict[time_ass][1][sensor]['data_format']:
                # if 'pygimli' in obs2map[time_ass]['data_format']:
                # if 'pygimli' in self.dict_obs[time_ass]['ERT']['data_format']:
                    data2add.append(np.array(items_dict[time_ass][1][sensor]['data']['rhoa']))
                else:
                    data2add.append(items_dict[time_ass][1][sensor]['data']['resist'].to_numpy())
                data2add= np.hstack(data2add)
            else:
                data2add.append(items_dict[time_ass][1][sensor]['data'])
                
            return data2add
                    
       
        if time_ass is None:
            time_ass = self.count_DA_cycle

        data = []
        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------       
        items_dict = list(self.dict_obs.items())
        for sensor in items_dict[time_ass][1].keys():  
            if 'all' in list_assimilated_obs:
                # obskey2map.append(sensor)
                data2add = extract_data(sensor, time_ass, data)
                data.append(data2add)
            else:
                for l in list_assimilated_obs:
                    if match:
                        if l == sensor:     
                            # obskey2map.append(sensor)
                            data2add = extract_data(sensor, time_ass, data)
                            data.append(data2add)
                    else:
                        if l in sensor:  
                            # obskey2map.append(sensor)
                            data2add = extract_data(sensor, time_ass, data)
                            data.append(data2add)
            # len(data2add)      
        return np.hstack(data), obskey2map
        
        
    def _DA_analysis(self, prediction, 
                     DA_type='enkf_Evensen2009_Sakov',
                     list_update_parm=['St. var.'],
                     list_assimilated_obs='all',
                     ens_valid=[]):
        """
        THIS SHOULD BE MOVED TO DA CLASS

        Analysis ensemble using DA
        print(list_assimilated_obs)


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
        rejected_ens : list
            list of ensemble to reject
        Returns
        -------
        None.

        Note: for now only implemented for EnkF DA
        """

        update_key = 'ini_perturbation'
        if self.count_DA_cycle > 0:
            update_key = 'update_nb' + str(self.count_DA_cycle)

        # find the method to map for this time step
        # --------------------------------------------------------
        obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)

        # prepare states
        # ---------------------------------------------------------------------
        ensemble_psi, ensemble_sw, ens_size, sim_size = self._read_state_ensemble()

        # ---------------------------------------------------------------------
        # prepare data
        # ---------------------------------------------------------------------
        # select ONLY data to assimilate
        # --------------------------------------------------------------------- 
        data, _ = self._get_data2assimilate(list_assimilated_obs)
        print('data size: ' + str(len(data)))

        # check size of data_cov
        # ---------------------------------------------------------------------
        # if len(self.stacked_data_cov)/len(self.dict_obs.items()) != len(data):
        # if 'ERT' in sensors:
        #     if np.shape(self.stacked_data_cov[self.count_DA_cycle])[0] != len(data[0]):
        #         raise ValueError('need to compute data covariance')
        # else:
        if np.shape(self.stacked_data_cov[self.count_DA_cycle])[0] != len(data):
            raise ValueError('need to compute data covariance')
                
        data_cov =self.stacked_data_cov[self.count_DA_cycle]
    
        
        
        # # filter non-valid ensemble before Analysis
        # # ---------------------------------------------------------------------
        prediction_valid = prediction
        ensemble_psi_valid = ensemble_psi[:,ens_valid]            
        ensemble_sw_valid = ensemble_sw[:,ens_valid]            
        # print('prediction valid size: ' + str(len(prediction_valid[:,0])))

        
        
        # prepare parameters
        # ---------------------------------------------------------------------
        # When updating the states only, the elements of X are the
        # pressure heads at each node of the finite element grid, while
        # the state augmentation technique is used when also updat-
        # ing the parameters
        param = self.transform_parameters(list_update_parm,update_key)
        print('parm size: ' + str(len(param)))

        # if update_key == 'ini_perturbation':
        #     param = self.spatialize_parameters(list_update_parm,param)

        param_valid = []
        if len(param)>0:          
            param_valid = param[:, ens_valid]  

            
        # run Analysis
        # ---------------------------------------------------------------------
                
        result_analysis = cathy_DA.run_analysis(DA_type,
                                                data,
                                                data_cov,
                                                param_valid,
                                                list_update_parm,
                                                [ensemble_psi_valid,ensemble_sw_valid],
                                                prediction_valid,
                                                alpha=self.damping,
                                                )
                     
            
        # plot ensemble covariance matrices and changes (only for ENKF)
        # ---------------------------------------------------------------------
        if len(result_analysis)>2:
            [A, Amean, dA, 
             dD, MeasAvg, S, 
             COV, B, dAS, 
             analysis, 
             analysis_param]  = result_analysis
            
            try:
                plt_CT.show_DA_process_ens(ensemble_psi_valid,
                                            data,COV,
                                            dD,dAS,B,analysis, 
                                            savefig=True, 
                                            savename= os.path.join(self.workdir,
                                                                   self.project_name,
                                                                   'DA_Ensemble',
                                                                   'DA_Matrices_t' 
                                                                   +str(self.count_DA_cycle)
                                                                  ),
                                            # label_sensor=str([*self.dict_obs[0]])
                                            label_sensor=str([*obskey2map])
                                            ) 
            except:
                self.console.print("[b]Impossible to plot matrices[/b]")
                
                
        # particule filter
        # ----------------
        else:
            [analysis, 
             analysis_param]  = result_analysis
            

        # Back-transformation of the parameters
        # ---------------------------------------------------------------------
        # for pp in enumerate(list_update_parm[:]):
        #     if 'St. var.' in pp[1]:
        #         pass 
        #     else:  
            
        # analysis_param_valid = []
        # if len(param)>0:          
        #     analysis_param_valid = analysis_param[ens_valid]  
            
        self.transform_parameters(list_update_parm,
                                  param=np.hstack(analysis_param),
                                  back=True)
                            
                # if self.dict_parm_pert[pp[1]]['transf_type'] == 'log':
                #     print('back log transformation')
                #     analysis_param = np.exp(analysis_param)

                    
                    
        # return ensemble_psi, Analysis, AnalysisParam, Data, data_cov
        # ---------------------------------------------------------------------
        return(ensemble_psi, ensemble_sw, data,  analysis, 
                 analysis_param)



                            
    def _obs_key_select(self,list_assimilated_obs):
        ''' find data to map from dictionnary of observations '''
        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------
        obskey2map = [] # data type 2 map
        obs2map = []    # data dict 2 map            
        items_dict = list(self.dict_obs.items())
        for sensor in items_dict[self.count_DA_cycle][1].keys():  
            if 'all' in list_assimilated_obs:
                obskey2map.append(sensor)
                obs2map.append(items_dict[self.count_DA_cycle][1][sensor])
            else:
                for l in list_assimilated_obs:
                    if l in sensor:           
                        obskey2map.append(sensor)
                        obs2map.append(items_dict[self.count_DA_cycle][1][sensor])
        return obskey2map,obs2map    

            
    
    
    def map_states2Observations(self,
                                 list_assimilated_obs='all', 
                                 parallel=False,
                                 default_state = 'psi',
                                 verbose = False,
                                 **kwargs):
        '''
        Translate (map) the state values (pressure head or saturation) to observations (or predicted) values
        
        
        In other words: apply H mapping operator to convert model to predicted value Hx
        
        In practice, for each assimilation time, loop over the observation 
        dictionnary and infer observation type. Stack Hx.
        Hx = [Hx0,Hx1,Hx_nb_of_observation]

        .. note::
        
            - For ERT data, H is Archie law, SWC --> ER0 --> ERapp
            - For SWC data, H is a multiplicator (porosity)
            - For tensiometer data, no mapping needed or VGP model if model state used is saturation
            - For discharge, no mapping needed
            - For scale data: no mapping needed


        .. note::
        
            - For ERT data each ERT mesh node is an observation

        THIS SHOULD BE MOVED TO DA CLASS

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
        **kwargts : TYPE
            DESCRIPTION.

        Returns
        -------
        Hx : np.array
            Ensemble of the simulated (predicted) observations.
        '''
        
        
        if parallel:
            path_fwd_CATHY_list = [] # list of ensemble path of cathy output
            for ens_nb in self.ens_valid: # loop over ensemble files
                path_fwd_CATHY_list.append(os.path.join(self.workdir,
                                        self.project_name, 
                                        'DA_Ensemble/cathy_' + str(ens_nb+1)))



            
        Hx_ens = [] # matrice of predicted observation for each ensemble realisation
        for ens_nb in self.ens_valid: # loop over ensemble files
            # print('ens nb:' + str(ens_nb))

            path_fwd_CATHY = os.path.join(self.workdir, self.project_name,
                                  './DA_Ensemble/cathy_' + str(ens_nb+1))
                
            df_psi = self.read_outputs(filename='psi',
                                       path=os.path.join(path_fwd_CATHY,
                                                         self.output_dirname)
                                       )
            df_sw = self.read_outputs(filename='sw',
                                      path=os.path.join(path_fwd_CATHY,
                                                        self.output_dirname)
                                      )
                                      
            # infer soil parameters properties
            # ---------------------------------
            porosity = self.soil_SPP['SPP'][:, 4][0]
    
            # find data to map with dictionnary of observations
            # --------------------------------------------
            obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)
            state = [df_psi[-1],df_sw[-1]]
            
            
            Hx_stacked= [] # stacked predicted observation
            # Loop over observations to map
            # ---------------------------------------------------------------------
            for i, obs_key in enumerate(obskey2map):
                
                # print('obs_key nb:' + obs_key)
                
                if 'tensio' in obs_key:
                    # case 1: pressure head assimilation (Hx_PH)
                    # -------------------------------------------------------------
                    Hx_PH = state[0][obs2map[i]['mesh_nodes']]
                    Hx_stacked.append(Hx_PH)
    
                if 'swc' in obs_key:
                    # case 2: sw assimilation (Hx_SW)
                    # --------------------------------------------------------------------
                    Hx_SW = state[1][obs2map[i]['mesh_nodes']] * porosity
                    Hx_stacked.append(Hx_SW)
                    # note: the value of the porosity can be unique or not depending on the soil physical properties defined
    
                if 'scale' in obs_key:
                    #Atmpot-vf (9) : Potential atmospheric forcing (rain +ve / evap -ve) as a volumetric flux [L^3/T]
                    #Atmpot-v (10) : Potential atmospheric forcing volume [L^3] (See parm input file for units)
                    #Atmpot-r (11) : Potential atmospheric forcing rate [L/T]
                    #Atmpot-d (12) : Potential atmospheric forcing depth [L]
                    #Atmact-vf(13) : Actual infiltration (+ve) or exfiltration (-ve) at atmospheric BC nodes as a volumetric flux [L^3/T]
                    #Atmact-v (14) : Actual infiltration (+ve) or exfiltration (-ve) volume [L^3]
                    #Atmact-r (15) : Actual infiltration (+ve) or exfiltration (-ve) rate [L/T]
                    #Atmact-d (16) : Actual infiltration (+ve) or exfiltration (-ve) depth [L]
                    
                    print('Not yet implemented')
                    
                    df_dtcoupling = self.read_outputs(filename='dtcoupling',
                                              path=os.path.join(path_fwd_CATHY,
                                                                self.output_dirname)
                                              )
                    
                    Hx_scale = state[1][obs2map[i]['mesh_nodes']] * porosity
                    Hx_stacked.append(Hx_scale)
                                       


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
                    if parallel == False:
                        Hx_ERT = self._map_ERT(state,
                                                path_fwd_CATHY,
                                                ens_nb,
                                                **kwargs)
                        if 'pygimli' in obs2map[i]['data_format']:
                            Hx_stacked.append(Hx_ERT['rhoa'])
                        else:
                            Hx_stacked.append(Hx_ERT['resist'])
                        Hx_stacked = np.hstack(Hx_stacked)

            if len(Hx_stacked)>2:
                np.hstack(Hx_stacked)
            Hx_ens.append(Hx_stacked)
        
        if len(Hx_ens)>2:
            Hx_ens= np.hstack(Hx_ens)

        
        # special case of ERT // during sequential assimilation
        # ---------------------------------------------------------------------
        for i, obs_key in enumerate(obskey2map):
            if 'ERT' in obs_key:
                if parallel:
                    Hx_ens_ERT = self._map_ERT_parallel(path_fwd_CATHY_list,
                                                        savefig = True,
                                                        DA_cnb = self.count_DA_cycle,
                                                        )    
                    if len(Hx_ens)>0:
                        Hx_ens = np.vstack([Hx_ens,Hx_ens_ERT])
                        # np.shape(Hx_ens)
                    else:
                        Hx_ens = Hx_ens_ERT

        return Hx_ens #meas_size * ens_size 
    
    
    def _add_2_ensemble_Archie(self, df_Archie_2add):
        '''
        Store in a dataframe Archie relationship for all ensembles and all assimilation times

        Parameters
        ----------
        df_Archie_2add : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        
        if hasattr(self,'Archie') == False:
            self.Archie =  pd.DataFrame(columns=['time',
                                                 'ens_nb', 
                                                 'sw',
                                                 'ER_converted',
                                                 'OL'])
      
        # print(self.Archie['time'].max())
        # print(self.Archie['ens_nb'].max())
        self.Archie = pd.concat([self.Archie,df_Archie_2add])
        # print(self.Archie['time'].max())
        # print(self.Archie['ens_nb'].max())



        
    def _add_2_ensemble_Hx(self, Hx, Hx_2add):
        '''
        Store in an array predicted value for all ensembles and all assimilation times
        '''
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
        # N_col = len(ens_valid)
        
        # Ensemble matrix X of M rows and N  + fill with psi values
        # --------------------------------------------------------------------
        psi = np.zeros([M_rows, N_col])*1e99
        sw = np.zeros([M_rows, N_col])*1e99
        for j in range(self.NENS):
        # for j in range(len(ens_valid)):
            
            try:
                df_psi = out_CT.read_psi(os.path.join(self.workdir, self.project_name,
                                                   "DA_Ensemble/cathy_" + str(j+1),
                                                   'output/psi'))
                psi[:, j] = df_psi[-1, :]
            except:
                pass
            
            try:
                df_sw = out_CT.read_sw(os.path.join(self.workdir, self.project_name,
                                                   "DA_Ensemble/cathy_" + str(j+1),
                                                   'output/sw'))
                sw[:, j] = df_sw[-1, :]
            except:
                pass
        # check if there is still zeros
        if np.count_nonzero(psi==1e99) != 0:
            print('!!!unconsistent filled X, missing values!!!')
            # sys.exit()
            
        return psi, sw, N_col, M_rows 



    def _create_subfolders_ensemble(self, NENS):
        '''
        THIS SHOULD BE MOVED TO DA CLASS
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
            raise ValueError(":worried_face: [b]processor exe not found[/b]")
            # self.console.print(":worried_face: [b]processor exe not found[/b]")

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

        
        
    def _DA_df(self, state=[None,None], state_analysis=None, rejected_ens=[], **kwargs):
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
                     'psi_bef_update', 
                     'sw_bef_update_',
                     'analysis',
                     'OL',# Open loop boolean
                     'rejected'] 


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

            data_df_root = [t_ass*np.ones(len(state[0])), 
                            ens_nb*np.ones(len(state[0])), 
                            state[0],
                            state[1],
                            state[0],
                            True*np.ones(len(state[0])),
                            False*np.ones(len(state[0]))]
            df_DA_ti = pd.DataFrame(np.transpose(data_df_root),
                                  columns=cols_root)
                
        else:
 
            for n in range(self.NENS):               
                
                data_df_root = [self.count_atmbc_cycle*np.ones(len(state[0])), 
                                n*np.ones(len(state[0])), 
                                state[0][:,n],
                                state[1][:,n],
                                state_analysis[:,n],
                                False*np.ones(len(state[0])),
                                int(rejected_ens[n])**np.ones(len(state[0]))
                                ]
                df_DA_ti_ni = pd.DataFrame(np.transpose(data_df_root),
                                      columns=cols_root)
        
                df_DA_ti= pd.concat([df_DA_ti, df_DA_ti_ni], axis=0)   

                       
        self.df_DA= pd.concat([self.df_DA, df_DA_ti], axis=0, ignore_index=True)


        pass




    def _update_input_ensemble(self, NENS, ENS_times,
                               parm_pert=[],
                               update_parm_list=['St. var.'],
                               analysis=[],
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

        self.console.print(":sponge: [b]update input ensemble[/b]")

        

        # resample atmbc file for the given DA window
        # ---------------------------------------------------------------------
        self.selec_atmbc_window(NENS, ENS_times)

        # ---------------------------------------------------------------------
        
        # analysis = []
        # if 'analysis' in kwargs:
        #     analysis = kwargs['analysis']
            
        self.update_ENS_files(parm_pert, 
                              update_parm_list=update_parm_list,
                              cycle_nb=self.count_atmbc_cycle,
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
        '''

        if len(self.grid3d) == 0:
            self.run_processor(IPRT1=3, verbose=False, DAFLAG=0)
            self.grid3d = in_CT.read_grid3d(
                os.path.join(self.workdir, self.project_name))

        # read full simulation atmbc and filter time window
        # ----------------------------------------------------------------------
        df_atmbc, HSPATM, IETO = in_CT.read_atmbc(os.path.join(
                                                                self.workdir, 
                                                                self.project_name,
                                                                'input', 'atmbc'
                                                               ), 
                                                  grid=self.grid3d
                                                  )
        
        if self.count_atmbc_cycle is not None:
            try:
                time_window_atmbc = [
                                        df_atmbc.iloc[self.count_atmbc_cycle]['time'], 
                                        df_atmbc.iloc[self.count_atmbc_cycle+1]['time']
                                    ]
            except:
                pass
        else:
                time_window_atmbc = [ENS_times[0], ENS_times[1]]

        df_atmbc_window = df_atmbc[(df_atmbc['time'] >= time_window_atmbc[0]) &
                      (df_atmbc['time'] <= time_window_atmbc[1])]


        for ens_nb in range(NENS):
            os.chdir(os.path.join(self.workdir, self.project_name,
                                  './DA_Ensemble/cathy_' + str(ens_nb+1)))            
            diff_time = time_window_atmbc[1] - time_window_atmbc[0]
            self.update_parm(TIMPRTi=[0,diff_time],TMAX=diff_time,
                              filename=os.path.join(os.getcwd(), 'input/parm'),
                              backup=True)
            self.console.print(":warning: [b]Making the assumption that atmbc are homogeneous![/b]")          
            VALUE = []
            for t in df_atmbc_window['time'].unique():
                # VALUE.append(df_atmbc_window[df_atmbc_window['time']==t]['value'].mean())
                VALUE.append(df_atmbc_window[df_atmbc_window['time']==t]['value'])
            if len(VALUE)>0:
                self.update_atmbc(HSPATM=1,IETO=0,
                                  time=[0,diff_time], 
                                  VALUE=[VALUE[0]],
                                  filename=os.path.join(os.getcwd(), 'input/atmbc'))
            else:
                pass
        
        pass


    def update_ENS_files(self, parm_pert, update_parm_list, **kwargs):
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        Update by overwriting ensemble files (usually after analysis step or initially to build the ensemble)
        - update state
        - update model parameters: 
            - initial conditions
            - Archie
            - Hydraulic conductivity ks
            - Feddes parameters
            - Atmbc boundary conditions
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
                        df_psi = self.read_outputs(filename='psi',
                                                  path=os.path.join(os.getcwd(),
                                                                    self.output_dirname))
                        self.update_ic(INDP=1, IPOND=0,
                                       pressure_head_ini=analysis[:,ens_nb],
                                       filename=os.path.join(os.getcwd(), 'input/ic'),
                                       backup=True)
            else:
                if kwargs['cycle_nb'] > 0:
                    raise ValueError('no state variable update - use last iteration as initial conditions ?')
                    df_psi = self.read_outputs(filename='psi',
                                              path=os.path.join(os.getcwd(),
                                                                self.output_dirname))
                    if kwargs['cycle_nb'] > 0:
                            self.update_ic(INDP=1, IPOND=0,
                                           pressure_head_ini=df_psi[-1, :],
                                           filename=os.path.join(os.getcwd(), 'input/ic'),
                                           backup=False)
                   

        # loop over dict of perturbated variable
        # ----------------------------------------------------------------------
        FeddesParam_mat_ens = [] # matrice of ensemble for Feddes parameters
        Archie_parms_mat_ens = [] # matrice of ensemble for Archie parameters
        VG_parms_mat_ens = [] # matrice of ensemble for VG parameters
        PERMX_het_ens = [] # matrice of ensemble for hydraulic conductivity parameters
        Archie_p_names = ['porosity', 'rFluid_Archie','a_Archie','m_Archie','n_Archie','pert_sigma_Archie']
        VG_p_possible_names = ['n_VG', 'thetar_VG', 'alpha_VG','VGPSATCELL_VG']
        VG_p_possible_names_positions_in_soil_table = [5,6,7,7]
        
        for parm_i, key in enumerate(update_parm_list):  # loop over perturbated variables dict
            key_root = re.split('(\d+)', key)  
            if len(key_root)==1:
                key_root.append('0')

            # ------------------------------------------------------------------
            for ens_nb in range(self.NENS):
                # change directory according to ensmble file nb
                os.chdir(os.path.join(self.workdir, self.project_name,
                                      './DA_Ensemble/cathy_' + str(ens_nb+1)))

                # atmbc update
                # --------------------------------------------------------------
                if key == 'atmbc':
                    print('hietograph perturbated not yet implemented')
                    self.update_atmbc(verbose=True)

                # ic update (i.e. water table position update)
                # --------------------------------------------------------------
                elif key.casefold() in 'ic'.casefold():      
                    if kwargs['cycle_nb']==0:
                        self.update_ic(INDP=0, IPOND=0,
                                       pressure_head_ini=parm_pert[key
                                           ]['ini_perturbation'][ens_nb],
                                       filename=os.path.join(os.getcwd(), 'input/ic'),
                                       backup=True)

                # kss update
                # --------------------------------------------------------------
                elif key_root[0].casefold() in 'ks'.casefold():
                    if len(PERMX_het_ens)==0:   
                        SPP =  self.soil_SPP['SPP_map'] 
                        PERMX_het_ens = np.ones([len(SPP['PERMX']),self.NENS])
                        PERMX_het_ens[:] = np.nan
                    
                    PERMX_het_ens[int(key_root[1]),ens_nb]=parm_pert[key][update_key][ens_nb]
                    
                    SPP_ensi = []
                    for es in range(self.NENS):
                        spd = dict()
                        spd['PERMX'] = list(PERMX_het_ens[:,ens_nb])
                        spd['PERMY'] = list(PERMX_het_ens[:,ens_nb])
                        spd['PERMZ'] = list(PERMX_het_ens[:,ens_nb])
                        SPP_ensi.append(spd)
                        
                    SPP.update(SPP_ensi[ens_nb])
                    
                    # print(PERMX_het_ens[:,ens_nb])
                    if np.isnan(PERMX_het_ens[:,ens_nb]).any()==False:
                        self.update_soil(SPP=SPP,
                                          FP=self.soil_FP['FP_map'],
                                          verbose=True,
                                          filename=os.path.join(os.getcwd(), 'input/soil')
                                          )  # specify filename path
                # VG parameters update
                # --------------------------------------------------------------
                elif key_root[0] in VG_p_possible_names: 
                    
                    # equivalence_CATHY = {
                    #                      'thetar_VG':'VGRMCCELL',
                    #                      'alpha_VG':'-1/VGPSATCELL',
                    #                      'n_VG':'VGNCELL',
                    #                      }
                    
                    # equivalence_CATHY
                    # retention curves parameters VGN, VGRMC, and VGPSAT
                    # - 'VGNCELL' (NSTR, NZONE): van Genuchten curve exponent  = n
                    # - 'VGRMCCELL' (NSTR, NZONE): residual moisture content = \thetaR
                    # - 'VGPSATCELL' (NSTR, NZONE): van Genuchten curve exponent --> 
                    #                               VGPSAT == -1/alpha (with alpha expressed in [L-1]);

                    SPP =  self.soil_SPP['SPP_map'] 
                    if len(VG_parms_mat_ens)==0:                       
                        VG_parms_mat = self.soil_SPP['SPP'] 
                        np.shape(VG_parms_mat)
                        if len(SPP['PERMX'])<2:
                            VG_parms_mat_ens = np.matlib.repmat(VG_parms_mat[0], self.NENS, 1)
                        else:
                            # VG_parms_mat_ens = np.repeat(VG_parms_mat, self.NENS, axis=0)
                            print('Non homogeneous soil VG properties not implemented')

                    
                    if len(SPP['PERMX'])<2:
                        # for k in parm_pert.keys():
                        for k in VG_p_possible_names:
                            match = [parm_incr for parm_incr in parm_pert.keys() if k in parm_incr]
                            if len(match)>0:
                                # print(match)
                                # print(k)
                                idVG = VG_p_possible_names_positions_in_soil_table[VG_p_possible_names.index(k)]
                                # print(idVG)
                                VG_parms_mat_ens[:,idVG] = parm_pert[match[0]][update_key]
                            else:
                                pass
                    else:
                        # VG_parms_mat_ens[:,:,idVG] = np.matlib.repmat(parm_pert[key][update_key],len(VG_parms_mat),1).T
                        print('Non homogeneous soil VG properties not implemented')
                   
                    VG_parms_ensi = []
                    for es in range(self.NENS):
                        fed = dict()
                        for i, f in  enumerate(list(SPP.keys())):
                            if len(SPP['PERMX'])<2:
                                if f == 'alpha_VG': # convert to CATHY VGPSATCELL == -1/alpha
                                    fed[f] = -1/VG_parms_mat_ens[es,i]
                                else:
                                    fed[f] = [VG_parms_mat_ens[es,i]]
                            else:
                                print('Non homogeneous soil VG properties not implemented')
                                # fed[f] = list(VG_parms_mat_ens[es,:,i])
                        VG_parms_ensi.append(fed)
                    
                    self.update_soil(SPP=VG_parms_ensi[ens_nb],
                                     FP=self.soil_FP['FP_map'],
                                     verbose=False,
                                     filename=os.path.join(os.getcwd(), 'input/soil'))  # specify filename path

                # FeddesParam update
                # --------------------------------------------------------------
                elif key_root[0] in ['PCANA', 'PCREF', 'PCWLT', 'ZROOT', 'PZ', 'OMGC']:                   
                    SPP =  self.soil_SPP['SPP_map'] 
                                        
                    if len(FeddesParam_mat_ens)==0:                       
                        FeddesParam_mat = self.soil_FP['FP']
                        FeddesParam_mat_ens = np.repeat(FeddesParam_mat, self.NENS, axis=0)
                        FeddesParam_mat_ens = np.reshape(FeddesParam_mat_ens,(self.NENS,self.cathyH["MAXVEG"],6))

                    idFeddes = ['PCANA', 'PCREF', 'PCWLT', 'ZROOT', 'PZ', 'OMGC'].index(key_root[0])
                    
                    
                    # if self.cathyH["MAXVEG"]<2: # Case where only 1 vegetation type
                    FeddesParam_mat_ens[:,int(key_root[1]),idFeddes] = parm_pert[key][update_key]
                    # FeddesParam_mat_ens[:,:,idFeddes] = parm_pert[key][update_key]                       
                    # else: # Case where only 1 vegetation type

                    
                    FeddesParam_ensi = []
                    for es in range(self.NENS):
                        fed = dict()
                        for i, f in  enumerate(['PCANA', 'PCREF', 'PCWLT', 'ZROOT', 'PZ', 'OMGC']):
                            fed[f] = list(FeddesParam_mat_ens[es,:,i])
                        FeddesParam_ensi.append(fed)
                    self.update_soil(FP=FeddesParam_ensi[ens_nb],
                                     SPP=SPP,
                                     verbose=False,
                                     filename=os.path.join(os.getcwd(), 'input/soil'))  # specify filename path

                # Archie_p update
                # --------------------------------------------------------------
                elif key_root[0] in Archie_p_names:
                    # self.Archie_parms = {'rFluid':rFluid, 'a':a, 'm':m, 'n':n, 'pert_sigma':pert_sigma}
                    idArchie = Archie_p_names.index(key_root[0])
                    if len(Archie_parms_mat_ens)==0:
                        for p in Archie_p_names:
                            if len(self.Archie_parms[p]) != self.NENS:
                                self.Archie_parms[p]=list(self.Archie_parms[p]*np.ones(self.NENS))
                        Archie_parms_mat_ens = np.zeros([len(Archie_p_names),self.NENS]) 
                        Archie_parms_mat_ens[:] = np.nan
                    Archie_parms_mat_ens[idArchie,ens_nb] = parm_pert[key][update_key][ens_nb]    
                                        
                    if np.isnan(Archie_parms_mat_ens[idArchie,:]).any()==False:
                        self.Archie_parms[key_root[0]]=list(Archie_parms_mat_ens[idArchie,:])
                        
                    # print(self.Archie_parms)

                else:
                    # variable perturbated to ommit: state, Archie param
                    var_pert_to_ommit = ['St. var.']

                    if not key in var_pert_to_ommit:
                        raise  ValueError('key:' + str(key) + ' to update is not existing!')
            

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
            

        if hasattr(self, 'df_performance') is False:
            # initiate if simulation just started             
            # -------------------------------------------------------------
            self.df_performance = pd.DataFrame()   
            self.RMSE_avg_stacked = []
            self.RMSE_sensor_stacked = []
                    
            
            # NAs default     
            # -------------------------------------------------------------
            # RMSE_avg = np.nan # root mean square error (RMSE)
            # NMRMSE = np.nan # time-averaged normalized root mror: unexpected indent
            
           
        
        # average differences over the number of ensemble
        # ALL OBSERVATIONS
        # ------------------------------
        
        
        all_Obs_diff_mat = np.zeros(np.shape(prediction))
        for i in range(len(data)):
            for j in range(np.shape(prediction)[1]): # Loop over ensemble collumns
                all_Obs_diff_mat[i,j] = abs(data[i]-prediction[i,j])
                # all_Obs_diff_mat[:,j] = abs(data[i]-prediction[i,j])
        
        all_Obs_diff_avg = np.nansum(all_Obs_diff_mat,axis=1)*(1/len(prediction))

        # average differences over all the observations
        all_Obs_RMSE_avg_ti = np.sum(all_Obs_diff_avg,axis=0)*(1/len(data))
        
            
        
        # # compute metrics for each observation variable
        # # ------------------------------------------
        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------
        obs2eval_key = [] # data type 2 eval
        obs2eval = []    # data dict 2 eval

        _ , obs2eval_key = self._get_data2assimilate(list_assimilated_obs)
        
        
        start_line_obs = 0
        for s, name_sensor in enumerate(obs2eval_key): # case where there is other observations than ERT
            obs2eval, _ = self._get_data2assimilate([name_sensor], match=True)

            # number of observation at a given time
            # for ERT number of observation = is the number of grid cells
            n_obs = len(obs2eval)
            prediction2eval = prediction[start_line_obs:start_line_obs+n_obs]
            obs2eval_diff_mat = np.zeros(np.shape(prediction2eval))
                        
            for i in range(len(obs2eval)):
                for j in range(np.shape(prediction2eval)[1]): # Loop over ensemble collumns
                    obs2eval_diff_mat[i,j] = abs(obs2eval[i]-prediction2eval[i,j])
    
            obs2eval_diff_avg = np.nansum(obs2eval_diff_mat,axis=1)*(1/len(prediction2eval))

            start_line_obs = start_line_obs + n_obs
            
            if 'ERT' in name_sensor:
                RMSE_sensor_ti = np.sum(obs2eval_diff_avg,axis=0)*(1/len(data))
                # Plot here scatter Transfer resistance scatter plot between the observed and simulated data
                # plt.scatter(obs2eval[i],prediction2eval[i,j])
            else:
                RMSE_sensor_ti = obs2eval_diff_avg
    
    
            # compute normalised RMSE
            #---------------------------------------------------------
            
            self.RMSE_avg_stacked.append(all_Obs_RMSE_avg_ti) 
            self.RMSE_sensor_stacked.append(RMSE_sensor_ti) 
        
            # if hasattr(self, 'RMSE') is False:
            if t_obs == 0:
                NMRMSE_sensor_ti = RMSE_sensor_ti
                NMRMSE_avg_ti = all_Obs_RMSE_avg_ti
            else:
                NMRMSE_sensor_ti = (1/(t_obs+1))*np.sum(self.RMSE_sensor_stacked)
                NMRMSE_avg_ti = (1/(t_obs+1))*np.sum(self.RMSE_avg_stacked)
    

    
    
            # root names for the collumns name                
            # -------------------------------------------------------------
            cols_root = ['time', 
                         'ObsType',
                         'RMSE'+ name_sensor,
                         'RMSE_avg',
                         'NMRMSE'+ name_sensor,
                         'NMRMSE_avg',
                         'OL']        
            
            # root data                
            # -------------------------------------------------------------
            data_df_root = [[t_obs], 
                            [name_sensor],
                            [RMSE_sensor_ti],
                            [all_Obs_RMSE_avg_ti],
                            [NMRMSE_sensor_ti],
                            [NMRMSE_avg_ti],
                            [OL_bool]]
            
            
            df_performance_ti = pd.DataFrame(np.array(data_df_root,dtype=object).T,
                                  columns=cols_root)
            
    
            # concatenate with main RMSE dataframe
            # -------------------------------------------------------------
            self.df_performance= pd.concat([self.df_performance, df_performance_ti], 
                                       axis=0, 
                                       ignore_index=True)
            

            

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
        datetime =  None
        if 'datetime' in kwargs:
            datetime = kwargs['datetime']


        # discharge type read
        # ---------------------------------------------------------------------
        if data_type == 'discharge':
            df = in_meas.read_discharge(filename)
            units = '$m^{3}/s$'


        # weight type read
        # ---------------------------------------------------------------------
        if data_type == 'scale':
            # infer ET from weight data
            if isinstance(filename, str):
                df = in_meas.read_scale(filename)
            else:
                df = filename
                filename = None

            units = '$mm/s$'
            obs_cov_type = None 



        # discharge type read
        # ---------------------------------------------------------------------
        elif data_type == 'swc':
            if isinstance(filename, str):
                df = in_meas.read_swc(filename)
            else:
                df = filename
                filename = None
            units = '$m^{3}/m^{3}$'
            obs_cov_type = None 
            # point sensors --> covariance between the sensors is formed later

            # if 'date' in df.keys:
            #     date_yyyymmdd_hhmmss = df['date']
            #     dict_obs_2add.update(
            #          date_yyyymmdd_hhmmss = date_yyyymmdd_hhmmss
            #          )
                
                
        # tensiometer type read
        # ---------------------------------------------------------------------
        elif data_type == 'tensiometer':
            if isinstance(filename, str):
                df = in_meas.read_tensiometers(filename)
            else:
                df = filename
                filename = None

            units = '$kPa$'
            # convert in m
            # -------------------
            
            df = utils_CT.kPa2m(df)
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

            df, dict_ERT = in_meas.read_ERT(filename,data_format)
            units = '$\Omega$'
            
            if 'elecs' in dict_ERT.keys():
                elecs = dict_ERT['elecs']

                
            dict_obs_2add.update(
                     elecs = elecs
                     )


        # no file specified (raise error)
        # ---------------------------------------------------------------------
        else:
            print('no file specified')

            
        dict_obs_2add.update(
                            filename = filename,
                            data_type =  data_type,
                            units= units,  # units
                            data=  df,
                            data_err=  data_err,
                            mesh_nodes =  mesh_nodes,
                            assimilation_times=  tA,
                            datetime =  datetime
                            )
        
        
        # add optionnal data metadata to the dictionnary
        # ---------------------------------------------------------------------
        if 'meta' in kwargs:
            meta = kwargs['meta']
            dict_obs_2add.update(meta)
            
        
            
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
            
            # resipy format
            try:
                err_arr = np.tile(dict_obs['data_err'],len(dict_obs['data']))
                data_cov = np.diag(err_arr)
                
            # pygimli format
            except:
                err_arr = np.tile(dict_obs['data_err'],len(dict_obs['data']['a']))
                data_cov = np.diag(err_arr)

            
            


        return data_cov 


    def resynchronise_times(self,data_measure,atmbc_df):
        ''' old key is elapsed time in second from the first observation, 
        while new key is from the first atmbc time
        '''
        for d in range(len(data_measure.keys())):
            items_dict = list(self.dict_obs.items())
            # print(items_dict)           
            # list(items_dict[d][1].keys())
            for sensor in list(items_dict[d][1].keys()):
                new_key = (atmbc_df[atmbc_df['sensor_name']==sensor]['diff'].dt.total_seconds().to_numpy()[d])
                old_key = list(data_measure.keys())[d]
                self.dict_obs[old_key][sensor]['assimilation_times'] = new_key
                self.dict_obs[new_key] = self.dict_obs.pop(old_key)
        pass
    
    # self.dict_obs.keys()

    def dictObs_2pd(self):
        '''dict of observation to dataframe of observation'''
        df_obs = pd.DataFrame.from_dict(self.dict_obs).stack().to_frame()
        df_obs = pd.DataFrame(df_obs[0].values.T.tolist(), 
                                index=df_obs.index)
        df_obs.index.names= ['sensorNameidx','assimilation time']
        return df_obs

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
            if type(dict_props[dp]) == float: # homogeneous properties
                self.update_mesh_vtk(prop=dp, prop_value=np.ones(nnodes)*dict_props[dp])
            elif type(dict_props[dp]) == list: # homogeneous properties
                if len(dict_props[dp]) == 1:
                    self.update_mesh_vtk(prop=dp, prop_value=np.ones(nnodes)*dict_props[dp])


        
        pass


    def map_prop_veg(self,dict_props):
               
        if hasattr(self,'veg_map')==False:
            warnings.warn('no known existing vegetation map.')
            pass
        
        else:
            map_veg_dict = {}
            update_map_veg = {}
          
            for d in dict_props.keys():
                map_veg = np.zeros(np.shape(self.veg_map))
                for i, value in enumerate(dict_props[d]):
                    map_veg[self.veg_map == i+1]=value
                update_map_veg[d] = map_veg
            
            return update_map_veg
        
        
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
        if filename == 'mbeconv':
            df=out_CT.read_mbeconv(path)
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

    def find_nearest_node(self, node_coords, grid3d=[]):
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

        if len(grid3d)==0:
            grid3d=in_CT.read_grid3d(self.project_name)
 

        closest_idx=[]
        closest=[]
        
        for i, nc in enumerate(node_coords):
            # euclidean distance
            d=((grid3d['nodes_idxyz'][:, 1] - nc[0]) ** 2 +
                   (grid3d['nodes_idxyz'][:, 2] - nc[1]) ** 2 +
                   (abs(grid3d['nodes_idxyz'][:, 3]) - abs(nc[2])) ** 2
                   ) ** 0.5

            closest_idx.append(np.argmin(d))
            closest.append(grid3d['nodes_idxyz'][closest_idx[i], 1:])

            threshold=5e-1
            if d[np.argmin(d)] > threshold:
                self.console.print(":warning: [b]No node close to the required points![/b]")                
                print(d[np.argmin(d)])
                print(grid3d['nodes_idxyz'][closest_idx[i], 1:])
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





    def backup_simu(self):
        '''
        Save a copy of the simulation for reuse within python
        
        Returns
        -------
        project_filename.pkl

        '''
        
        with open(os.path.join(self.workdir, self.project_name, self.project_name + '.pkl'), 'wb') as f:
            pickle.dump(CATHY,f)
        
        f.close()
        
   
    def backup_results_DA(self,meta_DA=[]):
        '''
        Save minimal dataframes of the simulation for result visualisation within python
        
        Returns
        -------
        project_filename_df.pkl

        '''
        
        
        with open(os.path.join(self.workdir, self.project_name, self.project_name + '_df.pkl'), 'wb') as f:
            pickle.dump(meta_DA,f)
            pickle.dump(self.dict_parm_pert,f)
            pickle.dump(self.df_DA,f)
            pickle.dump(self.dict_obs,f)
            
            
            # pickle.dump(
                
            #     'data': Data: Sensors: 64 data: 1400, nonzero entries: ['a', 'b', 'k', 'm', 'n', 'r', 'rhoa', 'valid']
                
                
                
            #     )
            
            if hasattr(self,'df_performance'):
                pickle.dump(self.df_performance,f)

            
            if hasattr(self,'Archie'):
                pickle.dump(self.Archie,f)
                

        
        f.close()
        

    def load_pickle_backup(self,filename=''):
        
        if len(filename)==0:
            filename = os.path.join(self.workdir, self.project_name, self.project_name + '_df.pkl')
    
        # print(filename)
        backup_list = []
        all_names = [
                     'meta_DA',
                     'dict_parm_pert',
                     'df_DA',
                     'dict_obs',
                     'df_performance',
                     'Archie',
                     ]
        names = []
        i=0
        with open(filename,"rb") as f:
            while True:
                try:
                    backup_list.append(pickle.load(f))
                    names.append(all_names[i])
                    i += 1
                except EOFError:
                    break
        f.close()
        
        dict_backup = {}
        for i, n in enumerate(names):
            dict_backup[n]=backup_list[i]
            
        # return df_DA, dict_parm_pert, df_performance, dict_obs
        return dict_backup
