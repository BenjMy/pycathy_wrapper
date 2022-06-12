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
        recompile : TYPE, optional
            DESCRIPTION. The default is True.
        runProcess : TYPE, optional
            DESCRIPTION. The default is True.
        verbose : TYPE, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.
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
        """
        This method receives arguments for the Data ASSIMILATION methods. 
        MUST BE IMPLEMENTED BY CHILD CLASSES.
        """
        raise NotImplementedError





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
            self.hapin[hapin[i]] = tmp_param_value[i]

        # print('--'*30)
        # print(self.hapin)

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
            for t, v in zip(time, v_atmbc):

                if verbose == True:
                    print(t, v)
                atmbcfile.write("{:.3e}".format(t) + "\t" + "time" + "\n")
                # atmbcfile.close()

                if isinstance(v, float) | isinstance(v, int):
                    atmbcfile.write("{:.3e}".format(v) + "\t" + "VALUE" + "\n")
                    # atmbcfile.write(str(v) + "\t" + "VALUE" + "\n")
                else:
                    print(t, v)
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
        self._write_SOIL_file(SoilPhysProp, FeddesParam)
        
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

        if backup:
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
    
        print(filename)
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
