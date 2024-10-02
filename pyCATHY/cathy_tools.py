"""
Main class controlling the wrapper.

This class is managing three main components: 
    
1. Reading/writing inputs and output files 

We used a generic formalation for the function's names based on the CATHY files. 

Example:
-------
    In order to update the soil file:
        `update_soil(arguments)`
    In order to update the atmbc file:
        `update_atmbc(arguments)`
        
    Note that the name of the arguments are similar to the name of the variable to update
    In order to update the parm file with a new minimum time step DTMIN:
        `update_parm(DTMIN=1e3)`
        
The update function is composed of three main actions:
    - Set the defaut parameters (it is actually not reading the input file to set the defaut parameters)
    - replace values by the new ones introduced via the function argument
    - write the new file (it overwrites the existing file)
    
Note also that updating an input file may affect the CATHY.H control file. 
When needed, the code takes care of updating values of the CATHY.H file retroactivelly.

Remenber to update all your prepo files before calling `run_preprocessor()`, and
all you input files before calling `run_processor()` 

2. Compiling the fortran files
3. Running executable

Once all the files were updated, the `run_preprocessor()` or `run_processor()` 
are taking care of recompiling the source files via bash cmd.
   
"""

import os
import sys
import warnings
from pathlib import Path

import numpy as np
from git import Repo  # In order to fetch the file directly from the repo

import pyCATHY.meshtools as mt
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.plotters import cathy_plots as plt_CT

warnings.filterwarnings("ignore", category=DeprecationWarning)

import glob
import pickle
import re
import shutil
import subprocess
import time  # measure time of computation
from collections import OrderedDict
from os import listdir
from os.path import isfile, join

import pandas as pd
import pyvista as pv  # read .vtk files
import rich.console
from rich import print

import matplotlib.pyplot as plt
from scipy.interpolate import Rbf

# -----------------------------------------------------------------------------
def subprocess_run_multi(pathexe_list):
    """
    Run multiple exe files in parallel.

    Multiprocessor functions outside main CATHY object
    # https://stackoverflow.com/questions/44144584/typeerror-cant-pickle-thread-lock-objects
    """

    print(f"x= {pathexe_list}, PID = {os.getpid()}")
    os.chdir(pathexe_list)
    callexe = "./" + "cathy"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p = subprocess.run(
            [callexe],
            text=True,
            # capture_output=True
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    return p


def make_console(verbose):
    """
    Start up the :class:`rich.console.Console` instance we'll use.

    Parameters
    ----------
    verbose : bool
        Whether or not to print status messages to stderr.
    """
    return rich.console.Console(stderr=True, quiet=not verbose)


# -----------------------------------------------------------------------------
class CATHY:
    """
    Main CATHY object.

    When instantiated it creates the tree project directories with 'prj_name' as root folder.
    The src files are fetched from the online repository if not existing
    (note that it is possible to call a specific version).

    """

    def __init__(
        self,
        dirName=None,
        prj_name="my_cathy_prj",
        notebook=False,
        version="1.0.0",
        verbose=True,
        **kwargs,
    ):
        """
        Create CATHY object.

        ..note:
            All variables in CAPITAL LETTERS use the semantic from the CATHY legacy fortran codes,
            while the others are created exclusively for the wrapper.

        ..note:
            All file variable are self dictionnary objects; example: soil file is contained in self.soil

        """
        self.t0 = time.time()  # executation time estimate
        self.console = make_console(verbose)
        self.console.print(":checkered_flag: [b]Initiate CATHY object[/b]")
        self.verbose = False

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
        # os.chdir(self.workdir)

        self.project_name = prj_name

        if not os.path.exists(os.path.join(self.workdir, self.project_name)):
            os.makedirs(os.path.join(self.workdir, self.project_name), exist_ok=True)

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

        # THIS SHOULD BE REMOVED #
        # dict related to Data assimilation
        # ---------------------------------------------------------------------
        self.DAFLAG = False  # Flag to trigger data assimilation process
        # (self.dict_obs is assimiliatiom times size)
        self.dict_obs = OrderedDict()  # dictionnary containing all the observation data
        # (self.stacked_data_cov is data_size * assimiliatiom times size)
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
                            os.path.join(self.workdir, self.project_name, "src")
                        )
                        shutil.rmtree(
                            os.path.join(self.workdir, self.project_name, "input")
                        )
                        shutil.rmtree(
                            os.path.join(self.workdir, self.project_name, "prepro")
                        )

        # fetch src files if not existing from Gitbucket repo
        # ---------------------------------------------------------------------
        if not os.path.exists(os.path.join(self.workdir, self.project_name, "src")):
            self.console.print(":worried_face: [b]src files not found[/b]")
            self.console.print("working directory is:" + str((self.workdir)))

            if version == "1.0.0":
                try:
                    # Repo.clone_from(
                    #    "https://bitbucket.org/cathy1_0/cathy.git",
                    #    os.path.join(self.workdir, self.project_name, "tmp_src"),
                    #    branch="master",
                    #)
                    
                    Repo.clone_from(
                        "https://github.com/CATHY-Org/CATHY_src_org.git",
                        os.path.join(self.workdir, self.project_name, "tmp_src"),
                        branch="main",
                    )
                    
                    
                    self.console.print(":inbox_tray: [b]Fetch cathy src files[/b]")
                    shutil.move(
                        os.path.join(self.workdir, self.project_name, "tmp_src/src"),
                        os.path.join(self.workdir, self.project_name, "src"),
                    )
                    
                    pathsrc = os.path.join(
                        self.workdir, self.project_name, "tmp_src/runs/weilletal/"
                    )

                    onlyfiles = [
                        f for f in listdir(pathsrc) if isfile(join(pathsrc, f))
                    ]

                    for (
                        file
                    ) in (
                        onlyfiles
                    ):  
                        shutil.move(
                            os.path.join(pathsrc, file),
                            os.path.join(self.workdir, self.project_name, file),
                        )
                        
                except:
                    print("no internet connection to fetch the files")
                    sys.exit()
                    pass
                
            if version == "G. Manoli":
                print("fetch cathy G. Manoli src files")
                path_manoli = "/home/ben/Documents/CATHY/CathyGitbucket/Test_Gabriele/1_Gabriele_Piante_NON_modificato/CATHY_RWU_ABL_1D/"
                shutil.copytree(
                    path_manoli, os.path.join(self.workdir, self.project_name, "src")
                )

        if not os.path.exists(os.path.join(self.workdir, self.project_name, "prepro")):
            self.console.print(
                ":inbox_tray: [b]Fetch cathy prepro src files[/b]"
            )
            shutil.move(
                os.path.join(
                    self.workdir,
                    self.project_name,
                    "tmp_src/runs/weilletal/prepro",
                ),
                os.path.join(self.workdir, self.project_name, "prepro"),
            )
            # os.remove(os.path.join(
            #      self.workdir, self.project_name, "prepro") + '/dem' )
        if not os.path.exists(os.path.join(self.workdir, self.project_name, "input")):
            self.console.print(":inbox_tray: [b]Fetch cathy input files[/b]")
            shutil.move(
                os.path.join(
                    self.workdir,
                    self.project_name,
                    "tmp_src/runs/weilletal/input",
                ),
                os.path.join(self.workdir, self.project_name, "input"),
            )

        if not os.path.exists(os.path.join(self.workdir, self.project_name, "output")):
                self.create_output(output_dirname="output")

                shutil.rmtree(
                    os.path.join(self.workdir, self.project_name, "tmp_src")
                )
                os.remove(os.path.join(self.workdir, self.project_name, "readme.txt"))

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
                            os.path.join(self.workdir, self.project_name, "output")
                        )
                        shutil.rmtree(
                            os.path.join(self.workdir, self.project_name, "vtk")
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

        self.console.print(":cooking: [b]gfortran compilation[/b]")
        for loopi in range(2):  # run it twice (to avoid the first error)
            os.chdir(os.path.join(self.workdir, self.project_name, "prepro/src/"))

            # clean all files compiled
            for file in glob.glob("*.o"):
                os.remove(file)

            # if self.notebook==False:
            bashCommand = "gfortran -O -o pycppp mpar.f90 mbbio.f90 wbb_sr.f90 csort.f90 qsort.f90 depit.f90 cca.f90 smean.f90 dsf.f90 facet.f90 hg.f90 mrbb_sr.f90 bb2shp_sr.f90 shape.f90 dbase.f90 streamer.f90 cppp.f90"
            try:
                p = os.system(bashCommand + "> /dev/null 2>&1")
                # run it twice (to avoid the first error)
                p = os.system(bashCommand + "> /dev/null 2>&1")

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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            p = subprocess.run(
                [bashcmd],
                text=True,
                input=my_data,
                # capture_output=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
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
            idoutlet = np.where(self.DEM == min(np.unique(self.DEM)))
            self.DEM[idoutlet[0], idoutlet[1]] = max(np.unique(self.DEM))

        if hasattr(self, 'DEM') is False:
            
            DEM_mat, DEM_header = in_CT.read_dem(
                                        os.path.join(self.workdir, self.project_name, "prepro/dem"),
                                        os.path.join(self.workdir, self.project_name, "prepro/dtm_13.val"),
                                    )
            self.DEM = DEM_mat
        
        # if hasattr(self, 'hapin') is False:

            # self.update_cathyH(
            #                     ROWMAX=self.hapin['N'],
            #                     COLMAX=self.hapin['M']
            #                     )
            
        # if np.shape(self.DEM)[1]==self.hapin['N']:
        #     self.console.rule(''':warning: !transposing DEM dimension!
        #                           This should be avoided in a future version
        #                           :warning:''',
        #                             style="yellow")
        #     self.update_prepo_inputs(N=self.hapin['M'],
        #                               M=self.hapin['N'])


        # self.update_zone()
        os.chdir(self.workdir)

        pass

    def recompileSrc(self, verbose=False):
        """
        Other option is self.run_processor(runProcess=False)

        """
        ti = time.time()  # executation time estimate
        self.console.print(
            ":hammer_and_wrench: [b] Recompile src files[/b] ["
            + str(int(abs(ti - self.t0)))
            + "s]"
        )
        
        
        # gfortran_executable = 'gfortran.exe'

        # # Execute the where command to find the path to gfortran
        # where_command = ['where', gfortran_executable]
        # result = subprocess.run(where_command, stdout=subprocess.PIPE, shell=True, check=True)
        
        # # Get the path to gfortran from the output of the where command
        # gfortran_path = result.stdout.decode('utf-8').strip().split('\r\n')[0]

        # clean all files previously compiled
        for file in glob.glob("*.o"):
            os.remove(file)
        # list all the fortran files to compile and compile
        for file in glob.glob("*.f"):
            
            # filepath = Path(self.workdir) / self.project_name / 'src' / str(file)
            filepath = str(file)
            bashCommand = "gfortran -c " + str(filepath)
            # gfortran -c *.f
            # gfortran *.o -L\MinGW\lib -llapack -lblas -o cathy
            # print(bashCommand)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                process = subprocess.Popen(
                    bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                # if verbose==True:
                #     output, error = process.communicate()
                output, error = process.communicate()
        # list all the fortran compiled files to compile and run
        files = ""
        for file in glob.glob("*.o"):
            # files += " " + os.path.join(self.workdir, self.project_name, 'prepro/src/' + str(file))
            files += " " + str(file)
        bashCommand = "gfortran" + files + " -llapack -lblas -o " + self.processor_name
        ti = time.time()  # executation time estimate
        self.console.print(
            ":cooking: [b]gfortran compilation[/b] ["
            + str(int(abs(ti - self.t0)))
            + "s]"
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            process = subprocess.Popen(
                bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            output, error = process.communicate()
            
        print(error)


        try:
            shutil.move(
                os.path.join(
                    self.workdir, self.project_name, "src", self.processor_name
                ),
                os.path.join(self.workdir, self.project_name, self.processor_name),
            )
        except:
            self.console.print(":pensive_face: [b]Cannot find the new processsor[/b]")

        pass

    def run_processor(self, recompile=True, runProcess=True, verbose=False, **kwargs):
        """
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
        """

        # set parallel flag for mutirun (used for DA ensemble)
        # --------------------------------------------------------------------
        parallel = False
        if "parallel" in kwargs:
            parallel = kwargs["parallel"]

        # update parm and cathyH
        # --------------------------------------------------------------------
        # VERY VERY IMPORTANT NEVER COMMENT !

        self.check_DEM_versus_inputs() # to uncomment
        # if len(kwargs)>0:
        self.update_parm(**kwargs) # to uncomment
        self.update_cathyH(**kwargs) ### to uncomment
        # self.cathyH
        if recompile:
            # recompile
            # --------------------------------------------------------------------
            os.chdir(os.path.join(self.workdir, self.project_name, "src"))
            self.recompileSrc()

        if "path_CATHY_folder" in kwargs:
            os.chdir(kwargs["path_CATHY_folder"])
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
        if runProcess:
            
            # t0 = time.time()  # executation time estimate
            self.console.print(":athletic_shoe: [b]Run processor[/b]")
            callexe = "./" + self.processor_name

            # case of simple simulation
            # ----------------------------------------------------------------
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                
                p = subprocess.run(
                    [callexe],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                
                # process = subprocess.Popen(
                #     [callexe], stdout=subprocess.PIPE, stderr=subprocess.PIPE
                # )
                # # if verbose==True:
                # #     output, error = process.communicate()
                # output, error = process.communicate()
                
                # Redirecting output to /dev/null or NUL
                # output_redirection = "> /dev/null 2>&1" if os.name != 'nt' else "> NUL 2>&1"
                
                # # Execute the command without displaying the shell
                # os.system(f"{callexe} {output_redirection}")
                # result = os.popen(f"{callexe} {output_redirection}").read()


                
                # p = subprocess.run(
                #     [callexe],
                #     stdout=subprocess.DEVNULL, 
                #     stderr=subprocess.DEVNULL,
                # )
                
                # p = subprocess.Popen([callexe], 
                #                         )
                
                # p = subprocess.Popen([callexe], 
                #                         stdout=subprocess.PIPE, 
                #                         stderr=subprocess.PIPE
                #                         )

                # warnings.simplefilter("ignore")
                # p = subprocess.Popen([callexe], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                # p = subprocess.run(
                #     [callexe],
                #     text=True,
                #     capture_output=True,
                #     check=True
                #     )
                
            # print('verbose='+ str(verbose))
            if verbose:
                # stdout, stderr = p.communicate()  # This will block until the process finishes
                # print(result)
                print(p.stdout)
                # print(p.stderr)
            os.chdir(os.path.join(self.workdir))
            
                        
            try:
                self.grid3d = out_CT.read_grid3d(
                    os.path.join(self.workdir, self.project_name, 'output', 'grid3d')
                )
            except:
                self.run_processor(IPRT1=3)
                
            # if not "nnod3" in self.grid3d.keys():
            #     self.run_processor(IPRT1=3)
                
            
            # computation time
            # ----------------------------------------------------------------
            # t1 = time.time()
            # self.total_computation_time = t1 - t0

        return

    def create_output(self, output_dirname="output"):
        """
        Create output directories

        Parameters
        ----------
        output_dirname : str, optional
            Name of the dir. The default is "output".

        Returns
        -------
        None.

        """

        # create project_name/output folders
        if not os.path.exists(os.path.join(self.project_name, output_dirname)):
            os.mkdir(os.path.join(self.workdir, self.project_name, output_dirname))

        if not os.path.exists(os.path.join(self.project_name, "vtk")):
            os.mkdir(os.path.join(self.workdir, self.project_name, "vtk"))

        # compute the processor (make clean make) using the Makefile
        pass

    def display_time_run(self):

        return self.total

    #%% Checkings

    def check_DEM_versus_inputs(self):
        """
        Check shape consistency between attributes and DEM

        """

        try:
            if hasattr(self, "veg_map"):
                if np.shape(self.veg_map) != np.shape(self.DEM):
                    raise ValueError(
                        "Inconsistent shapes between vegetation map and DEM - need to update veg first"
                    )
        except:
            self.update_veg_map(indice_veg=np.ones(np.shape(self.DEM)))
            pass

    #%%
    # --------------------------------------------------------------------------- #
    # update files
    # --------------------------------------------------------------------------- #

    def update_cathyH(self, verbose=False, **kwargs):
        """
        Update CathyH file input based on **kwargs parameters.

        Parameters ---------- **kwargs : all variable parameters included into cathyH. For instance
        ROWMAX # maximum NROW, with NROW = number of rows in the DEM.

        Returns ------- type Write a modified cathyH file based on new parameters parsed

        """

        indent = ""
        for kk, value in kwargs.items():
            if kk == "indent":
                indent = value

        if verbose:
            self.console.print(
                indent + ":arrows_counterclockwise: [b]Update cathyH files[/b]"
            )

        CATHYH_file = open(
            os.path.join(self.workdir, self.project_name, "src", "CATHY.H"), "r"
        )
        Lines0_109 = CATHYH_file.readlines()
        CATHYH_file.close()

        # check if hapin and parm exist
        # create them if not existing
        if hasattr(self, "hapin") is False:
            self.update_prepo_inputs()
        if hasattr(self, "dem_parameters") is False:
            self.update_dem_parameters()
        if hasattr(self, "veg_map") is False:
            self.update_veg_map()
            
        if "NR" not in self.parm:
            self.console.rule(
                ":warning: warning messages above :warning:", style="yellow"
            )
            self.console.print(
                r"""
                            The parm dictionnary is empty
                            Falling back to defaults to update CATHYH
                            This can have consequences !!
                            """,
                style="yellow",
            )
            self.console.rule("", style="yellow")
            self.update_parm()

        DEMRES = 1                
        ROWMAX = self.hapin["M"]
        COLMAX = self.hapin["N"]
        
        
        if len(self.cathyH) == 0:

            self.cathyH = {
                # "ROWMAX": self.hapin["N"],  # maximum NROW, with NROW = number of rows in the DEM
                # "COLMAX": self.hapin["M"],  # maximum NCOL, with NCOL = number of columns in the DEM
                "ROWMAX": ROWMAX,  # maximum NROW, with NROW = number of rows in the DEM
                "COLMAX": COLMAX,  # maximum NCOL, with NCOL = number of columns in the DEM
                # 'COARSE': ,
                "MAXCEL": int(self.hapin["N"]) * int(self.hapin["M"]),
                "MAXRES": 1,
                "DEMRES": DEMRES,
                "NODMAX": (int(self.hapin["N"]) / DEMRES + 1) * (int(self.hapin["M"]) / DEMRES + 1),
                "NTRMAX": 2* (int(self.hapin["N"]) * int(self.hapin["M"]))/ (DEMRES * DEMRES),
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
                "MAXZON": self.dem_parameters["nzone"],  # maximum NZONE, with NZONE = number of material types in the porous medium
                "MAXTRM": self.hapin["N"]*self.hapin["M"]*self.dem_parameters["nstr"]*30,
                "MAXIT": 30,
                "NRMAX": self.parm["NR"]+1,
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
                "MAXVEG": self.MAXVEG, #
                # "MAXVEG": 1, #
            }

        for kk, value in kwargs.items():
            if kk in self.cathyH.keys():
                self.cathyH[kk] = value
        
        self.cathyH['MAXCEL'] = int(self.cathyH["ROWMAX"]) * int(self.cathyH["COLMAX"])
        self.cathyH['NODMAX'] = int((int(self.cathyH["ROWMAX"]) / DEMRES + 1) * (int(self.cathyH["COLMAX"]) / DEMRES + 1))
        self.cathyH['NTRMAX'] =  int((2*self.cathyH['MAXCEL'])/ (DEMRES * DEMRES))
        self.cathyH['MAXTRM'] = self.cathyH["ROWMAX"]*self.cathyH["COLMAX"]*self.dem_parameters["nstr"]*30      
        self.cathyH['NFACEMAX'] = self.cathyH['NODMAX']*3
        # self.cathyH['MAXTRM'] = 1599000
        # MAXTRM=1599000
        
        
        
        
        # cathyH_laC = {
        #     # "ROWMAX": 247,  # maximum NROW, with NROW = number of rows in the DEM
        #     # "COLMAX": 221,  # maximum NCOL, with NCOL = number of columns in the DEM
        #     # # 'COARSE': ,
        #     # "MAXRES": 1,
        #     # "DEMRES": 1,
        #     "NODMAX": 28440,
        #     "NTRMAX": 2*27848,
        #     "MAXTRM": 2364075,
        #     "NRMAX": 1,
        # }
        # for k in cathyH_laC:
        #     self.cathyH[k] = cathyH_laC[k]
            
        
        # ---------------------------------------------------------------------
        # write new cathy H
        with open(
            os.path.join(self.workdir, 
                         self.project_name, 
                         "src",
                         "CATHY.H"
                         ), "w+") as CATHYH_file:
            for i, l in enumerate(Lines0_109):
                if i < 109:
                    CATHYH_file.write(l)

            CATHYH_file.write(
                "      PARAMETER (ROWMAX={},COLMAX={},DEMRES={})\n".format(
                    self.cathyH["ROWMAX"], self.cathyH["COLMAX"], self.cathyH["DEMRES"]
                )
            )
            CATHYH_file.write("      PARAMETER (MAXCEL=ROWMAX*COLMAX,MAXRES=1)\n")
            # CATHYH_file.write(
            #     "      PARAMETER (NODMAX=(ROWMAX/DEMRES+1)*(COLMAX/DEMRES+1))\n"
            # )
            CATHYH_file.write(
                "      PARAMETER (NODMAX={})\n".format(self.cathyH["NODMAX"])
                )
            # CATHYH_file.write(
            #     "      PARAMETER (NTRMAX=2*MAXCEL/(DEMRES*DEMRES))\n".format()
            # )
            CATHYH_file.write(
                "      PARAMETER (NTRMAX={})\n".format(self.cathyH["NTRMAX"])
                )
            CATHYH_file.write(
                "      PARAMETER (NP2MAX=1,MAXSTR={})\n".format(self.cathyH["MAXSTR"])
            )
            CATHYH_file.write(
                # "      PARAMETER (NFACEMAX={})\n".format(self.cathyH['NFACEMAX']))
                "      PARAMETER (NFACEMAX={})\n".format(self.cathyH['NFACEMAX']))
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
            CATHYH_file.write("      PARAMETER (NIAUXMAX=NMAX + MAXTRM + 1)\n".format())
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

    def update_hapin(self, Lines, hapin, tmp_lnb, tmp_param_value):

        L = []
        tmp_param_value_new = []

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
                                line = xnew[0] + "=              " + str(xnew[1]) + "\n"
                        tmp_param_value_new.append(value)

                count += 1  # count line nb
            L.append(line)

        # # write the new hap.in file

        hap_file = open(
            os.path.join(self.workdir, self.project_name, "prepro/hap.in"), "w+"
        )
        hap_file.writelines(L)
        hap_file.close()

    def update_prepo_inputs(self, DEM=None, verbose=False, show=False, **kwargs):
        """
        Update default prepro inputs i.e. hap.in and dtm_13.val files based on kwargs

        Parameters
        ----------
        DEM : TYPE, optional
            DESCRIPTION. The default is None.
        verbose : TYPE, optional
            DESCRIPTION. The default is False.
        show : TYPE, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        """

        # print("update hap.in")
        self.console.print(":arrows_counterclockwise: [b]Update hap.in file[/b]")

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
                list_tmp = tmp_param_value[i].split("/t")
                tmp_param_value[i] = [float(x) for x in list_tmp]

            self.hapin[hapin[i]] = tmp_param_value[i]

        # if self.hapin["dr"] != self.hapin["delta_x"]:
        # if  self.hapin["dr"] % self.hapin["delta_x"] != 0:
        #     print("adapt rivulet param to dem resolution")
        #     self.hapin["dr"] = self.hapin["delta_x"]
            # sys.exit()

        for key, value in kwargs.items():
            self.hapin[key] = value

        self.update_hapin(Lines, hapin, tmp_lnb, tmp_param_value)


        self.update_dem_parameters(**kwargs)
        # dtm_13.val
        # If we start with a DEM file ("dtm_13.val") for an already delineated
        # catchment (i.e., a "catchment" DEM file instead of a "full" DEM file), then
        # only the last step in the above procedure (i.e., run "cppp" just once) is
        # needed (make sure that "Boundary channel construction" is set to 0 in
        # "hap.in").

        if DEM is not None:
            self.DEM = DEM
            self.hapin["N"] = np.shape(DEM)[1]
            self.hapin["M"] = np.shape(DEM)[0]

            self.update_hapin(Lines,
                              hapin, 
                              tmp_lnb, 
                              tmp_param_value
                              )

            # check presence of the outlet
            if len(np.unique(DEM)) == 1:
                print("Error: outlet not defined")
                DEM[0, 0] = 0
                
            self.console.print(":arrows_counterclockwise: [b]Update dtm_13 file[/b]")
            self.update_dem(DEM)
            self.update_zone()

            with open(
                os.path.join(self.workdir, self.project_name, "prepro/dtm_13.val"), "w+"
            ) as f:
                # use exponential notation
                np.savetxt(f, DEM, fmt="%1.4e")
                # np.shape(DEM)

            self.update_cathyH(
                ROWMAX=self.hapin["M"],
                COLMAX=self.hapin["N"],
                MAXSTR=self.dem_parameters['nstr'],
                DEMRES=min([self.hapin["delta_x"], self.hapin["delta_y"]])
            )

        self.update_dem_parameters(**kwargs)

        if show == True:
            plt_CT.show_dem(DEM, self.hapin)

        pass

    def update_dem_parameters(self, **kwargs):
        """
        Update DEM parameters - if no args reset to default values


        ivert
        =0 each layer will be parallel to the surface, including the base of the 3‐d grid.
            ZRATIO is applied to each vertical cross section.
        =1 base of the 3‐d grid will be flat.
            ZRATIO is applied to each vertical cross section.
        =2 base of the 3‐d grid will be flat, as will the NSTR‐1 horizontal cross sections above it.
            ZRATIO is applied only to the vertical cross section having the lowest elevation.
        =3 for each cell of the dem a single depth value is read in file input IIN60 (basement).
            ZRATIO is applied to each vertical cross section.

        =4 the first NSTR‐1 layers from the surface will be parallel to the surface and the base of the 3‐d grid will be flat.
            ZRATIO is applied only to the vertical cross section having the lowest elevation.

        isp
        =0 for flat surface layer (only one Z value is read in, and is replicated to all surface nodes);
            otherwise surface layer is not flat (Z values read in for each surface node);
            (for ISP=0, IVERT=0, 1, and 2 yield the same 3‐d mesh, given the same values of BASE and ZRATIO).

        base
        Value which defines the thickness or base of the 3‐d mesh.
        For IVERT=0, BASE is subtracted from each surface elevation value,
            so that each vertical cross section will be of thickness BASE,
            and the base of the 3‐d mesh will be parallel to the surface.
        For IVERT=1 or 2, BASE is subtracted from the lowest surface elevation value,
            say ZMIN, so that each vertical cross section will be of thickness (Z ‐ ZMIN) + BASE,
            where Z is the surface elevation for that cross section.
            The base of the 3‐d mesh will thus be flat.

        """

        self.console.print(
            ":arrows_counterclockwise: [b]Update dem_parameters file [/b]"
        )

        if hasattr(self, "dem_parameters") == False:
            dem_parameters_dict = in_CT.read_dem_parameters(
                os.path.join(
                    self.workdir,
                    self.project_name,
                    "input",
                    "dem_parameters",
                )
            )

            self.dem_parameters = dem_parameters_dict

        if 'zratio(i).i=1.nstr' in self.dem_parameters:
            self.dem_parameters['zratio(i),i=1,nstr'] = self.dem_parameters.pop('zratio(i).i=1.nstr')

        # create dictionnary from kwargs
        for keykwargs, value in kwargs.items():
            if keykwargs == "zratio":
                key = "zratio(i),i=1,nstr"
                value = list(value)
                if (abs(float(sum(value))-float(1)))>1e-3:
                    self.console.rule(
                        ":warning: warning messages above :warning:", style="yellow"
                    )
                    self.console.print(
                        r"The sum of all the layers is not equal to 1 but to {}".format(
                            np.round(sum(value), 2)
                        ),
                        style="yellow",
                    )
                    self.console.rule("", style="yellow")
            else:
                key = keykwargs

            try:
                self.dem_parameters[key]
                self.dem_parameters[key] = value
            except:
                pass

        dem_parameters_tmp = {}
        for key, value in self.dem_parameters.items():
            if isinstance(value, list):
                strlst = "\t".join(str(e) for e in value)
                self.dem_parameters[key] = strlst

        #         if (key +'_list') in self.dem_parameters.keys() is False:
        #             dem_parameters_tmp[key +'_list'] = value


        if 'zratio' in kwargs:
            self.dem_parameters['nstr'] = len(kwargs['zratio'])
        # self.dem_parameters = self.dem_parameters | dem_parameters_tmp
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
                    # print(str(list(self.dem_parameters.values())[counth]))
                    dem_parametersfile.write(
                        str(list(self.dem_parameters.values())[counth]) + "\t" + "\n"
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
                        str(list(self.dem_parameters.keys())[counth]) + "\t" + "\n"
                    )
                    counth += 1

        dem_parametersfile.close()

        pass

    def update_lakes_map(self, lakes_map=[]):
        """
        Update lakes_map file

        Parameters
        ----------
        dem : np.array([]), optional
            2d numpy array describing the zone for x and y direction.
            The array dim must be equal to DEM dim.
            The default is []."""

        with open(
            os.path.join(self.workdir, self.project_name, "prepro/lakes_map"), "w+"
        ) as lakesfile:

            lakesfile.write("north:     0" + "\n")
            lakesfile.write("south:     " + str(self.hapin["yllcorner"]) + "\n")
            lakesfile.write("east:     0" + "\n")
            lakesfile.write("west:     " + str(self.hapin["xllcorner"]) + "\n")

            lakesfile.write("rows:     " + str(self.hapin["N"]) + "\n")
            lakesfile.write("cols:     " + str(self.hapin["M"]) + "\n")
            np.savetxt(lakesfile, lakes_map, fmt="%i")
            # np.shape(DEM)
        lakesfile.close()
        
        pass
                
            
    def update_dem(self, DEM=[]):
        """
        Update zone file

        Parameters
        ----------
        dem : np.array([]), optional
            2d numpy array describing the zone for x and y direction.
            The array dim must be equal to DEM dim.
            The default is []."""

        with open(
            os.path.join(self.workdir, self.project_name, "prepro/dem"), "w+"
        ) as demfile:

            demfile.write("north:     0" + "\n")
            demfile.write("south:     " + str(self.hapin["yllcorner"]) + "\n")
            demfile.write("east:     0" + "\n")
            demfile.write("west:     " + str(self.hapin["xllcorner"]) + "\n")

            demfile.write("rows:     " + str(self.hapin["N"]) + "\n")
            demfile.write("cols:     " + str(self.hapin["M"]) + "\n")
            np.savetxt(demfile, DEM, fmt="%1.4e")
            # np.shape(DEM)
        demfile.close()
        
        with open(
            os.path.join(self.workdir, self.project_name, "prepro/dtm_13.val"), "w+"
        ) as f:
            np.savetxt(f, DEM, fmt="%1.4e")
                
        self.DEM = DEM
        pass

    def update_zone(self, zone=[]):
        """
        Update zone file

        Parameters
        ----------
        zone : np.array([]), optional
            2d numpy array describing the zone for x and y direction.
            The array dim must be equal to DEM dim.
            The smallest numbering is 1 (not 0!)
            If zone is empty, reset to homogeneous zone = 1
            The default is [].

        Notes
        -----
        .. note::
            - Create object zone variable
            - Update dem_parameters file (number of zones).
            - Update parm file.
            - Update CATHYH file.
            - Update vtk mesh file.

        """
        self.console.print(":arrows_counterclockwise: [b]update zone file [/b]")

        with open(
            os.path.join(self.workdir, self.project_name, "prepro/zone"), "w+"
        ) as zonefile:

            zonefile.write("north:     0" + "\n")
            zonefile.write("south:     " + str(self.hapin["yllcorner"]) + "\n")
            zonefile.write("east:     0" + "\n")
            zonefile.write("west:     " + str(self.hapin["xllcorner"]) + "\n")
            zonefile.write("rows:     " + str(self.hapin["M"]) + "\n")
            zonefile.write("cols:     " + str(self.hapin["N"]) + "\n")
                        
            if len(zone) == 0:
                zone = np.c_[np.ones([self.hapin["M"], self.hapin["N"]])]
                np.savetxt(zonefile, zone, fmt="%i")
            else:
                # if np.shape(zone)== :
                np.savetxt(zonefile, zone, fmt="%i")

        zonefile.close()

        self.zone = zone

        # update number of zone in the dem parameter file
        self.update_dem_parameters(nzone=len(np.unique(zone)))
        self.update_parm()
        self.update_cathyH(MAXZON=len(np.unique(zone)))

    # %% INPUT FILES

    def update_parm(self, **kwargs):
        """
        Update parm file
        .. note:
            - create object zone parm
            - Updated CATHYH file (NPRT and NUMVP).
        """

        self.console.print(":arrows_counterclockwise: [b]Update parm file [/b]")
        warnings_parm = []


        if len(self.parm) == 0:
            #%%
            dict_parm = in_CT.read_parm(os.path.join(self.workdir, 
                                                     self.project_name, 
                                                     "input",
                                                     "parm"
                                                     )
                                        )
            #%%

            self.parm = dict_parm
            
            if '(TIMPRT(I)I=1NPRT)' in self.parm:
                self.parm['(TIMPRT(I),I=1,NPRT)'] = self.parm.pop('(TIMPRT(I)I=1NPRT)')


        # create self.parm dictionnary from kwargs
        # --------------------------------------------------------------------
        for kk, value in kwargs.items():
            if kk in self.parm.keys():
                if self.verbose == True:
                    print(f"key: {kk} | value: {value}")

            # times of interest TIMPRTi
            # ----------------------------------------------------------------
            if kk == "TIMPRTi":
                key = "(TIMPRT(I),I=1,NPRT)"

                self.parm[key] = list(value)

            # points of interest NODVP
            # ----------------------------------------------------------------
            elif kk == "NODVP":
                key = "(NODVP(I),I=1,NUMVP)"
                self.parm[key] = value

                # check if consistency between node of interest and
                # number of nodes of interest
                # ------------------------------------------------------------
                if len(value) != self.parm["NUMVP"]:
                    warnings_parm.append(
                        "Adjusting NUMVP with respect to NODVP requested" + "\n"
                    )
                    self.parm["NUMVP"] = len(value)

            # points of interest NR
            # ----------------------------------------------------------------
            elif kk == "ID_NR":
                key = "CONTR(I),I=1,NR"
                self.parm[key] = value

                # check if consistency between node of interest and
                # number of nodes of interest
                # ------------------------------------------------------------
                if len(value) != self.parm["NR"]:
                    warnings_parm.append(
                        "Adjusting NR with respect to CONTR requested" + "\n"
                    )
                    self.parm["NR"] = len(value)
                    
            # points of interest NR
            # ----------------------------------------------------------------
            elif kk == "ID_QOUT":
                key = "(ID_QOUT(I),I=1,NUM_QOUT)"
                self.parm[key] = value

                # check if consistency between node of interest and
                # number of nodes of interest
                # ------------------------------------------------------------
                if len(value) != self.parm["NUM_QOUT"]:
                    warnings_parm.append(
                        "Adjusting NUM_QOUT with respect to ID_QOUT requested" + "\n"
                    )
                    self.parm["NUM_QOUT"] = len(value)


            # other type of kwargs
            # ----------------------------------------------------------------
            else:
                if kk in self.parm.keys():
                    self.parm[kk] = value

                    # check if consistency between times of interest and TMAX
                    # ------------------------------------------------------------
                    if "(TIMPRT(I),I=1,NPRT)" in kk:
                        if self.parm["TMAX"] < max(self.parm["(TIMPRT(I),I=1,NPRT)"]):
                            self.parm["TMAX"] = max(self.parm["(TIMPRT(I),I=1,NPRT)"])

        # check if consistency between times of interest and
        # number of times of interest
        # ------------------------------------------------------------
        if len(self.parm["(TIMPRT(I),I=1,NPRT)"]) != self.parm["NPRT"]:
            warnings_parm.append(
                "Adjusting NPRT with respect to time of interests requested" + "\n"
            )
            self.parm["NPRT"] = len(self.parm["(TIMPRT(I),I=1,NPRT)"])

        # check if consistency between DELTAT, DTMIN, and DTMAX
        # ------------------------------------------------------------

        if min(np.diff(self.parm["(TIMPRT(I),I=1,NPRT)"])) <= 0:
            raise ValueError("Vtk time steps should be monotonically increasing")

        if self.parm["DTMIN"] > min(np.diff(self.parm["(TIMPRT(I),I=1,NPRT)"])
        # ) or self.parm["DELTAT"] > 1e-1 * min(
        #     np.diff(self.parm["(TIMPRT(I),I=1,NPRT)"])
        ):
            warnings_parm.append(
                "Adjusting DTMIN with respect to time of interests requested" + "\n"
            )
            self.parm["DTMIN"] = min(np.diff(self.parm["(TIMPRT(I),I=1,NPRT)"])) / 1e2

        if self.parm["DELTAT"] < self.parm["DTMIN"]:
            warnings_parm.append(
                "adjusting 2*DTMIN == DELTAT: " + str(self.parm["DTMIN"]) + "\n"
            )
            self.parm["DTMIN"] = self.parm["DELTAT"] / 2

        if self.parm["DELTAT"] >= self.parm["DTMAX"]:
            warnings_parm.append("Adjusting DTMAX == 2*DELTAT" + "\n")
            self.parm["DTMAX"] = 1e2 * self.parm["DELTAT"]

        if self.parm["TMAX"] < max(self.parm["(TIMPRT(I),I=1,NPRT)"]):
            warnings_parm.append(
                "Adjusting TMAX with respect to time of interests requested" + "\n"
            )
            self.parm["TMAX"] = max(self.parm["(TIMPRT(I),I=1,NPRT)"])

        if len(warnings_parm) > 0:
            self.console.rule(
                ":warning: warning messages above :warning:", style="yellow"
            )
            for ww in warnings_parm:
                self.console.print(warnings_parm, style="yellow")
            self.console.rule("", style="yellow")

        # transform array args to list
        # --------------------------------------------------------------------
        for kk, value in self.parm.items():
            if "numpy.array" in str(type(value)):
                value = list(value)
                self.parm[kk] = value

        # update CATHYH
        # --------------------------------------------------------------------
        self.update_cathyH(
            MAXPRT=self.parm["NPRT"],
            MAXVP=self.parm["NUMVP"],
            indent="           :arrow_right_hook:",
            **kwargs,
        )

        # write parm file
        # --------------------------------------------------------------------
        if "filename" in kwargs:
            file2write = kwargs["filename"]
        else:
            file2write = os.path.join(
                self.workdir, self.project_name, self.input_dirname, "parm"
            )

        # backup file during DA scheme cycle
        # --------------------------------------------------------------------
        backup = False
        if "backup" in kwargs:
            backup = kwargs["backup"]

        if backup == True:
            if self.count_DA_cycle is not None:
                dst_dir = file2write + str(self.count_DA_cycle)
                shutil.copy(file2write, dst_dir)

        self._write_parm_file(file2write)

        pass

    def _write_parm_file(self, file2write):
        """
        Overwrite existing parm file

        Returns
        -------
        New overwritten file.

        """        
        if int(self.parm['NR'])==0:
            header_fmt_parm = [3, 3, 2, 4, 4, 3, 3, 2, 4, 3, 3, 4, 4, 4, 2, 1, 1, 1]
            try:
                del self.parm["CONTR(I),I=1,NR"]
            except:
                pass
        else:
            header_fmt_parm = [3, 3, 2, 4, 4, 3, 3, 2, 4, 3, 3, 4, 4, 4, 2, 1, 1, 1, 1]
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
                left = str(list(self.parm.values())[counth : counth + h])
                left = left.strip("[]").replace(",", "")
                left = left.strip("[]").replace("[", "")

                # right = keys
                # ------------------------------------------------------------
                right = str(list(self.parm.keys())[counth : counth + h])
                right = right.strip("[]").replace("',", "")
                right = right.replace("'", "")

                # add left + right
                # ------------------------------------------------------------
                line = left + "\t" + right + "\n"
                counth += h
                parmfile.write(str(line))

        parmfile.close()

        pass

    def update_ic(self, INDP=2, IPOND=0, WTPOSITION=0, verbose=False, **kwargs):
        """
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
            
            
        kwargs : 
            - pressure_head_ini (int): uniform value of pressure head to assign (only if INDP=0)

        Returns
        -------
        new ic file written/overwritten.

        """

        if "shellprint_update" in kwargs:
            if kwargs["shellprint_update"] is False:
                pass
        else:
            self.console.print(":arrows_counterclockwise: [b]Update ic[/b]")

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
            "WTPOSITION": WTPOSITION,
        }

        # write ic file
        # --------------------------------------------------------------------

        if "filename" in kwargs:
            filename = kwargs["filename"]
        else:
            filename = os.path.join(
                self.workdir, self.project_name, self.input_dirname, "ic"
            )

        backup = False
        if "backup" in kwargs:
            backup = kwargs["backup"]

        if backup == True:
            if self.count_DA_cycle is not None:
                dst_dir = filename + str(self.count_DA_cycle - 1)
                shutil.copy(filename, dst_dir)

        with open(filename, "w+") as icfile:

            if INDP == 0:
                icfile.write(str(INDP) + "\t" + str(IPOND) + "\t" + "INDP \n")
                if "pressure_head_ini" in kwargs:
                    np.savetxt(icfile, [kwargs["pressure_head_ini"]], fmt="%1.3e")
                    
                    # try:
                    self.map_prop2mesh({"ic": kwargs["pressure_head_ini"]})
                    # except:
                    #     print('Impossible to map pressure head to mesh')

                else:
                    raise ValueError("Missing initial pressure value")
            elif INDP == 1:
                pressure_head_ini = kwargs["pressure_head_ini"]
                icfile.write(
                    str(INDP)
                    + "\t"
                    + str(IPOND)
                    + "\t"
                    + "INDP"
                    + "\t"
                    + "IPOND"
                    + "\n"
                )
                np.savetxt(icfile, pressure_head_ini, fmt="%1.3e")
                self.map_prop2mesh({"ic": kwargs["pressure_head_ini"]})
                
            elif INDP==2:
                icfile.write(
                    str(INDP)
                    + "\t"
                    + str(IPOND)
                    + "\t"
                    + "INDP"
                    + "\t"
                    + "IPOND"
                    + "\n"
                )

            elif INDP in [3, 4]:
                icfile.write(
                    str(INDP)
                    + "\t"
                    + str(IPOND)
                    + "\t"
                    + "INDP"
                    + "\t"
                    + "IPOND"
                    + "\n"
                )
                icfile.write(str(WTPOSITION) + "\t" + "WTPOSITION" + "\n")

        icfile.close()

        pass

    def update_atmbc(
        self,
        HSPATM=1,
        IETO=0,
        time=None,
        VALUE=[None, None],
        netValue=[],
        show=False,
        verbose=False,
        **kwargs,
    ):
        """
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
        IETO : int, optional
            - =0 for linear interpolation of the atmospheric boundary condition inputs between different different ATMTIM
            - otherwise the inputs are assigned as a piecewise constant function (ietograph).
            The default is 0.
        time : list
            ATMBC Times in seconds. The default is None.
        VALUE : list
            List of array. The default is [Precipitation, EvapoTranspiration].
        show : bool, optional
            Plot atmbc. The default is False.
        verbose : bool, optional
            Display. The default is False.
        **kwargs

        Returns
        -------

        ..note:

                - Update parm file (NPRT).
                - Update CATHYH file (MAXPRT).

        """

        if "shellprint_update" in kwargs:
            if kwargs["shellprint_update"] is False:
                pass
        else:
            self.console.print(":arrows_counterclockwise: [b]Update atmbc[/b]")

        # type(netValue)
        if len(netValue) > 0:
            v_atmbc = netValue
        else:
            # atmbc are homoegenous
            # -----------------------------------------------------------------
            if HSPATM == 1:
                if len(VALUE) == 2:  # take the difference between Precipitation and EvapoTranspiration
                    v_atmbc = VALUE[0] - abs(VALUE[1])
                else: # Assume it is already the net 
                    v_atmbc = VALUE
                    print('Assuming it is already the net')

            else:
                print('spatial v_atmbc must be provided as net between Precipitation and EvapoTranspiration')


        # set default parameters
        # --------------------------------------------------------------------
        self.atmbc = {"HSPATM": HSPATM, "IETO": IETO, "time": time, "VALUE": v_atmbc}

        # len(self.atmbc['time'])
        # overwrite existing input atmbc file
        # --------------------------------------------------------------------

        if "filename" in kwargs:
            filename = kwargs["filename"]
        else:
            filename = os.path.join(
                self.workdir, self.project_name, self.input_dirname, "atmbc"
            )

        backup = True
        if "backup" in kwargs:
            backup = kwargs["backup"]

        if backup == True:
            if self.count_DA_cycle is not None:
                dst_dir = filename + str(self.count_DA_cycle - 1)
                shutil.copy(filename, dst_dir)

        with open(filename, "w+") as atmbcfile:
            atmbcfile.write(
                str(HSPATM) + "\t" + str(IETO) + "\t" + "HSPATM" + "\t" + "IETO" + "\n"
            )

            # atmbc are homogeneous
            # -----------------------------------------------------------------
            if HSPATM == 1:
                for t, v in zip(time, v_atmbc):
                    if verbose:
                        print(t, v)
                    atmbcfile.write("{:.0f}".format(t) + "\t" + "time" + "\n")
                    # atmbcfile.close()
                    if isinstance(v, float) | isinstance(v, int):
                        atmbcfile.write("{:.3e}".format(v) + "\t" + "VALUE" + "\n")
                    else:
                        atmbcfile.write("{:.3e}".format(float(v)) + "\t" + "VALUE" + "\n")
                    # atmbcfile.write(str(v) + "\t" + "VALUE" + "\n")
                        
            # atmbc are heteregeneous
            # -----------------------------------------------------------------
            else:
                for t, v in zip(time, v_atmbc):
                    atmbcfile.write("{:.0f}".format(t) + "\t" + "time" + "\n")
                    np.savetxt(atmbcfile, v, fmt="%.3e")

        atmbcfile.close()

        self.update_parm(TIMPRTi=[time[0],time[-1]], NPRT=2, TMAX=max(time))
        # self.update_parm(TIMPRTi=time, NPRT=len(time), TMAX=max(time))

        # don't need to update if sequential DA as the cathy.exe is already created
        # ---------------------------------------------------------------------
        if "omit_cathyH" not in kwargs:
            self.update_cathyH(MAXPRT=len(time))

        if show:
            if HSPATM ==0:
                print('impossible to plot for non homogeneous atmbc')
            #     # sys.exit()
            else:
            # x_units = "sec"
            # for key, value in kwargs.items():
            #     if key == "x_units":
            #         x_units = value
                if len(netValue) > 0:
                    plt_CT.show_atmbc(time, v_atmbc, **kwargs)
                else:
                    plt_CT.show_atmbc(time, VALUE, **kwargs)
        pass

    def update_nansfdirbc(
        self,
        time=[],
        NDIR=0,
        NDIRC=0,
        NQ3=None,
        no_flow=False,
        pressure_head=[],
        mesh_bc_df=None,
    ):
        """
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
        pressure_head : TYPE, optional
            Specify a value of node pressure head to impose as Dirichlet boundary condition.
            The default is [].


        Returns
        -------
        None.

        """

        # check that the mesh exist
        # --------------------------------------------------------------------
        try:
            self.grid3d = out_CT.read_grid3d(os.path.join(self.workdir,
                                                          self.project_name, 
                                                          'output', 'grid3d')
                                             )
        except OSError:
            print("grid3d missing - need to run the processor with IPRT1=3 first")

        # check that mesh_bound_cond_df exist
        # --------------------------------------------------------------------
        if not hasattr(self, "mesh_bound_cond_df"):
            if len(time)>25:
                print('Nb of times is too big to handle bc condition in the df')
                # time = [0]
            # self.init_boundary_conditions(
            #                                 "nansfdirbc",
            #                                 time,
            #                                 NDIR=NDIR,
            #                                 NDIRC=NDIRC,
            #                                 pressure_head=pressure_head,
            #                                 no_flow=no_flow,
            #                             )
            # self.assign_mesh_bc_df("nansfdirbc", time, 
            #                        pressure_head=pressure_head,
            #                        no_flow=no_flow
            #                        )

        else:         
            # if len(time)>25:
            #     print('Nb of times is too big to handle bc condition in the df')
            #     time = 0
            # self.update_mesh_boundary_cond(mesh_bc_df)
            self.assign_mesh_bc_df("nansfdirbc", time, 
                                   no_flow=no_flow
                                   )
        # apply BC
        # --------------------------------------------------------------------
        # if no_flow:  # Dirichlet  == 0
        #     print("shortcut set_BC_laterals mesh dataframe")
        #     # self.set_BC_laterals(time=time, BC_type='Dirichlet', val=0)
        # else:
        #     raise ValueError(
        #         "Non homogeneous Dirichlet Boundary conditions Not yet implemented"
        #     )

        self.update_parm()
        self.update_cathyH()

        with open(
            os.path.join(
                self.workdir, self.project_name, self.input_dirname, "nansfdirbc"
            ),
            "w+",
        ) as nansfdirbcfile:

            # self.mesh_bound_cond_df

            # To simulate the no-flow boundaries conditions for the bottom and
            # vertical sides of the domain --> NDIR and NDIRC equal to zero
            # -------------------------------------------------------------
            if no_flow:  # Dirichlet
                if len(time) == 0:
                    time = self.atmbc["time"]
                for tt in time:
                    # nansfdirbcfile.write(str(tt) + "\t" + "time" + "\n")
                    nansfdirbcfile.write("{:.0f}".format(tt) + "\t" + "time" + "\n")

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
            else:
                if len(time) == 0:
                    time = self.atmbc["time"]
                for tt in time:
                    # nansfdirbcfile.write(str(tt) + "\t" + "time" + "\n")
                    nansfdirbcfile.write("{:.d}".format(tt) + "\t" + "time" + "\n")

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
                raise ValueError(
                    "Non homogeneous Dirichlet Boundary conditions Not yet implemented"
                )

        nansfdirbcfile.close()

        self.update_parm()
        self.update_cathyH()


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

    def update_nansfneubc(
        self, time=[], 
        ZERO=0, NQ=0, CONTQ=[],
        fixed_flux=None, 
        no_flow=False,
        **kwargs,
    ):
        """
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

        """


        whereBC = None
        val = fixed_flux
        for k in kwargs:
            # {'xmin_bound': 1.2e-07}
            val = kwargs.get(k)
            if type(k)==str:
                whereBC = k
                    
        # check that the mesh existNeuma
        # --------------------------------------------------------------------
        try:
            self.grid3d = out_CT.read_grid3d(os.path.join(self.workdir,
                                                          self.project_name, 
                                                          'output', 'grid3d')
                                             )
        except OSError:
            print("grid3d missing - need to run the processor with IPRT1=3 first")

        # check that mesh_bound_cond_df exist
        # --------------------------------------------------------------------
        if hasattr(self, "mesh_bound_cond_df")==False:
            if len(time)>25:
                print('Nb of times is too big to handle bc condition in the df')
            else:
                self.create_mesh_bounds_df("nansfneubc", 
                                           self.grid3d['mesh3d_nodes'], 
                                           time
                                           )
                # self.init_boundary_conditions("nansfneubc", time, NQ=NQ, ZERO=ZERO)
                # self.assign_mesh_bc_df("nansfneubc", 
                #                        time, 
                #                        no_flow=no_flow
                #                        )
        else:
            if no_flow:
                self.assign_mesh_bc_df("nansfneubc", 
                                       time, 
                                       no_flow=no_flow
                                       )
                self.update_cathyH(NQMAX=1)

            elif whereBC is not None:
                NQ = np.sum(self.mesh_bound_cond_df.set_index('time').loc[0,whereBC]==True)
                CONTQ = list(np.where(self.mesh_bound_cond_df[whereBC]==True)[0])
                nodes_id_CONTQ = list(self.mesh_bound_cond_df.loc[CONTQ,'id_node'].unique())
                # len(nodes_id_CONTQ)
                self.update_cathyH(NQMAX=NQ)

            else:
                raise ValueError('Not yet implemented')
                
            for tt in time:
                self.update_mesh_boundary_cond(
                                                tt,
                                                BC_name='nansfneubc',
                                                BC_val=val,
                                                nodesId=CONTQ
                                               )
        # read existing input nansfneubc file
        # --------------------------------------------------------------------
        with open(os.path.join(self.workdir, 
                               self.project_name, 
                               self.input_dirname, 
                               "nansfneubc"
                               ),
                  "w+",
                  ) as nansfneubcfile:
            if no_flow:  # Neumanm
                print("shortcut set_BC_laterals mesh dataframe")   
                if len(time) == 0:
                    time = self.atmbc["time"]
                for tt in time:
                    nansfneubcfile.write("{:.0f}".format(tt) + "\t" + "time" + "\n")
                    nansfneubcfile.write(
                        str(ZERO) + "\t" + str(NQ) + "\t" + "ZERO" + "\t" + "NQ" + "\n"
                    )
            elif whereBC:
                print("Non homogeneous Neumanm Boundary conditions requested")
                for tt in self.mesh_bound_cond_df.time.unique():
                    nansfneubcfile.write("{:.0f}".format(tt) + "\t" + "time" + "\n")
                    nansfneubcfile.write(
                                        str(ZERO) + "\t" + str(NQ) + "\t" + "ZERO" + "\t" + "NQ" + "\n"
                                        )
                    # write node numbers where Neumann BC applied
                    # --------------------------------------------
                    np.savetxt(nansfneubcfile, nodes_id_CONTQ, fmt="%i")
                    # write values of flow rate (m3/s) for each node numbers where Neumann BC applied
                    # --------------------------------------------
                    values_nansfneubc = self.mesh_bound_cond_df.set_index('time').loc[tt,'nansfneubc'].values
                    values_nansfneubc = list(values_nansfneubc[~np.isnan(values_nansfneubc)])
                    np.savetxt(nansfneubcfile, values_nansfneubc, fmt="%.3e")
                    
                    
        nansfneubcfile.close()    


        pass

    def update_sfbc(self, time=[], sfbc=False, no_flow=False):
        """
        Seepage face boundary conditions at time t

        Parameters
        ----------
        time : np.array([]), optional
            Absolute simulation time in sec.
            The default is [].

        Returns
        -------
        None.

        """

        # check that the mesh exist
        # --------------------------------------------------------------------
        try:
            self.grid3d = out_CT.read_grid3d(os.path.join(self.workdir,
                                                          self.project_name,
                                                          'output', 'grid3d')
                                             )
        except OSError:
            print("grid3d missing - need to run the processor with IPRT1=3 first")

        # check that mesh_bound_cond_df exist
        # --------------------------------------------------------------------
        if hasattr(self, "mesh_bound_cond_df") is False:
            if len(time)>25:
                print('Nb of times is too big to handle bc condition in the df')
                # time = [0]
            # self.init_boundary_conditions("sfbc", time)
            # self.assign_mesh_bc_df('sfbc', time, 
            #                        no_flow=no_flow
            #                        )
        else:
            # if len(time)>25:
            #     print('Nb of times is too big to handle bc condition in the df')
            #     time = 0
            self.assign_mesh_bc_df('sfbc', time, 
                                   no_flow=no_flow
                                   )
        #     self.update_mesh_boundary_cond()

        # apply BC
        # --------------------------------------------------------------------
        # if no_flow:
        #     # if len(time) == 0:
        #     #     time = self.atmbc["time"]
        #     # self.set_BC_laterals(time=time, BC_type='sfbc', val=0)
        #     print("shortcut set_BC_laterals mesh dataframe")

        # else:
        #     raise ValueError(
        #         "Non homogeneous Neumanm Boundary conditions Not yet implemented"
        #     )

        with open(
            os.path.join(self.workdir, self.project_name, self.input_dirname, "sfbc"),
            "w+",
        ) as sfbcfile:

            if len(time) == 0:
                time = self.atmbc["time"]
            for tt in time:
                sfbcfile.write("{:.0f}".format(tt) + "\n")
                sfbcfile.write("0" + "\n")

        sfbcfile.close()

        pass

    def update_soil(
        self,
        IVGHU=[],
        FP=[],
        FP_map=[],
        SPP=[],
        SPP_map=[],
        zone3d=[],
        show=False,
        **kwargs,
    ):
        """
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
            - 'PCREF': float, field capacity --> h3 (≈ -4m)
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


        FP : Dict, optional
            Dictionnary containing the SPP properies per zones (root map, indice of vegetation)

        SPP : np.array, optional
            Array of size (nb of layers*nb of zones)*8 (8: Physical properties)

        SPP_map : Dict or Pd.DataFrame, optional
            Dictionnary containing the SPP properies per zones
            
            
            Example for 2 zones:
                {
                    'PERMX': [0.000188, 9.4e-05],
                    'PERMY': [0.000188, 0.000188],
                    'PERMZ': [0.000188, 0.000188],
                    'ELSTOR': [1e-05, 1e-05],
                    'POROS': [0.55, 0.55],
                    'VGNCELL': [1.46, 1.46],
                    'VGRMCCELL': [0.15, 0.15],
                    'VGPSATCELL': [0.03125, 0.03125]
                 }
                
            
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

            Example for 1 zone:

                             PERMX     PERMY     PERMZ  ...  VGNCELL  VGRMCCELL  VGPSATCELL
                str zone                                ...
                0   0     0.000188  0.000188  0.000188  ...     1.46       0.15     0.03125
                1   0     0.000188  0.000188  0.000188  ...     1.46       0.15     0.03125
                2   0     0.000188  0.000188  0.000188  ...     1.46       0.15     0.03125
                3   0     0.000188  0.000188  0.000188  ...     1.46       0.15     0.03125
                4   0     0.000188  0.000188  0.000188  ...     1.46       0.15     0.03125
                    

        ..note::

            At the file level, this is an example of 3 layers and  2 zones; The inner reading cycle is by zone, and the outer one by layers, i.e.:

            PERMX_z1_str1 PERMY_z1_str1  PERMZ_z1_str1  ELSTOR_z1_str1 POROS_z1_str1 VGNCELL_z1_str1 VGRMCCELL__z1_str1 VGPSATCELL__z1_str1
            PERMX_z2_str1 PERMY_z2_str1  PERMZ_z2_str1  ELSTOR_z2_str1 POROS_z2_str1 VGNCELL_z2_str1 VGRMCCELL__z2_str1 VGPSATCELL__z2_str1
            PERMX_z1_str2 PERMY_z1_str2  PERMZ_z1_str2  ELSTOR_z1_str2 POROS_z1_str2 VGNCELL_z1_str2 VGRMCCELL__z1_str2 VGPSATCELL__z1_str2
            PERMX_z2_str2 PERMY_z2_str2  PERMZ_z2_str2  ELSTOR_z2_str2 POROS_z2_str2 VGNCELL_z2_str2 VGRMCCELL__z2_str2 VGPSATCELL__z2_str2
            PERMX_z1_str3 PERMY_z1_str3  PERMZ_z1_str3  ELSTOR_z1_str3 POROS_z1_str3 VGNCELL_z1_str3 VGRMCCELL__z1_str3 VGPSATCELL__z1_str3
            PERMX_z2_str3 PERMY_z2_str3  PERMZ_z2_str3  ELSTOR_z2_str3 POROS_z2_str3 VGNCELL_z2_str3 VGRMCCELL__z2_str3  VGPSATCELL__z2_str3

            Make sure all 8 parameters for each layer/zone are on the same line.
            Also, make sure NZONE = 2 in dem_parameters and MAXZON in CATHY.H is larger than or equal to NZONE.

        Returns
        -------
        - update parm file
        - update CATHY.H file
        - update mesh vtk file

        """
        if "shellprint_update" in kwargs:
            if kwargs["shellprint_update"] is False:
                pass
        else:
            self.console.print(":arrows_counterclockwise: [b]Update soil[/b]")

        
        # set default parameters if SPP and/or FP args are not existing yet
        # --------------------------------------------------------------------
        if len(self.soil) == 0:
            self.set_SOIL_defaults()
            try:
                df = self.read_inputs('soil', 
                                      MAXVEG=self.MAXVEG
                                      )
                if len(FP_map)==0:
                    # FP_map = df[1].reset_index(drop=True).to_dict(orient='list')
                    FP_map = df[1] #.reset_index(drop=True).to_dict(orient='list')
            except:
                pass

        if len(SPP_map) == 0:
            SPP_map = self.set_SOIL_defaults(SPP_map_default=True)

        if len(FP_map) == 0:
            FP_map = self.set_SOIL_defaults(FP_map_default=True)
            
        # elif len(FP_map) == 0:
        #     # FP_map = self.soil_FP["FP_map"]
        #     FP_map = self.set_SOIL_defaults(FP_map_default=True)

        # check size of the heteregeneity
        # -----------------------------------
        if len(zone3d) == 0:
            print("homogeneous soil")           
            # check size of soil properties map versus nb of zones/ nb of layers
            # --------------------------------------------------------------------
            if isinstance(SPP_map["PERMX"], float):
                if self.dem_parameters["nzone"] != 1:
                    raise ValueError("Wrong number of zones")
        
        else:
            self.zone3d=zone3d
            self.dem_parameters['nzone'] = np.size(self.zone3d[0])
            self.update_dem_parameters()
            self.update_cathyH(MAXZON=self.dem_parameters['nzone'])

        # read function arguments kwargs and udpate soil and parm files
        # --------------------------------------------------------------------
        for kk, value in kwargs.items():
            self.soil[kk] = value
            self.parm[kk] = value

        # loop over Feddes parameters
        # --------------------------------------------------------------------
        for fp in FP_map:  # loop over fedded parameterssoil_het_dim
            self.soil[fp] = FP_map[fp]

        # check consistency between parameters
        # --------------------------------------------------------------------
        if isinstance(SPP_map["PERMX"], float):
            for k in SPP_map:
                SPP_map[k] = [SPP_map[k]]
        
        if isinstance(SPP_map, dict):
            SPP_map_dict = SPP_map
            if hasattr(self, 'zone3d'):
                num_rows, num_cols = np.shape(self.zone3d)[0],(np.shape(self.zone3d)[1] * np.shape(self.zone3d)[2])
                df_SPP_map = self.init_soil_SPP_map_df(num_cols, num_rows)
        
                for i, layersi_zones in enumerate(self.zone3d):
                    layersi_zones = np.hstack(layersi_zones)
                    unique_zones = np.unique(layersi_zones)
                    layersi_zones_prop = np.zeros(np.shape(layersi_zones))
                    for c in SPP_map_dict.keys():
                        for zi in unique_zones:
                            layersi_zones_prop[layersi_zones == zi] = SPP_map_dict[c][int(zi-1)]
        
                        df_SPP_map.loc[(slice(None), i + 1), c] = layersi_zones_prop
            else:
                nzones, nstr = self.dem_parameters["nzone"], self.dem_parameters["nstr"]
                df_SPP_map = self.init_soil_SPP_map_df(nzones, nstr)
        
                for key, values in SPP_map.items():
                    if len(values) == 1:
                        values_all_layers = values*nstr
                        self.console.rule(
                            ":warning: Assuming that soil is homogeneous with depth :warning:", style="yellow"
                        )
                    df_SPP_map[key] = values_all_layers
        
            SPP_map = df_SPP_map

        if (SPP_map["VGRMCCELL"] >= SPP_map["POROS"]).any():
            self.console.rule(
                ":warning: residual water content is > porosity :warning:", style="yellow"
            )

        # create prepro inputs if not existing (containing info about the DEM)
        # --------------------------------------------------------------------
        if hasattr(self, "dem_parameters") is False:
            self.update_prepo_inputs()

        # Soil Physical Properties strat by strat
        # --------------------------------------------------------------------
        self.soil_SPP = {}
        self.soil_SPP["SPP_map"] = SPP_map  # mapping with respect to zones
        if len(SPP) > 0:
            self.soil_SPP["SPP"] = SPP  # matrice with respect to zones
        else:
            SoilPhysProp = SPP_map.values
            self.soil_SPP["SPP"] = SoilPhysProp  # matrice with respect to zones/layers
            
        # Vegetation properties (PCANA,PCREF,PCWLT,ZROOT,PZ,OMGC)
        # --------------------------------------------------------------------
        FeddesParam = self._prepare_SOIL_vegetation_tb(FP_map)
        self.soil_FP = {}
        self.soil_FP["FP"] = FeddesParam
        self.soil_FP["FP_map"] = FP_map  # mapping with respect to zones
        
        map_veg = np.ones(np.shape(self.veg_map))
        for i, value in enumerate(FP_map['PCANA']):
            map_veg[self.veg_map == i + 1] = i + 1
        self.update_veg_map(map_veg)


        if show:
            update_map_veg = self.map_prop_veg(FP_map)
            fig, ax = plt_CT.dem_plot_2d_top(update_map_veg, label="all")
            fig.savefig(
                os.path.join(self.workdir, self.project_name, "map_veg.png"), dpi=400
            )

        # write soil file
        # --------------------------------------------------------------------
        # print('write soil')
        self._write_SOIL_file(self.soil_SPP["SPP"], 
                              FeddesParam, 
                              **kwargs
                              )

        # map SPP to the mesh
        # --------------------------------------------------------------------
        
        dem, dem_header = self.read_inputs('dem')

        if len(zone3d) > 0:
            
            saveMeshPath = os.path.join(self.workdir, self.project_name,
                                        'vtk',
                                        self.project_name + '.vtk'
                                        )
            if 'saveMeshPath' in kwargs:
                saveMeshPath = kwargs.pop('saveMeshPath')
                
            mt.add_markers2mesh(
                                zone3d,
                                dem,
                                self.mesh_pv_attributes,
                                self.dem_parameters,
                                self.hapin,
                                self.grid3d,
                                to_nodes=False,
                                show=show,
                                saveMeshPath=saveMeshPath,
                            )

            # mt.add_markers2mesh(
            #                         zone3d,
            #                         self.mesh_pv_attributes,
            #                         self.dem_parameters,
            #                         self.hapin,
            #                         to_nodes=True
            #                     )

            # for spp in SPP_map:
            #     self.map_prop_2mesh_markers(spp, SPP_map[spp], to_nodes=False)

        pass
    
    def _get_soil_SPP_columnsNames(self):
        
        # Define the column names
        columns = ['PERMX', 'PERMY', 'PERMZ','ELSTOR',
                   'POROS',
                   'VGNCELL', 'VGRMCCELL', 'VGPSATCELL'
                   ]
        return columns
    
    def _get_soil_FP_columnsNames(self):
            
        # Define the column names
        columns = ['PCANA', 'PCREF', 'PCWLT',
                   'ZROOT',
                   'PZ','OMGC'
                   ]
        return columns
    
    def init_soil_FP_map_df(self, nveg):

        columns = self._get_soil_FP_columnsNames()
        FP_map = pd.DataFrame(index=np.arange(1,nveg+1), 
                              columns=columns
                              )
        FP_map.index.name = 'Veg nb'
        return FP_map

    def init_soil_SPP_map_df(self, nzones, nstr):
        
        columns = self._get_soil_SPP_columnsNames()
        
        # Generate the lists of integers for zones and strings
        zones = [ni + 1 for ni in range(nzones)]
        strings = [ni + 1 for ni in range(nstr)]
        
        # Create multi-level index
        multi_index = pd.MultiIndex.from_product([zones, strings], 
                                                 names=['zone', 'str']
                                                 )
        # Create an empty DataFrame with multi-level index and specified columns
        SPP_map = pd.DataFrame(index=multi_index, columns=columns)
        
        return SPP_map
            
            
    def set_SOIL_defaults(self, 
                          FP_map_default=False, 
                          SPP_map_default=False
                          ):

        self.soil = {
            "PMIN": -5.0,
            "IPEAT": 0,
            "SCF": 1.0,  # here we assume that all soil is covered by the vegetation
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

        if FP_map_default:

            # FP_map = {
            #     # Feddes parameters default values
            #     "PCANA": [0.0],
            #     "PCREF": [-4.0],
            #     "PCWLT": [-150],
            #     "ZROOT": [1.0],
            #     "PZ": [1.0],
            #     "OMGC": [1.0],
            # }
            PCANA= 0.0
            PCREF=-4.0
            PCWLT=-150
            ZROOT=1
            PZ=1
            OMGC=1
            # self.soil.update(FP)
            nveg = len(np.unique(self.veg_map))
            FP_map = self.init_soil_FP_map_df(nveg)

            for c in FP_map.columns:
                FP_map[c] = eval(c)
                
            return FP_map

        # # set Soil Physical Properties defaults parameters
        # # --------------------------------------------------------------------

        if SPP_map_default:

            PERMX = PERMY = PERMZ = 1.88e-04
            ELSTOR = 1.00e-05
            POROS = 0.55
            VGNCELL = 1.46
            VGRMCCELL = 0.15
            VGPSATCELL = 0.03125

            # Replace these values with your actual number of zones and strings
            nzones = self.dem_parameters["nzone"]
            nstr = self.dem_parameters["nstr"]
                    
            SPP_map = self.init_soil_SPP_map_df(nzones,nstr)
            
            for c in SPP_map.columns:
                SPP_map[c] = eval(c)
                
            return SPP_map

        pass

    def _prepare_SPP_tb(self, SPP, zone3d):
        """
        prepare SOIL Physical Properties table

        Parameters
        ----------
        SPP : TYPE
            DESCRIPTION.

        Returns
        -------
        np.array describing the SoilPhysProp with rows corresponding to the layer.

        """

        # check number of zones
        if self.dem_parameters["nzone"] > 1 or len(zone3d) > 0:
            if len(SPP["PERMX"]) <= 1:
                for i, spp in enumerate(SPP):
                    SPP[spp] = SPP[spp] * np.ones(self.dem_parameters["nzone"])
            else:
                pass

            # check size of the heteregeneity of SPP
            # ----------------------------------------
            # 1d --> uniform
            # 2d --> lateral variations due to zones defined in surface
            # 3d --> lateral + vertical variations due to zones and strates

            SoilPhysProp = []

            if len(zone3d) == 0:
                # loop over strates
                # -----------------------------------------------------------------
                for istr in range(self.dem_parameters["nstr"]):
                    #  loop over zones (defined in the zone file)
                    # --------------------------------------------------------------
                    LayeriZonei = np.zeros([self.dem_parameters["nzone"], 8])
                    for izone in range(self.dem_parameters["nzone"]):
                        for i, spp in enumerate(SPP):
                            LayeriZonei[izone, i] = SPP[spp][izone]
                    SoilPhysProp.append(LayeriZonei)

            else:
                zones3d_def = np.arange(1, self.hapin["M"] * self.hapin["N"] + 1)
                zones3d_def = np.reshape(zones3d_def, [self.hapin["M"], self.hapin["N"]])
                self.update_zone(zones3d_def)
                # loop over strates
                # -----------------------------------------------------------------
                for istr in range(self.dem_parameters["nstr"]):
                    #  loop over zones (defined in the zone file)
                    # --------------------------------------------------------------
                    LayeriZonei = np.zeros([self.dem_parameters["nzone"], 8])
                    for izone in range(self.hapin["N"] * self.hapin["M"]):
                        for i, spp in enumerate(SPP):
                            flag_zone = int(np.ravel(zone3d[istr])[izone])
                            LayeriZonei[izone, i] = SPP[spp][flag_zone - 1]
                    SoilPhysProp.append(LayeriZonei)
            SoilPhysProp = np.vstack(SoilPhysProp)
            # np.shape(SoilPhysProp)

        # case if there is only one zone in the mesh
        # -------------------------------------------------------------
        else:
            if len(SPP["PERMX"]) <= 1:
                izoneSoil = []
                for spp in SPP:
                    izoneSoil.append(SPP[spp])
                izoneSoil = np.hstack(izoneSoil)
                SoilPhysProp = np.tile(izoneSoil, (self.dem_parameters["nstr"], 1))

            # case where it is heterogeneous in depth (z)
            # ------------------------------------------------------
            else:
                izoneSoil_per_layer = []
                for stri in range(self.dem_parameters["nstr"]):
                    izoneSoil = []
                    for spp in SPP:
                        izoneSoil.append(SPP[spp][stri])
                    izoneSoil_per_layer.append(np.hstack(izoneSoil))
                SoilPhysProp = np.vstack(izoneSoil_per_layer)

        return SoilPhysProp

    def _prepare_SOIL_vegetation_tb(self, FP_map):
        """
        _prepare_SOIL_vegetation_tb

        Parameters
        ----------
        FP_map : dict
            dict containing Feddes parameters.
            - 'PCANA': anaerobiosis point
            - 'PCREF': field capacity
            - 'PCWLT': wilting point
            - 'ZROOT': float root depth
            - 'PZ': float, can be negative, defines the root profile
            - 'OMGC': float, compensatory mechanism
            For details, see http://dx.doi.org/10.1002/2015WR017139
        Returns
        -------
        FeddesParam: numpy array
            table or array describing Feddes parameters for a given DEM




        """
        # Vegetation properties (PCANA,PCREF,PCWLT,ZROOT,PZ,OMGC)
        # --------------------------------------------------------------------
        # Check if root_map file exist and is updated
        # -------------------------------------------
        if hasattr(self, "veg_map") is False:
            if len(FP_map[list(FP_map)[0]]) == 1:
                self.update_veg_map()
            else:
                raise ValueError(
                    "Found multiple values of Feddes zones"
                    + "but vegetation map is not defined"
                )

        
        # Check vegetation heterogeneity dimension
        # ----------------------------------------
        if self.cathyH["MAXVEG"] != len(FP_map[list(FP_map)[0]]):
            raise ValueError(
                "Wrong number of vegetations: PCANA size is "
                + str(len(FP_map[list(FP_map)[0]]))
                + " while MAXVEG is "
                + str(self.cathyH["MAXVEG"])
            )

        # check number of vegetation
        # --------------------------------------------------------------------
        if self.cathyH["MAXVEG"] > 1:
            FeddesParam = np.zeros([self.cathyH["MAXVEG"], 6])
            for iveg in range(self.cathyH["MAXVEG"]):  # loop over veg zones within a strate
                izoneVeg_tmp = []
                for sfp in FP_map:
                    izoneVeg_tmp.append(FP_map[sfp].loc[iveg+1])

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
        """
        _write_SOIL_file

        Parameters
        ----------
        SoilPhysProp : Soil physical properties
            Numpy array of Soil physical properties.
        FeddesParam : Feddes Parameters
            Dictionnatry of Feddes Parameters
        """

        # backup file during DA scheme cycle
        # --------------------------------------------------------------------
        backup = True
        if "backup" in kwargs:
            backup = kwargs["backup"]

            
        # number of side header for each row
        header_fmt_soil = [1, 2, 2, 6, 1, 5, 1, 2, 3]

        # open soil file
        # --------------------------------------------------------------------
        if self.DAFLAG:
            soil_filepath = os.path.join(os.getcwd(), self.input_dirname, "soil")
        else:
            soil_filepath = os.path.join(
                self.workdir, self.project_name, self.input_dirname, "soil"
            )

        if "filename" in kwargs:
            soil_filepath = os.path.join(kwargs["filename"]) #, "soil")
            
            
        # print('v'*13)
        # print(backup)
        # print('backup soil')
        # print(self.count_DA_cycle)
        # print(soil_filepath)

        if backup:
            if self.count_DA_cycle is not None:
                dst_dir = soil_filepath + str(self.count_DA_cycle)
                shutil.copy(soil_filepath, dst_dir)
                # print('backup soil ' + dst_dir)

        with open(os.path.join(soil_filepath), "w+") as soilfile:

            counth = 0  # count header index

            # Write line by line according to header format
            # ----------------------------------------------------------------
            for i, h in enumerate(header_fmt_soil):

                # left = values
                # ------------------------------------------------------------
                left = right = []
                left = str(list(self.soil.values())[counth : counth + h])
                left = left.strip("[]").replace(",", "")

                # right = keys
                # ------------------------------------------------------------
                right = str(list(self.soil.keys())[counth : counth + h])
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
                "PERMX PERMY  PERMZ  ELSTOR POROS,VGNCELL,VGRMCCELL,VGPSATCELL" + "\n"
            )

        soilfile.close()

    def update_veg_map(self, indice_veg=None, show=False, **kwargs):
        """
        Contains the raster map describing which type of vegetation every cell belongs to.


        Parameters
        ----------
        indice_veg : raster, optional
            DESCRIPTION. The default is 1.
        show : bool, optional
            Plot the vegetation raster. The default is False.
        Returns
        -------
        indice_veg : raster
            New vegetation distribution raster.

        """
        
        if indice_veg is None:
            indice_veg, str_hd_rootmap = in_CT.read_root_map(os.path.join(
                                                                             self.workdir, 
                                                                             self.project_name, 
                                                                             self.input_dirname, 
                                                                             "root_map"
                                                                             )
                                                            )
        self.veg_map = indice_veg

        if hasattr(self, "hapin") is False:
            self.update_prepo_inputs()

        with open(
            os.path.join(
                self.workdir, self.project_name, self.input_dirname, "root_map"
            ),
            "w+",
        ) as rootmapfile:
            rootmapfile.write("north:     0" + "\n")
            rootmapfile.write("south:     " + str(self.hapin["yllcorner"]) + "\n")
            rootmapfile.write("east:     0" + "\n")
            rootmapfile.write("west:     " + str(self.hapin["xllcorner"]) + "\n")
            rootmapfile.write("rows:     " + str(self.hapin["M"]) + "\n")
            rootmapfile.write("cols:     " + str(self.hapin["N"]) + "\n")

            if isinstance(indice_veg, int):
                indice_veg = (
                                np.c_[np.ones([int(self.hapin["M"]), int(self.hapin["N"])])]*indice_veg
                             )
                np.savetxt(rootmapfile, indice_veg, fmt="%i")
            else:
                np.savetxt(rootmapfile, indice_veg, fmt="%i")

        rootmapfile.close()

        # exclude vegetation label from number of vegetation if is it outside the DEM domain
        # i.e if DEM values are negative
        # ---------------------------------------------------------------------------------
        if len(np.unique(indice_veg))>1:
            exclude_veg = self._check_outside_DEM(indice_veg)
            self.MAXVEG = len(np.unique(indice_veg)) - exclude_veg
            self.MAXVEG = len(np.unique(indice_veg))# - exclude_veg

            if exclude_veg>0:
                print('excluding outside DEM')
                print('MAXVEG='+ str(self.MAXVEG))
            
        else:
            self.MAXVEG = len(np.unique(indice_veg))


        self.update_cathyH(MAXVEG=self.MAXVEG) # to uncomment

        if show:
            ax = plt_CT.show_indice_veg(self.veg_map, **kwargs)
            return indice_veg, ax
        return indice_veg


    def _check_outside_DEM(self,raster2check):
        
        if hasattr(self,'DEM') is False:
            DEM_mat, DEM_header = in_CT.read_dem(
                                        os.path.join(self.workdir, self.project_name, "prepro/dem"),
                                        os.path.join(self.workdir, self.project_name, "prepro/dtm_13.val"),
                                    )
            
            self.DEM = DEM_mat
        # exclude vegetation label from number of vegetation if is it outside the DEM domain
        # i.e if DEM values are negative

        exclude_out_ind = 0
        if np.min(self.DEM)<0:
        # if len(raster2check[self.DEM<0]):
            exclude_out_ind = 1
            
        return exclude_out_ind
    
        
    #%% Add inputs and outputs attributes to the mesh

    def init_boundary_conditions(self, BCtypName, time, **kwargs):
        """
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

        """
        # self.update_nansfdirbc()
        # self.update_nansfneubc()
        # self.update_sfbc()

        try:
            self.console.print(
                ":orange_square: [b] init boundary condition dataframe [/b]"
            )

            # self.create_mesh_vtk()
            if not "nnod3" in self.grid3d.keys():
                self.run_processor(IPRT1=3)
            if hasattr(self, "mesh_bound_cond_df") is False:
                self.create_mesh_bounds_df(
                                            BCtypName, 
                                            self.grid3d["mesh3d_nodes"], 
                                            time, 
                                            **kwargs
                )
        except:
            raise ValueError("cannot init boundary conditions dataframe")
            pass

        pass

    def check_for_inconsistent_BC(self):
        pass

    def set_BC_laterals(self, time, BC_type="", val=0):
        """
        Set all sides expect surface one
        """

        nnodes = len(
            self.mesh_bound_cond_df[self.mesh_bound_cond_df["time (s)"] == 0]["id_node"]
        )
        for tt in time:
            BC_bool_name = []
            BC_bool_val = []

            for id_node in range(nnodes):
                if self.mesh_bound_cond_df["bound"].loc[int(id_node)] == True:
                    BC_bool_name.append(BC_type)
                    BC_bool_val.append(0)
                else:
                    BC_bool_name.append(None)
                    BC_bool_val.append(val)

            self.update_mesh_boundary_cond(
                time=tt, BC_bool_name=BC_bool_name, BC_bool_val=BC_bool_val
            )

        self.check_for_inconsistent_BC()

        # self.update_mesh_vtk(BC_bool_name,BC_bool_val)

        pass

       

    def get_outer_nodes(self, x, y, z):
        x_min, x_max = np.min(x), np.max(x)
        y_min, y_max = np.min(y), np.max(y)
        z_min, z_max = np.min(z), np.max(z)
        outer_mask = np.logical_or.reduce((x == x_min, x == x_max, y == y_min, y == y_max, z == z_min, z == z_max))
        outer_nodes = np.column_stack((x[outer_mask], y[outer_mask], z[outer_mask]))
        return outer_nodes, outer_mask
        
        
    def create_mesh_bounds_df(self, BCtypName, grid3d, times, **kwargs):
        """
        Create a dataframe with flag for different boundary condtions assigned to each nodes

        Parameters
        ----------
        grid3d : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        mesh_bound_np = np.c_[np.arange(0,len(grid3d)),grid3d]
        self.mesh_bound_cond_df = pd.DataFrame(mesh_bound_np)
        self.mesh_bound_cond_df = self.mesh_bound_cond_df.rename(
                                                                    columns={
                                                                        0: "id_node",
                                                                        1: "x",
                                                                        2: "y",
                                                                        3: "z",
                                                                    }
                                                                )
        self.mesh_bound_cond_df["id_node"] = self.mesh_bound_cond_df["id_node"].astype(int)
        self.mesh_bound_cond_df["time"] = 0
      
        top_id = np.arange(0, self.grid3d['nnod'], 1)
        top_mask = np.zeros(len(self.grid3d['mesh3d_nodes']), dtype=bool)
        top_mask[top_id.astype(int)] = True

        bot_id = np.arange(self.grid3d['nnod3']-self.grid3d['nnod'], 
                           self.grid3d['nnod3'],
                           1
                           )
        bot_mask = np.zeros(len(self.grid3d['mesh3d_nodes']), dtype=bool)
        bot_mask[bot_id.astype(int)] = True
        
        x_threshold = self.grid3d['mesh3d_nodes'][:,0].min()
        xmin_side_mask =  self.mesh_bound_cond_df['x'] <= x_threshold
        
        x_threshold = self.grid3d['mesh3d_nodes'][:,0].max()
        xmax_side_mask =  self.mesh_bound_cond_df['x'] >= x_threshold
        
        y_threshold = self.grid3d['mesh3d_nodes'][:,1].min()
        ymin_side_mask =  self.mesh_bound_cond_df['y'] <= y_threshold
        
        y_threshold = self.grid3d['mesh3d_nodes'][:,1].max()
        ymax_side_mask =  self.mesh_bound_cond_df['y'] >= y_threshold
        
        all_sides_mask = np.c_[xmin_side_mask,
                                xmax_side_mask,
                                ymin_side_mask,
                                ymax_side_mask,
                                ]
        all_sides_mask = np.any(all_sides_mask, axis=1)

      
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
       
        # # Scatter plot in 3D
        # sc = ax.scatter(self.grid3d['mesh3d_nodes'][:, 0],
        #                 self.grid3d['mesh3d_nodes'][:, 1],
        #                 self.grid3d['mesh3d_nodes'][:, 2],
        #                 # c=self.grid3d['mesh3d_nodes'][:, 2],  # Coloring by z-axis values
        #                 c='k',  # Coloring by z-axis values
        #                 # cmap='viridis'
        #                 )  # Optional: Change color map
        # # Scatter plot in 3D
        # sc = ax.scatter(self.grid3d['mesh3d_nodes'][all_sides_mask, 0],
        #                 self.grid3d['mesh3d_nodes'][all_sides_mask, 1],
        #                 self.grid3d['mesh3d_nodes'][all_sides_mask, 2],
        #                 # c=self.grid3d['mesh3d_nodes'][left_side_mask, 2],  # Coloring by z-axis values
        #                 c='r',  # Coloring by z-axis values
        #                 # cmap='viridis'
        #                 )  # Optional: Change color map
            


        # Step 4: Fill the dataframe with flag for outer nodes
        # -------------------------------------------------------------------
        # print('Fill the dataframe with flag for outer nodes')
        self.mesh_bound_cond_df["noflow_bound"] = np.ones(len(self.grid3d['mesh3d_nodes']))*True
        self.mesh_bound_cond_df["top_bound"] = top_mask
        self.mesh_bound_cond_df["bot_bound"] = bot_mask
        self.mesh_bound_cond_df["xmin_bound"] = xmin_side_mask
        self.mesh_bound_cond_df["ymin_bound"] = ymin_side_mask
        self.mesh_bound_cond_df["xmax_bound"] = xmax_side_mask
        self.mesh_bound_cond_df["ymax_bound"] = ymax_side_mask
        self.mesh_bound_cond_df["all_sides"] = all_sides_mask
               
        # Step 4: Replicate df for all given times
        # -------------------------------------------------------------------   
        if len(times) > 1:
            mesh_bound_cond_df_withtimes = pd.concat(
                [self.mesh_bound_cond_df.assign(time=ti) for ti in times]
            ).reset_index(drop=True)
            self.mesh_bound_cond_df = mesh_bound_cond_df_withtimes             
                
            

    def assign_mesh_bc_df(self, BCtypName, times=0, **kwargs):

        
        # step 5 add flag for each type of BC
        # ---------------------------------------
        # specified pressure
        if "nansfdirbc" in BCtypName:
            if "no_flow" in kwargs:
                self.mesh_bound_cond_df[BCtypName] = -99
                mask = self.mesh_bound_cond_df["noflow_bound"] == True
                self.mesh_bound_cond_df.loc[mask, BCtypName] = 0
        # specified flux
        if "nansfneubc" in BCtypName:
            if "no_flow" in kwargs:
                self.mesh_bound_cond_df[BCtypName] = -99
                mask = self.mesh_bound_cond_df["noflow_bound"] == True
                self.mesh_bound_cond_df.loc[mask, BCtypName] = 0
        # specified seepage
        if "sfbc" in BCtypName:
            if "no_flow" in kwargs:
                self.mesh_bound_cond_df[BCtypName] = -99
                mask = self.mesh_bound_cond_df["noflow_bound"] == True
                self.mesh_bound_cond_df.loc[mask, BCtypName] = 0
    
            # print(
            #     "SKip time dependence init boundary condition dataframe - consequences (?)"
            # )

        pass

    def update_mesh_boundary_cond(
        self,
        time,
        BC_name=[],
        BC_val=[],
        nodesId=[],
    ):
        """
        update_mesh_bounds

        Parameters
        ----------
        bound_type : str
            Neumann or Dirichlet.
        bound_bool : TYPE
            Boolean for bound cond.
        """
        self.console.print(":sponge: [b]update boundary condition dataframe[/b]")
        
        if 'BC_type' not in self.mesh_bound_cond_df.columns:
            self.mesh_bound_cond_df.loc[
                                        self.mesh_bound_cond_df["time"] == time, 
                                        BC_name
                                        ] = None
        self.mesh_bound_cond_df.loc[
                                    (self.mesh_bound_cond_df["time"] == time) & 
                                    (self.mesh_bound_cond_df["id_node"].isin(nodesId)),
                                    BC_name
                                    ] = BC_val        
        pass

    def create_mesh_vtkris3d_vtk9(self):
        """
        Create mesh for vtk format version 9


        # vtk DataFile Version 9.0
        3D Unstructured Grid of Linear Triangles
        ASCII
        DATASET STRUCTURED_GRID
        FIELD FieldData  1
        TIME 1 1 double
                   0.00000
        POINTS     7056 float
        TETRA (5,NT)       - element connectivities in 3-d mesh (TETRA(5,I)
        # C                        indicates material type for 3-d element I)
        """

        if not "nnod3" in self.grid3d.keys():
            self.run_processor(IPRT1=3)

        with open(
            os.path.join(self.workdir, self.project_name, "vtk/mesh_tmp.vtk"), "w+"
        ) as vtkmesh:
            vtkmesh.write("# vtk DataFile Version 2.0\n")
            vtkmesh.write("3D Unstructured Grid of Linear Triangles\n")
            vtkmesh.write("ASCII\n")
            vtkmesh.write("DATASET UNSTRUCTURED_GRID\n")
            vtkmesh.write("FIELD FieldData  1\n")
            vtkmesh.write("TIME 1 1 double\n")
            vtkmesh.write("           0.00000\n")
            vtkmesh.write(
                "POINTS " + "{:3.0f}".format(self.grid3d["nnod3"]) + " float\n"
            )
            np.savetxt(vtkmesh, self.grid3d["mesh3d_nodes"], fmt="%1.6e")
            len(self.grid3d["mesh3d_nodes"])
            ntetra = len(self.grid3d["mesh_tetra"])
            numb = ntetra * 5
            mesh_tretra_m = self.grid3d["mesh_tetra"].T[0:4] - 1
            # np.shape(mesh_tretra_m)
            tetra_mesh = np.vstack([4 * np.ones(ntetra), mesh_tretra_m]).T
            vtkmesh.write(
                "CELLS " + "{:d}".format(ntetra) + "\t" + "{:d}".format(numb) + "\n"
            )
            np.savetxt(vtkmesh, tetra_mesh, fmt="%d", delimiter="\t")
            vtkmesh.write("CELL_TYPES " + "{:d}".format(ntetra) + "\n")
            np.savetxt(vtkmesh, 10 * np.ones(ntetra).T, fmt="%d")
            vtkmesh.close()

        pass

    def create_mesh_vtkris3d_vtk2(self, verbose):
        """
        Create mesh for vtk format version 2


        # vtk DataFile Version 2.0
        3D Unstructured Grid of Linear Triangles
        ASCII
        DATASET STRUCTURED_GRID
        FIELD FieldData  1
        TIME 1 1 double
                   0.00000
        POINTS     7056 float
        TETRA (5,NT)       - element connectivities in 3-d mesh (TETRA(5,I)
        # C                        indicates material type for 3-d element I)
        """

        with open(
            os.path.join(self.workdir, self.project_name, "vtk/mesh_tmp.vtk"), "w+"
        ) as vtkmesh:
            vtkmesh.write("# vtk DataFile Version 2.0\n")
            vtkmesh.write("3D Unstructured Grid of Linear Triangles\n")
            vtkmesh.write("ASCII\n")
            vtkmesh.write("DATASET UNSTRUCTURED_GRID\n")
            vtkmesh.write("FIELD FieldData  1\n")
            vtkmesh.write("TIME 1 1 double\n")
            vtkmesh.write("           0.00000\n")
            vtkmesh.write(
                "POINTS " + "{:3.0f}".format(self.grid3d["nnod3"]) + " float\n"
            )
            np.savetxt(vtkmesh, self.grid3d["mesh3d_nodes"], fmt="%1.6e")
            len(self.grid3d["mesh3d_nodes"])
            ntetra = len(self.grid3d["mesh_tetra"])
            numb = ntetra * 5
            mesh_tretra_m = self.grid3d["mesh_tetra"].T[0:4] - 1
            # np.shape(mesh_tretra_m)
            tetra_mesh = np.vstack([4 * np.ones(ntetra), mesh_tretra_m]).T
            vtkmesh.write(
                "CELLS " + "{:d}".format(ntetra) + "\t" + "{:d}".format(numb) + "\n"
            )
            np.savetxt(vtkmesh, tetra_mesh, fmt="%d", delimiter="\t")
            vtkmesh.write("CELL_TYPES " + "{:d}".format(ntetra) + "\n")
            np.savetxt(vtkmesh, 10 * np.ones(ntetra).T, fmt="%d")
            vtkmesh.close()

        pass

    def create_mesh_vtk(self, verbose=False):
        """
        Create custum mesh
        THIS SHOULD BE MOVED TO MESHTOOLS
        """

        # if not 'nnod3' in self.grid3d.keys():
        self.run_preprocessor(verbose=verbose)
        self.run_processor(IPRT1=3, verbose=verbose)
        
        self.create_mesh_vtkris3d_vtk2(verbose)
        # self.create_mesh_vtkris3d_vtk9()
        self.mesh_pv_attributes = pv.read(
            os.path.join(self.workdir, self.project_name, "vtk/mesh_tmp.vtk")
        )
        self.mesh_pv_attributes.save(
            os.path.join(
                self.workdir, self.project_name, "vtk/", self.project_name + ".vtk"
            ),
            binary=False,
        )

        # THIS IS A TEMPORARY IMPLEMENTATION OF THE PYVISTA MESH
        #  https://docs.pyvista.org/examples/00-load/create-structured-surface.html
        pass

    def update_mesh_vtk(self, prop="", prop_value=[], replaceVTK=True, **kwargs):
        """
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

        """
        self.mesh_pv_attributes.add_field_data(prop_value, prop)

        if replaceVTK:
            self.mesh_pv_attributes.save(
                os.path.join(
                    self.workdir, self.project_name, "vtk/", self.project_name + ".vtk"
                ),
                binary=False,
            )

        pass

    def map_prop_2mesh_markers(self, 
                               prop_name, 
                               prop_map, 
                               zones_markers_3d=None,
                               to_nodes=False,
                               **kwargs):
        """
        Map a physical property to the CATHY mesh nodes/cells.
        The mapping length should be equal to the mesh node markers 
        unique value length.
        If no markers (i.e zones_markers_3d = None) are defined 
        for the mesh nodes then each layers 
        is associated with a unique marker

        Parameters
        ----------
        prop_name : str
            property name i.e. ic, POROS, ... .
        prop_map : list
            Values of the property. 
        to_nodes : bool, optional
            Map to the mesh nodes. The default is False.



        
        Returns
        -------
        pv.Mesh
            Updated pyvista mesh with new property.
        """
        
        if hasattr(self, "mesh_pv_attributes") == False:
            self.create_mesh_vtk()
        
        saveMeshPath = os.path.join(
                                    self.workdir, 
                                    self.project_name, 
                                    "vtk/", 
                                    self.project_name + ".vtk"
                                    )
        if 'saveMeshPath' in kwargs:
            saveMeshPath = kwargs.pop('saveMeshPath')
        # if 'node_markers_new' in self.mesh_pv_attributes == False:
            
        if zones_markers_3d is None:
            zones_markers_3d = []
            for l in range(self.dem_parameters['nstr']):
                zones_markers_3d.append(np.ones([self.hapin["M"], self.hapin["N"]])*l)
            
            mt.add_markers2mesh(
                                zones_markers_3d,
                                self.DEM,
                                self.mesh_pv_attributes,
                                self.dem_parameters,
                                self.hapin,
                                self.grid3d,
                                to_nodes=False,
                                show=False,
                                saveMeshPath = saveMeshPath
                            )
            
        
        if to_nodes:
            prop_mesh_nodes = np.zeros(len(self.mesh_pv_attributes["node_markers_new"]))
            for m in range(len(prop_map)):
                prop_mesh_nodes[
                                self.mesh_pv_attributes["node_markers_new"] == m
                                ] = prop_map[m]
            self.mesh_pv_attributes[prop_name] = prop_mesh_nodes
            self.mesh_pv_attributes.save(saveMeshPath,
                                            binary=False,
                                            )
            return list(prop_mesh_nodes)

        else:
            prop_mesh_cells = np.zeros(len(self.mesh_pv_attributes["cell_markers"]))
            for m in range(len(prop_map)):
                prop_mesh_cells[
                    self.mesh_pv_attributes["cell_markers"] == m
                ] = prop_map[m]
            self.mesh_pv_attributes[prop_name] = prop_mesh_cells
            self.mesh_pv_attributes.save(saveMeshPath,
                                            binary=False,
                                            )
            self.mesh_pv_attributes.set_active_scalars(prop_name)
            prop_mesh_nodes = self.mesh_pv_attributes.cell_data_to_point_data()

            return prop_mesh_cells, prop_mesh_nodes[prop_name]

    def map_prop2mesh(self, dict_props):
        """
        Add a given physical property to the CATHY mesh
        """
        if hasattr(self, "mesh_pv_attributes") == False:
            self.create_mesh_vtk()

        for dp in dict_props.keys():
            if ~isinstance(dict_props[dp], list):
                print(
                    "Single value detected for "
                    + str(dp)
                    + " ==> assumming it homogeneous"
                )
                self.update_mesh_vtk(
                    prop=dp,
                    prop_value=np.ones(len(self.mesh_pv_attributes.points))
                    * dict_props[dp],
                )
            else:
                self.update_mesh_vtk(prop=dp, prop_value=dict_props[dp])
        pass

    def map_prop_veg(self, dict_props):
        if hasattr(self, "veg_map") == False:
            warnings.warn("no known existing vegetation map.")
            pass

        else:
            map_veg_dict = {}
            update_map_veg = {}

            for d in dict_props.keys():
                map_veg = np.zeros(np.shape(self.veg_map))
                for i, value in enumerate(dict_props[d]):
                    map_veg[self.veg_map == i+1] = value
                update_map_veg[d] = map_veg

            return update_map_veg

    def map_prop2zone(self, dict_props, prop):
        if hasattr(self, "zone") == False:
            warnings.warn("no known existing zones.")
            pass
        else:
            prop_zones = np.zeros(np.shape(self.zone))
            for z in range(len(np.unique(self.zone))):
                prop_zones[self.zone == z + 1] = dict_props[prop][z]
            return prop_zones

    # ------------------------------------------------------------------------
    #%% Plot call
    # ------------------------------------------------------------------------
    def show(self, prop="hgsfdet", **kwargs):
        """
        Call and parse to cathy.plotter from the main CATHY class

        Parameters
        ----------
        prop : str
            property to plot.
        **kwargs : kwargs

        Returns
        -------
        None.

        """
        try:    
            df = self.read_outputs(filename=prop)
        except:
            pass
        if prop == "hgsfdet":
            plt_CT.show_hgsfdet(df, **kwargs)
        elif prop == "hgraph":
            plt_CT.show_hgraph(df, **kwargs)
        elif prop == "cumflowvol":
            plt_CT.show_COCumflowvol(df, **kwargs)
        elif prop == "dtcoupling":
            plt_CT.show_dtcoupling(df, **kwargs)
        elif prop == "wtdepth":
            plt_CT.show_wtdepth(df, **kwargs)
        elif prop == "spatialET":
            df_fort777 = out_CT.read_fort777(os.path.join(self.workdir,
                                                          self.project_name,
                                                          'fort.777'),
                                              )
            cmap = plt_CT.show_spatialET(df_fort777, **kwargs)
            return cmap
            
        elif prop == "WTD": # water table depth
        
            xyz = self.read_outputs('xyz')
            df_psi = self.read_outputs('psi')
            grid3d = self.read_outputs('grid3d')
            nstr = self.dem_parameters['nstr']+1
            nnod = int(grid3d['nnod'])
            NPRT = np.shape(df_psi.values)[0]
            XYZsurface= xyz[['x','y','z']].iloc[0:nnod].to_numpy()    
            WT, FLAG = self.infer_WTD_from_psi(df_psi.values,nnod,nstr,NPRT,xyz,XYZsurface)
            cmap = plt_CT.plot_WTD(XYZsurface,WT,**kwargs)
            return cmap
            
            
        else:
            print("no proxy to plot")
        # elif filename == 'psi':
        #     plt_CT.show_psi(path)
        # elif filename == 'sw':
        #     plt_CT.show_sw(path)
        pass

    def show_bc(self, BCtypName=None, time=0, ax=None, **kwargs):
        """Show bc"""
        
        if BCtypName is None:
            fig = plt.figure()
            # Plot nansfdirbc
            # ------------------------------------
            ax = fig.add_subplot(1, 3, 1, projection='3d')
            cmap = plt_CT.plot_mesh_bounds('nansfdirbc', 
                                    self.mesh_bound_cond_df, 
                                    time, 
                                    ax
                                    )
            # Plot nansfneubc
            # ------------------------------------
            ax = fig.add_subplot(1, 3, 2, projection='3d')
            cmap = plt_CT.plot_mesh_bounds('nansfneubc', 
                                    self.mesh_bound_cond_df, 
                                    time, 
                                    ax
                                    )           
            # Plot sfbc
            # ------------------------------------
            ax = fig.add_subplot(1, 3, 3, projection='3d')
            cmap = plt_CT.plot_mesh_bounds('sfbc', 
                                    self.mesh_bound_cond_df, 
                                    time, 
                                    ax
                                    )    
            plt.tight_layout()
        else:
            cmap = plt_CT.plot_mesh_bounds(BCtypName, 
                                    self.mesh_bound_cond_df, 
                                    time, 
                                    ax
                                    )           

        pass

    def show_input(self, prop="hgsfdet", ax=None, **kwargs):
        """
        Call and parse to cathy.plotter from the main CATHY class

        Parameters
        ----------
        prop : str
            property to plot.
        **kwargs : kwargs

        Returns
        -------
        None.

        """
        df = self.read_inputs(filename=prop, **kwargs)
        if prop == "atmbc":
            plt_CT.show_atmbc(df["time"], df["value"], ax=ax, **kwargs)
        elif prop == "root_map":
            plt_CT.show_indice_veg(df[0], ax=ax, **kwargs)
        elif prop == "dem":
            if hasattr(self, "hapin") is False:
                self.update_prepo_inputs()
            hapin = self.hapin
            plt_CT.show_dem(df[0], hapin, ax=ax, **kwargs)
        elif prop == "zone":
            plt_CT.show_zone(df[0], ax=ax)           
        elif prop == "soil":
            SPP_colname = self._get_soil_SPP_columnsNames()
            FP_colname = self._get_soil_FP_columnsNames()

            yprop = "PERMX"
            if "yprop" in kwargs:
                yprop = kwargs["yprop"]
                    
            if yprop in SPP_colname:
                # in 2 dimensions
                # -------------
                zone_mat = in_CT.read_zone(
                    os.path.join(self.workdir, self.project_name, "prepro/zone")
                )
    
                layer_nb = 0
                if "layer_nb" in kwargs:
                    layer_nb = kwargs["layer_nb"]
    

    
                soil_map_prop = zone_mat[0]
                
                exclude_zone = self._check_outside_DEM(zone_mat[0])
                NZONES = len(np.unique(zone_mat[0])) - exclude_zone
    
                if NZONES-1>1:
                    for z in range(NZONES):
                        soil_map_prop[zone_mat[0] == z+1] = df[0][yprop].xs(
                                                                    (z+1,layer_nb)
                                                                    )
                else:
                    soil_map_prop[zone_mat[0] == 1] = df[0][yprop].xs((1, 
                                                                       layer_nb)
                                                                      )
    
                cmap = plt_CT.show_soil(soil_map_prop, ax=ax,
                                 **kwargs)
                return cmap
                
            elif yprop in FP_colname:
                FP_map_prop = np.copy(self.veg_map)
                for vegi in range(len(np.unique(FP_map_prop))):
                    FP_map_prop[FP_map_prop == vegi+1] = df[1][yprop][vegi+1]
                    
                cmap = plt_CT.show_soil(FP_map_prop, ax=ax,
                                 **kwargs)
                return cmap

        elif ("dtm_" in prop) | ("lakes_map" in prop):
            raster_mat, header_raster = in_CT.read_raster(
                os.path.join(self.workdir, self.project_name, "prepro/" + prop)
            )
            if hasattr(self, "hapin") is False:
                self.update_prepo_inputs()
            hapin = self.hapin
            plt_CT.show_raster(
                raster_mat, header_raster, prop=prop, hapin=hapin, ax=ax, **kwargs
            )

        else:
            print("no proxy to plot")

        pass

    #%% Read outputs/inputs
    # ------------------------------------------------------------------------
    def read_outputs(self, filename, **kwargs):
        """
        Read CATHY format output file


        Parameters
        ----------
        filename : str
            name of the output file to read.

        Returns
        -------
        A dataframe or a dict describing file data/entries.

        """
        path = os.path.join(self.workdir, self.project_name, "output", filename)
        if "path" in kwargs:
            path = kwargs["path"] + "/" + filename

        if filename == "vp":
            df = out_CT.read_vp(path)
            return df
        elif filename == "hgraph":
            df = out_CT.read_hgraph(path)
            return df
        elif filename == "dtcoupling":
            df = out_CT.read_dtcoupling(path)
            return df
        elif filename == "hgsfdet":
            df = out_CT.read_hgsfdet(path)
            return df
        elif filename == "psi":
            df = out_CT.read_psi(path)
            return df
        elif filename == "sw":
            df = out_CT.read_sw(path)
            return df
        elif filename == "mbeconv":
            df = out_CT.read_mbeconv(path)
            return df
        elif filename == "cumflowvol":
            df = out_CT.read_cumflowvol(path)
            return df
        elif filename == "wtdepth":
            df = out_CT.read_wtdepth(path)
        elif filename == "xyz":
            df = out_CT.read_xyz(path)
            return df
        elif filename == "grid3d":
            df = out_CT.read_grid3d(path)
            return df
        else:
            print("no file specified")
        pass

    def read_inputs(self, filename, **kwargs):
        """
        Read CATHY format input file


        Parameters
        ----------
        filename : str
            name of the input CATHY file to read.

        Returns
        -------
        A dataframe or a dict describing file data/entries.

        """
        if filename == "atmbc":
            
            if len(self.grid3d) == 0:
                # self.run_processor(IPRT1=3, DAFLAG=0)
                self.grid3d = out_CT.read_grid3d(os.path.join(self.workdir,
                                                              self.project_name, 
                                                              'output', 'grid3d')
                                                 )
            df, HSPATM, IETO, = in_CT.read_atmbc(
                os.path.join(self.workdir, self.project_name, "input", filename),
                grid=self.grid3d
            )
            self.atmbc = {"HSPATM": HSPATM, 
                          "IETO": IETO, 
                          "time": df['time'].unique(), 
                          "VALUE": df['value'],
                          "atmbc_df": df
                          }

            return df
        elif filename == "dem":
            df = in_CT.read_dem(
                os.path.join(self.workdir, self.project_name, "prepro/dem"),
                os.path.join(self.workdir, self.project_name, "prepro/dtm_13.val"),
            )
            return df
        elif filename == "root_map":
            df = in_CT.read_root_map(
                os.path.join(self.workdir, self.project_name, "input", filename)
            )
            return df
        elif filename == "soil":
            dem_parm = in_CT.read_dem_parameters(
                os.path.join(self.workdir, self.project_name, "input", "dem_parameters")
            )
            # self.dem_parameters
            # self.update_dem_parameters()
            

            if 'MAXVEG' in kwargs:
                MAXVEG = kwargs['MAXVEG']
            else:
                MAXVEG = self.MAXVEG
                
            df = in_CT.read_soil(
                os.path.join(self.workdir, self.project_name, "input", filename),
                dem_parm,
                MAXVEG=MAXVEG,
            )
            # df[0]
            return df
        elif filename == "zone":
            df = in_CT.read_zone(
                os.path.join(self.workdir, self.project_name, "prepro/zone")
            )
            return df
        elif ("dtm_" in filename) | ("lakes_map" in filename):
            raster_mat, header_raster = in_CT.read_raster(
                os.path.join(self.workdir, self.project_name, "prepro/" + filename)
            )
            return raster_mat, header_raster
        else:
            print("unknown file requested")
        pass

    # -------------------------------------------------------------------#
    # %% utils
    # -------------------------------------------------------------------#

    
    def infer_WTD_from_psi(self,psi,nnod,nstr,NPRT,xyz,XYZsurface):
        # Vertical profiles in rows and layers in columns,
        # define topography Z
        Z = np.zeros([nnod,nstr])
        for l in range(nstr):
            for nn in range(nnod):
                Z[nn,l] = xyz['z'].iloc[l*nnod+nn]
        
        Z = np.fliplr(Z)
        FLAG = np.zeros([nnod,NPRT]);
        WT=[]
        for nprti in range(NPRT): # loop over NPRT
           
            vpPSI = np.zeros([nnod,nstr])
            for l in range(nstr):
                for nn in range(nnod):
                    vpPSI[nn,l] = psi[nprti][l*nnod+nn]
            
            vpPSI = np.fliplr(vpPSI)
            # % Z0 contains for every vertical profile the height of the watertable,
            Z0 = XYZsurface[:,2].copy()
        
            for ni in range(nnod): ## loop over all mesh nodes 
                for nstri in range(nstr-1): ## loop over mesh layers 
                    # print(nstri)
                    if vpPSI[ni, nstri] > 0 and vpPSI[ni, nstri+1] < 0 and FLAG[ni, nprti] == 0:
                        # Watertable, interpolate linearly (Z=rc*VPPSI+z0)
                        rc = (Z[ni, nstri] - Z[ni, nstri+1]) / (vpPSI[ni, nstri] - vpPSI[ni, nstri+1])
                        Z0[ni] = Z[ni, nstri] - rc * vpPSI[ni, nstri]
                        FLAG[ni, nprti] = 1
                    elif vpPSI[ni, nstri] > 0 and vpPSI[ni, nstri+1] < 0 and FLAG[ni, nprti] == 1:
                            # Second watertable
                        FLAG[ni, nprti] = 2
                    elif nstri == nstr - 2 and vpPSI[ni, nstri+1] >= 0 and FLAG[ni, nprti] == 0:
                        # Watertable not encountered and nodes still saturated
                        FLAG[ni, nprti] = 3                
                    elif nstri == nstr -2 and FLAG[ni, nprti] == 0:
                        # Watertable not encountered, and not fully saturated profile. Set watertable to lowest node
                        Z0[ni] = Z[ni, nstri - nstr + 2]
                        FLAG[ni, nprti] = 4
            
            Z0 = XYZsurface[:,2]-Z0;
            WT.append(Z0)
        
        WT = np.vstack(WT)
        
        return WT, FLAG


    def find_nearest_node(self, node_coords, grid3d=[]):
        """
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

        """
        if np.array(node_coords).ndim <= 1:
            node_coords = [node_coords]
        if len(grid3d) == 0:
            grid3d = out_CT.read_grid3d(os.path.join(self.workdir,
                                                          self.project_name, 
                                                          'output', 'grid3d')
                                             )

        closest_idx = []
        closest = []
        for i, nc in enumerate(node_coords):
            # euclidean distance
            if len(nc)==3:
                d = (
                    (grid3d["mesh3d_nodes"][:, 0] - nc[0]) ** 2
                    + (grid3d["mesh3d_nodes"][:, 1] - nc[1]) ** 2
                    # + (abs(grid3d["mesh3d_nodes"][:, 2]) - abs(nc[2])) ** 2
                    + (grid3d["mesh3d_nodes"][:, 2] - nc[2]) ** 2
                ) ** 0.5
            else:
                d = (
                    (grid3d["mesh3d_nodes"][:, 0] - nc[0]) ** 2
                    + (grid3d["mesh3d_nodes"][:, 1] - nc[1]) ** 2
                ) ** 0.5
            closest_idx.append(np.argmin(d))
            closest.append(grid3d["mesh3d_nodes"][closest_idx[i],:])
            threshold = 5e-1
            if d[np.argmin(d)] > threshold:
                self.console.print(
                    ":warning: [b]No node close to the required points![/b]"
                )

        return closest_idx, closest

    def rich_display(self, title="Star Wars Movies", **kwargs):
        """
        Describe the variable state and fate during the simulation with a rich table

        Returns
        -------
        None.

        """
        self.console.print(eval("self." + str(title)))
        pass

    def backup_simu(self):
        """
        Save a copy of the simulation for reuse within python

        Returns
        -------
        project_filename.pkl

        """
        with open(
            os.path.join(self.workdir, self.project_name, self.project_name + ".pkl"),
            "wb",
        ) as f:
            pickle.dump(CATHY, f)
        f.close()
        pass

                
    def backup_results_DA(self, meta_DA=[]):
        """
        Save minimal dataframes of the simulation for result visualisation within python

        Returns
        -------
        project_filename_df.pkl

        """
        with open(
            os.path.join(
                self.workdir, self.project_name, self.project_name + "_df.pkl"
            ),
            "wb",
        ) as f:
            pickle.dump(meta_DA, f)
            pickle.dump(self.dict_parm_pert, f)
            pickle.dump(self.df_DA, f)
            pickle.dump(self.dict_obs, f)
            if hasattr(self, "df_performance"):
                pickle.dump(self.df_performance, f)
            if hasattr(self, "df_Archie"):
                pickle.dump(self.df_Archie, f)
        f.close()

    def load_pickle_backup(self,filename=""):
        if len(filename) == 0:
            filename = os.path.join(
                self.workdir, 
                self.project_name, 
                self.project_name + "_df.pkl"
            )
            # filename = os.path.join(
            #     self.workdir, self.project_name, str(idsimu) + "_df.pkl"
            # )
        backup_list = []
        all_names = [
            "meta_DA",
            "dict_parm_pert",
            "df_DA",
            "dict_obs",
            "df_performance",
            "df_Archie",
        ]
        names = []
        i = 0
        # pickle_off = open(filename, "rb")
        # df = pd.read_pickle(pickle_off)

        with open(filename, "rb") as f:
                for i in range(len(all_names)):
                    try:
                        df = pd.read_pickle(f)
                        backup_list.append(df)
                        # backup_list.append(pickle.load(f))
                        names.append(all_names[i])
                        # i += 1
                    except:
                        pass
        f.close()
        dict_backup = {}
        for i, n in enumerate(names):
            dict_backup[n] = backup_list[i]
        return dict_backup
