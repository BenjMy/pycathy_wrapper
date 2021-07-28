# -*- coding: utf-8 -*-
"""
@author: bmary
"""
from __future__ import print_function
import numpy as np
import os,warnings
import subprocess
# import git
import meshtools as mt

class CATHY(object):
    '''Main CATHY object.'''
    def __init__(self,dirName):
        '''Create CATHY object.
        '''
        print('init CATHY object')
        self.dirName = dirName
        
        
        self.processor_name = 'cathy'
        self.project_name = 'my_cathy_prj'
        
        if not os.path.exists(self.project_name):
            os.mkdir(self.project_name, mode=0o777)
            
        if not os.path.exists(self.project_name + '/src'):
            os.mkdir(self.project_name + '/src', mode=0o777)
            #git.Git(elf.project_name + '/src').clone("git@bitbucket.org:cathy1_0/cathy.git")

            


        #check if src files are existing

        self.mesh_gmsh = [] # gmsh meshfile

        # check if all the folders/files are included
        # prepro
        # input
        # output (empty folder)
        # cathy.fnames and cathy.exe


    def init(self):
        '''Make CATHY input files.
        '''
        self.prepare_makefile()
        self.create_cathy_exe(exe_name=self.processor_name)
        self.create_output(output_dirname='output')

    def prepare_makefile(self, job_type=0):

        # check flags for job type  r
        # check if all inputs files are present

        # DAFLAG - flag for the choice of the data assimilation scheme:
        #           = 0 nudging (if NUDN=0, no data assimilation)
        #           = 1 EnKF with Evensen's algorithm (Ocean Dynamics, 2004)
        #           = 2 EnKF with parameters update
        #           = 3 Particle filtering (SIR algorithm)
        #           = 4 Particle filtering (SIR algorithm) with parameters update


        # Data assimilation = 2
        # Plant model = 1


        # mandatory inputs common to all type of simulation
        self.set_atmbc()
        # copy into project_name/inputs

        # inputs for Data assimilation

        # inputs for plant model


        return

    def create_cathy_exe(self):
        """Compile cathy inputs files and create cathy.exe

        Returns
        -------
        type
            Description of returned object.

        """

        # check if all files/folder are included
        # 'Makefile' # make makefile equivalent to cmd gfortran -c ... gfortran *.o ...
        # -------------

        # 'cathy.fnames'
        # -------------
        # Some of these inputs are automatically generated during the pre-processing (Table 2),
        # while others are new input files that should be updated with appropriate parameters for the specific case study
        # The unit number associated to each input and output file (e.g. unit IIN1) are used by CATHY in reading and writing statements
        # use the backslash (\) in quotes in Windows, the forward slash (/) are used in Mac OS X.
        #'cathy_main.f'
        # -------------


        #'datin.f'
        # -------------

        #'CATHY.H'
        # -------------

        # move to the directory where the source FORTRAN files are contained (cathy_main.f)

        # edit the openio.f file by replacing the forward slash (/) with back slash (\)
        # for Windows and the CATHY.H file to properly define the dimensions of the specific problem (dimension of the arrays).


        # compute the processor (make clean make) using the Makefile
        #make clean
        #make Makefile


        return

    def run_processor(self):
        """ Run cathy.exe

        Returns
        -------
        type
            Description of returned object.

        """
        #self.init()
        # check if input folder

        # check if output folder

        # check if vtk folder

        # check location of the exe (should be in project_name folder)
        callexe_path = os.path.join(self.dirName, self.project_name, self.processor_name)
        print(self.processor_name)

        #subprocess.call(["gfortran","-o","output.fx","source.f90"])#create
        #subprocess.call(["output.fx"])                             #execute

        try:
            subprocess.call([callexe_path, '-ARG'], shell=True)
            subprocess.Popen(["cathyEnv/bin/python"])
        except:
            pass

        return

    def create_output(self,output_dirname='output'):

        # create project_name/output folders
        if not os.path.exists(os.path.join(self.project_name, output_dirname)):
            os.mkdir(os.path.join(self.project_name, output_dirname), mode=0o777)

        if not os.path.exists(os.path.join(self.project_name, 'vtk')):
            os.mkdir(os.path.join(self.project_name, 'vtk'), mode=0o777)
            
            
        # compute the processor (make clean make) using the Makefile
        return





    #%% INPUT FILES

    def create_inputs(self,input_dirname='input'):

        # create project_name/output folders
        if not os.path.exists(os.path.join(self.project_name, input_dirname)):
            os.mkdir(os.path.join(self.project_name, input_dirname), mode=0o777)
        #self.create_parm()
        self.create_atmbc()
        #self.create_nansfdirbc()
        # compute the processor (make clean make) using the Makefile
        return


    def create_parm(IPRT1,DAFLAG,
                    ISIMGR, PONDH_MIN,
                    KSLOPE, TOLKSL,
                    PKRL,   PKRR,   PSEL,   PSER,
                    PDSE1L, PDSE1R, PDSE2L, PDSE2R,
                    ISFONE ,ISFCVG, DUPUIT,
                    TETAF,  LUMP,   IOPT,
                    NLRELX, OMEGA,
                    L2NORM, TOLUNS, TOLSWI,  ERNLMX,
                    ITUNS, ITUNS1, ITUNS2,
                    ISOLV,  ITMXCG, TOLCG,
                    DELTAT, DTMIN,  DTMAX,  TMAX,
                    DTMAGA, DTMAGM, DTREDS, DTREDM,
                    IPRT,   NPRT,   TIMPRT ,I):
        
        with open('parm', 'w+') as parmfile:

            parmfile.write(str(NDIR) + "\t" + str(NDIRC) + "\t" 
                                + 'NDIR' + "\t" + 'NDIRC' + "\n")
            
            
            parmfile.close()
            
        return
    
    def create_ic(INDP,WTHEIGHT,IPOND):
        
        return
        
        
    def create_atmbc(HSPATM=0,IETO=0):
        """Short summary.

        Parameters
        ----------
        HSPATM : type
            =0 for spatially variable atmospheric boundary condition
                inputs; blank or =9999 if unit IIN6 input is to be ignored;
                otherwise atmospheric BC's are homogeneous in space.
        IETO : type
            =0 for linear interpolation of the atmospheric boundary condition
                inputs between different ATMTIM; otherwise the inputs are assigned
                as a piecewise constant function (ietograph).

        Returns
        -------
        type
            Description of returned object.

        """
        
        with open('atmbc', 'w+') as atmbcfile:
            atmbcfile.write(str(HSPATM) + "\t" + str(IETO) + "\t" 
                            + 'HSPATM' + "\t" + 'IETO' + "\n")
            atmbcfile.close()        
        
        return
    
    def create_nansfdirbc(HSPATM=0,IETO=0, noflow=False):
        
        with open('nansfdirbc', 'w+') as nansfdirbcfile:

            if noflow:
                nansfdirbcfile.write(str(NDIR) + "\t" + str(NDIRC) + "\t" 
                                + 'NDIR' + "\t" + 'NDIRC' + "\n")

        return
       
    def soil(PMIN=-999,IPEAT=0,SCF=1.0,**kwargs):
             
             # SMCREF,SMCWLT, SMCANA,
             # IVGHU,
             # CBETA0, THETA0, CANG,
             # VGN, VGM, VGRMC, VGPSAT,
             # HUN, HUA, HUB, HUALFA, HUBETA, HUGAMA, HUPSIA, HUSWR,
             # BCBETA, BCRMC, BCPSAT,
             # ZONE_VERT,
             # STR_ZON)
        
        with open('soil', 'w+') as soilfile:
            soilfile.write(str(NDIR) + "\t" + str(NDIRC) + "\t" 
                                 + 'NDIR' + "\t" + 'NDIRC' + "\n")           
                         
        pass
    
    def root_map():
        
        pass
    
    def plant():
        # plant parameters only exist for CATHYv Manoli
        

        pass
    
    #%% Meshtool functions 

    def create_3dmesh_CATHY(self,gmsh_mesh=[],
                            NZONE=[],NSTR=[],N1=[],
                            NNOD=None, NTRI=None,
                            ZRATIO=[],Z1=[],
                            IVERT=0, ISP=0, BASE=4,
                            debug=False):
        
        mesh_dict, con_matrix, phys_entity = mt.mshParse(gmsh_mesh)

        # insert all the elements of the superficial mesh i.e.
        # the headers, the 4th columns desribing the triangles and the nodes coordinates


        if NNOD is None:
            print('NNOD')

        print(ZRATIO)

        with open('grid', 'w+') as gridfile:
            gridfile.write("\t"  + str(NZONE) + "\t" + str(NSTR) + "\t" + str(N1) + "\t"
                            + 'NZONE' + "\t" + 'NSTR' + "\t" + 'N1' + "\n")
            gridfile.write("\t"  + str(NNOD) + "\t" + str(NTRI) + "\t"
                            + 'NNOD' + "\t" + 'NTRI' + "\n")
            gridfile.write("\t"  + str(IVERT) + "\t" + str(ISP) + "\t" + str(BASE) + "\t"
                            + 'IVERT' + "\t" + 'ISP' + "\t" + 'BASE (m) ' + "\n")
            #gridfile.write("\t" + '\t'.join(ZRATIO[1:]) +
            #                + 'ZRATIO' + "\n")
            np.savetxt('grid', list(ZRATIO))
            gridfile.write('ZRATIO'+ "\n")

            gridfile.write("\t"  + str(Z1) + "\t"
                            + 'Z(1) m' + "\n")

            #np.savetxt('grid', np.transpose([mesh_dict['node_x'], mesh_dict['node_y']]),fmt='%1.2f')
            np.savetxt('grid', np.transpose([mesh_dict['node_x'], mesh_dict['node_y']]),fmt='%1.2f')
            gridfile.close()

        for k, v in mesh_dict.items():    # TypeError: 'list' object is not callable
            print(k)
        #print(mesh_dict['elm_id'])
        #gridfile.write("durin's day\n")
        # copy mesh attribute to 'grid' files


    #%% DATA ASSIMILATION

    def create_archie():
        pass
    
    def create_data_ass():
        
        pass
        
    def create_elec_nodes():

        pass


    #%% ERT DATA
