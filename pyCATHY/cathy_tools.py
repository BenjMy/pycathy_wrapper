# -*- coding: utf-8 -*-
"""
@author: bmary
"""
from __future__ import print_function
import numpy as np
import os,warnings
import subprocess
import glob, os
import shutil

import git
from git import Repo

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

        # fetch src files if not existing
        try:
            Repo.clone_from('https://bitbucket.org/cathy1_0/cathy.git', 
                            os.path.join(self.dirName,self.project_name),
                             branch='master')
            print('fetch cathy src files')
            shutil.rmtree(os.path.join(self.dirName,self.project_name,'runs'))
        except:
            pass
                



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




        #'datin.f'
        # -------------

        #'CATHY.H'
        # -------------


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
        # check location of the exe (should be in project_name folder)
        os.chdir(os.path.join(self.dirName, self.project_name, 'src'))

        for file in glob.glob("*.f"):
            bashCommand = "gfortran -c " + str(file)
             
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()

        files = ""
        for file in glob.glob("*.o"):
            files += " " + str(file) 
        
        bashCommand = "gfortran" + files + " -llapack -lblas -o " + self.processor_name
      
            
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

        # move to the directory where the source FORTRAN files are contained (cathy_main.f)
        shutil.move(os.path.join(self.dirName, self.project_name, 'src',self.processor_name),
                    os.path.join(self.dirName, self.project_name, self.processor_name))

        print('move ' + self.processor_name + ' into ' + self.project_name + ' directory')
        callexe_path = os.path.join(self.dirName, self.project_name, self.processor_name)

        # 'cathy.fnames'
        # -------------
        # Some of these inputs are automatically generated during the pre-processing (Table 2),
        # while others are new input files that should be updated with appropriate parameters for the specific case study
        # The unit number associated to each input and output file (e.g. unit IIN1) are used by CATHY in reading and writing statements
        # use the backslash (\) in quotes in Windows, the forward slash (/) are used in Mac OS X.
        #'cathy_main.f'
        # -------------
        try:
            os.path.exists(os.path.join(self.dirName, self.project_name,'cathy.fnames'))
        except OSError:
            print('cathy.fnames missing')                               
            sys.exit()

        # check if input folder
        # if not os.path.exists(os.path.join(self.dirName, self.project_name,'input'):

        # check if output folder

        # check if vtk folder
        
        
        #try:
        #    subprocess.call([callexe_path, '-ARG'], shell=True)
        #    subprocess.Popen(["cathyEnv/bin/python"])
        #except:
        #    pass

        return

    def create_output(self,output_dirname='output'):
        """Short summary.

        Parameters
        ----------
        output_dirname : type
            Description of parameter `output_dirname`.

        Returns
        -------
        type
            Description of returned object.

        """

        # create project_name/output folders
        if not os.path.exists(os.path.join(self.project_name, output_dirname)):
            os.mkdir(os.path.join(self.project_name, output_dirname), mode=0o777)

        if not os.path.exists(os.path.join(self.project_name, 'vtk')):
            os.mkdir(os.path.join(self.project_name, 'vtk'), mode=0o777)


        # compute the processor (make clean make) using the Makefile
        return





    #%% INPUT FILES

    def create_inputs(self,input_dirname='input'):
        """Short summary.

        Parameters
        ----------
        input_dirname : type
            Description of parameter `input_dirname`.

        Returns
        -------
        type
            Description of returned object.

        """

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
        """Short summary.

        Parameters
        ----------
        INDP : int
            Description of parameter `INDP`.
        WTHEIGHT : type
            Description of parameter `WTHEIGHT`.
        IPOND : type
            Description of parameter `IPOND`.

        Returns
        -------
        type
            Description of returned object.

        """

        return


    def create_atmbc(HSPATM=0,IETO=0,TIME=None,VALUE=None):
        """Atmospheric forcing term (atmbc - IIN6).

        Parameters
        ----------
        HSPATM : int
            =0 for spatially variable atmospheric boundary condition
                inputs; blank or =9999 if unit IIN6 input is to be ignored;
                otherwise atmospheric BC's are homogeneous in space.
        IETO : int
            =0 for linear interpolation of the atmospheric boundary condition
                inputs between different ATMTIM; otherwise the inputs are assigned
                as a piecewise constant function (ietograph).
        TIME : np.ndarray
            Time at current time level (seconds)
        VALUE : np.ndarray
            (meters per seconds)
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

    def create_nansfdirbc(NDIR=0,NDIRC=0, NQ3=None):
        """Boundary conditions (nansfdirbc - IIN8, nansfneubc - IIN9, sfbc - IIN7)
        The boundary conditions are defined in the nansfdirbc (Dirichlet),
        nansfneubc (Neumann), and sfbc (seepage face) files.
        To simulate the no-flow boundaries conditions for the bottom and
        vertical sides of the domain it is necessary to set NDIR and NDIRC
        equal to zero.
        To simulate different boundary conditions, it is necessary to indicate
        the number of selected nodes through NDIR or NDIRC, then to specify the
        node ID’s that you want to consider and eventually the value of pressure head
         or flux that you want to assign.

        Parameters
        ----------
        NDIR : int
            Number of non-atmospheric, non‐seepage face Dirichlet nodes in 2-d mesh.
            The BC's assigned to these surface nodes are replicated vertically (compare NDIRC)
        NDIRC : int
            Number of 'fixed' non-atmospheric,
            non-seepage face Dirichlet nodes in 3‐d mesh
            ('fixed' in the sense that these BC's are not replicated to other nodes ‐ compare NDIR)
        NQ3 : int
            Number of non-atmospheric, non‐seepage face Neumann nodes in 3‐d mesh.

        Returns
        -------
        type
            Description of returned object.

        """

        with open('nansfdirbc', 'w+') as nansfdirbcfile:

            if noflow:
                nansfdirbcfile.write(str(NDIR) + "\t" + str(NDIRC) + "\t"
                                + 'NDIR' + "\t" + 'NDIRC' + "\n")

        return

    def soil(PMIN=-999,IPEAT=0,SCF=1.0,**kwargs):
        """Soil parameters (soil - IIN4).
        The porous media properties are defined in the soil file.
        The first thing that must be decides is the type of relationship
        to describe the hydraulic characteristics of the unsaturated soil (i.e. retention curves).
        This can be done through the choice of the parameter IVGHU amongst the several options.

        Parameters
        ----------
        PMIN : int
            air dry' pressure head value (for switching control of atmospheric boundary conditions during evaporation)
            [m sec]
        IPEAT : int
            Flag for peat soil deformation
            =0 constant porosity (in unsaturated soil)
            =1 consider porosity variations with water saturation
        SCF : int
            soil cover fraction (fraction of soil covered in vegetation)
        **kwargs : type
            Description of parameter `**kwargs`.

        Returns
        -------
        type
            Description of returned object.

        """

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
        """Contains the raster map of the root zone depth as shown in the figure below.

        Returns
        -------
        type
            Description of returned object.

        """

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
        """Short summary.

        Parameters
        ----------
        gmsh_mesh : type
            Description of parameter `gmsh_mesh`.
        NZONE : type
            # of material types in the porous medium.
        NSTR : type
            The number of vertical layers.
        N1 : type
            The maximum number of element connections to a node.
        NNOD : type
            Description of parameter `NNOD`.
        NTRI : type
            Description of parameter `NTRI`.
        ZRATIO : type
            The thickness of vertical layers or the fraction of total grid height
            that each layer is to occupy (ZRATIO (1) is for the surface‐most layer.
            ZRATIO values must sum to 1.).
        Z1 : type
            Description of parameter `Z1`.
        IVERT : type
            =0 each layer will be parallel to the surface, including the base of the 3‐d grid.
            `ZRATIO` is applied to each vertical cross section.
            =1 base of the 3‐d grid will be flat, and `ZRATIO` is applied to each vertical cross section
            =2 base of the 3‐d grid will be flat, as will the NSTR‐1 horizontal cross sections above it.
            `ZRATIO` is applied only to the vertical cross section having the lowest elevation.
            =3 for each cell of the dem a single depth value is read in file input IIN60 (basement).
            `ZRATIO` is applied to each vertical cross section.
            =4 the first NSTR‐1 layers from the surface will be parallel to the surface and the base of the 3‐d grid will be flat.
            `ZRATIO` is applied only to the vertical cross section having the lowest elevation.
        ISP : type
            =0 for flat surface layer (only one Z value is read in, and is replicated to all surface nodes);
            otherwise surface layer is not flat (Z values read in for each surface node);
            (for ISP=0, IVERT=0, 1, and 2 yield the same 3‐d mesh, given the same values of BASE and ZRATIO).
        BASE : type
            Value which defines the thickness or base of the 3‐d mesh.
            For `IVERT`=0, BASE is subtracted from each surface elevation value,
            so that each vertical cross section will be of thickness BASE,
            and the base of the 3‐d mesh will be parallel to the surface.
            For IVERT=1 or 2, BASE is subtracted from the lowest surface elevation value,
            say ZMIN, so that each vertical cross section will be of thickness (Z ‐ ZMIN) + BASE,
            where Z is the surface elevation for that cross section.
            The base of the 3‐d mesh will thus be flat.
        debug : type
            Description of parameter `debug`.

        Returns
        -------
        type
            Description of returned object.

        """

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
