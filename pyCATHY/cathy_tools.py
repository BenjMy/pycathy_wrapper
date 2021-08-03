# -*- coding: utf-8 -*-
"""
@author: bmary
"""
from __future__ import print_function
import sys
import numpy as np
import os,warnings
import subprocess
import glob, os
from os import listdir
from os.path import isfile, join
import shutil

import git
from git import Repo

#import resipy
#from resipy import Project

import meshtools as mt

class CATHY(object):
    '''Main CATHY object.'''
    def __init__(self,dirName,prjName='my_cathy_prj',**kwargs):
        '''Create CATHY object.
        '''
        print('init CATHY object')
        self.workdir = os.path.join(os.getcwd() , dirName)
        os.chdir(self.workdir)


        self.processor_name = 'cathy'
        self.project_name = prjName
        self.input_dirname = 'input'
        self.output_dirname = 'output'


        # self.time = []
        # infitration
        self.drippers = []

        # ERT
        self.elecs = []



        for key,value in kwargs.items():
            if key == 'clear_src':
                if value == True:
                    if os.path.exists(os.path.join(self.workdir,self.project_name,'src')):
                        print('clear src files')
                        shutil.rmtree(os.path.join(self.workdir,self.project_name,'src'))
                        shutil.rmtree(os.path.join(self.workdir,self.project_name,'input'))
                        shutil.rmtree(os.path.join(self.workdir,self.project_name,'tmp_src'))
                        shutil.rmtree(os.path.join(self.workdir,self.project_name,'prepro'))
                

        if not os.path.exists(os.path.join(self.project_name)):
            os.makedirs(os.path.join(self.project_name),exist_ok=True)

        # if not os.path.exists(os.path.join(self.project_name,'prepro')):
        #     os.makedirs(os.path.join(self.project_name,'prepro'),exist_ok=True)

        # fetch src files if not existing
        if not os.path.exists(os.path.join(self.project_name,'src')):
            print('src files not found')
            # try:
            Repo.clone_from('https://bitbucket.org/cathy1_0/cathy.git',
                            os.path.join(self.workdir,self.project_name,'tmp_src'),
                             branch='master')
            print('fetch cathy src files')
            shutil.move(os.path.join(self.workdir,self.project_name,'tmp_src/src'),
                        os.path.join(self.workdir,self.project_name,'src'))

            print('fetch cathy prepro src files')
            shutil.move(os.path.join(self.workdir,self.project_name,'tmp_src/runs/weilletal/prepro'),
                        os.path.join(self.workdir,self.project_name,'prepro'))

            print('fetch cathy input files')
            shutil.move(os.path.join(self.workdir,self.project_name,'tmp_src/runs/weilletal/input'),
                        os.path.join(self.workdir,self.project_name,'input'))


            pathsrc = os.path.join(os.getcwd(),self.project_name,'tmp_src/runs/weilletal/')


            onlyfiles = [f for f in listdir(pathsrc) if isfile(join(pathsrc, f))]

            for file in onlyfiles: # You could shorten this to one line, but it runs on a bit.
                shutil.move(os.path.join(pathsrc, file),
                            os.path.join(self.project_name,file))

        

                #shutil.rmtree(os.path.join(self.project_name,'tmp_src/runs'))

                # shutil.move(os.path.join(self.dirName,self.project_name,'tmp_src/README'),
            #             os.path.join(self.dirName,self.project_name,'README'))
            # except:
            #     pass

        for key,value in kwargs.items():
            if key == 'clear_outputs':
                if value == True:
                        if not os.path.exists(os.path.join(self.workdir,self.project_name,'output')):
                            self.create_output(output_dirname='output')
                        else:
                            shutil.rmtree(os.path.join(self.workdir,self.project_name,'output'))
                            shutil.rmtree(os.path.join(self.workdir,self.project_name,'vtk'))
                            self.create_output(output_dirname='output')


                    
                    
                    

        #check if src files are existing

        self.mesh_gmsh = [] # gmsh meshfile

        # check if all the folders/files are included
        # prepro
        # input
        # output (empty folder)
        # cathy.fnames and cathy.exe


    # def init(self):
    #     '''Make CATHY input files.
    #     '''
        # self.prepare_makefile()
        # self.create_cathy_exe(exe_name=self.processor_name)
        # self.create_output(output_dirname='output')



    def update_prepo_inputs(self,DEM=None,verbose=False,**kwargs):
        """ Update default prepro inputs i.e. hap.in and dtm_13.val files based on kwargs

        Parameters
        ----------
        DEM : type
            Description of parameter `DEM`.
        verbose : type
            Description of parameter `verbose`.
        **kwargs : type
            Description of parameter `**kwargs`.

        Returns
        -------
        type
            Description of returned object.

        """


        #%% hap.in
        
        print('update hap.in')

        structural_parameter = ['delta_x','delta_y','N','M','N_celle',
                               'xllcorner','yllcorner']

        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! need to update rivulet_parameter')


        hap_file = open(os.path.join(self.project_name , 'prepro/hap.in'), 'r')
        Lines = hap_file.readlines()

        # read current values in hap.in

        tmp_param_value = []
        tmp_lnb = [] # save line nb where values are existing
        count=0
        for line in Lines:
            x = line.strip().split("=")
            if len(x)>1: # skip lines with no =
                try:
                    tmp_param_value.append(float(x[1]))
                    tmp_lnb.append(count)
                    count += 1
                except:
                    pass
        # iterate over the lines
        # check from kwargs if values are changed
        # save new line with new value if needed

        L=[]
        tmp_param_value_new=[]
        count=0
        for i, line in enumerate(Lines):
            xnew = line.strip().split("=")
            if len(xnew)>1: # skip lines with no =
                for key,value in kwargs.items():
                    if count<len(structural_parameter): # only consider structural parameters
                        if key==structural_parameter[tmp_lnb[count]]:
                            print(f'key: {key} | value: {value}')
                            if value!=tmp_param_value[count]:
                                xnew[1]=value
                                line= xnew[0] + '=              ' + str(xnew[1])  + '\n'
                        tmp_param_value_new.append(value)

                count += 1 # count line nb
            L.append(line)

        # write the new hap.in file
        hap_file.close()

        hap_file = open(os.path.join(self.project_name , 'prepro/hap.in'), 'w+')
        hap_file.writelines(L)
        hap_file.close()


        self.hapin = {}

        for i in range(len(structural_parameter)):
            key = structural_parameter[i]
            self.hapin[key] = tmp_param_value_new[i]

        #%% dtm_13.val
        # If we start with a DEM file ("dtm_13.val") for an already delineated
        # catchment (i.e., a "catchment" DEM file instead of a "full" DEM file), then
        # only the last step in the above procedure (i.e., run "cppp" just once) is
        # needed (make sure that "Boundary channel construction" is set to 0 in
        # "hap.in").


        if DEM is not None:
            self.DEM = DEM

            # check consistency with hapin

            # check presence of the outlet
            if len(np.unique(DEM))==1:
                 print("Error: outlet not defined")
                 DEM_withOutlet = DEM
                 DEM_withOutlet[0,0]=0
                 
            print('update dtm_13.val')
            print(os.getcwd())
            print(DEM_withOutlet)
            with open(os.path.join(self.project_name,'prepro/dtm_13.val'), 'w+') as f:
                np.savetxt(f, DEM_withOutlet, fmt='%1.4e')   # use exponential notation


        self.update_dem_parameters(**kwargs)
        # self.update_transp(**kwargs)



        return

    def update_dem_parameters(self,**kwargs):
        """Short summary.

        Parameters
        ----------
        **kwargs : type
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

        Returns
        -------
        type
            Description of returned object.

        """
        
        print('update dem paramaters')
        # set default parameters
        
        ltmp= [0.002,	0.004,	0.006,	0.008,	0.01,	0.01,	0.02, 0.02,
               0.05,	0.05,	0.1,	0.1,	0.2, 0.2, 0.22]
        
        self.dem_parameters =	{
                      "delta_x": self.hapin['delta_x'],
                      "delta_y": self.hapin['delta_y'],
                      "factor": 1.0e+0,
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

        for keykwargs,value in kwargs.items():
           if keykwargs=='zratio':
                key = 'zratio(i),i=1,nstr'
                if sum(value)!=1:
                    print('sum is not equal to 1 -->' + str(sum(value)))
           else:
                key = keykwargs

           try:
               self.dem_parameters[key]
               self.dem_parameters[key]=value
           except:
               pass



        # write file
        header_fmt = [1,1,1,1,3,3,1]

        for key,value in self.dem_parameters.items():
            if isinstance(value, list):
                strlst = "\t".join(str(e) for e in value)
                print(strlst)
                self.dem_parameters[key]=strlst
                
                
                
        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'dem_parameters'), 'w+') as dem_parametersfile:
            
            # for key,value in self.dem_parameters.items():
            #     dem_parametersfile.write(str(value) + "\n")
                
                counth=0
                for h in header_fmt:
                    if h==3:
                        dem_parametersfile.write(str(list(self.dem_parameters.values())[counth])+ "\t" + 
                                       str(list(self.dem_parameters.values())[counth+1]) + "\t" + 
                                       str(list(self.dem_parameters.values())[counth+2]) + "\t" + "\n")
                        counth += 3
                    if h==1:
                        dem_parametersfile.write(str(list(self.dem_parameters.values())[counth])+ "\t" +"\n")
                        counth += 1
                        
                        
            # for key,value in self.dem_parameters.items():
            #     dem_parametersfile.write(str(value) + "\n")
                
                counth=0
                for h in header_fmt:
                    if h==3:
                        dem_parametersfile.write(str(list(self.dem_parameters.keys())[counth])+ "\t" + 
                                       str(list(self.dem_parameters.keys())[counth+1]) + "\t" + 
                                       str(list(self.dem_parameters.keys())[counth+2]) + "\t" + "\n")
                        counth += 3
                    if h==1:
                        dem_parametersfile.write(str(list(self.dem_parameters.keys())[counth])+ "\t" +"\n")
                        counth += 1


        # with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'dem_parameters'), 'w+') as dem_parametersfile:
        #     for key,value in self.dem_parameters.items():
        #         dem_parametersfile.write(str(value) + "\n")
        #     for key,value in self.dem_parameters.items():
        #         dem_parametersfile.write(key + "\n")

        dem_parametersfile.close()


    def run_preprocessor(self,KeepOutlet=True,verbose=False,**kwargs):
        """Run cppp.exe

        Returns
        -------
        type
            Running the executable file has allowed to generate a complete
            set of files describing physiographic features of the drainage system, as shown in Table 2.

        """

        os.chdir(os.path.join(self.project_name, 'prepro/src/'))

        #clean all files compiled
        for file in glob.glob("*.o"):
            os.remove(file)


        bashCommand = 'gfortran -O -o cppp mpar.f90 mbbio.f90 wbb_sr.f90 csort.f90 qsort.f90 depit.f90 cca.f90 smean.f90 dsf.f90 facet.f90 hg.f90 mrbb_sr.f90 bb2shp_sr.f90 shape.f90 dbase.f90 streamer.f90 cppp.f90'

        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

        if verbose:
            output, error = process.communicate()



        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

        if verbose:
            output, error = process.communicate()



        os.chdir(self.workdir)

        # move to the directory where the source FORTRAN files are contained (cathy_main.f)
        shutil.move(os.path.join(self.workdir,self.project_name, 'prepro/src/cppp'),
                    os.path.join(self.workdir,self.project_name, 'prepro/cppp'))


        print('run preprocessor')
        os.chdir(os.path.join(self.workdir, self.project_name, 'prepro'))

        bashcmd = './cppp'
        try:
             p = subprocess.Popen([bashcmd],
                                  stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  bufsize=1,
                                  universal_newlines=True)
        except Exception as e:
             print("Error: "+str(e))
             sys.exit(-1)

        p.stdin.write('2\n')
        p.stdin.flush()
        p.stdin.write('0\n')
        p.stdin.flush()
        p.stdin.write('1\n')
        p.stdin.flush()


        if  KeepOutlet==False:
            print('remove outlet')
            with open(os.path.join(self.workdir,self.project_name,'prepro/dtm_13.val'), 'w+') as f:
                    np.savetxt(f, self.DEM , fmt='%1.4e')   # use exponential

        os.chdir(self.workdir)


        return

    def run_processor(self,verbose=False,**kwargs):
        """ Run cathy.exe

        Returns
        -------
        type
            Description of returned object.

        """

        # # call subroutines from kwargs here
        #         for key,value in kwargs.items():
        #             if count<len(structural_parameter):
        #                 if key==structural_parameter[tmp_lnb[count]]:
        #                     print(f'key: {key} | value: {value}')

        # print('kwargs.items()')
        # print(kwargs.items())
            
            
        if len(kwargs.items())>0:
            print('kwargs.items()')
            print(kwargs)
        
            self.update_parm(**kwargs)
            
        # for kk in kwargs.items():
        #     print(kk)
        #     self.update_parm(kk)
        
        
        # check location of the exe (should be in project_name folder)
        print(os.getcwd())
        os.chdir(os.path.join(self.workdir,self.project_name,'src'))



        #clean all files compiled
        for file in glob.glob("*.o"):
            os.remove(file)

        for file in glob.glob("*.f"):
            bashCommand = "gfortran -c " + str(file)

            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            if verbose:
                output, error = process.communicate()

        files = ""
        for file in glob.glob("*.o"):
            files += " " + str(file)

        bashCommand = "gfortran" + files + " -llapack -lblas -o " + self.processor_name
        # print(bashCommand)


        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

        if verbose:
            output, error = process.communicate()

        # move to the directory where the source FORTRAN files are contained (cathy_main.f)
        # and OVERWRITE IF ALREADY EXISTING (full path spacified)
        os.chdir(os.path.join(self.workdir))
        shutil.move(os.path.join(self.project_name, 'src',self.processor_name),
                    os.path.join(self.project_name, self.processor_name))

        print('move ' + self.processor_name + ' into ' + self.project_name + ' directory')

        


        # 'cathy.fnames'
        # -------------
        # Some of these inputs are automatically generated during the pre-processing (Table 2),
        # while others are new input files that should be updated with appropriate parameters for the specific case study
        # The unit number associated to each input and output file (e.g. unit IIN1) are used by CATHY in reading and writing statements
        # use the backslash (\) in quotes in Windows, the forward slash (/) are used in Mac OS X.
        #'cathy_main.f'
        # -------------

        fnames_file = open(os.path.join(self.project_name , 'cathy.fnames'), 'r')
        Lines = fnames_file.readlines()
        
        # count=0
        # for line in Lines:
        #     x = line.strip().split("unit")
        #     if len(x)>1:
        #         # check if file exist
        #         # print(str(x[0]))
        #         # x[0]).replace('', '')
        #         # x0 = x[0].replace(' '' ', '')
        #         # print(str(x0))
        #         if not os.path.exists(x[0]):
        #             print(str(x[0]))





        os.chdir(os.path.join(self.project_name))

        # try:
        #     os.path.exists(os.path.join('cathy.fnames'))
        # except OSError:
        #     print('cathy.fnames missing')
        #     sys.exit()

        # check if input folder
        try:
            os.path.exists(os.path.join('inputs'))
        except OSError:
            print('input folder missing')
            sys.exit()

        # check if output folder
        if not os.path.exists('output'):
            os.makedirs('output',exist_ok=True)


        # check if vtk folder
        if not os.path.exists('vtk'):
            os.makedirs('vtk',exist_ok=True)


        print('run processor')
        print(os.getcwd())

        callexe = './' + self.processor_name
        process = subprocess.Popen(callexe.split(), stdout=subprocess.PIPE)
        if verbose:
            output, error = process.communicate()
        os.chdir(os.path.join(self.workdir))

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
            os.mkdir(os.path.join(self.workdir,self.project_name, output_dirname))

        if not os.path.exists(os.path.join(self.project_name, 'vtk')):
            os.mkdir(os.path.join(self.workdir,self.project_name, 'vtk'))


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


    def update_parm(self,**kwargs):


        # set default parameters
        self.parm =	{
                      "IPRT1": 2,
                      "NCOUT": 0,
                      "TRAFLAG": 1,
                      "ISIMGR": 2,
                      "PONDH_MIN": 0.00,
                      "VELREC": 0,
                      "KSLOPE": 0,
                      "TOLKSL": 0.01,
                      "PKRL": -3.0,
                      "PKRR": -1.0,
                      "PSEL": -3.0,
                      "PSER": -1.0,
                      "PDSE1L": -3.0,
                      "PDSE1R": -1.0,
                      "PDSE2L": -3.0,
                      "PDSE2R": -1.0,
                      "ISFONE": 0,
                      "ISFCVG": 0,
                      "DUPUIT": 0,
                      "TETAF": 1.,
                      "LUMP": 1,
                      "IOPT": 1,
                      "NLRELX": 0,
                      "OMEGA": 0.8,
                      "L2NORM": 0,
                      "TOLUNS": 1.0e-4,
                      "TOLSWI": 1.0e+30,
                      "ERNLMX": 1.0e+30,
                      "ITUNS":  10,
                      "ITUNS1": 5,
                      "ITUNS2": 7,
                      "ISOLV":  10,
                      "ITMXCG": 5,
                      "TOLCG": 7,
                      "DELTAT":  0.01,
                      "DTMIN": .00001,
                      "DTMAX": 10.,
                      "TMAX":  3600.0, #Time at end of simulation (TMAX is set to 0.0 for steady state problem)
                      "DTMAGA":  0.0,
                      "DTMAGM": 1.1,
                      "DTREDS": 0.0,
                      "DTREDM": .5,
                      "IPRT": 4,
                      "VTKF":7 ,
                      "NPRT": 3,
                      "(TIMPRT(I),I=1,NPRT)": [1800.,3600.,7200.],
                      "NUMVP": 1,
                      "(NODVP(I),I=1,NUMVP)": [441],
                      "NR": 0,
                      "NUM_QOUT": 0,
                      "(ID_QOUT(I),I=1,NUM_QOUT)": [441]
                      }

        # create dictionnary from kwargs

        for keykwargs,value in kwargs.items():
            print(f'keykwargs: {keykwargs} | value: {value}')
            if keykwargs=='TIMPRTi':
                key = '(TIMPRT(I),I=1,NPRT)'
                self.parm[key]=value
            else:
                self.parm[keykwargs]=value



        # write file
        header_fmt = [3,3,2,4,4,3,3,2,4,3,3,4,4,4,2,1,2]
        counth=0

        for key,value in self.parm.items():
            if isinstance(value, list):
                strlst = '\n '.join(str(e) for e in value)
                self.parm[key]=strlst
            if key == 'NUMVP':
                self.parm[key]= str(value) + " \n"



        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'parm'), 'w+') as parmfile:
            for h in header_fmt:
                if h==4:
                    parmfile.write(str(list(self.parm.values())[counth])+ "\t" + str(list(self.parm.values())[counth+1]) +
                                   "\t" + str(list(self.parm.values())[counth+2]) + "\t" + str(list(self.parm.values())[counth+3]) +
                                   "\t" + str(list(self.parm.keys())[counth]) + "\t" + str(list(self.parm.keys())[counth+1]) +
                                   "\t" + str(list(self.parm.keys())[counth+2]) + "\t" + str(list(self.parm.keys())[counth+3]) + "\n")
                    counth += 4
                if h==3:
                    parmfile.write(str(list(self.parm.values())[counth])+ "\t" + str(list(self.parm.values())[counth+1]) + "\t" + str(list(self.parm.values())[counth+2]) + "\t" +
                                        str(list(self.parm.keys())[counth]) + "\t" + str(list(self.parm.keys())[counth+1]) + "\t" + str(list(self.parm.keys())[counth+2]) + "\n")
                    counth += 3

                if h==2:
                    parmfile.write(str(list(self.parm.values())[counth])+ "\t" + str(list(self.parm.values())[counth+1]) + "\t" +
                                        str(list(self.parm.keys())[counth])+ "\t" + str(list(self.parm.keys())[counth+1]) + "\n")
                    counth += 2
                if h==1:
                    parmfile.write(str(list(self.parm.values())[counth])+ "\t"  +
                                        str(list(self.parm.keys())[counth])+ "\n")
                    counth += 1

        parmfile.close()

        pass





    def update_ic(self,INDP=2,IPOND=0,WTHEIGHT=[]):
        """Short summary.

        Parameters
        ----------
        INDP : int
            Flag for pressure head initial conditions (all nodes)
            =0 for input of uniform initial conditions (one value read in)
            =1 for input of non-uniform IC's (one value read in for each node)
            =2 for calculation of fully saturated vertical hydrostatic equilibrium IC's (calculated in subroutine ICVHE). In the case of IPOND>0, the fully saturated hydrostatic IC is calculated (in subroutine ICVHEPOND) starting from the ponding head values at the surface nodes, rather than surface pressure heads of 0.
            =3 for calculation of partially saturated vertical hydrostatic equilibrium IC's (calculated in subroutine ICVHWT) with the water table height (relative to the base of the 3‐d grid) given by parameter WTHEIGHT 
            =4 for calculation of partially saturated vertical hydrostatic equilibrium IC's (calculated in subroutine ICVDWT) with the water table depth (relative to the surface of the 3‐d grid) given by parameter WTHEIGHT 
        WTHEIGHT : type
            For the case INDP=3, specifies the initial water table height relative to the base of the 3‐d grid
        IPOND : type
            Flag for ponding head initial conditions (surface nodes)
            =0 no input of ponding head initial conditions; otherwise (IPOND = 1 or 2) ponding head initial conditions are read into PONDNOD, and, where PONDNOD > 0, these values are used to update the surface node values in PTIMEP read in according to the previous INDP flag
            =1 uniform ponding head initial conditions (one value read in)
            =2 non-uniform ponding head initial conditions (one value read in for each node)

        Returns
        -------
        type
            Description of returned object.

        """

        # set default parameters
        self.ic =	{
                      "INDP": INDP,
                      "WTHEIGHT": WTHEIGHT,
                      "IPOND": IPOND
                      }
        
        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'ic'), 'w+') as icfile:
            icfile.write(str(INDP) + "\t" + str(IPOND) + "\t"
                            + 'INDP' + "\t" + 'IPOND' + "\n")
            icfile.write(str(WTHEIGHT) + "\t" +  'WTHEIGHT' + "\n")                         
        icfile.close()

        pass



    def update_atmbc(self, HSPATM=0,IETO=0,TIME=None,VALUE=None):
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
        
        # set default parameters
        self.atmbc =	{
                      "HSPATM": HSPATM,
                      "IETO": IETO,
                      "TIME": TIME,
                      "VALUE": VALUE
                     }

        # C     Write atmbc file
        #       write(32,*) '1  1  HSPATM IETO'
        #          write(32,*) 0.D0, 'TIME'
        #          write(32,*) 0.D0
        #          write(32,*) 1e+20, 'TIME'
        #          write(32,*) 0.D0
        # c      do i=1,nt
        # c         write(32,*) time(i), 'TIME'
        # c         write(32,*) atmbc(i)
        # c      enddo

        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'atmbc'), 'w+') as atmbcfile:
            atmbcfile.write(str(HSPATM) + "\t" + str(IETO) + "\t"
                            + 'HSPATM' + "\t" + 'IETO' + "\n")
            for t,v in zip(TIME,VALUE):
                print(t,v)
                atmbcfile.write(str(t) + 
                                 "\t" + 'TIME' + "\n")
                atmbcfile.write(str(v) + 
                                 "\t" + 'VALUE' + "\n")
                 
        atmbcfile.close()
        pass

    def update_nansfdirbc(NDIR=0,NDIRC=0, NQ3=None):
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
        
        self.read_grid3d()

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
              
      
        # set default parameters
        # self.nansfdirbc =	{
        #               "NDIR": HSPATM,
        #               "NDIRC": IETO,
        #              }

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

    def read_grid3d(self):
      
         grid3d
         hap_file = open(os.path.join(self.project_name , 'output/grid3d'), 'r')
         Lines = hap_file.readlines()
            
         #    c     Read grid3d
         #  read(28,*) nnod,nnod3,nel
         #  do i=1,nel
         #     read(28,*) 
         #  enddo
         #  a=0
         #  do i=1,nnod3
         #     read(28,*) x(i),y(i),z(i)
         #     if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
         # 1       (y(i).eq.5))then
         #     a=a+1
         #     endif
         #  enddo
      
    
         return
      

    def create_3dmesh_CATHY(self,gmsh_mesh=[],
                            NZONE=[],NSTR=[],N1=[],
                            NNOD=None, NTRI=None,
                            ZRATIO=[],Z1=[],
                            IVERT=0, ISP=0, BASE=4,
                            debug=False):
        """Create 3d mesh (grid file) from gmsh file.

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
        k = Project(typ='R2') # create a Project object in a working directory (can also set using k.setwd())
        k.importMesh(gmsh_mesh)
        k.showMesh()
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





    # -------------------------------------------------------------------#
    #%% Infitration DATA


    def create_infitration(self,dirfiles):
        self.set_drippers(dirfiles)


        pass


    def set_drippers(self,dirfiles,drip_pos='drippers.txt'):

        print(os.getcwd())
        if isinstance(drip_pos, str):
            self.drippers = np.loadtxt(os.path.join(self.project_name,dirfiles,drip_pos),skiprows=1,delimiter=',')
        else:
            self.drippers = drip_pos

        # check drippers position against DEM
        self.hapin
        mesh_x_max = float(self.hapin['xllcorner']) + float(self.hapin['delta_x'])*float(self.hapin['N'])
        mesh_y_max = float(self.hapin['yllcorner']) + float(self.hapin['delta_y'])*float(self.hapin['M'])

        mesh_x=[]
        for xx in range(int(self.hapin['N'])):
            mesh_x.append(float(self.hapin['xllcorner']) + float(self.hapin['delta_x'])*xx)

        mesh_y=[]
        for yy in range(int(self.hapin['M'])):
            mesh_y.append(float(self.hapin['yllcorner']) + float(self.hapin['delta_y'])*yy)


        print(mesh_x)
        print(mesh_y)



        print(mesh_x_max)
        print(max(self.drippers[:,0]))

        if mesh_x_max < max(self.drippers[:,0]):
            print('Error: max mesh_x=' + str(mesh_x_max) + '; max dripper x pos=' + str(max(self.drippers[:,0])))

        if mesh_y_max< max(self.drippers[:,1]):
            print('Error: max mesh_x=' + str(mesh_y_max) + '; max dripper y pos=' + str(max(self.drippers[:,1])))

       # for dd in self.drippers:
       #     dd==

       #  (A==B).all()




        self.drippers_nodes = []


        pass

    # -------------------------------------------------------------------#
    #%% ERT DATA

    def set_elecs(self,filename='elecs.csv'):

        self.elecs = np.loadtxt(filename,skiprows=1,delimiter=',')



        pass

    # -------------------------------------------------------------------#
    #%% DATA ASSIMILATION

    def create_archie():
        pass

    def create_data_ass():

        pass

    def create_elec_nodes():

        pass
