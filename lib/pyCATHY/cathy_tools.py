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

#import meshtools as mt

from pyCATHY import plot_tools as pltCT
import plot_tools as pltCT

# from pyCATHY import cathy_utils as utilsCT



class CATHY(object):
    '''Main CATHY object.'''
    
    def __init__(self,dirName,prjName='my_cathy_prj', notebook=False, **kwargs):
        '''Create CATHY object.
        '''
        
        print('init CATHY object')
        self.notebook = notebook # flag if the script is run in a notebook

        self.workdir = os.path.join(os.getcwd() , dirName)
        
        if not os.path.exists(os.path.join(self.workdir)):
            os.makedirs(self.workdir,exist_ok=True)
        os.chdir(self.workdir)

        self.project_name = prjName

        if not os.path.exists(os.path.join(self.workdir,self.project_name)):
            os.makedirs(os.path.join(self.workdir,self.project_name),exist_ok=True)
            
            
        self.processor_name = 'cathy'
        self.input_dirname = 'input'
        self.output_dirname = 'output'


        # inputs
        self.parm = {} # dict of parm input parameters
        self.soil = {} # dict of soil input parameters
        self.ic = {} # dict of ic input parameters

        
        # self.time = []
        # infitration
        self.drippers = []

        # ERT
        self.elecs = []

        # meshing
        self.mesh_gmsh = [] # gmsh meshfile



        for key,value in kwargs.items():
            
            if key == 'clear_src': # clear src files 
                if value == True:
                    if os.path.exists(os.path.join(self.workdir,self.project_name,'src')):
                        print('clear src files')
                        shutil.rmtree(os.path.join(self.workdir,self.project_name,'src'))
                        shutil.rmtree(os.path.join(self.workdir,self.project_name,'input'))
                        shutil.rmtree(os.path.join(self.workdir,self.project_name,'tmp_src'))
                        shutil.rmtree(os.path.join(self.workdir,self.project_name,'prepro'))
                
        
        # fetch src files if not existing from Gitbucket repo
        if not os.path.exists(os.path.join(self.project_name,'src')):
            print('src files not found')
            print(self.workdir)
            try:
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
            except:
                pass

        for key,value in kwargs.items():
            if key == 'clear_outputs':
                if value == True:
                        if not os.path.exists(os.path.join(self.workdir,self.project_name,'output')):
                            self.create_output(output_dirname='output')
                        else:
                            shutil.rmtree(os.path.join(self.workdir,self.project_name,'output'))
                            shutil.rmtree(os.path.join(self.workdir,self.project_name,'vtk'))
                            self.create_output(output_dirname='output')


        pass


    def run_preprocessor(self,KeepOutlet=True,verbose=False,**kwargs):
        """Run cppp.exe

        1. updates cathy parameters prepo
        2. recompile de source files prepo
        3. run the preprocessor using bash cmd
        
        Returns
        -------
        type
            Running the executable file has allowed to generate a complete
            set of files describing physiographic features of the drainage system, as shown in Table 2.

        """

        for loopi in range(2):                 # run it twice (to avoid the first error)
            os.chdir(os.path.join(self.workdir,self.project_name, 'prepro/src/'))
    
            #clean all files compiled
            for file in glob.glob("*.o"):
                os.remove(file)
    
            if self.notebook==False:
                bashCommand = 'gfortran -O -o pycppp mpar.f90 mbbio.f90 wbb_sr.f90 csort.f90 qsort.f90 depit.f90 cca.f90 smean.f90 dsf.f90 facet.f90 hg.f90 mrbb_sr.f90 bb2shp_sr.f90 shape.f90 dbase.f90 streamer.f90 cppp.f90'
                
                try:
                    # bashCommand = 'gcc -O -o pycppp mpar.f90 mbbio.f90 wbb_sr.f90 csort.f90 qsort.f90 depit.f90 cca.f90 smean.f90 dsf.f90 facet.f90 hg.f90 mrbb_sr.f90 bb2shp_sr.f90 shape.f90 dbase.f90 streamer.f90 cppp.f90'
                    
                    p = os.system(bashCommand)
                    p = os.system(bashCommand) # run it twice (to avoid the first error)

                    # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                    # process.communicate()
                    # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                    # process.communicate()

                    # if verbose:
                    #     output, error = process.communicate()
            
                    # process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            
                    # if verbose:
                    #     output, error = process.communicate()

                        
                except:
                    print('bash cmd not recognized')
                    


        
        
                os.chdir(self.workdir)
    
        #try:
        # move to the directory where the source FORTRAN files are contained (cathy_main.f)
        shutil.move(os.path.join(self.workdir,self.project_name, 'prepro/src/pycppp'),
                    os.path.join(self.workdir,self.project_name, 'prepro/pycppp'))

    
        print('run preprocessor')
        os.chdir(os.path.join(self.workdir, self.project_name, 'prepro'))

        bashcmd = './pycppp'             
        my_data = "2\n0\n1\n" # user input data
        p = subprocess.run([bashcmd], text=True, input=my_data, capture_output=True)
     
        # try:
        #     p = subprocess.run([bashcmd], text=True, input=my_data, capture_output=True)
        # except:
        #     shutil.move(os.path.join(self.workdir,self.project_name, 'prepro/src/cppp'),
        #     os.path.join(self.workdir,self.project_name, 'prepro/cppp'))
        #     p = subprocess.run([bashcmd], text=True, input=my_data, capture_output=True)

        if verbose:
            print(p.stdout)
            print(p.stderr)
        
        if  "catchment with more than one outlet cell!" in p.stdout:
            print("catchment with more than one outlet cell!")
            sys.exit()


        if  KeepOutlet==False:
            print('remove outlet')
            self.DEM[0]=self.DEM[1]
            with open(os.path.join(self.workdir,self.project_name,'prepro/dtm_13.val'), 'w+') as f:
                    np.savetxt(f, self.DEM , fmt='%1.4e')   # use exponential

        os.chdir(self.workdir)


        # -------------------------------------------- #
        # run processor only to build the 3d mesh
        # self.run_processor(verbose=True,IPRT1=2,TRAFLAG=0)
        # self.read_grid3d()
        # -------------------------------------------- #
 
        return

    def recompileSrc(self,verbose=True):

        os.chdir(os.path.join(self.workdir,self.project_name,'src'))

        print('recompileSrc')

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
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

        if verbose:
            output, error = process.communicate()

        # move to the directory where the source FORTRAN files are contained (cathy_main.f)
        # and OVERWRITE IF ALREADY EXISTING (full path spacified)
        
        os.chdir(os.path.join(self.workdir))
        shutil.move(os.path.join(self.project_name, 'src',self.processor_name),
                    os.path.join(self.project_name, self.processor_name))

        print('move ' + self.processor_name + ' into ' + self.project_name + ' directory')

        os.chdir(os.path.join(self.project_name))
        pass
    
        
    def run_processor(self,recompile=True,verbose=False,runProcess=True,**kwargs):
        """ Run cathy.exe
        
        1. updates cathy parameters
        2. recompile de source files (set False for notebook)
        3. run the processor using bash cmd

        Returns
        -------
        type
            Description of returned object.

        """

         
        import time

        t0 = time.time()


   
        if len(kwargs.items())>0:
            
            if verbose == True:
                print(kwargs)
        
            self.update_parm(**kwargs)
            self.update_cathyH(**kwargs)
            
       
        # 'CATHY.H'
        # --------
        # if isinstance(CATHY,'hapin'):
        #     self.update_cathyH()

        # check location of the exe (should be in project_name folder)
        os.chdir(os.path.join(self.workdir,self.project_name,'src'))

        if recompile==True:
            
            # self.recompileSrc()
            print('recompile src files')
            
            # gfortran -c *.f 
            # gfortran *.o -L\MinGW\lib -llapack -lblas -o cathy


            #clean all files compiled
            for file in glob.glob("*.o"):
                os.remove(file)
    
            for file in glob.glob("*.f"):
                bashCommand = "gfortran -c " + str(file)
    
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()
    
            files = ""
            for file in glob.glob("*.o"):
                files += " " + str(file)
    
            bashCommand = "gfortran" + files + " -llapack -lblas -o " + self.processor_name
            # print(bashCommand)
    
    
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE) 
            output, error = process.communicate()
    
            # move to the directory where the source FORTRAN files are contained (cathy_main.f)
            # and OVERWRITE IF ALREADY EXISTING (full path spacified)
            os.chdir(os.path.join(self.workdir))
            
            try:
                shutil.move(os.path.join(self.workdir, self.project_name, 'src', self.processor_name),
                            os.path.join(self.workdir, self.project_name, self.processor_name))
        
                print('move ' + self.processor_name + ' into ' + self.project_name + ' directory')
            except:
                print('cannot find the new processsor:' + str(self.processor_name))
    
            


        # 'cathy.fnames'
        # -------------
        # Some of these inputs are automatically generated during the pre-processing (Table 2),
        # while others are new input files that should be updated with appropriate parameters for the specific case study
        # The unit number associated to each input and output file (e.g. unit IIN1) are used by CATHY in reading and writing statements
        # use the backslash (\) in quotes in Windows, the forward slash (/) are used in Mac OS X.
        #'cathy_main.f'
        # -------------


        # fnames_file = open(os.path.join(self.project_name , 'cathy.fnames'), 'r')
        # Lines = fnames_file.readlines()
        
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
        
        os.chdir(os.path.join(self.workdir,self.project_name))

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

        if runProcess==True:

            print('run processor')
            # print(os.getcwd())
    
            callexe = './' + self.processor_name
            # process = subprocess.Popen(callexe.split(), stdout=subprocess.PIPE)
            # if verbose:
            #     output, error = process.communicate()
            p = subprocess.run([callexe], text=True, capture_output=True)
            
            if verbose == True:
                print(p.stdout)
                print(p.stderr)
    
            os.chdir(os.path.join(self.workdir))
    
            #try:
            #    subprocess.call([callexe_path, '-ARG'], shell=True)
            #    subprocess.Popen(["cathyEnv/bin/python"])
            #except:
            #    pass
        
            # print('timer')

        t1 = time.time()
        self.total = t1-t0
        

        
        return

    def time_run(self):
        
        print(self.total)
        print(self.cathyH['MAXCEL'])
        print(self.cathyH['MAXZON'])
        print(self.cathyH['MAXSTR'])
        # print(self.cathyH['TMAX'])
        
        pass



    def update_cathyH(self,verbose=False,**kwargs):
        """Short summary.

        Parameters
        ----------
        **kwargs : type
            ROWMAX : type
                'ROWMAX': self.hapin['N'], # maximum NROW, with NROW = number of rows in the DEM.

        Returns
        -------
        type
            Description of returned object.

        """
        
        if verbose == True:
            print('update_cathyH')
            print('-.-'*10)
        CATHYH_file = open(os.path.join(self.workdir,self.project_name ,'src' ,'CATHY.H'), 'r')
        Lines0_109 = CATHYH_file.readlines()
        CATHYH_file.close()
        
        DEMRES = 1
    
        self.cathyH = {
                        'ROWMAX': self.hapin['M'], # maximum NROW, with NROW = number of rows in the DEM
                        'COLMAX': self.hapin['N'], # maximum NCOL, with NCOL = number of columns in the DEM
                        # 'COARSE': ,
                        'MAXCEL': int(self.hapin['N'])*int(self.hapin['M']),
                        'MAXRES': 1,
                        'DEMRES': DEMRES,
                        'NODMAX': (int(self.hapin['N'])/DEMRES+1)*(int(self.hapin['M'])/DEMRES+1),
                        'NTRMAX': 2*(int(self.hapin['N'])*int(self.hapin['M']))/(DEMRES*DEMRES),
                        'NP2MAX': 1,
                        'MAXSTR': self.dem_parameters['nstr'],
                        'NPMAX': 1,
                        'NQMAX': 1,
                        'NSFMAX': 1,
                        'NNSFMX': 1,
                        'MAXNUDN': 1,
                        'MAXNUDT': 1,
                        'MAXNUDC': 1,
                        # 'MAXNENS': ,
                        'MAXZON': self.dem_parameters['nzone'], #maximum NZONE, with NZONE = number of material types in the porous medium
                        'MAXTRM': 52111,
                        'MAXIT': 30,
                        'NRMAX': self.parm['NR'],           
                        'MAXPRT': self.parm['NPRT'],#maximum NPRT (ref. parm file), with NPRT = number of time values for detailed output
                        'MAXVP': int(self.parm['NUMVP']), #maximum NUMVP (ref. parm input file), NUMVP = number of surface nodes for vertical profile output
                        'N1MAX' : self.dem_parameters['n1'],    #maximum N1 (it is good to have N1 ≤ 20), N1 = maximum number of element connections to a node     
                        'MAXBOT': 1,
                        'INTBOT': 1,
                        'NTAB': 100,
                        'MAXQOUT': 1,
                        'MAXENNOUT': 52560,   #Is related to data assimilation (DA) (to do not considered if you do not have DA)          
                        'MAXVEG': 1       

                    }       
        # self.run_processor(verbose=True,IPRT1=3,TRAFLAG=0)
        # self.read_grid3d()
        # # print('nnod' + str(self.nnod))
        # # print('NODMAX' + str(self.cathyH['NODMAX']))
        # print('N1MAX' + str(self.cathyH['N1MAX']))
        # print('MAXVP' + str(self.cathyH['MAXVP']))
        # sys.exit()
        
        
        # create dictionnary from kwargs
        for keykwargs,value in kwargs.items():
            if verbose == True:
                print(f'modified: {keykwargs} | value: {value}')
            self.cathyH[keykwargs]=value
           
        with open(os.path.join(self.workdir,self.project_name ,'src' ,'CATHY.H'), 'w+') as CATHYH_file:
             for i, l in enumerate(Lines0_109):
                 if i < 109:
                     CATHYH_file.write(l)

             CATHYH_file.write('      PARAMETER (ROWMAX={},COLMAX={},DEMRES={})\n'.format(self.cathyH['ROWMAX'], self.cathyH['COLMAX'], self.cathyH['DEMRES']))
             CATHYH_file.write('      PARAMETER (MAXCEL=ROWMAX*COLMAX,MAXRES=1)\n')
             CATHYH_file.write('      PARAMETER (NODMAX=(ROWMAX/DEMRES+1)*(COLMAX/DEMRES+1))\n')
             CATHYH_file.write('      PARAMETER (NTRMAX=2*MAXCEL/(DEMRES*DEMRES))\n'.format())
             CATHYH_file.write('      PARAMETER (NP2MAX=1,MAXSTR={})\n'.format(self.cathyH['MAXSTR']))
             CATHYH_file.write('      PARAMETER (NFACEMAX=74000)\n'.format())
             CATHYH_file.write('      PARAMETER (NMAX=NODMAX*(MAXSTR + 1),NTEMAX=3*NTRMAX*MAXSTR)\n'.format())
             CATHYH_file.write('      PARAMETER (NPMAX={},NPMAX_TRA=1,NQMAX={},NSFMAX={})\n'.format(self.cathyH['NPMAX'],
                                                                                                    self.cathyH['NQMAX'],
                                                                                                    self.cathyH['NSFMAX']))
             CATHYH_file.write('      PARAMETER (NNSFMX={},MAXDIR=NODMAX+NPMAX+NSFMAX*NNSFMX)\n'.format(self.cathyH['NNSFMX']))
             CATHYH_file.write('      PARAMETER (MAXNUDN={},MAXNUDT={},MAXNUDC={})\n'.format(self.cathyH['MAXNUDN'], self.cathyH['MAXNUDT'], 
                                                                                             self.cathyH['MAXNUDC']))
             CATHYH_file.write('      PARAMETER (MAXZON={},MAXTRM={},MAXIT={},MAXVEG={})\n'.format(self.cathyH['MAXZON'],self.cathyH['MAXTRM'],
                                                                                                   self.cathyH['MAXIT'],self.cathyH['MAXVEG']))
             CATHYH_file.write('      PARAMETER (NRMAX={},MAXPRT={},MAXVP={})\n'.format(self.cathyH['NRMAX'],
                                                                                        self.cathyH['MAXPRT'],
                                                                                        self.cathyH['MAXVP']))
             CATHYH_file.write('      PARAMETER (N1MAX={},NTPMAX=N1MAX*NMAX)\n'.format(self.cathyH['N1MAX']))
             CATHYH_file.write('      PARAMETER (MAXBOT={},INTBOT={},MAXQOUT={})\n'.format(self.cathyH['MAXBOT'],
                                                                                           self.cathyH['INTBOT'],
                                                                                           self.cathyH['MAXQOUT']))
             CATHYH_file.write('cxcx  PARAMETER (NIAUXMAX=NFACEMAX + MAXTRM + 1)\n'.format())
             CATHYH_file.write('cxcx  PARAMETER (NRAUXMAX=5*NFACEMAX + MAXTRM,NQMAX_TRA=NODMAX)\n'.format())
             CATHYH_file.write('      PARAMETER (NIAUXMAX=NMAX + MAXTRM + 1)\n'.format())
             CATHYH_file.write('      PARAMETER (NRAUXMAX=5*NMAX + MAXTRM,NQMAX_TRA=NODMAX)\n'.format())
             CATHYH_file.write('      PARAMETER (MAXVTKPRT=9)\n'.format())
             CATHYH_file.write('      PARAMETER (MAXFCONTONODE=100,MAXLKP=3)\n')
      
        CATHYH_file.close()
        pass

            
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

    def update_prepo_inputs(self,DEM=None,verbose=False,show=False,**kwargs):
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

        terrain_parameter = ['pt','imethod','lambda','CC_threshold','ndcf','nchc','A_threshold',
                                'ASk_threshold','kas','DN_threshold','local_slope_t',
                                'p_outflow_vo','bcc','cqm','cqg']
        
        # 'RIVULET NETWORK PARAMETERS (HYDRAULIC GEOMETRY OF THE SINGLE RIVULET)'
        rivulet_parameter = ['dr','As_rf','(Qsf_rf,w_rf)','(Wsf_rf,b1_rf,b2_rf)',
                               '(kSsf_rf,y1_rf,y2_rf)','Qsi_rf']
        

        channel_parameter = ['As_cf','(Qsf_cf,w_cf)','(Wsf_cf,b1_cf,b2_cf)','(kSsf_cf,y1_cf,y2_cf)',
                               'Qsi_cf']

        # idhapin = np.ones(22) + [1,1,2,3,3,1] #+ [1,2,3,3,1]
        hapin = structural_parameter + terrain_parameter + rivulet_parameter + channel_parameter
        
        hap_file = open(os.path.join(self.workdir, self.project_name , 'prepro/hap.in'), 'r')
        Lines = hap_file.readlines()
        hap_file.close()

        # read current values in hap.in


        
        tmp_param_value = []
        tmp_lnb = [] # save line nb where values are existing
        count=0
        for line in Lines:
            x = line.strip().split("=")
            if len(x)>1: # skip lines with no =
                xs= " ".join(x[1].split())
                xsr = xs.replace(" ", ",")
                l = xsr.split(',')
                if isinstance(l,list):
                    l = "/t".join(l)
                tmp_param_value.append(l)
                tmp_lnb.append(count)
                count += 1


        self.hapin = {}

        for i in range(len(tmp_param_value)):
            self.hapin[hapin[i]] = tmp_param_value[i]
            
        # print('--'*30)
        # print(self.hapin)

        
        if self.hapin['dr'] != self.hapin['delta_x']:
            print('adapt rivulet param to dem resolution')
            self.hapin['dr'] = self.hapin['delta_x']

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

        L=[]
        tmp_param_value_new=[]
        for key,value in kwargs.items():
            self.hapin[key] = value

        # # print('--'*30)
        # print(self.hapin)

        count=0
        for i, line in enumerate(Lines):
            xnew = line.strip().split("=")
            if len(xnew)>1: # skip lines with no =
                for key,value in self.hapin.items():
                    if count<len(hapin): # only consider structural parameters
                        # print(count)
                        if key==hapin[tmp_lnb[count]]:
                            # print(f'key: {key} | value: {value}')
                            if value!=tmp_param_value[count]:
                                if isinstance(value,list):
                                    value_str = "/t".join(value)
                                    xnew[1]=value_str
                                else:
                                    xnew[1]=value
                                line= xnew[0] + '=              ' + str(xnew[1])  + '\n'
                                # print(line)
                        tmp_param_value_new.append(value)

                count += 1 # count line nb
            L.append(line)

        # # write the new hap.in file

        hap_file = open(os.path.join(self.workdir, self.project_name , 'prepro/hap.in'), 'w+')
        hap_file.writelines(L)
        hap_file.close()

        # print(self.hapin)

        # self.hapin = {}

        # for i in range(len(structural_parameter)):
        #     key = structural_parameter[i]
        #     self.hapin[key] = tmp_param_value_new[i]

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
                  self.DEM = DEM_withOutlet

            print('update dtm_13.val')
            with open(os.path.join(self.project_name,'prepro/dtm_13.val'), 'w+') as f:
                np.savetxt(f, self.DEM, fmt='%1.4e')   # use exponential notation


        self.update_dem_parameters(**kwargs)
        # self.update_transp(**kwargs)
        
        if show==True:
            print('show')
            # pltCT.test(self.workdir, self.project_name)
            
        pass

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
        if hasattr(self, 'dem_parameters') == False:
            
            print('cannot finc dem paramters')

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
           print(keykwargs,value)
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


        for key,value in self.dem_parameters.items():
            if isinstance(value, list):
                strlst = "\t".join(str(e) for e in value)
                print(strlst)
                self.dem_parameters[key]=strlst
                

        # write file
        header_fmt = [1,1,1,1,3,3,1]
        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'dem_parameters'), 'w+') as dem_parametersfile:
            
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



        dem_parametersfile.close()
        
        pass


    #%% Ouput/INPUT FILES


    def update_zone(self,zone_xyz=[]):
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
        with open(os.path.join(self.workdir , self.project_name, 'prepro/zone'), 'w+') as zonefile:

            zonefile.write('north:     0' + "\n")         
            zonefile.write('south:     0' + "\n")         
            zonefile.write('east:     0' + "\n")         
            zonefile.write('west:     0' + "\n")         
            zonefile.write('rows:     ' + str(self.hapin['M']) + "\n")         
            zonefile.write('cols:     ' + str(self.hapin['N']) + "\n") 
            if len(zone_xyz)==0:
                zone_xyz = np.c_[np.ones([self.hapin['M'],self.hapin['N']])]
                np.savetxt(zonefile,zone_xyz,fmt='%i')            
            else:
                # if np.shape(zone_xyz)== :
                np.savetxt(zonefile,zone_xyz,fmt='%i')            

        zonefile.close()
        
        # update number of zone in the dem parameter file
        self.update_dem_parameters(nzone=len(np.unique(zone_xyz)))
        self.update_parm()
        self.update_cathyH(MAXZON=len(np.unique(zone_xyz)))


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


    def update_parm(self,verbose=False,**kwargs):


        # set default parameters
        self.parm =	{
                      "IPRT1": 2,
                      "NCOUT": 0,
                      "TRAFLAG": 1,
                      "ISIMGR": 2, #Flag for type of simulation and type of surface grid
                      "PONDH_MIN": 0.00, #Minimum ponding head
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
                      "DTMIN": .00001, #Minimum FLOW3D time step size allowed
                      "DTMAX": 10., #Maximum FLOW3D time step size allowed
                      "TMAX":  3600.0, #Time at end of simulation (TMAX is set to 0.0 for steady state problem)
                      "DTMAGA":  0.0,
                      "DTMAGM": 1.1,
                      "DTREDS": 0.0,
                      "DTREDM": .5,
                      "IPRT": 4,
                      "VTKF":7,
                      "NPRT": 3,
                      "(TIMPRT(I),I=1,NPRT)": [1800.,3600.,7200.],
                      "NUMVP": 1,
                      "(NODVP(I),I=1,NUMVP)": [441],
                      "NR": 0,
                      "NUM_QOUT": 0,
                      "(ID_QOUT(I),I=1,NUM_QOUT)": [441]
                      }


        #%%
        
        # DAFLAG:
        # Flag for the choice of the data assimilation scheme:
        # = 0 nudging (if NUDN=0, no data assimilation)
        # = 1 EnKF with Evensen's algorithm (Ocean Dynamics, 2004)
        # = 2 EnKF with parameters update
        # = 3 Particle filtering (SIR algorithm)
        # = 4 Particle filtering (SIR algorithm) with parameters update



        # create dictionnary from kwargs

        for keykwargs,value in kwargs.items():
            if verbose == True:
                print(f'keykwargs: {keykwargs} | value: {value}')
            if keykwargs=='TIMPRTi':
                key = '(TIMPRT(I),I=1,NPRT)'
                self.parm[key]=value
                if len(value) != self.parm['NPRT']:
                    self.parm['NPRT']=len(value)

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
        
        self.update_cathyH(MAXPRT=self.parm['NPRT'])

        pass





    def update_ic(self,INDP=2,IPOND=0,WTPOSITION=[]):
        """Short summary.

        Parameters
        ----------
        INDP : int
            Flag for pressure head initial conditions (all nodes)
            - =0 for input of uniform initial conditions (one value read in)
            - =1 for input of non-uniform IC's (one value read in for each node)
            - =2 for calculation of fully saturated vertical hydrostatic equilibrium IC's (calculated in subroutine ICVHE). In the case of IPOND>0, the fully saturated hydrostatic IC is calculated (in subroutine ICVHEPOND) starting from the ponding head values at the surface nodes, rather than surface pressure heads of 0.
            - =3 for calculation of partially saturated vertical hydrostatic equilibrium IC's (calculated in subroutine ICVHWT) with the water table height (relative to the base of the 3‐d grid) given by parameter WTHEIGHT 
            - =4 for calculation of partially saturated vertical hydrostatic equilibrium IC's (calculated in subroutine ICVDWT) with the water table depth (relative to the surface of the 3‐d grid) given by parameter WTHEIGHT 
        WTPOSITION : type
            For the case INDP=3, specifies the initial water table height relative to the base of the 3‐d grid
        IPOND : type
            Flag for ponding head initial conditions (surface nodes)
            - =0 no input of ponding head initial conditions; otherwise (IPOND = 1 or 2) ponding head initial conditions are read into PONDNOD, and, where PONDNOD > 0, these values are used to update the surface node values in PTIMEP read in according to the previous INDP flag
            - =1 uniform ponding head initial conditions (one value read in)
            - =2 non-uniform ponding head initial conditions (one value read in for each node)

        Returns
        -------
        type
            Description of returned object.

        """

        # check value of WTPOSITION
        # if WTPOSITION>0:
        #     print('WTPOSITION must be negative - water table height relative to the base of the 3‐d grid')
        #     sys.exit() 
            
            
        # set default parameters
        self.ic =	{
                      "INDP": INDP,
                      "WTPOSITION": WTPOSITION, #For the case INDP=3, specifies the initial water table 
                                                  #height relative to the base of the 3‐d grid
                      "IPOND": IPOND
                      }
        
        
        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'ic'), 'w+') as icfile:
            icfile.write(str(INDP) + "\t" + str(IPOND) + "\t"
                            + 'INDP' + "\t" + 'IPOND' + "\n")
            icfile.write(str(WTPOSITION) + "\t" +  'WTPOSITION' + "\n")                         
        icfile.close()

        pass



    def update_atmbc(self, HSPATM=0,IETO=0,TIME=None,VALUE=[None,None],show=False,verbose=False):
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
        vdiff = VALUE[0]-VALUE[1]

        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'atmbc'), 'w+') as atmbcfile:
            atmbcfile.write(str(HSPATM) + "\t" + str(IETO) + "\t"
                            + 'HSPATM' + "\t" + 'IETO' + "\n")
            
            for t,v in zip(TIME,vdiff):
                
                if verbose == True:
                    print(t,v)
                atmbcfile.write(str(t) + 
                                  "\t" + 'TIME' + "\n")
                atmbcfile.write(str(v) + 
                                  "\t" + 'VALUE' + "\n")
                 
        atmbcfile.close()
               
        
        self.update_parm(NPRT=len(TIME))
        self.update_cathyH(MAXPRT=len(TIME))
        
        if show == True:
            # if HSPATM !=0:
            #     print('impossible to plot for non homogeneous atmbc')
            #     # sys.exit()
            # else:
            pltCT.atmbc_inputs_plot(TIME,VALUE)
            
        
        pass
    

    def update_nansfdirbc(self,TIME=[], NDIR=0,NDIRC=0, NQ3=None, noflow=True, bound_xyz=None, pressureHead=[]):
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
         
         !!! update_nansfdirbc use the grid3d 

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
        
        try:
            self.read_grid3d()
        except OSError:
            print('grid3d missing - need to run the processor with IPRT1=2 first')
            # self.run_processor(verbose=True,IPRT1=2,TRAFLAG=0)
            # self.read_grid3d()
            # # sys.exit()            

     
        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'nansfdirbc'), 'w+') as nansfdirbcfile:

            if noflow:    
                if len(TIME)==0:
                    TIME = self.atmbc['TIME']
                for tt in TIME:
                    nansfdirbcfile.write(str(tt) + "\t" + 'TIME' + "\n")
                    nansfdirbcfile.write(str(NDIR) + "\t" + str(NDIRC) + "\t"
                                    + 'NDIR' + "\t" + 'NDIRC' + "\n")
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

    def update_nansfneubc(self,TIME=[], NQ=0,ZERO=0):

        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'nansfneubc'), 'w+') as nansfneubcfile:

            if len(TIME)==0:
                TIME = self.atmbc['TIME']
            for tt in TIME:
                nansfneubcfile.write(str(tt) + "\t" + 'TIME' + "\n")
                nansfneubcfile.write(str(ZERO) + "\t" + str(NQ) + "\t"
                                + 'ZERO' + "\t" + 'NQ' + "\n")
            
        nansfneubcfile.close()

        # set default parameters
        # self.nansfneubc =	{
        #               "NDIR": HSPATM,
        #               "NDIRC": IETO,
        #              }
        
        pass

    def update_sfbc(self,TIME=[]):

        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'sfbc'), 'w+') as sfbcfile:

            if len(TIME)==0:
                TIME = self.atmbc['TIME']
            for tt in TIME:
                sfbcfile.write(str(tt) + "\n")
                sfbcfile.write('0' + "\n")
            
        sfbcfile.close()

        # set default parameters
        # self.nansfneubc =	{
        #               "NDIR": HSPATM,
        #               "NDIRC": IETO,
        #              }
        
        pass    

    def update_soil(self, FP= [], SPP=[], verbose=False,                   
                    **kwargs):
        """Soil parameters (soil - IIN4).
        The porous media properties are defined in the soil file.
        The first thing that must be decides is the type of relationship
        to describe the hydraulic characteristics of the unsaturated soil (i.e. retention curves).
        This can be done through the choice of the parameter IVGHU amongst the several options.

        Parameters
        ----------
           
        FP: Feddes Parameters [PCANA PCREF PCWLT ZROOT PZ OMGC] : [float float float float float float]
            'PCANA': float
                anaerobiosis point
            'PCREF': float
                field capacity
            'PCWLT': float
                wilting point
            'ZROOT': float
                root depth
            'PZ': float
                ??
            'OMGC': float
                ??
            For details, see http://dx.doi.org/10.1002/2015WR017139

                  
        SPP : Soil Physical Properties
            'PERMX' (NSTR, NZONE): saturated hydraulic conductivity - xx
            'PERMY' (NSTR, NZONE): saturated hydraulic conductivity - yy
            'PERMZ' (NSTR, NZONE): saturated hydraulic conductivity - zz
            'ELSTOR' (NSTR, NZONE): specific storage
            'POROS' (NSTR, NZONE): porosity (moisture content at saturation)
            
            Parameters for van Genuchten and extended van Genuchten moisture curves
            'VGNCELL': 
            'VGRMCCELL': residual moisture content
            'VGPSATCELL': saturated water content
                
        **kwargs : type
            PMIN : int
                air dry' pressure head value (for switching control of atmospheric boundary conditions during evaporation)
                [m sec]
            IPEAT : int
                Flag for peat soil deformation
                =0 constant porosity (in unsaturated soil)
                =1 consider porosity variations with water saturation
            SCF : int
                soil cover fraction (fraction of soil covered in vegetation)
                
            
        Returns
        -------
        type
            write the soil file.

        """
        
        if len(SPP)==0:#set defaults parameters
            PERMX = PERMY = PERMZ =1.88E-04
            ELSTOR = 1.00E-05
            POROS = 0.55
            VGNCELL = 1.46
            VGRMCCELL = 0.15
            VGPSATCELL =  0.03125
                
            SPP = {'PERMX':PERMX,'PERMY':PERMY,'PERMZ':PERMZ,
                    'ELSTOR':ELSTOR,'POROS':POROS,
                    'VGNCELL':VGNCELL,'VGRMCCELL':VGRMCCELL,'VGPSATCELL':VGPSATCELL}

        # set default parameters         
        self.soil =	{
                      "PMIN":-5.0,
                      
                      "IPEAT": 0,
                      "SCF":1.0,
                      
                      "CBETA0": 0.4,
                      "CANG": 0.225,
                      
                      # Feddes parameters default values
                       "PCANA":[0.0],
                       "PCREF":[-4.0],
                       "PCWLT":[-150],
                       "ZROOT":[1.0],
                       "PZ":[1.0],
                       "OMGC":[1.0],
            
                      
                      "IVGHU": 0,
                      
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
                      "BCPSAT": -0.345
                }

        # if len(FP)==0:#set defaults parameters
           
        #     FP = {'PCANA':self.soil['PCANA'],'PCREF':self.soil['PCREF'],'PCWLT':self.soil['PCWLT'],
        #                       'ZROOT':self.soil['ZROOT'],'PZ':self.soil['PZ'],
        #                       'OMGC':self.soil['OMGC']}
            
    

        # read kwargs
        for keykwargs,value in kwargs.items():
            if verbose == True:
                print(f'key kwargs: {keykwargs} | value: {value}')
            self.soil[keykwargs]=value
            self.parm[keykwargs]=value
        
        for fp in FP: # loop over fedded parameters
            if verbose == True:
                print(fp,FP[fp])
            self.soil[fp]=FP[fp]


        
        if hasattr(self,'dem_parameters') is False:
            self.update_prepo_inputs()
        
            
        #%% Soil physical properties       
        # Physical prop strat by strat
        if self.dem_parameters['nzone']>1:
            SoilPhysProp = np.ones([self.dem_parameters['nzone']*self.dem_parameters['nstr'],8])
            k=0
            for istr in range(self.dem_parameters['nstr']): # loop over strates
                izoneSoil = np.zeros([self.dem_parameters['nzone'],8])
                for izone in range(self.dem_parameters['nzone']): # loop over zones within a strate
                    izoneSoil_tmp = []
                    for spp in SPP:
                        izoneSoil_tmp.append(SPP[spp][izone])
                    izoneSoil_tmp = np.hstack(izoneSoil_tmp)
                    izoneSoil[izone,:]= izoneSoil_tmp
                    ki = k 
                    ke = k + self.dem_parameters['nzone']
                    SoilPhysProp[ki:ke,:]=izoneSoil
                    # SoilPhysProp[self.dem_parameters['nzone']*istr*izone:self.dem_parameters['nzone']*istr*izone+self.dem_parameters['nzone'],:]=izoneSoil
                k+= self.dem_parameters['nzone']
        else:
            izoneSoil=[]
            for spp in SPP:
                izoneSoil.append(SPP[spp])
            # SoilPhysProp = np.ones([self.dem_parameters['nstr'],8])*izoneSoil
            izoneSoil = np.hstack(izoneSoil)
            SoilPhysProp = np.tile(izoneSoil,(self.dem_parameters['nstr'],1))

        #%%  Vegetation properties
        
        # PCANA,PCREF,PCWLT,ZROOT,PZ,OMGC
        # Read only if IVGHU=0
        print(self.cathyH['MAXVEG'])

        if self.cathyH['MAXVEG']>1:
            print(self.cathyH['MAXVEG'])
            FeddesParam = np.zeros([self.cathyH['MAXVEG'],6])
            for iveg in range(self.cathyH['MAXVEG']): # loop over veg zones within a strate
                izoneVeg_tmp = []
                for sfp in FP:
                    izoneVeg_tmp.append(FP[sfp][iveg])
                    # if iveg==1:
                    #     izoneVeg_tmp.append('\t PCANA PCREF PCWLT ZROOT PZ OMGC')
                    
                izoneVeg_tmp = np.hstack(izoneVeg_tmp)
                FeddesParam[iveg,:]= izoneVeg_tmp
        else:
            FeddesParam= np.c_[self.soil['PCANA'],self.soil['PCREF'],
                            self.soil['PCWLT'],self.soil['ZROOT'],
                            self.soil['PZ'],self.soil['OMGC']]
            



        #%% write soil file
        counth=0
        header_fmt_soil = [1,2,2,6,1,5,1,2,3]

        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'soil'), 'w+') as soilfile:
            for i, h in enumerate(header_fmt_soil):
                left=right=[]
                left = str(list(self.soil.values())[counth:counth+h])            
                left = left.strip('[]').replace(",", "")
                
                right = str(list(self.soil.keys())[counth:counth+h])
                right = right.strip('[]').replace(",", "")
                right = right.replace("'", "")
                
                if i==3: # Feddes parameters
                    np.savetxt(soilfile,FeddesParam,fmt='%1.3e')
                    counth += h
                else:    
                    line = left + '\t' + right + "\n"
                    counth += h
                    soilfile.write(str(line))
                    
            np.savetxt(soilfile,SoilPhysProp,fmt='%1.3e')
            soilfile.write('PERMX    PERMY    PERMZ    ELSTOR   POROS,VGNCELL,VGRMCCELL,VGPSATCELL' + "\n")
                
        soilfile.close()
        pass
    
    

    def update_root_map(self, root_depth=1, show=False):
        """Contains the raster map describing which type of vegetation every cell belongs to.

        Returns
        -------
        type
            Description of returned object.

        """
        if isinstance(root_depth,int):
            root_depth = float(root_depth)

        if hasattr(self,'hapin') is False:
            self.update_prepo_inputs()
            
            
        print('update root map')
        with open(os.path.join(self.workdir , self.project_name, self.input_dirname, 'root_map'), 'w+') as rootmapfile:
 
            rootmapfile.write('north:     0' + "\n")         
            rootmapfile.write('south:     0' + "\n")         
            rootmapfile.write('east:     0' + "\n")         
            rootmapfile.write('west:     0' + "\n")         
            rootmapfile.write('rows:     ' + str(self.hapin['M']) + "\n")         
            rootmapfile.write('cols:     ' + str(self.hapin['N']) + "\n")  
            
            if isinstance(root_depth,float):
                # if  root_depth>self.dem_parameters['base']:
                #     print('max root mesh > max mesh depth')
                #     sys.exit()
                root_depth = np.c_[np.ones([int(self.hapin['M']),
                                                      int(self.hapin['N'])])]*root_depth
                np.savetxt(rootmapfile,root_depth,
                           fmt='%1.2e')
            else:
                # if np.shape(zone_xyz)== :
                # if  max(max(root_depth))>self.dem_parameters['base']:
                    # print('max root mesh > max mesh depth')
                    # sys.exit()
                np.savetxt(rootmapfile,root_depth,fmt='%1.2e')    
                
        rootmapfile.close()
        
        self.update_cathyH(MAXVEG=len(np.unique(root_depth)))
        
        if show is not None:
            pltCT.rootMap_plot(root_depth)
        
        pass

    def plant():
        # plant parameters only exist for CATHYv Manoli


        pass

    #%% Meshtool functions

    def read_grid3d(self):
      
         
         print('reading grid3d') 
         grid3d_file = open(os.path.join(self.project_name , 'output/grid3d'), 'r')
         # Lines = grid3d_file.readlines()
         self.nnod,self.nnod3,self.nel = np.loadtxt(grid3d_file,max_rows=1)
         grid3d_file.close()

         grid3d_file = open(os.path.join(self.project_name , 'output/grid3d'), 'r')
         mesh_tetra = np.loadtxt(grid3d_file,skiprows=1,max_rows=int(self.nel)-1)
         grid3d_file.close()


         # mesh_tetra = np.loadtxt(grid3d_file,skiprows=1+int(self.nel)-1, max_rows=1+self.nel+self.nnod3-1)
         grid3d_file = open(os.path.join(self.project_name , 'output/grid3d'), 'r')
         mesh3d_nodes = np.loadtxt(grid3d_file,skiprows=1+int(self.nel), max_rows=1+int(self.nel)+int(self.nnod3)-1)
         grid3d_file.close()


         xyz_file = open(os.path.join(self.project_name , 'output/xyz'), 'r')
         nodes_idxyz = np.loadtxt(xyz_file,skiprows=1)
         xyz_file.close()
         
         
         # self.xmesh = mesh_tetra[:,0]
         # self.ymesh = mesh_tetra[:,1]
         # self.zmesh = mesh_tetra[:,2]     
         # return mesh3d_nodes
         
         self.grid = {
                     'nnod': self.nnod, # number of surface nodes
                     'nnod3': self.nnod3, # number of volume nodes
                     'nel': self.nel,
                     'mesh3d_nodes': mesh3d_nodes,
                     'mesh_tetra': mesh_tetra,                     
                     'nodes_idxyz': nodes_idxyz,                     
                     }

         return self.grid
      

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
