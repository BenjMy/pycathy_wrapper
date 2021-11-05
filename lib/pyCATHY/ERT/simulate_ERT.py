"""Class managing ERT data simulation and inversion + petrophysical transformation
"""


import matplotlib.pyplot as plt
import numpy as np
from resipy import R2 # geophysics tools
import pyvista as pv # if not installed : pip install with conda install pyvista (in a conda terminal)
import numpy as np
import os 


        
def create_ERT_survey(pathERT,elecsXYZ,sequence,mesh, **kwargs):
    
    #https://hkex.gitlab.io/resipy/api.html

    # os.chdir(pathERT) 
    
    isExist = os.path.exists(pathERT)

    if not isExist:
      
      # Create a new directory because it does not exist 
      os.makedirs(pathERT) 
    
    ERT = R2(pathERT + 'ERT_fwdmodel', typ='R3t')
    ERT.setTitle('Rhizo_synth')
    
    

    ERT.setElec(np.c_[elecsXYZ[:,0],elecsXYZ[:,2],elecsXYZ[:,1],elecsXYZ[:,3]])
    
    ERT.importMesh(mesh)
    
    if 'res0' in kwargs:
        ERT.setRefModel(kwargs['res0'])

    ERT.addRegion(np.array([[0.1,0.1],[0.1,0.4],[0.3,0.4],[0.3,0.1],[0.1,0.1]]), 500, iplot=False)
    ERT.importSequence(sequence)
    # -----------------------------------------------
    
    return ERT

def fwd_ERT_survey(ERT,noise,show=False):

    # ----------#
    # fwd modelling
    # ---------------------------------------------------------#

    ERT.forward(noise=0.05, iplot=show) # forward modelling with 5 % noise added to the output



def invert_ERT_survey(ERT,show=False):

    # ---------------------------------------------------------#
    # inversion
    # ---------------------------------------------------------#
    ERT.param['num_xy_poly'] = 0
    ERT.param['zmin'] = -np.inf
    ERT.param['zmax'] = np.inf
    ERT.param['data_type'] = 1 # using log of resistitivy
    ERT.err = False # if we want to use the error from the error models fitted before
    ERT.param['a_wgt'] = 0.001
    ERT.param['b_wgt'] = 0.05
    
    ERT.invert()
    
    if show==True:
        ERT.showResults(index=0) # show the initial model
        ERT.showResults(index=1, sens=False) # show the inverted model
    


def Archie_ERT(ERT,rFluid,porosity,m,n):
    '''
    Archie transformation
    (only possible on inverted values)

    Parameters
    ----------
    ERT : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    
    # (rho/rFluid/porosity^-m)^(-1/n)
    rFluid=1
    porosity=0.55
    m = 1
    n = 1
    
    str_formula = f"(x['Resistivity']* {rFluid} * {porosity} **(-m))**(-1/ {n})"
    # print(str_formula)
       
    
    # only possible on inverted mesh (meshResults)
    ERT.computeAttribute(str_formula, name='saturation')
    
    pass




# def ERT_init(ARCHIE,NERT):
#     '''
    

#     Parameters
#     ----------
#     ARCHIE : TYPE
#         DESCRIPTION.
#     NERT : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     None.

#     '''


#     OPEN(18,FILE='input/archie',status='old')
    
#     READ(18,*) ARCHIE(1)
#     READ(18,*) ARCHIE(2)
#     READ(18,*) ARCHIE(3)
#     READ(18,*) ARCHIE(4)
#     READ(18,*) NERT
#     CLOSE(18)

#     return



# def update_archie(A_AR, CW_AR, M_AR, N_AR):
#     '''
    

#     Parameters
#     ----------
#     A_AR : TYPE
#         DESCRIPTION.
#     CW_AR : TYPE
#         DESCRIPTION.
#     M_AR : TYPE
#         DESCRIPTION.
#     N_AR : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     A_AR : TYPE
#         DESCRIPTION.
#     CW_AR : TYPE
#         DESCRIPTION.
#     M_AR : TYPE
#         DESCRIPTION.
#     N_AR : TYPE
#         DESCRIPTION.

#     '''
    
#     A_AR=ARCHIE(1)
#     CW_AR=ARCHIE(2)
#     M_AR=ARCHIE(3)
#     N_AR=ARCHIE(4)
      
      
#     return A_AR, CW_AR, M_AR, N_AR



# def ERT_measure(NLKP,N,NT,IVGHU,TETRA,TP,
#                 ENPNEW,ENPNODI,ENSNODI,ARCHIE,EN_ERT,
#                 ENPT,NERT,NTRI,NENS,ELECTRC_NODI,
#                 PNEW,SNODI,PNODI,SW,CKRW,SWE,CKRWE):
#     '''
#     Reads pressure values, calculates saturation at the assimilation
#     time and calculates the electrical conductivity of the soil by the
#     Archie law = ELECTRC_NODI (for the 2D cross section only).
#     ELETRC_NODI is calculated for the expanded 2D grid -> input for
#     SAT2D

#     Parameters
#     ----------
#     NLKP : TYPE
#         # of lookup values in the moisture curve tables for
#         IVGHU = -1. NLKP values of PCAP, SATC, and KRWC will be
#         read in for each zone (inner input loop) and for each layer
#         (outer input loop). NLKP must be at least 3, and for each
#         set of NLKP values, the PCAP values must be in ascending
#         order (e.g., -10.0, -9.9, -9.8, ..., -0.04, -0.02, 0.0)..
#     N : TYPE
#         DESCRIPTION.
#     NT : TYPE
#         DESCRIPTION.
#     IVGHU : TYPE
#         =-1 for moisture curve lookup table
#         =0 for van Genuchten moisture curves
#         =1 for extended van Genuchten moisture curves
#         =2 for moisture curves from Huyakorn et al (WRR 20(8) 1984,
#         WRR 22(13) 1986) with Kr=Se**n conductivity relationship
#         =3 for moisture curves from Huyakorn et al (WRR 20(8) 1984,
#         WRR 22(13) 1986) with conductivity relationship from
#         Table 3 of 1984 paper (log_10 Kr(Se) curve).
#     TETRA : TYPE
#         DESCRIPTION.
#     TP : TYPE
#         DESCRIPTION.
#     ENPNEW : TYPE
#         DESCRIPTION.
#     ENPNODI : TYPE
#         DESCRIPTION.
#     ENSNODI : TYPE
#         DESCRIPTION.
#     ARCHIE : TYPE
#         - parameters for the application of Archie low:
#         C_el= 1/a C_w PNODI^m SW^n
#         - ARCHIE(1) = a
#         - ARCHIE(2) = C_w
#         - ARCHIE(3) = m
#         - ARCHIE(4) = n
#     EN_ERT : TYPE
#         ERT measures associated with the ensemble.
#     ENPT : TYPE
#         DESCRIPTION.
#     NERT : TYPE
#         number of ERT observations (measurements) in the EnKF streamflow or ERT observation values .
#     NTRI : TYPE
#         DESCRIPTION.
#     NENS : TYPE
#         DESCRIPTION.
#     ELECTRC_NODI : TYPE
#         DESCRIPTION.
#     PNEW : TYPE
#         DESCRIPTION.
#     SNODI : TYPE
#         DESCRIPTION.
#     PNODI : TYPE
#         DESCRIPTION.
#     SW : TYPE
#         DESCRIPTION.
#     CKRW : TYPE
#         DESCRIPTION.
#     SWE : TYPE
#         DESCRIPTION.
#     CKRWE : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     None.

#     '''


#     OPEN(18,FILE='input/archie',status='old')
    
#     READ(18,*) ARCHIE(1)
#     READ(18,*) ARCHIE(2)
#     READ(18,*) ARCHIE(3)
#     READ(18,*) ARCHIE(4)
#     READ(18,*) NERT
#     CLOSE(18)

#     return






#       DO NRE=1,NENS
#          OPEN (19,file='output/el_conductivity',status='unknown',
#      1          form='unformatted')
#          J=ENPT(NRE)
#          DO I=1,N
#             PNEW(I)=ENPNEW(I,J)
#             SNODI(I)=ENSNODI(I,J)
#             PNODI(I)=ENPNODI(I,J)
#          END DO
#          CALL CHPARM(N,SNODI,PNODI,IVGHU)
#          CALL CHVELO(NLKP,N,NT,NTRI,IVGHU,TETRA,TP,
#      1               PNEW,SNODI,PNODI,SW,CKRW,SWE,CKRWE)
#          I=1
#          CONT=0
#          CONT1=0
#          DO K=1,N
# c              CONT1=CONT1+1
# c              CONT=CONT+1
# c              IF (CONT1 .EQ.20) I=I+2
# c              IF (CONT.EQ.204) I=I+10
# c              ELECTRC_NODI(K)=I
# c              IF (CONT1 .EQ.51) THEN
# c                 I  =I-2
# c                 CONT1=0
# c              END IF
# c              IF (CONT .EQ.612) THEN
# c                 I  =I-10
# c                 CONT=0
# c              END IF
#              ELECTRC_NODI(K)=
#      1          1/A_AR*CW_AR*PNODI(K)**M_AR*SW(K)**N_AR
#              WRITE(19)ELECTRC_NODI(k)
#          END DO
#          CLOSE(19)
#          call system('./prot.sh')
#          open(18,file='risp/ert.risp',status='unknown')
#          DO I=1,NERT
#              read(18,*) EN_ERT(I,J)
#          END DO
#          close(18)
#       END DO 
#       do i=1,nert      
#              write(777,*) i,(en_ert(i,j),j=1,nens)
#       end do
# C
#       RETURN
#       END
      
      


    




# src/ert_measure.f:C**************************  ERT_MEASURE *******************************
# src/ert_measure.f:      SUBROUTINE ERT_MEASURE(NLKP,N,NT,IVGHU,TETRA,TP,
# src/datin.f:         READ(IIN40,*) ERT_FLAG
# src/cathy_main.f:C                ENKFNOD(I,2) = 3 -->  measures of ERT (geolettriche)
# src/cathy_main.f:C   ERT_FLAG - flag for assimilation of ERT measures
# src/cathy_main.f:C          = 0 no ERT measures assimilated
# src/cathy_main.f:C          = 1 ERT measures assimilated 
# src/cathy_main.f:c          = 2 ERT measures as output in the detailed output times 
# src/cathy_main.f:C   NERT   - number of ERT observations (measurements) in the EnKF
# src/cathy_main.f:C                        streamflow or ERT observation values 
# src/cathy_main.f:C   EN_ERT(NERT,NENS)  - ERT measures associated with the ensemble
# src/cathy_main.f:C If ERT measures are assimilated read the parameters for archie's low
# src/cathy_main.f:                 CALL ERT_INIT(ARCHIE,NERT)
# src/cathy_main.f:                  CALL ERT_MEASURE(NLKP,N,NT,IVGHU,TETRA,TP,
# src/cathy_main.f:C  associated to the ensemble are computed in the subroutine ERT_MEASURE
# src/cathy_main.f:                  CALL ERT_MEASURE(NLKP,N,NT,IVGHU,TETRA,TP,
# src/cathy_main.f:                  CALL ERT_MEASURE(NLKP,N,NT,IVGHU,TETRA,TP,
# src/ert_init.f:      SUBROUTINE ERT_INIT(ARCHIE,NERT)
# src/ert_init.f:c      CALL ERT_MESH2D(LENGTH_X,LENGTH_Y,NCOL3D-1,NROW3D-1,TRIANG2D)
# src/ert_init.f:cc      CALL ERT_COOR2D(LENGTH_X,LENGTH_Y,NCOL3D-1,NROW3D-1,

                                       

#

# src/updatesir.f:     3                     DSMEASMAX,EN_ERT)
# src/updatesir.f:      INTEGER CONT_ERT
# src/updatesir.f:      REAL*8  SUMW,SUMPSI,SUMTHETA,SUMQ,D,SUMW1,SUM_ERT
# src/updatesir.f:      REAL*8  EN_ERT(MAXNUDN,MAXNENS)
# src/updatesir.f:      write(IOUT56,101)'NRE','SUMPSI','SUMTHETA','SUMQ','SUM_ERT',
# src/updatesir.f:         SUM_ERT=0.0d0
# src/updatesir.f:         CONT_ERT=0
# src/updatesir.f:               CONT_ERT=CONT_ERT+1
# src/updatesir.f:c               EN_ERT(CONT_ERT,J)=EN_ERT(CONT_ERT,J)/1000
# src/updatesir.f:               SUM_ERT=SUM_ERT+(ENKFVAL(CTR,I)-
# src/updatesir.f:     &                EN_ERT(CONT_ERT,J))**2/SIG2Y
# src/updatesir.f:         WSIRN(J)=WSIR(J)*exp(-0.5d0*(SUMPSI+SUMTHETA+SUMQ+SUM_ERT))
# src/updatesir.f:         WRITE(IOUT56,102)J,SUMPSI,SUMTHETA,SUMQ,SUM_ERT,WSIRN(J)
# src/ert_measure.f:C**************************  ERT_MEASURE *******************************
# src/ert_measure.f:      SUBROUTINE ERT_MEASURE(NLKP,N,NT,IVGHU,TETRA,TP,
# src/ert_measure.f:     1                       ENPNEW,ENPNODI,ENSNODI,ARCHIE,EN_ERT,
# src/ert_measure.f:     2                       ENPT,NERT,NTRI,NENS,ELECTRC_NODI,
# src/ert_measure.f:      INTEGER  NLKP,N,NT,IVGHU,NERT,NTRI,NENS
# src/ert_measure.f:      REAL*8   ARCHIE(4),EN_ERT(MAXNUDN,*),ELECTRC_NODI(*)
# src/ert_measure.f:         DO I=1,NERT
# src/ert_measure.f:             read(18,*) EN_ERT(I,J)
# src/pert_soil.f:C************************ PERT_SOIL **********************************
# src/pert_soil.f:      SUBROUTINE PERT_SOIL(NENS,N,NSTR,NZONE,NROW,NCOL,PERMX,PERMY,
# src/pert_ic.f:C************************ PERT_IC **************************************
# src/pert_ic.f:      SUBROUTINE PERT_IC(NENS,N,TETAF,PTIMEP,DSIC,ENPTIMEP,ENPNEW,
# src/pert_ic.f:      WRITE(IOUT56,103) 'NRE','PERTURB'
# src/newton.f:      SUBROUTINE NEWTON(CPUVEC,N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN,
# src/newton.f:      INTEGER  N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN
# src/newton.f:         NITERT=NITERT+NITER
# src/grdsys.f:     1                  NSTR,IVERT,LUMP,IMAX,
# src/grdsys.f:      INTEGER   NSTR,IVERT,LUMP,IMAX
# src/grdsys.f:      CALL GEN3D(N,NT,NNOD,NSTR,NTRI,IPRT1,IVERT,
# src/datin.f:      SUBROUTINE DATIN(ISIMGR,IVERT,ISP,WTPOSITION,BASE,ZRATIO, 
# src/datin.f:     J                 DAFLAG,ERT_FLAG,NENS,NOBS,DSRETC,
# src/datin.f:      INTEGER   ISIMGR,IVERT,ISP
# src/datin.f:      INTEGER   DAFLAG,ERT_FLAG,NENS,NOBS,ENKFT,NENSMIN,NEFFMIN
# src/datin.f:         READ(IIN2,*) IVERT,ISP,BASE
# src/datin.f:         WRITE(IOUT2,1350) IVERT,ISP,BASE
# src/datin.f:            READ(IIN11,*) IVERT,ISP,BASE
# src/datin.f:            IF (IVERT.EQ.3) THEN
# src/datin.f:            WRITE(IOUT40,1350) IVERT,ISP,BASE
# src/datin.f:         READ(IIN40,*) ERT_FLAG
# src/datin.f: 1350 FORMAT(/,5X,'IVERT  (TYPE OF VERTICAL DISCRETIZATION)= ',I6,
# src/picard.f:      SUBROUTINE PICARD(CPUVEC,N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN,
# src/picard.f:      INTEGER  N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN
# src/picard.f:      NITERT=NITERT+NITER
# src/gen3d.f:      SUBROUTINE GEN3D(N,NT,NNOD,NSTR,NTRI,IPRT1,IVERT,
# src/gen3d.f:      INTEGER   N,NT,NNOD,NSTR,NTRI,IPRT1,IVERT
# src/gen3d.f:            IF (IVERT .EQ. 0) THEN
# src/gen3d.f:               IF (IVERT .EQ. 1) THEN
# src/gen3d.f:               ELSE IF (IVERT .EQ. 2) THEN
# src/gen3d.f:               ELSE IF (IVERT .EQ. 3) THEN
# src/gen3d.f:	       ELSE IF (IVERT .EQ. 4) THEN
# src/inital.f:     2                  INDP,NITERT,ITLIN,ITRTOT,ITER,NSTEP,
# src/inital.f:      INTEGER  NITERT,ITLIN,ITRTOT,ITER,NSTEP
# src/inital.f:      CALL INIT1(NITERT,ITLIN,ITRTOT,ITER,NSTEP,
# src/flow3d.f:     2                  ITER,ITUNS,NITER,NITERT,ITLIN,ITRTOT,
# src/flow3d.f:      INTEGER  ITER,ITUNS,NITER,NITERT,ITLIN,ITRTOT
# src/flow3d.f:         CALL PICARD(CPUVEC,N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN,
# src/flow3d.f:         CALL NEWTON(CPUVEC,N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN,
# src/init1.f:      SUBROUTINE INIT1(NITERT,ITLIN,ITRTOT,ITER,NSTEP,
# src/init1.f:      INTEGER  NITERT,ITLIN,ITRTOT,ITER,NSTEP
# src/init1.f:      NITERT=0
# src/nudone.f:     1       /,5X,'NUDRZ  (VERTICAL RADIUS OF INFLUENCE)   = ',1PE15.5,
# Binary file src/updatesir.o matches
# Binary file src/datin.o matches
# src/cathy_main.f:C                ENKFNOD(I,2) = 3 -->  measures of ERT (geolettriche)
# src/cathy_main.f:C   Aug/06: Added IVERT=4 option. The base of the 3D grid is flat,
# src/cathy_main.f:C           IVERT=3. MC
# src/cathy_main.f:C   IVERT  - =0 each layer will be parallel to the surface, including
# src/cathy_main.f:C            (for ISP=0, IVERT=0, 1, and 2 yield the same 3-d mesh, 
# src/cathy_main.f:C   NITERT - number of iterations for the linear solver at each 
# src/cathy_main.f:C   ERT_FLAG - flag for assimilation of ERT measures
# src/cathy_main.f:C          = 0 no ERT measures assimilated
# src/cathy_main.f:C          = 1 ERT measures assimilated 
# src/cathy_main.f:c          = 2 ERT measures as output in the detailed output times 
# src/cathy_main.f:C   NERT   - number of ERT observations (measurements) in the EnKF
# src/cathy_main.f:C            (NERT < NOBS )
# src/cathy_main.f:C              mesh. For IVERT=0, BASE is subtracted from each surface
# src/cathy_main.f:C              be parallel to the surface. For IVERT=1 or 2, BASE is 
# src/cathy_main.f:C                        only if IVERT=3
# src/cathy_main.f:C                        is to occupy (see also description of IVERT). 
# src/cathy_main.f:C                        streamflow or ERT observation values 
# src/cathy_main.f:C   EN_ERT(NERT,NENS)  - ERT measures associated with the ensemble
# src/cathy_main.f:      INTEGER KPRT,IPEAT,IPRT,IPRT1,ISIMGR,IVERT,ISP,INDP
# src/cathy_main.f:      INTEGER ITMXCG,NITER,NITERT,ITLIN,KLSFAI,IMAX,MAXITER
# src/cathy_main.f:      INTEGER MINBOT,NDZ,NUDFLAG,DAFLAG,WFLAG,ERT_FLAG
# src/cathy_main.f:      INTEGER NEFFMIN,NSTEP1,ITRTOT1,NERT,NROUT
# src/cathy_main.f:      REAL*8  EN_ERT(MAXNUDN,MAXNENS)
# src/cathy_main.f:      CALL DATIN(ISIMGR,IVERT,ISP,WTPOSITION,BASE,ZRATIO, 
# src/cathy_main.f:     J           DAFLAG,ERT_FLAG,NENS,NOBS,DSRETC,
# src/cathy_main.f:     1        NSTR,IVERT,LUMP,IMAX,
# src/cathy_main.f:     2        INDP,NITERT,ITLIN,ITRTOT,ITER,NSTEP,
# src/cathy_main.f:         CALL PERT_SOIL(NENS,N,NSTR,NZONE,NROW,NCOL,PERMX0,PERMY0,
# src/cathy_main.f:            CALL PERT_IC(NENS,N,TETAF,PTIMEP,DSIC,ENPTIMEP,ENPNEW,
# src/cathy_main.f:C If ERT measures are assimilated read the parameters for archie's low
# src/cathy_main.f:            IF (ERT_FLAG.GE.1) THEN
# src/cathy_main.f:                 CALL ERT_INIT(ARCHIE,NERT)
# src/cathy_main.f:     1                     DFLOAT(NITERT)/DFLOAT(ITER),
# src/cathy_main.f:         NITERT=0
# src/cathy_main.f:     2        ITER,ITUNS,NITER,NITERT,ITLIN,ITRTOT,
# src/cathy_main.f:     1           HTIDIR,HTINEU,ITER,NITERT,KBACKT,KBACK,NNOD,NSF,
# src/cathy_main.f:     1        DFLOAT(NITERT)/DFLOAT(ITER),
# src/cathy_main.f:               IF (ERT_FLAG.EQ.2) THEN
# src/cathy_main.f:                  CALL ERT_MEASURE(NLKP,N,NT,IVGHU,TETRA,TP,
# src/cathy_main.f:     1                             ENPNEW,ENPNODI,ENSNODI,ARCHIE,EN_ERT,
# src/cathy_main.f:     2                             ENPT,NERT,NTRI,NENS,SCR,
# src/cathy_main.f:         CALL TIMNXT(NSTEP,ITER,NITERT,KBACKT,NSURFT,CPUSUB,CPUSURF)
# src/cathy_main.f:     1                        DFLOAT(NITERT)/DFLOAT(ITER),
# src/cathy_main.f:C  associated to the ensemble are computed in the subroutine ERT_MEASURE
# src/cathy_main.f:            IF (ERT_FLAG.EQ.1) THEN
# src/cathy_main.f:                  CALL ERT_MEASURE(NLKP,N,NT,IVGHU,TETRA,TP,
# src/cathy_main.f:     1                             ENPNEW,ENPNODI,ENSNODI,ARCHIE,EN_ERT,
# src/cathy_main.f:     2                             ENPT,NERT,NTRI,NENS,SCR,
# src/cathy_main.f:     3                DSMEASMAX,EN_ERT)
# src/cathy_main.f:     1        DFLOAT(NITERT)/DFLOAT(ITER),
# src/cathy_main.f:              IF (ERT_FLAG.EQ.2) THEN
# src/cathy_main.f:                  CALL ERT_MEASURE(NLKP,N,NT,IVGHU,TETRA,TP,
# src/cathy_main.f:     1                             ENPNEW,ENPNODI,ENSNODI,ARCHIE,EN_ERT,
# src/cathy_main.f:     2                             ENPT,NERT,NTRI,NENS,SCR,
# src/cathy_main.f:c            CALL TIMNXT(ENNSTEP(J),ITER,NITERT,KBACKT,NSURFT,CPUSUB,
# Binary file src/cathy matches
# src/bkstep.f:     1                  IETO,HTIDIR,HTINEU,ITER,NITERT,KBACKT,KBACK,
# src/bkstep.f:      INTEGER  ITER,NITERT,KBACKT,KBACK
# src/bkstep.f:      NITERT=0
# src/ert_init.f:      SUBROUTINE ERT_INIT(ARCHIE,NERT)
# src/ert_init.f:      INTEGER NERT
# src/ert_init.f:      READ(18,*) NERT
# src/ert_init.f:c      INTEGER NNOD2D,NTRI2D,N12D,NZONE2D,NERT
# src/ert_init.f:c      CALL ERT_MESH2D(LENGTH_X,LENGTH_Y,NCOL3D-1,NROW3D-1,TRIANG2D)
# src/ert_init.f:cc      CALL ERT_COOR2D(LENGTH_X,LENGTH_Y,NCOL3D-1,NROW3D-1,
# src/flow3d.f.old:     2                  ITER,ITUNS,NITER,NITERT,ITLIN,ITRTOT,
# src/flow3d.f.old:      INTEGER  ITER,ITUNS,NITER,NITERT,ITLIN,ITRTOT
# src/flow3d.f.old:         CALL PICARD(CPUVEC,N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN,
# src/flow3d.f.old:         CALL NEWTON(CPUVEC,N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN,
# Binary file src/pert_ic.o matches
# src/timnxt.f:      SUBROUTINE TIMNXT(NSTEP,ITER,NITERT,KBACKT,NSURFT,CPUSUB,CPUSURF)
# src/timnxt.f:      INTEGER   NSTEP,ITER,NITERT,KBACKT,NSURFT
# src/timnxt.f:      NITERT=0
# Binary file src/nudone.o matches


# elec.dat :
# contiene le nuove posizioni degli elettrodi nella griglia grande
# indicando il numero di elettrodi totali e, per ciascun elettrodo, il numero
# di nodo ed il corrispettivo numero di elettrodo;


#%% Archie

# src/ert_measure.f:C     Archie law = ELECTRC_NODI (for the 2D cross section only).
# src/cathy_main.f:C   ARCHIE(4)          - parameters for the application of Archie low:

# src/ert_measure.f:     1                       ENPNEW,ENPNODI,ENSNODI,ARCHIE,EN_ERT,
# src/ert_measure.f:      REAL*8   ARCHIE(4),EN_ERT(MAXNUDN,*),ELECTRC_NODI(*)
# src/ert_measure.f:      A_AR=ARCHIE(1)
# src/ert_measure.f:      CW_AR=ARCHIE(2)
# src/ert_measure.f:      M_AR=ARCHIE(3)
# src/ert_measure.f:      N_AR=ARCHIE(4)
# src/cathy_main.f:C   ARCHIE(4)          - parameters for the application of Archie low:
# src/cathy_main.f:C                      - ARCHIE(1) = a
# src/cathy_main.f:C                      - ARCHIE(2) = C_w
# src/cathy_main.f:C                      - ARCHIE(3) = m
# src/cathy_main.f:C                      - ARCHIE(4) = n
# src/cathy_main.f:      REAL*8  ARCHIE(4)
# src/cathy_main.f:                 CALL ERT_INIT(ARCHIE,NERT)
# src/cathy_main.f:     1                             ENPNEW,ENPNODI,ENSNODI,ARCHIE,EN_ERT,
# src/cathy_main.f:     1                             ENPNEW,ENPNODI,ENSNODI,ARCHIE,EN_ERT,
# src/cathy_main.f:     1                             ENPNEW,ENPNODI,ENSNODI,ARCHIE,EN_ERT,
# src/ert_init.f:      SUBROUTINE ERT_INIT(ARCHIE,NERT)
# src/ert_init.f:      REAL*8 ARCHIE(4)
# src/ert_init.f:      READ(18,*) ARCHIE(1)
# src/ert_init.f:      READ(18,*) ARCHIE(2)
# src/ert_init.f:      READ(18,*) ARCHIE(3)
# src/ert_init.f:      READ(18,*) ARCHIE(4)
# src/ert_init.f:c      REAL*8 ARCHIE(4),ZRATIO(NSTR),BASE


# -------------------------------------------------------------------#
#%% ERT DATA

# 2 files
# ERT measure
# ERT init