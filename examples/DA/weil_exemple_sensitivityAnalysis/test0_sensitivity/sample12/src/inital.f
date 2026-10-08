C
C**************************  INITAL ************************************
C
C  further input and initialization for ICs, BCs, nudging, parameters,
C  and misc counters, flags, and arrays
C
C***********************************************************************
C
      SUBROUTINE INITAL(NNOD,NSTR,N,NT,NTRI,NP,NQ,NSF,NDIR,NDIRC,NNEU,
     1                  NNEUC,HTIDIR,HTINEU,HSPATM,HTIATM,IETO,
     2                  IPRT1,IPOND,INDP,NITERT,ITLIN,ITRTOT,ITER,NSTEP,
     3                  NSURFT,NSURFT_T,NSURFT_TB,
     4                  KPRT,ISFCVG,KSFCV,KSFCVT,KBACKT,KBACK,KLSFAI,
     5                  IOPT,IPEAT,IVGHU,KSLOPE,HGFLAG,SFFLAG,
     6                  TETRA,TP,NODDIR,CONTP,CONTQ,IFATM,IFATMP,
     7                  NSFNUM,NSFNOD,SFEX,SFEXP,SFEXIT,DUPUIT,
     8                  NUDN,NUDCTR,NUDT,NUDC,NUDG,NUDFLAG,WFLAG,
     9                  NUDTIM,NUDRXY,NUDRZ,NUDTAU,
     A                  NUDX,NUDY,NUDZ,NUDEPS,NUDTET,NUDCUM,
     B                  DTGMIN,SFCHEK,SURF,DEM,PONDING,PONDP,
     C                  WTPOSITION,DELTAT,TMAX,TETAF,DTMIN,DTMAX,
     D                  DTSMAL,DTBIG,DTAVG,TIME,TIMEP,PONDH_MIN,
     E                  AREATOT,VOLTOT,RMAX,PEL,
     F                  X,Y,Z,XC,YC,ZC,VOLU,ARENOD,POROS,ELSTOR,
     G                  VGNCELL,VGRMCCELL,VGPSATCELL,SNODI,PNODI,DEF,
     H                  PORE,INDE,PNEW,PTIMEP,POLD,PTNEW,PTOLD,
     I                  SFQP,Q,QPOLD,QTIM,QINP,PRESC,PTIM,PINP,
     J                  ATMPOT,ATMACT,ATMOLD,ATMTIM,ATMINP,
     K                  OVFLNOD,OVFLP,PONDNOD,
     L                  NCELL,NCELNL,NUMRES,
     M                  CPUSUB,CPUSUB_T,CPUSURF,CPUSURF_T,CPUNL,CPUVEC,
     N                  ANP,ANQ,ACONTP,ACONTQ,
     O                  SFV,SFVNUM,SFVNOD,SFVTIM,QTRANIE,
     P                  IFSF,IFSFP,PUNTDIRSF_NODE,SEL,SW,SWP,
     R                  PUNTDIRFLOW_NODE,PUNTNEUFLOW_NODE,
     S                  CONDIR_NODE,CONNEU_NODE,
     T                  LEL,KEL,LAMBDA,LAMBDANODI,KD,KDNODI,
     V                  TRIANG)
     
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,j
      INTEGER  NNOD,NSTR,N,NT,NTRI,NSF,DUPUIT
      INTEGER  HTIDIR,HTINEU,HSPATM,HTIATM,IETO,IPRT1,IPOND,INDP
      INTEGER  NITERT,ITLIN,ITRTOT,ITER,NSTEP
      INTEGER  NSURFT,NSURFT_T,NSURFT_TB
      INTEGER  KPRT,ISFCVG,KSFCV,KSFCVT,KBACKT,KBACK,KLSFAI
      INTEGER  IOPT,IPEAT,IVGHU,KSLOPE
      INTEGER  NUDN,NUDCTR,NUDT,NUDC,NUDFLAG,WFLAG
      integer  NCELL,NCELNL,NUMRES
      INTEGER  ANP,ANQ
      INTEGER  NDIR(3),NDIRC(3),NP(3),NQ(3),NNEU(3),NNEUC(3)
      INTEGER  ACONTP(*),ACONTQ(*),TRIANG(4,*)
      INTEGER  HGFLAG(9),SFFLAG(5)
      INTEGER  TETRA(5,*),TP(*),NODDIR(*),CONTP(3,*),CONTQ(3,*)
      INTEGER  IFATM(*),IFATMP(*),IFSF(*),IFSFP(*),PUNTDIRSF_NODE(*)
      INTEGER  PUNTNEUFLOW_NODE(*),PUNTDIRFLOW_NODE(*) 
      INTEGER  NSFNUM(*),NSFNOD(NSFMAX,*)
      INTEGER  SFEX(NSFMAX,*),SFEXP(NSFMAX,*),SFEXIT(NSFMAX,*)
      INTEGER  NUDTET(*)
      integer  sfv(2),sfvnum(2,nsfmax),sfvnod(2,nsfmax,nnsfmx)
      LOGICAL  DTGMIN,SFCHEK,SURF,DEM,PONDING,PONDP
      REAL     CPUSUB,CPUSUB_T,CPUSURF,CPUSURF_T,CPUNL,CPUVEC(*)
      REAL*8   NUDG
      REAL*8   WTPOSITION,DELTAT,TMAX,TETAF,DTMIN,DTMAX
      REAL*8   DTSMAL,DTBIG,DTAVG,TIME,TIMEP,PONDH_MIN
      REAL*8   AREATOT,VOLTOT,RMAX
      REAL*8   NUDTIM(*),NUDRXY(*),NUDRZ(*),NUDTAU(MAXNUDT,*)
      REAL*8   NUDX(*),NUDY(*),NUDZ(*),NUDEPS(*),NUDCUM(*)
      REAL*8   X(*),Y(*),Z(*),XC(*),YC(*),ZC(*)
      REAL*8   VOLU(*),ARENOD(*),PEL(*)
      REAL*8   POROS(MAXSTR,*),ELSTOR(MAXSTR,*)
      REAL*8   VGNCELL(MAXSTR,*),VGRMCCELL(MAXSTR,*)
      REAL*8   VGPSATCELL(MAXSTR,*)
      REAL*8   LAMBDA(MAXSTR,*),KD(MAXSTR,*)
      REAL*8   SNODI(*),PNODI(*),PORE(*),INDE(*),DEF(*)
      REAL*8   LAMBDANODI(*),KDNODI(*)
      REAL*8   PNEW(*),PTIMEP(*),POLD(*),PTNEW(*),PTOLD(*)
      REAL*8   SFQP(NSFMAX,*),Q(*),QPOLD(*),QTIM(*),QINP(3,*)
      REAL*8   PRESC(*),PTIM(*),PINP(3,*)
      REAL*8   ATMPOT(*),ATMACT(*),ATMOLD(*),ATMTIM(*),ATMINP(3,*)
      REAL*8   OVFLNOD(*),OVFLP(*),PONDNOD(*)
      REAL*8   SEL(*),SW(*),SWP(*),LEL(*),KEL(*)
      REAL*8   CONDIR_NODE(*),CONNEU_NODE(*)
      real*8   sfvtim(2)
      REAL*8   QTRANIE(NMAX)
      INCLUDE 'MB_HGRAPH.H'
      INCLUDE 'SOILCHAR.H'
      INCLUDE 'SURFWATER.H'
      INCLUDE 'RIVERNETWORK.H'
      INCLUDE 'IOUNITS.H'
C
C  calculate fully (INDP=2) or partially (INDP=3) saturated
C  vertical hydrostatic equilibrium IC's
C
      IF (INDP .EQ. 2) THEN
         CALL ICVHE(NNOD,NSTR,N,IPRT1,IPOND,PTIMEP,PONDNOD,Z)
      ELSE IF (INDP .EQ. 3) THEN
         CALL ICVHWT(NNOD,NSTR,N,IPRT1,IPOND,WTPOSITION,
     1               PTIMEP,PONDNOD,Z)
      ELSE IF (INDP.EQ.4) THEN
         CALL ICVDWT(NNOD,NSTR,N,IPRT1,IPOND,WTPOSITION,
     1               PTIMEP,PONDNOD,Z)
      END IF
C
C  set up TMAX, TETAF, and time step sizes for steady state problem
C
      IF (DELTAT .GE. 1.0D+15) THEN
         TMAX=0.0D0
         TETAF=1.0D0
         DELTAT=RMAX
         DTMIN=RMAX
         DTMAX=RMAX
      END IF
C
C  initialization of various counters, flags, and arrays
C
      CALL INIT0R(NT,PEL)      
      CALL INIT0R(NT,SEL)
      CALL INIT1(NITERT,ITLIN,ITRTOT,ITER,NSTEP,
     1           DELTAT,DTMAX,DTMIN,DTGMIN,TIMEP,TIME,
     2           RMAX,DTSMAL,DTBIG,DTAVG,
     3           KPRT,KSFCV,KSFCVT,KBACKT,KBACK,KLSFAI,
     4           N,TP,SNODI,PNODI,PORE,INDE,PNEW,PTIMEP,DEF,
     5           NP,NODDIR,CONTP,QPOLD,
     6           HGFLAG,SFFLAG,
     7           CPUSUB,CPUSUB_T,CPUSURF,CPUSURF_T,CPUNL,CPUVEC)
C
C  calculate total catchment area and volume
C
      AREATOT = 0.0D0
      DO I=1,NNOD
         AREATOT = AREATOT + ARENOD(I)
      END DO
      VOLTOT = 0.0D0
      DO I=1,NT
         VOLTOT = VOLTOT + VOLU(I)
      END DO
      IF (DEM) THEN
         WRITE(IOUT2,1000) AREATOT,NCELL*DELTA_X*DELTA_Y,VOLTOT
      ELSE
         WRITE(IOUT2,1020) AREATOT,VOLTOT
      END IF
      IF (SURF) WRITE(IOUT43,1100) NNOD,NCELL,AREATOT
      WRITE(IOUT43,1110)
      
      CALL INIT0R(NMAX,SW)
      CALL INIT0R(NMAX,SWP)
C
C  initialization for SURF_ROUTE (for the subsurface-surface coupling case)
C
      CALL INIT_SR(NNOD,NCELNL,NUMRES,NSURFT,NSURFT_T,NSURFT_TB,
     1     OVFLNOD,OVFLP)
C
C  read and initialize non-atmospheric, non-seepage face Dirichlet
C  boundary conditions parameters and arrays for first time step
C  (times 0.0 and DELTAT)
C
C     IF (NP(3) .NE. 0) THEN
      CALL BCONE('   NATM, NSF DIRICHLET',IIN8,IOUT2,IOUT19,N,NNOD,
     1     NPMAX,IPRT1,NDIR,NDIRC,NP,NSTR,HTIDIR,TIME,
     2     DELTAT,PTIM,PINP,CONTP,PRESC,ANP,ACONTP,
     3     PUNTDIRFLOW_NODE,CONDIR_NODE)

C     WRITE(222,*) TIME                                       
C     WRITE(222,*) ANP                                        
C     WRITE(222,*) (ACONTP(I),I=1,ANP)                        
C     WRITE(222,*) (PRESC(I),I=1,ANP)
C     END IF
C
C  read and initialize non-atmospheric, non-seepage face Neumann
C  boundary conditions parameters and arrays for first time step
C  (times 0.0 and DELTAT)
C
C     IF (NQ(3) .NE. 0) THEN
C     CALL INIT0I(3,ZERO)
C     CALL INIT0I(3,ZEROC)
      CALL BCONE('     NATM, NSF NEUMANN',IIN9,IOUT2,IOUT20,N,NNOD,
     1              NQMAX,IPRT1,NNEU,NNEUC,NQ,NSTR,HTINEU,TIME,
     2              DELTAT,QTIM,QINP,CONTQ,Q,ANQ,ACONTQ,
     3              PUNTNEUFLOW_NODE,CONNEU_NODE)

C     WRITE(223,*) TIME                                       
C     WRITE(223,*) ANP                                        
C     WRITE(223,*) (ACONTP(I),I=1,ANP)                        
C     WRITE(223,*) (PRESC(I),I=1,ANP)
C     END IF
C  
C  Initialize QPOLD and QTRANIE, since init1.f doesn't anymore
C
      CALL INIT0R(ANP,QPOLD)
      CALL INIT0R(N,QTRANIE)
C     WRITE(224,*) TIME                                       
C     WRITE(224,*) ANP                                        
C     WRITE(224,*) (ACONTP(I),I=1,ANP)                        
C     WRITE(224,*) (PRESC(I),I=1,ANP)
C
C  if coupled SURF_ROUTE-FLOW3D simulation but no initial ponding
C  (IPOND = 0) set the flag PONDING to FALSE (so that SURF_ROUTE
C  is not called for the first time step)
C
      IF (IPOND .NE. 0) THEN
         PONDING = .TRUE.
      ELSE
         PONDING = .FALSE.
      END IF
      PONDP = PONDING
C
C unit IIN7 initial input: seepage face boundary conditions
C
      call sfvone(iin7,iterm,iout2,iprt1,n,nsfmax,nnsfmx,
     1     time,nsf,nsfnum,nsfnod,sfv,sfvnum,
     2     sfvnod,sfvtim,ifsf,ifsfp,z)
      
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            SFQP(I,J)=0.0D0
         END DO
      END DO
      IF (ISFCVG .EQ. 1 .AND. NSF .GT. 0) THEN
         SFCHEK=.TRUE.
      ELSE
         SFCHEK=.FALSE.
      END IF
c
C  calculate seepage face exit points based on initial pressure heads
C
      CALL SFINIT(N,NSF,NSFNUM,NSFNOD,SFEX,SFEXP,SFEXIT,DUPUIT,
     1            PTIMEP,PNEW,SFFLAG,DELTAT,Z,
     1            PUNTDIRSF_NODE)
C
C  read and initialize atmospheric boundary condition parameters and
C  arrays for first time step (times 0.0 and DELTAT)
C
      CALL ATMONE(NNOD,HSPATM,HTIATM,IETO,TIME,DELTAT,PONDH_MIN,
     1     IFATM,IFATMP,ARENOD,ATMPOT,ATMACT,ATMOLD,ATMTIM,
     2     ATMINP,PNEW,PTIMEP,ANP,ANQ,ACONTP,
     3     ACONTQ,NSF,NSFNUM,NSFNOD)
C
C  read and initialize nudging parameters and arrays
C
      CALL NUDONE(N,NT,TETRA,X,Y,Z,XC,YC,ZC,VOLU,
     1            NUDN,NUDCTR,NUDT,NUDC,NUDG,NUDFLAG,WFLAG,
     2            NUDTIM,NUDRXY,NUDRZ,NUDTAU,
     3            NUDX,NUDY,NUDZ,NUDEPS,NUDTET,NUDCUM)
c     nudging is not yet supported for Newton iteration
      IF (IOPT .NE. 1  .AND.  NUDN .NE. 0) THEN
         WRITE(IOUT2,1300)
         CALL CLOSIO
         STOP
      END IF
C     WRITE(227,*) TIME                                       
C     WRITE(227,*) ANP                                        
C     WRITE(227,*) (ACONTP(I),I=1,ANP)                        
C     WRITE(227,*) (PRESC(I),I=1,ANP)
C
C  initialize pressure heads for first time step
C
      CALL VCOPYR(N,POLD,PNEW)
      CALL WEIGHT(N,TETAF,PNEW,PTIMEP,PTNEW)
      CALL VCOPYR(N,PTOLD,PTNEW)
C     WRITE(228,*) TIME                                       
C     WRITE(228,*) ANP                                        
C     WRITE(228,*) (ACONTP(I),I=1,ANP)                        
C     WRITE(228,*) (PRESC(I),I=1,ANP)
C
C  initialize mass balance and hydrograph parameters for first time
C  level (time 0.0)
C
      CALL MBINIT(SFFLWP,AACTP,ADINP,ADOUTP,NDINP,NDOUTP,
     1            ANINP,ANOUTP,NNINP,NNOUTP,
     2            NUDINP,NUDOUTP,VNUDTOT,
     3            NNOD,IFATMP,ATMOLD,ANQ,Q,
     4            VAPOT_T,VAACT_T,VSFTOT,VNDTOT,VNNTOT,VTOT,
     5            STORE2,CDSTOR,CVIN,CVOUT,CERRAS,CAERAS)
C     WRITE(229,*) TIME                                       
C     WRITE(229,*) ANP                                        
C     WRITE(229,*) (ACONTP(I),I=1,ANP)                        
C     WRITE(229,*) (PRESC(I),I=1,ANP)
C
C  distribute on a nodal basis parameters which are input on an
C  element basis (POROS, ELSTOR)
C
        CALL TPNODI(N,NT,NTRI,TETRA,TRIANG,PNODI,SNODI,TP,IPEAT,
     1            POROS,ELSTOR,PORE,INDE,DEF,PEL,SEL,
     2            VGNCELL,VGRMCCELL,VGPSATCELL,
     3            LEL,LAMBDA,LAMBDANODI,KEL,KD,KDNODI)

C     WRITE(230,*) TIME                                       
C     WRITE(230,*) ANP                                        
C     WRITE(230,*) (ACONTP(I),I=1,ANP)                        
C     WRITE(230,*) (PRESC(I),I=1,ANP)
C
C  calculate miscellaneous van Genuchten, Huyakorn, and
C  Brooks-Corey moisture curve parameters, including VGPNOT for
C  van Genuchten and extended van Genuchten cases and BCPORM
C  for Brooks-Corey case
C
      IF (IVGHU .NE. -1) CALL CHPARM(N,SNODI,PNODI,IVGHU)
C     WRITE(231,*) TIME                                       
C     WRITE(231,*) ANP                                        
C     WRITE(231,*) (ACONTP(I),I=1,ANP)                        
C     WRITE(231,*) (PRESC(I),I=1,ANP)
C
C  calculate localized tangent slope approximations of moisture
C  curve derivatives for the case KSLOPE=4
C
      IF (IVGHU .NE. -1) THEN
         IF (KSLOPE .EQ. 4) THEN
            IF (IOPT .EQ. 1) THEN
               CALL CHTANP(N,SNODI,PNODI,IVGHU)
            ELSE
               CALL CHTANN(N,SNODI,PNODI,IVGHU)
            END IF
         END IF
      END IF
C     WRITE(232,*) TIME                                       
C     WRITE(232,*) ANP                                        
C     WRITE(232,*) (ACONTP(I),I=1,ANP)                        
C     WRITE(232,*) (PRESC(I),I=1,ANP)
C
      RETURN
 1000 FORMAT(/,5X,'AREATOT (TOTAL CATCHMENT SURFACE AREA)  = ',1PE15.5,
     1       /,5X,'(SURFACE AREA FROM NCELL*DELTAX*DELTAY) = ',1PE15.5,
     2       /,5X,'VOLTOT  (TOT VOL OF DISCRETIZED CATCH.) = ',1PE15.5)
 1020 FORMAT(/,5X,'AREATOT (TOTAL CATCHMENT SURFACE AREA)  = ',1PE15.5,
     1       /,5X,'VOLTOT  (TOT VOL OF DISCRETIZED CATCH.) = ',1PE15.5)
 1100 FORMAT('#',17X,' ***** Surface vs subsurface diagnostics ***** ',
     1       /,'#','NNOD    (# of surface nodes)            = ',I6,
     2       /,'#','NCELL   (# of DEM cells)                = ',I6,
     3       /,'#','AREATOT (total catchment surface area)  = ',1PE15.5)
 1110 FORMAT(  '#','Step      (1) : Time step',
     1       /,'#','Deltat    (2) : Time step size',
     2       /,'#','Time      (3) : See parm input file for units',
     3       /,'#','Back      (4) : # of back-stepping occurrences',
     4       /,'#','NL-l      (5) : # of nonlinear iterations for the',
     5            ' successful (last) time step',
     6       /,'#','NL-a      (6) : # of nonlinear iterations for the',
     7            ' successful time step and any back-steps (= NL-last',
     8            ' + ITUNS*Back)',
     9       /,'#','Sdt-l     (7) : # of time steps in the surface',
     A            ' routing module for the successful (last)',
     B            ' subsurface module time step',
     C       /,'#','Sdt-a     (8) : # of time steps in the surface',
     D            ' routing module for the successful subsurface',
     E            ' module time step and any back-steps',
     F       /,'#','Atmpot-vf (9) : Potential atmospheric forcing',
     G            ' (rain +ve / evap -ve) as a volumetric flux [L^3/T]',
     H       /,'#','Atmpot-v (10) : Potential atmospheric forcing',
     I            ' volume [L^3] (See parm input file for units)',
     J       /,'#','Atmpot-r (11) : Potential atmospheric forcing',
     K            ' rate [L/T]',
     L       /,'#','Atmpot-d (12) : Potential atmospheric forcing',
     M            ' depth [L]',
     N       /,'#','Atmact-vf(13) : Actual infiltration (+ve) or',
     O            ' exfiltration (-ve) at atmospheric BC nodes as a',
     P            ' volumetric flux [L^3/T]',
     Q       /,'#','Atmact-v (14) : Actual infiltration (+ve) or',
     R            ' exfiltration (-ve) volume [L^3]',
     S       /,'#','Atmact-r (15) : Actual infiltration (+ve) or',
     T            ' exfiltration (-ve) rate [L/T]',
     U       /,'#','Atmact-d (16) : Actual infiltration (+ve) or',
     V            ' exfiltration (-ve) depth [L]',
     W       /,'#','Horton   (17) : Fraction of the surface nodes',
     X            ' that are saturated or ponded due to Horton',
     Y            ' infiltration excess (Note: based on PNEW and',
     Z            ' not on IFATM)',
     A       /,'#','Dunne    (18) : Fraction of the surface nodes',
     B            ' that are saturated or ponded due to Dunne',
     C            ' saturation excess (see previous note)',
     D       /,'#','Ponded   (19) : Fraction of the surface nodes',
     E            ' that are ponded (PNEW > PONDH_MIN)',
     F       /,'#','Satur    (20) : Fraction of the surface nodes',
     G            ' that are saturated or ponded (PNEW > 0)',
     H       /,'#','CPU-sub  (21) : CPU seconds for the subsurface',
     I            ' flow module',
     J       /,'#','CPU-surf (22) : CPU seconds for the surface',
     K            ' routing module',
     1       /,'#   (1)       (2)       (3)   (4)   (5)   (6)   (7)',
     2         '   (8)        (9)       (10)       (11)       (12)',
     3         '       (13)       (14)       (15)       (16)',
     4         '   (17)   (18)   (19)   (20)       (21)       (22)',
     5       /,'#  Step    Deltat      Time  Back  NL-l  NL-a Sdt-l',
     6         ' Sdt-a  Atmpot-vf   Atmpot-v   Atmpot-r   Atmpot-d',
     7         '  Atmact-vf   Atmact-v   Atmact-r   Atmact-d',
     8         ' Horton  Dunne Ponded  Satur   CPU-sub  CPU-surf')
 1300 FORMAT(/,5X,'INVALID IOPT AND NUDN COMBINATION (nudging is not ',
     1            'yet supported for ',
     2       /,5X,'the Newton scheme)')
      END
