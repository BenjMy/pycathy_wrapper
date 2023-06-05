C
C**************************  FLOW3D ************************************
C
C  subsurface three-dimensional variably saturated flow solver
C
C***********************************************************************
C
      SUBROUTINE FLOW3D(IOPT,IPRT1,NLRELX,L2NORM,N,NT,NTRI,NTERM,NSTR,
     1                  NDIR,NDIRC,NP,NQ,HSPATM,HTIATM,HTIDIR,HTINEU,
     2                  ITER,ITUNS,NITER,NITERT,ITLIN,ITRTOT,
     3                  ISOLV,ITMXCG,IPEAT,IVGHU,KSLOPE,ISFONE,KSFCV,
     4                  KSFCVT,IBOT,MINBOT,NNOD,NUMDIR,NSF,KLSFAI,
     5                  NLKP,NUDN,NUDC,NUDCTR,NUDFLAG,WFLAG,
     6                  TP,TETRA,TETJA,IA,JA,TOPOL,CONTP,CONTQ,IVOL,
     7                  INSYM,IFATM,IFATMP,NODDIR,NSFNUM,NSFNOD,
     8                  SFEX,SFEXIT,SFEXP,SFFLAG,NUDTET,DUPUIT,
     9                  SURF,PONDING,
     A                  DTGMIN,LSFAIL,ERRGMX,NORMCV,ITAGEN,SFCHEK,
     B                  KSFZER,NOBACK,BCKSTP,CPUVEC,CPUNL,TIME,TIMEP,
     C                  DELTAT,DTMIN,DTREDM,DTREDS,TETAF,OMEGA,OMEGAP,
     D                  RMAX,RMIN,TOLUNS,TOLSWI,TOLCG,ERNLMX,
     E                  PONDH_MIN,NUDG,
     F                  PONDNOD,OVFLNOD,PERMX,PERMY,PERMZ,INDE0,
     G                  SNODI,PNODI,PEL,PORE,INDE,PRESC,PTIM,PINP,DEF,
     H                  Q,QTIM,QINP,ATMTIM,ATMINP,AI,BI,CI,DI,LMASS,
     I                  ARENOD,VOLNOD,VOLU,VOLUR,
     J                  ETAI,ETAE,ET1,ET2,ET1E,DETAI,DETAIE,SW,SWTIMEP,
     K                  SWE,CKRW,CKRWE,DCKRW,DCKRWE,SENODI,SEELT,
     L                  POLD,PDIFF,PTOLD,PNEW,PTNEW,PTIMEP,
     M                  TNOTI,XT5,LHSP,COEF1,COEF2,COEF3,COEF4,SCR,
     N                  RNSYM,ATMACT,ATMPOT,ATMOLD,LHSATM,LHSSF,
     O                  QPNEW,QPOLD,SFQ,SFQP,BKCFLOW_NODE,X,Y,Z,
     P                  NUDTIM,NUDRXY,NUDRZ,NUDX,NUDY,NUDZ,NUDCUM,
     Q                  NUDEPS,NUDDIF,NUDSMC,NUDNOD,NUDVAL,NUDTAU,
     R                  PUNTDIRSF_NODE,QTRANIE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  IKMAX,KSF
      INTEGER  I,II,J
      INTEGER  IOPT,IPRT1,NLRELX,L2NORM,N,NT,NTRI,NTERM,NSTR
      INTEGER  NDIR,NDIRC,NP,NQ,HSPATM,HTIATM,HTIDIR,HTINEU
      INTEGER  ITER,ITUNS,NITER,NITERT,ITLIN,ITRTOT
      INTEGER  ISOLV,ITMXCG,IPEAT,IVGHU,KSLOPE,ISFONE,KSFCV,KSFCVT
      INTEGER  IBOT,MINBOT,NNOD,NUMDIR,NSF,KLSFAI,DUPUIT
      INTEGER  NLKP,NUDN,NUDC,NUDCTR,NUDFLAG,WFLAG
      INTEGER  TETRA(5,*),TETJA(4,4,*),IA(*),JA(*),TOPOL(*)
      INTEGER  CONTP(*),CONTQ(*),IVOL(*),TP(*)
      INTEGER  INSYM(*),IFATM(*),IFATMP(*),NODDIR(*)
      INTEGER  NSFNUM(*),NSFNOD(NSFMAX,*)
      INTEGER  SFEX(NSFMAX,*),SFEXIT(NSFMAX,*),SFEXP(NSFMAX,*)
      INTEGER  SFFLAG(*),NUDTET(*),PUNTDIRSF_NODE(*)
      LOGICAL  SURF,PONDING,DTGMIN,LSFAIL,ERRGMX,NORMCV,ITAGEN,SFCHEK
      LOGICAL  KSFZER,NOBACK,BCKSTP
      REAL     CPUV9,CPUNL1
      REAL     CPUVEC(*),CPUNL
      REAL*8   PIKMAX,PINF,PL2,FINF,FL2
      REAL*8   TIME,TIMEP
      REAL*8   DELTAT,DTMIN,DTREDM,DTREDS,TETAF,OMEGA,OMEGAP
      REAL*8   RMAX,RMIN,TOLUNS,TOLSWI,TOLCG,ERNLMX
      REAL*8   PONDH_MIN,NUDG
      REAL*8   PONDNOD(*),OVFLNOD(*)
      REAL*8   PERMX(MAXSTR,*),PERMY(MAXSTR,*),PERMZ(MAXSTR,*)
      REAL*8   SNODI(*),PNODI(*),PEL(*),PRESC(*),PTIM(*),PINP(3,*)
      REAL*8   PORE(*),INDE(*),DEF(*),INDE0(*)
      REAL*8   Q(*),QTIM(*),QINP(3,*),QTRANIE(*)
      REAL*8   ATMTIM(*),ATMINP(3,*)
      REAL*8   AI(4,*),BI(4,*),CI(4,*),DI(4,*),LMASS(4,4)
      REAL*8   ARENOD(*),VOLNOD(*),VOLU(*),VOLUR(*)
      REAL*8   ETAI(*),ETAE(*),DETAI(*),DETAIE(*),SW(*),SWE(*)
      REAL*8   ET1(*),ET2(*),ET1E(*),SWTIMEP(*)
      REAL*8   CKRW(*),CKRWE(*),DCKRW(*),DCKRWE(*),SENODI(*),SEELT(*)
      REAL*8   POLD(*),PDIFF(*),PTOLD(*),PNEW(*),PTNEW(*),PTIMEP(*)
      REAL*8   TNOTI(*),XT5(*),LHSP(*),BKCFLOW_NODE(*)
      REAL*8   COEF1(*),COEF2(*),COEF3(*),COEF4(*),SCR(*),RNSYM(*)
      REAL*8   ATMACT(*),ATMPOT(*),ATMOLD(*),LHSATM(*),LHSSF(NSFMAX,*)
      REAL*8   QPNEW(*),QPOLD(*),SFQ(NSFMAX,*),SFQP(NSFMAX,*)
      REAL*8   X(*),Y(*),Z(*)
      REAL*8   NUDTIM(*),NUDRXY(*),NUDRZ(*)
      REAL*8   NUDX(*),NUDY(*),NUDZ(*),NUDCUM(*)
      REAL*8   NUDEPS(*),NUDDIF(*),NUDSMC(*),NUDNOD(*)
      REAL*8   NUDVAL(MAXNUDC,*),NUDTAU(MAXNUDT,*)
      REAL*8   SWNEW(NMAX)
      INCLUDE 'NORMVL.H'
      INCLUDE 'MB_HGRAPH.H'
      INCLUDE 'IOUNITS.H'
C           
C----------------------------------------------START NONLINEAR ITERATION
C           
      CALL TIM(CPUNL1,1)
      WRITE(IOUT2,1050)
      WRITE(ITERM,1020)
  200 CONTINUE
      IF (IOPT .EQ. 1) THEN
         CALL PICARD(CPUVEC,N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN,
     1               ITMXCG,TOLCG,TIME,DELTAT,TETAF,IPEAT,IVGHU,
     2               NLKP,NUDN,NUDC,NUDCTR,NUDG,NUDFLAG,WFLAG,
     3               KSLOPE,TP,TETRA,TETJA,JA,TOPOL,
     4               PERMX,PERMY,PERMZ,SNODI,PNODI,PEL,CONTP,CONTQ,
     5               PORE,INDE,PRESC,Q,AI,BI,CI,DI,LMASS,X,Y,Z,
     6               INDE0,VOLNOD,VOLU,VOLUR,IVOL,RMAX,RMIN,
     7               ETAI,ETAE,ET1,ET2,ET1E,SW,SWE,CKRW,CKRWE,
     C               SENODI,SEELT,SWNEW,SWTIMEP,
     8               POLD,PDIFF,PTOLD,PNEW,DEF,
     9               PTNEW,PTIMEP,TNOTI,XT5,LHSP,COEF1,COEF2,COEF4,
     A               SCR,NNOD,IFATM,ATMACT,ATMOLD,LHSATM,NUMDIR,
     B               NODDIR,NSF,NSFNUM,NSFNOD,SFEX,LHSSF,DUPUIT,
     C               LSFAIL,KLSFAI,
     D               NUDTET,NUDTIM,NUDRXY,NUDRZ,NUDX,NUDY,NUDZ,
     E               NUDEPS,NUDDIF,NUDSMC,NUDNOD,NUDVAL,NUDTAU,QTRANIE)
      ELSE  
         CALL NEWTON(CPUVEC,N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN,
     1               ISOLV,ITMXCG,TOLCG,TIME,DELTAT,TETAF,IVGHU,
     2               NLKP,KSLOPE,TP,TETRA,TETJA,IA,JA,TOPOL,
     3               PERMX,PERMY,PERMZ,Z,
     4               SNODI,PNODI,CONTP,CONTQ,PRESC,Q,BI,CI,DI,
     5               LMASS,VOLU,VOLUR,IVOL,RMAX,RMIN,ETAI,ETAE,
     6               DETAI,DETAIE,SW,SWE,CKRW,CKRWE,DCKRW,DCKRWE,
     7               POLD,PDIFF,PTOLD,PNEW,PTNEW,PTIMEP,
     8               TNOTI,XT5,LHSP,COEF1,COEF2,COEF3,SCR,
     9               IBOT,MINBOT,INSYM,RNSYM,NNOD,IFATM,
     A               ATMACT,ATMOLD,LHSATM,NUMDIR,NODDIR,NSF,
     B               NSFNUM,NSFNOD,SFEX,LHSSF,LSFAIL,KLSFAI,DUPUIT)
      END IF
C     WRITE(5000,*)ITER,INDE(13)
C
C  before checking for convergence of the nonlinear scheme, we:
C  (i)   back-calculate fluxes at Dirichlet nodes
C  (ii)  compute (but don't output) mass balance errors
C  (iii) apply nonlinear relaxation (if required)
C  (iv)  calculate the nonlinear convergence and residual error norms
C  (v)   check for switching of atmospheric boundary conditions
C        (depending on the setting of TOLSWI)
C  (vi)  calculate new position of the exit point along each seepage
C        face and check for seepage face exit point convergence
C  Steps (i) and (ii) are done first because pressure head and boundary
C  flux values may be modified in steps (iii), (v), and (vi),
C  and we want to use the unmodified values in the mass balance
C  calculations (since these are the values obtained after solving
C  the system of equations). We consider the modified values however to
C  be the solution for the current iteration or time level, in the sense
C  that these are the values which are saved (in POLD, PTOLD, PTIMEP,
C  ATMOLD, SFQP, etc) before we pass on to the next iteration or
C  time level.
C  Mass balance errors are only output at the end of the nonlinear
C  iterative procedure.
C  It is necessary to perform step (vi) before checking for convergence
C  of the nonlinear scheme since in the case SFCHEK=TRUE the
C  seepage face exits points must also converge for the nonlinear
C  iterations to converge.
C
      CALL TIM(CPUV9,1)
      CALL MASBAL(IOPT,N,NP,NQ,NNOD,NSF,IPEAT,IVGHU,NUDC,
     1            NLKP,NT,NTRI,TETAF,TIME,DELTAT,
     2            TETRA,TP,TOPOL,JA,CONTP,IFATM,
     3            SFEX,SFFLAG,NSFNUM,NSFNOD,
     4            PNEW,PTIMEP,PDIFF,SNODI,PNODI,PORE,INDE,
     5            SENODI,ETAI,SEELT,ETAE,Z,VOLNOD,QPNEW,QPOLD,Q,
     6            ATMACT,ATMOLD,COEF1,XT5,SCR,SWNEW,SWTIMEP,
     7            BKCFLOW_NODE,SFQ,SFQP,NUDNOD,NUDCUM)
C  
C  apply nonlinear relaxation (if required)
C        
      IF (NLRELX .NE. 0) THEN
         IF (NLRELX .EQ. 1) THEN
            CALL RELAX(N,OMEGA,PNEW,POLD)
         ELSE
            CALL RELXOM(N,ITER,OMEGA,OMEGAP,PIKMXV,PNEW,POLD)
            CALL RELAX(N,OMEGA,PNEW,POLD)
         END IF   
      END IF      
C
C  check nonlinear convergence and switching of atmospheric
C  and seepage face boundary conditions
C
      CALL CONVER(ITER,N,NNOD,NSF,KSFCV,KSFCVT,NITER,ISFONE,
     1            L2NORM,IKMAX,KSF,NSFNUM,NSFNOD,SFEX,SFEXIT,
     2            SFFLAG,IFATM,SURF,PONDING,KSFZER,TOLSWI,
     3            TIME,DELTAT,PIKMAX,PINF,PL2,FINF,FL2,PONDH_MIN,
     4            PNEW,POLD,TNOTI,ATMACT,ATMPOT,SFQ,QTRANIE,
     5            ARENOD,PONDNOD,OVFLNOD,DUPUIT,Z,PUNTDIRSF_NODE)
      CALL TIM(CPUV9,2)
      CPUVEC(9)=CPUVEC(9)+CPUV9
C     
C------------------------------CHECK FOR CONVERGENCE OF NONLINEAR SCHEME
C     
C  case   (i) LSFAIL = FALSE and ERRGMX = FALSE and NORMCV = FALSE and
C             ITAGEN = TRUE
C             Iterate again: update pressure heads for next iteration.
C  case  (ii) LSFAIL = FALSE and ERRGMX = FALSE and NORMCV = TRUE
C             Convergence achieved: output mass balance errors;
C             input (if necessary), interpolate, and check for switching
C             of atmospheric boundary conditions; update pressure
C             heads for next time level.
C  case (iii) (LSFAIL = TRUE or ERRGMX = TRUE or (NORMCV = FALSE and
C             ITAGEN = FALSE)) and DTGMIN = TRUE
C             Convergence failed - back-step:
C             (see BKSTEP routine - reset values to previous
C             time level; reduce time step size; interpolate
C             non-atmospheric, non-seepage face Dirichlet and Neumann
C             boundary conditions; interpolate and check
C             for switching of atmospheric boundary conditions; re-start
C             the time step using smaller time step size).
C  case  (iv) (LSFAIL = TRUE or ERRGMX = TRUE or (NORMCV = FALSE and
C             ITAGEN = FALSE)) and DTGMIN = FALSE
C             Convergence failed - no back-stepping possible: output
C             mass balance errors; calculate hydrograph fluxes;
C             end the simulation.
C  Note: in the case SFCHEK=TRUE an additional condition for convergence
C  is that the seepage face exit points converged (i.e. KSFZER=TRUE).
C        
      BCKSTP=.FALSE.
      IF (ITER .LT. ITUNS) THEN
         ITAGEN=.TRUE.
      ELSE  
         ITAGEN=.FALSE.
      END IF      
      IF (PL2 .GE. ERNLMX  .OR.  PINF .GE. ERNLMX  .OR.
     1    FL2 .GE. ERNLMX  .OR.  FINF .GE. ERNLMX) THEN
         ERRGMX=.TRUE.
      ELSE 
         ERRGMX=.FALSE.
      END IF
      IF (L2NORM .EQ. 0) THEN
         IF (PINF .LE. TOLUNS) THEN
            NORMCV=.TRUE.
         ELSE 
            NORMCV=.FALSE.
         END IF
      ELSE    
         IF (PL2 .LE. TOLUNS) THEN
            NORMCV=.TRUE.
         ELSE 
            NORMCV=.FALSE.
         END IF
      END IF  
      IF (     ( (.NOT. LSFAIL) .AND. (.NOT. ERRGMX) .AND.
     1           (.NOT. NORMCV) .AND. ITAGEN )
     2    .OR. ( (.NOT. LSFAIL) .AND. (.NOT. ERRGMX) .AND.
     3           ITAGEN .AND. SFCHEK .AND. (.NOT. KSFZER) ) ) THEN  
C        
C  case (i): iterate again
C             
c         CALL VCOPYI(NSF,SFEXIT,SFEX) 
        DO I=1,NSF
           DO J=1,NSFNUM(I)
               SFEXIT(I,J)=SFEX(I,J)
            END DO
         END DO
         CALL VCOPYR(N,POLD,PNEW)
         CALL VCOPYR(N,PTOLD,PTNEW)
         CALL WEIGHT(N,TETAF,PNEW,PTIMEP,PTNEW)
         ITER=ITER+1
         GO TO 200
      END IF
C        
C  cases (ii), (iii), (iv): exit iteration loop
C        
      ITRTOT=ITRTOT+ITER
      CALL TIM(CPUNL1,2)
      CPUNL=CPUNL+CPUNL1
      WRITE(IOUT2,1020)
      DO II=1,ITER
         WRITE(IOUT2,1030) II,PL2V(II),PIKMXV(II),ITUMAX(II),PCURRV(II),
     1                     PPREVV(II),FL2V(II),FINFV(II)
      END DO
      IF (     ( (.NOT. SFCHEK) .AND. (.NOT. LSFAIL) .AND.
     1           (.NOT. ERRGMX) .AND. NORMCV )
     2    .OR. (       SFCHEK .AND. (.NOT. LSFAIL) .AND.
     3           (.NOT. ERRGMX) .AND. NORMCV .AND. KSFZER ) ) THEN
C        
C  case (ii): convergence achieved
C        
         NOBACK=.FALSE.  
         WRITE(IOUT2,1090) ITER
         WRITE(ITERM,1090) ITER
      ELSE
         WRITE(IOUT2,1040) ITER
         WRITE(ITERM,1040) ITER
         IF (DTGMIN) THEN
C    
C  case (iii): convergence failed - back-step
C  
            WRITE(IOUT2,1095)
            BCKSTP=.TRUE.
            NOBACK=.FALSE.  
         ELSE                
C    
C  case (iv): convergence failed - no back-stepping possible
C
            NOBACK=.TRUE.    
            WRITE(IOUT2,1100)
         END IF
      END IF
C
      RETURN
 1020 FORMAT(/,'                     NONLINEAR CONVERGENCE BEHAVIOR ',
     1       /,' iter- convergence error norms  node    PNEW at',
     2         '    POLD at  residual error norms',
     3       /,' ation         PL2      PIKMAX IKMAX      IKMAX',
     4         '      IKMAX        FL2       FINF')
 1030 FORMAT(I6,2(1PE12.4),I6,2(1PE11.2),2(1PE11.3)) 
 1040 FORMAT(1X,'CONVERGENCE NOT ACHIEVED IN ',I4,' ITERATIONS')
 1050 FORMAT(/,'      LINEAR SOLVER CONVERGENCE BEHAVIOR',
     1       /,' iter       residual  real residual')
 1090 FORMAT(1X,'CONVERGENCE ACHIEVED IN ',I4,' ITERATIONS')  
 1095 FORMAT(1X,'>>> BACK-STEPPING...')
 1100 FORMAT(1X,'>>> NO BACK-STEPPING POSSIBLE -- SIMULATION WILL ',
     1          'TERMINATE')
      END
