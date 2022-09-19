C
C**************************  PICARD ************************************
C
C  Picard iteration for the nonlinear equation
C
C***********************************************************************
C
      SUBROUTINE PICARD(CPUVEC,N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN,
     1                  ITMXCG,TOLCG,TIME,DELTAT,TETAF,IPEAT,IVGHU,
     2                  NLKP,NUDN,NUDC,NUDCTR,NUDG,NUDFLAG,WFLAG,
     3                  KSLOPE,TP,TETRA,TETJA,JA,TOPOL,
     4                  PERMX,PERMY,PERMZ,SNODI,PNODI,PEL,CONTP,CONTQ,
     5                  PORE,INDE,PRESC,Q,AI,BI,CI,DI,LMASS,X,Y,Z,
     6                  INDE0,VOLNOD,VOLU,VOLUR,IVOL,RMAX,RMIN,
     7                  ETAI,ETAE,ET1,ET2,ET1E,SW,SWE,CKRW,CKRWE,
     C                  SENODI,SEELT,SWNEW,SWTIMEP,
     8                  POLD,PDIFF,PTOLD,PNEW,DEF,
     9                  PTNEW,PTIMEP,TNOTI,XT5,LHSP,COEF1,COEF2,COEF4,
     A                  SCR,NNOD,IFATM,ATMACT,ATMOLD,LHSATM,NUMDIR,
     B                  NODDIR,NSF,NSFNUM,NSFNOD,SFEX,LHSSF,DUPUIT,
     C                  LSFAIL,KLSFAI,
     D                  NUDTET,NUDTIM,NUDRXY,NUDRZ,NUDX,NUDY,NUDZ,
     E                  NUDEPS,NUDDIF,NUDSMC,NUDNOD,NUDVAL,NUDTAU,
     F                  QTRANIE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  K,I
      INTEGER  N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN
      INTEGER  ITMXCG,IPEAT,IVGHU
      INTEGER  NLKP,NUDN,NUDC,NUDCTR,NUDFLAG,WFLAG,KSLOPE
      INTEGER  NNOD,NUMDIR,NSF,DUPUIT,KLSFAI
      INTEGER  TP(*),TETRA(5,*),TETJA(4,4,*),JA(*),TOPOL(*)
      INTEGER  CONTP(*),CONTQ(*),IVOL(*),IFATM(*),NODDIR(*)
      INTEGER  NSFNUM(*),NSFNOD(NSFMAX,*),SFEX(NSFMAX,*),NUDTET(*)
      LOGICAL  LSFAIL
      REAL     CPUL1
      REAL     CPUVEC(*)
      REAL*8   TOLCG,TIME,DELTAT,TETAF,NUDG,RMAX,RMIN
      REAL*8   PERMX(MAXSTR,*),PERMY(MAXSTR,*),PERMZ(MAXSTR,*)
      REAL*8   SNODI(*),PNODI(*),PEL(*),PORE(*),INDE(*),PRESC(*)
      REAL*8   Q(*),AI(4,*),BI(4,*),CI(4,*),DI(4,*)
      REAL*8   LMASS(4,4),X(*),Y(*),Z(*)
      REAL*8   INDE0(*),VOLNOD(*),VOLU(*),VOLUR(*)
      REAL*8   ETAI(*),ETAE(*),SW(*),SWE(*),CKRW(*),CKRWE(*)
      REAL*8   ET1(*),ET2(*),ET1E(*)
      REAL*8   SENODI(*),SEELT(*)
      REAL*8   POLD(*),PDIFF(*),PTOLD(*),PNEW(*),DEF(*)
      REAL*8   PTNEW(*),PTIMEP(*),QTRANIE(NMAX)
c CARLOTTA SWNEW saturation at current iteration
c          SWTIMEP at the previous nonlinear iteration
      REAL*8   SWNEW(NMAX),SWTIMEP(NMAX)
c CARLOTTA
      REAL*8   TNOTI(*),XT5(*),LHSP(*)
      REAL*8   COEF1(*),COEF2(*),COEF4(*),COEF5(MAXTRM)
      REAL*8   SCR(*),ATMACT(*),ATMOLD(*),LHSATM(*),LHSSF(NSFMAX,*)
      REAL*8   NUDTIM(*),NUDRXY(*),NUDRZ(*)
      REAL*8   NUDX(*),NUDY(*),NUDZ(*)
      REAL*8   NUDEPS(*),NUDDIF(*),NUDSMC(*),NUDNOD(*)
      REAL*8   NUDVAL(MAXNUDC,*),NUDTAU(MAXNUDT,*)
      INCLUDE 'IOUNITS.H'
C
C  calculate soil moisture characteristics needed for Picard scheme
C  (and also SW which is needed for nudging)
C
      CALL TIM(CPUL1,1)
      IF (IPEAT .EQ. 1) then
          CALL PEATCH(NLKP,N,NNOD,NT,NTRI,IVGHU,TETRA,TP,
     1                PTNEW,Z,DEF,INDE,INDE0,PORE,
     2                SNODI,PNODI,SW,CKRW,SENODI,ETAI,
     3                CKRWE,SEELT,ETAE)
      ELSE
c CARLOTTA lo sto facendo solo in questo caso
          CALL PICUNS(NLKP,N,NT,KSLOPE,IVGHU,NTRI,TP,TETRA,
     1                PTNEW,PTOLD,PNEW,PTIMEP,SWNEW,SWTIMEP,
     2                SNODI,PNODI,SW,CKRW,ETAI,
     3                ET1,ET2,ET1E,SWE,CKRWE,ETAE,PEL)
      END IF
 
      CALL TIM(CPUL1,2)
      CPUVEC(1)=CPUVEC(1)+CPUL1
C
C  initialize the arrays which will contain the global stiffness
C  and mass matrices
C
      CALL TIM(CPUL1,1)
      CALL INIT0R(NTERM,COEF1)
      CALL INIT0R(NTERM,COEF2)
      CALL INIT0R(NTERM,COEF4)
      CALL INIT0R(NTERM,COEF5)
      CALL TIM(CPUL1,2)
      CPUVEC(2)=CPUVEC(2)+CPUL1
C
C  assemble global stiffness and mass matrices from the local
C  contributions of each element
C
      CALL TIM(CPUL1,1)
      CALL ASSPIC(NT,NTRI,TETRA,TETJA,LMASS,COEF1,COEF2,
     1            COEF4,ET1E,PEL,
     1            PERMX,PERMY,PERMZ,CKRWE,ETAE,
     2            VOLU,VOLUR,BI,CI,DI)
      CALL TIM(CPUL1,2)
      CPUVEC(3)=CPUVEC(3)+CPUL1
C
C  assemble RHS vector, without contributions of the unsaturated
C  zone gravitational term, the boundary conditions, and the
C  nudging term
C
      CALL TIM(CPUL1,1)
      CALL RHSPIC(N,JA,TOPOL,DELTAT,PTNEW,PNEW,PTIMEP,
     1            SWNEW,SWTIMEP,COEF1,COEF2,COEF4,TNOTI)
      CALL TIM(CPUL1,2)
      CPUVEC(4)=CPUVEC(4)+CPUL1
C
C  assemble the global LHS system matrix from the stiffness and mass
C  matrices
C
      CALL TIM(CPUL1,1)
      CALL CFMATP(N,NTERM,TETAF,DELTAT,COEF1,COEF2,COEF4,COEF5,
     1            TOPOL,ET2)
      CALL TIM(CPUL1,2)
      CPUVEC(5)=CPUVEC(5)+CPUL1
C
C  add contribution of the unsaturated zone gravitational term
C  to the RHS vector
C
      CALL TIM(CPUL1,1)
      CALL RHSGRV(NT,NTRI,TETRA,TNOTI,DI,PERMZ,IVOL,CKRWE)
C
C  add nudging contribution to the RHS vector
C
      IF (IPEAT .EQ. 0) THEN
         CALL NUDPIC(N,NT,NUDN,NUDC,NUDCTR,TIME,NUDG,NUDFLAG,WFLAG,
     1               TETRA,NUDTET,IVOL,PTNEW,
     2               X,Y,Z,VOLUR,VOLU,SW,PNODI,TNOTI,AI,BI,CI,DI,
     3               NUDTIM,NUDRXY,NUDRZ,NUDX,NUDY,NUDZ,NUDEPS,
     4               NUDDIF,NUDSMC,NUDNOD,NUDVAL,NUDTAU)
      ELSE
         CALL NUDPIC(N,NT,NUDN,NUDC,NUDCTR,TIME,NUDG,NUDFLAG,WFLAG,
     1               TETRA,NUDTET,IVOL,PTNEW,
     2               X,Y,Z,VOLUR,VOLU,SW,PORE,TNOTI,AI,BI,CI,DI,
     3               NUDTIM,NUDRXY,NUDRZ,NUDX,NUDY,NUDZ,NUDEPS,
     4               NUDDIF,NUDSMC,NUDNOD,NUDVAL,NUDTAU)
      END IF
C
C  save a copy of RHS vector before imposing Dirichlet boundary
C  conditions (needed for back-calculation of fluxes)
C
      CALL VCOPYR(N,XT5,TNOTI)
C
C  impose boundary conditions, and save diagonal elements of
C  LHS system matrix corresponding to Dirichlet nodes
C
      CALL BCPIC(NP,NQ,CONTP,CONTQ,LHSP,Q,NNOD,N,
     1           LHSATM,IFATM,ATMACT,ATMOLD,TOPOL,COEF1,TNOTI,
     2           TETAF,RMAX,NUMDIR,NODDIR,
     3           NSF,NSFNUM,NSFNOD,SFEX,LHSSF,QTRANIE)
      CALL TIM(CPUL1,2)
      CPUVEC(6)=CPUVEC(6)+CPUL1
C
C  solve the linear system of equations and calculate the
C  residual error in the linear solution
C 
      CALL TIM(CPUL1,1)
      CALL SYMSLV(IOUT2,N,NTERM,NUMDIR,NITER,ITMXCG,TOLCG,RMIN,
     1            NODDIR,JA,TOPOL,PDIFF,TNOTI,COEF1,COEF2,SCR)
      ITLIN=ITLIN+NITER
      NITERT=NITERT+NITER
C
C  set flag if the linear solver failed
C
      IF (NITER .GE. ITMXCG) THEN
         LSFAIL=.TRUE.
         KLSFAI=KLSFAI + 1
      ELSE
         LSFAIL=.FALSE.
      END IF
      CALL TIM(CPUL1,2)
      CPUVEC(7)=CPUVEC(7)+CPUL1
C
C  extract pressure head solution PNEW from the difference
C  solution PDIFF
C
      CALL TIM(CPUL1,1)
      DO K=1,N
         PNEW(K)=PNEW(K)+PDIFF(K)
      END DO
C
C  restore diagonal elements of LHS system matrix corresponding to
C  Dirichlet nodes.
C  Also, set solution at Dirichlet nodes to the prescribed
C  values. This is done since the solved solution at
C  Dirichlet nodes may not be exactly equal to the prescribed
C  values, due to inaccuracies and roundoff errors which could
C  arise from the way we imposed Dirichlet conditions (by multiplying
C  diagonal terms by a 'large' number).
C
      CALL SHLPIC(NP,TOPOL,CONTP,PRESC,LHSP,COEF1,PNEW,POLD,
     1            NSF,NSFNUM,NSFNOD,SFEX,LHSSF,DUPUIT,Z,
     2            NNOD,IFATM,LHSATM)
      CALL TIM(CPUL1,2)
      CPUVEC(8)=CPUVEC(8)+CPUL1
C
      RETURN
      END
