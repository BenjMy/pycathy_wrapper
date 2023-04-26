C
C**************************  NEWTON ************************************
C
C  Newton iteration for the nonlinear equation
C
C***********************************************************************
C
      SUBROUTINE NEWTON(CPUVEC,N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN,
     1                  ISOLV,ITMXCG,TOLCG,TIME,DELTAT,TETAF,IVGHU,
     2                  NLKP,KSLOPE,TP,TETRA,TETJA,IA,JA,TOPOL,
     3                  PERMX,PERMY,PERMZ,Z,
     4                  SNODI,PNODI,CONTP,CONTQ,PRESC,Q,BI,CI,DI,
     5                  LMASS,VOLU,VOLUR,IVOL,RMAX,RMIN,ETAI,ETAE,
     6                  DETAI,DETAIE,SW,SWE,CKRW,CKRWE,DCKRW,DCKRWE,
     7                  POLD,PDIFF,PTOLD,PNEW,PTNEW,PTIMEP,
     8                  TNOTI,XT5,LHSP,COEF1,COEF2,COEF3,SCR,
     9                  IBOT,MINBOT,INSYM,RNSYM,NNOD,IFATM,
     A                  ATMACT,ATMOLD,LHSATM,NUMDIR,NODDIR,NSF,
     B                  NSFNUM,NSFNOD,SFEX,LHSSF,LSFAIL,KLSFAI,DUPUIT)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  K,IERSYM
      INTEGER  N,NT,NTRI,NTERM,NP,NQ,NITER,NITERT,ITLIN
      INTEGER  ISOLV,ITMXCG,IVGHU
      INTEGER  NLKP,KSLOPE,IBOT,MINBOT,NNOD,NUMDIR,NSF,KLSFAI,DUPUIT
      INTEGER  TP(*),TETRA(5,*),TETJA(4,4,*),IA(*),JA(*),TOPOL(*)
      INTEGER  CONTP(*),CONTQ(*),IVOL(*)
      INTEGER  INSYM(*),IFATM(*),NODDIR(*)
      INTEGER  NSFNUM(*),NSFNOD(NSFMAX,*),SFEX(NSFMAX,*)
      LOGICAL  LSFAIL
      REAL     CPUL1
      REAL     CPUVEC(*)
      REAL*8   TOLCG,TIME,DELTAT,TETAF,RMAX,RMIN
      REAL*8   PERMX(MAXSTR,*),PERMY(MAXSTR,*),PERMZ(MAXSTR,*),Z(*)
      REAL*8   SNODI(*),PNODI(*),PRESC(*)
      REAL*8   Q(*),BI(4,*),CI(4,*),DI(4,*)
      REAL*8   LMASS(4,4),VOLU(*),VOLUR(*)
      REAL*8   ETAI(*),ETAE(*)
      REAL*8   DETAI(*),DETAIE(*),SW(*),SWE(*)
      REAL*8   CKRW(*),CKRWE(*),DCKRW(*),DCKRWE(*)
      REAL*8   POLD(*),PDIFF(*),PTOLD(*),PNEW(*),PTNEW(*),PTIMEP(*)
      REAL*8   TNOTI(*),XT5(*),LHSP(*)
      REAL*8   COEF1(*),COEF2(*),COEF3(*),SCR(*),RNSYM(*)
      REAL*8   ATMACT(*),ATMOLD(*),LHSATM(*),LHSSF(NSFMAX,*)
      INCLUDE 'IOUNITS.H'
C
C  calculate soil moisture characteristics needed for Newton scheme
C
      CALL TIM(CPUL1,1)
      CALL NEWUNS(NLKP,N,NT,KSLOPE,IVGHU,NTRI,TP,TETRA,
     1            PTNEW,PTOLD,SNODI,PNODI,
     2            SW,CKRW,ETAI,DCKRW,DETAI,
     3            SWE,CKRWE,ETAE,DCKRWE,DETAIE)
      CALL TIM(CPUL1,2)
      CPUVEC(1)=CPUVEC(1)+CPUL1
C
C  initialize the arrays which will contain the global stiffness,
C  mass, and Jacobian (derivative terms) matrices
C
      CALL TIM(CPUL1,1)
      CALL INIT0R(NTERM,COEF1)
      CALL INIT0R(NTERM,COEF2)
      CALL INIT0R(NTERM,COEF3)
      CALL TIM(CPUL1,2)
      CPUVEC(2)=CPUVEC(2)+CPUL1
C
C  assemble global stiffness and mass matrices from the local
C  contributions of each element, and also the derivative term
C  components of the Jacobian
C
      CALL TIM(CPUL1,1)
      CALL ASSNEW(NT,NTRI,TETRA,TETJA,LMASS,COEF1,COEF2,
     1            COEF3,PERMX,PERMY,PERMZ,CKRWE,ETAE,
     2            PTNEW,PTIMEP,PNEW,DCKRW,DETAI,TETAF,DELTAT,
     3            VOLU,VOLUR,IVOL,BI,CI,DI)
      CALL TIM(CPUL1,2)
      CPUVEC(3)=CPUVEC(3)+CPUL1
C
C  assemble RHS vector, without contribution of the unsaturated
C  zone gravitational term and without the boundary conditions
C
      CALL TIM(CPUL1,1)
      CALL RHSNEW(N,JA,TOPOL,DELTAT,PTNEW,PNEW,
     1            PTIMEP,COEF1,COEF2,TNOTI)
      CALL TIM(CPUL1,2)
      CPUVEC(4)=CPUVEC(4)+CPUL1
C
C  assemble the global LHS system matrix (the Jacobian) from the
C  stiffness and mass matrices and the derivative term components
C  of the Jacobian
C
      CALL TIM(CPUL1,1)
      CALL CFMATN(NTERM,TETAF,DELTAT,COEF1,COEF2,COEF3)
      CALL TIM(CPUL1,2)
      CPUVEC(5)=CPUVEC(5)+CPUL1
C
C  add contribution of the unsaturated zone gravitational term
C  to the RHS vector
C
      CALL TIM(CPUL1,1)
      CALL RHSGRV(NT,NTRI,TETRA,TNOTI,DI,PERMZ,IVOL,CKRWE)
C
C  save a copy of RHS vector before imposing Dirichlet boundary
C  conditions (needed for back-calculation of fluxes)
C
      CALL VCOPYR(N,XT5,TNOTI)
C
C  impose boundary conditions, and
C  save diagonal elements of LHS system matrix corresponding to
C  Dirichlet nodes
C
      CALL BCNEW(NP,NQ,CONTP,CONTQ,LHSP,Q,NNOD,LHSATM,
     1           IFATM,ATMACT,ATMOLD,JA,TOPOL,COEF1,TNOTI,
     2           TETAF,RMAX,NUMDIR,NODDIR,
     3           NSF,NSFNUM,NSFNOD,SFEX,LHSSF)
      CALL TIM(CPUL1,2)
      CPUVEC(6)=CPUVEC(6)+CPUL1
C
C  solve the linear system of equations and calculate the
C  residual error in the linear solution
C 
      CALL TIM(CPUL1,1)
      CALL NSYSLV(ISOLV,IOUT2,N,NTERM,NUMDIR,NODDIR,ITMXCG,IBOT,MINBOT,
     1            MAXBOT,IERSYM,NITER,IA,JA,TOPOL,INSYM,TOLCG,RMIN,
     2            COEF1,COEF2,COEF3,SCR,RNSYM,PDIFF,TNOTI)
      IF (ISOLV .NE. 3) THEN
         ITLIN=ITLIN+NITER
         NITERT=NITERT+NITER
      END IF
C
C  set flag if the linear solver failed
C
      IF (ISOLV .EQ. 3) THEN
         IF (IERSYM .NE. 0) THEN
            LSFAIL=.TRUE.
            KLSFAI=KLSFAI + 1
         ELSE
            LSFAIL=.FALSE.
         END IF
      ELSE
         IF (NITER .GE. ITMXCG) THEN
            LSFAIL=.TRUE.
            KLSFAI=KLSFAI + 1
         ELSE
            LSFAIL=.FALSE.
         END IF
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
C  Dirichlet nodes may not be exactly equal to the presribed
C  values, due to inaccuracies and roundoff errors which could
C  arise from the way we imposed Dirichlet conditions (by multiplying
C  diagonal terms by a 'large' number).
C
      CALL SHLNEW(NP,TOPOL,JA,CONTP,PRESC,LHSP,COEF1,PNEW,POLD,
     1            NSF,NSFNUM,NSFNOD,SFEX,LHSSF,DUPUIT,Z,
     2            NNOD,IFATM,LHSATM)
      CALL TIM(CPUL1,2)
      CPUVEC(8)=CPUVEC(8)+CPUL1
C
      RETURN
      END
