C
C**************************  SHLSYM ************************************
C
C  restore diagonal elements of LHS system matrix corresponding to
C  Dirichlet nodes. 
C  Also, set solution at Dirichlet nodes to the prescribed 
C  values. This is done since the solved solution at
C  Dirichlet nodes may not be exactly equal to the presribed
C  values, due to inaccuracies and roundoff errors which could 
C  arise from the way we imposed Dirichlet conditions (by multiplying
C  diagonal terms by a 'large' number).
C  (symmetric storage).
C
C***********************************************************************
C
      SUBROUTINE SHLSYM(NP,TOPOL,CONTP,PRESC,LHSP,COEF1,PNEW)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J,IND,INOD
      INTEGER   NP
      INTEGER   TOPOL(*),CONTP(*)
      REAL*8    PRESC(*),LHSP(*),COEF1(*),PNEW(*)
C
      DO I=1,NP
         J=CONTP(I)
         PNEW(J)=PRESC(I)
         IND=TOPOL(J)
         COEF1(IND)=LHSP(I)
      END DO
C
      RETURN
      END
