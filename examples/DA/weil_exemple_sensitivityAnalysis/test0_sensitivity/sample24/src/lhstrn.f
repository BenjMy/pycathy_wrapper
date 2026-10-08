C
C**************************  LHSTRN ************************************
C
C  save diagonal elements of LHS system matrix corresponding to 
C  Dirichlet nodes
C
C***********************************************************************
C
      SUBROUTINE LHSTRN(NPC,NNPC,TOPOLC,JAC,LHSC,COEF1C)
C
      IMPLICIT  NONE
      INTEGER   K,J,IND
      INTEGER   NPC
      INTEGER   NNPC(*),TOPOLC(*),JAC(*)
      REAL*8    LHSC(*),COEF1C(*)
C
      DO K=1,NPC
         J=NNPC(K)
         IND=TOPOLC(J)
         LHSC(K)=COEF1C(IND)
      END DO
C
      RETURN
      END
