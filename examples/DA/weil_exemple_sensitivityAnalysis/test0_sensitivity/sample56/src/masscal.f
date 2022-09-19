C
C**************************  STORCAL ***********************************
C
C  calculate volume of water in the subsurface
C
C***********************************************************************
C
      SUBROUTINE MASSCAL(NT,MASST,CNEW_EL,SW_EL,PEL,VOLU)
C
      IMPLICIT  NONE
      INTEGER   I
      INTEGER   NT
      REAL*8    MASST
      REAL*8    CNEW_EL(*),SW_EL(*),PEL(*),VOLU(*)
C
      MASST = 0.0D0
      DO I=1,NT
         MASST = MASST + CNEW_EL(I)*SW_EL(I)*VOLU(I)*PEL(I)
      END DO
C
      RETURN
      END
