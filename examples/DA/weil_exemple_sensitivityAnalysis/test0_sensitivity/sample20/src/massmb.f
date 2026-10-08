C
C**************************  MASSMB ************************************
C
C  calculate mass of change in storage
C  between the current time level and the previous time level.
C  DSMASS > 0 for net increase in mass.
C
C***********************************************************************
C
      SUBROUTINE MASSMB(NT,DSMASS,CNEW_EL,COLDOLD,SW_EL,
     1            SWP_EL,PEL,VOLU)
C
      IMPLICIT NONE
      INTEGER  I
      INTEGER  NT
      REAL*8   DSMASS
      REAL*8   CNEW_EL(*),COLDOLD(*),SW_EL(*),SWP_EL(*)
      REAL*8   PEL(*),VOLU(*)
C
      DSMASS=0.0D0
      DO I=1,NT
         DSMASS=DSMASS + VOLU(I)*PEL(I)*(CNEW_EL(I)*SW_EL(I)
     1                 -COLDOLD(I)*SWP_EL(I))
      END DO
C
      RETURN
      END
