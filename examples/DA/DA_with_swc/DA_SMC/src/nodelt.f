C
C**************************  NODELT ************************************
C
C  average nodal input array VNOD to obtain values at each 
C  element (output array VELT)
C
C***********************************************************************
C
      SUBROUTINE NODELT(NT,TETRA,VNOD,VELT)
C
      IMPLICIT  NONE
      INTEGER   J,IEL,INOD
      INTEGER   NT
      INTEGER   TETRA(5,*)
      REAL*8    R4
      REAL*8    VNOD(*),VELT(*)
C
      R4=1.0D0/4.0D0
      DO IEL=1,NT
         VELT(IEL)=0.0D0
         DO J=1,4
            INOD=TETRA(J,IEL)
            VELT(IEL)=VELT(IEL) + VNOD(INOD)
         END DO
         VELT(IEL)=VELT(IEL)*R4
      END DO
C
      RETURN
      END
