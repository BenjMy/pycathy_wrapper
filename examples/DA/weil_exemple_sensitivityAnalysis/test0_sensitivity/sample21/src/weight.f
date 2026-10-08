C
C**************************  WEIGHT ************************************
C
C  weighting parameter TETA is used to calculate the weighted
C  values VWGHT from the values at the current (VCURR) and
C  previous (VPREV) time or iteration levels
C
C***********************************************************************
C
      SUBROUTINE WEIGHT(N,TETA,VCURR,VPREV,VWGHT)
C
      IMPLICIT  NONE
      INTEGER   K
      INTEGER   N
      REAL*8    TETA1
      REAL*8    TETA
      REAL*8    VCURR(*),VPREV(*),VWGHT(*)
C
      TETA1=1.0D0 - TETA
      DO K=1,N
         VWGHT(K)=TETA*VCURR(K) + TETA1*VPREV(K)
      END DO
C
      RETURN
      END
