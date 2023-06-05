C
C**************************  ICVDWT ************************************
C
C  calculate partially saturated vertical hydrostatic equilibrium IC's
C
C***********************************************************************
C
      SUBROUTINE ICVDWT(NNOD,NSTR,N,IPRT1,IPOND,WTPOSITION,
     1                  PTIMEP,PONDNOD,Z)
C
      IMPLICIT  NONE
      INTEGER   I,K,KK
      INTEGER   NNOD,NSTR,N,IPRT1,IPOND
      REAL*8    WTPOSITION
      REAL*8    ZERO,ONE
      REAL*8    PTIMEP(N),PONDNOD(NNOD),Z(N)
      INCLUDE  'IOUNITS.H'
      PARAMETER (ZERO = 0.0, ONE = 1.0)
C
      DO I=1,NNOD
         DO K=0,NSTR
            KK=K*NNOD + I
            PTIMEP(KK)= Z(I) - Z(KK) - WTPOSITION
         END DO
      END DO
      IF (IPRT1 .GE. 1) WRITE(IOUT2,1000) (K,PTIMEP(K),K=1,N)
C
      RETURN
 1000 FORMAT(/,1X,' INITIAL PRESSURE HEAD (HYDROSTATIC, PART. SAT.)',/,
     1       (4(I6,2X,1PE11.3)))
      END
