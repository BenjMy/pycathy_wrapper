C
C**************************  ICVHE  ************************************
C
C  calculate fully saturated vertical hydrostatic equilibrium IC's
C
C***********************************************************************
C
      SUBROUTINE ICVHE(NNOD,NSTR,N,IPRT1,IPOND,PTIMEP,PONDNOD,Z)
C
      IMPLICIT  NONE
      INTEGER   I,K,KK
      INTEGER   NNOD,NSTR,N,IPRT1,IPOND
      REAL*8    ZERO,ONE,TOTHED,MULT
      REAL*8    PTIMEP(N),PONDNOD(NNOD),Z(N)
      INCLUDE  'IOUNITS.H'
      PARAMETER (ZERO = 0.0, ONE = 1.0)
C
      IF (IPOND .EQ. 0) THEN
         MULT = ZERO
      ELSE
         MULT = ONE
      END IF
      DO I=1,NNOD
         TOTHED=Z(I) + MULT*PONDNOD(I)
         DO K=0,NSTR
            KK=K*NNOD + I
            PTIMEP(KK)=TOTHED - Z(KK)
         END DO
      END DO
      IF (IPRT1 .GE. 1) WRITE(IOUT2,1000) (K,PTIMEP(K),K=1,N)
C
      RETURN
 1000 FORMAT(/,1X,' INITIAL PRESSURE HEAD (HYDROSTATIC)',/,
     1       (4(I6,2X,1PE11.3)))
      END
