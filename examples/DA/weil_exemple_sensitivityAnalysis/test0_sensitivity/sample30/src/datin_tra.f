C
C**************************  DATIN_TRA  ************************************
C
C  initialization of parameters for transport > READ in the transp input file
C
C***************************************************************************
C
      SUBROUTINE DATIN_TRA(NSTR,NZONE,NTRI,NTRI3,NADV,NADVFL,LUMPC,
     1     CDIFFUS,FLAG_TEMP,FLAG_INTERP,FLAG_LIMITER,CFLINP,TETAC,
     2     DIFFUS,GAMMAS,IPARM,ALFAL,ALFAT,KD,LAMBDA,MIXPART,REACFLAG)
C     
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J
      INTEGER   NSTR,NZONE,NTRI
      INTEGER   NTRI3,NADV,NADVFL,LUMPC,CDIFFUS
      INTEGER   IPRTCG,IPREC,IMAX,ISOL,IEXIT
      INTEGER   FLAG_TEMP,FLAG_INTERP,FLAG_LIMITER
      INTEGER   IPARM(10)
      REAL*8    CFLINP,TETAC,DIFFUS,GAMMAS 
      REAL*8    ALFAL(MAXSTR,MAXZON),ALFAT(MAXSTR,MAXZON)
      REAL*8    KD(MAXSTR,MAXZON),LAMBDA(MAXSTR,MAXZON),MIXPART
      LOGICAL   REACFLAG
C
      INCLUDE 'IOUNITS.H'
C      
      REACFLAG=.FALSE.
      READ(IIN61,*) TETAC,LUMPC,CDIFFUS
      READ(IIN61,*) IPRTCG,IPREC,IMAX,ISOL,IEXIT
      DO I=1,10
         IPARM(I)=0
      END DO
      IPARM(1) = 98
      IPARM(2) = IPRTCG
      IPARM(3) = IPREC
      IPARM(6) = IMAX
      IPARM(7) = ISOL
      IPARM(10) = IEXIT
      READ(IIN61,*) FLAG_TEMP,FLAG_INTERP,FLAG_LIMITER
      READ(IIN61,*) CFLINP,NADV,NADVFL         
      READ(IIN61,*) MIXPART         
      DO I=1,NSTR
         DO J=1,NZONE
            READ(IIN61,*) ALFAL(I,J),ALFAT(I,J),KD(I,J),LAMBDA(I,J)
            IF (KD(I,J).NE.0) THEN
               REACFLAG=.TRUE.
            END IF
            IF (LAMBDA(I,J).NE.0) THEN
               REACFLAG=.TRUE.
            END IF
         END DO
      END DO
      READ(IIN61,*) DIFFUS, GAMMAS         
      NTRI3 = 3 * NTRI
C     
      
      RETURN
C

      END
