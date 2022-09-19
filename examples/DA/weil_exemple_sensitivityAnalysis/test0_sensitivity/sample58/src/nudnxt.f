C
C**************************  NUDNXT ************************************
C
C  input (if necessary) new nudging observation data and update
C  NUDC and NUDCTR
C
C***********************************************************************
C
      SUBROUTINE NUDNXT(NUDN,NUDT,NUDC,NUDCTR,WFLAG,
     1                  TIME,NUDTIM,NUDTAU,NUDVAL)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,J,JJ,K,KK
      INTEGER  NUDN,NUDT,NUDC,NUDCTR,WFLAG
      REAL*8   TIMOBS
      REAL*8   TIME
      REAL*8   NUDTIM(*),NUDTAU(MAXNUDT,*),NUDVAL(MAXNUDC,*)
      INCLUDE 'IOUNITS.H'
C
      IF (NUDN .EQ. 0) GO TO 900
      IF (NUDCTR .EQ. 0) GO TO 500
      JJ=NUDC
      DO K=1,NUDC
         KK=NUDCTR-NUDC+K
         IF (WFLAG .EQ. 0) THEN
            TIMOBS=NUDTIM(KK) + NUDTAU(KK,1)
         ELSE
            TIMOBS=NUDTIM(KK) + 3.0d0*DLOG(10.0d0)*NUDTAU(KK,1)
         END IF
         IF (TIME .GT. TIMOBS) THEN
            DO I=2,JJ
               DO J=1,NUDN
                  NUDVAL(I-1,J)=NUDVAL(I,J)
               END DO
            END DO
            JJ=JJ-1
         END IF
      END DO
      NUDC=JJ
  500 CONTINUE
      IF (NUDCTR .EQ. NUDT) GO TO 900
      IF (WFLAG .EQ. 0) THEN
         TIMOBS=NUDTIM(NUDCTR+1) - NUDTAU(NUDCTR+1,1)
      ELSE
         TIMOBS=NUDTIM(NUDCTR+1)
      END IF
      IF (TIME .LT. TIMOBS) GO TO 900
      NUDCTR=NUDCTR+1
      NUDC=NUDC+1
      READ(IIN50,*) (NUDVAL(NUDC,I),I=1,NUDN)
C
  900 RETURN
      END
