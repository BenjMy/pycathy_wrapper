C
C**************************  ATMNXT_TRA ************************************
C
C  input (if necessary) and interpolate atmospheric boundary 
C  conditionsFOR TRANSPORT for next time level
C
C***********************************************************************
C
      SUBROUTINE ATMNXT_TRA(NNOD,CAUSPATM,CAUTIATM,IETOCAU,TIME,
     1                  ATMCONC,CONCTIM,CONCINP,DELTAT)
C
      IMPLICIT  NONE
      INCLUDE   'CATHY.H'
      INTEGER   I,J
      INTEGER   NNOD,CAUSPATM,CAUTIATM,IETOCAU
      REAL*8    TIMEIN,SLOPE
      REAL*8    CONCTIM(*),CONCINP(3,*)
      REAL*8    ATMCONC(*)
      REAL*8    TIME,DELTAT
      INCLUDE  'IOUNITS.H'

C
      IF (CAUTIATM .NE. 0) GO TO 800
  200 IF (TIME .LE. CONCTIM(3)) GO TO 300
      CONCTIM(1)=CONCTIM(2)
      CONCTIM(2)=CONCTIM(3)
      DO I=1,NNOD
         CONCINP(1,I)=CONCINP(2,I)
         CONCINP(2,I)=CONCINP(3,I)
      END DO
      READ(IIN63,*,END=700) TIMEIN
      CONCTIM(3)=TIMEIN
      IF (CAUSPATM .EQ. 0) THEN
         READ(IIN63,*) (CONCINP(3,I),I=1,NNOD)
      ELSE
         READ(IIN63,*) CONCINP(3,1)
         DO I=2,NNOD
            CONCINP(3,I)=CONCINP(3,1)
         END DO
      END IF
      GO TO 200
  300 IF (CONCTIM(3) .GT. CONCTIM(2)) THEN
         DO I=1,NNOD
            SLOPE=(CONCINP(3,I)- CONCINP(2,I))/(CONCTIM(3) - CONCTIM(2))
            IF (IETOCAU.NE.0) SLOPE=0.0D0
            ATMCONC(I)=(CONCINP(2,I) + SLOPE*(TIME - CONCTIM(2)))
         END DO
      ELSE
         DO I=1,NNOD
            ATMCONC(I)=CONCINP(3,I)
         END DO
      END IF
      GO TO 800
C     
  700 CAUTIATM=1
      GO TO 300
      
C
  800 RETURN
      END
