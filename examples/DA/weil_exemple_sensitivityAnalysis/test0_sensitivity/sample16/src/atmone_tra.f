C
C**************************  ATMONE_TRA ************************************
C
C  read and initialize atmospheric BC FOR TRANSPORT 
C   for first time step (times 0 and DELTAT),
C
C***********************************************************************
C
      SUBROUTINE ATMONE_TRA(NNOD,CAUSPATM,CAUTIATM,IETOCAU,TIME,DELTAT,
     1     ATMCONC,CONCOLD,CONCTIM,CONCINP)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J,INOD
      INTEGER   NNOD,CAUSPATM,CAUTIATM,IETOCAU
      REAL*8    TIMEIN,SLOPE
      REAL*8    TIME,DELTAT
      REAL*8    ATMCONC(*),CONCOLD(*)
      REAL*8    CONCTIM(*),CONCINP(3,*)
      INCLUDE  'IOUNITS.H'
C
C INITIALISATION
      CAUTIATM=0
      CONCTIM(1)=0.0D0
      CONCTIM(2)=0.0D0
      DO I=1,NNOD
         CONCINP(1,I)=0.0D0
         CONCINP(2,I)=0.0D0
         ATMCONC(I)=0.0D0
         CONCOLD(I)=0.0D0
      END DO
C
      READ(IIN63,*,END=800) CAUSPATM,IETOCAU
C
      IF (CAUSPATM .EQ. 9999) THEN
         CAUTIATM=1
         GO TO 800
      END IF
C
      READ(IIN63,*,END=600) TIMEIN
C
      IF (DELTAT .GE. 1.0D+15) TIMEIN=0.0D0
      CONCTIM(3)=TIMEIN
C
      IF (CAUSPATM .EQ. 0) THEN
         READ(IIN63,*) (CONCINP(3,I),I=1,NNOD)
      ELSE
         READ(IIN63,*) CONCINP(3,1)
         DO I=2,NNOD
            CONCINP(3,I)=CONCINP(3,1)
         END DO
      END IF
C
      IF (CONCTIM(3) .LE. 0.0D0) THEN
         DO I=1,NNOD
            CONCOLD(I)=CONCINP(3,I)
         END DO
      END IF
C
      IF (DELTAT .GE. 1.0D+15) GO TO 300
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
 500  CONTINUE
      GO TO 800
       
C
C
  600 CAUTIATM=1
      GO TO 800
  700 CAUTIATM=1
      GO TO 300
C
  800 RETURN
 1100 FORMAT(/,5X,'HSPATM(0 SPAT. VAR. ATM BC, ELSE HOMO.) = ',I6)
      END
