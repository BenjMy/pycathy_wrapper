C
C**************************  ATMONE ************************************
C
C  read and initialize atmospheric boundary condition parameters and 
C  arrays for first time step (times 0 and DELTAT), and check for
C  switching of atmospheric boundary conditions for this first time step
C
C***********************************************************************
C
      SUBROUTINE ATMONE(NNOD,HSPATM,HTIATM,IETO,TIME,DELTAT,
     1                  PONDH_MIN,IFATM,IFATMP,ARENOD,ATMPOT,ATMACT,
     2                  ATMOLD,ATMTIM,ATMINP,PNEW,PTIMEP,
     3                  NP,NQ,CONTP,CONTQ,NSF,NSFNUM,NSFNOD)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   I,J,INOD
      INTEGER   NNOD,HSPATM,HTIATM,IETO,NSF
      INTEGER   NP,NQ
      INTEGER   CONTP(*),CONTQ(*),IFATM(*),IFATMP(*)
      INTEGER   NSFNUM(*),NSFNOD(NSFMAX,*)
      REAL*8    TIMEIN,SLOPE,GASDEV
      REAL*8    TIME,DELTAT,PONDH_MIN,DSATM
      REAL*8    ARENOD(*),ATMPOT(*),ATMACT(*),ATMOLD(*)
      REAL*8    ATMTIM(*),ATMINP(3,*),PNEW(*),PTIMEP(*)
      INCLUDE  'IOUNITS.H'
      INCLUDE  'SOILCHAR.H'
C     
      HTIATM=0
      ATMTIM(1)=0.0D0
      ATMTIM(2)=0.0D0
      DO I=1,NNOD
         IFATM(I)=0
         IFATMP(I)=0
         ATMINP(1,I)=0.0D0
         ATMINP(2,I)=0.0D0
         ATMPOT(I)=0.0D0
         ATMACT(I)=0.0D0
         ATMOLD(I)=0.0D0
      END DO
      READ(IIN6,*,END=800) HSPATM,IETO
      IF (HSPATM .EQ. 9999) THEN
         HTIATM=1
         DO I=1,NNOD
            IFATM(I)=-1
            IFATMP(I)=-1
         END DO
         GO TO 800
      END IF
      WRITE(IOUT2,1100) HSPATM
      READ(IIN6,*,END=600) TIMEIN
      IF (DELTAT .GE. 1.0D+15) TIMEIN=0.0D0
      ATMTIM(3)=TIMEIN
      IF (HSPATM .EQ. 0) THEN
         READ(IIN6,*) (ATMINP(3,I),I=1,NNOD)
      ELSE
         READ(IIN6,*) ATMINP(3,1)
         DO I=2,NNOD
            ATMINP(3,I)=ATMINP(3,1)
         END DO
      END IF
      IF (ATMTIM(3) .LE. 0.0D0) THEN
         DO I=1,NNOD
            ATMOLD(I)=ATMINP(3,I)*ARENOD(I)
         END DO
      END IF
      IF (DELTAT .GE. 1.0D+15) GO TO 300
  200 IF (TIME .LE. ATMTIM(3)) GO TO 300
      ATMTIM(1)=ATMTIM(2)
      ATMTIM(2)=ATMTIM(3)
      DO I=1,NNOD
         ATMINP(1,I)=ATMINP(2,I)
         ATMINP(2,I)=ATMINP(3,I)
      END DO
      READ(IIN6,*,END=700) TIMEIN
      ATMTIM(3)=TIMEIN
      IF (HSPATM .EQ. 0) THEN
         READ(IIN6,*) (ATMINP(3,I),I=1,NNOD)
      ELSE
         READ(IIN6,*) ATMINP(3,1)
         DO I=2,NNOD
            ATMINP(3,I)=ATMINP(3,1)
         END DO
      END IF
      GO TO 200
  300 CONTINUE
      IF (ATMTIM(3) .GT. ATMTIM(2)) THEN
         DO I=1,NNOD
            SLOPE=(ATMINP(3,I) - ATMINP(2,I))/(ATMTIM(3) - ATMTIM(2))
            IF (IETO.NE.0) SLOPE=0.0D0
            ATMPOT(I)=(ATMINP(2,I) + SLOPE*(TIME - ATMTIM(2)))*ARENOD(I)
         END DO
      ELSE
         DO I=1,NNOD
            ATMPOT(I)=ATMINP(3,I)*ARENOD(I)
         END DO
      END IF
  400 DO I=1,NP
         INOD=CONTP(I)
         IF (INOD .LE. NNOD) THEN
            IFATM(INOD) =-1
            IFATMP(INOD)=-1
         END IF
      END DO
      DO I=1,NQ
         INOD=CONTQ(I)
         IF (INOD .LE. NNOD) THEN
            IFATM(INOD) =-1
            IFATMP(INOD)=-1
         END IF
      END DO
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            INOD=NSFNOD(I,J)
            IF (INOD .LE. NNOD) THEN
               IFATM(INOD) =-1
               IFATMP(INOD)=-1
            END IF
         END DO
      END DO
      DO 500 I=1,NNOD
         IF (IFATM(I) .EQ. -1) GO TO 500
         IF (PNEW(I) .GE. PONDH_MIN) THEN
            IFATM(I) = 2
            IFATMP(I) = 2
            GO TO 500
         END IF
         IF (PNEW(I) .GE. 0.0D0  .AND.  ATMPOT(I) .GT. 0.0D0) THEN
            IFATM(I)=1
         END IF
         IF (PTIMEP(I) .GE. 0.0D0  .AND. ATMOLD(I) .GT. 0.0D0) THEN
            IFATMP(I)=1
         END IF
         IF (PNEW(I) .LE. PMIN  .AND.  ATMPOT(I) .LT. 0.0D0) THEN
            PNEW(I)=PMIN
            IFATM(I)=1
         END IF
         IF (PTIMEP(I) .LE. PMIN  .AND.  ATMOLD(I) .LT. 0.0D0) THEN
            PTIMEP(I)=PMIN
            IFATMP(I)=1
         END IF
  500 CONTINUE
C
C  there are no back-calculated fluxes for first time step so we use 0.0
C
      DO I=1,NNOD
         IF (IFATM(I) .EQ. 0) THEN
            IF (ATMPOT(I).GE.0.0d0) THEN
                ATMACT(I)=ATMPOT(I)
            ELSE
                ATMACT(I)=(1.0d0-SCF)*ATMPOT(I)
            END IF
         ELSE
            ATMACT(I)=0.0D0
         END IF
         IF (IFATMP(I) .EQ. 1 .OR. IFATMP(I) .EQ. 2) THEN
            ATMOLD(I)=0.0D0
         END IF
      END DO
      GO TO 800
C
  600 HTIATM=1
      GO TO 400
  700 HTIATM=1
      GO TO 300
C
  800 RETURN
 1100 FORMAT(/,5X,'HSPATM(0 SPAT. VAR. ATM BC, ELSE HOMO.) = ',I6)
      END
