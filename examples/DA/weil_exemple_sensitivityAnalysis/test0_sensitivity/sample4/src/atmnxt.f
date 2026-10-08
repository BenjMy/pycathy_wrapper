C
C**************************  ATMNXT ************************************
C
C  input (if necessary) and interpolate atmospheric boundary 
C  conditions for next time level
C
C***********************************************************************
C
      SUBROUTINE ATMNXT(NNOD,HSPATM,HTIATM,IETO,TIME,IFATM,AREA,
     1                  ATMPOT,ATMACT,ATMTIM,ATMINP,DELTAT,
     2                  NP,NQ,CONTP,CONTQ,NSF,NSFNUM,NSFNOD,SCF)
C
      IMPLICIT  NONE
      INCLUDE   'CATHY.H'
      INTEGER   I,J,INOD
      INTEGER   NNOD,HSPATM,HTIATM,IETO
      INTEGER   NP,NQ,NSF
      INTEGER   IFATM(*)
      INTEGER   CONTP(*),CONTQ(*)
      INTEGER   NSFNUM(*),NSFNOD(NSFMAX,*)
      REAL*8    TIMEIN,SLOPE,DSATM,ALPHA,SCF
      REAL*8    TIME,GASDEV,DELTAT
      REAL*8    AREA(*),ATMPOT(*),ATMACT(*),ATMTIM(*),ATMINP(3,*)
      INCLUDE  'IOUNITS.H'
C
      IF (HTIATM .NE. 0) GO TO 800
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
            ATMPOT(I)=(ATMINP(2,I) + SLOPE*(TIME - ATMTIM(2)))*AREA(I)
            IF (IFATM(I) .EQ. 0) THEN
                IF (ATMPOT(I).GE.0.0d0) THEN
                    ATMACT(I)=ATMPOT(I)
                ELSE
                    ATMACT(I)=(1.0d0-SCF)*ATMPOT(I)
                END IF
            END IF
         END DO
      ELSE
         DO I=1,NNOD
            ATMPOT(I)=ATMINP(3,I)*AREA(I)
            IF (IFATM(I) .EQ. 0) THEN
                IF (ATMPOT(I).GE.0.0d0) THEN
                    ATMACT(I)=ATMPOT(I)
                ELSE
                    ATMACT(I)=(1.0d0-SCF)*ATMPOT(I)
                END IF
            END IF
         END DO
      END IF
      GO TO 800
C
  700 HTIATM=1
      GO TO 300
      
C
  800 CONTINUE
C
C  resets IFATM to exclude non atmospheric surface nodes
c  in case some boundary nodes have changed in BCNXT
C
      DO I=1,NP
         INOD=CONTP(I)
         IF (INOD .LE. NNOD) THEN
            IFATM(INOD) =-1
         END IF
      END DO
      DO I=1,NQ
         INOD=CONTQ(I)
         IF (INOD .LE. NNOD) THEN
            IFATM(INOD) =-1
         END IF
      END DO
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            INOD=NSFNOD(I,J)
            IF (INOD .LE. NNOD) THEN
               IFATM(INOD) =-1
            END IF
         END DO
      END DO
c
      RETURN
      END
