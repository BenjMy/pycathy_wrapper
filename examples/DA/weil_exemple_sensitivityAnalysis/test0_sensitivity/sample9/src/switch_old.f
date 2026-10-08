C
C**************************  SWITCH_OLD ********************************
C
C  switching control of atmospheric boundary conditions
C  (old version, for the uncoupled, subsurface flow only case; should
C  be superceded by the new ADRSTN routine for both coupled and
C  uncoupled cases, but for the time being we keep the old version of
C  SWITCH active as well)
C
C***********************************************************************
C
      SUBROUTINE SWITCH_OLD(NNOD,IFATM,ATMACT,ATMPOT,PNEW)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  NNOD
      INTEGER  IFATM(*)
      REAL*8   ATMACT(*),ATMPOT(*),PNEW(*)
      INCLUDE 'SOILCHAR.H'
C
      DO 500 I=1,NNOD
         IF (IFATM(I) .EQ. -1) GO TO 500
         IF (IFATM(I) .EQ. 1   .AND.  PNEW(I) .GE. 0.0D0) THEN
            IF (ATMPOT(I).LT.0.0D0  .OR.  ATMACT(I).GT.ATMPOT(I)) THEN
               IFATM(I)=0
               ATMACT(I)=ATMPOT(I)
               GO TO 500
            END IF
         END IF
         IF (IFATM(I) .EQ. 1   .AND.  PNEW(I) .LE. PMIN) THEN
            IF (ATMPOT(I).GT.0.0D0  .OR.  ATMACT(I).LT.ATMPOT(I)) THEN
               IFATM(I)=0
               ATMACT(I)=ATMPOT(I)
               GO TO 500
            END IF
         END IF
         IF (IFATM(I) .EQ. 0   .AND.  PNEW(I) .GE. 0.0D0
     1                         .AND.  ATMPOT(I) .GE. 0.0D0) THEN
            IFATM(I)=1
            PNEW(I)=0.0D0
            GO TO 500
         END IF
         IF (IFATM(I) .EQ. 0   .AND.  PNEW(I) .LE. PMIN
     1                         .AND.  ATMPOT(I) .LT. 0.0D0) THEN
            IFATM(I)=1
            PNEW(I)=PMIN
            GO TO 500
         END IF
  500 CONTINUE
C
      RETURN
      END
