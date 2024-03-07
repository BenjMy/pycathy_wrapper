C
C**************************  PONDUPD ***********************************
C
C  update surface BCs for FLOW3D to reflect ponding situation as
C  calculated in SURF_ROUTE, and update PONDING flag
C
C***********************************************************************
C
      SUBROUTINE PONDUPD(NNOD,PONDING,IFATM,PONDH_MIN,DELTAT,PONDNOD,
     1                   ARENOD,ATMPOT,ATMACT,PNEW)
C
      IMPLICIT NONE
      INTEGER  NNOD 
      INTEGER  I
      INTEGER  IFATM(*)
      LOGICAL  PONDING
      REAL*8   PONDH_MIN,DELTAT
      REAL*8   ZERO,ONE,DTR
      REAL*8   PONDNOD(*),ARENOD(*),ATMPOT(*),ATMACT(*)
      REAL*8   PNEW(*)
      PARAMETER (ZERO = 0.0D+00, ONE = 1.0D+00)
C    
      PONDING = .FALSE.
      DTR = ONE/DELTAT
      DO I = 1,NNOD
         IF(IFATM(I) .EQ. -1) GO TO 500
         IF (PONDNOD(I) .GE. PONDH_MIN) THEN
            IF (IFATM(I) .EQ. 1 .OR. IFATM(I) .EQ. 2) THEN
               PNEW(I) = PONDNOD(I)
               PONDING = .TRUE.
            ELSE IF (IFATM(I) .EQ. 0) THEN
                PONDING = .TRUE.
               ATMACT(I) = ATMPOT(I) + PONDNOD(I)*ARENOD(I)*DTR
!CMS                ATMACT(I) = ATMPOT(I)
            END IF
         END IF
!CMS         IF (I .EQ. 113 .OR. I .EQ. 116 .OR. I .EQ. 125 .OR. I
!CMS     1         .EQ. 128 .OR. I .EQ. 101 .OR. I .EQ. 102 .OR. I .EQ.
!CMS     2          103 .OR. I .EQ. 104 .OR. I .EQ. 124 .OR. I .EQ.
!CMS     3          125 .OR. I .EQ. 128 .OR. I .EQ. 129 .OR. I .EQ.
!CMS     4          114 .OR. I .EQ. 115 .OR. I .EQ. 116) THEN
!CMS                PONDING = .FALSE.
!CMS                ATMACT(I) = ATMPOT(I)
!CMS                PONDNOD(I) = ZERO
!CMS         END IF

 500     CONTINUE
      END DO
C
      RETURN  
      END
