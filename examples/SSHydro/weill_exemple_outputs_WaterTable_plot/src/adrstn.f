C
C**************************  ADRSTN ************************************
C
C  Special handling (BC switching) for case when surface is `air dry'
C  and the potential atmospheric flux is > 0 (rainfall). We force a
C  'switch' to Neumann BC. We also switch to Neumann BC if the surface
C  is 'air dry' and the actual flux has become greater in magnitude
C  than the potential (evaporation) rate.
C
C***********************************************************************
C
      SUBROUTINE ADRSTN(NNOD,IFATM,ATMPOT,ATMACT,PNEW)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I
      INTEGER  NNOD
      INTEGER  IFATM(NNOD)
      REAL*8   ATMACT(NNOD),ATMPOT(NNOD),PNEW(NNOD)
      INCLUDE 'SOILCHAR.H'
C
      DO 500 I=1,NNOD
         IF (PNEW(I) .LE. PMIN) THEN
            IF (ATMPOT(I) .GT. 0.0D0 .OR. ATMACT(I) .LT. ATMPOT(I)) THEN
               IFATM(I)=0
               ATMACT(I)=ATMPOT(I)
               PNEW(I)=PMIN
            END IF
         END IF
  500 CONTINUE
C
      RETURN
      END
