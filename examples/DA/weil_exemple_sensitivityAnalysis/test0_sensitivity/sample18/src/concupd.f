C
C**************************  CONCUPD ***********************************
C
C  update surface CONCENTRATION FOR COUPLING SURFACE AND SUBSURFACE TRANSPORT
C
C***********************************************************************
C
      SUBROUTINE CONCUPD(NNOD,IFATM,IFATMP,DELTAT,PONDNOD,
     1     ARENOD,ATMPOT,ATMACT,CONCNOD,ATMCONC,
     2     CONCNODUPD)
C
      IMPLICIT NONE
      INTEGER  NNOD,I
      INTEGER  IFATM(*),IFATMP(*)
      REAL*8   DELTAT
      REAL*8   PONDNOD(*),ARENOD(*),ATMPOT(*),ATMACT(*)
      REAL*8   CONCNOD(*),ATMCONC(*),CONCNODUPD(*)
      REAL*8   VOLTOT(NNOD),MASSURF(NNOD),MASSPOT(NNOD)
      INCLUDE  'IOUNITS.H'
C    
      
      DO I = 1,NNOD
         CONCNODUPD(I) = CONCNOD(I)
         IF (IFATMP(I) .EQ. 2 .AND. IFATM(I) .EQ. 2) THEN
C     RAINFALL/INFILTRATION CASE         
            IF ((ATMPOT(I) .GE. 0.) .AND. (ATMACT(I) .GE. 0.)) THEN
               VOLTOT(I) = PONDNOD(I)*ARENOD(I) + ATMPOT(I)*DELTAT
               MASSURF(I) =  PONDNOD(I)*ARENOD(I)*CONCNOD(I)
               MASSPOT(I) = ATMPOT(I)*DELTAT*ATMCONC(I)
C     
               CONCNODUPD(I) = (MASSURF(I)+MASSPOT(I))/VOLTOT(I)
C     EVAPORATION/INFILTRATION CASE         
            ELSEIF ((ATMPOT(I) .LT. 0.) .AND. (ATMACT(I) .GT. 0.)) THEN
               VOLTOT(I) = PONDNOD(I)*ARENOD(I) + ATMPOT(I)*DELTAT
               MASSURF(I) =  PONDNOD(I)*ARENOD(I)*CONCNOD(I)
C     
               CONCNODUPD(I) = MASSURF(I)/VOLTOT(I)
            END IF           
         END IF
      END DO      
C     
      RETURN  
      END
