C
C**************************  SOURCE_TRASUR ************************************
C
C Procedure that calculates the source term for the surface transport equation 
C after the resolution of the subsurface transport equation. The objective is to 
C calculate TRAFLNOD that is the source term at each node for the surface transport 
C equation 
C
C***********************************************************************
C
      SUBROUTINE SOURCE_TRASUR(NNOD,IFATM,IFATMP,IFCONC,DELTAT,ATMPOT,
     1     ATMACT,PONDNOD,CONCNOD,CNNEW,ATMCONC,TRAFLNOD,ARENOD,
     2     CONCNODUPD,PONDNODP)
C
      IMPLICIT  NONE
      INCLUDE 'CATHY.H'
      INTEGER   I,J
      INTEGER   NNOD
      INTEGER   IFATM(*),IFATMP(*),IFCONC(*)
      REAL*8    ATMCONC(*)
      REAL*8    DELTAT,INFPOT,ZERO,FLUTRA,ETACT
      REAL*8    PONDNOD(*),ATMPOT(*),ATMACT(*),PONDNODP(*)
      REAL*8    CONCNOD(*),CNNEW(*),TRAFLNOD(*)
      REAL*8    ARENOD(*),CONCNODUPD(*)
C      
      PARAMETER (ZERO = 0.0)
C
      DO I=1,NNOD
      TRAFLNOD(I) = 0.
         IF (IFCONC(I) .EQ. -1)  THEN
            GO TO 500
         END IF
c-----------------------------------         
c     PONDING situation at T+1
c----------------------------------- 
         IF (IFATM(I) .EQ. 2) THEN
            IF (ATMPOT(I) .GE. ZERO) THEN
               IF (ATMACT(I) .GE. ZERO) THEN
C-----------------------------------
C     rainfall, infiltration 
C-----------------------------------
                     TRAFLNOD(I)= ATMPOT(I)*ATMCONC(I) - 
     1                   ATMACT(I)*CONCNODUPD(I) 
                     GO TO 500
               ELSE
C     
C-----------------------------------
C     rainfall, exfiltration (==> return flow or seepage) 
C-----------------------------------
                  TRAFLNOD(I)= ATMPOT(I)*ATMCONC(I) - 
     1                 ATMACT(I)*CNNEW(I) 
C
                  GO TO 500
               END IF
            ELSE
               IF (ATMACT(I) .GE. ZERO) THEN
C     
C-----------------------------------
C     evaporation, infiltration 
C----------------------------------
                  TRAFLNOD(I) = -1.*ATMACT(I)*CONCNODUPD(I)
                     if ( TRAFLNOD(I).GT.0.0d0) THEN
                        END IF
                  GO TO 500
               ELSE
C     
C-----------------------------------
C     evaporation, exfiltration (==> return flow or seepage) 
C-----------------------------------
C     
                  TRAFLNOD(I) = -1.*ATMACT(I)*CNNEW(I)
                       if ( TRAFLNOD(I).GT.0.0d0) THEN
                        END IF
                  GO TO 500
               END IF
            END IF
         END IF
C     
c-----------------------------------         
c     NO PONDING situation at T+1
c-----------------------------------
         IF (IFATM(I) .EQ. 0) THEN
            IF (ATMPOT(I) .GE. ZERO) THEN
               IF (ATMACT(I) .GE. ZERO) then 
c     
c-----------------------------------
C     rainfall,infiltration
c-----------------------------------
C     
                  IF (IFATMP(I) .EQ. 0) THEN
                     TRAFLNOD(I) = 0.
                     GO TO 500
                  ELSE
                      TRAFLNOD(I) = ATMPOT(I)*ATMCONC(I) - 
     1                   ATMACT(I)*CONCNOD(I) 
                     GO TO 500
                  END IF
               ELSE 
                  GO TO 500
               END IF
            ELSE
               IF (ATMACT(I) .GE. ZERO) THEN
                  
c-----------------------------------
C     evaporation,infiltration
C-----------------------------------
                  TRAFLNOD(I)=-1.*ATMACT(I)*CONCNOD(I)
            if ( TRAFLNOD(I).GT.0.0d0) THEN
                        END IF
                  GO TO 500
               ELSE
c-----------------------------------
C     evaporation,exfiltration
C-----------------------------------
                  IF (IFATMP(I) .EQ. 0) THEN
                     TRAFLNOD(I)=0.
                      if ( TRAFLNOD(I).GT.0.0d0) THEN
                        END IF
                     GO TO 500
                  ELSE
C TO BE DETERMINED!!!!!!!
                     GO TO 500
                  END IF           
               END IF
            END IF
         END IF
 500     CONTINUE
      END DO
C     
      RETURN
 9000 FORMAT(  'SWITCH anomaly: unsaturated node ',I6,
     1         ' (PL = ',1PE12.5,')',
     2       /,'  generated where IFATM was 2 and under ',
     3         'rainfall + exfiltration conditions ',
     4       /,'  (ie, return flow/seepage)',
     5       /,'  (TIME=',1PE8.2,', DELTAT=',1PE8.2,')')
 9020 FORMAT(  'SWITCH anomaly: unsaturated node ',I6,
     1         ' (PL = ',1PE12.5,')',
     2       /,'  generated where IFATM was 1 and under ',
     3         'rainfall + exfiltration conditions ',
     4       /,'  (ie, return flow/seepage)',
     5       /,'  (TIME=',1PE8.2,', DELTAT=',1PE8.2,')')
 9040 FORMAT(  'SWITCH anomaly: air dry node ',I6,' with ',
     1         'Dirichlet BC yields infiltration flux ',
     2       /,'  (ATMACT = ',1PE12.5,')',
     3       /,'  (TIME=',1PE8.2,', DELTAT=',1PE8.2,')')
 9060 FORMAT(  'SWITCH anomaly: Neumann node ',I6,' under ',
     1         'rainfall conditions becomes air dry ',
     2       /,'  (ATMPOT = ',1PE12.5,')',
     3       /,'  (TIME=',1PE8.2,', DELTAT=',1PE8.2,')')
      END
