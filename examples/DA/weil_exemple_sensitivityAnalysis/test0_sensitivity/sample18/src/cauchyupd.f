C
C**************************  CAUCHYUPD ************************************
C
C procedure used to calculate the total Cauchy BC imposed at land surface 
C for subsurface transport. The total QCAUCHY BC is calculated from flow 
C results surface concentration and atmospheric cauchy BC for transport
C
C***********************************************************************
C
      SUBROUTINE CAUCHYUPD(NNOD,NQ_TRA,CONTQ_TRA,QCAUCHY,IFATM,IFCONC,
     1     ATMPOT,ATMACT,OVFLNOD,PNEW,CONCNOD,DELTAT,ATMCONC,
     2     ARENOD,CNNEW,IFATMP,PONDNOD,CONCNODUPD)
C
      IMPLICIT  NONE
      INCLUDE 'CATHY.H'
      INTEGER   I
      INTEGER   NQ_TRA,NNOD
      INTEGER   IFCONC(*),IFATM(*),IFATMP(*),CONTQ_TRA(*)
      REAL*8    QCAUCHY(*)
      REAL*8    ATMCONC(*)
      REAL*8    DELTAT,INFPOT,ZERO
      REAL*8    PNEW(*),ATMPOT(*),ATMACT(*)
      REAL*8    CONCNOD(*),ARENOD(*),OVFLNOD(*)
      REAL*8    CNNEW(*),CONCNODUPD(*),PONDNOD(*)
C      

      PARAMETER (ZERO = 0.0D0)
     
      NQ_TRA = 0
C     
      DO I=1,NNOD
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
C      rainfall, infiltration 
C-----------------------------------         
                  IF (IFATMP(I) .EQ. 0) THEN
                     NQ_TRA = NQ_TRA + 1
                     CONTQ_TRA(NQ_TRA) = I
                     QCAUCHY(NQ_TRA) = ATMACT(I) * ATMCONC(I) 
                     GO TO 500
                  ELSE
                     NQ_TRA = NQ_TRA + 1
                     CONTQ_TRA(NQ_TRA) = I
                     QCAUCHY(NQ_TRA) = ATMACT(I) * CONCNODUPD(I)   
                     GO TO 500
                  END IF
               ELSE
C-----------------------------------              
C     rainfall, exfiltration (==> return flow or seepage) 
C-----------------------------------              
C     ALWAYS DEFAULT BC ==> NEUMAN =0                  
                  GO TO 500
               END IF
            ELSE
               IF (ATMACT(I) .GE. ZERO) THEN
C     
C-----------------------------------
C     evaporation, infiltration 
C----------------------------------
                  NQ_TRA = NQ_TRA + 1
                  CONTQ_TRA(NQ_TRA) = I
                  QCAUCHY(NQ_TRA) = ATMACT(I) * CONCNODUPD(I) 
                  GO TO 500
               ELSE
C-----------------------------------
C     evaporation, exfiltration (==> return flow or seepage) 
C-----------------------------------
C     
C     ALWAYS DEFAULT BC ==> NEUMAN =0                  
                  GO TO 500
               END IF
            END IF
         END IF
C     
c-----------------------------------         
c    NO PONDING situation at T+1
c-----------------------------------  
         
C     
         IF (IFATM(I) .EQ. 0) THEN
            IF (ATMPOT(I) .GE. ZERO) THEN
               IF (ATMACT(I) .GE. ZERO) THEN
C     
c-----------------------------------         
C     rainfall,infiltration
c-----------------------------------
C     
                  IF (IFATMP(I) .EQ. 0) THEN
                     NQ_TRA = NQ_TRA + 1
                     CONTQ_TRA(NQ_TRA) = I 
                     QCAUCHY(NQ_TRA)=ATMPOT(I)*ATMCONC(I)
                     GO TO 500
                  ELSE
                     NQ_TRA = NQ_TRA + 1
                     CONTQ_TRA(NQ_TRA) = I
                     QCAUCHY(NQ_TRA)=(PONDNOD(I)*ARENOD(I)/DELTAT)*
     1                    CONCNOD(I) + ATMPOT(I)*ATMCONC(I)
                     GO TO 500
                  END IF
               ELSE
c-----------------------------------         
C     rainfall,exfiltration
c-----------------------------------                 
C     THEORETICALLY IMPOSSIBLE TO BE UNSATUREATED WITH RAINFALL 
C     AND EXFILTRATION                  
                  GO TO 500
               END IF
            ELSE
               IF (ATMACT(I) .GE. ZERO) THEN
c-----------------------------------
C     evaporation,infiltration
C-----------------------------------
                  NQ_TRA = NQ_TRA + 1
                  CONTQ_TRA(NQ_TRA) = I
                  QCAUCHY(NQ_TRA) = ATMACT(I)*CONCNOD(I)
               GO TO 500
            ELSE
c-----------------------------------
C     evaporation,exfiltration
C-----------------------------------               
               IF (IFATMP(I) .EQ. 0) THEN
                  NQ_TRA = NQ_TRA + 1
                  CONTQ_TRA(NQ_TRA) = I 
                  QCAUCHY(NQ_TRA)=0.
                  GO TO 500
               ELSE
C TO BE DETERMINED!!!!!!
                  GO TO 500
               END IF
            END IF
         END IF
      END IF
 500  CONTINUE         
cm    WRITE(222,9070)I,IFATM(I),IFATMP(I),ATMPOT(I),ATMACT(I),
cm   1            ATMCONC(I),CONCNODUPD(I)
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
c9070 FORMAT(I6,I6,I6,4(1PE12.5))
      END
