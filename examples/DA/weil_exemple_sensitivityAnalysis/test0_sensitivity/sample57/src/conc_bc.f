C
C**************************  CONC_BC ***********************************
C
C  update surface CONCENTRATION FOR COUPLING SURFACE AND SUBSURFACE TRANSPORT
C
C***********************************************************************
C
      SUBROUTINE CONC_BC(NNOD,IFATM,IFATMP,DELTAT,PONDNOD,
     1     PONDNODP,ARENOD,ATMPOT,ATMACT,CONCNOD,CONCNOD_OLD,ATMCONC,
     2     CNNEW,ANP_TRA,CONTP_TRA,PRESC_TRA,CONCNODBC,CONCNODUPD)
C
      IMPLICIT NONE
      INTEGER  NNOD,I
      INTEGER  IFATM(*),IFATMP(*),ANP_TRA,CONTP_TRA(*)
      REAL*8   DELTAT
      REAL*8   PONDNOD(*),PONDNODP(*),ARENOD(*),ATMPOT(*),ATMACT(*)
      REAL*8   CONCNOD(*),ATMCONC(*),CNNEW(*),PRESC_TRA(*)
      REAL*8   CONCNOD_OLD(*),CONCNODBC(*),CONCNODUPD(*)
      REAL*8   VOLTOT(NNOD),MASSURF(NNOD),MASSPOT(NNOD)
C    
C      
      DO I = 1,NNOD
         CONCNODBC(I) = CONCNOD(I)
         CONCNODUPD(I) = CONCNOD(I)
C -------------------------------
C  PONDING T
C -------------------------------
         IF ((IFATM(I) .EQ. 2).AND.(PONDNOD(I).NE.0.0d0)) THEN
C     RAINFALL / INFILTRATION CASE         
          IF ((ATMPOT(I) .GE. 0.0d0) .AND. (ATMACT(I) .GE. 0.0d0)) THEN
            
               VOLTOT(I) = PONDNOD(I) * ARENOD(I) + ATMPOT(I) * DELTAT
c            
               MASSURF(I) = PONDNOD(I) * ARENOD(I) * CONCNOD(I)
               MASSPOT(I) = ATMPOT(I) * DELTAT * ATMCONC(I)
C
               CONCNODBC(I) = (MASSURF(I) + MASSPOT(I)) / VOLTOT(I)
               CONCNODUPD(I) = CONCNODBC(I)   
               GO TO 500
          ELSEIF ((ATMPOT(I).GE.0.0d0).AND.(ATMACT(I).LT.0.0d0)) THEN
               VOLTOT(I) = PONDNOD(I) * ARENOD(I) + ATMPOT(I) * DELTAT -
     1              ATMACT(I) * DELTAT
c
               MASSURF(I) = PONDNOD(I) * ARENOD(I) * CONCNOD(I)
               MASSPOT(I) = ATMPOT(I) * DELTAT * ATMCONC(I)
     1                     - ATMACT(I) * CNNEW(I) * DELTAT
c
                CONCNODUPD(I) = (MASSURF(I) + MASSPOT(I)) / VOLTOT(I)    
                CONCNODBC(I) = CNNEW(I)
               GO TO 500
C     EVAPORATION / EXFILTRATION CASE         
           ELSEIF ((ATMPOT(I).LT.0.0d0).AND.(ATMACT(I).LT.0.0d0)) THEN
             write(*,*) 'exfiltration 2 !!!!!',I
               VOLTOT(I) = PONDNOD(I) * ARENOD(I) + ATMPOT(I) * DELTAT -
     1              ATMACT(I) * DELTAT
c
               MASSURF(I) = PONDNOD(I) * ARENOD(I) * CONCNOD(I)
               MASSPOT(I) = - ATMACT(I) * CNNEW(I) * DELTAT
c
                CONCNODUPD(I) = (MASSURF(I) + MASSPOT(I)) / VOLTOT(I)    
                CONCNODBC(I) = CNNEW(I)
               GO TO 500
C     EVAPORATION / INFILTRATION CASE         
           ELSEIF ((ATMPOT(I).LT.0.0d0).AND.(ATMACT(I).GE.0.0d0)) THEN
               VOLTOT(I) = PONDNOD(I) * ARENOD(I) + ATMPOT(I) * DELTAT
               CONCNODBC(I) =  PONDNOD(I) * ARENOD(I) * CONCNOD(I)
     1                          / VOLTOT(I)
               CONCNODUPD(I) = CONCNODBC(I)
               GO TO 500
           END IF   
C
C
C -------------------------------
C  NO PONDING T
C -------------------------------
        ELSEIF ((IFATM(I).NE.2).AND.(IFATMP(I).NE.2)) THEN
c          IF (IFATMP(I).EQ.0.0d0 THEN
            
C     RAINFALL/INFILTRATION CASE         
           IF ((ATMPOT(I).GT.0.0d0).AND.(ATMACT(I).GE.0.0d0)) THEN
              CONCNODBC(I) = ATMCONC(I)
              CONCNODUPD(I) = 0.0d0
              GO TO 500
C     EVAPORATION/INFILTRATION CASE         
           ELSEIF ((ATMPOT(I).LT.0.0d0).AND.(ATMACT(I).GE.0.0d0)) THEN
CM            WRITE(*,*) 'IMPOSSIBLE CASE AT THE NODE : ',I
              GO TO 500
C     EVAPORATION/EXFILTRATION CASE         
           ELSEIF ((ATMPOT(I).LT.0.0d0).AND.(ATMACT(I).LT.0.0d0)) THEN
CM             write(*,*) 'UNLIKELY CASE - EXTREME EVAPO'              
              GO TO 500
           END IF
C
C------ BUT PONDING at t-1
c     RAINFALL/INFILTRATION CASE 
c            ELSEIF ((IFATM(I).NE.2).AND.(IFATMP(I).NE.0.0d0) THEN
c            VOLTOT(I) = PONDNODP(I) * ARENOD(I) + ATMPOT(I) * DELTAT
c            
c            MASSURF(I) = PONDNODP(I) * ARENOD(I) * CONCNOD(I)
c            MASSPOT(I) = ATMPOT(I) * DELTAT * ATMCONC(I)
cC
c            CONCNODBC(I) = (MASSURF(I) + MASSPOT(I)) / VOLTOT(I)
c            CONCNODUPD(I) = 0.0d0
         END IF
500   CONTINUE  
      END DO      
C     

      IF (ANP_TRA.GT.0) THEN
          DO I=1,ANP_TRA
             CONCNODBC(CONTP_TRA(I))=PRESC_TRA(I)
          END DO
      END IF

C     
      RETURN  
      END
