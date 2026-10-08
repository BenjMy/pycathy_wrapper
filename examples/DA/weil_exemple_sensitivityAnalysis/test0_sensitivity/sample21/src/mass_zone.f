C
C**************************  MASS_ZONE *********************************
C
C  Compute the solute total mass under each zones at each time step
C
C***********************************************************************
C
      SUBROUTINE MASS_ZONE(NT,TETRA,NZONE,MASSZONE,CNEW,VOLU,PEL,SWNEW,
     1              VOLZONE,MASSSOLZONE,KEL)

C
      IMPLICIT NONE
c     INCLUDE 'CATHY.H'
      INTEGER  K,ZO,MTYPE
      INTEGER  NT,NZONE
      INTEGER  TETRA(5,*)
      REAL*8   MASSZONE(*),VOLU(*),CNEW(*),PEL(*),SWNEW(*)
      REAL*8   VOLZONE(*),MASSSOLZONE(*),KEL(*)
c      REAL*8   CNEW_justformass(NT),CNNEW(*)
c      INCLUDE 'IOUNITS.H'
c      INCLUDE 'SOILCHAR.H'
C
        CALL INIT0R(NZONE,MASSZONE)
        CALL INIT0R(NZONE,MASSSOLZONE)
        CALL INIT0R(NZONE,VOLZONE)
        DO ZO=1,NZONE
            DO K=1,NT
            MTYPE=TETRA(5,K)
                IF (MTYPE.EQ.ZO) THEN
                MASSZONE(ZO)=MASSZONE(ZO)+
     1                   CNEW(K)*VOLU(K)*PEL(K)*SWNEW(K)
                VOLZONE(ZO)=VOLZONE(ZO)+
     1                   VOLU(K)*PEL(K)*SWNEW(K)
                MASSSOLZONE(ZO)= MASSSOLZONE(ZO) +
     1                   KEL(K)*CNEW(K)*VOLU(K)*2*(1-PEL(K))
                END IF
            END DO
        END DO
            
          DO ZO=1,NZONE  
          IF ((MASSZONE(ZO).GT.0.0d0).AND.(MASSZONE(ZO).LE.1e-30)) THEN
            MASSZONE(ZO)=0.0d0
          END IF
          IF ((VOLZONE(ZO).GT.0.0d0).AND.
     1                      (VOLZONE(ZO).LE.1e-30)) THEN
            VOLZONE(ZO)=0.0d0
          END IF
          IF ((MASSSOLZONE(ZO).GT.0.0d0).AND.
     1                      (MASSSOLZONE(ZO).LE.1e-30)) THEN
            MASSSOLZONE(ZO)=0.0d0
          END IF
          END DO
      RETURN
      
      END
