C
C************************** MASBIL *************************************
c THIS ROUTINE PERFORMS MASSE BALANCE AT THE SURFACE AND GIVES
C THE CONCENTRATIONS
C***********************************************************************
C

      SUBROUTINE MASBIL(NCELL,TIPO_R,DELTA_X,DELTA_Y,
     1                  QMASS_IN_KK_SN,QMASS_IN_KKP1_SN,
     2                  QMASS_OUT_KK_SN_1,QMASS_OUT_KK_SN_2,
     3                  QMASS_OUT_KKP1_SN_1,QMASS_OUT_KKP1_SN_2,
     4                  MASS_KK_SN,MASS_KKP1_SN,SURFACE_CONC_SN,
     5                  CONC_KKP1_SN,VOLUME_KKP1_SN,DELTA_T,MIX_CORRECT,
     6                  surface_mix_sn) 

      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INCLUDE 'IOUNITS.H'

      INTEGER INDEX_C,NCELL,k
      INTEGER TIPO_R(*)
      REAL*8  CONC_KKP1
      REAL*8  QMASS_IN_KK,QMASS_IN_KKP1,QMASS_OUT_KK_1,QMASS_OUT_KK_2
      REAL*8  QMASS_OUT_KKP1_1,QMASS_OUT_KKP1_2
      REAL*8  DELTA_X,DELTA_Y,DELTA_M,DELTA_T,MASS_KK,MASS_KKP1
      REAL*8  SURFACE_CONC,VOLUME_KKP1
      REAL*8  QMASS_IN_KK_SN(*),QMASS_IN_KKP1_SN(*)
      REAL*8  QMASS_OUT_KK_SN_1(*), QMASS_OUT_KK_SN_2(*)
      REAL*8  QMASS_OUT_KKP1_SN_1(*),QMASS_OUT_KKP1_SN_2(*)
      REAL*8  CONC_KKP1_SN(*),VOLUME_KKP1_SN(*)
      REAL*8  SURFACE_CONC_SN(*)
      REAL*8  MASS_KKP1_SN(*),MASS_KK_SN(*),MIX_CORRECT(*)
      REAL*8  SURFACE_MIX_SN(*),SURFACE_MIX
C
        CALL INIT0R(ncell,MIX_CORRECT)
C
      DO INDEX_C=1,NCELL
 
         IF(TIPO_R(INDEX_C).EQ.0) THEN
            
            QMASS_IN_KK=QMASS_IN_KK_SN(INDEX_C)
            QMASS_IN_KKP1=QMASS_IN_KKP1_SN(INDEX_C)
            QMASS_OUT_KK_1=QMASS_OUT_KK_SN_1(INDEX_C)
            QMASS_OUT_KK_2=QMASS_OUT_KK_SN_2(INDEX_C)
            QMASS_OUT_KKP1_1=QMASS_OUT_KKP1_SN_1(INDEX_C)
            QMASS_OUT_KKP1_2=QMASS_OUT_KKP1_SN_2(INDEX_C)
            SURFACE_CONC=SURFACE_CONC_SN(INDEX_C)
            SURFACE_MIX=SURFACE_MIX_SN(INDEX_C)
            VOLUME_KKP1 = VOLUME_KKP1_SN(INDEX_C)
            
C   INITIAL VOLUME
          
            MASS_KK=MASS_KK_SN(INDEX_C)
          
C   DELTA_V CALCULATION
            
            DELTA_M= (QMASS_IN_KK+QMASS_IN_KKP1)/2*DELTA_T +
     &           SURFACE_CONC*DELTA_T -
     &           (QMASS_OUT_KK_1+QMASS_OUT_KK_2)/2*DELTA_T -
     &           (QMASS_OUT_KKP1_1+QMASS_OUT_KKP1_2)/2*DELTA_T
C
C     FINAL VOLUME 
          
          MASS_KKP1= MASS_KK + DELTA_M
          
C   CALCULATION OF WATER HEIGHTS
          IF(MASS_KKP1 .GE. 0.D0 .AND. VOLUME_KKP1 .GT. 1e-8) THEN
             CONC_KKP1=MASS_KKP1/VOLUME_KKP1
             IF (CONC_KKP1.GE.1.0e30) THEN
             write(*,*) INDEX_C,CONC_KKP1,MASS_KKP1,VOLUME_KKP1
             stop
             END IF
          ELSE
             MIX_CORRECT(INDEX_C)=MASS_KKP1
             MASS_KKP1=0.0d0
             CONC_KKP1=0.0D0
          END IF
          
          MASS_KKP1_SN(INDEX_C)=MASS_KKP1 
          CONC_KKP1_SN(INDEX_C)=CONC_KKP1
          
       END IF
      END DO
C
      RETURN
      END
