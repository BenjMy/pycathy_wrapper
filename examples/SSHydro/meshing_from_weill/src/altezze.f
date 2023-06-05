C
C************************** ALTEZZE *************************************
C This subroutines calculates the water heights by performing the water
C balance at the surface cells.
C Warning: in the calculation of delta_v it has been deleted the multiplicative
C factor delta_x*delta_y since the surface_water term has already been 
C calculated in m3/s by the flow3d module.
C***********************************************************************
C

      SUBROUTINE ALTEZZE(NCELL,TIPO_R,DELTA_X,DELTA_Y,
     1                  Q_IN_KK_SN,Q_IN_KKP1_SN,
     2                  Q_OUT_KK_SN_1,Q_OUT_KK_SN_2,
     3                  Q_OUT_KKP1_SN_1,Q_OUT_KKP1_SN_2,
     4                  VOLUME_KK_SN,VOLUME_KKP1_SN,SURFACE_WATER_SN,
     5                  H_WATER_KKP1_SN,DELTA_T,
     6                  H_POOL_KKP1_VEC,NUM_R,ELEVATION) 

      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INCLUDE 'IOUNITS.H'

      INTEGER INDEX_C,NCELL,AA
      INTEGER TIPO_R(*),NUM_R(*)
      REAL*8  H_WATER_KKP1
      REAL*8  Q_IN_KK,Q_IN_KKP1,Q_OUT_KK_1,Q_OUT_KK_2
      REAL*8  Q_OUT_KKP1_1,Q_OUT_KKP1_2
      REAL*8  DELTA_X,DELTA_Y,DELTA_V,DELTA_T,VOLUME_KK,VOLUME_KKP1
      REAL*8  SURFACE_WATER
      REAL*8  Q_IN_KK_SN(*),Q_IN_KKP1_SN(*)
      REAL*8  Q_OUT_KK_SN_1(*), Q_OUT_KK_SN_2(*)
      REAL*8  Q_OUT_KKP1_SN_1(*),Q_OUT_KKP1_SN_2(*)
      REAL*8  H_WATER_KKP1_SN(*)
      REAL*8  SURFACE_WATER_SN(*),ELEVATION(*)
      REAL*8  H_POOL_KKP1_VEC(*)
      REAL*8  VOLUME_KKP1_SN(*),VOLUME_KK_SN(*)
C
      DO INDEX_C=1,MAXCEL
 
       IF(TIPO_R(INDEX_C).EQ.0) THEN

          Q_IN_KK=Q_IN_KK_SN(INDEX_C)
          Q_IN_KKP1=Q_IN_KKP1_SN(INDEX_C)
          Q_OUT_KK_1=Q_OUT_KK_SN_1(INDEX_C)
          Q_OUT_KK_2=Q_OUT_KK_SN_2(INDEX_C)
          Q_OUT_KKP1_1=Q_OUT_KKP1_SN_1(INDEX_C)
          Q_OUT_KKP1_2=Q_OUT_KKP1_SN_2(INDEX_C)
          SURFACE_WATER=SURFACE_WATER_SN(INDEX_C)

C   INITIAL VOLUME

          VOLUME_KK=VOLUME_KK_SN(INDEX_C)

C   DELTA_V CALCULATION

          DELTA_V= (Q_IN_KK+Q_IN_KKP1)/2*DELTA_T +
     &              SURFACE_WATER*DELTA_T -
     &             (Q_OUT_KK_1+Q_OUT_KK_2)/2*DELTA_T -
     &             (Q_OUT_KKP1_1+Q_OUT_KKP1_2)/2*DELTA_T

C   FINAL VOLUME 
 
          VOLUME_KKP1=VOLUME_KK + DELTA_V
          VOLUME_KKP1_SN(INDEX_C)=VOLUME_KKP1

C   CALCULATION OF WATER HEIGHTS

          IF(VOLUME_KKP1 .GE. 0.D0) THEN
             H_WATER_KKP1=VOLUME_KKP1/(DELTA_X*DELTA_Y)
          ELSE
             VOLUME_KKP1_SN(INDEX_C)=0.0D0
             H_WATER_KKP1=0.0D0
          END IF

          H_WATER_KKP1_SN(INDEX_C)=H_WATER_KKP1
          
       
       ELSE IF(TIPO_R(INDEX_C).LT.1000) THEN
          
          H_WATER_KKP1_SN(INDEX_C)=H_POOL_KKP1_VEC(NUM_R(INDEX_C))-
     &                            ELEVATION(INDEX_C) 

          IF(H_WATER_KKP1_SN(INDEX_C).LT.0) THEN
              H_WATER_KKP1_SN(INDEX_C)=0.0D0
          END IF

       ELSE  
 
          AA=NUM_R(TIPO_R(INDEX_C)-1000)
          
          H_WATER_KKP1=H_POOL_KKP1_VEC(AA)-
     &                 ELEVATION(INDEX_C)

          IF(H_WATER_KKP1.GT.0.D0) THEN
             
             H_WATER_KKP1_SN(INDEX_C)=H_WATER_KKP1
             
          ELSE 

            H_WATER_KKP1_SN(INDEX_C)=0.0D0

          END IF

        END IF

      END DO

      RETURN
      END
