C
C************************** ROUTE_TRA **********************************
C  
C  CALCULATION  OF THE INFLOW RATE FOR THE DOWNSTREAM CELL AT T+1 AND T,
C  Q_INFLOW_KKP1_SN AND Q_INFLOW_KK_SN, AND OF THE LATERAL INFLOW OF THE 
C  CELL, Q_OVERLAND.
C
C***********************************************************************
      SUBROUTINE ROUTE_TRA(NROW,NCELL,TIPO_R,RESERVR,
     1     CELLS_R,N_HA,TIME,DELTAT)

      IMPLICIT NONE

      INCLUDE 'CATHY.H'
      INCLUDE 'RIVERNETWORK.H'
      INCLUDE 'SURFWATER.H'
      INCLUDE 'IOUNITS.H'
      INCLUDE 'TRANSPSURF.H'
      
      INTEGER NCELL,NROW
      INTEGER N_O_QUOTA,I,J,I_BASIN,JR
      INTEGER NUM_RR,INDEX_RR
      INTEGER I_RR,J_RR,TIPO_RR,CELLS_RR
      INTEGER III_1,JJJ_1,III_2,JJJ_2
      INTEGER I_CV_1,J_CV_1,I_CV_2,J_CV_2
      INTEGER I_BASIN_CV_1,I_BASIN_CV_2
      INTEGER P_OUTFLOW_1,P_OUTFLOW_2

      INTEGER TIPO_R(MAXCEL)
      INTEGER CELLS_R(MAXRES),N_HA(MAXRES)
      INTEGER RESERVR(MAXCEL)
      
      REAL*8 W_1,W_2
      REAL*8 TIME,DELTAT
      REAL*8 LOCAL_SLOPE_1,LOCAL_SLOPE_2,EPL_1,EPL_2
      REAL*8 WS1_SF_1,WS1_SF_2,KSS1_SF_1,KSS1_SF_2
      REAL*8 B1_SF,Y1_SF,NRC
      REAL*8 QMASS_IN_KK,QMASS_IN_KKP1,QMASS_OUT_KK,QMASS_OUT_KKP1
      REAL*8 SURFACE_CONC,QMASS_OVERLAND_1,QMASS_OVERLAND_2
      REAL*8 H_POOL_KK,H_POOL_KKP1,H_FONDO_RR
      REAL*8 CR,CU,PE,X,C1,AK
C
      CU_MAX=0.0D0
      
C CALCULATION OF THE INFLOW HYDROGRAPHS. THE CELLS ARE PROCESSED
C INTO DESCENDING ELEVATION ORDER

      DO N_O_QUOTA=1,NCELL
         I_BASIN=QOI_SN(N_O_QUOTA)
         JR=MOD(I_BASIN,NROW)
         IF (JR.NE.0) THEN
            J=JR
            I=(I_BASIN-J)/NROW+1
         ELSE
            J=NROW
            I=I_BASIN/NROW
         END IF
           
C FEATURES OF THE CURRENT CELL AS READ IN DATIN.F
     
         W_1=DTM_W_1(I,J)
         W_2=DTM_W_2(I,J)
         P_OUTFLOW_1=DTM_P_OUTFLOW_1(I,J)
         P_OUTFLOW_2=DTM_P_OUTFLOW_2(I,J)
         LOCAL_SLOPE_1=DTM_LOCAL_SLOPE_1(I,J)
         LOCAL_SLOPE_2=DTM_LOCAL_SLOPE_2(I,J)
         EPL_1=DTM_EPL_1(I,J)
         EPL_2=DTM_EPL_2(I,J)

         KSS1_SF_1=DTM_KSS1_SF_1(I,J)
         KSS1_SF_2=DTM_KSS1_SF_2(I,J)
         WS1_SF_1=DTM_WS1_SF_1(I,J)
         WS1_SF_2=DTM_WS1_SF_2(I,J)
         B1_SF=DTM_B1_SF(I,J)        
         Y1_SF=DTM_Y1_SF(I,J)
         NRC=DTM_NRC(I,J)
               
        

C LOCAL CONTRIBUTION TO SURFACE RUNOFF FROM FLOW3D,SURFACE_CONC

         SURFACE_CONC=SURFACE_CONC_SN(I_BASIN)/NRC
      
                 
         IF (W_1.NE.0) THEN
   
C CALCULATION OF LATERAL INFLOW RATE, QMASS_OVERLAND

            CR=1.0D0/EPL_1
            QMASS_OVERLAND_1=SURFACE_CONC*W_1*CR

C CALCULATION OF DOWNSTREAM CELL COORDINATES

            III_1=NINT((P_OUTFLOW_1-5)/3.)
            JJJ_1=P_OUTFLOW_1-5-3*III_1
            I_CV_1=I+III_1      !COORDINATES OF DOWNSTREAM CELL (CARDINAL DIRECTION)
            J_CV_1=J+JJJ_1      !COORDINATES OF DOWNSTREAM CELL (CARDINAL DIRECTION)
            I_BASIN_CV_1 = (I_CV_1 - 1)*NROW + J_CV_1 !COORDINATES OF DOWNSTREAM CELL (CARDINAL DIRECTION) FOR I_BASIN INDEX

            
C     INFLOW AND OUTFLOW HYDROGRAPHS FOR MC.F OR SERBATIO.F      

            QMASS_IN_KK=QMASS_IN_KK_SN(I_BASIN)*W_1/NRC
            QMASS_OUT_KK=QMASS_OUT_KK_SN_1(I_BASIN)/NRC
            QMASS_IN_KKP1=QMASS_IN_KKP1_SN(I_BASIN)*W_1/NRC

            
            IF (TIPO_R(I_BASIN).LT.1000) THEN
               IF (TIPO_R(I_BASIN).NE.0) THEN
                  INDEX_RR=I_BASIN
                  I_RR=I
                  J_RR=J
                  TIPO_RR=TIPO_R(I_BASIN)
                  CELLS_RR=CELLS_R(RESERVR(I_BASIN))
                  H_POOL_KK=H_POOL_KK_VEC(RESERVR(I_BASIN))
                  NUM_RR = RESERVR(I_BASIN)
                  H_FONDO_RR=H_FONDO(RESERVR(I_BASIN))
 
                  CALL SERBATOIO(INDEX_RR,I_RR,J_RR,TIPO_RR,QMASS_IN_KK,
     1                 QMASS_IN_KKP1,QMASS_OUT_KKP1,H_POOL_KKP1,
     2                 H_POOL_KK,DELTAT,SURFACE_CONC,
     3                 NUM_RR,N_HA,HRES,ARES,H_SOGLIA,
     4                 H_FONDO_RR) 
         
                  IF(QMASS_OUT_KKP1.LT.0.0D0) THEN
!     WRITE(IOUT41,*)'ATTENZIONE! PORTATE USCENTI NEGATIVE'
                     QMASS_OUT_KKP1=0.0D0
                  END IF

                  H_POOL_KKP1_VEC(RESERVR(I_BASIN))=H_POOL_KKP1
                  QMASS_OUT_KKP1_SN_1(I_BASIN)=QMASS_OUT_KKP1
                  
               ELSE             ! IF (TIPO_R(I_BASIN).EQ.0) THEN
 
C     MUSKINGUM-CUNGE
                  CALL MC(LOCAL_SLOPE_1,EPL_1,KSS1_SF_1,WS1_SF_1,B1_SF,
     1             Y1_SF,DELTAT,CU,PE,X,C1,AK,QMASS_IN_KK,QMASS_IN_KKP1,
     2                 QMASS_OUT_KK,QMASS_OVERLAND_1,QMASS_OUT_KKP1)
                  
                  IF(QMASS_OUT_KKP1.LT.0.0D0) THEN
!     CXCX          WRITE(IOUT41,*)'ATTENZIONE! PORTATE USCENTI NEGATIVE'
                     QMASS_OUT_KKP1=0.0D0
                  END IF

                  QMASS_OUT_KKP1_SN_1(I_BASIN)=QMASS_OUT_KKP1*NRC

                  IF (CU .GE. CU_MAX) THEN
                     CU_MAX=CU
                     AK_MAX=AK
                  END IF
                  
c                  WRITE(IOUT44,1500)TIME,I,J,CU,PE,X,C1

               END IF
           ELSE                  ! IF (TIPO_R(INDEX_C).GE.1000) THEN
              WRITE(6,*) 'INCOMPRENSIBILE CASE!'
              STOP
C             QMASS_OUT_KKP1=QMASS_IN_KKP1 + SURFACE_CONC             
C             QMASS_OUT_KKP1_SN_1(I_BASIN)=QMASS_OUT_KKP1
            END IF
            
C     SUPERIMPOSITION OF OUTFLOW HYDROGRAPHS FOR THE DOWNSTREAM CELLS
            IF (N_O_QUOTA .NE. NCELL) THEN
               QMASS_IN_KKP1_SN(I_BASIN_CV_1) = 
     1              QMASS_IN_KKP1_SN(I_BASIN_CV_1) + 
     1              QMASS_OUT_KKP1_SN_1(I_BASIN)
            END IF                 
         END IF                 !  IF (W_1.NE.0)

         IF (W_2.NE.0) THEN    
C CALCULATION OF LATERAL INFLOW RATE, QMASS_OVERLAND
            CR=1.0D0/EPL_2
            QMASS_OVERLAND_2=SURFACE_CONC*W_2*CR
C CALCULATION OF DOWNSTREAM CELL COORDINATES            
            III_2=NINT((P_OUTFLOW_2-5)/3.)
            JJJ_2=P_OUTFLOW_2-5-3*III_2
            I_CV_2=I+III_2      !COORDINATES OF DOWNSTREAM CELLS (DIAGONAL DIRECTION)
            J_CV_2=J+JJJ_2      !COORDINATES OF DOWNSTREAM CELLS (DIAGONAL DIRECTION)
            I_BASIN_CV_2 = (I_CV_2 - 1)*NROW + J_CV_2 !COORDINATES OF DOWNSTREAM CELL (DIAGONAL DIRECTION) FOR I_BASIN INDEX

            
C     INFLOW AND OUTFLOW HYDROGRAPHS FOR MC.F OR SERBATIO.F      
            QMASS_IN_KK=QMASS_IN_KK_SN(I_BASIN)*W_2/NRC
            QMASS_OUT_KK=QMASS_OUT_KK_SN_2(I_BASIN)/NRC
            QMASS_IN_KKP1=QMASS_IN_KKP1_SN(I_BASIN)*W_2/NRC     

            IF (TIPO_R(I_BASIN).LT.1000) THEN
              IF (TIPO_R(I_BASIN).NE.0) THEN
                 INDEX_RR=I_BASIN
                 I_RR=I
                 J_RR=J
                 TIPO_RR=TIPO_R(I_BASIN)
                 CELLS_RR=CELLS_R(RESERVR(I_BASIN))
                 H_POOL_KK=H_POOL_KK_VEC(RESERVR(I_BASIN))
                 NUM_RR = RESERVR(I_BASIN)
                 H_FONDO_RR=H_FONDO(RESERVR(I_BASIN))

                 CALL SERBATOIO(INDEX_RR,I_RR,J_RR,TIPO_RR,QMASS_IN_KK,
     1                 QMASS_IN_KKP1,QMASS_OUT_KKP1,H_POOL_KKP1,
     2                 H_POOL_KK,DELTAT,SURFACE_CONC,
     3                 NUM_RR,N_HA,HRES,ARES,H_SOGLIA,
     4                 H_FONDO_RR) 
        
                 IF(QMASS_OUT_KKP1.LT.0.0D0) THEN
c      WRITE(IOUT41,*)'ATTENZIONE! PORTATE USCENTI NEGATIVE'
                    QMASS_OUT_KKP1=0.0D0
                 END IF

                 H_POOL_KKP1_VEC(RESERVR(I_BASIN))=H_POOL_KKP1
                 QMASS_OUT_KKP1_SN_2(I_BASIN)=QMASS_OUT_KKP1

              ELSE             ! IF (TIPO_R(I_BASIN).EQ.0) THEN

C        MUSKINGUM-CUNGE
                 
                 CALL MC(LOCAL_SLOPE_2,EPL_2,KSS1_SF_2,WS1_SF_2,B1_SF,
     1             Y1_SF,DELTAT,CU,PE,X,C1,AK,QMASS_IN_KK,QMASS_IN_KKP1,
     2                QMASS_OUT_KK,QMASS_OVERLAND_2,QMASS_OUT_KKP1)

                  IF(QMASS_OUT_KKP1.LT.0.0D0) THEN
!     CXCX          WRITE(IOUT41,*)'ATTENZIONE! PORTATE USCENTI NEGATIVE'
                     QMASS_OUT_KKP1=0.0D0
                  END IF

              
                  QMASS_OUT_KKP1_SN_2(I_BASIN)=QMASS_OUT_KKP1*NRC
                  
                  IF (CU .GE. CU_MAX) THEN
                     CU_MAX=CU
                     AK_MAX=AK
                  END IF

c                  WRITE(IOUT44,1500)TIME,I,J,CU,PE,X,C1
                  
               END IF
            ELSE                ! IF (TIPO_R(I_BASIN).GE.1000) THEN
 
               WRITE(6,*) 'INCOMPRENSIBILE CASE!'
               STOP              
C               QMASS_OUT_KKP1=QMASS_IN_KKP1 + SURFACE_CONC             
C               QMASS_OUT_KKP1_SN_2(I_BASIN)=QMASS_OUT_KKP1
            END IF
            

C     SUPERIMPOSITION OF OUTFLOW HYDROGRAPHS FOR THE DOWNSTREAM CELLS
            QMASS_IN_KKP1_SN(I_BASIN_CV_2)= 
     1           QMASS_IN_KKP1_SN(I_BASIN_CV_2) + 
     1           QMASS_OUT_KKP1_SN_2(I_BASIN)
            
         END IF                 !  IF (W_2.NE.0)      

      
      END DO                    !  DO N_O_QUOTA=1,N_CELLE
    
     
      RETURN 
 1500 FORMAT(1PE16.8,2(2X,I4),4(2X,E15.7))
      END
