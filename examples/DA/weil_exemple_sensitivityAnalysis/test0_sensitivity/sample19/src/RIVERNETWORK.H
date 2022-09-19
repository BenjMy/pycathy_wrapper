C
C**************************  COMMON BLOCK RIVERNETWORK ********************
C            
C  integer and real scalars and arrays for surface river network and watershed 
C  characteristics
C            
C***********************************************************************
C
      
      INTEGER   QOI_SN(MAXCEL)
      INTEGER   DTM_P_OUTFLOW_1(COLMAX,ROWMAX)
      INTEGER   DTM_P_OUTFLOW_2(COLMAX,ROWMAX)
      REAL*8    NORTH,SOUTH,EAST,WEST
      REAL*8    DELTA_X,DELTA_Y
      REAL*8    FACTOR
      REAL*8    DTM_W_1(COLMAX,ROWMAX),DTM_W_2(COLMAX,ROWMAX)
      REAL*8    DTM_LOCAL_SLOPE_1(COLMAX,ROWMAX)
      REAL*8    DTM_LOCAL_SLOPE_2(COLMAX,ROWMAX)
      REAL*8    DTM_EPL_1(COLMAX,ROWMAX)
      REAL*8    DTM_EPL_2(COLMAX,ROWMAX)
      REAL*8    DTM_KSS1_SF_1(COLMAX,ROWMAX)
      REAL*8    DTM_KSS1_SF_2(COLMAX,ROWMAX)
      REAL*8    DTM_WS1_SF_1(COLMAX,ROWMAX)
      REAL*8    DTM_WS1_SF_2(COLMAX,ROWMAX)
      REAL*8    DTM_B1_SF(COLMAX,ROWMAX)
      REAL*8    DTM_Y1_SF(COLMAX,ROWMAX)
      REAL*8    DTM_NRC(COLMAX,ROWMAX)

      
      COMMON /RIVERNETWORK/
     1       NORTH,SOUTH,EAST,WEST,FACTOR,DELTA_X,DELTA_Y,
     2       DTM_W_1,DTM_W_2,
     3       DTM_LOCAL_SLOPE_1,DTM_LOCAL_SLOPE_2,
     4       DTM_EPL_1,DTM_EPL_2,
     5       DTM_KSS1_SF_1,DTM_KSS1_SF_2,
     6       DTM_WS1_SF_1,DTM_WS1_SF_2,DTM_B1_SF,DTM_Y1_SF,
     7       DTM_NRC,DTM_P_OUTFLOW_1,DTM_P_OUTFLOW_2,QOI_SN
