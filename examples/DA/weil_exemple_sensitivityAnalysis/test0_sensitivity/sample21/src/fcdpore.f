C
C*************************** FCDPORE ***********************************
C
C  calculate porosity derivative respect to saturation using Camporese
C  adaptation of Pyatt and John relation for peat volume changes
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FCDPORE(INDEi,SWi,D,PNODIi)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   DEDS
      REAL*8   INDEi,SWi,D,PNODIi
      REAL*8   CBETA,THETA0
      REAL*8   UNO
      PARAMETER (UNO=1.0D0)
      INCLUDE 'SOILCHAR.H'
C
      IF (SWi .LT. UNO) THEN
         CBETA=CBETA0+CANG*D 
         THETA0=PNODIi/(UNO-PNODIi)
         IF (CBETA .GT. UNO) THEN
            CBETA=UNO
         END IF 
         DEDS=(ABS(THETA0+UNO)**(UNO-CBETA)*CBETA*ABS(INDEi*SWi+UNO)**
     1        (CBETA-UNO))/(UNO-ABS(THETA0+UNO)**(UNO-CBETA)*
     2        CBETA*ABS(INDEi*SWi+UNO)**(CBETA-UNO)*SWi)
C        FCDPORE=(1.0D0/(1.0D0+INDEi))**2*DEDS
         FCDPORE=(1.0D0/(1.0D0+INDEi))*DEDS
      ELSE
         FCDPORE=0.0D0
      END IF
C
      RETURN
      END
