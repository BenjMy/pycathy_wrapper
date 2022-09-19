C
C*************************** FCINDE ************************************
C
C  calculate void ratio index using Camporese
C  adaptation of Pyatt and John relation for peat volume changes
C
C***********************************************************************
C
      DOUBLE PRECISION FUNCTION FCINDE(INDE,SW,D,PSI,PNODI,SNODI,PORE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      REAL*8   INDE,SW,D,PSI,PNODI,SNODI,PORE
      REAL*8   CBETA,THETA0
      REAL*8   UNO
      PARAMETER (UNO=1.0D0)
      INCLUDE 'SOILCHAR.H'
C
      IF (SW .LT. UNO) THEN
         CBETA=CBETA0+CANG*D 
         THETA0=PNODI/(UNO-PNODI)
         IF (CBETA .GT. UNO) THEN
            CBETA=UNO
         END IF
         FCINDE=ABS(THETA0+UNO)**(UNO-CBETA)*ABS(INDE*SW+UNO)**CBETA
     1          -UNO 
      ELSE
         PORE=UNO-(UNO-PNODI)*EXP(-SNODI*PSI)
         FCINDE=PORE/(UNO-PORE)
      END IF
C
      RETURN
      END
