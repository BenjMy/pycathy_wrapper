C
C**************************  COMMON BLOCK NORMVL ***********************
C            
C  integer and real arrays for calculation of residual and
C  difference norms
C            
C***********************************************************************
C
      INTEGER ITUMAX(MAXIT)
      REAL*8  PIKMXV(MAXIT),PCURRV(MAXIT),PPREVV(MAXIT)   
      REAL*8  PL2V(MAXIT),FINFV(MAXIT),FL2V(MAXIT) 
      COMMON  /NORMVL/PIKMXV,PCURRV,PPREVV,PL2V,FINFV,FL2V,ITUMAX
