C
C**************************  NUDONE ************************************
C
C  read and initialize nudging parameters and arrays.
C  Note that no nudging is done during the first time step, so the
C  observation values NUDVAL are only read in subroutine NUDNXT (the
C  implied assumption is that the first nudging observation time is at
C  least NUDTAU time units away from times 0.0 and DELTAT; in other
C  words "re-initialization" is presumably not needed so close to the
C  initial conditions).
C
C***********************************************************************
C
      SUBROUTINE NUDONE(N,NT,TETRA,X,Y,Z,XC,YC,ZC,VOLU,
     1                  NUDN,NUDCTR,NUDT,NUDC,NUDG,NUDFLAG,WFLAG,
     2                  NUDTIM,NUDRXY,NUDRZ,NUDTAU,
     3                  NUDX,NUDY,NUDZ,NUDEPS,NUDTET,NUDCUM)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  I,J
      INTEGER  N,NT
      INTEGER  NUDN,NUDCTR,NUDT,NUDC,NUDFLAG,WFLAG
      INTEGER  TETRA(5,*)
      INTEGER  NUDTET(*)
      REAL*8   NUDG
      REAL*8   X(*),Y(*),Z(*),XC(*),YC(*),ZC(*),VOLU(*)
      REAL*8   NUDTIM(*),NUDRXY(*),NUDRZ(*),NUDTAU(MAXNUDT,*)
      REAL*8   NUDX(*),NUDY(*),NUDZ(*),NUDEPS(*),NUDCUM(*)
      INCLUDE 'IOUNITS.H'
C
C  input the nudging parameters. Each observation point (and time) can
C  be assigned a different value of NUDRXY, NUDRZ, NUDEPS (and NUDTAU),
C  but for now only one value is read in for each of these parameters
C  and replicated.
C
      READ(IIN50,*) NUDN,NUDFLAG,WFLAG
      WRITE(IOUT2,1000) NUDN
      IF (NUDN .EQ. 0) THEN
         NUDT=0
         NUDG=0.0
         GO TO 900
      END IF
      READ(IIN50,*) NUDT,NUDG
      WRITE(IOUT2,1010) NUDT,NUDG
      READ(IIN50,*) (NUDTIM(I),I=1,NUDT)
      WRITE(IOUT2,1020) (NUDTIM(I),I=1,NUDT)
      DO I=1,NUDN
         READ(IIN50,*) NUDX(I),NUDY(I),NUDZ(I),NUDEPS(I)
      END DO
CM    READ(IIN50,*) NUDRXY(1),NUDRZ(1),NUDTAU(1,1),NUDEPS(1)
      READ(IIN50,*) NUDRXY(1),NUDRZ(1),NUDTAU(1,1)
      WRITE(IOUT2,1030) NUDRXY(1),NUDRZ(1),NUDTAU(1,1),NUDEPS(1)
      DO I=2,NUDN
         NUDRXY(I)=NUDRXY(1)
         NUDRZ(I)=NUDRZ(1)
CM       NUDEPS(I)=NUDEPS(1)
         NUDTAU(1,I)=NUDTAU(1,1)
         DO J=2,NUDT
            NUDTAU(J,I)=NUDTAU(1,1)
         END DO
      END DO
      DO J=2,NUDT
         NUDTAU(J,1)=NUDTAU(1,1)
      END DO
C
C  locate the nudging observation points within our grid  
C
      CALL NUDLOCATE(NT,NUDN,TETRA,NUDTET,
     1               NUDX,NUDY,NUDZ,X,Y,Z,VOLU)
      WRITE(IOUT2,1040)
      DO I=1,NUDN
         J=NUDTET(I)
         WRITE(IOUT2,1050) I,NUDX(I),NUDY(I),NUDZ(I),J,XC(J),YC(J),ZC(J)
      END DO
C
C  initialize NUDCTR, NUDC, and NUDCUM
C
      NUDCTR=0
C
  900 CONTINUE
      NUDC=0
      CALL INIT0R(N,NUDCUM)
C
      RETURN
 1000 FORMAT(/,5X,'NUDN   (# OF NUDGING OBSERVATION POINTS)= ',I6)
 1010 FORMAT(  5X,'NUDT   (# OF NUDGING OBSERVATION TIMES) = ',I6,
     1       /,5X,'NUDG   (NUDGING FACTOR "G")             = ',1PE15.5)
 1020 FORMAT(  5X,'OBSERVATION TIMES',/,(5(2X,1PE12.4)))
 1030 FORMAT(  5X,'NUDRXY (HORIZONTAL RADIUS OF INFLUENCE) = ',1PE15.5,
     1       /,5X,'NUDRZ  (VERTICAL RADIUS OF INFLUENCE)   = ',1PE15.5,
     2       /,5X,'NUDTAU (HALF PERIOD OF TIME INFL WINDOW)= ',1PE15.5,
     3       /,5X,'NUDEPS (OBS.DATA QUALITY FACTOR "EPS")  = ',1PE15.5)
 1040 FORMAT(' OBS.POINT   X-COORD   Y-COORD   Z-COORD  NUDTET ',
     1       ' XCnudtet  YCnudtet  ZCnudtet')
 1050 FORMAT(5X,I5,3(1PE10.2),2X,I6,3(1PE10.2))
      END
