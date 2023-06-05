C**********************************************************************
C                            SUBROUTINE NEUMANN
C compute vertically duplicated BCs and free-drainage fluxes
C**********************************************************************
      SUBROUTINE NEUMANN(TIME,NNEU,NNEUC,NNOD,NSTR,CKRW,KZNOD,
     1                   ARENOD,ACONTQ,QINP,Q)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER I,K,NN
      INTEGER NNOD,NSTR,STR
      INTEGER NNEU(3),NNEUC(3),ACONTQ(*)
      REAL*8  CKRW(*),ARENOD(*),Q(*),TIME,SUM1,SUM2
      REAL*8  KZNOD(*),QINP(3,*)
C
      K = 0
      SUM1=0.0d0
      SUM2=0.0d0
      IF (NNEU(2) .LT. 0) THEN
          DO I=1,NNOD
             K = K + 1
             NN = ACONTQ(K)
             Q(K)=-1.0d0*ARENOD(I)*CKRW(NN)*KZNOD(NN)
             SUM1=SUM1+Q(K)
          END DO
          DO I = 1,NNEUC(2)
             K = K + 1
             NN = ACONTQ(K)
             STR = MIN(INT(REAL(NN-1)/NNOD)+1,NSTR)
             Q(K) = QINP(2,K)
             SUM2=SUM2+Q(K)
          END DO
cm        write(999,*)TIME,SUM1,SUM2
      END IF
C
      RETURN
      END
