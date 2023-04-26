      SUBROUTINE AXBSYM(N,NTERM,TOPOL,JA,COEF1,XVEC,BVEC)
C
C  BVEC:=A.XVEC
C
      IMPLICIT NONE
      INTEGER N,NTERM
      INTEGER I,K,M,MM
cxcx  INTEGER JA(NTERM),TOPOL(N+1)
cxcx  REAL*8  COEF1(NTERM),XVEC(N),BVEC(N)
      INTEGER JA(*),TOPOL(*)
      REAL*8  COEF1(*),XVEC(*),BVEC(*)
C
      DO K=1,N
         BVEC(K)=0.0D0
      END DO
      DO K=1,N
         M=TOPOL(K)
         MM=TOPOL(K+1)-1
         BVEC(K)=BVEC(K)+COEF1(M)*XVEC(JA(M))
         DO I=M+1,MM
            BVEC(K)=BVEC(K)+COEF1(I)*XVEC(JA(I))
            BVEC(JA(I))=BVEC(JA(I))+COEF1(I)*XVEC(K)
         END DO
      END DO
C
      RETURN
      END
