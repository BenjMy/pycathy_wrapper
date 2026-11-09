C     blas_lite.f : the only symbols the CATHY sources need from a system BLAS/LAPACK
C     (checked by linking CATHY_src_WithIrrigation code/ without any library: undefined
C     references were exactly DCOPY, DGETRF and DGETRS; DAXPY, DDOT and DNRM2 are defined
C     in solscal-extended.f). DGETRF/DGETRS are only used for a 4 x 4 solve (iperplane.f).
C     Plain reference algorithms with the standard interfaces, so that no system BLAS/LAPACK
C     (-lblas -llapack) is needed to link the solver.
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
      INTEGER N,INCX,INCY,I,IX,IY
      DOUBLE PRECISION DX(*),DY(*)
      IF (N.LE.0) RETURN
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DY(IY) = DX(IX)
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
      END

C     LU factorisation with partial pivoting (unblocked, same result layout as LAPACK DGETRF)
      SUBROUTINE DGETRF(M,N,A,LDA,IPIV,INFO)
      INTEGER M,N,LDA,INFO,IPIV(*),I,J,K,JP
      DOUBLE PRECISION A(LDA,*),AMAX,T,PIV
      INFO = 0
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
      DO 40 J = 1,MIN(M,N)
         JP = J
         AMAX = ABS(A(J,J))
         DO 10 I = J+1,M
            IF (ABS(A(I,J)).GT.AMAX) THEN
               AMAX = ABS(A(I,J))
               JP = I
            END IF
   10    CONTINUE
         IPIV(J) = JP
         IF (A(JP,J).NE.0.0D0) THEN
            IF (JP.NE.J) THEN
               DO 20 K = 1,N
                  T = A(J,K)
                  A(J,K) = A(JP,K)
                  A(JP,K) = T
   20          CONTINUE
            END IF
            PIV = 1.0D0/A(J,J)
            DO 30 I = J+1,M
               A(I,J) = A(I,J)*PIV
   30       CONTINUE
         ELSE IF (INFO.EQ.0) THEN
            INFO = J
         END IF
         IF (J.LT.MIN(M,N)) THEN
            DO 38 K = J+1,N
               DO 36 I = J+1,M
                  A(I,K) = A(I,K) - A(I,J)*A(J,K)
   36          CONTINUE
   38       CONTINUE
         END IF
   40 CONTINUE
      RETURN
      END

C     Solve A*X = B ('N') or A**T*X = B ('T','C') with the factors of DGETRF
      SUBROUTINE DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
      CHARACTER TRANS
      INTEGER N,NRHS,LDA,LDB,INFO,IPIV(*),I,J,K,IP
      DOUBLE PRECISION A(LDA,*),B(LDB,*),T
      LOGICAL NOTRAN
      INFO = 0
      NOTRAN = (TRANS.EQ.'N' .OR. TRANS.EQ.'n')
      IF (.NOT.NOTRAN .AND. TRANS.NE.'T' .AND. TRANS.NE.'t'
     &    .AND. TRANS.NE.'C' .AND. TRANS.NE.'c') THEN
         INFO = -1
         RETURN
      END IF
      IF (N.EQ.0 .OR. NRHS.EQ.0) RETURN
      IF (NOTRAN) THEN
         DO 10 J = 1,N
            IP = IPIV(J)
            IF (IP.NE.J) THEN
               DO 5 K = 1,NRHS
                  T = B(J,K)
                  B(J,K) = B(IP,K)
                  B(IP,K) = T
    5          CONTINUE
            END IF
   10    CONTINUE
         DO 30 K = 1,NRHS
            DO 20 J = 1,N
               IF (B(J,K).NE.0.0D0) THEN
                  DO 15 I = J+1,N
                     B(I,K) = B(I,K) - B(J,K)*A(I,J)
   15             CONTINUE
               END IF
   20       CONTINUE
            DO 28 J = N,1,-1
               IF (B(J,K).NE.0.0D0) THEN
                  B(J,K) = B(J,K)/A(J,J)
                  DO 25 I = 1,J-1
                     B(I,K) = B(I,K) - B(J,K)*A(I,J)
   25             CONTINUE
               END IF
   28       CONTINUE
   30    CONTINUE
      ELSE
         DO 60 K = 1,NRHS
            DO 40 J = 1,N
               T = B(J,K)
               DO 35 I = 1,J-1
                  T = T - A(I,J)*B(I,K)
   35          CONTINUE
               B(J,K) = T/A(J,J)
   40       CONTINUE
            DO 50 J = N,1,-1
               T = B(J,K)
               DO 45 I = J+1,N
                  T = T - A(I,J)*B(I,K)
   45          CONTINUE
               B(J,K) = T
   50       CONTINUE
   60    CONTINUE
         DO 80 J = N,1,-1
            IP = IPIV(J)
            IF (IP.NE.J) THEN
               DO 70 K = 1,NRHS
                  T = B(J,K)
                  B(J,K) = B(IP,K)
                  B(IP,K) = T
   70          CONTINUE
            END IF
   80    CONTINUE
      END IF
      RETURN
      END
