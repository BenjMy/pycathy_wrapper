C
C**************************  STRPIC ************************************
C
C  routine per l'analisi topologica del reticolo nel caso simmetrico
C
C***********************************************************************
C
      SUBROUTINE STRPIC(N,NTERM,TETRA,JA,TOPOL,NT,N1,IMAX)
C
      IMPLICIT    NONE
      INTEGER     I,J,K,L,M,III,KK,KKK,MM,MCONTR
      INTEGER     N,N1,NT,NTERM,IMAX
      INTEGER     I1(4),I2(4)
      INTEGER     TOPOL(*),JA(*),TETRA(5,*)
      INCLUDE    'IOUNITS.H'
C
C  fissa gli elementi diagonali a distanza costante N1
C        
      NTERM=N1*N
      CALL INIT0I(NTERM,JA)
      JA(1)=1
      DO K=1,NTERM-1
         IF (K/N1*N1 .EQ. K) JA(K+1)=K/N1+1
      END DO
C
C  analizza tutti i tetraedri ( NT )
C
      DO 400 K=1,NT
         DO I=1,4
            I2(I)=TETRA(I,K)
         END DO
C
C  ordina i nodi dell'elemento in senso crescente 
C
         DO J=1,4
            KKK=1
            KK=I2(1)
            DO I=2,4
               IF(I2(I).LT.KK) THEN
                  KK=I2(I)
                  KKK=I
               END IF
            END DO
            I1(J)=KK
            I2(KKK)=IMAX
         END DO
C  
C  genera il vettore JA
C
         DO 6 I=1,3
            J=I+1
            DO 7 L=J,4
               M=N1*(I1(I)-1)+L-J+2
               MCONTR=N1*I1(I)
    9          IF(I1(L)-JA(M))8,7,11
   11          IF(JA(M).EQ.0)GO TO 10
               M=M+1
               IF(M-MCONTR.GE.0)GO TO 99
               GO TO 9
   10          JA(M)=I1(L)
               GO TO 7
    8          MM=M
   18          MM=MM+1
               IF(MM-MCONTR.GE.0) GO TO 99
               IF(JA(MM))18,13,18
   13          JA(MM)=JA(MM-1)
               MM=MM-1
               IF(MM.GT.M)GO TO 13
               JA(M)=I1(L)
    7       CONTINUE
    6    CONTINUE
  400 CONTINUE
C
C  costruisce il vettore TOPOL
C
      TOPOL(1)=1
      M=1
      J=1
      DO 20 K=1,NTERM,N1
         DO 21 I=1,N1
            IF(JA(K+I-1).EQ.0)GO TO 21
            M=M+1
   21    CONTINUE
         J=J+1
         TOPOL(J)=M
   20 CONTINUE
C
C  compatta il vettore JA eliminando gli zeri
C
      M=0
      DO 14 K=1,NTERM
         IF(JA(K).EQ.0) GO TO 14
         M=M+1
         JA(M)=JA(K)
   14 CONTINUE
      NTERM=M
C
      RETURN
   99 III=MCONTR/N1
      WRITE(IOUT2,101)III,K
  101 FORMAT(1X,'RIGA = ',I6,' TETRAEDRO = ',I6,'  AUMENTARE N1 ')
      CALL CLOSIO
      STOP
      END
