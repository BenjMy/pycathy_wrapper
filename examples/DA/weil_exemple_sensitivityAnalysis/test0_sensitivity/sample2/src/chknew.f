C
C**************************  CHKNEW ************************************
C
C  routine per controllare la correttezza della rappresentazione
C  della matrice in forma compatta nel caso non simmetrico
C  calcola il numero di elementi diagonali nulli della matrice
C  ritornandolo in NDZ
C
C***********************************************************************
C
      SUBROUTINE CHKNEW(N,NTERM,TOPOL,JA,NDZ,IER)
C
      IMPLICIT  NONE
      INTEGER   I,K,IS,IE,JS,JE,IEM1
      INTEGER   N,NTERM,NDZ
      INTEGER   TOPOL(*),IER(*),JA(*)
      INCLUDE  'IOUNITS.H'
C
      IER(1)=0
C
C  controlla il numero delle righe
C
      IF (N-2) 11,1,1
C
C  controlla il numero totale di elementi della matrice
C
    1 IF(TOPOL(N+1)-1-NTERM) 12,2,12
C
C  controlla che la topologia della matrice sia corretta .
C  controlla che la lunghezza di ciascuna riga sia corretta .
C  calcola NDZ ,cioe' il numero di elementi diagonali nulli .
C
    2 NDZ=0
      DO 10 I=1,N
         IS=TOPOL(I)
         IE=TOPOL(I+1)-1
         IF(IE-IS) 3,7,5
    3    IF(IE-IS+1) 13,4,24
    4    NDZ=NDZ+1
         GO TO 10
    5    IEM1=IE-1
         DO 6 K=IS,IEM1
C
C  salta se un indice di colonna e' minore od uguale allo
C  indice che nell'ordine lo ha preceduto
C
            IF(JA(K+1)-JA(K)) 14,14,6
    6    CONTINUE
    7    JS=JA(IS)
         JE=JA(IE)
         DO 66 K=IS,IE
         IF(JA(K)-I)66,39,66
   66    CONTINUE
         NDZ=NDZ+1
C
C  salta se il primo indice di colonna della riga I non
C  e' maggiore di 1 .
C
   39    IF(JS-1) 15,9,9
C
C  salta se l'ultimo indice di colonna della riga i e' maggiore di N .
C
    9    IF(JE-N)10,10,16
   10 CONTINUE
      GO TO 20
C
C  imposta gli indicatori di errore .
C
   11 IER(1)=1000
      GO TO 20
   12 IER(1)=2000
      GO TO 20
   13 IER(1)=7000
      GO TO 19
   14 IER(1)=6000
      GO TO 18
   15 IER(1)=4000
      GO TO 17
   16 IER(1)=5000
   17 IER(3)=JS
      IER(4)=JE
      GO TO 19
   18 IER(2)=K
      IER(3)=JA(K)
      IER(4)=JA(K+1)
   19 IER(5)=IS
      IER(6)=IE
      IER(7)=I
   20 IF(IER(1))22,23,22
   22 WRITE(IOUT2,100) IER(1)
   23 RETURN
C
C  questo errore non dovrebbe accadere .
C
   24 IER(1)=9061
      GO TO 19
  100 FORMAT( ' ROUTINE CHKNEW  CODICE DI ERRORE =',I5)
      END
