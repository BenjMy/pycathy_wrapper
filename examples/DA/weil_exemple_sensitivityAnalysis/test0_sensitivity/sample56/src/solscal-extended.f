c          parametri aggiunti nelle nuove subroutines
c          ----------------------------------------------
c ncoef == maximum dimension of solution vector
c
c lfil    = integer. The fill-in parameter. Each row of L and
c           each row of U will have a maximum of lfil elements
c           in addition to their original number of nonzero elements.
c           Thus storage can be determined beforehand.
c           lfil must be .ge. 0.
c
c tolilu== tolerance for ilut preconditioning                          
c                                                                     
c meth  == integer - type of preconditioning
c        = 1  ilu(0)
c        = 2  milu
c        = 3  ilut
c
c im    == size of krylov subspace:  should not exceed 50 in this      
c          version (can be reset by changing parameter command for    
c          kmax below)                                                 
c
c coefl, jal,                                                        
c ial   == the  matrix  L in compressed sparse row format:             
c       coefl (1:nl)  = nonzero elements of L stored row-wise in order 
c       jal(1:nl) = corresponding column indices.                      
c       ial(1:n+1) = pointer to beginning of each row in coefl and jal 
c         here nl = number of nonzero elements in L = ial(n+1)-1     

c coefu, jau,                                                          
c iau   == the  matrix  U in compressed sparse row format:             
c       coefu (1:nu)  = nonzero elements of U stored row-wise in order 
c       jau(1:nu) = corresponding column indices.                      
c       iau(1:n+1) = pointer to beginning of each row in coefu and jau 
c         here nu = number of nonzero elements in U = iau(n+1)-1       
c
c nwk     = integer. The minimum length of arrays alu and jlu
c
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju      = integer array of length n containing the pointers to
c           the beginning of each row of U in the matrix alu,jlu.
c
c work arrays:                                                         
c=============                                                         
c vv    == work matrix  (real) of dimensions (ncoef, im)                                     *
c iv    == work matrix  (integer) of dimensions (ncoef, 3)                                     *
c
c ac, jac,                                                             
c iac   == the  matrix complete LU in compressed sparse row format:             *
c       ac (1:nc)  = nonzero elements of LU stored row-wise in order   
c       jac(1:nc) = corresponding column indices.                      
c  topolc(1:nc+1) = pointer to beginning of each row in coefu and jac  
c         here nc = number of nonzero elements in LU = iac(n+1)-1       
c-----------------------------------------------------------------------
c----------------------------------------------------------------------c
      subroutine amux (n, x, y, a,ja,ia)
      real*8  x(*), y(*), a(*), t
      integer ja(*), ia(*)
c-----------------------------------------------------------------------
c         A times a vector
c-----------------------------------------------------------------------
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=Ax
c
c-----------------------------------------------------------------------
c     --- main loop
      do 100 i = 1,n
         k1 = ia(i)
         k2 = ia(i+1) -1
         t = 0.0d0
         do 99 k=k1, k2
            t = t + a(k)*x(ja(k))
 99      continue
         y(i) = t
 100  continue
      return
c---------end-of-amux---------------------------------------------------
c-----------------------------------------------------------------------
      end
c  prodotto fra la matrice coeff e il vettore vett, risultato in r
      subroutine aperbn(n2,topol,ja,coeff,vet,r)
      implicit real*8 (a-h,o-z)
      double precision coeff(*) , vet(*) , r(*)
      integer   ja(*) , topol(*)
      do k = 1 , n2
        r(k) = 0
        do i = topol(k) , topol(k+1)-1
           r(k) = r(k) + coeff(i) * vet(ja(i))
        end do
      end do
      return
      end
C     *** AVINNS *************************************************** ***
C     *                                                                *
C     *   SOLUZIONE MEDIANTE LE SOSTITUZIONI IN AVANTI E               *
C     *   ALL'INDIETRO DI A*X=B CON A MATRICE GENERALE SPARSA          *
C     *   NON SINGOLARE, FATTORIZZATA NEL PRODOTTO L*U ,               *
C     *   IMPIEGANGO PERMUTAZIONI DI RIGHE E DI COLONNE                *
C     *                                                                *
C     ******************************************************************
C
      SUBROUTINE AVINNS(N,IPR,IPC,IC,IU,JC,DI,C,B,X)
      IMPLICIT  REAL*8 (A-H,O-Z)
      DIMENSION IC(1),IU(1),DI(1),JC(1),C(1),B(1),X(1)
      DIMENSION IPR(1),IPC(1)
C
      NM1  =N-1
C
C  PARTE LA SOSTITUZIONE IN AVANTI NELLA  RIGA IPR(1) E COLONNA IPC(1)
C
      I1   =IPR(1)
      I2   =IPC(1)
      X(I2)=B(I1)*DI(1)
C
C     RISOLVE L*Y=B, Y E' MEMORIZZATA IN X
C
      DO 3 I=2,N
C
C    PROCESSA RIGHE E COLONNE NELLA SUCCESSIONE IN CUI ERANO I PIVOTS
C
         I1   =IPR(I)
         I2   =IPC(I)
         IL   =IC(I)
         IH   =IU(I)-1
         SUM  =B(I1)
         IF (IL-IH) 1,1,3
    1    DO 2 J=IL,IH
            JCI  =JC(J)
            JCI  =IPC(JCI)
            AUX  =C(J)
    2       SUM  =SUM-AUX*X(JCI)
    3    X(I2)=SUM*DI(I)
C
C     PARTE LA SOSTITUZIONE ALL'INDIETRO
C     RISOLVE U*X=Y   LA SOLUZIONE VIENE RESTITUITA NEL VETTORE X
C
      DO 6 I=1,NM1
         I1   =N-I
         I2   =IPC(I1)
         SUM  =X(I2)
         IL   =IU(I1)
         IH   =IC(I1+1)-1
         IF (IL-IH) 4,4,6
    4    DO 5 J=IL,IH
            JCI  =JC(J)
            JCI  =IPC(JCI)
            AUX  =C(J)
    5       SUM  =SUM-AUX*X(JCI)
    6    X(I2)=SUM
C
      RETURN
      END
       SUBROUTINE DAXPY (N,T,X,INDX,Y,INDY)
       real*8 x(1), y(1), t
c-------------------------------------------------------------------
c does the following operation
c y <--- y + t * x ,   (replace by the blas routine daxpy )
c indx and indy are supposed to be one here
c-------------------------------------------------------------------
       do 1 k=1,n
       y(k) = y(k) + x(k)*t
1      continue
       return
c*******   end of daxpy  *******************************************
        end
c---------------------------------------------------
       DOUBLE PRECISION FUNCTION DDOT(N,X,IX,Y,IY)
       real*8 x(1), y(1), t
c-------------------------------------------------------------------
c computes the inner product t=(x,y) -- replace by blas routine..
c-------------------------------------------------------------------
       t = 0.0d0
               do 10 j=1,n
               t = t + x(j)*y(j)
10     continue
       ddot=t
       return
c*******   end of ddot   *******************************************
       end
c  calcola i reciproci  degi elementi diagonali di
c  una matrice memorizzata nella
c  forma compressed sparse row format (non symmetric case)
      subroutine diagn(n,topol,ia,ja,coef1,d)
      implicit none
      integer n,topol(*),ia(*),ja(*),k,nterm
      real*8 coef1(*),d(*)
      nterm=topol(n+1)-1
      do 10 k=1,nterm
          if(ia(k).eq.ja(k)) then
          d(ia(k))=1./coef1(K)
          endif
10      continue
        return
        end
C   EMULA LA SUBROUTINE SNRM2 DI SCILIB
C   CALCOLA LA NORMA EUCLIDEA DI UN VETTORE SX
C   DI SOLITO INCX=1
      DOUBLE PRECISION FUNCTION DNRM2(N,SX,INCX)
      REAL*8 SX(1),S
      S=0.
      DO 10 I=1,N
      IX=(I-1)*INCX+1
      S=S+SX(IX)*SX(IX)
10    CONTINUE
      DNRM2 =DSQRT(S)
      RETURN
      END
c-----------------------------------------------------------------------
c   subr. che esegue il passaggio da Doolittle (incomleta)
c   a Crout (incompleta)
      subroutine dolcro(n,ac,jac,topolc,dl)
      implicit none
      real*8 ac(*),dl(*)
      integer n,jac(*),topolc(*),i,j
      do 150 i=1,n
      do 250 j=topolc(i),topolc(i+1)-1
       if(jac(j).lt.i) then
         ac(j)=ac(j)*dl(jac(j))
       endif
       if(jac(j).gt.i) then
        ac(j)=ac(j)/dl(i)
        endif
250     continue
150     continue
        return
        end
C   DATA UNA MATRICE MEMORIZZATA IN CSR FORMAT
C    CALCOLA I PUNTATORI ALLA TRIANGOLARE BASSA E ALTA
C    INOLTRE ESTRAE L ED U
      SUBROUTINE ESTRLU(N,JA,TOPOL,IAL,JAL,IAU,JAU,COEF2,COEFL,
     &                  COEFU)
      IMPLICIT  REAL*8 (A-H,O-Z)
      REAL*8 COEF2(1),COEFL(1),COEFU(1)
      INTEGER JA(1),TOPOL(1),IAL(1),JAL(1),IAU(1),JAU(1)
      IAL(1)=1
      I1=1
      DO 10 I=1,N
      DO 20 J=TOPOL(I),TOPOL(I+1)-1
      IF(JA(J).LE.I) THEN
      JAL(I1)=JA(J)
      COEFL(I1)=COEF2(J)
      I1=I1+1
      ENDIF
20    CONTINUE
      IAL(I+1)=I1
10    CONTINUE
C  TRIANGOLARE ALTA
      IAU(1)=1
      I2=1
      DO 30 I=1,N
      DO 40 J=TOPOL(I),TOPOL(I+1)-1
      IF(JA(J).GE.I) THEN
      JAU(I2)=JA(J)
      COEFU(I2)=COEF2(J)
      I2=I2+1
      ENDIF
40    CONTINUE
      IAU(I+1)=I2
30    CONTINUE
C  CAMBIO GLI ELEMENTI DIAGONALI DI COEFU
      DO 100 I=1,N
      DO 200 J=IAU(I),IAU(I+1)-1
      IF(JAU(J).EQ.I) COEFU(J)=1.
200   CONTINUE
100   CONTINUE
      RETURN
      END
C     ***  FANUNS  **** ************************************************
C     *                                                                *
C     *   FATTORIZZAZIONE NUMERICA DI UNA GENERALE MATRICE SPARSA      *
C     *   NON SINGOLARE USANDO PERMUTAZIONI DI RIGHE E COLONNE         *
C     *                                                                *
C     ******************************************************************
C
      SUBROUTINE FANUNS(N,IPR,IPCI,IA,JA,A,IC,IU,JC,DI,C,IER,X)
      IMPLICIT  REAL*8 (A-H,O-Z)
      DIMENSION IA(1),JA(1),A(1),IC(1),JC(1),C(1),IU(1),DI(1),X(1)
      DIMENSION IPR(1),IPCI(1)
C
      FUZZ =0.D0
      INERR=1
      IF (IER+12345) 2,1,2
    1 INERR=2
    2 IER  =0
      IE   =0
C
C  IL SEGUENTE CICLO CALCOLA I VALORI NUMERICI DELLA RIGA I DI L E DI U
C
      DO 14 I=1,N
         IL   =IE+1
         IE   =IC(I+1)-1
         IH   =IU(I)-1
C
C   PROCESSA LE RIGE NELLA SEQUENZA IPR(1),...,IPR(N) INVECE DI 1,...,N
C
         I1   =IPR(I)
         IAL  =IA(I1)
         IAH  =IA(I1+1)-1
C
C     INIZIALIZZA LO SPAZIO DI LAVORO X CON ZERI NELLE POSIZIONI
C     RISULTANTI (RIGA I DI JC)
C
         X(I) =0.D0
         DO 3 K=IL,IE
            JCI  =JC(K)
    3       X(JCI)=0.D0
C
C     COPIA LA RIGA I DI A, IN ORDINE PERMUTATO, IN X
C
         DO 4 K=IAL,IAH
            JAP  =JA(K)
            JAP  =IPCI(JAP)
    4       X(JAP)=A(K)
C
C     SALTA SE LA RIGA I DI L CONTIENE SOLO L'ELEMENTO DIAGONALE
C
         IF (IL-IH) 5,5,9
C
C     CALCOLA I VALORI NELLA RIGA I DELLA MATRICE TRIANGOLARE BASSA L
C
    5    DO 8 J=IL,IH
            JCI  =JC(J)
            XJI  =X(JCI)
            C(J) =XJI
            JL   =IU(JCI)
            JH   =IC(JCI+1)-1
            IF (JL-JH) 6,6,8
    6       DO 7 K=JL,JH
               JCI  =JC(K)
    7          X(JCI)=X(JCI)-XJI*C(K)
    8       CONTINUE
C
C     MEMORIZZA IL RECIPROCO DELL'ELEMENTO DIAGONALE L(II)
C
    9    XD   =X(I)
         IF (XD-FUZZ) 10,10,11
   10    IF (XD+FUZZ) 11,18,18
   11    XD   =1.D0/XD
         DI(I)=XD
         JL   =IU(I)
C
C     SALTA SE LA RIGA I DI U CONTIENE SOLO L'ELEM. DIAGONALE U(II)=1
C
         IF (JL-IE) 12,12,14
C
C     CALCOLA I VALORI NELLA RIGA I DELLA TRIANGOLARE ALTA U
C
   12    DO 13 J=JL,IE
            JCI  =JC(J)
   13       C(J) =XD*X(JCI)
   14    CONTINUE
   15 RETURN
C
C    IL VALORE NUMERICO DEL PIVOT E' QUASI ZERO
C
   18 IER  =5000
      N    =I-1
      GOTO 15
      END
C     *** FASINS *******************************************************
C     *                                                                *
C     *   FATTORIZZAZIONE SIMBOLICA DI UNA GENERALE MATRICE            *
C     *   SPARSA NON SINGOLARE USANDO PERMUTAZIONI DI RIGHE E DI       *
C     *   COLONNE                                                      *
C     *                                                                *
C     ******************************************************************
C
      SUBROUTINE FASINS(N,IBOT,IPR,IPCI,IA,JA,IC,IU,JC,IK,IER,IX,IP)
      DIMENSION IA(*),JA(*),IC(*),IU(*),JC(*),IX(*),IP(*)
      DIMENSION IPR(*),IPCI(*)
C
      INERR=1
      IF (IER+12345) 2,1,2
    1 INERR=2
    2 IER  =0
      INX  =1
      IK   =0
C
C   IL SEGUENTE CICLO CALCOLA LE POSIZIONI DEGLI ELEMENTI
C   NELL' I-ESIMA RIGA DI L E DI U RISPETTIVAMENTE
C
      IXPS =1
      DO 30 I=1,N
         IC(I)=INX
         IDC  =0
C
C   PROCESSA LE RIGHE NELLA SUCCESSIONE  IPR(1),...,IPR(N)
C   INVECE CHE IN QUELLA 1, 2,....,N
         I1   =IPR(I)
         IAN  =IA(I1+1)-IA(I1)
         IAB  =IA(I1)-1
         IF (IAN) 35,35,3
C
C   COPIA GLI INDICI DI COLONNA PER GLI ELEMENTI NELLA RIGA
C   I1=IPR(I) DA JA IN IX, LI PERMUTA, E LI ORDINA IN ORDINE
C   ASCENDENTE
C   INIZIALIZZA IL VETTORE DI LAVORO IP CON I PUNTATORI 2,3,...,N-1,1
C
    3    DO 4 K=1,IAN
            K1   =K+IAB
            JAP  =JA(K1)
    4       IX(K)=IPCI(JAP)
         CALL SORTNS(IX,IP,1,IAN,IER)
         DO 5 K=1,IAN
    5       IP(K)=K+1
         IP(IAN)=1
         IN   =IAN+1
    6    IX1  =IX(IXPS)
C
C     SALTA AD ALTRE ISTRUZIONI A SECONDA CHE LA POSIZIONE CORRENTE
C     DI PARTENZA IN IX SIA PRIMA, SULLA, O DOPO (ERRORE) LA DIAGONALE
C
         IF (IX1-I) 7,22,37
    7    IL   =IU(IX1)
         IH   =IC(IX1+1)-1
C
C     SALTA SE LA RIGA IX1 DI U CONTIENE SOLO L'ELEMENTO DIAGONALE
C
         IF (IL-IH) 8,8,21
    8    IK   =IK+IH-IL+1
         IXPB =IXPS
         IXP  =IP(IXPB)
         IX1  =IX(IXP)
    9    IC1  =JC(IL)
   10    IF (IX1-IC1) 11,13,15
C
C     L'INDICE DI COLONNA DI IX E' MINORE DELL'INDICE DI COLONNA DI U
C
   11    IXPB =IXP
         IXP  =IP(IXPB)
C
C     SALTA SE IX E' FINITO PRIMA DI U ALTRIMENTI
C     CONFRONTA IL PROSSIMO INDICE DI COLONNA DI IX CON L'INDICE
C     DI COLONNA DI U
C
         IF (IXP-1) 38,17,12
   12    IX1  =IX(IXP)
         GOTO 10
C
C     L'INDICE DI COLONNA DI IX=INDICE DI COLONNA DI U
C     SALTA SE E' RAGGIUNTA LA FINE DELLA RIGA IN U
C
   13    IF (IL-IH) 14,21,38
   14    IL   =IL+1
         IC1  =JC(IL)
         GOTO 11
C
C    L'INDICE DI COLONNA DI IX E' MAGGIORE DELL'INDICE DI COLONNA
C    DI U, AGGIUNGE ELEMENTI DI U A IX IN ORDINE CRESCENTE; INOLTRE
C    AGGIORNA IP
C
   15    IX(IN)=JC(IL)
         IP(IN)=IXP
         IP(IXPB)=IN
C
C     SALTA SE LA RIGA IN U E' FINITA PRIMA DI IX
C
         IF (IL-IH) 16,20,38
   16    IL   =IL+1
         IXPB =IN
         IN   =IN+1
         GOTO 9
C
C     COPIA IL RESTO DI U IN IX
C
   17    IX(IN)=JC(IL)
         IP(IXPB)=IN
         IF (IL-IH) 18,19,38
   18    IL   =IL+1
         IXPB =IN
         IN   =IN+1
         GOTO 17
C
C    AGGIORNA L'INDICATORE DI FINE DI IP E FA PARTIRE L'INDICATORE
C    IXPS
C
   19    IP(IN)=1
   20    IN   =IN+1
   21    IXPS =IP(IXPS)
         IDC  =IDC+1
         GOTO 6
C
C     DIAGONALE RAGGIUNTA
C
   22    IXPS =1
C
C     LO SPAZIO DISPONIBILE PER JC E' SUPERATO SE LA PRIMA POSIZIONE
C     IN  JC+ NO. DI ELEMENTI DI IX E' MAGGIORE DELLA DIMENSIONE
C     SPECIFICATA DI JC
C
         IF (INX+IN-IBOT) 24,24,23
   23    IER  =3000
         GOTO 36
C
C     COPIA IX IN JC
C
   24    IF (IDC) 38,27,25
   25    IH   =INX+IDC-1
         DO 26 K=INX,IH
            JC(K)=IX(IXPS)
   26       IXPS =IP(IXPS)
         INX  =IH+1
   27    IXPS =IP(IXPS)
         IF (IXPS-1) 38,29,28
   28    JC(INX)=IX(IXPS)
         INX  =INX+1
         GOTO 27
C
C     IL PUNTATORE ALL'ELEMENTO DIAGONALE = AL PUNTATORE AL PRIMO
C     ELEMENTO NELLA RIGA + NO. DI ELEMENTI A SINISTRA DELLA DIAGONALE
C
   29    IU(I)=IC(I)+IDC
         IK   =IK+INX-IU(I)+1
   30    CONTINUE
      IC(N+1)=INX
   31 RETURN
C
C     LA LUNGHEZZA DELLA RIGA I DI A E' ZERO, A E' SINGOLARE
C
   35 IER  =1000
   36 IP(1)=IA(I1)
      IP(2)=IA(I1+1)-1
      N    =I-1
      GOTO 31
C
C    IL PIVOTING ESEGUITO IN ACCORDO A IPR PORTA A UNO ZERO SIMBOLICO
C    LA MATRICE PUO' ESSERE SINGOLARE
C
   37 IER  =2000
      GOTO 36
C
C     ERRORE CAUSATO DA SBAGLI NEI DATI DI INGRESSO
C
   38 IER  =4000
      GOTO 36
      END
C
C********************  GCRK  *******************************************
C    GCRK  :    GCR(K)
C        ITK e' il numero di direzioni da salvare
C
      SUBROUTINE GCRK(N,POTEN,TNOTI,X,BP,R,P,IA,JA,COEF1,
     1                COEF2,RES,NTERM,TOPOL,imax,tol,scr1,scr2,
     2                ITK,PPP,BPJ,BETAI,NP,CONTP,NITER,ERR)
      implicit real*8 (a-h,o-z)
      real*8  COEF1(*),COEF2(*)
      REAL*8  POTEN(*),TNOTI(*),X(*),BP(*),R(*),P(*)
      REAL*8  A,BB,C,ALFA,BETA,RES(*),scr1(*),scr2(*)
      INTEGER IA(*),JA(*)
      INTEGER TOPOL(*),CONTP(*)
      REAL*8 A1
      REAL*8 PPP(n,*),BPJ(n,*),BETAI(*),DD,AB
C
      NITER=1
      XLUNG=0.
      DO 200 K=1,N
      X(K)=TNOTI(K)
      P(K)=0.
200   R(K)=0.
      DO I=1,NP
         IN=CONTP(I)
         X(IN)=0.D0
      END DO
      DO I=1,N
         XLUNG=XLUNG+X(I)*X(I)
      END DO
      DO 201 K=1,NTERM
201   R(IA(K))=R(IA(K))+COEF1(K)*POTEN(JA(K))
      DO 202 K=1,N
      R(K)=TNOTI(K)-R(K)
202   CONTINUE
C    CALCOLO L**(-1)*R0
      CALL PRODL(P,R,N,IA,JA,COEF2,NTERM,TOPOL)
      DO 203 K=1,N
       RES(K)=P(K)
       R(K)=0.0
203   CONTINUE
C    CALCOLO Y0=U*X0
      CALL PRU(R,POTEN,N,IA,JA,COEF2,NTERM,TOPOL)
       DO 204 K=1,N
       POTEN(K)=R(K)
204   CONTINUE
16    A=0.
      DO 302 K=1,N
      BP(K)=0.
      X(K)=0.
302   CONTINUE
      CALL PRODBH(BP,P,N,IA,JA,COEF1,COEF2,NTERM,TOPOL,scr1,scr2)
      A=0.
      BB=0.
      DO 303 K=1,N
      A=A+RES(K)*BP(K)
      BB=BB+BP(K)*BP(K)
303   CONTINUE
      IM=MOD(NITER,ITK)
      IF(IM.GT.0) JTK=IM
      IF(IM.EQ.0) JTK=ITK
      DO 515 I=1,N
      PPP(I,JTK)=P(I)
      BPJ(I,JTK)=BP(I)
515   CONTINUE
      ALFA=A/BB
      DO 304 K=1,N
      R(K)=0.
      RES(K)=RES(K)-ALFA*BP(K)
      POTEN(K)=POTEN(K)+ALFA*P(K)
304   CONTINUE
      CALL PRODBH(R,RES,N,IA,JA,COEF1,COEF2,NTERM,TOPOL,scr1,scr2)
      DO 545 J=1,JTK
      C=0.
      DD=0.
      DO 305 K=1,N
      C=C+R(K)*BPJ(K,J)
      DD=DD+BPJ(K,J)*BPJ(K,J)
305   CONTINUE
      BETAI(J)=-C/DD
545   CONTINUE
      DO 306 K=1,N
      AB=0.
      DO 386 J=1,JTK
      AB=AB+BETAI(J)*PPP(K,J)
386   CONTINUE
      P(K)=RES(K)+AB
306   CONTINUE
      A=0.
      DO K=1,N
         A=A+(RES(K))**2
      END DO
      IF(XLUNG.GT.0.D0) THEN
         ERR=DSQRT(A/XLUNG)
      ELSE
         ERR=DSQRT(A/N)
      END IF
      IF ( NITER .LT. IMAX  .AND.  ERR .GT. TOL ) THEN
         NITER=NITER+1
         GO TO 16
      END IF
      DO 520 K=1,N
      R(K)=POTEN(K)
520   POTEN(K)=0.
C     TROVO X=U**(-1)*Y
      CALL PRODU(POTEN,R,N,IA,JA,COEF2,NTERM,TOPOL)
      RETURN
      END
C
C********************  GCRKUT  *******************************************
C    GCRK  :    GCR(K)
C        ITK e' il numero di direzioni da salvare
C
      SUBROUTINE GCRKUT(N,POTEN,TNOTI,X,BP,R,P,JA,COEF1,
     1                RES,NTERM,TOPOL,imax,tol,scr1,scr2,
     2                ITK,PPP,BPJ,BETAI,NP,CONTP,NITER,ERR,
     3     ial,jal,coefl,iau,jau,coefu,dl)                         
      implicit real*8 (a-h,o-z)
      real*8 coefl(*),coefu(*),dl(*)
      integer ial(*),jal(*),iau(*),jau(*)
      real*8  COEF1(*)
      REAL*8  POTEN(*),TNOTI(*),X(*),BP(*),R(*),P(*)
      REAL*8  A,BB,C,ALFA,BETA,RES(*),scr1(*),scr2(*)
      INTEGER JA(*)
      INTEGER TOPOL(*),CONTP(*)
      REAL*8 A1
      REAL*8 PPP(n,*),BPJ(n,*),BETAI(*),DD,AB
C
      write(99,*) 'iter, residual'
      NITER=1
      XLUNG=0.
      DO 200 K=1,N
      X(K)=TNOTI(K)
      P(K)=0.
200   R(K)=0.
      DO I=1,NP
         IN=CONTP(I)
         X(IN)=0.D0
      END DO
      DO I=1,N
         XLUNG=XLUNG+X(I)*X(I)
      END DO
ccc     DO 201 K=1,NTERM
ccc201   R(IA(K))=R(IA(K))+COEF1(K)*POTEN(JA(K))
       call aperbn(n,topol,ja,coef1,poten,r)
      DO 202 K=1,N
      R(K)=TNOTI(K)-R(K)
202   CONTINUE
C    CALCOLO L**(-1)*R0
cc      CALL PRODL(P,R,N,IA,JA,COEF2,NTERM,TOPOL)
       call yprdly(p,r,n,ial,jal,coefl,dl)
      DO 203 K=1,N
       RES(K)=P(K)
       R(K)=0.0
203   CONTINUE
C    CALCOLO Y0=U*X0
cc      CALL PRU(R,POTEN,N,IA,JA,COEF2,NTERM,TOPOL)
        call aperbn(n,iau,jau,coefu,poten,r)
       DO 204 K=1,N
       POTEN(K)=R(K)
204   CONTINUE
16    A=0.
      DO 302 K=1,N
      BP(K)=0.
      X(K)=0.
302   CONTINUE
cc      CALL PRODBH(BP,P,N,IA,JA,COEF1,COEF2,NTERM,TOPOL,scr1,scr2)
        call yprdhg(bp,p,n,ial,jal,coefl,iau,jau,coefu,ja,topol,
     1              coef1,dl,scr1,scr2)
      A=0.
      BB=0.
      DO 303 K=1,N
      A=A+RES(K)*BP(K)
      BB=BB+BP(K)*BP(K)
303   CONTINUE
      IM=MOD(NITER,ITK)
      IF(IM.GT.0) JTK=IM
      IF(IM.EQ.0) JTK=ITK
      DO 515 I=1,N
      PPP(I,JTK)=P(I)
      BPJ(I,JTK)=BP(I)
515   CONTINUE
      ALFA=A/BB
      DO 304 K=1,N
      R(K)=0.
      RES(K)=RES(K)-ALFA*BP(K)
      POTEN(K)=POTEN(K)+ALFA*P(K)
304   CONTINUE
cc      CALL PRODBH(R,RES,N,IA,JA,COEF1,COEF2,NTERM,TOPOL,scr1,scr2)
        call yprdhg(r,res,n,ial,jal,coefl,iau,jau,coefu,ja,topol,
     1              coef1,dl,scr1,scr2)
      DO 545 J=1,JTK
      C=0.
      DD=0.
      DO 305 K=1,N
      C=C+R(K)*BPJ(K,J)
      DD=DD+BPJ(K,J)*BPJ(K,J)
305   CONTINUE
      BETAI(J)=-C/DD
545   CONTINUE
      DO 306 K=1,N
      AB=0.
      DO 386 J=1,JTK
      AB=AB+BETAI(J)*PPP(K,J)
386   CONTINUE
      P(K)=RES(K)+AB
306   CONTINUE
      A=0.
      DO K=1,N
         A=A+(RES(K))**2
      END DO
      IF(XLUNG.GT.0.D0) THEN
         ERR=DSQRT(A/XLUNG)
      ELSE
         ERR=DSQRT(A/N)
      END IF
      IF ( NITER .LT. IMAX  .AND.  ERR .GT. TOL ) THEN
         write(99,*) niter,err
         NITER=NITER+1
         GO TO 16
      END IF
      DO 520 K=1,N
      R(K)=POTEN(K)
520   POTEN(K)=0.
C     TROVO X=U**(-1)*Y
cc      CALL PRODU(POTEN,R,N,IA,JA,COEF2,NTERM,TOPOL)
       call yprduy(poten,r,n,iau,jau,coefu)
      RETURN
      END
C  ---------------------------------
C                GCSTAD
C       VERSIONE CON PRECONDIZIONAMENTO DIAGONALE
C         PER ELABORATORE SCALARE
       SUBROUTINE GCSTAD(N,POTEN,TNOTI,PJM,VJM,P,V,S,T,R,RT,
     1 TOPOL,JA,COEF1,
     2 IMAX,TOL,SCR1,NP,CONTP,NITER,ERR,DD)
        IMPLICIT NONE
        INTEGER JA(*),TOPOL(*)
        INTEGER I, N,IMAX,NITER,K
        REAL*8 POTEN(*),TNOTI(*),PJM(*),VJM(*),P(*),V(*)
        REAL*8 S(*),T(*),R(*),RT(*),BETA
        REAL*8 COEF1(*),SCR1(*),DD(*)
        REAL*8 ROJ,ROJM1,ALFJ,ALFJM1,ETAJ,SIGMA,AA,BB,XLUNG,AZ,ERR
        INTEGER NP,CONTP(*),IN,J
        REAL*8 TOL
C-----------------------
      write(99,*) 'iter, residual'
        niter=1
        xlung=0.0
        do 200 i=1,n
        pjm(i)=tnoti(i)
200     continue
        do i=1,np
         in=contp(i)
         pjm(in)=0.d0
        end do
         do i=1,n
           xlung=xlung+pjm(i)*pjm(i)
         end do
c  calcolo A*x0
        call aperbn(n,topol,ja,coef1,poten,scr1)
c   calcolo r0=b-A*x0
        do 205 k=1,n
         scr1(k)=tnoti(k)-scr1(k)
205     continue
c   calcolo L**(-1)*r0
         do 100 k=1,n
           R(K)=SCR1(K)*DD(K)
100       continue
c  calcolo y0=U*x0=I*x0
          do 96 i=1,n
          RT(I)=R(I)
          pjm(i)=0.d0
          vjm(i)=0.d0
96        continue
          rojm1=1.D0
          alfjm1=1.
          etaj=1.
16        continue
          roj=0.d0
          do 105 i=1,n
          roj=roj+rt(i)*r(i)
105       continue
          beta=(roj/rojm1)*(alfjm1/etaj)
           do 110 i=1,n
            p(i)=r(i)+beta*(pjm(i)-etaj*vjm(i))
110        continue
c    Vj=B*pj
         call aperbn(n,topol,ja,coef1,p,v)
        DO 185 I=1,N
        v(I)=v(I)*dd(i)
185     continue
c---------------
        sigma=0.d0
        do 120 i=1,n
         sigma=sigma+rt(i)*v(i)
120     continue
        alfj=roj/sigma
         do 130 i=1,n
          s(i)=r(i)-alfj*v(i)
130     continue
c    tj=Bsj
         call aperbn(n,topol,ja,coef1,s,t)
       do 195 i=1,n
       t(I)=t(I)*dd(I)
195    continue
c-------------------
       aa=0.d0
       bb=0.d0
       do 140 i=1,n
        aa=aa+s(i)*t(i)
        bb=bb+t(i)*t(i)
140    continue
       etaj=aa/bb
       do 150 i=1,n
        poten(i)=poten(i)+alfj*p(i)+etaj*s(i)
        r(i)=s(i)-etaj*t(i)
150    continue
          do i=1,np
           in=contp(i)
           r(in)=0.d0
            end do
        aa=0.d0
        do 170 i=1,n
         aa=aa+r(i)*r(i)
170     continue
      IF(XLUNG.GT.0.D0) THEN
         ERR=DSQRT(AA/XLUNG)
      ELSE
         ERR=DSQRT(AA/N)
      END IF
c
c  check if the norm is smaller than the tolerance
c
          if((niter.le.imax).and.(err.gt.tol)) then
             write(99,*) niter,err
             niter=niter+1
             rojm1=roj
             alfjm1=alfj
              do 180 i=1,n
              pjm(i)=p(i)
              vjm(i)=v(i)
180           continue
              goto 16
            ENDIF
           return
           end
C  ---------------------------------
C                GCSTAN
C       VERSIONE NON PRECONDIZIONATA PER ELABORATORE SCALARE
       SUBROUTINE GCSTAN(N,POTEN,TNOTI,PJM,VJM,P,V,S,T,R,RT,
     1 TOPOL,JA,COEF1,
     2 IMAX,TOL,SCR1,NP,CONTP,NITER,ERR)
        IMPLICIT NONE
        INTEGER JA(*),TOPOL(*)
        INTEGER I, N,IMAX,NITER,K
        REAL*8 POTEN(*),TNOTI(*),PJM(*),VJM(*),P(*),V(*)
        REAL*8 S(*),T(*),R(*),RT(*),BETA
        REAL*8 COEF1(*),SCR1(*)
        REAL*8 ROJ,ROJM1,ALFJ,ALFJM1,ETAJ,SIGMA,AA,BB,XLUNG,AZ,ERR
        INTEGER NP,CONTP(*),IN
        REAL*8 TOL
      write(99,*) 'iter, residual'
        niter=1
        xlung=0.0
        do 200 i=1,n
        pjm(i)=tnoti(i)
200     continue
        do i=1,np
         in=contp(i)
         pjm(in)=0.d0
        end do
         do i=1,n
           xlung=xlung+pjm(i)*pjm(i)
         end do
c  calcolo A*x0
        call aperbn(n,topol,ja,coef1,poten,scr1)
c   calcolo r0=b-A*x0
        do 205 k=1,n
         scr1(k)=tnoti(k)-scr1(k)
205     continue
c   calcolo L**(-1)*r0
         do 100 k=1,n
           R(K)=SCR1(K)
          rt(k)=r(k)
100       continue
c  calcolo y0=U*x0=I*x0
          do 96 i=1,n
          pjm(i)=0.d0
          vjm(i)=0.d0
96        continue
          rojm1=1.D0
          alfjm1=1.
          etaj=1.
16        continue
          roj=0.d0
          do 105 i=1,n
          roj=roj+rt(i)*r(i)
105       continue
          beta=(roj/rojm1)*(alfjm1/etaj)
           do 110 i=1,n
            p(i)=r(i)+beta*(pjm(i)-etaj*vjm(i))
110        continue
c    Vj=B*pj
         call aperbn(n,topol,ja,coef1,p,v)
        sigma=0.d0
        do 120 i=1,n
         sigma=sigma+rt(i)*v(i)
120     continue
        alfj=roj/sigma
         do 130 i=1,n
          s(i)=r(i)-alfj*v(i)
130     continue
c    tj=Bsj
         call aperbn(n,topol,ja,coef1,s,t)
       aa=0.d0
       bb=0.d0
       do 140 i=1,n
        aa=aa+s(i)*t(i)
        bb=bb+t(i)*t(i)
140    continue
       etaj=aa/bb
       do 150 i=1,n
        poten(i)=poten(i)+alfj*p(i)+etaj*s(i)
        r(i)=s(i)-etaj*t(i)
150    continue
          do i=1,np
           in=contp(i)
           r(in)=0.d0
            end do
        aa=0.d0
        do 170 i=1,n
         aa=aa+r(i)*r(i)
170     continue
      IF(XLUNG.GT.0.D0) THEN
         ERR=DSQRT(AA/XLUNG)
      ELSE
         ERR=DSQRT(AA/N)
      END IF
c
c  check if the norm is smaller than the tolerance
c
          if((niter.le.imax).and.(err.gt.tol)) then
             write(99,*) niter,err
             niter=niter+1
             rojm1=roj
             alfjm1=alfj
              do 180 i=1,n
              pjm(i)=p(i)
              vjm(i)=v(i)
180           continue
              goto 16
            ENDIF
           return
           end
C                GCSTAB
C       VERSIONE PRECONZIZIONATA PER ELABORATORE SCALARE
       SUBROUTINE GCSTAS(N,POTEN,TNOTI,PJM,VJM,P,V,S,T,R,RT,
     1 TOPOL,IA,JA,COEF1,COEF2,NTERM,
     2 IMAX,TOL,SCR1,SCR2,NP,CONTP,NITER,ERR)
        IMPLICIT NONE
        INTEGER JA(*),TOPOL(*)
        INTEGER I, N,IMAX,NITER,K
        REAL*8 POTEN(*),TNOTI(*),PJM(*),VJM(*),P(*),V(*)
        REAL*8 S(*),T(*),R(*),RT(*),BETA
        REAL*8 COEF1(*),SCR1(*),SCR2(*)
        REAL*8 ROJ,ROJM1,ALFJ,ALFJM1,ETAJ,SIGMA,AA,BB,XLUNG,AZ,ERR
        INTEGER IA(*),NP,CONTP(*),NTERM,IN
        REAL*8 TOL,COEF2(*)
      write(99,*) 'iter, residual'
        niter=1
        xlung=0.0
        do 200 i=1,n
        pjm(i)=tnoti(i)
200     continue
        do i=1,np
         in=contp(i)
         pjm(in)=0.d0
        end do
         do i=1,n
           xlung=xlung+pjm(i)*pjm(i)
         end do
c  calcolo A*x0
        call aperbn(n,topol,ja,coef1,poten,scr1)
c   calcolo r0=b-A*x0
        do 205 k=1,n
         scr1(k)=tnoti(k)-scr1(k)
205     continue
c   calcolo L**(-1)*r0
         call prodl(r,scr1,n,ia,ja,coef2,nterm,topol)
         do 100 k=1,n
          rt(k)=r(k)
100       continue
c  calcolo y0=U*x0
          call pru(scr1,poten,n,ia,ja,coef2,nterm,topol)
           do 94  k=1,n
            poten(k)=scr1(k)
94        continue
          do 96 i=1,n
          pjm(i)=0.d0
          vjm(i)=0.d0
96        continue
          rojm1=1.D0
          alfjm1=1.
          etaj=1.
16        continue
          roj=0.d0
          do 105 i=1,n
          roj=roj+rt(i)*r(i)
105       continue
          beta=(roj/rojm1)*(alfjm1/etaj)
           do 110 i=1,n
            p(i)=r(i)+beta*(pjm(i)-etaj*vjm(i))
110        continue
c    Vj=B*pj
      call prodbh(v,p,n,ia,ja,coef1,coef2,nterm,topol,scr1,scr2)
C       STOP
        sigma=0.d0
        do 120 i=1,n
         sigma=sigma+rt(i)*v(i)
120     continue
        alfj=roj/sigma
         do 130 i=1,n
          s(i)=r(i)-alfj*v(i)
130     continue
c    tj=Bsj
      call prodbh(t,s,n,ia,ja,coef1,coef2,nterm,topol,scr1,scr2)
       aa=0.d0
       bb=0.d0
       do 140 i=1,n
        aa=aa+s(i)*t(i)
        bb=bb+t(i)*t(i)
140    continue
       etaj=aa/bb
       do 150 i=1,n
        poten(i)=poten(i)+alfj*p(i)+etaj*s(i)
        r(i)=s(i)-etaj*t(i)
150    continue
          do i=1,np
           in=contp(i)
           r(in)=0.d0
            end do
        aa=0.d0
        do 170 i=1,n
         aa=aa+r(i)*r(i)
170     continue
      IF(XLUNG.GT.0.D0) THEN
         ERR=DSQRT(AA/XLUNG)
      ELSE
         ERR=DSQRT(AA/N)
      END IF
c
c  check if the norm is smaller than the tolerance
c
          if((niter.le.imax).and.(err.gt.tol)) then
             write(99,*) niter,err
ccccciter             write(6,138)niter,err
             niter=niter+1
             rojm1=roj
             alfjm1=alfj
              do 180 i=1,n
              pjm(i)=p(i)
              vjm(i)=v(i)
180           continue
              goto 16
            ENDIF
ccccciter       write(6,138)niter,err
138     format(i5,3x,'err=',e12.5)
           do 520 i=1,n
             r(i)=poten(i)
             poten(i)=0.d0
520      continue
c  trovo x=U**(-1)y
          call produ(poten,r,n,ia,ja,coef2,nterm,topol)
           return
           end
C                GCSTAB
C       VERSIONE PRECONZIZIONATA PER ELABORATORE SCALARE
       SUBROUTINE GCSTUT(N,POTEN,TNOTI,PJM,VJM,P,V,S,T,R,RT,
     1 TOPOL,JA,COEF1,NTERM,
     2 IMAX,TOL,SCR1,SCR2,NP,CONTP,NITER,ERR,
     3     ial,jal,coefl,iau,jau,coefu,dl)                         
        IMPLICIT NONE
      real*8 coefl(*),coefu(*),dl(*)
      integer ial(*),jal(*),iau(*),jau(*)
        INTEGER JA(*),TOPOL(*)
        INTEGER I, N,IMAX,NITER,K
        REAL*8 POTEN(*),TNOTI(*),PJM(*),VJM(*),P(*),V(*)
        REAL*8 S(*),T(*),R(*),RT(*),BETA
        REAL*8 COEF1(*),SCR1(*),SCR2(*)
        REAL*8 ROJ,ROJM1,ALFJ,ALFJM1,ETAJ,SIGMA,AA,BB,XLUNG,AZ,ERR
        INTEGER NP,CONTP(*),NTERM,IN
        REAL*8 TOL
      write(99,*) 'iter, residual'
        niter=1
        xlung=0.0
        do 200 i=1,n
        pjm(i)=tnoti(i)
200     continue
        do i=1,np
         in=contp(i)
         pjm(in)=0.d0
        end do
         do i=1,n
           xlung=xlung+pjm(i)*pjm(i)
         end do
c  calcolo A*x0
        call aperbn(n,topol,ja,coef1,poten,scr1)
c   calcolo r0=b-A*x0
        do 205 k=1,n
         scr1(k)=tnoti(k)-scr1(k)
205     continue
c   calcolo L**(-1)*r0
cc        call prodl(r,scr1,n,ia,ja,coef2,nterm,topol)
           call yprdly(r,scr1,n,ial,jal,coefl,dl)
         do 100 k=1,n
          rt(k)=r(k)
100       continue
c  calcolo y0=U*x0
cc          call pru(scr1,poten,n,ia,ja,coef2,nterm,topol)
           call aperbn(n,iau,jau,coefu,poten,scr1)
           do 94  k=1,n
            poten(k)=scr1(k)
94        continue
          do 96 i=1,n
          pjm(i)=0.d0
          vjm(i)=0.d0
96        continue
          rojm1=1.D0
          alfjm1=1.
          etaj=1.
16        continue
          roj=0.d0
          do 105 i=1,n
          roj=roj+rt(i)*r(i)
105       continue
          beta=(roj/rojm1)*(alfjm1/etaj)
           do 110 i=1,n
            p(i)=r(i)+beta*(pjm(i)-etaj*vjm(i))
110        continue
c    Vj=B*pj
ccc     call prodbh(v,p,n,ia,ja,coef1,coef2,nterm,topol,scr1,scr2)
        call yprdhg(v,p,n,ial,jal,coefl,iau,jau,coefu,ja,topol,
     1              coef1,dl,scr1,scr2)
        sigma=0.d0
        do 120 i=1,n
         sigma=sigma+rt(i)*v(i)
120     continue
        alfj=roj/sigma
         do 130 i=1,n
          s(i)=r(i)-alfj*v(i)
130     continue
c    tj=Bsj
cc      call prodbh(t,s,n,ia,ja,coef1,coef2,nterm,topol,scr1,scr2)
        call yprdhg(t,s,n,ial,jal,coefl,iau,jau,coefu,ja,topol,
     1              coef1,dl,scr1,scr2)
       aa=0.d0
       bb=0.d0
       do 140 i=1,n
        aa=aa+s(i)*t(i)
        bb=bb+t(i)*t(i)
140    continue
       etaj=aa/bb
       do 150 i=1,n
        poten(i)=poten(i)+alfj*p(i)+etaj*s(i)
        r(i)=s(i)-etaj*t(i)
150    continue
          do i=1,np
           in=contp(i)
           r(in)=0.d0
            end do
        aa=0.d0
        do 170 i=1,n
         aa=aa+r(i)*r(i)
170     continue
      IF(XLUNG.GT.0.D0) THEN
         ERR=DSQRT(AA/XLUNG)
      ELSE
         ERR=DSQRT(AA/N)
      END IF
c
c  check if the norm is smaller than the tolerance
c
          if((niter.le.imax).and.(err.gt.tol)) then
             write(99,*) niter,err
ccccciter             write(6,138)niter,err
             niter=niter+1
             rojm1=roj
             alfjm1=alfj
              do 180 i=1,n
              pjm(i)=p(i)
              vjm(i)=v(i)
180           continue
              goto 16
            ENDIF
ccccciter       write(6,138)niter,err
138     format(i5,3x,'err=',e12.5)
           do 520 i=1,n
             r(i)=poten(i)
             poten(i)=0.d0
520      continue
c  trovo x=U**(-1)y
c          call produ(poten,r,n,ia,ja,coef2,nterm,topol)
           call yprduy(poten,r,n,iau,jau,coefu)
           return
           end
C**********************    GRADME  **********************************
      SUBROUTINE GRADDP(N,POTEN,scr,RES,B,R,P,TOPOL,JA,COEF1,
     *COEF2,NTERM,TNOTI,tol,imax,np,contp,NITER,ERR)
c
c  conjugate gradient method
c
c  TNOTI contains the RHS and is never changed
c  R and SCR are scratch vectors that are initialized to TNOTI
c
c  CHANGES from Giuseppe verions:
c               included IMAX and TOL as exchange parameters
c               included NP and CONTP as exchange parameters
c               the norm of the residual is calculated excluding the
c                     Dirichlet boundary conditions
c               the test for convergence is made on the maximum norm AND
c                     the average norm of residual being less than TOL
c
      implicit real*8 (a-h,o-z)
      REAL*8 COEF1(*),COEF2(*)
      REAL*8 POTEN(*),scr(*),RES(*),B(*),R(*),P(*),TNOTI(*)
      INTEGER*4 JA(*),contp(*)
      INTEGER*4 TOPOL(*)
c
      A=0.
      NITER=1
      DO 200 K=1,N
      scr(k)=tnoti(k)
200   R(K)=TNOTI(K)
      DO I=1,NP
         IN=CONTP(I)
         SCR(IN)=0.D0
      END DO
      XLUNG=0.D0
      DO I=1,N
         XLUNG=XLUNG+SCR(I)*SCR(I)
      END DO
      DO I=1,NP
         IN=CONTP(I)
         SCR(IN)=TNOTI(IN)
      END DO
      DO 205 K=1,N
      I=TOPOL(K)
      J=TOPOL(K+1)-1
      DO 205 M=I,J
      R(K)=R(K)-COEF1(M)*POTEN(JA(M))
      res(k)=r(k)
      IF(M.NE.I) then
         R(JA(M))=R(JA(M))-COEF1(M)*POTEN(K)
         res(ja(m))=r(ja(m))
      endif
205   CONTINUE
c
c  takes the dirichelet conditions away from residual
c
      do 202 k=1,np
         res(contp(k))=0.
 202  continue
c
c  calculates norm of residual
c
      a=0.
      do 203 i=1,n
         a=a+res(i)*res(i)
 203  continue
      a=dsqrt(a/n)
      CALL PRODDP(P,R,N,TOPOL,JA,COEF2,NTERM)
16    A=0.
      DO 302 K=1,N
302   B(K)=0.
      DO 10 K=1,N
      I=TOPOL(K)
      J=TOPOL(K+1)-1
      DO 10 M=I,J
      B(K)=B(K)+COEF1(M)*P(JA(M))
      IF(M.NE.I) B(JA(M))=B(JA(M))+COEF1(M)*P(K)
10    CONTINUE
      A=0.
      BB=0.
      DO 303 K=1,N
      A=A+P(K)*R(K)
303   BB=BB+P(K)*B(K)
      ALFA=A/BB
      DO 304 K=1,N
      R(K)=R(K)-ALFA*B(K)
304   POTEN(K)=POTEN(K)+ALFA*P(K)
      CALL PRODDP(scr,R,N,TOPOL,JA,COEF2,NTERM)
      C=0.
      DO 305 K=1,N
305   C=C+B(K)*scr(K)
      BETA=-C/BB
      DO 306 K=1,N
306   P(K)=scr(K)+BETA*P(K)
c
c  takes the dirichelet conditions away from residual
c
      do 7 k=1,n
         res(k)=r(k)
   7  continue
      do 8 k=1,np
         res(contp(k))=0.
   8  continue
c
c  calculates norm of residual
c
      A=0.
      DO K=1,N
         A=A+(RES(K))**2
      END DO
      IF(XLUNG.GT.0.D0) THEN
         ERR=DSQRT(A/XLUNG)
      ELSE
         ERR=DSQRT(A/N)
      END IF
c
c  check if the norm is smaller than the tolerance
c
      if(err.gt.tol .AND. NITER.LT.IMAX) then
         NITER=NITER+1
         GO TO 16
      endif
      RETURN
      END
c
c****************************  GRAMRB  *******************
c
C    GRAMRB E' LA NUOVA VERSIONE DI GRAD1 CON IL MEDODO
C    PRECONDITIONED MINIMUM RESIDUAL METHOD  (MR)
C
      SUBROUTINE GRAMRB(N,POTEN,TNOTI,X,BP,R,P,IA,JA,COEF1,
     1                  COEF2,RES,NTERM,TOPOL,imax,tol,scr1,
     2                  scr2,NP,CONTP,NITER,ERR)
      implicit real*8 (a-h,o-z)
      real*8 COEF1(*),COEF2(*)
      REAL*8 POTEN(*),TNOTI(*),X(*),BP(*),R(*),P(*)
      REAL*8 RES(*),scr1(*),scr2(*)
      INTEGER*4 IA(*),JA(*)
      INTEGER*4 TOPOL(*),CONTP(*)
c
      write(99,*) 'iter, residual'
      NITER=1
      XLUNG=0.
      DO 200 K=1,N
      X(K)=TNOTI(K)
      P(K)=0.
200   R(K)=0.
      DO I=1,NP
         IN=CONTP(I)
         X(IN)=0.D0
      END DO
      DO I=1,N
         XLUNG=XLUNG+X(I)*X(I)
      END DO
      DO 201 K=1,NTERM
201   R(IA(K))=R(IA(K))+COEF1(K)*POTEN(JA(K))
      DO 202 K=1,N
      R(K)=TNOTI(K)-R(K)
202   CONTINUE
c
C    CALCOLO L**(-1)*R0
c
      CALL PRODL(P,R,N,IA,JA,COEF2,NTERM,TOPOL)
      DO 203 K=1,N
       RES(K)=P(K)
       R(K)=0.0
203   CONTINUE
c
C    CALCOLO Y0=U*X0
c
      CALL PRU(R,POTEN,N,IA,JA,COEF2,NTERM,TOPOL)
       DO 204 K=1,N
       POTEN(K)=R(K)
204   CONTINUE
16    A=0.
      DO 302 K=1,N
      R(K)=0.
      X(K)=0.
302   CONTINUE
      CALL PRODBH(R,RES,N,IA,JA,COEF1,COEF2,NTERM,TOPOL,scr1,scr2)
      A=0.
      BB=0.
      DO 303 K=1,N
      A=A+RES(K)*R(K)
      BB=BB+R(K)*R(K)
303   CONTINUE
      ALFA=A/BB
      DO 304 K=1,N
      POTEN(K)=POTEN(K)+ALFA*RES(K)
      RES(K)=RES(K)-ALFA*R(K)
304   CONTINUE
      DO I=1,NP
         IN=CONTP(I)
         RES(IN)=0.D0
      END DO
      A=0.
      DO K=1,N
         A=A+(RES(K))**2
      END DO
      IF(XLUNG.GT.0.D0) THEN
         ERR=DSQRT(A/XLUNG)
      ELSE
         ERR=DSQRT(A/N)
      END IF
      IF (NITER .LT. IMAX  .AND.  ERR .GT. TOL) THEN
         write(99,*) niter,err
         NITER=NITER+1
         GO TO 16
      END IF
      DO K=1,N
         R(K)=POTEN(K)
         POTEN(K)=0.0
      END DO
c
C     TROVO X=U**(-1)*Y
c
      CALL PRODU(POTEN,R,N,IA,JA,COEF2,NTERM,TOPOL)
      RETURN
      END
C
c****************************  GRAMUT  *******************
c
C    PRECONDITIONED MINIMUM RESIDUAL METHOD  (MR)
C
      SUBROUTINE GRAMUT(N,POTEN,TNOTI,X,BP,R,P,JA,COEF1,
     1                  RES,NTERM,TOPOL,imax,tol,scr1,
     2                  scr2,NP,CONTP,NITER,ERR,
     3     ial,jal,coefl,iau,jau,coefu,dl)                         
      implicit real*8 (a-h,o-z)
      real*8 coefl(*),coefu(*),dl(*)
      integer ial(*),jal(*),iau(*),jau(*)
      real*8 COEF1(*)
      REAL*8 POTEN(*),TNOTI(*),X(*),BP(*),R(*),P(*)
      REAL*8 RES(*),scr1(*),scr2(*)
      INTEGER*4 JA(*)
      INTEGER*4 TOPOL(*),CONTP(*)
c
      write(99,*) 'iter, residual'
      NITER=1
      XLUNG=0.
      DO 200 K=1,N
      X(K)=TNOTI(K)
      P(K)=0.
200   R(K)=0.
      DO I=1,NP
         IN=CONTP(I)
         X(IN)=0.D0
      END DO
      DO I=1,N
         XLUNG=XLUNG+X(I)*X(I)
      END DO
cc     DO 201 K=1,NTERM
cc      R(IA(K))=R(IA(K))+COEF1(K)*POTEN(JA(K))
cc201   continue
       call aperbn(n,topol,ja,coef1,poten,r)
      DO 202 K=1,N
      R(K)=TNOTI(K)-R(K)
202   CONTINUE
c
C    CALCOLO L**(-1)*R0
c
cc    CALL PRODL(P,R,N,IA,JA,COEF2,NTERM,TOPOL)
      call yprdly(p,r,n,ial,jal,coefl,dl)
      DO 203 K=1,N
       RES(K)=P(K)
       R(K)=0.0
203   CONTINUE
c
C    CALCOLO Y0=U*X0
c
cc      CALL PRU(R,POTEN,N,IA,JA,COEF2,NTERM,TOPOL)
       call aperbn(n,iau,jau,coefu,poten,r)
       DO 204 K=1,N
       POTEN(K)=R(K)
204   CONTINUE
16    A=0.
      DO 302 K=1,N
      R(K)=0.
      X(K)=0.
302   CONTINUE
cc    CALL PRODBH(R,RES,N,IA,JA,COEF1,COEF2,NTERM,TOPOL,scr1,scr2)
        call yprdhg(r,res,n,ial,jal,coefl,iau,jau,coefu,ja,topol,
     1              coef1,dl,scr1,scr2)
      A=0.
      BB=0.
      DO 303 K=1,N
      A=A+RES(K)*R(K)
      BB=BB+R(K)*R(K)
303   CONTINUE
      ALFA=A/BB
      DO 304 K=1,N
      POTEN(K)=POTEN(K)+ALFA*RES(K)
      RES(K)=RES(K)-ALFA*R(K)
304   CONTINUE
      DO I=1,NP
         IN=CONTP(I)
         RES(IN)=0.D0
      END DO
      A=0.
      DO K=1,N
         A=A+(RES(K))**2
      END DO
      IF(XLUNG.GT.0.D0) THEN
         ERR=DSQRT(A/XLUNG)
      ELSE
         ERR=DSQRT(A/N)
      END IF
      IF (NITER .LT. IMAX  .AND.  ERR .GT. TOL) THEN
         write(99,*) niter,err
         NITER=NITER+1
         GO TO 16
      END IF
      DO K=1,N
         R(K)=POTEN(K)
         POTEN(K)=0.0
      END DO
c
C     TROVO X=U**(-1)*Y
c
cc      CALL PRODU(POTEN,R,N,IA,JA,COEF2,NTERM,TOPOL)
       call yprduy(poten,r,n,iau,jau,coefu)
      RETURN
      END
c----------------------------------------------------------------------
        subroutine ilu0 (n, a, ja, ia, alu, jlu, ju, iw, ierr)
        implicit real*8 (a-h,o-z)
        real*8 a(*), alu(*)
        integer ja(*), ia(*), ju(*), jlu(*), iw(*)
c------------------ right preconditioner ------------------------------*
c                    ***   ilu(0) preconditioner.   ***                *
c----------------------------------------------------------------------*
c Note that this has been coded in such a way that it can be used
c with pgmres. Normally, since the data structure of a, ja, ia is
c the same as that of a, ja, ia, savings can be made. In fact with
c some definitions (not correct for general sparse matrices) all we
c need in addition to a, ja, ia is an additional diagonal.
c Ilu0 is not recommended for serious problems. It is only provided
c here for comparison purposes.
c-----------------------------------------------------------------------
c
c on entry:
c---------
c n       = dimension of matrix
c a, ja,
c ia      = original matrix in compressed sparse row storage.
c
c on return:
c-----------
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju	  = pointer to the diagonal elements in alu, jlu.
c	
c ierr	  = integer indicating error code on return
c	     ierr = 0 --> normal return
c	     ierr = k --> code encountered a zero pivot at step k.
c work arrays:
c-------------
c iw	    = integer work array of length n.
c------------
c IMPORTANT
c-----------
c it is assumed that the the elements in the input matrix are stored
c    in such a way that in each row the lower part comes first and
c    then the upper part. To get the correct ILU factorization, it is
c    also necessary to have the elements of L sorted by increasing
c    column number. It may therefore be necessary to sort the
c    elements of a, ja, ia prior to calling ilu0. This can be
c    achieved by transposing the matrix twice using csrcsc.
c
c-----------------------------------------------------------------------
        ju0 = n+2
        jlu(1) = ju0
c
c initialize work vector to zero's
c
        do  31 i=1, n
           iw(i) = 0
 31     continue
c
c main loop
c
       do 500 ii = 1, n
           js = ju0
c
c generating row number ii of L and U.
c
           do 100 j=ia(ii),ia(ii+1)-1
c
c     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
c
              jcol = ja(j)
              if (jcol .eq. ii) then
                 alu(ii) = a(j)
                 iw(jcol) = ii
                 ju(ii)  = ju0
              else
                 alu(ju0) = a(j)
                 jlu(ju0) = ja(j)
                 iw(jcol) = ju0
                 ju0 = ju0+1
              endif
 100       continue
           jlu(ii+1) = ju0
           jf = ju0-1
           jm = ju(ii)-1
c
c exit if diagonal element is reached.
c
           do 150 j=js, jm
              jrow = jlu(j)
              tl = alu(j)*alu(jrow)
              alu(j) = tl
c
c perform  linear combination
c
              do 140 jj = ju(jrow), jlu(jrow+1)-1
                 jw = iw(jlu(jj))
                 if (jw .ne. 0) alu(jw) = alu(jw) - tl*alu(jj)
 140          continue
 150       continue
c
c     invert  and store diagonal element.
c
           if (alu(ii) .eq. 0.0d0) goto 600
           alu(ii) = 1.0d0/alu(ii)
c
c     reset pointer iw to zero
c
           iw(ii) = 0
           do 201 i = js, jf
 201          iw(jlu(i)) = 0
 500       continue
           ierr = 0
           return
c
c     zero pivot :
c
 600       ierr = ii
c
           return
c------- end of ilu0 ---------------------------------------------------
c-----------------------------------------------------------------------
           end
c-----------------------------------------------------------------------
       subroutine ilut (n,a,ja,ia,lfil,tol,alu,jlu,ju,iwk,
     *                  wu,wl,jr,jwl,jwu,ierr)
c-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       real*8 a(*), alu(*), wu(n), wl(n), tol
       integer ja(*),ia(n+1),jlu(*),ju(n),jr(n), jwu(n),
     *      jwl(n), n, lfil, iwk, ierr
c----------------------------------------------------------------------*
c                      *** ILUT preconditioner ***                     *
c                      ---------------------------                     *
c      incomplete LU factorization with dual truncation mechanism      *
c      VERSION 2 : sorting  done for both L and U.                     *
c                                                                      *
c----------------------------------------------------------------------*
c---- coded by Youcef Saad May, 5, 1990. ------------------------------*
c---- Dual drop-off strategy works as follows.                         *
c                                                                      *
c     1) Theresholding in L and U as set by tol. Any element whose size*
c        is less than some tolerance (relative to the norm of current  *
c        row in u) is dropped.                                         *
c                                                                      *
c     2) Keeping only the largest lfil+il(i) elements in the i-th row  *
c        of L and the largest lfil+iu(i) elements in the i-th row of   *
c        U where il(i), iu(i) are the original number of nonzero       *
c        elements of the L-part and the U-part of the i-th row of A    *
c                                                                      *
c Flexibility: one can use tol=0 to get a strategy based on keeping the*
c largest elements in each row of L and U. Taking tol .ne. 0 but lfil=n*
c will give the usual threshold strategy (however, fill-in is then     *
c impredictible).                                                      *
c                                                                      *
c----------------------------------------------------------------------*
c PARAMETERS
c-----------
c
c on entry:
c==========
c n       = integer. The dimension of the matrix A.
c
c a,ja,ia = matrix stored in Compressed Sparse Row format.
c
c lfil    = integer. The fill-in parameter. Each row of L and
c           each row of U will have a maximum of lfil elements
c           in addition to their original number of nonzero elements.
c           Thus storage can be determined beforehand.
c           lfil must be .ge. 0.
c
c iwk     = integer. The minimum length of arrays alu and jlu
c
c On return:
c===========
c
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju      = integer array of length n containing the pointers to
c           the beginning of each row of U in the matrix alu,jlu.
c
c ierr    = integer. Error message with the following meaning.
c           ierr  = 0    --> successful return.
c           ierr .gt. 0  --> zero pivot encountered at step number ierr.
c           ierr  = -1   --> Error. input matrix may be wrong.
c                            (The elimination process has generated a
c                            row in L or U whose length is .gt.  n.)
c           ierr  = -2   --> The matrix L overflows the array al.
c           ierr  = -3   --> The matrix U overflows the array alu.
c           ierr  = -4   --> Illegal value for lfil.
c           ierr  = -5   --> zero row encountered.
c
c work arrays:
c=============
c jr,jwu,jwl 	  = integer work arrays of length n.
c wu, wl          = real work arrays of length n+1, and n resp.
c
c Notes:
c ------
c A must have all nonzero diagonal elements.
c-----------------------------------------------------------------------
        if (lfil .lt. 0) goto 998
c-------------------------------
c initialize ju0 (points to next element to be added to alu,jlu)
c and pointer.
c
        ju0 = n+2
        jlu(1) = ju0
c
c  integer double pointer array.
c
        do 1 j=1, n
         jr(j)  = 0
 1           continue
c-----------------------------------------------------------------------
c  beginning of main loop.
c-----------------------------------------------------------------------
           do 500 ii = 1, n
           j1 = ia(ii)
           j2 = ia(ii+1) - 1
           tnorm = 0.0d0
           do 501 k=j1,j2
              tnorm = tnorm+abs(a(k))
 501          continue
              if (tnorm .eq. 0.0) goto 999
              tnorm = tnorm/real(j2-j1+1)
c
c--- unpack L-part and U-part of row of A in arrays wl, wu --
c
       lenu = 1
       lenl = 0
       jwu(1) = ii
       wu(1) = 0.0
       jr(ii) = 1
c
       do 170  j = j1, j2
           k = ja(j)
           t = a(j)
           if (abs(t) .lt. tol*tnorm .and. k .ne. ii) goto 170
           if (k .lt. ii) then
              lenl = lenl+1
              jwl(lenl) = k
              wl(lenl) = t
              jr(k) = lenl
           else if (k .eq. ii) then
              wu(1) = t
           else
              lenu = lenu+1
              jwu(lenu) = k
              wu(lenu) = t
              jr(k) = lenu
           endif
 170    continue
        tnorm = tnorm/real(j2-j1+1)
        lenl0 = lenl
        lenu0 = lenu
        jj = 0
        nl = 0
c-------------------------------------------------------------------
c---------------------- eliminate previous rows --------------------
c-------------------------------------------------------------------
 150    jj = jj+1
        if (jj .gt. lenl) goto 160
c-------------------------------------------------------------------
c in order to do the elimination in the correct order we need to
c exchange the current row number with the one that has
c smallest column number, among jj,jj+1,...,lenl.
c-------------------------------------------------------------------
        jrow = jwl(jj)
        k = jj
c
c determine smallest column index
c
        do 151 j=jj+1,lenl
           if (jwl(j) .lt. jrow) then
              jrow = jwl(j)
              k = j
           endif
 151    continue
c
c exchange in jwl
c
       if (k .ne. jj) then
           j = jwl(jj)
           jwl(jj) = jwl(k)
           jwl(k) = j
c
c exchange in jr
c
           jr(jrow) = jj
           jr(j) = k
c
c exchange in wl
c
           s = wl(jj)
           wl(jj) = wl(k)
           wl(k) = s
        endif
c
        if (jrow .ge. ii) goto 160
c---------get the multiplier for row to be eliminated: jrow
        fact = wl(jj)*alu(jrow)
c zero out element in row by setting jr(jrow) = 0
        jr(jrow) = 0
        if (abs(fact)*wu(n+2-jrow) .le. tol*tnorm) goto 150
c-------------------------------------------------------------------
c------------ combine current row and row jrow ---------------------
c-------------------------------------------------------------------
        do 203 k = ju(jrow), jlu(jrow+1)-1
           s = fact*alu(k)
           j = jlu(k)
           jpos = jr(j)
c
c if fill-in element is small then disregard:
c
           if (abs(s) .lt. tol*tnorm .and. jpos .eq. 0) goto 203
           if (j .ge. ii) then
c
c     dealing with upper part.
c
              if (jpos .eq. 0) then
c     this is a fill-in element
                 lenu = lenu+1
                 if (lenu .gt. n) goto 995
                 jwu(lenu) = j
                 jr(j) = lenu
                 wu(lenu) = - s
              else
c     no fill-in element --
                 wu(jpos) = wu(jpos) - s
              endif
           else
c
c     dealing with lower part.
c
              if (jpos .eq. 0) then
c     this is a fill-in element
                 lenl = lenl+1
                 if (lenl .gt. n) goto 995
                 jwl(lenl) = j
                 jr(j) = lenl
                 wl(lenl) = - s
              else
c     no fill-in element --
                 wl(jpos) = wl(jpos) - s
              endif
           endif
 203    continue
        nl = nl+1
        wl(nl) = fact
        jwl(nl)  = jrow
        goto 150
c----------------------------------------------------------
c------------ update l-matrix -----------------------------
c----------------------------------------------------------
 160    len = min0(nl,lenl0+lfil)
c 160    len = min0(nl,lfil)
 
        call qsplit (wl,jwl,nl,len)
c
        do 204 k=1, len
           if (ju0 .gt. iwk) goto 996
           alu(ju0) =  wl(k)
           jlu(ju0) =  jwl(k)
c          ju0 = ju0+1
          ju0=ju0+1
204       continue
c
c  save pointer to beginning of row ii of U
c
        ju(ii) = ju0
c
c  reset double-pointer jr to zero (L-part - except first
c  jj-1 elements which have already been reset)
        do 306 k= jj, lenl
              jr(jwl(k)) = 0
 306     continue
        len = min0(lenu,lenu0+lfil)
c        len = min0(lenu,lfil)
         call qsplit (wu(2), jwu(2), lenu-1,len)
c----------------------------------------------------------
c------------ update u-matrix -----------------------------
c----------------------------------------------------------
        t = 0.0d0
        if (len + ju0 .gt. iwk) goto 997
        do 302 k=2, len
           jlu(ju0) = jwu(k)
           alu(ju0) = wu(k)
           t = t + abs(wu(k) )
           ju0 = ju0+1
 302   continue
c
c     save norm in wu (backwards). Norm is in fact average abs value
c
        wu(n+2-ii) = t / real(len+1)
c
c     store inverse of diagonal element of u
c
        if (wu(1) .eq. 0.0) wu(1) = (0.0001 + tol)*tnorm
c
        alu(ii) = 1.0d0/ wu(1)
c
c     update pointer to beginning of next row of U.
c
       jlu(ii+1) = ju0
c
c     reset double-pointer jr to zero (U-part)
c
       do 308 k=1, lenu
           jr(jwu(k)) = 0
 308    continue
c-----------------------------------------------------------------------
c     end main loop
c-----------------------------------------------------------------------
 500    continue
        ierr = 0
        return
c
c     zero pivot :
c
c 900    ierr = ii
c        return
c
c     incomprehensible error. Matrix must be wrong.
c
 995    ierr = -1
        return
c
c     insufficient storage in L.
c
 996    ierr = -2
        return
c
c     insufficient storage in U.
c
 997    ierr = -3
        return
c
c     illegal lfil entered.
c
 998    ierr = -4
        return
c
c     zero row encountered
c
 999    ierr = -5
        return
c---------------- end of ilut  -----------------------------------------
c-----------------------------------------------------------------------
        end
c      METH=1     ILU(0)
c      METH=2     MILU
c      METH=3     ILUT
C*****************************************************************
      SUBROUTINE ILUT1(meth,N,NMAX,A,JA,TOPOL,LFIL,TOLILU,ALU,JLU,JU,
     1NWK,VV,IW,IERR)
      IMPLICIT NONE
      INTEGER NMAX
      integer meth,lfil,nwk,im,k,i,nterm,ierr,n,imax
      INTEGER TOPOL(*),JA(*)
      integer jlu(*),ju(*),iw(n,3)
      real*8 a(*),alu(*)
      REAL*8 VV(N,2)
      real*8 tolilu
C-----------------------------------
C
      If (meth.eq.1) then
       call ilu0(n,a,ja,topol,alu,jlu,ju,iw,ierr)
      endif
      if (meth.eq.2) then
       call milu0(n,a,ja,topol,alu,jlu,ju,iw,ierr)
        endif
       if(meth.eq.3) then
        call ilut(n,a,ja,topol,lfil,tolilu,alu,jlu,ju,nwk,vv,
     1      vv(1,2),iw,iw(1,2),iw(1,3),ierr)
        endif
       if(ierr.ne.0) then
       WRITE(6,*)' IERR=',IERR
       endif
           return
            end
c      METH=1     ILU(0)
c      METH=2     MILU
c      METH=3     ILUT
C*****************************************************************
      subroutine ilutme(meth,n,nmax,a,ja,topol,lfil,tolilu,alu,jlu,
     1ju,nwk,vv,iw,ierr,ac,jac,topolc,ial,jal,coefl,iau,jau,coefu,dl)
      IMPLICIT NONE
      INTEGER NMAX
      integer meth,lfil,nwk,im,k,i,nterm,ierr,n,imax
      INTEGER TOPOL(*),JA(*)
      integer jlu(*),ju(*),iw(n,3)
      real*8 a(*),alu(*)
      REAL*8 VV(N,2)
      real*8 tolilu
C--------------------------------------
       REAL*8 AC(*),COEFL(*),COEFU(*),DL(*)
       INTEGER JAC(*),TOPOLC(*),JAL(*),JAU(*)
       INTEGER IAL(*),IAU(*),J
C-----------------------------------
C
      If (meth.eq.1) then
       call ilu0(n,a,ja,topol,alu,jlu,ju,iw,ierr)
      endif
      if (meth.eq.2) then
       call milu0(n,a,ja,topol,alu,jlu,ju,iw,ierr)
        endif
       if(meth.eq.3) then
        call ilut(n,a,ja,topol,lfil,tolilu,alu,jlu,ju,nwk,vv,
     1      vv(1,2),iw,iw(1,2),iw(1,3),ierr)
        endif
        if(ierr.ne.0) then
       WRITE(6,*)' IERR=',IERR
       endif
C        WRITE(6,119)(JU(I),I=1,N+1)
119      format(' vett ju'/(16i5))
C         WRITE(6,121)(JLU(I),ALU(I),I=1,JU(N))
c------------------------
c   ripristino gli elementi diagonali di l
c   perche' in ILUT c'erano i reciproci
       do 700 i=1,n
       ALU(I)=1./ALU(I)
       DL(I)=ALU(I)
700    continue
C        WRITE(6,119)(JU(I),I=1,N+1)
C         WRITE(6,121)(JLU(I),ALU(I),I=1,JU(N))
121      format('  jlu - alu'/ (4(i5,e12.5)))
c   Da  MSR a CSR
        CALL MSRCSR(N,ALU,JLU,AC,JAC,TOPOLC,VV)
c  passaggio da Doolittle (incompleta) a Crout (incompleta)
        call dolcro(n,ac,jac,topolc,dl)
C         WRITE(6,127)(JAC(I),AC(I),I=1,17)
127      FORMAT(' JAC -  AC'/ (4(I5,E12.5)))
c  Estrae L  ed U e pone uguali a 1 gli elementi diagonali di U
      call estrlu(n,jac,topolc,ial,jal,iau,jau,ac,coefl,coefu)
c---------------
        CALL ORDJA(N,IAL,JAL,COEFL)
        CALL ORDJA(N,IAU,JAU,COEFU)
C        WRITE(6,97)(IAL(I),I=1,N)
97       FORMAT(' IAL'/(10I5))
C        WRITE(6,98)(JAL(I),COEFL(I),I=1,IAL(N+1)-1)
98       FORMAT(' JAL - COEFL'/(4(I5,1X,E12.5)))
C-----------------------
C        WRITE(6,397)(IAU(I),I=1,N,50)
397       FORMAT(' IAU'/(10I8))
C        WRITE(6,198)(JAU(I),COEFU(I),I=1,IAU(N+1)-1,200)
198       FORMAT(' JAU - COEFU'/(4(I5,1X,E12.5)))
C          WRITE(6,713)(DL(I),I=1,N)
713        FORMAT(' DL'/(5E12.5))
           return
            end
c
c  Conjugate gradient method with CHOLESKY preconditioner
c
c      DOUBLE PRECISION
c
C**********************   INCLME  ******************************
      SUBROUTINE INCLDP(iout,JA,TOPOL,COEFR,LRIGA,N,NTERM)
c
c  IOUT is used for the output device
c
      implicit real*8 (a-h,o-z)
c
      REAL*8 LRIGA(*),COEFR(*)
      INTEGER*4 JA(*)
      INTEGER*4 TOPOL(*)
      integer iout
      DO 10 K=1,NTERM
10    LRIGA(K)=0.
      DO 1 KK=1,N
      K=TOPOL(KK)
      A=COEFR(K)-LRIGA(K)
      IF(A.LE.0.0d0) WRITE(iout,100) KK ,A
      IF(A.LE.0.0d0) WRITE(iout,101) LRIGA(TOPOL(KK-1))
C     IF(A.LE.0.0d0) A=2
      IF(A.LE.0.0d0) A=(LRIGA(TOPOL(KK-1)))**2
      LRIGA(K)=DSQRT(A)
      IF(KK.EQ.N) GO TO 1
      I=TOPOL(KK)+1
      J=TOPOL(KK+1)-1
      IF(J.LT.I) GO TO 1
      DO 2 K1=I,J
2     LRIGA(K1)=(COEFR(K1)-LRIGA(K1))/LRIGA(K)
C     IF(I.EQ.J) GO TO 1
      DO 4 K2=I,J
      I1=K2+1
      J1=TOPOL(JA(K2))
      LRIGA(J1)=LRIGA(J1)+LRIGA(K2)**2
      IF(K2.EQ.J) GO TO 4
      J1=J1+1
      IF(J1.GE.TOPOL(JA(K2)+1)) GO TO 4
      GO TO 7
6     I1=I1+1
      IF(I1.GT.J) GO TO 4
      GO TO 7
3     J1=J1+1
      IF(J1.GE.TOPOL(JA(K2)+1)) GO TO 4
7     IF(JA(J1)-JA(I1))3,5,6
5     LRIGA(J1)=LRIGA(J1)+LRIGA(K2)*LRIGA(I1)
      GO TO 6
4     CONTINUE
1     CONTINUE
12    RETURN
100   FORMAT(1X,'ELEMENTO DIAGONALE DI L NULLO,I,J =',I5,2X,E16.5)
101   FORMAT(1X,'ELEMENTO DIAGONALE PRECEDENTE =  ',E16.8)
      END
c
c************************  INCLU  ***************************
c
      SUBROUTINE INCLU(iout,IA,JA,TOPOL,COEFR,
     1                 COEFC,LRIGA,LCOLON,N,NTERM)
      implicit real*8 (a-h,o-z)
      real*8 LRIGA(*),LCOLON(*),COEFR(*),COEFC(*)
      INTEGER*4 IA(*),JA(*)
      INTEGER*4 TOPOL(*)
c
      I=1
      TOPOL(1)=1
      DO 53 K=2,NTERM
      IF(IA(K).GT.IA(K-1)) I=I+1
      IF(IA(K).GT.IA(K-1)) TOPOL(I)=K
53    CONTINUE
      DO 50 K=1,NTERM
      I=TOPOL(JA(K))
51    IF(JA(I).LT.IA(K)) I=I+1
      IF(JA(I).EQ.IA(K)) GO TO 52
      GO TO 51
52    COEFC(I)=COEFR(K)
50    CONTINUE
      K=1
1     LCOLON(K)=COEFC(K)
      LRIGA(K)=COEFR(K)/COEFR(1)
      IF(K.EQ.1) LRIGA(1)=COEFR(1)
      K=K+1
      IF(IA(K).EQ.IA(K-1)) GO TO 1
11    I=K-1
3     J=TOPOL(JA(I))
5     IF(JA(J).LT.IA(I)) GO TO 4
      LRIGA(J)=LCOLON(I)
      LCOLON(J)=LRIGA(I)
      GO TO 2
4     J=J+1
      GO TO 5
2     I=I-1
      IF(JA(I).GT.IA(I)) GO TO 3
7     IF(IA(K).LE.JA(K)) GO TO 10
      K=K+1
      GO TO 7
10    KK=K
6     J=TOPOL(JA(K))
      I=TOPOL(IA(K))
      A=0.
      B=0.
9     IF(JA(I).GE.IA(K).OR.JA(J).GE.IA(K)) GO TO 8
      IF(JA(I)-JA(J)) 20,21,22
20    I=I+1
      GO TO 9
21    A=A+LCOLON(I)*LRIGA(J)
      IF(I.NE.J)B=B+LRIGA(I)*LCOLON(J)
      I=I+1
      J=J+1
      IF(I.GT.NTERM)GO TO 8
      GO TO 9
22    J=J+1
      GO TO 9
8     LCOLON(K)=COEFC(K)-A
      IF(K.EQ.KK) GO TO 13
14    IF(IA(K).LT.JA(K)) LRIGA(K)=(COEFR(K)-B)/LRIGA(KK)
      IF(IA(K).EQ.JA(K)) LRIGA(K)=LCOLON(K)
      K=K+1
      IF(K.GT.NTERM) GO TO 12
      IF(IA(K).EQ.IA(K-1)) GO TO 6
      GO TO 11
13    IF(LCOLON(K).NE.0.0d0)GO TO 14
      WRITE(iout,100) IA(KK)
      KKK=KK
15    KKK=KKK-1
      IF(IA(KKK).EQ.JA(KKK)) LCOLON(K)=LCOLON(KKK)
      IF(IA(KKK).EQ.JA(KKK)) GO TO 14
      GO TO 15
12    RETURN
100   FORMAT(1X,'ELEMENTO DIAGONALE DI L NULLO,I,J =',I5)
      END
C
C (17/03/91 -- REMOVAL OF 'WRITE' STATEMENTS -- C. PANICONI)
C
C     *** INPENS *******************************************************
C     *                                                                *
C     *   CALCOLA IL VETTORE PERMUTAZIONE INVERSA O SE OPT='2' ,       *
C     *   I VETTORI TRASPOSTI DEL DATO E LE PERMUTAZIONI INVERSE       *
C     *                                                                *
C     ******************************************************************
C
      SUBROUTINE INPENS(IP,N,IOPT,INP,IER)
      DIMENSION IP(1),INP(1)
C
      INDER=IER
      IER  =0
      IF (N) 1,1,2
    1 IER  =1000
      GOTO 11
    2 DO 3 J=1,N
    3    INP(J)=0
      DO 7 J=1,N
         IIP  =IP(J)
         IF (IIP) 6,6,4
    4    IF (N-IIP) 6,5,5
    5    IF (INP(IIP)) 6,7,6
    6    IER  =2000
         GOTO 11
    7    INP(IIP)=J
      IF (IOPT-2) 11,8,11
    8 DO 9 J=1,N
         IJ   =IP(J)
         NJ   =INP(J)
         IP(NJ)=IJ
    9    INP(IJ)=NJ
   11 RETURN
      END
C     ***************************
      SUBROUTINE  IORDY(N,VET,AV)
C     ***********************************************************
C     SUBROUTINE ORDINAMENTO DI UN VETTORE INTERO E IN MANIERA
C     COERENTE UNVETTORE REAL AV
C     VARIABILI DI SCAMBIO  VET,N : INOUT
C     ***********************************************************
      INTEGER VET(1),VEL
      REAL*8 AV(1)
C     WRITE(6,28)(VET(I),I=1,N)
28    FORMAT(' VETT '/(10(I5)))
      DO 80 L=2,N
      DO 64 I=L-1,1,-1
      IF(VET(I).GT.VET(L)) GOTO 64
      M=I+2
      GOTO 68
   64 CONTINUE
      M=2
   68 IF(M.GT.L)GOTO 80
      VEL=VET(L)
      A1=AV(L)
      DO73 J=L,M,-1
      VET(J)=VET(J-1)
      AV(J)=AV(J-1)
   73 CONTINUE
      VET(M-1)=VEL
      AV(M-1)=A1
   80 CONTINUE
C     WRITE(6,78)(VET(I),I=1,N)
78    FORMAT(' VETT ORDINATO'/(10(I5)))
C     WRITE(6,79)(AV(I),I=1,N)
79    FORMAT(' VETTORE REAL'/(6(2X,F10.2)))
      RETURN
      END
c-----------------------------------------------------------------------
          subroutine lusol0 (n, y, x, alu, jlu, ju)
          real*8 x(n), y(n), alu(*)
          integer n, jlu(*), ju(*)
c-----------------------------------------------------------------------
c
c performs a forward followed by a backward solve
c for LU matrix as produced by  ILUT
c
c-----------------------------------------------------------------------
c local variables
c
        integer i,k
c
c forward solve
c
        do 40 i = 1, n
           x(i) = y(i)
           do 41 k=jlu(i),ju(i)-1
              x(i) = x(i) - alu(k)* x(jlu(k))
 41        continue
 40     continue
c
c     backward solve.
c
       do 90 i = n, 1, -1
          do 91 k=ju(i),jlu(i+1)-1
              x(i) = x(i) - alu(k)*x(jlu(k))
 91    continue
           x(i) = alu(i)*x(i)
 90     continue
c
        return
c----------------end of lusol0 -----------------------------------------
c-----------------------------------------------------------------------
        end
c----------------------------------------------------------------------
           subroutine milu0 (n, a, ja, ia, alu, jlu, ju, iw, ierr)
           implicit real*8 (a-h,o-z)
           real*8 a(*), alu(*)
           integer ja(*), ia(*), ju(*), jlu(*), iw(*)
c----------------------------------------------------------------------*
c                *** simple milu(0) preconditioner. ***                *
c----------------------------------------------------------------------*
c Note that this has been coded in such a way that it can be used
c with pgmres. Normally, since the data structure of a, ja, ia is
c the same as that of a, ja, ia, savings can be made. In fact with
c some definitions (not correct for general sparse matrices) all we
c need in addition to a, ja, ia is an additional diagonal.
c Ilu0 is not recommended for serious problems. It is only provided
c here for comparison purposes.
c-----------------------------------------------------------------------
c
c on entry:
c----------
c n       = dimension of matrix
c a, ja,
c ia      = original matrix in compressed sparse row storage.
c
c on return:
c----------
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju	  = pointer to the diagonal elements in alu, jlu.
c	
c ierr	  = integer indicating error code on return
c	     ierr = 0 --> normal return
c	     ierr = k --> code encountered a zero pivot at step k.
c work arrays:
c-------------
c iw	    = integer work array of length n.
c------------
c Note (IMPORTANT):
c-----------
C it is assumed that the the elements in the input matrix are ordered
c    in such a way that in each row the lower part comes first and
c    then the upper part. To get the correct ILU factorization, it is
c    also necessary to have the elements of L ordered by increasing
c    column number. It may therefore be necessary to sort the
c    elements of a, ja, ia prior to calling milu0. This can be
c    achieved by transposing the matrix twice using csrcsc.
c-----------------------------------------------------------
          ju0 = n+2
          jlu(1) = ju0
c initialize work vector to zero's
          do 31 i=1, n
 31           iw(i) = 0
c
c-------------- MAIN LOOP ----------------------------------
c
          do 500 ii = 1, n
           js = ju0
c
c generating row number ii or L and U.
c
           do 100 j=ia(ii),ia(ii+1)-1
c
c copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
c
              jcol = ja(j)
              if (jcol .eq. ii) then
                 alu(ii) = a(j)
                 iw(jcol) = ii
                 ju(ii)  = ju0
              else
                 alu(ju0) = a(j)
                 jlu(ju0) = ja(j)
                 iw(jcol) = ju0
                 ju0 = ju0+1
              endif
 100       continue
           jlu(ii+1) = ju0
           jf = ju0-1
           jm = ju(ii)-1
c s accumulates fill-in values
           s = 0.0d0
           do 150 j=js, jm
              jrow = jlu(j)
              tl = alu(j)*alu(jrow)
              alu(j) = tl
c-----------------------perform linear combination --------
              do 140 jj = ju(jrow), jlu(jrow+1)-1
                 jw = iw(jlu(jj))
                 if (jw .ne. 0) then
                       alu(jw) = alu(jw) - tl*alu(jj)
                    else
                       s = s + tl*alu(jj)
                    endif
 140          continue
 150       continue
c----------------------- invert and store diagonal element.
           alu(ii) = alu(ii)-s
           if (alu(ii) .eq. 0.0d0) goto 600
           alu(ii) = 1.0d0/alu(ii)
c----------------------- reset pointer iw to zero
           iw(ii) = 0
           do 201 i = js, jf
 201          iw(jlu(i)) = 0
 500       continue
           ierr = 0
           return
c     zero pivot :
 600       ierr = ii
           return
c------- end of milu0 --------------------------------------------------
c-----------------------------------------------------------------------
           end
c-----------------------------------------------------------------------
      subroutine msrcsr (n,a,ja,ao,jao,iao,wk)
      real*8 a(*),ao(*),wk(n)
      integer ja(*),jao(*),iao(n+1)
c-----------------------------------------------------------------------
c       Modified - Sparse Row  to   Compressed Sparse Row
c
c-----------------------------------------------------------------------
c converts a compressed matrix using a separated diagonal
c (modified sparse row format) in the Compressed Sparse Row
c format.
c does not check for zero elements in the diagonal.
c
c
c-----  correzione G. PINI   28/9/93
c on entry :
c---------
c n        = row dimension of matrix
c a, ja  = sparse matrix in msr sparse storage format
c            see routine csrmsr for details
c
c on return :
c-----------
c
c a0, ja0, ia0 = matrix in csr format. note that the
c            algorithm is in place: ao, jao can be the same
c            as a, ja, in which case it will be overwritten on it
c            upon return.
c
c              here nnz = number of nonzero elements+1
c work arrays:
c------------
c wk    = real work array of length n
c
c notes:
c------- In place algorithm (see a, ja, ia).
c-----------------------------------------------------------------------
      logical added
      do 1 i=1,n
         wk(i) = a(i)
 1    continue
      iao(1) = 1
      iptr = 1
c---------
      do 500 ii=1,n
         added = .false.
         idiag = iptr + (ja(ii+1)-ja(ii))
         do 100 k=ja(ii),ja(ii+1)-1
            j = ja(k)
            if (j .lt. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr+1
            elseif (added) then
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr+1
            else
c add diag element - only reserve a position for it.
               idiag = iptr
               iptr = iptr+1
               added = .true.
c     then other element
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr+1
            endif
 100     continue
         ao(idiag) = wk(ii)
         jao(idiag) = ii
         if (.not. added) iptr = iptr+1
         iao(ii+1) = iptr
 500  continue
      return
c------------ end of subroutine csrmsr ---------------------------------
c-----------------------------------------------------------------------
      end
C
C
C**************************  NSYSLT ************************************
C
C  Solve a nonsymmetric linear system of equations using a method
C  of minimum residuals (nonsymmetric conjugate gradients) or
C  a nonsymmetric direct solver (note: COEF2, COEF3, and SCR
C  in the conjugate gradients routines are used as scratch vectors)
c---------------  precondizionatore
c  METH = 1   ILU(0)
c       = 2   MILU
c       = 3   ILUT
c----------------  solutore
C  ISOLV=-6   CGMRES ( precondizionto con K-1)
C       =-5   BCGSTAB ( con precondizionamento diagonale)
C       =-4   BCGSTAB ( non precondizionato)
C       =-3   TFQMR ( con precondizionamento diagonale)
C       =-2   TFQMR ( non precondizionato)
C       =-1   TFQMR (precondizionato con K-1)
C       =0    BCGSTAB (precondizionato con K-1)
C       =1    minimum residuals (GRAMRB).
C       =2    GCRK(5)
C       =3    IBM's NONSYM
C  Calculate the residual from the conjugate gradients or NONSYM
C  solution of the nonsymmetric linear system.
C
C  Note: to get the iteration number and the residuals on unit 6,
C         uncomment all the 'WRITE(6,...' statements
C
c    employes: ial,jal,coefl,dl,iau,jau,coefu
C***********************************************************************
C
      SUBROUTINE NSYSLT(ISOLV,IOUT,N,NTERM,NUMDIR,NODDIR,ITMXCG,
     1                  IBOT,MINBOT,MAXBOT,IERSYM,NITER,
     2                  IA,JA,TOPOL,INSYM,TOLCG,RMIN,COEF1,COEF2,
     3                  COEF3,SCR,RNSYM,PNEW,TNOTI,
     4                  ncoef,ial,jal,coefl,iau,jau,coefu,
     5 meth,lfil,tolilu,alu,jlu,ju,nwk,vv,iw,ac,jac,topolc,maxlu,im)
C
      IMPLICIT  NONE
      integer maxlu,ncoef,meth,ierr,im
      integer lfil,jlu(*),ju(*),nwk,iw(ncoef,*)
      integer jac(*),topolc(*)
      real*8 ac(*)
      real*8 tolilu,alu(*),vv(ncoef,*)
      INTEGER   I,J,K,NFASIN
      INTEGER   IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11
      INTEGER   IM1,IM2,IM3
      INTEGER   INS1,INS2,INS3,INS4,INS5,INS6,INS7
      INTEGER   ISOLV,IOUT,N,NTERM,NUMDIR,ITMXCG
      INTEGER   IBOT,MINBOT,MAXBOT,IERSYM,NITER
      INTEGER   NODDIR(*),IA(*),JA(*),TOPOL(*),INSYM(*)
      REAL*8    ERR,RES,RHSN,RESINI
      REAL*8    TOLCG,RMIN,RMAX
      REAL*8    COEF1(*),COEF2(*),COEF3(*),SCR(*),RNSYM(*)
      REAL*8    PNEW(*),TNOTI(*)
      integer   ial(*),jal(*),iau(*),jau(*),in12
      real*8    coefl(*),coefu(*)
      integer its
      LOGICAL   FATSYM
      COMMON /STAT/ITS
      DATA      FATSYM/.FALSE./
      RMAX=RNSYM(1)
C
CCCC      write(6,*)'      NSYSLT  (ilut  - gmres -  pnew=0)    '
C
         IN1=1
         IN2=IN1+N
         IN3=IN2+N
         IN4=IN3+N
         IN5=IN4+N
         IN6=IN5+N
         IN7=IN6+N
         IN8=IN7+N
         IN9=IN8+N
         IN10=IN9+N
         IN11=IN10+N
         in12=in11+n
         IM1=1
         IM2=IM1+5*N
         IM3=IM2+5*N
C

      IF (ISOLV .LE. 2) THEN
C
C  calculates the norm of the RHS vector
C
         RHSN=0.0D0
         DO I =1,N
            RHSN=RHSN+TNOTI(I)*TNOTI(I)
         END DO
         IF (RHSN .EQ. 0.0D0) THEN
             DO I=1,N
               PNEW(I)=0.0D0
             END DO
             RETURN
         END IF
         if(isolv.eq.-6) then
            CALL ilut1(meth,n,ncoef,coef1,ja,topol,lfil,tolilu,alu,jlu,
     1                 ju,nwk,vv,iw,ierr)
            do i=1,n
               scr(i)=tnoti(i)
               scr(in2+i-1)=tnoti(i)
            end do
c  x0=K**(-1)b
            call lusol0(n,scr(in2),pnew,alu,jlu,ju)

c-----------------------
            CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN3),N,
     1                  NUMDIR,NODDIR,RESINI,RMIN)

CCCC          write(6,143) res
            call pgmres(n,im,scr,pnew,vv,tolcg,itmxcg,99,coef1,
     1                 ja,topol,alu,jlu,ju,ierr,err)
C
C  calculate the residual
C
            CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1                  NUMDIR,NODDIR,RES,RMIN)
            WRITE(IOUT,1000) its,RESINI,ERR,RES
            NITER=its
            return
         endif
         if(isolv.ge.-1) then
            CALL ilutme(meth,n,ncoef,coef1,ja,topol,lfil,tolilu,alu,jlu,
     1                 ju,nwk,vv,iw,ierr,ac,jac,topolc,
     2                 ial,jal,coefl,iau,jau,coefu,scr(in12))
CCC            write(iout,*)'n. el L=', ial(n+1)-1,'   n. el. U=',iau(n)
            if(iau(n).gt.maxlu) then
               write(iout,*)' aumentare il valore maxlu'
               stop
            endif
C
c  x0=K^-1b   di solito
              call yprdly(scr,tnoti,n,ial,jal,coefl,scr(in12))
              call yprduy(pnew,scr,n,iau,jau,coefu)
         endif
         if((isolv.eq.-5).or.(isolv.eq.-3)) then
            call diagn(n,topol,ia,ja,coef1,scr(in11))
            call prdy(pnew,tnoti,n,scr(in11))
         endif
         if( (isolv.eq.-4).or.(isolv.eq.-2)) then
            call setone(numdir,noddir,topol,ja,coef1,tnoti,rmax)
            do i=1,n
               pnew(I)=tnoti(i)
            end do
         endif
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1               NUMDIR,NODDIR,RESINI,RMIN)
CCC          write(6,143) res
CCC 143      format(' initial residual=',e12.5)
         IF(ISOLV.EQ.-5) THEN
            CALL GCSTAD(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN9),
     3                  NUMDIR,NODDIR,NITER,ERR,SCR(IN11))
         ELSE IF(ISOLV.EQ.-4) THEN
            CALL GCSTAN(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN9),
     3                  NUMDIR,NODDIR,NITER,ERR)
         ELSE IF(ISOLV.EQ.-3) THEN
            CALL TFQWD(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                 SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                 SCR(IN9),TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN10),
     3                 NUMDIR,NODDIR,NITER,ERR,SCR(IN11))
         ELSE IF(ISOLV.EQ.-2) THEN
            CALL TFQWN(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                 SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                 SCR(IN9),TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN10),
     3                 NUMDIR,NODDIR,NITER,ERR)
         ELSE IF(ISOLV.EQ.-1) THEN
            CALL TFQWut(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                SCR(IN9),TOPOL,JA,COEF1,NTERM,ITMXCG,
     3                TOLCG,SCR(IN10),SCR(IN11),NUMDIR,NODDIR,NITER,
     4                ERR,
     5                ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF(ISOLV.EQ.0) THEN
            CALL GCSTUT(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,NTERM,
     3                  ITMXCG,TOLCG,SCR(IN9),SCR(IN10),
     4                  NUMDIR,NODDIR,NITER,ERR,
     5                  ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF (ISOLV .EQ. 1) THEN
            CALL GRAMUT(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),JA,COEF1,SCR(IN5),NTERM,TOPOL,
     2                  ITMXCG,TOLCG,SCR(IN6),SCR(IN7),NUMDIR,NODDIR,
     3                  NITER,ERR,
     4                  ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF(ISOLV.EQ.2) THEN
            CALL GCRKUT(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                SCR(IN4),JA,COEF1,SCR(IN5),NTERM,TOPOL,
     2                ITMXCG,TOLCG,SCR(IN6),SCR(IN7),5,
     3                COEF3(IM1),COEF3(IM2),COEF3(IM3),NUMDIR,NODDIR,
     4                NITER,ERR,
     5                ial,jal,coefl,iau,jau,coefu,scr(in12))
         ENDIF
C
C  calculate the residual
C
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1               NUMDIR,NODDIR,RES,RMIN)
         WRITE(IOUT,1020) NITER,RESINI,ERR,RES,IAL(N+1)-1,IAU(N)
      ELSE
         INS1=1
         INS2=INS1+N
         INS3=INS2+N
         INS4=INS3+N+1
         INS5=INS4+N
         INS6=INS5+N
         INS7=INS6+N
C
C  in using the direct solver (NONSYM) to solve the system, we
C  need to perform symbolic factorization, and we do so only once
C
         IF (.NOT. FATSYM) THEN
            FATSYM=.TRUE.
            DO I=1,N
               J=N+I
               INSYM(I)=I
               INSYM(J)=I
            END DO
            IBOT=MAXBOT
            MINBOT=IBOT
            NFASIN=N
            CALL FASINS(NFASIN,IBOT,INSYM(INS1),INSYM(INS2),TOPOL,JA,
     1                  INSYM(INS3),INSYM(INS4),INSYM(INS7),K,IERSYM,
     2                  INSYM(INS5),INSYM(INS6))
            IF (IERSYM .NE. 0) THEN
               WRITE(IOUT,1300) IERSYM,NFASIN,N
               RETURN
            END IF
            MINBOT=INSYM(INS3+N)+N
         END IF
         CALL FANUNS(N,INSYM(INS1),INSYM(INS2),TOPOL,JA,COEF1,
     1               INSYM(INS3),INSYM(INS4),INSYM(INS7),SCR(INS1),
     2               RNSYM,IERSYM,SCR(INS2))
         IF (IERSYM .NE. 0) THEN
            WRITE(IOUT,1310) IERSYM
            RETURN
         END IF
         CALL AVINNS(N,INSYM(INS1),INSYM(INS2),INSYM(INS3),INSYM(INS4),
     1               INSYM(INS7),SCR(INS1),RNSYM,TNOTI,PNEW)
C
C  calculate the residual
C
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(INS1),N,
     1               NUMDIR,NODDIR,RES,RMIN)
         WRITE(IOUT,1010) RESINI,RES
      END IF
C
      RETURN
 1000 FORMAT(1X,I4,3(1PE15.6),'  <<NONSYMMETRIC SOLVER>>')
 1010 FORMAT(5X,1PE15.6,15X,1PE15.6,'  <<NONSYMMETRIC SOLVER>>')
 1020 FORMAT(1X,I4,3(1PE15.6),2I10,
     1       '  (RESINI,ERR,RES,N.EL. L,N.EL. U)')
 1300 FORMAT(//,1X,' ERROR FROM NONSYM ROUTINE FASINS, CODE ',I4,
     1        /,6X,I7,' NODES PROCESSED SO FAR OUT OF ',I7)
 1310 FORMAT(//,1X,' ERROR FROM NONSYM ROUTINE FANUNS, CODE ',I4)
      END
C
C**************************  NSYSLU ************************************
C  (funziona solo con matrici aventi pattern simmetrico)
c-------------------------------------------------------------
C  Solve a nonsymmetric linear system of equations using a method
C  of minimum residuals (nonsymmetric conjugate gradients) or
C  a nonsymmetric direct solver (note: COEF2, COEF3, and SCR
C  in the conjugate gradients routines are used as scratch vectors)
C  ISOLV=-5   BCGSTAB ( con precondizionamento diagonale)
C       =-4   BCGSTAB ( non precondizionato)
C       =-3   TFQMR ( con precondizionamento diagonale)
C       =-2   TFQMR ( non precondizionato)
C       =-1   TFQMR (precondizionato con K-1)
C       =0    BCGSTAB (precondizionato con K-1)
C       =1    minimum residuals (GRAMRB).
C       =2    GCRK(5)
C       =3    IBM's NONSYM
C  Calculate the residual from the conjugate gradients or NONSYM
C  solution of the nonsymmetric linear system.
C
C  Note: to get the iteration number and the residuals on unit 6,
C         uncomment all the 'WRITE(6,...' statements
C
c    employes: ial,jal,coefl,dl,iau,jau,coefu
C***********************************************************************
C
      SUBROUTINE NSYSLU(ISOLV,IOUT,N,NTERM,NUMDIR,NODDIR,ITMXCG,
     1                  IBOT,MINBOT,MAXBOT,IERSYM,NITER,
     2                  IA,JA,TOPOL,INSYM,TOLCG,RMIN,COEF1,COEF2,
     3                  COEF3,SCR,RNSYM,PNEW,TNOTI,
     4                  ial,jal,coefl,iau,jau,coefu)
C
      IMPLICIT  NONE
      INTEGER   I,J,K,NFASIN
      INTEGER   IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11
      INTEGER   IM1,IM2,IM3
      INTEGER   INS1,INS2,INS3,INS4,INS5,INS6,INS7
      INTEGER   ISOLV,IOUT,N,NTERM,NUMDIR,ITMXCG
      INTEGER   IBOT,MINBOT,MAXBOT,IERSYM,NITER
      INTEGER   NODDIR(*),IA(*),JA(*),TOPOL(*),INSYM(*)
      REAL*8    ERR,RES,RHSN,RESINI
      REAL*8    TOLCG,RMIN,RMAX
      REAL*8    COEF1(*),COEF2(*),COEF3(*),SCR(*),RNSYM(*)
      REAL*8    PNEW(*),TNOTI(*)
      integer   ial(*),jal(*),iau(*),jau(*),in12
      real*8    coefl(*),coefu(*)
C
      LOGICAL   FATSYM
      DATA      FATSYM/.FALSE./
      RMAX=RNSYM(1)
C
CCC      write(6,*)'      NSYSLU  (inclu  - ial,jal,coefl,..) '
C
      IF (ISOLV .LE. 2) THEN
C
C  calculates the norm of the RHS vector
C
         RHSN=0.0D0
         DO I =1,N
            RHSN=RHSN+TNOTI(I)*TNOTI(I)
         END DO
         IF (RHSN .EQ. 0.0D0) THEN
             DO I=1,N
               PNEW(I)=0.0D0
             END DO
             RETURN
         END IF
         IN1=1
         IN2=IN1+N
         IN3=IN2+N
         IN4=IN3+N
         IN5=IN4+N
         IN6=IN5+N
         IN7=IN6+N
         IN8=IN7+N
         IN9=IN8+N
         IN10=IN9+N
         IN11=IN10+N
         in12=in11+n
         IM1=1
         IM2=IM1+5*N
         IM3=IM2+5*N
         if(isolv.ge.-1) then
            CALL INCLU(IOUT,IA,JA,TOPOL,COEF1,COEF3,COEF2,SCR,N,NTERM)
            call prepy(n,ia,ja,topol,ial,jal,iau,jau,coef2,
     1                 coefl,coefu,scr(in12))
            call yprdly(scr,tnoti,n,ial,jal,coefl,scr(in12))
            call yprduy(pnew,scr,n,iau,jau,coefu)
         endif
         if((isolv.eq.-5).or.(isolv.eq.-3)) then
            call diagn(n,topol,ia,ja,coef1,scr(in11))
            call prdy(pnew,tnoti,n,scr(in11))
         endif
         if( (isolv.eq.-4).or.(isolv.eq.-2)) then
            call setone(numdir,noddir,topol,ja,coef1,tnoti,rmax)
            do i=1,n
               pnew(I)=tnoti(i)
            end do
         endif
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1               NUMDIR,NODDIR,RESINI,RMIN)
          write(99,*) 0,resini
         IF(ISOLV.EQ.-5) THEN
            CALL GCSTAD(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN9),
     3                  NUMDIR,NODDIR,NITER,ERR,SCR(IN11))
         ELSE IF(ISOLV.EQ.-4) THEN
            CALL GCSTAN(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN9),
     3                  NUMDIR,NODDIR,NITER,ERR)
         ELSE IF(ISOLV.EQ.-3) THEN
            CALL TFQWD(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                 SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                 SCR(IN9),TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN10),
     3                 NUMDIR,NODDIR,NITER,ERR,SCR(IN11))
         ELSE IF(ISOLV.EQ.-2) THEN
            CALL TFQWN(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                 SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                 SCR(IN9),TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN10),
     3                 NUMDIR,NODDIR,NITER,ERR)
         ELSE IF(ISOLV.EQ.-1) THEN
            CALL TFQWut(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                SCR(IN9),TOPOL,JA,COEF1,NTERM,ITMXCG,
     3                TOLCG,SCR(IN10),SCR(IN11),NUMDIR,NODDIR,NITER,
     4                ERR,
     5                ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF(ISOLV.EQ.0) THEN
            CALL GCSTut(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,NTERM,
     3                  ITMXCG,TOLCG,SCR(IN9),SCR(IN10),
     4                  NUMDIR,NODDIR,NITER,ERR,
     5                  ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF (ISOLV .EQ. 1) THEN
            CALL GRAMut(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),JA,COEF1,SCR(IN5),NTERM,TOPOL,
     2                  ITMXCG,TOLCG,SCR(IN6),SCR(IN7),NUMDIR,NODDIR,
     3                  NITER,ERR,
     4                  ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF(ISOLV.EQ.2) THEN
            CALL GCRKut(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                SCR(IN4),JA,COEF1,SCR(IN5),NTERM,TOPOL,
     2                ITMXCG,TOLCG,SCR(IN6),SCR(IN7),5,
     3                COEF3(IM1),COEF3(IM2),COEF3(IM3),NUMDIR,NODDIR,
     4                NITER,ERR,
     5                ial,jal,coefl,iau,jau,coefu,scr(in12))
         ENDIF
C
C  calculate the residual
C
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1               NUMDIR,NODDIR,RES,RMIN)
         WRITE(IOUT,1000) NITER,RESINI,ERR,RES
      ELSE
         INS1=1
         INS2=INS1+N
         INS3=INS2+N
         INS4=INS3+N+1
         INS5=INS4+N
         INS6=INS5+N
         INS7=INS6+N
C
C  in using the direct solver (NONSYM) to solve the system, we
C  need to perform symbolic factorization, and we do so only once
C
         IF (.NOT. FATSYM) THEN
            FATSYM=.TRUE.
            DO I=1,N
               J=N+I
               INSYM(I)=I
               INSYM(J)=I
            END DO
            IBOT=MAXBOT
            MINBOT=IBOT
            NFASIN=N
            CALL FASINS(NFASIN,IBOT,INSYM(INS1),INSYM(INS2),TOPOL,JA,
     1                  INSYM(INS3),INSYM(INS4),INSYM(INS7),K,IERSYM,
     2                  INSYM(INS5),INSYM(INS6))
            IF (IERSYM .NE. 0) THEN
               WRITE(IOUT,1300) IERSYM,NFASIN,N
               RETURN
            END IF
            MINBOT=INSYM(INS3+N)+N
         END IF
         CALL FANUNS(N,INSYM(INS1),INSYM(INS2),TOPOL,JA,COEF1,
     1               INSYM(INS3),INSYM(INS4),INSYM(INS7),SCR(INS1),
     2               RNSYM,IERSYM,SCR(INS2))
         IF (IERSYM .NE. 0) THEN
            WRITE(IOUT,1310) IERSYM
            RETURN
         END IF
         CALL AVINNS(N,INSYM(INS1),INSYM(INS2),INSYM(INS3),INSYM(INS4),
     1               INSYM(INS7),SCR(INS1),RNSYM,TNOTI,PNEW)
C
C  calculate the residual
C
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(INS1),N,
     1               NUMDIR,NODDIR,RES,RMIN)
         WRITE(IOUT,1010) RES
      END IF
C
      RETURN
 1000 FORMAT(1X,I4,3(1PE15.6),'  <<NONSYMMETRIC SOLVER>>')
 1010 FORMAT(20X,1PE15.6,'  <<NONSYMMETRIC SOLVER>>')
 1300 FORMAT(//,1X,' ERROR FROM NONSYM ROUTINE FASINS, CODE ',I4,
     1        /,6X,I7,' NODES PROCESSED SO FAR OUT OF ',I7)
 1310 FORMAT(//,1X,' ERROR FROM NONSYM ROUTINE FANUNS, CODE ',I4)
      END
C
C**************************  NSYSLV ************************************
C
C  Solve a nonsymmetric linear system of equations using a method
C  of minimum residuals (nonsymmetric conjugate gradients) or
C  a nonsymmetric direct solver (note: COEF2, COEF3, and SCR
C  in the conjugate gradients routines are used as scratch vectors)
C  ISOLV=-5   BCGSTAB ( con precondizionamento diagonale)
C       =-4   BCGSTAB ( non precondizionato)
C       =-3   TFQMR ( con precondizionamento diagonale)
C       =-2   TFQMR ( non precondizionato)
C       =-1   TFQMR (precondizionato con K-1)
C       =0    BCGSTAB (precondizionato con K-1)
C       =1    minimum residuals (GRAMRB).
C       =2    GCRK(5)
C       =3    IBM's NONSYM
C  Calculate the residual from the conjugate gradients or NONSYM
C  solution of the nonsymmetric linear system.
C
C  Note: to get the iteration number and the residuals on unit 6,
C         uncomment all the 'WRITE(6,...' statements
C
C***********************************************************************
C
      SUBROUTINE NSYSLV(ISOLV,IOUT,N,NTERM,NUMDIR,NODDIR,ITMXCG,
     1                  IBOT,MINBOT,MAXBOT,IERSYM,NITER,
     2                  IA,JA,TOPOL,INSYM,TOLCG,RMIN,COEF1,COEF2,
     3                  COEF3,SCR,RNSYM,PNEW,TNOTI)
C
      IMPLICIT  NONE
      INTEGER   I,J,K,NFASIN
      INTEGER   IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11
      INTEGER   IM1,IM2,IM3
      INTEGER   INS1,INS2,INS3,INS4,INS5,INS6,INS7
      INTEGER   ISOLV,IOUT,N,NTERM,NUMDIR,ITMXCG
      INTEGER   IBOT,MINBOT,MAXBOT,IERSYM,NITER
      INTEGER   NODDIR(*),IA(*),JA(*),TOPOL(*),INSYM(*)
      REAL*8    ERR,RES,RHSN,RESINI
      REAL*8    TOLCG,RMIN,RMAX
      REAL*8    COEF1(*),COEF2(*),COEF3(*),SCR(*),RNSYM(*)
      REAL*8    PNEW(*),TNOTI(*)
      LOGICAL   FATSYM
      DATA      FATSYM/.FALSE./
      RMAX=RNSYM(1)
C
CCCC      write(*,*)'  NSYSLV   (inclu  - ia,ja,...)  '
C
      IF (ISOLV .LE. 2) THEN
C
C  calculates the norm of the RHS vector
C
         RHSN=0.0D0
         DO I =1,N
            RHSN=RHSN+TNOTI(I)*TNOTI(I)
         END DO
         IF (RHSN .EQ. 0.0D0) THEN
             DO I=1,N
               PNEW(I)=0.0D0
             END DO
             RETURN
         END IF
         IN1=1
         IN2=IN1+N
         IN3=IN2+N
         IN4=IN3+N
         IN5=IN4+N
         IN6=IN5+N
         IN7=IN6+N
         IN8=IN7+N
         IN9=IN8+N
         IN10=IN9+N
         IN11=IN10+N
         IM1=1
         IM2=IM1+5*N
         IM3=IM2+5*N
         if(isolv.ge.-1) then
            CALL INCLU(IOUT,IA,JA,TOPOL,COEF1,COEF3,COEF2,SCR,N,NTERM)
            CALL PROD1(PNEW,TNOTI,N,IA,JA,COEF2,NTERM,TOPOL)
         endif
         if((isolv.eq.-5).or.(isolv.eq.-3)) then
            call diagn(n,topol,ia,ja,coef1,scr(in11))
            call prdy(pnew,tnoti,n,scr(in11))
         endif
         if( (isolv.eq.-4).or.(isolv.eq.-2)) then
            call setone(numdir,noddir,topol,ja,coef1,tnoti,rmax)
            do i=1,n
               pnew(I)=tnoti(i)
            end do
         endif
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1               NUMDIR,NODDIR,RESINI,RMIN)
          write(99,*) 0,resini
CCCC          write(6,*)' residuo iniziale',res
         IF(ISOLV.EQ.-5) THEN
            CALL GCSTAD(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN9),
     3                  NUMDIR,NODDIR,NITER,ERR,SCR(IN11))
         ELSE IF(ISOLV.EQ.-4) THEN
            CALL GCSTAN(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN9),
     3                  NUMDIR,NODDIR,NITER,ERR)
         ELSE IF(ISOLV.EQ.-3) THEN
            CALL TFQWD(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                 SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                 SCR(IN9),TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN10),
     3                 NUMDIR,NODDIR,NITER,ERR,SCR(IN11))
         ELSE IF(ISOLV.EQ.-2) THEN
            CALL TFQWN(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                 SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                 SCR(IN9),TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN10),
     3                 NUMDIR,NODDIR,NITER,ERR)
         ELSE IF(ISOLV.EQ.-1) THEN
            CALL TFQW(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                SCR(IN9),TOPOL,IA,JA,COEF1,COEF2,NTERM,ITMXCG,
     3                TOLCG,SCR(IN10),SCR(IN11),NUMDIR,NODDIR,NITER,
     4                ERR)
         ELSE IF(ISOLV.EQ.0) THEN
             CALL GCSTAS(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,IA,JA,COEF1,COEF2,NTERM,
     3                  ITMXCG,TOLCG,SCR(IN9),SCR(IN10),
     4                  NUMDIR,NODDIR,NITER,ERR)
         ELSE IF (ISOLV .EQ. 1) THEN
            CALL GRAMRB(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),IA,JA,COEF1,COEF2,SCR(IN5),NTERM,TOPOL,
     2                  ITMXCG,TOLCG,SCR(IN6),SCR(IN7),NUMDIR,NODDIR,
     3                  NITER,ERR)
         ELSE IF(ISOLV.EQ.2) THEN
            CALL GCRK(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                SCR(IN4),IA,JA,COEF1,COEF2,SCR(IN5),NTERM,TOPOL,
     2                ITMXCG,TOLCG,SCR(IN6),SCR(IN7),5,
     3                COEF3(IM1),COEF3(IM2),COEF3(IM3),NUMDIR,NODDIR,
     4                NITER,ERR)
         ENDIF
C
C  calculate the residual
C
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1               NUMDIR,NODDIR,RES,RMIN)
         WRITE(IOUT,1000) NITER,RESINI,ERR,RES
      ELSE
         INS1=1
         INS2=INS1+N
         INS3=INS2+N
         INS4=INS3+N+1
         INS5=INS4+N
         INS6=INS5+N
         INS7=INS6+N
C
C  in using the direct solver (NONSYM) to solve the system, we
C  need to perform symbolic factorization, and we do so only once
C
         IF (.NOT. FATSYM) THEN
            FATSYM=.TRUE.
            DO I=1,N
               J=N+I
               INSYM(I)=I
               INSYM(J)=I
            END DO
            IBOT=MAXBOT
            MINBOT=IBOT
            NFASIN=N
            CALL FASINS(NFASIN,IBOT,INSYM(INS1),INSYM(INS2),TOPOL,JA,
     1                  INSYM(INS3),INSYM(INS4),INSYM(INS7),K,IERSYM,
     2                  INSYM(INS5),INSYM(INS6))
            IF (IERSYM .NE. 0) THEN
               WRITE(IOUT,1300) IERSYM,NFASIN,N
               RETURN
            END IF
            MINBOT=INSYM(INS3+N)+N
         END IF
         CALL FANUNS(N,INSYM(INS1),INSYM(INS2),TOPOL,JA,COEF1,
     1               INSYM(INS3),INSYM(INS4),INSYM(INS7),SCR(INS1),
     2               RNSYM,IERSYM,SCR(INS2))
         IF (IERSYM .NE. 0) THEN
            WRITE(IOUT,1310) IERSYM
            RETURN
         END IF
         CALL AVINNS(N,INSYM(INS1),INSYM(INS2),INSYM(INS3),INSYM(INS4),
     1               INSYM(INS7),SCR(INS1),RNSYM,TNOTI,PNEW)
C
C  calculate the residual
C
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(INS1),N,
     1               NUMDIR,NODDIR,RES,RMIN)
         WRITE(IOUT,1010) RESINI,RES
      END IF
C
      RETURN
 1000 FORMAT(1X,I4,3(1PE15.6),'  <<NONSYMMETRIC SOLVER>>')
 1010 FORMAT(5X,1PE15.6,15X,1PE15.6,'  <<NONSYMMETRIC SOLVER>>')
 1300 FORMAT(//,1X,' ERROR FROM NONSYM ROUTINE FASINS, CODE ',I4,
     1        /,6X,I7,' NODES PROCESSED SO FAR OUT OF ',I7)
 1310 FORMAT(//,1X,' ERROR FROM NONSYM ROUTINE FANUNS, CODE ',I4)
      END
C
C
C**************************  NSYST ************************************
C
C  Solve a nonsymmetric linear system of equations using a method
C  of minimum residuals (nonsymmetric conjugate gradients) or
C  a nonsymmetric direct solver (note: COEF2, COEF3, and SCR
C  in the conjugate gradients routines are used as scratch vectors)
c---------------  precondizionatore
c  METH = 1   ILU(0)
c       = 2   MILU
c       = 3   ILUT
c----------------  solutore
C  ISOLV=-6   CGMRES ( precondizionto con K-1)
C       =-5   BCGSTAB ( con precondizionamento diagonale)
C       =-4   BCGSTAB ( non precondizionato)
C       =-3   TFQMR ( con precondizionamento diagonale)
C       =-2   TFQMR ( non precondizionato)
C       =-1   TFQMR (precondizionato con K-1)
C       =0    BCGSTAB (precondizionato con K-1)
C       =1    minimum residuals (GRAMRB).
C       =2    GCRK(5)
C       =3    IBM's NONSYM
C  Calculate the residual from the conjugate gradients or NONSYM
C  solution of the nonsymmetric linear system.
C
C  Note: to get the iteration number and the residuals on unit 6,
C         uncomment all the 'WRITE(6,...' statements
C
c    employes: ial,jal,coefl,dl,iau,jau,coefu
C***********************************************************************
C
      SUBROUTINE NSYST(ISOLV,IOUT,N,NTERM,NUMDIR,NODDIR,ITMXCG,
     1                  IBOT,MINBOT,MAXBOT,IERSYM,NITER,
     2                  IA,JA,TOPOL,INSYM,TOLCG,RMIN,COEF1,COEF2,
     3                  COEF3,SCR,RNSYM,PNEW,TNOTI,
     4                  ncoef,ial,jal,coefl,iau,jau,coefu,
     5                  meth,lfil,tolilu,nwk,maxlu,im,maxdir)
C
      IMPLICIT  NONE
      integer lfil,nwk,maxdir
      integer maxlu,ncoef,meth,ierr,im
      real*8 tolilu
      INTEGER   I,J,K,NFASIN
      INTEGER   IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11
      INTEGER   IM1,IM2,IM3
      INTEGER   INS1,INS2,INS3,INS4,INS5,INS6,INS7
      integer it1,it2,it3,is1,is2,is3,is4,is5
      INTEGER   ISOLV,IOUT,N,NTERM,NUMDIR,ITMXCG
      INTEGER   IBOT,MINBOT,MAXBOT,IERSYM,NITER
      INTEGER   NODDIR(*),IA(*),JA(*),TOPOL(*),INSYM(*)
      REAL*8    ERR,RES,RHSN,RESINI
      REAL*8    TOLCG,RMIN,RMAX
      REAL*8    COEF1(*),COEF2(*),COEF3(*),SCR(*),RNSYM(*)
      REAL*8    PNEW(*),TNOTI(*)
      integer   ial(*),jal(*),iau(*),jau(*),in12
      real*8    coefl(*),coefu(*)
      integer its
      LOGICAL   FATSYM
      COMMON /STAT/ITS
      DATA      FATSYM/.FALSE./
      RMAX=RNSYM(1)
C
CCCC      write(6,*)'      NSYSLT  (ilut  - gmres -  pnew=K^-1 b )    '
C
         IN1=1
         IN2=IN1+N
         IN3=IN2+N
         IN4=IN3+N
         IN5=IN4+N
         IN6=IN5+N
         IN7=IN6+N
         IN8=IN7+N
         IN9=IN8+N
         IN10=IN9+N
         IN11=IN10+N
         in12=in11+n
         IM1=1
         IM2=IM1+5*N
         IM3=IM2+5*N
C
         it1=in12+n
         it2=it1+2*maxlu
         it3=it2+n*maxdir
         is1=1
         is2=is1+2*maxlu
         is3=is2+n
         is4=is3+3*n
         is5=is4+2*maxlu
cc    alu=scr(it1), vv=scr(it2), al=scr(it3)
cc    jlu=insym(is1), ju=insym(is2), iw=insym(is3),
cc    jac=insym(is4), topolc=insym(is5)
c
      IF (ISOLV .LE. 2) THEN
C
C  calculates the norm of the RHS vector
C
         RHSN=0.0D0
         DO I =1,N
            RHSN=RHSN+TNOTI(I)*TNOTI(I)
         END DO
         IF (RHSN .EQ. 0.0D0) THEN
             DO I=1,N
               PNEW(I)=0.0D0
             END DO
             RETURN
         END IF
         if(isolv.eq.-6) then
          CALL ilut1(meth,n,ncoef,coef1,ja,topol,lfil,tolilu,scr(it1)
     1     ,insym(is1),insym(is2),nwk,scr(it2),insym(is3),ierr)
            do i=1,n
               scr(i)=tnoti(i)
               scr(in2+i-1)=tnoti(i)
            end do
c  x0=K**(-1)b
           call lusol0(n,scr(in2),pnew,scr(it1),insym(is1),insym(is2))


c-----------------------
            CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN3),N,
     1                  NUMDIR,NODDIR,RESINI,RMIN)
          write(99,*) 0,resini
CCCC          write(6,143) res
            call pgmres(n,im,scr,pnew,scr(it2),tolcg,itmxcg,99,coef1,
     1            ja,topol,scr(it1),insym(is1),insym(is2),ierr,err)
C
C  calculate the residual
C
            CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1                  NUMDIR,NODDIR,RES,RMIN)
            WRITE(IOUT,1000) its,RESINI,ERR,RES
            return
         endif
         if(isolv.ge.-1) then
            CALL ilutme(meth,n,ncoef,coef1,ja,topol,lfil,tolilu,
     1               scr(it1),insym(is1),insym(is2),nwk,scr(it2),
     2               insym(is3),ierr,scr(it3),insym(is4),insym(is5),
     3               ial,jal,coefl,iau,jau,coefu,scr(in12))
CCC            write(iout,*)'n. el L=', ial(n+1)-1,'   n. el. U=',iau(n)
            if(iau(n).gt.maxlu) then
               write(iout,*)' aumentare il valore maxlu'
               stop
            endif
C
c  x0=K^-1b   di solito
              call yprdly(scr,tnoti,n,ial,jal,coefl,scr(in12))
              call yprduy(pnew,scr,n,iau,jau,coefu)
         endif
         if((isolv.eq.-5).or.(isolv.eq.-3)) then
            call diagn(n,topol,ia,ja,coef1,scr(in11))
            call prdy(pnew,tnoti,n,scr(in11))
         endif
         if( (isolv.eq.-4).or.(isolv.eq.-2)) then
            call setone(numdir,noddir,topol,ja,coef1,tnoti,rmax)
            do i=1,n
               pnew(I)=tnoti(i)
            end do
         endif
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1               NUMDIR,NODDIR,RESINI,RMIN)
CCC          write(6,143) res
CCC 143      format(' initial residual=',e12.5)
         IF(ISOLV.EQ.-5) THEN
            CALL GCSTAD(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN9),
     3                  NUMDIR,NODDIR,NITER,ERR,SCR(IN11))
         ELSE IF(ISOLV.EQ.-4) THEN
            CALL GCSTAN(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN9),
     3                  NUMDIR,NODDIR,NITER,ERR)
         ELSE IF(ISOLV.EQ.-3) THEN
            CALL TFQWD(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                 SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                 SCR(IN9),TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN10),
     3                 NUMDIR,NODDIR,NITER,ERR,SCR(IN11))
         ELSE IF(ISOLV.EQ.-2) THEN
            CALL TFQWN(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                 SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                 SCR(IN9),TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN10),
     3                 NUMDIR,NODDIR,NITER,ERR)
         ELSE IF(ISOLV.EQ.-1) THEN
            CALL TFQWut(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                SCR(IN9),TOPOL,JA,COEF1,NTERM,ITMXCG,
     3                TOLCG,SCR(IN10),SCR(IN11),NUMDIR,NODDIR,NITER,
     4                ERR,
     5                ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF(ISOLV.EQ.0) THEN
            CALL GCSTUT(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,NTERM,
     3                  ITMXCG,TOLCG,SCR(IN9),SCR(IN10),
     4                  NUMDIR,NODDIR,NITER,ERR,
     5                  ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF (ISOLV .EQ. 1) THEN
            CALL GRAMUT(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),JA,COEF1,SCR(IN5),NTERM,TOPOL,
     2                  ITMXCG,TOLCG,SCR(IN6),SCR(IN7),NUMDIR,NODDIR,
     3                  NITER,ERR,
     4                  ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF(ISOLV.EQ.2) THEN
            CALL GCRKUT(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                SCR(IN4),JA,COEF1,SCR(IN5),NTERM,TOPOL,
     2                ITMXCG,TOLCG,SCR(IN6),SCR(IN7),5,
     3                COEF3(IM1),COEF3(IM2),COEF3(IM3),NUMDIR,NODDIR,
     4                NITER,ERR,
     5                ial,jal,coefl,iau,jau,coefu,scr(in12))
         ENDIF
C
C  calculate the residual
C
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1               NUMDIR,NODDIR,RES,RMIN)
         WRITE(IOUT,1020) NITER,RESINI,ERR,RES,IAL(N+1)-1,IAU(N)
      ELSE
         INS1=1
         INS2=INS1+N
         INS3=INS2+N
         INS4=INS3+N+1
         INS5=INS4+N
         INS6=INS5+N
         INS7=INS6+N
C
C  in using the direct solver (NONSYM) to solve the system, we
C  need to perform symbolic factorization, and we do so only once
C
         IF (.NOT. FATSYM) THEN
            FATSYM=.TRUE.
            DO I=1,N
               J=N+I
               INSYM(I)=I
               INSYM(J)=I
            END DO
            IBOT=MAXBOT
            MINBOT=IBOT
            NFASIN=N
            CALL FASINS(NFASIN,IBOT,INSYM(INS1),INSYM(INS2),TOPOL,JA,
     1                  INSYM(INS3),INSYM(INS4),INSYM(INS7),K,IERSYM,
     2                  INSYM(INS5),INSYM(INS6))
            IF (IERSYM .NE. 0) THEN
               WRITE(IOUT,1300) IERSYM,NFASIN,N
               RETURN
            END IF
            MINBOT=INSYM(INS3+N)+N
         END IF
         CALL FANUNS(N,INSYM(INS1),INSYM(INS2),TOPOL,JA,COEF1,
     1               INSYM(INS3),INSYM(INS4),INSYM(INS7),SCR(INS1),
     2               RNSYM,IERSYM,SCR(INS2))
         IF (IERSYM .NE. 0) THEN
            WRITE(IOUT,1310) IERSYM
            RETURN
         END IF
         CALL AVINNS(N,INSYM(INS1),INSYM(INS2),INSYM(INS3),INSYM(INS4),
     1               INSYM(INS7),SCR(INS1),RNSYM,TNOTI,PNEW)
C
C  calculate the residual
C
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(INS1),N,
     1               NUMDIR,NODDIR,RES,RMIN)
         WRITE(IOUT,1010) RESINI,RES
      END IF
C
      RETURN
 1000 FORMAT(1X,I4,3(1PE15.6),'  <<NONSYMMETRIC SOLVER>>')
 1010 FORMAT(5X,1PE15.6,15X,1PE15.6,'  <<NONSYMMETRIC SOLVER>>')
 1020 FORMAT(1X,I4,3(1PE15.6),2I10,
     1       '  (RESINI,ERR,RES,N.EL. L,N.EL. U)')
 1300 FORMAT(//,1X,' ERROR FROM NONSYM ROUTINE FASINS, CODE ',I4,
     1        /,6X,I7,' NODES PROCESSED SO FAR OUT OF ',I7)
 1310 FORMAT(//,1X,' ERROR FROM NONSYM ROUTINE FANUNS, CODE ',I4)
      END
C
C**************************  NSYT ************************************
C
C  Solve a nonsymmetric linear system of equations using a method
C  of minimum residuals (nonsymmetric conjugate gradients) or
C  a nonsymmetric direct solver (note: COEF2, COEF3, and SCR
C  in the conjugate gradients routines are used as scratch vectors)
c---------------  precondizionatore
c  METH = 1   ILU(0)
c       = 2   MILU
c       = 3   ILUT
c----------------  solutore
C  ISOLV=-5   BCGSTAB ( con precondizionamento diagonale)
C       =-4   BCGSTAB ( non precondizionato)
C       =-3   TFQMR ( con precondizionamento diagonale)
C       =-2   TFQMR ( non precondizionato)
C       =-1   TFQMR (precondizionato con K-1)
C       =0    BCGSTAB (precondizionato con K-1)
C       =1    minimum residuals (GRAMRB).
C       =2    GCRK(5)
C       =3    IBM's NONSYM
C  Calculate the residual from the conjugate gradients or NONSYM
C  solution of the nonsymmetric linear system.
C
C  Note: to get the iteration number and the residuals on unit 6,
C         uncomment all the 'WRITE(6,...' statements
C
c    employes: ial,jal,coefl,dl,iau,jau,coefu
C***********************************************************************
C
      SUBROUTINE NSYT(ISOLV,IOUT,N,NTERM,NUMDIR,NODDIR,ITMXCG,
     1                  IBOT,MINBOT,MAXBOT,IERSYM,NITER,
     2                  IA,JA,TOPOL,INSYM,TOLCG,RMIN,COEF1,COEF2,
     3                  COEF3,SCR,RNSYM,PNEW,TNOTI,
     4                  ncoef,ial,jal,coefl,iau,jau,coefu,
     5 meth,lfil,tolilu,alu,jlu,ju,nwk,vv,iw,ac,jac,topolc,maxlu)
C
      IMPLICIT  NONE
      integer maxlu,ncoef,meth,ierr
      integer lfil,jlu(*),ju(*),nwk,iw(ncoef,*)
      integer jac(*),topolc(*)
      real*8 ac(*)
      real*8 tolilu,alu(*),vv(ncoef,*)
      INTEGER   I,J,K,NFASIN
      INTEGER   IN1,IN2,IN3,IN4,IN5,IN6,IN7,IN8,IN9,IN10,IN11
      INTEGER   IM1,IM2,IM3
      INTEGER   INS1,INS2,INS3,INS4,INS5,INS6,INS7
      INTEGER   ISOLV,IOUT,N,NTERM,NUMDIR,ITMXCG
      INTEGER   IBOT,MINBOT,MAXBOT,IERSYM,NITER
      INTEGER   NODDIR(*),IA(*),JA(*),TOPOL(*),INSYM(*)
      REAL*8    ERR,RES,RHSN,RESINI
      REAL*8    TOLCG,RMIN,RMAX
      REAL*8    COEF1(*),COEF2(*),COEF3(*),SCR(*),RNSYM(*)
      REAL*8    PNEW(*),TNOTI(*)
      integer   ial(*),jal(*),iau(*),jau(*),in12
      real*8    coefl(*),coefu(*)
      LOGICAL   FATSYM
      DATA      FATSYM/.FALSE./
      RMAX=RNSYM(1)
CCC      write(6,*)'       NSYT  (ilut  -  pnew = K**(-1)b)    '
         IN1=1
         IN2=IN1+N
         IN3=IN2+N
         IN4=IN3+N
         IN5=IN4+N
         IN6=IN5+N
         IN7=IN6+N
         IN8=IN7+N
         IN9=IN8+N
         IN10=IN9+N
         IN11=IN10+N
         in12=in11+n
         IM1=1
         IM2=IM1+5*N
         IM3=IM2+5*N
C
      IF (ISOLV .LE. 2) THEN
C
C  calculates the norm of the RHS vector
C
         RHSN=0.0D0
         DO I =1,N
            RHSN=RHSN+TNOTI(I)*TNOTI(I)
         END DO
         IF (RHSN .EQ. 0.0D0) THEN
             DO I=1,N
               PNEW(I)=0.0D0
             END DO
             RETURN
         END IF
         if(isolv.ge.-1) then
            CALL ilutme(meth,n,ncoef,coef1,ja,topol,lfil,tolilu,alu,jlu,
     1                  ju,nwk,vv,iw,ierr,ac,jac,topolc,
     2                  ial,jal,coefl,iau,jau,coefu,scr(in12))
CCC           write(6,*)'n. el L=', ial(n+1)-1,'   n. el. U=',iau(n)
            if(iau(n).gt.maxlu) then
               write(iout,*)' aumentare il valore maxlu'
               stop
            endif
            call yprdly(scr,tnoti,n,ial,jal,coefl,scr(in12))
            call yprduy(pnew,scr,n,iau,jau,coefu)
         endif
         if((isolv.eq.-5).or.(isolv.eq.-3)) then
            call diagn(n,topol,ia,ja,coef1,scr(in11))
            call prdy(pnew,tnoti,n,scr(in11))
         endif
         if( (isolv.eq.-4).or.(isolv.eq.-2)) then
            call setone(numdir,noddir,topol,ja,coef1,tnoti,rmax)
            write(*,*)  'rmax'
            write(*,'(5e15.6)') rmax
            write(*,*)  'tnoti'
            write(*,'(5e15.6)') (tnoti(i),i=1,n)
            do i=1,n
               pnew(I)=tnoti(i)
            end do
         endif
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1               NUMDIR,NODDIR,RESINI,RMIN)
CCC          write(6,143) res
CCC 143      format(' initial residual=',e12.5)
         IF(ISOLV.EQ.-5) THEN
            CALL GCSTAD(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN9),
     3                  NUMDIR,NODDIR,NITER,ERR,SCR(IN11))
         ELSE IF(ISOLV.EQ.-4) THEN
            CALL GCSTAN(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN9),
     3                  NUMDIR,NODDIR,NITER,ERR)
         ELSE IF(ISOLV.EQ.-3) THEN
            CALL TFQWD(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                 SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                 SCR(IN9),TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN10),
     3                 NUMDIR,NODDIR,NITER,ERR,SCR(IN11))
         ELSE IF(ISOLV.EQ.-2) THEN
            CALL TFQWN(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                 SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                 SCR(IN9),TOPOL,JA,COEF1,ITMXCG,TOLCG,SCR(IN10),
     3                 NUMDIR,NODDIR,NITER,ERR)
         ELSE IF(ISOLV.EQ.-1) THEN
            CALL TFQWut(IOUT,N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                SCR(IN9),TOPOL,JA,COEF1,NTERM,ITMXCG,
     3                TOLCG,SCR(IN10),SCR(IN11),NUMDIR,NODDIR,NITER,
     4                ERR,
     5                ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF(ISOLV.EQ.0) THEN
            CALL GCSTut(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),SCR(IN5),SCR(IN6),SCR(IN7),SCR(IN8),
     2                  TOPOL,JA,COEF1,NTERM,
     3                  ITMXCG,TOLCG,SCR(IN9),SCR(IN10),
     4                  NUMDIR,NODDIR,NITER,ERR,
     5                  ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF (ISOLV .EQ. 1) THEN
            CALL GRAMut(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                  SCR(IN4),JA,COEF1,SCR(IN5),NTERM,TOPOL,
     2                  ITMXCG,TOLCG,SCR(IN6),SCR(IN7),NUMDIR,NODDIR,
     3                  NITER,ERR,
     4                  ial,jal,coefl,iau,jau,coefu,scr(in12))
         ELSE IF(ISOLV.EQ.2) THEN
            CALL GCRKut(N,PNEW,TNOTI,SCR(IN1),SCR(IN2),SCR(IN3),
     1                SCR(IN4),JA,COEF1,SCR(IN5),NTERM,TOPOL,
     2                ITMXCG,TOLCG,SCR(IN6),SCR(IN7),5,
     3                COEF3(IM1),COEF3(IM2),COEF3(IM3),NUMDIR,NODDIR,
     4                NITER,ERR,
     5                ial,jal,coefl,iau,jau,coefu,scr(in12))
         ENDIF
C
C  calculate the residual
C
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1               NUMDIR,NODDIR,RES,RMIN)
         WRITE(IOUT,1020) NITER,RESINI,ERR,RES,IAL(N+1)-1,IAU(N)
      ELSE
         INS1=1
         INS2=INS1+N
         INS3=INS2+N
         INS4=INS3+N+1
         INS5=INS4+N
         INS6=INS5+N
         INS7=INS6+N
C
C  in using the direct solver (NONSYM) to solve the system, we
C  need to perform symbolic factorization, and we do so only once
C
         IF (.NOT. FATSYM) THEN
            FATSYM=.TRUE.
            DO I=1,N
               J=N+I
               INSYM(I)=I
               INSYM(J)=I
            END DO
            IBOT=MAXBOT
            MINBOT=IBOT
            NFASIN=N
            CALL FASINS(NFASIN,IBOT,INSYM(INS1),INSYM(INS2),TOPOL,JA,
     1                  INSYM(INS3),INSYM(INS4),INSYM(INS7),K,IERSYM,
     2                  INSYM(INS5),INSYM(INS6))
            IF (IERSYM .NE. 0) THEN
               WRITE(IOUT,1300) IERSYM,NFASIN,N
               RETURN
            END IF
            MINBOT=INSYM(INS3+N)+N
         END IF
         CALL FANUNS(N,INSYM(INS1),INSYM(INS2),TOPOL,JA,COEF1,
     1               INSYM(INS3),INSYM(INS4),INSYM(INS7),SCR(INS1),
     2               RNSYM,IERSYM,SCR(INS2))
         IF (IERSYM .NE. 0) THEN
            WRITE(IOUT,1310) IERSYM
            RETURN
         END IF
         CALL AVINNS(N,INSYM(INS1),INSYM(INS2),INSYM(INS3),INSYM(INS4),
     1               INSYM(INS7),SCR(INS1),RNSYM,TNOTI,PNEW)
C
C  calculate the residual
C
         CALL RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(INS1),N,
     1               NUMDIR,NODDIR,RES,RMIN)
         WRITE(IOUT,1010) RES
      END IF
C
      RETURN
 1000 FORMAT(1X,I4,2(1PE15.6),'  <<NONSYMMETRIC SOLVER>>')
 1010 FORMAT(20X,1PE15.6,'  <<NONSYMMETRIC SOLVER>>')
 1020 FORMAT(1X,I4,3(1PE15.6),2I10,
     1       '  (RESINI,ERR,RES,N.EL. L,N.EL. U)')
 1300 FORMAT(//,1X,' ERROR FROM NONSYM ROUTINE FASINS, CODE ',I4,
     1        /,6X,I7,' NODES PROCESSED SO FAR OUT OF ',I7)
 1310 FORMAT(//,1X,' ERROR FROM NONSYM ROUTINE FANUNS, CODE ',I4)
      END
C  MA MATRICE E' MEMORIZZATA IN CSR FORMAT.
c  La presente routine ordina in modo crescente, dentro ciascuna riga,
C IL VETTORE  JA, E IN MANIERA COERENTE VIENE ORDINATO IL VETT. A
       subroutine ordja(n,topol,ja,a)
       implicit none
       real*8 a(*),av(200)
       integer n,topol(*),ja(*),jv(200),i1,i2,i,j
       do 10 i=1,n
       i1=0
       do 20 J=topol(i),topol(i+1)-1
       i1=i1+1
       jv(i1)=ja(j)
       av(i1)=a(j)
20     continue
c  ordina jv in maniera crescente e in maniera coerente av
       CALL IORDY(I1,JV,AV)
       i2=0
       do 30 j=topol(i),topol(i+1)-1
       i2=i2+1
       ja(j)=jv(i2)
       a(j)=av(i2)
30     continue
10     continue
       return
       end
c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c                   ITERATIVE SOLVERS MODULE                           c
c----------------------------------------------------------------------c
c contents: (May 7, 1990)  -- revised ilut Feb, 11, 1992               c
c----------                                                            c
c                                                                      c
c pgmres  : preconditioned GMRES solver                                c
c ilut    : incomplete LU factorization with dual truncation strategy  c
c ilu0    : simple ILU(0) preconditioning                              c
c milu0   : MILU(0) preconditioning                                    c
c amux    : SPARSKIT routine for performing matrix-by-vector products  c
c lusol0  : forward followed by backward triangular solve (Precond.)   c
c qsplit  : quick split routine used by ilut to sort out the k largest c
c           elements in absolute value                                 c
c                                                                      c
c Note: all preconditioners are preprocessors to pgmres.               c
c usage: call preconditioner then call pgmres.                         c
c                                                                      c
c----------------------------------------------------------------------c
       subroutine pgmres (n, im, rhs, sol, vv, eps, maxits, iout,
     *                    aa, ja, ia, alu, jlu, ju, ierr, RO)
c-----------------------------------------------------------------------
       implicit real*8 (a-h,o-z)
       integer n, im, maxits, iout, ierr, ja(*), ia(n+1), jlu(*), ju(n)
       real*8 vv(n,*), rhs(n), sol(n), aa(*), alu(*), eps
       COMMON /STAT/ITS
c----------------------------------------------------------------------*
c                                                                      *
c                 *** ILUT - Preconditioned GMRES ***                  *
c                                                                      *
c----------------------------------------------------------------------*
c This is a simple version of the ILUT preconditioned GMRES algorithm. *
c The ILUT preconditioner uses a dual strategy for dropping elements   *
c instead  of the usual level of-fill-in approach. See details in ILUT *
c subroutine documentation. PGMRES uses the L and U matrices generated *
c from the subroutine ILUT to precondition the GMRES algorithm.        *
c The preconditioning is applied to the right. The stopping criterion  *
c utilized is based simply on reducing the residual norm by epsilon.   *
c This preconditioning is more reliable than ilu0 but requires more    *
c storage. It seems to be much less prone to difficulties related to   *
c strong nonsymmetries in the matrix. We recommend using a nonzero tol *
c (tol=.005 or .001 usually give good results) in ILUT. Use a large    *
c lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the    *
c more reliable the code is. Efficiency may also be much improved.     *
c Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as *
c Gaussian elimination without pivoting.                               *
c                                                                      *
c ILU(0) and MILU(0) are also provided for comparison purposes         *
c USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and *
c then call pgmres.                                                    *
c----------------------------------------------------------------------*
c Coded by Y. Saad - This version dated May, 7, 1990.                  *
c----------------------------------------------------------------------*
c parameters                                                           *
c-----------                                                           *
c on entry:                                                            *
c==========                                                            *
c                                                                      *
c n     == integer. The dimension of the matrix.                       *
c im    == size of krylov subspace:  should not exceed 50 in this      *
c          version (can be reset by changing parameter command for     *
c          kmax below)                                                 *
c rhs   == real vector of length n containing the right hand side.     *
c          Destroyed on return.                                        *
c sol   == real vector of length n containing an initial guess to the  *
c          solution on input. approximate solution on output           *
c eps   == tolerance for stopping criterion. process is stopped        *
c          as soon as ( ||.|| is the euclidean norm):                  *
c          || current residual||/||initial residual|| <= eps           *
c maxits== maximum number of iterations allowed                        *
c iout  == output unit number number for printing intermediate results *
c          if (iout .le. 0) nothing is printed out.                    *
c                                                                      *
c aa, ja,                                                              *
c ia    == the input matrix in compressed sparse row format:           *
c          aa(1:nnz)  = nonzero elements of A stored row-wise in order *
c          ja(1:nnz) = corresponding column indices.                   *
c          ia(1:n+1) = pointer to beginning of each row in aa and ja.  *
c          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)  *
c                                                                      *
c alu,jlu== A matrix stored in Modified Sparse Row format containing   *
c           the L and U factors, as computed by subroutine ilut.       *
c                                                                      *
c ju     == integer array of length n containing the pointers to       *
c           the beginning of each row of U in alu, jlu as computed     *
c           by subroutine ILUT.                                        *
c                                                                      *
c on return:                                                           *
c==========                                                            *
c sol   == contains an approximate solution (upon successful return).  *
c ierr  == integer. Error message with the following meaning.          *
c          ierr = 0 --> successful return.                             *
c          ierr = 1 --> convergence not achieved in itmax iterations.  *
c          ierr =-1 --> the initial guess seems to be the exact        *
c                       solution (initial residual computed was zero)  *
c                                                                      *
c----------------------------------------------------------------------*
c                                                                      *
c work arrays:                                                         *
c=============                                                         *
c vv    == work array of length  n x (im+1) (used to store the Arnoli  *
c          basis)                                                      *
c----------------------------------------------------------------------*
c subroutines called :                                                 *
c amux   : SPARSKIT routine to do the matrix by vector multiplication  *
c          delivers y=Ax, given x  -- see SPARSKIT/BLASSM/amux         *
c lusol0 : combined forward and backward solves (Preconditioning ope.) *
c BLAS1  routines.                                                     *
c----------------------------------------------------------------------*
       parameter (kmax=50)
       real*8 hh(kmax+1,kmax), c(kmax), s(kmax), rs(kmax+1),t
c-------------------------------------------------------------
c arnoldi size should not exceed kmax=50 in this version..
c to reset modify paramter kmax accordingly.
c-------------------------------------------------------------
       data epsmac/1.d-16/
       n1 = n + 1
       its = 0
c-------------------------------------------------------------
c outer loop starts here..
c-------------- compute initial residual vector --------------
       call amux (n, sol, vv, aa, ja, ia)
       do 21 j=1,n
          vv(j,1) = rhs(j) - vv(j,1)
 21    continue
c-------------------------------------------------------------
 20    ro = dnrm2(n, vv, 1)
       if (iout .gt. 0 .and. its .eq. 0)
     *      write(iout, 199) its, ro
       if (ro .eq. 0.0d0) goto 999
       t = 1.0d0/ ro
       do 210 j=1, n
          vv(j,1) = vv(j,1)*t
 210   continue
       if (its .eq. 0) eps1=eps*ro
c     ** initialize 1-st term  of rhs of hessenberg system..
       rs(1) = ro
       i = 0
 4     i=i+1
       its = its + 1
       i1 = i + 1
       call lusol0 (n, vv(1,i), rhs, alu, jlu, ju)
       call amux (n, rhs, vv(1,i1), aa, ja, ia)
c-----------------------------------------
c     modified gram - schmidt...
c-----------------------------------------
       do 55 j=1, i
          t = ddot(n, vv(1,j),1,vv(1,i1),1)
          hh(j,i) = t
          call daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1)
 55    continue
       t = dnrm2(n, vv(1,i1), 1)
       hh(i1,i) = t
       if ( t .eq. 0.0d0) goto 58
       t = 1.0d0/t
       do 57  k=1,n
          vv(k,i1) = vv(k,i1)*t
 57    continue
c
c     done with modified gram schimd and arnoldi step..
c     now  update factorization of hh
c
 58    if (i .eq. 1) goto 121
c--------perfrom previous transformations  on i-th column of h
       do 66 k=2,i
          k1 = k-1
          t = hh(k1,i)
          hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
          hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
 66    continue
 121   gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
c
c     if gamma is zero then any small value will do...
c     will affect only residual estimate
c
       if (gam .eq. 0.0d0) gam = epsmac
c
c     get  next plane rotation
c
       c(i) = hh(i,i)/gam
       s(i) = hh(i1,i)/gam
       rs(i1) = -s(i)*rs(i)
       rs(i) =  c(i)*rs(i)
c
c     detrermine residual norm and test for convergence-
c
       hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
       ro = abs(rs(i1))
 131   format(1h ,2e14.4)
       if (iout .gt. 0)
     *      write(iout, 199) its, ro
       if (i .lt. im .and. (ro .gt. eps1))  goto 4
c
c     now compute solution. first solve upper triangular system.
c
       rs(i) = rs(i)/hh(i,i)
       do 30 ii=2,i
          k=i-ii+1
          k1 = k+1
          t=rs(k)
          do 40 j=k1,i
             t = t-hh(k,j)*rs(j)
 40       continue
          rs(k) = t/hh(k,k)
 30    continue
c
c     form linear combination of v(*,i)'s to get solution
c
       t = rs(1)
       do 15 k=1, n
          rhs(k) = vv(k,1)*t
 15    continue
       do 16 j=2, i
          t = rs(j)
          do 161 k=1, n
             rhs(k) = rhs(k)+t*vv(k,j)
 161      continue
 16    continue
c
c     call preconditioner.
c
       call lusol0 (n, rhs, rhs, alu, jlu, ju)
       do 17 k=1, n
          sol(k) = sol(k) + rhs(k)
 17    continue
c
c     restart outer loop  when necessary
c
       if (ro .le. eps1) goto 990
       if (its .gt. maxits) goto 991
c
c     else compute residual vector and continue..
c
       do 24 j=1,i
          jj = i1-j+1
          rs(jj-1) = -s(jj-1)*rs(jj)
          rs(jj) = c(jj-1)*rs(jj)
 24    continue
       do 25  j=1,i1
          t = rs(j)
          if (j .eq. 1)  t = t-1.0d0
          call daxpy (n, t, vv(1,j), 1,  vv, 1)
 25    continue
 199   format('   its =', i4, ' res. norm =', d20.6)
c     restart outer loop.
       goto 20
 990   ierr = 0
       return
 991   ierr = 1
       return
 999   continue
       ierr = -1
       return
c-----------------end of pgmres ---------------------------------------
c-----------------------------------------------------------------------
       end
c  prdy esegue il precondizionamento diagonale
C   PRDY dovrebbe ESEGUireE IL PRODOTTO D**(-1)*B
c si suppone tuttavia che in d siano gia' memorizzati i reciproci
c  degli elementi diagonali
      SUBROUTINE PRDY(X,B,N,D)
      implicit none
      REAL*8 X(1),B(1),D(1)
      integer n,i
      DO 10 I=1,N
      X(I)=B(I)*D(I)
10    CONTINUE
      RETURN
      END
C    ***********     ALTRE ROUTINES DI SERVIZIO
C    CALCOLA I PUNTATORI ALLA TRIANGOLARE BASSA E ALTA
C    MEMORIZZA IN DL GLI ELEMNTI DIAGONALI DI COEFL
      SUBROUTINE PREPY(N,IA,JA,TOPOL,IAL,JAL,IAU,JAU,COEF2,COEFL,
     &                  COEFU,DL)
      IMPLICIT  REAL*8 (A-H,O-Z)
      REAL*8 COEF2(1),COEFL(1),COEFU(1),DL(1)
      INTEGER IA(1),JA(1),TOPOL(1),IAL(1),JAL(1),IAU(1),JAU(1)
      IAL(1)=1
      I1=1
      DO 10 I=1,N
      DO 20 J=TOPOL(I),TOPOL(I+1)-1
      IF(JA(J).LE.I) THEN
      JAL(I1)=JA(J)
      COEFL(I1)=COEF2(J)
      I1=I1+1
      ENDIF
20    CONTINUE
      IAL(I+1)=I1
10    CONTINUE
C  TRIANGOLARE ALTA
      IAU(1)=1
      I2=1
      DO 30 I=1,N
      DO 40 J=TOPOL(I),TOPOL(I+1)-1
      IF(JA(J).GE.I) THEN
      JAU(I2)=JA(J)
      COEFU(I2)=COEF2(J)
      I2=I2+1
      ENDIF
40    CONTINUE
      IAU(I+1)=I2
30    CONTINUE
C  CAMBIO GLI ELEMENTI DIAGONALI DI COEFU
      DO 100 I=1,N
      DO 200 J=IAU(I),IAU(I+1)-1
      IF(JAU(J).EQ.I) COEFU(J)=1.
200   CONTINUE
100   CONTINUE
C     MEMORIZZO IN DL GLI ELEMTI DIAGONALI DI COEFL
      DO 300 I=1,N
      DO 400 J=IAL(I),IAL(I+1)-1
      IF(JAL(J).EQ.I) DL(I)=COEFL(J)
400   CONTINUE
300   CONTINUE
      RETURN
      END
c
c******************************  PRL  ***********************
c
C     CALCOLO IL PRODOTTO L*X=Y
c
      SUBROUTINE PRL(Y,X,N,IA,JA,COEF2,NTERM,TOPOL)
      implicit real*8 (a-h,o-z)
      real*8 COEF2(*)
      INTEGER*4 JA(*),IA(*)
      INTEGER*4 TOPOL(*)
      REAL*8 X(*),Y(*)
c
      DO 45 I=1,N
45    Y(I)=0.
       DO 6 K=1,NTERM
       IF(IA(K).LT.JA(K)) GOTO 6
       Y(IA(K))=Y(IA(K))+COEF2(K)*X(JA(K))
6      CONTINUE
       RETURN
       END
c
c******************  PROD1  **************************
c
      SUBROUTINE PROD1(X,B,N,IA,JA,COEF2,NTERM,TOPOL)
      implicit real*8 (a-h,o-z)
      real*8 COEF2(*)
      INTEGER*4 IA(*),JA(*)
      INTEGER*4 TOPOL(*)
      REAL*8 X(*),B(*)
c
5     X(1)=B(1)/COEF2(1)
      A=0.
      DO 6 K=2,NTERM
      IF(JA(K)-IA(K))11,12,6
11    A=A+COEF2(K)*X(JA(K))
      GO TO 6
12    X(IA(K))=(B(IA(K))-A)/COEF2(K)
      A=0.
6     CONTINUE
      A=0.
      M=NTERM-TOPOL(N)+2
      DO 7 K=M,NTERM
      I=NTERM-K+1
      IF(IA(I)-JA(I)) 1,2,7
1     A=A+COEF2(I)*X(JA(I))
      GO TO 7
2     X(IA(I))=X(IA(I))-A
      A=0.
7     CONTINUE
8     RETURN
      END
c
c*********************  PRODBH  ********************
c
C     POSTO B=L**(-1)*A*U**(-1), ESEGUE IL PRODOTTO DI B PER UN VETTORE
C     H E LO METTE IN Z. USA DUE VETTORI DI SERVIZIO X E Y CHE SI
C     POTERBBERO RIDURRE A UNO SOLO
c
      SUBROUTINE PRODBH(Z,H,N,IA,JA,COEF1,COEF2,NTERM,TOPOL,X,Y)
      implicit real*8 (a-h,o-z)
      real*8 COEF1(*),COEF2(*)
      REAL*8 Z(*),H(*)
      REAL*8 X(*),Y(*)
      INTEGER*4 IA(*),JA(*)
      INTEGER*4 TOPOL(*)
c
C     ESEGUO IL PRODOTTO U**(-1)*H=X
c
      CALL PRODU(X,H,N,IA,JA,COEF2,NTERM,TOPOL)
c
C     ESEGUO IL PRODOTTO A*X=Y
c
      DO 10 K=1,N
      Z(K)=0.
10    Y(K)=0.
      DO 20 K=1,NTERM
      Y(IA(K))=Y(IA(K))+COEF1(K)*X(JA(K))
20    CONTINUE
c
C     ESEGUO IL PRODOTTO L**(-1)*Y=Z
c
      CALL PRODL(Z,Y,N,IA,JA,COEF2,NTERM,TOPOL)
      RETURN
      END
C********************    PRODME     *********************************
      SUBROUTINE PRODDP(X,B,N,TOPOL,JA,COEF2,NTERM)
      implicit real*8 (a-h,o-z)
      REAL*8 COEF2(*)
      INTEGER*4 JA(*)
      INTEGER*4 TOPOL(*)
      REAL*8 X(*),B(*),A
      DO 10 K=1,N
10    X(K)=0.
20    DO 26 K=1,N
      I=TOPOL(K)
      J=TOPOL(K+1)-1
      DO 26 M=I,J
      IF(M.EQ.I) GO TO 22
      X(JA(M))=X(JA(M))+COEF2(M)*X(K)
      GO TO 26
22    X(K)=(B(K)-X(K))/COEF2(M)
26    CONTINUE
      A=0.
      DO 27 K=1,N
      I=TOPOL(N-K+1)
      J=TOPOL(N-K+2)-1
      DO 27 M=I,J
      MM=J-M+I
      IF(MM.EQ.I) GO TO 23
      A=A+COEF2(MM)*X(JA(MM))
      GO TO 27
23    X(N-K+1)=(X(N-K+1)-A)/COEF2(MM)
      A=0.
27    CONTINUE
8     RETURN
      END
c
c********************* PRODL  ***********************
c
C     ESEGUE IL PRODOTTO L**(-1)*Y=Z RISOLVENDO IL SISTEMA L*Z=X
C
      SUBROUTINE PRODL(Z,Y,N,IA,JA,COEF2,NTERM,TOPOL)
      implicit real*8 (a-h,o-z)
      real*8 COEF2(*)
      REAL*8 Z(*),Y(*)
      INTEGER*4 IA(*),JA(*)
      INTEGER*4 TOPOL(*)
c
5     Z(1)=Y(1)/COEF2(1)
      A=0.
      DO 6 K=2,NTERM
      IF(JA(K)-IA(K)) 11,12,6
11    A=A+COEF2(K)*Z(JA(K))
      GOTO 6
12    Z(IA(K))=(Y(IA(K))-A)/COEF2(K)
      A=0.
6     CONTINUE
      RETURN
      END
c
c*******************************  PRODU  **********************
c
C     ESEGUE IL PRODOTTO U**(-1)*H=X RISOLVENDO IL SISTEMA U*X=H
C
      SUBROUTINE PRODU(X,H,N,IA,JA,COEF2,NTERM,TOPOL)
      implicit real*8 (a-h,o-z)
      real*8 COEF2(*)
      REAL*8 X(*),H(*)
      INTEGER*4 IA(*),JA(*)
      INTEGER*4 TOPOL(*)
      X(N)=H(N)
      A=0.
      M=NTERM-TOPOL(N)+2
      DO 7 K=M,NTERM
      I=NTERM-K+1
      IF(IA(I)-JA(I))1,2,7
1     A=A+COEF2(I)*X(JA(I))
      GOTO 7
2     X(IA(I))=H(IA(I))-A
      A=0.
7     CONTINUE
      RETURN
      END
c
c********************  PRU  *********************************
c
C     CALCOLO IL PRODOTTO U*X=Y
      SUBROUTINE PRU(Y,X,N,IA,JA,COEF2,NTERM,TOPOL)
      implicit real*8 (a-h,o-z)
      real*8 COEF2(*)
      INTEGER*4 JA(*),IA(*)
      INTEGER*4 TOPOL(*)
      REAL*8 X(*),Y(*)
c
       DO 6 K=1,NTERM
       IF(IA(K).GT.JA(K)) GOTO 6
       IF(IA(K).EQ.JA(K)) Y(IA(K))=X(IA(K))
       IF(IA(K).EQ.JA(K)) GOTO 6
       Y(IA(K))=Y(IA(K))+COEF2(K)*X(JA(K))
6      CONTINUE
       RETURN
       END
c-----------------------------------------------------------------------
        subroutine qsplit  (a, ind, n, ncut)
        real*8 a(n)
        integer ind(n), n, ncut
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c     a(i) .le. a(ncut) for i .le. ncut and
c     a(i) .ge. a(ncut) for i .ge. ncut
c     ind(1:n) is an integer array which permuted in the same way as a(*).
c-----------------------------------------------------------------------
        real*8 tmp, abskey
        integer itmp, first, last
c-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
 1      mid = first
        abskey = abs(a(mid))
        do 2 j=first+1, last
           if (abs(a(j)) .gt. abskey) then
              mid = mid+1
c     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
c
c     interchange
c
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
c
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
c
c     test for while loop
c
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
c----------------end-of-qsplit------------------------------------------
c-----------------------------------------------------------------------
        end
C
C**************************  RESNSY ************************************
C
C  calculate the residual from the conjugate gradients or NONSYM
C  solution of a nonsymmetric linear system
C
C***********************************************************************
C
      SUBROUTINE RESNSY(TOPOL,JA,COEF1,TNOTI,PNEW,RESV,N,
     1                  NUMDIR,NODDIR,RES,RMIN)
C
      IMPLICIT  NONE
      INTEGER   I,K,M,MM
      INTEGER   N,NUMDIR
      INTEGER   TOPOL(*),JA(*),NODDIR(*)
      REAL*8    XLUNG
      REAL*8    RES,RMIN
      REAL*8    COEF1(*),TNOTI(*),PNEW(*),RESV(*)
C
      DO K=1,N
         RESV(K)=TNOTI(K)
      END DO
      DO K=1,NUMDIR
         RESV(NODDIR(K))=0.0D0
      END DO
      XLUNG=0.0D0
      DO K=1,N
         XLUNG=XLUNG+RESV(K)*RESV(K)
      END DO
      DO K=1,N
         M=TOPOL(K)
         MM=TOPOL(K+1)-1
         DO I=M,MM
            RESV(K)=RESV(K)-COEF1(I)*PNEW(JA(I))
         END DO
      END DO
C
C  inizializza a 0 la variabile RES = residuo quadratico medio
C
      DO K=1,NUMDIR
         RESV(NODDIR(K))=0.0D0
      END DO
      RES=0.0D0
C
C  calcola il residuo quadratico medio
C
      DO K=1,N
         RES=RES + RESV(K)*RESV(K)
      END DO
      IF(XLUNG .LE. RMIN) THEN
         RES=DSQRT(RES/N)
      ELSE
         RES=DSQRT(RES/XLUNG)
      ENDIF
C
      RETURN
      END
C
C**************************  RESSYM ************************************
C
C  calculate the residual from the conjugate gradients
C  solution of a symmetric linear system
C
C***********************************************************************
C
      SUBROUTINE RESSYM(TOPOL,JA,COEF1,TNOTI,PNEW,RESV,N,
     1                  NUMDIR,NODDIR,RES,RMIN)
C
      IMPLICIT  NONE
      INTEGER   I,K,M,MM
      INTEGER   N,NUMDIR
      INTEGER   TOPOL(*),JA(*),NODDIR(*)
      REAL*8    XLUNG
      REAL*8    RES,RMIN
      REAL*8    COEF1(*),TNOTI(*),PNEW(*),RESV(*)
C
      DO K=1,N
         RESV(K)=TNOTI(K)
      END DO
      DO K=1,NUMDIR
         RESV(NODDIR(K))=0.0D0
      END DO
      XLUNG=0.0D0
      DO K=1,N
         XLUNG=XLUNG+RESV(K)*RESV(K)
      END DO
      DO K=1,N
         M=TOPOL(K)
         MM=TOPOL(K+1)-1
         DO I=M,MM
            IF(I.GT.M) RESV(JA(I))=RESV(JA(I))-COEF1(I)*PNEW(K)
            RESV(K)=RESV(K)-COEF1(I)*PNEW(JA(I))
         END DO
      END DO
C
C  inizializza a 0 la variabile RES = residuo quadratico medio
C
      DO K=1,NUMDIR
         RESV(NODDIR(K))=0.0D0
      END DO
      RES=0.0D0
C
C  calcola il residuo quadratico medio
C
      DO K=1,N
         RES=RES + RESV(K)*RESV(K)
      END DO
      IF(XLUNG .LE. RMIN) THEN
         RES=DSQRT(RES/N)
      ELSE
         RES=DSQRT(RES/XLUNG)
      END IF
C
      RETURN
      END
      SUBROUTINE SETONE(NUMDIR,NODDIR,TOPOL,JA,COEF1,TNOTI,RMAX)
C
      IMPLICIT  NONE
      INTEGER   NUMDIR
      INTEGER   I,K,IND,IND1
      INTEGER   NODDIR(*),TOPOL(*),JA(*)
      REAL*8    RMAX,RMINV
      REAL*8    COEF1(*),TNOTI(*)
C
      RMINV=1.0D0/RMAX
      DO I=1,NUMDIR
         IND=TOPOL(NODDIR(I))
         IND1=TOPOL(NODDIR(I)+1)-1
         DO K=IND,IND1
            COEF1(K)=0.0D0
            IF( JA(K) .EQ. NODDIR(I)) COEF1(K)=1.0D0
         END DO
         TNOTI(NODDIR(I))=TNOTI(NODDIR(I))*RMINV
      END DO
C
      RETURN
      END
C     *** SORTNS *******************************************************
C     *                                                                *
C     *   TALE ROUTINE ORDINA IN ORDINE ASCENDENTE GLI ELEMENTI DEL    *
C     *   VETTORE DI INTERI IA A PARTIRE DA IA(II) FINO A IA(JJ).      *
C     *   IL VETTORE R E' AGGIORNATO IN ACCORDO ALLA SUCCESSIONE DEGLI *
C     *   ELEMENTI DI IA A PARTIRE DALLA POSIZIONE II FINO A JJ        *
C     *   LA DIMENSIONE K DEI VETTORI INTERNI IU E IL CONSENTE DI      *
C     *   ORDINARE FINO A 2**(K+1)-1 ELEMENTI                          *
C     *   QUI K=16, COSI SORTNS PUO' ORDINARE 131071 ELEMETI           *
C     *                                                                *
C     ******************************************************************
C
      SUBROUTINE SORTNS(IA,R,II,JJ,IER)
      DIMENSION IA(1),IU(16),IL(16)
      INTEGER R(1)
C
      INERR=1
      IF (IER+12345) 2,1,2
    1 INERR=2
    2 IER  =0
      IF (II) 4,4,3
    3 IF (II-JJ) 5,4,4
    4 IER  =1000
      GOTO 29
C
C     INIZIALIZZA M= NO. DI SEGMENTI
C     I = INIZIO CORRENTE DEL SEGMENTO INFERIORE
C     J = FINE CORRENTE DEL SEMENTO SUPERIORE
C
    5 M    =1
      I    =II
      J    =JJ
C
C   NESSUNA AZIONE E' ESEGUITA PER I SEGMENTI DI LUNGHEZZA NON POSITIVA
C
    6 IF (I-J) 7,19,19
C
C     L'INDICE CORRENTE DEL SEGMENTO INFERIORE PARTE DA I
C     IT = VALORE  (LOCAZIONE) CHE SEPARA I SEGMENTI
C     SCEGLIE  IT PER RIORDINARE IA(I), IA(J), E IA(I+J/2)
C     SE NECESSARIO, IN MODO TALE CHE
C      IA(I) .LE. IA(I+J/2) .LE. IA(J) E SI PONE  IT = IA(I+J/2)
C
    7 K    =I
      IJ   =(J+I)/2
      IT   =IA(IJ)
      RT   =R(IJ)
C
C     DAPPRIMA SI SCAMBIA LA LOCAZIONE INIZIALE DEL SEGMENTO INFERIORE
C     E LA LOCAZIONE DI SEPARAZIONE SE NECESSARIO
C
      IF (IA(I)-IT) 9,9,8
    8 IA(IJ)=IA(I)
      R(IJ)=R(I)
      IA(I)=IT
      R(I) =RT
      IT   =IA(IJ)
      RT   =R(IJ)
C
C     L'INDICE CORRENTE L DEL SEGMENTO SUPERIORE PARTE DA J
C
    9 L    =J
C
C     ALLORA SCAMBIA L'ULTIMA LOCAZIONE DEL SEGMENTO SUPERIORE
C     E LA LOCAZIONE DI SEPARAZIONE SE NECESSARIO
C
      IF (IA(J)-IT) 10,13,13
   10 IA(IJ)=IA(J)
      R(IJ)=R(J)
      IA(J)=IT
      R(J) =RT
      IT   =IA(IJ)
      RT   =R(IJ)
C
C    INFINE SCAMBIA ANCORA LA LOCAZIONE INIZIALE NEL SEGMENTO INFERIORE
C    E LA LOCAZIONE DI SEPARAZIONE SE NECESSARIO PER OTTENERE
C    IA(I) .LE. IA(JJ) .LE. IA(J) E LA LOCAZIONE DI
C    SEPARAZIONE  IT=IA(JJ
C
      IF (IA(I)-IT) 13,13,11
   11 IA(IJ)=IA(I)
      R(IJ)=R(I)
      IA(I)=IT
      R(I) =RT
      IT   =IA(IJ)
      RT   =R(IJ)
      GOTO 13
C
C     SCAMBIA LA LOCAZIONE DELLA PARTE INFERIORE DEL SEGMENTO CHE
C     E' MAGGIORE DELLA LOCAZIONE DI SEPARAZIONE, CON LA LOCAZIONE
C     DELLA PARTE SUPERIORE DEL SEGMENTO CHE E' MINORE O UGUALE DELLA
C     LOCAZIONE DI SEPARAZIONE E CONTINUA A PROCESSARE QUESTO SEGMENTO
C
   12 ITT  =IA(L)
      RTT  =R(L)
      IA(L)=IA(K)
      R(L) =R(K)
      IA(K)=ITT
      R(K) =RTT
C
C     VA IN GIU' NEL SEGMENTO SUPERIORE TANTO QUANTO POSSIBILE
C    IN MODO CHE LE LOCAZIONI SIANO MAGGIORI
C    DELLE ENTRATE DI SEPARAZIONE
C
   13 L    =L-1
      IF (IA(L)-IT) 14,14,13
C
C     VA IN SU' NEL SEGMENTO INFERIORE TANTO QUANTO POSSIBILE
C    IN MODO CHE LE LOCAZIONI SIANO MINORI O UGUALI
C    DELLE ENTRATE DI SEPARAZIONE
C
   14 K    =K+1
      IF (IA(K)-IT) 14,15,15
C
C     SALTA SE ENTRAMBI I SEGMENTI HAMNO AL PIU' UNA LOCAZIONE IN
C     COMUNE
C
   15 IF (K-L) 12,12,16
C
C     SOVRAPPOSIZIONE DI ENTRAMBI I SEGMENTI SIGNIFICA COMPLETAMENTO
C    DELLA SEPARAZIONE, CONTINUA CON IL SEGMENTO PIU' CORTO E
C    MEMORIZZA IL SEGMENTO PIU' LUNGO
C
   16 IF (L-I-J+K) 18,18,17
C
C     L-I > J-K, CIOE' IL SEGMENTO INFERIORE E' PIU' LUNGO
C     PONI I=K E CONTINUA CON IL SEGMENTO SUPERIORE
C     INIZIO DELLA MEMORIZZAZIONE (IL(M)=I) E FINE (IU(M)=L)
C     DEL SEGMENTO SUPERIORE PER FUTURI PROCESSI
C
   17 IL(M)=I
      IU(M)=L
      I    =K
      M    =M+1
      GOTO 21
C
C     L-I E' UGUALE O MINORE DI J-K, CIOE' IL SEGMENTO SUPERIORE E'
C     PIU' LUNGO. PONI J=L E CONTINUA COL SEMENTO INFERIORE.
C     INIZIO DELLA MEMORIZZAZIONE (IL(M)=K) E FINE (IU(M)=J)
C     DEL SEGMENTO INFERIORE PER FUTURI PROCESSI
C
   18 IL(M)=K
      IU(M)=J
      J    =L
      M    =M+1
      GOTO 21
C
C     RIDUCE DI 1 IL NUMERO DEI SEGMENTI DA PROCESSARE
C     RETURN SE NON CI SONO PIU' SEGMENTI DA PROCESSARE
C
   19 M    =M-1
      IF (M) 20,29,20
C
C     PASSA AL PROSSIMO SEGMENTO DA  PROCESSARE
C     CONTINUA IL PROCESSO DI SEGMENTAZION PURCHE' LA LUNGHEZZA DEI
C     SEGMENTI SIA MAGGIORE DI 11
C
   20 I    =IL(M)
      J    =IU(M)
   21 IF (J-I-11) 22,7,7
   22 IF (I-II) 24,6,24
C
C     ORDINA GLI ELEMENTI DENTRO AL SEGMENTO PER MEZZO DI UN
C     INTERSCAMBIO DI COPPIE ADIACENTI SE IT NON PARTE DA II.
C
   23 I    =I+1
   24 IF (I-J) 25,19,25
   25 IT   =IA(I+1)
      RT   =R(I+1)
      IF (IA(I)-IT) 23,23,26
   26 K    =I
   27 IA(K+1)=IA(K)
      R(K+1)=R(K)
      K    =K-1
      IF (IT-IA(K)) 27,28,28
   28 IA(K+1)=IT
      R(K+1)=RT
      GOTO 23
   29 RETURN
      END
C   EMULA LA FUNCTION SPDOT DI SCILIB
C
      REAL*8 FUNCTION SPDOT(N,SY,INDEX,SX)
C
      IMPLICIT  REAL*8 (A-H,O-Z)
      REAL*8 SY(1),SX(1),S
      INTEGER INDEX(1)
      S=0.
      DO 10 I=1,N
      S=S+SY(INDEX(I))*SX(I)
10    CONTINUE
      SPDOT=S
      RETURN
      END
C
C   solutore nel caso scalare
C
C
C************************** SYMSLV ************************************
C
C  Solve a symmetric linear system of equations using a standard
C  conjugate gradients method (note: COEF2 and SCR in the conjugate
C  gradients routines are used as scratch vectors).
C  Calculate the residual from the conjugate gradients
C  solution of the symmetric linear system.
C
C***********************************************************************
C
      SUBROUTINE SYMSLV(IOUT,N,NTERM,NUMDIR,NITER,ITMXCG,TOLCG,RMIN,
     1                  NODDIR,JA,TOPOL,PNEW,TNOTI,COEF1,COEF2,SCR)
C
      IMPLICIT  NONE
      INTEGER   IN1,IN2,IN3,IN4,IN5
      INTEGER   IOUT,N,NTERM,NUMDIR,NITER,ITMXCG
      INTEGER   NODDIR(*),JA(*),TOPOL(*)
      REAL*8    ERR,RES
      REAL*8    TOLCG,RMIN
      REAL*8    PNEW(*),TNOTI(*),COEF1(*),COEF2(*),SCR(*)
C
      IN1=1
      IN2=IN1+N
      IN3=IN2+N
      IN4=IN3+N
      IN5=IN4+N
      CALL INCLDP(IOUT,JA,TOPOL,COEF1,COEF2,N,NTERM)
      CALL PRODDP(PNEW,TNOTI,N,TOPOL,JA,COEF2,NTERM)
      CALL GRADDP(N,PNEW,SCR(IN1),SCR(IN2),SCR(IN3),SCR(IN4),SCR(IN5),
     1            TOPOL,JA,COEF1,COEF2,NTERM,TNOTI,TOLCG,ITMXCG,
     2            NUMDIR,NODDIR,NITER,ERR)
C
C  calculate the residual
C
      CALL RESSYM(TOPOL,JA,COEF1,TNOTI,PNEW,SCR(IN1),N,
     1            NUMDIR,NODDIR,RES,RMIN)
      WRITE(IOUT,1000) NITER,ERR,RES
C
      RETURN
 1000 FORMAT(1X,I4,2(1PE15.6),'  <<SYMMETRIC SOLVER>>')
      END
C***********************************************************
c  versione per SOLSCAL  su risc
C   Transpose-Free QMR Algorithm  for non -Hermitian linear
C    SYSTEMS    (FREUND)  (ALGORITMO 5.2 CON PESI)
      SUBROUTINE TFQW(IOUT,N,X,TNOTI,P,U,RC,R,V,DM,Q,RT,AU,
     1                            TOPOL,ia,JA,COEF1,coef2,
     2nterm,IMAX,TOL,SCR1,SCR2,np,contp,niter,err)
      implicit none
      INTEGER N,TOPOL(*),JA(*),IMAX,IOUT
      real*8 x(*),tnoti(*),p(*),u(*),rc(*),r(*),v(*),dm(*)
      REAL*8 Q(*),RT(*),AU(*)
      REAL*8 COEF1(*),SCR1(*),SCR2(*)
      real*8 tol,ta0,tam1,tam,ron1,ron,tetm1,tetm
      real*8 sn,sn1,sim1,alf1,betn,omp1,etm,etm1,taf,cm
      real*8 coef2(*)
      integer ia(*),np,contp(*),nterm,in
      INTEGER NITER,MIT,I,MITER
      real*8 em1,tr,err
      write(99,*) 'iter, residual'
c  Calcolo Ax0
      call aperbn(n,topol,ja,coef1,x,scr1)
c  Calcolo r0=b-Ax0
      do 10 i=1,n
        scr1(i)=tnoti(i)-scr1(i)
10    continue
c  Calcolo L**(-1)r0
       call prodl(r,scr1,n,ia,ja,coef2,nterm,topol)
c   Calcolo y0=Ux0
       call pru(scr1,x,n,ia,ja,coef2,nterm,topol)
       do 15 i=1,n
       x(i)=scr1(i)
       p(i)=r(i)
       u(i)=r(i)
       rc(i)=r(i)
       dm(i)=0.
c      RT(I)=TNOTI(I)
       RT(I)=r(i)
15     continue
c  calcolo v0=Ay1
      call prodbh(v,p,n,ia,ja,coef1,coef2,nterm,topol,scr1,scr2)
c---------
      if(np.eq.0)  goto 49
       do 265 i=1,np
         in=contp(i)
         rc(in)=0.0
265    continue
49     continue
c------
      ron1=0.
      tetm1=0.
      etm1=0.
      sn1=0.
      do 20 i=1,n
      sn1=sn1+rc(i)**2
      ron1=ron1+rt(i)*r(i)
20    continue
      sn1=dsqrt(sn1)
      ta0=sn1
      tam1=sn1
      if(ron1.eq.0) then
        write(iout,*)' r0=',ron1,' cambiare r tilde'
         stop
      endif
      niter=0
      miter=0
30    continue
      niter=niter+1
      sim1=0.
      do 35 i=1,n
      sim1=sim1+rt(i)*v(i)
35    continue
      alf1=ron1/sim1
      do 40 i=1,n
      q(i)=u(i)-alf1*v(i)
40    continue
      do 50 i=1,n
      v(i)=u(i)+q(i)
50    continue
c   calcolo B(u_n-1 +q)
      call prodbh(au,v,n,ia,ja,coef1,coef2,nterm,topol,scr1,scr2)
      do 55 i=1,n
       rc(i)=rc(i)-alf1*au(i)
55    continue
c  calcolo  ||rn||
c---------
      if(np.eq.0)  goto 59
       do 365 i=1,np
         in=contp(i)
         rc(in)=0.0
365    continue
59     continue
c------
      sn=0.
      do 60 I=1,n
      sn=sn+rc(i)**2
60    continue
      sn=dsqrt(sn)
      mit=0
70    mit=mit+1
      miter=miter+1
      if(mit.eq.1) then
       omp1=dsqrt(sn1*sn)
      else
       omp1=sn
      endif
      tetm=omp1/tam1
      cm=1./dsqrt(1+tetm**2)
      tam=tam1*tetm*cm
      etm=cm**2*alf1
      taf=tetm1**2*etm1/alf1
       if(mit.eq.1) then
        do 100 i=1,n
        dm(i)=u(i)+taf*dm(i)
100     continue
       else
        do 110 i=1,n
        dm(i)=q(I)+taf*dm(i)
110    continue
       endif
       do 120 i=1,n
       x(i)=x(i)+etm*dm(i)
120    continue
c     aggiornamento
       tetm1=tetm
       tam1=tam
       etm1=etm
       em1=miter+1
c      tr=(ta0*tol)/dsqrt(em1)
       err=dsqrt(em1)*tam/ta0
       if((niter.ge.imax).or.(err.le.tol))goto 800
          write(99,*) niter,err
ccccciter             write(6,138)niter,err
138     format(i5,3x,'err=',e12.5)
        if(mit.eq.1) goto 70
c   esecuzione del punto c)
        ron=0.
        do 200 i=1,n
        ron=ron+rt(i)*rc(i)
200     continue
        betn=ron/ron1
        do 210 i=1,n
        u(i)=rc(i)+betn*q(i)
210     continue
        do 220 i=1,n
        p(i)=u(i)+betn*(q(i)+betn*p(i))
220     continue
      call prodbh(v,p,n,ia,ja,coef1,coef2,nterm,topol,scr1,scr2)
c   aggiornamento
         sn1=sn
         ron1=ron
         goto 30
800      continue
         do 520 i=1,n
         scr1(i)=x(i)
         x(i)=0.0
520      continue
c  Trovo x=U**(-1)y
         call produ(x,scr1,n,ia,ja,coef2,nterm,topol)
           return
           end
c***********************************************************
C--------------------------------------------------
C         PRECONDIZIONAMENTO DIAGONALE
C   Transpose-Free QMR Algorithm  for non -Hermitian linear
C    SYSTEMS    (FREUND)  (ALGORITMO 5.2 CON PESI)
      SUBROUTINE TFQWD(IOUT,N,X,TNOTI,P,U,RC,R,V,DM,Q,RT,AU,
     1TOPOL,JA,COEF1,IMAX,TOL,SCR1,
     2 NP,CONTP,NITER,ERR,DD)
      implicit none
      INTEGER N,TOPOL(*),JA(*),IMAX,J
      real*8 x(*),tnoti(*),p(*),u(*),rc(*),r(*),v(*),dm(*)
      REAL*8 Q(*),RT(*),AU(*),DD(*),SCR1(*)
      REAL*8 COEF1(*)
      real*8 tol,ta0,tam1,tam,ron1,ron,tetm1,tetm
      real*8 sn,sn1,sim1,alf1,betn,omp1,etm,etm1,taf,cm
      REAL*8 SS,EM1,TR,ERR
      INTEGER NITER,MITER,MIT,I,K
      INTEGER NP, CONTP(*),IN,IOUT
C-------------------
      write(99,*) 'iter, residual'
c  Calcolo Ax0
      call aperbn(n,topol,ja,coef1,x,scr1)
c  Calcolo r0=b-Ax0
      do 10 i=1,n
        scr1(i)=tnoti(i)-scr1(i)
10    continue
c  Calcolo L**(-1)r0
       DO 7 I=1,N
          R(I)=SCR1(I) *DD(I)
7      CONTINUE
C   CALCOLO Y0=Ux0=Ix0
       do 15 i=1,n
       p(i)=r(i)
       u(i)=r(i)
       rc(i)=r(i)
       dm(i)=0.
       RT(I)=R(I)
15     continue
C  CALCOLO V0=BY1
      CALL APERBN(N,TOPOL,JA,COEF1,P,V)
        do 650 i=1,n
         v(I)=v(I)*dd(I)
650     continue
c---------
      if(np.eq.0)  goto 49
       do 265 i=1,np
         in=contp(i)
         rc(in)=0.0
265    continue
49     continue
c------
      ron1=0.
      tetm1=0.
      etm1=0.
      sn1=0.
      do 20 i=1,n
      sn1=sn1+rc(i)**2
      ron1=ron1+rt(i)*r(i)
20    continue
      sn1=dsqrt(sn1)
      ta0=sn1
      tam1=sn1
      if(ron1.eq.0) then
        WRITE(IOUT,*)' R0=',RON1,' CAMBIARE R TILDE'
         stop
      endif
      niter=0
      miter=0
30    continue
      niter=niter+1
      sim1=0.
      do 35 i=1,n
      sim1=sim1+rt(i)*v(i)
35    continue
      alf1=ron1/sim1
      do 40 i=1,n
      q(i)=u(i)-alf1*v(i)
40    continue
      do 50 i=1,n
      V(I)=U(I)+Q(I)
50    continue
c   calcolo B(u_n-1 +q)
      CALL APERBN(N,TOPOL,JA,COEF1,V,AU)
        do 750 i=1,n
         au(I)=au(I)*dd(i)
750    continue
      do 55 i=1,n
       rc(i)=rc(i)-alf1*au(i)
55    continue
c  calcolo  ||rn||
c---------
      if(np.eq.0)  goto 59
       do 365 i=1,np
         in=contp(i)
         rc(in)=0.0
365    continue
59     continue
c------
      sn=0.
      do 60 I=1,n
      sn=sn+rc(i)**2
60    continue
      sn=dsqrt(sn)
      mit=0
70    mit=mit+1
      miter=miter+1
      if(mit.eq.1) then
       omp1=dsqrt(sn1*sn)
      else
       omp1=sn
      endif
      tetm=omp1/tam1
      cm=1./dsqrt(1+tetm**2)
      tam=tam1*tetm*cm
      etm=cm**2*alf1
      taf=tetm1**2*etm1/alf1
       if(mit.eq.1) then
        do 100 i=1,n
        dm(i)=u(i)+taf*dm(i)
100     continue
       else
        do 110 i=1,n
        dm(i)=q(I)+taf*dm(i)
110    continue
       endif
       do 120 i=1,n
       x(i)=x(i)+etm*dm(i)
120    continue
c     aggiornamento
       tetm1=tetm
       tam1=tam
       etm1=etm
       em1=miter+1
C      TR=(TA0*TOL)/DSQRT(EM1)
       ERR=DSQRT(EM1)*TAM/TA0
       IF((NITER.GE.IMAX).OR.(ERR.LE.TOL)) GOTO 800
        if(mit.eq.1) goto 70
c   esecuzione del punto c)
        ron=0.
        do 200 i=1,n
        ron=ron+rt(i)*rc(i)
200     continue
        betn=ron/ron1
        do 210 i=1,n
        u(i)=rc(i)+betn*q(i)
210     continue
        do 220 i=1,n
        p(i)=u(i)+betn*(q(i)+betn*p(i))
220     continue
      CALL APERBN(N,TOPOL,JA,COEF1,P,V)
         do 850 i=1,n
          v(I)=v(i)*dd(i)
850      continue
c   aggiornamento
         sn1=sn
         ron1=ron
         write(99,*) niter,err
         goto 30
800      continue
           return
           end
c********************************************************
C   SENZA PRECONDIZIONAMENTO
C   Transpose-Free QMR Algorithm  for non -Hermitian linear
C    SYSTEMS    (FREUND)  (ALGORITMO 5.2 CON PESI)
      SUBROUTINE TFQWN(IOUT,N,X,TNOTI,P,U,RC,R,V,DM,Q,RT,AU,
     1TOPOL,JA,COEF1,IMAX,TOL,SCR1,
     2 NP,CONTP,NITER,ERR)
      implicit none
      INTEGER N,TOPOL(*),JA(*),IMAX
      real*8 x(*),tnoti(*),p(*),u(*),rc(*),r(*),v(*),dm(*)
      REAL*8 Q(*),RT(*),AU(*),SCR1(*)
      REAL*8 COEF1(*)
      real*8 tol,ta0,tam1,tam,ron1,ron,tetm1,tetm
      real*8 sn,sn1,sim1,alf1,betn,omp1,etm,etm1,taf,cm
      REAL*8 SS,EM1,TR,ERR
      INTEGER NITER,MITER,MIT,I,K
      INTEGER NP, CONTP(*),IN,IOUT
      write(99,*) 'iter, residual'
c  Calcolo Ax0
      call aperbn(n,topol,ja,coef1,x,scr1)
c  Calcolo r0=b-Ax0
      do 10 i=1,n
        scr1(i)=tnoti(i)-scr1(i)
10    continue
c  Calcolo L**(-1)r0
       DO 7 I=1,N
          R(I)=SCR1(I)
7      CONTINUE
c   Calcolo y0=Ux0
       do 15 i=1,n
       p(i)=r(i)
       u(i)=r(i)
       rc(i)=r(i)
       dm(i)=0.
       RT(I)=R(I)
15     continue
C  CALCOLO V0=BY1
      CALL APERBN(N,TOPOL,JA,COEF1,P,V)
c---------
      if(np.eq.0)  goto 49
       do 265 i=1,np
         in=contp(i)
         rc(in)=0.0
265    continue
49     continue
c------
      ron1=0.
      tetm1=0.
      etm1=0.
      sn1=0.
      do 20 i=1,n
      sn1=sn1+rc(i)**2
      ron1=ron1+rt(i)*r(i)
20    continue
      sn1=dsqrt(sn1)
      ta0=sn1
      tam1=sn1
      if(ron1.eq.0) then
        WRITE(IOUT,*)' R0=',RON1,' CAMBIARE R TILDE'
         stop
      endif
      niter=0
      miter=0
30    continue
      niter=niter+1
      sim1=0.
      do 35 i=1,n
      sim1=sim1+rt(i)*v(i)
35    continue
      alf1=ron1/sim1
      do 40 i=1,n
      q(i)=u(i)-alf1*v(i)
40    continue
      do 50 i=1,n
      V(I)=U(I)+Q(I)
50    continue
c   calcolo B(u_n-1 +q)
      CALL APERBN(N,TOPOL,JA,COEF1,V,AU)
      do 55 i=1,n
       rc(i)=rc(i)-alf1*au(i)
55    continue
c  calcolo  ||rn||
c---------
      if(np.eq.0)  goto 59
       do 365 i=1,np
         in=contp(i)
         rc(in)=0.0
365    continue
59     continue
c------
      sn=0.
      do 60 I=1,n
      sn=sn+rc(i)**2
60    continue
      sn=dsqrt(sn)
      mit=0
70    mit=mit+1
      miter=miter+1
      if(mit.eq.1) then
       omp1=dsqrt(sn1*sn)
      else
       omp1=sn
      endif
      tetm=omp1/tam1
      cm=1./dsqrt(1+tetm**2)
      tam=tam1*tetm*cm
      etm=cm**2*alf1
      taf=tetm1**2*etm1/alf1
       if(mit.eq.1) then
        do 100 i=1,n
        dm(i)=u(i)+taf*dm(i)
100     continue
       else
        do 110 i=1,n
        dm(i)=q(I)+taf*dm(i)
110    continue
       endif
       do 120 i=1,n
       x(i)=x(i)+etm*dm(i)
120    continue
c     aggiornamento
       tetm1=tetm
       tam1=tam
       etm1=etm
       em1=miter+1
C      TR=(TA0*TOL)/DSQRT(EM1)
       ERR=DSQRT(EM1)*TAM/TA0
       IF((NITER.GE.IMAX).OR.(ERR.LE.TOL)) GOTO 800
        if(mit.eq.1) goto 70
c   esecuzione del punto c)
        ron=0.
        do 200 i=1,n
        ron=ron+rt(i)*rc(i)
200     continue
        betn=ron/ron1
        do 210 i=1,n
        u(i)=rc(i)+betn*q(i)
210     continue
        do 220 i=1,n
        p(i)=u(i)+betn*(q(i)+betn*p(i))
220     continue
      CALL APERBN(N,TOPOL,JA,COEF1,P,V)
c   aggiornamento
         sn1=sn
         ron1=ron
         write(99,*) niter,err
         goto 30
800      continue
           return
           end
c
c  Conjugate gradient method with CHOLESKY preconditioner
c
c      DOUBLE PRECISION
c
C***********************************************************
c  versione per SOLSCAL  su risc
C   Transpose-Free QMR Algorithm  for non -Hermitian linear
C    SYSTEMS    (FREUND)  (ALGORITMO 5.2 CON PESI)
      SUBROUTINE TFQWut(IOUT,N,X,TNOTI,P,U,RC,R,V,DM,Q,RT,AU,
     1                            TOPOL,JA,COEF1,
     2     nterm,IMAX,TOL,SCR1,SCR2,np,contp,niter,err,
     3     ial,jal,coefl,iau,jau,coefu,dl)
      implicit none
      real*8 coefl(*),coefu(*),dl(*)
      integer ial(*),jal(*),iau(*),jau(*)
      INTEGER N,TOPOL(*),JA(*),IMAX,IOUT
      real*8 x(*),tnoti(*),p(*),u(*),rc(*),r(*),v(*),dm(*)
      REAL*8 Q(*),RT(*),AU(*)
      REAL*8 COEF1(*),SCR1(*),SCR2(*)
      real*8 tol,ta0,tam1,tam,ron1,ron,tetm1,tetm
      real*8 sn,sn1,sim1,alf1,betn,omp1,etm,etm1,taf,cm
      integer np,contp(*),nterm,in
      INTEGER NITER,MIT,I,MITER
      real*8 em1,tr,err
      write(99,*) 'iter, residual'
c  Calcolo Ax0
      call aperbn(n,topol,ja,coef1,x,scr1)
c  Calcolo r0=b-Ax0
      do 10 i=1,n
        scr1(i)=tnoti(i)-scr1(i)
10    continue
c  Calcolo L**(-1)r0
cc       call prodl(r,scr1,n,ia,ja,coef2,nterm,topol)
        call yprdly(r,scr1,n,ial,jal,coefl,dl)
c   Calcolo y0=Ux0
cc       call pru(scr1,x,n,ia,ja,coef2,nterm,topol)
       call  aperbn(n,iau,jau,coefu,x,scr1)
       do 15 i=1,n
       x(i)=scr1(i)
       p(i)=r(i)
       u(i)=r(i)
       rc(i)=r(i)
       dm(i)=0.
c      RT(I)=TNOTI(I)
       RT(I)=r(i)
15     continue
c  calcolo v0=Ay1
cc      call prodbh(v,p,n,ia,ja,coef1,coef2,nterm,topol,scr1,scr2)
        call yprdhg(v,p,n,ial,jal,coefl,iau,jau,coefu,ja,topol,
     1              coef1,dl,scr1,scr2)
c---------
      if(np.eq.0)  goto 49
       do 265 i=1,np
         in=contp(i)
         rc(in)=0.0
265    continue
49     continue
c------
      ron1=0.
      tetm1=0.
      etm1=0.
      sn1=0.
      do 20 i=1,n
      sn1=sn1+rc(i)**2
      ron1=ron1+rt(i)*r(i)
20    continue
      sn1=dsqrt(sn1)
      ta0=sn1
      tam1=sn1
      if(ron1.eq.0) then
        write(iout,*)' r0=',ron1,' cambiare r tilde'
         stop
      endif
      niter=0
      miter=0
30    continue
      niter=niter+1
      sim1=0.
      do 35 i=1,n
      sim1=sim1+rt(i)*v(i)
35    continue
      alf1=ron1/sim1
      do 40 i=1,n
      q(i)=u(i)-alf1*v(i)
40    continue
      do 50 i=1,n
      v(i)=u(i)+q(i)
50    continue
c   calcolo B(u_n-1 +q)
cc      call prodbh(au,v,n,ia,ja,coef1,coef2,nterm,topol,scr1,scr2)
        call yprdhg(au,v,n,ial,jal,coefl,iau,jau,coefu,ja,topol,
     1              coef1,dl,scr1,scr2)
      do 55 i=1,n
       rc(i)=rc(i)-alf1*au(i)
55    continue
c  calcolo  ||rn||
c---------
      if(np.eq.0)  goto 59
       do 365 i=1,np
         in=contp(i)
         rc(in)=0.0
365    continue
59     continue
c------
      sn=0.
      do 60 I=1,n
      sn=sn+rc(i)**2
60    continue
      sn=dsqrt(sn)
      mit=0
70    mit=mit+1
      miter=miter+1
      if(mit.eq.1) then
       omp1=dsqrt(sn1*sn)
      else
       omp1=sn
      endif
      tetm=omp1/tam1
      cm=1./dsqrt(1+tetm**2)
      tam=tam1*tetm*cm
      etm=cm**2*alf1
      taf=tetm1**2*etm1/alf1
       if(mit.eq.1) then
        do 100 i=1,n
        dm(i)=u(i)+taf*dm(i)
100     continue
       else
        do 110 i=1,n
        dm(i)=q(I)+taf*dm(i)
110    continue
       endif
       do 120 i=1,n
       x(i)=x(i)+etm*dm(i)
120    continue
c     aggiornamento
       tetm1=tetm
       tam1=tam
       etm1=etm
       em1=miter+1
c      tr=(ta0*tol)/dsqrt(em1)
       err=dsqrt(em1)*tam/ta0
       if((niter.ge.imax).or.(err.le.tol))goto 800
ccccciter             write(6,138)niter,err
138     format(i5,3x,'err=',e12.5)
        if(mit.eq.1) goto 70
c   esecuzione del punto c)
        ron=0.
        do 200 i=1,n
        ron=ron+rt(i)*rc(i)
200     continue
        betn=ron/ron1
        do 210 i=1,n
        u(i)=rc(i)+betn*q(i)
210     continue
        do 220 i=1,n
        p(i)=u(i)+betn*(q(i)+betn*p(i))
220     continue
cc     call prodbh(v,p,n,ia,ja,coef1,coef2,nterm,topol,scr1,scr2)
        call yprdhg(v,p,n,ial,jal,coefl,iau,jau,coefu,ja,topol,
     1              coef1,dl,scr1,scr2)
c   aggiornamento
         sn1=sn
         ron1=ron
         write(99,*) niter,err
         goto 30
800      continue
         do 520 i=1,n
         scr1(i)=x(i)
         x(i)=0.0
520      continue
c  Trovo x=U**(-1)y
cc         call produ(x,scr1,n,ia,ja,coef2,nterm,topol)
           call yprduy(x,scr1,n,iau,jau,coefu)
           return
           end

C                                                                               
      SUBROUTINE YPRDHg(Z,H,N,IAL,JAL,COEFL,IAU,JAU,COEFU,                      
     *JA,TOPOL,COEF1, DL,x,y)                                                   
C     POSTO B=L**(-1)*A*U**(-1), ESEGUE IL PRODOTTO DI B PER UN VETTORE         
C     H E LO METTE IN Z. USA DUE VETTORI DI SERVIZIO X E Y CHE SI               
C     POTERBBERO RIDURRE A UNO SOLO                                             
        implicit none                                                           
      REAL*8 Z(1),H(1),COEFL(1),COEFU(1),DL(1),COEF1(*)                         
      REAL*8 X(*),Y(*)                                                          
      INTEGER N,K,IAL(1),JAL(1),IAU(1),JAU(1),TOPOL(1),JA(1)                    
C     ESEGUO IL PRODOTTO U**(-1)*H=X                                            
C      WRITE(6,150)(IAU(K),K=1,50)                                              
150    FORMAT(' IAU'/ (6I20))                                                   
      CALL YPRDUY(X,H,N,IAU,JAU,COEFU)                                          
C     ESEGUO IL PRODOTTO A*X=Y                                                  
      DO 10 K=1,N                                                               
      Z(K)=0.                                                                   
10    Y(K)=0.                                                                   
C     CALL PRODY(N,M,A,JAA,X,W,Y)                                               
      CALL APERBN(N,TOPOL,JA,COEF1,X,Y)                                         
C     WRITE(6,73)(Y(K),K=1,N)                                                   
73    FORMAT(5E12.5)                                                            
C     ESEGUO IL PRODOTTO L**(-1)*Y=Z                                            
      CALL YPRDLY(Z,Y,N,IAL,JAL,COEFL,DL)                                       
      RETURN                                                                    
      END                                                                       
C     ******************    ROUTINES DI SERVIZIO
C     ESEGUE IL PRODOTTO L**(-1)*B=X RISOLVENDO IL SISTEMA L*X=B
C     CON LE  FORWARD DIFFERENCES
C
      SUBROUTINE YPRDLY(X,B,N,IAL,JAL,COEFL,DL)
        implicit none
      real*8 X(1),B(1),COEFL(1),DL(1),s
      INTEGER IAL(1),JAL(1)
      INTEGER I,IL,NK,K,ILK,N
      X(1)=B(1)/COEFL(1)
      DO 10 I=2,N
      IL=IAL(I)
C     IN=IAL(I+1)-2
C     NK=IN-IL+1
      NK=IAL(I+1)-IAL(I)-1
cc    S=sPDOT(NK,X,JAL(IL),COEFL(IL))
      s=0.d0
      do 15 k=1,nk
      ilk=il+k-1
      s=s+x(jal(ilk))*coefl(ilk)
15    continue
      X(I)=(B(I)-S)/DL(I)
c      write(6,81)i,jal(il),coefl(il),jal(il+1),coefl(il+1)
81     format(2i5,2x,e12.5,1x,i6,e12.5)
c     write(6,85)i,il,nk,s,b(i),dl(i),x(i)
85    format( 3i5,4(1x,e12.5))
10    CONTINUE
      RETURN
      END
C     ESEGUE IL PRODOTTO U**(-1)*X=Y RISOLVENDO IL SISTEMA U*Y=X
C
      SUBROUTINE YPRDUY(Y,X,N,IAU,JAU,COEFU)
        implicit none
      INTEGER I,JN,NK,K,ILK,N
      real*8 Y(1),X(1),COEFU(1),s
      INTEGER IAU(1),JAU(1)
      Y(N)=X(N)/COEFU(IAU(N))
      DO 60 I=N-1,1,-1
C     IN=IAU(I+1)-1
      JN=IAU(I)+1
C     NK=IN-JN+1
      NK=IAU(I+1)-1-IAU(I)
c     S = SPDOT(NK,Y,JAU(JN),COEFU(JN))
      s=0.d0
      do 15 k=1,nk
      ilk=jn+k-1
      s=s+y(jau(ilk))*coefu(ilk)
15    continue
      Y(I)=X(I)-S
60    CONTINUE
      RETURN
      END
