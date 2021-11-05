C
C**************************  ERT_MEASURE *******************************
C
C     Reads pressure values, calculates saturation at the assimilation
C     time and calculates the electrical conductivity of the soil by the
C     Archie law = ELECTRC_NODI (for the 2D cross section only).
C     ELETRC_NODI is calculated for the expanded 2D grid -> input for
C     SAT2D
C 
C***********************************************************************
C
      SUBROUTINE ERT_MEASURE(NLKP,N,NT,IVGHU,TETRA,TP,
     1                       ENPNEW,ENPNODI,ENSNODI,ARCHIE,EN_ERT,
     2                       ENPT,NERT,NTRI,NENS,ELECTRC_NODI,
     3                       PNEW,SNODI,PNODI,SW,CKRW,SWE,CKRWE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  NRE,I,J,K,CONT,CONT1
      INTEGER  NLKP,N,NT,IVGHU,NERT,NTRI,NENS
      INTEGER  TETRA(5,*),TP(*),ENPT(*)
      REAL*8   A_AR,CW_AR,M_AR,N_AR
      REAL*8   ENPNEW(NMAX,*),ENPNODI(NMAX,*),ENSNODI(NMAX,*)
      REAL*8   ARCHIE(4),EN_ERT(MAXNUDN,*),ELECTRC_NODI(*)
      REAL*8   PNEW(*),SNODI(*),PNODI(*),SW(*),CKRW(*)
      REAL*8   SWE(*),CKRWE(*)
C
      A_AR=ARCHIE(1)
      CW_AR=ARCHIE(2)
      M_AR=ARCHIE(3)
      N_AR=ARCHIE(4)
      DO NRE=1,NENS
         OPEN (19,file='output/el_conductivity',status='unknown',
     1          form='unformatted')
         J=ENPT(NRE)
         DO I=1,N
            PNEW(I)=ENPNEW(I,J)
            SNODI(I)=ENSNODI(I,J)
            PNODI(I)=ENPNODI(I,J)
         END DO
         CALL CHPARM(N,SNODI,PNODI,IVGHU)
         CALL CHVELO(NLKP,N,NT,NTRI,IVGHU,TETRA,TP,
     1               PNEW,SNODI,PNODI,SW,CKRW,SWE,CKRWE)
         I=1
         CONT=0
         CONT1=0
         DO K=1,N
c              CONT1=CONT1+1
c              CONT=CONT+1
c              IF (CONT1 .EQ.20) I=I+2
c              IF (CONT.EQ.204) I=I+10
c              ELECTRC_NODI(K)=I
c              IF (CONT1 .EQ.51) THEN
c                 I  =I-2
c                 CONT1=0
c              END IF
c              IF (CONT .EQ.612) THEN
c                 I  =I-10
c                 CONT=0
c              END IF
             ELECTRC_NODI(K)=
     1          1/A_AR*CW_AR*PNODI(K)**M_AR*SW(K)**N_AR
             WRITE(19)ELECTRC_NODI(k)
         END DO
         CLOSE(19)
         call system('./prot.sh')
         open(18,file='risp/ert.risp',status='unknown')
         DO I=1,NERT
             read(18,*) EN_ERT(I,J)
         END DO
         close(18)
      END DO 
      do i=1,nert      
             write(777,*) i,(en_ert(i,j),j=1,nens)
      end do
C
      RETURN
      END
