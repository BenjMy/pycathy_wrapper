C
C**************************  MASS_BALANCE ******************************
C
C  Compute the total solute mass in the domain at actual time
C
C***********************************************************************
C
c      SUBROUTINE MASS_BALANCE(N,NT,TETRA,TIME,DELTAT,CNEW_EL,COLD_EL,
c     1                        COLDOLD,CNEW_NODE,VOLU,VOLNOD,SW_EL,
c     2                        SWP_EL,SW_NODE,PEL,PNODI,MASS_TOT,
c     3                        MASS_BILAN,MASS_OUT_SF,NSF,NSFNOD,SFQ,
c     4                        NSFNUM,MASS_OUT_DIR,ANP,ACONTP,QPNEW)
      SUBROUTINE MASS_BALANCE(N,NT,TIME,DELTAT,CNEW_EL,COLD_EL,
     1                        COLDOLD,VOLU,VOLNOD,SW_EL,
     2                        SWP_EL)
C
      IMPLICIT  NONE
c     INCLUDE   'CATHY.H'
c     INCLUDE  'IOUNITS.H'
      INTEGER   I,J
      INTEGER   NT,N
c      INTEGER   TETRA(5,*),NSF,NSFNOD(NSFMAX,*),NSFNUM(*)
c      INTEGER   ANP,ACONTP(*)
      REAL*8    TIME,DELTAT
      REAL*8    CNEW_EL(*),COLD_EL(*),COLDOLD(*)
      REAL*8    VOLU(*),VOLNOD(*),SW_EL(*),SWP_EL(*)
cxcx  REAL*8    PEL(NT)
c      REAL*8    SFQ(NSFMAX,*)
      REAL*8    MASS_TOT,MASS_BILAN,MASS_OUT_SF,MASS_OUT_DIR
c      REAL*8    QPNEW(*)
      REAL*8    MASS_TOT_unsat,MASS_TOT_sat
C
      MASS_TOT=0.0d0
      MASS_BILAN=0.0d0
      MASS_TOT_sat=0.0d0
      MASS_TOT_unsat=0.0d0
      DO i=1,NT
cxcx  PEL(I)=0.5
cxcx  MASS_TOT=MASS_TOT+(CNEW_EL(I)*VOLU(I)*SW_EL(I)*PEL(I))
      MASS_TOT=MASS_TOT+(CNEW_EL(I)*VOLU(I)*SW_EL(I)*0.5)
cxcx  MASS_BILAN=MASS_BILAN+
cxcx 1           (CNEW_EL(I)*VOLU(I)*SW_EL(I)*PEL(I)-
cxcx 2           COLDOLD(I)*VOLU(I)*SWP_EL(I)*PEL(I))
      MASS_BILAN=MASS_BILAN+
     1           (CNEW_EL(I)*VOLU(I)*SW_EL(I)*0.5-
     2           COLDOLD(I)*VOLU(I)*SWP_EL(I)*0.5)
        IF (SW_EL(I).EQ.1) THEN
cxcx    MASS_TOT_sat=MASS_TOT_sat+(CNEW_EL(I)*VOLU(I)*PEL(I))
        MASS_TOT_sat=MASS_TOT_sat+(CNEW_EL(I)*VOLU(I)*0.5)
        ELSE
cxcx    MASS_TOT_unsat=MASS_TOT_unsat
cxcx 1             +(CNEW_EL(I)*VOLU(I)*SW_EL(I)*PEL(I))
        MASS_TOT_unsat=MASS_TOT_unsat
     1             +(CNEW_EL(I)*VOLU(I)*SW_EL(I)*0.5)
        END IF
      END DO
      write(*,*) 'MASS_TOT',MASS_TOT
      write(*,*) 'MASS_BILAN',MASS_BILAN
      write(*,*) 'MASS_TOT_sat',MASS_TOT_sat
      write(*,*) 'MASS_TOT_unsat',MASS_TOT_unsat
C
C MASS BALANCE FROM SEE PAGE FACES
c      IF (NSF.NE.0) THEN
c         CALL ELTNOD_TRA(N,NT,TETRA,CNEW_EL,CNEW_NODE,PEL,PNODI,
c     1              VOLU,VOLNOD,SW_NODE,SW_EL)
C
c         MASS_OUT_SF=0.0D0
c         DO I=1,NSF
c           DO j=1,NSFNUM(I)
c           MASS_OUT_SF=MASS_OUT_SF+(SFQ(I,J)*CNEW_NODE(NSFNOD(I,J))*
c     1                       DELTAT)
c           END DO
c         END DO
c      END IF
C
C MASS BALANCE FROM DIRICHLET NODES
C
c      IF (ANP.NE.0) THEN
c      CALL ELTNOD_TRA(N,NT,TETRA,CNEW_EL,CNEW_NODE,PEL,PNODI,
c     1              VOLU,VOLNOD,SW_NODE,SW_EL)
C
c         MASS_OUT_DIR=0.0D0
c         DO I=1,ANP
c            IF (QPNEW(I) .LT. 0.0D0) THEN
c            MASS_OUT_DIR=MASS_OUT_DIR+(QPNEW(I)
c     1                      *CNEW_NODE(ACONTP(I))*DELTAT)
c            END IF
c        END DO
c      END IF
C
C OUTPUT WRITING
C            write(IOUT55,*) TIME,MASS_TOT,MASS_OUT_SF,MASS_OUT_DIR
C            write(*,*) 'TIME            SUBSURFACE SOLUTE MASS',
C     1     '           SF OUTPUT MASS        DIR OUTPUT MASS'
C            write(*,*) TIME,MASS_TOT,MASS_OUT_SF,MASS_OUT_DIR
C
C
      RETURN
 
        END
