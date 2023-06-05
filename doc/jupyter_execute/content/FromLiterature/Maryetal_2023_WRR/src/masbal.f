C
C**************************  MASBAL ************************************
C
C  Calculate mass balance for the run
C
C***********************************************************************
C
      SUBROUTINE MASBAL(IOPT,N,NP,NQ,NNOD,NSF,IPEAT,IVGHU,NUDC,
     1                  NLKP,NT,NTRI,TETAF,TIME,DELTAT,
     2                  TETRA,TP,TOPOL,JA,CONTP,IFATM,
     3                  SFEX,SFFLAG,NSFNUM,NSFNOD,
     4                  PNEW,PTIMEP,PDIFF,SNODI,PNODI,PORE,INDE,
     5                  SENODI,ETAI,SEELT,ETAE,Z,VOLNOD,QPNEW,QPOLD,Q,
     6                  ATMACT,ATMOLD,COEF1,XT5,SCR,SWNEW,SWTIMEP,
     7                  BKCFLOW_NODE,SFQ,SFQP,NUDNOD,NUDCUM)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  IOPT,N,NP,NQ,NNOD,NSF,IPEAT,IVGHU,NUDC
      INTEGER  NLKP,NT,NTRI
      INTEGER  TETRA(5,*),TP(*),TOPOL(*),JA(*)
      INTEGER  CONTP(*),IFATM(*),SFEX(NSFMAX,*),SFFLAG(*)
      INTEGER  NSFNUM(*),NSFNOD(NSFMAX,*)
      REAL*8   DELTAM
      REAL*8   TETAF,TIME,DELTAT
      REAL*8   PNEW(*),PTIMEP(*),PDIFF(*),SWNEW(*),SWTIMEP(*)
      REAL*8   SNODI(*),PNODI(*),PORE(*),INDE(*)
      REAL*8   SENODI(*),ETAI(*),SEELT(*),ETAE(*),Z(*),VOLNOD(*)
      REAL*8   QPNEW(*),QPOLD(*),Q(*),BKCFLOW_NODE(*)
      REAL*8   ATMACT(*),ATMOLD(*),COEF1(*),XT5(*),SCR(*)
      REAL*8   SFQ(NSFMAX,*),SFQP(NSFMAX,*),NUDNOD(*),NUDCUM(*)
      INCLUDE 'MB_HGRAPH.H'
C
C  back-calculate fluxes at all Dirichlet nodes
C
      IF (IOPT .EQ. 1) THEN
         CALL BKPIC(N,NP,NNOD,TETAF,TOPOL,JA,CONTP,IFATM,PDIFF,
     1              QPNEW,QPOLD,ATMACT,ATMOLD,COEF1,XT5,SCR,
     2              NSF,NSFNUM,NSFNOD,SFQ,SFQP,SFEX)
      ELSE
         CALL BKNEW(NP,NNOD,TETAF,TOPOL,JA,CONTP,IFATM,PDIFF,
     1              QPNEW,QPOLD,ATMACT,ATMOLD,COEF1,XT5,
     2              NSF,NSFNUM,NSFNOD,SFQ,SFQP,SFEX)
      END IF
C  CARLOTTA I am saving the back calculated nodal flow BKCFLOW_NODE
      CALL VCOPYR(N,BKCFLOW_NODE,SCR)
C  
C  calculate general storage term needed for mass balance calculations.
C  We need the storage term at the previous and current time levels. 
C  Storage values for the previous time level are computed and stored
C  in SCR, which is then passed to subroutine STORMB.
C
      IF  (IPEAT .EQ. 0) THEN
          CALL CHMASS(NLKP,N,NT,NTRI,IVGHU,TETRA,TP,
     1                PTIMEP,SNODI,PNODI,SCR,ETAE)
          CALL CHMASS(NLKP,N,NT,NTRI,IVGHU,TETRA,TP,
     1                PNEW,SNODI,PNODI,ETAI,ETAE)
      ELSE
          CALL CHMASPT(NLKP,N,NNOD,NT,NTRI,IVGHU,TETRA,TP,
     1                 PTIMEP,SNODI,PNODI,PORE,INDE,Z,
     2                 SENODI,SCR,SEELT,ETAE)
          CALL CHMASPT(NLKP,N,NNOD,NT,NTRI,IVGHU,TETRA,TP,
     1                 PNEW,SNODI,PNODI,PORE,INDE,Z,
     2                 SENODI,ETAI,SEELT,ETAE)
      END IF
C  
C  calculate mass balance errors 
C  
      IF (DELTAT .LT. 1.0D+15) THEN
         DELTAM=0.5D0*DELTAT 
CM       CALL STORMB(PNEW,PTIMEP,N,DSTORE,VOLNOD,SCR,ETAI)
         CALL STORMB(PNEW,PTIMEP,SWNEW,SWTIMEP,N,DSTORE,VOLNOD,SNODI,
     1               PNODI)
      ELSE
         DELTAM=1.0D0
         DSTORE=0.0D0  
      END IF
      CALL FLUXMB(N,NP,NQ,NNOD,NSF,NUDC,QPNEW,Q,IFATM,ATMACT,
     1            NDIN,NDOUT,NNIN,NNOUT,ADIN,ADOUT,ANIN,ANOUT, 
     2            NUDIN,NUDOUT,NUDNOD,NUDCUM,
     3            NSFNUM,NSFNOD,SFEX,SFQ,SFFLW,SFFLAG,
     4            TIME,DELTAT)
      VADIN =(ADIN  + ADINP) *DELTAM
      VADOUT=(ADOUT + ADOUTP)*DELTAM
      VNDIN =(NDIN  + NDINP) *DELTAM
      VNDOUT=(NDOUT + NDOUTP)*DELTAM
      VANIN =(ANIN  + ANINP) *DELTAM
      VANOUT=(ANOUT + ANOUTP)*DELTAM
      VNNIN =(NNIN  + NNINP) *DELTAM
      VNNOUT=(NNOUT + NNOUTP)*DELTAM
      VSFFLW=(SFFLW + SFFLWP)*DELTAM
      VNUDIN =(NUDIN  + NUDINP) * DELTAM
      VNUDOUT =(NUDOUT  + NUDOUTP) * DELTAM
      VIN =VADIN  + VNDIN  + VANIN  + VNNIN + VNUDIN
      VOUT=VADOUT + VNDOUT + VANOUT + VNNOUT + VSFFLW + VNUDOUT
      ERRAS=VIN + VOUT - DSTORE
      IF (DELTAT .GE. 1.0D+15) THEN
         IF (VIN .NE. 0.0D0) THEN
            ERREL=100.0D0*ERRAS/VIN
         ELSE
            ERREL=0.0D0
         END IF
      ELSE 
         IF ((VIN + VOUT) .NE. 0.0D0) THEN
            ERREL=100.0D0*ERRAS/(VIN + VOUT)
         ELSE
            ERREL=0.0D0
         END IF
      END IF 
C
      RETURN
      END
