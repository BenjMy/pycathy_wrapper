C
C**************************  STORMB ************************************
C
C  calculate volume of change in storage
C  between the current time level and the previous time level.
C  DSTORE > 0 for net increase in storage.
C
C***********************************************************************
C
      SUBROUTINE STORMB(PNEW,PTIMEP,SWNEW,SWTIMEP,N,DSTORE,VOLNOD,SNODI,
     1                  PNODI)
C
      IMPLICIT NONE
      INTEGER  K
      INTEGER  N
      REAL*8   DSTORE
      REAL*8   PNEW(*),PTIMEP(*),VOLNOD(*),SNODI(*),PNODI(*)
      REAL*8   SWNEW(*),SWTIMEP(*)
C
      DSTORE=0.0D0
      DO K=1,N
CM       DSTORE=DSTORE + VOLNOD(K)*(ETAI(K) + ETAIP(K))*0.5D0*
CM   1                             (PNEW(K) - PTIMEP(K))
         DSTORE=DSTORE + VOLNOD(K)*(SNODI(K)*(SWNEW(K)+SWTIMEP(K))
     1       *0.5*(PNEW(K)-PTIMEP(K))+PNODI(K)*(SWNEW(K)-SWTIMEP(K)))
      END DO
C
      RETURN
      END
