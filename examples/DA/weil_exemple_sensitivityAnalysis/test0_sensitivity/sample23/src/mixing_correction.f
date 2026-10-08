C
C************************** mixing_correction **************************
C
C Correction after the application of the mixing module in case of no
C respect of mass conservation. 
C -> MIX_CORRECT(MAXCEL) (filled out in SURF_FLOWTRA as debt/suplus mass) 
C is applied to correct first layer subsurface concentration 
C LG
C
C***********************************************************************
C

      SUBROUTINE mixing_correction(NROW,NCOL,INDEX,INDEX_WITH_LAKES,
     1                     N,NT,NNOD,TETRA,nstr,NCELL,CELL,MIX_CORRECT,
     2                     CNEW,CNNEW,PEL,PNODI,
     3                     VOLU,VOLNOD,SWNEW,SW,TP,COLD,CNOLD)

      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER   NROW,NCOL,IROW,ICOL,I,J,JR,K,ii,inod
      INTEGER   NCELL,cell(5,*),NNOD
      INTEGER   INDEX(ROWMAX,*),INDEX_WITH_LAKES(ROWMAX,*)
      INTEGER   I_BASIN_SURF,I_BASIN_SUB
      INTEGER   N,NT,NSTR
      INTEGER   TETRA(5,*),TP(*)
      REAL*8    CNEW(*),CNNEW(*),volu(*),volnod(*),pel(*),pnodi(*)
      real*8    sw(*),swnew(*),cold(*),CNOLD(*)
cm    real*8    MIX_CORRECT_node(nnod),MIX_CORRECT(*)
cm    real*8    CNNEW_NEW(NNOD),MIX_CORRECT_sub(ncell)
      real*8    MIX_CORRECT_node(NODMAX),MIX_CORRECT(*)
      real*8    CNNEW_NEW(NODMAX),MIX_CORRECT_sub(MAXCEL)


      CALL INIT0R(NNOD,CNNEW_new)
      CALL INIT0R(NCELL,MIX_CORRECT_sub)
      CALL INIT0R(NNOD,MIX_CORRECT_node)
C
      
C STEP 1 : from index_cell to cell
      DO ICOL=1,NCOL
         DO IROW=1,NROW
            I_BASIN_SUB=NCOL*(NROW-IROW)+ICOL
            I_BASIN_SURF=(ICOL-1)*NROW + IROW

            JR=MOD(I_BASIN_SUB,NCOL)
            IF (JR.NE.0) THEN
               J=JR
               I=(I_BASIN_SUB-J)/NCOL+1
            ELSE
               J=NCOL
               I=I_BASIN_SUB/NCOL
            END IF
            
            IF((INDEX(I,J).NE.0).AND.
     &           (INDEX_WITH_LAKES(I,J).NE.0))THEN
               
               MIX_CORRECT_sub(I_BASIN_SUB) = MIX_CORRECT(I_BASIN_SURF)
            END IF
         END DO
      END DO


C STEP 2 : from cell to node

      DO K=1,ncell
         DO II=1,4
            INOD=CELL(II,K)
            MIX_CORRECT_node(INOD)= MIX_CORRECT_node(INOD) + 
     1              MIX_CORRECT_sub(k)/4
         END DO
      END DO
C
C STEP 3 : CORRECTION SUBSURFACE ON THE NODES (CNNEW)
C
        DO I=1,NNOD
C       
        CNNEW_new(I)=(CNNEW(I)*VOLNOD(I)*PNODI(I)*SW(I)+
     1             MIX_CORRECT_node(i)) / (VOLNOD(I)*PNODI(I)*SW(I))
        END DO

C STEP 4 : CNNEW -> CNEW

      CALL nodetotetra_tra(N,NT,TETRA,nstr,CNEW,CNNEW_NEW,PEL,PNODI,
     1           VOLU,VOLNOD,SWNEW,SW,TP,NNOD,COLD,CNNEW)

      RETURN
      END
      
