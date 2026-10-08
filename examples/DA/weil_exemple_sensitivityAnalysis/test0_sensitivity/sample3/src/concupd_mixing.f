C
C**************************  CONCUPD_MIXING ****************************
C
C Update surface CONCENTRATION FOR MIXING PONDING WATER AND FIRST 
C LAYER -> ONLY FOR PONDED NODES (by LG)
C
C***********************************************************************
C
      SUBROUTINE CONCUPD_MIXING(NNOD,N,NT,TETRA,PONDNOD,VOLNOD,PNODI,SW,
     1     ARENOD,CONCNODUPD,CNNEW,CNEW,MIXPART,TP,PEL,VOLU,SWNEW,
     2     DELTAT,SOURCE_MIXING,ATMACT,TRAFLNOD_FLOW,OVFLNOD,PONDNODP)
C
      IMPLICIT NONE
      INTEGER  NNOD,I,NT,TETRA(5,NT),ITRIA,J,JJ,tp(*),J1,J2,J3,J4,N
      REAL*8   PONDNOD(*),ARENOD(*),VOLNOD(*),PNODI(*),SW(*)
      REAL*8   PEL(*),SWNEW(*),VOLU(*),TIME,DELTAT
      REAL*8   CONCNODUPD(*),CNNEW(*),MIXPART,CNEW(*),CNNEW_stock(nnod)
      REAL*8   VOLSURF(NNOD),VOLSUB(NNOD),MASSSURF(NNOD),MASSSUB(NNOD)
      REAL*8   mass_node(NT),SOURCE_MIXING(nnod),CONCNODUPD_new(nnod)
      real*8   test,futurepond(nnod),TRAFLNOD_FLOW(*),OVFLNOD(*)
      REAL*8   ATMACT(*),PONDNODP(*)
C    
      CALL INIT0R(NNOD,CONCNODUPD_new)
      CALL VCOPYR(NNOD,CNNEW_stock,cnnew)
      
Conditions : ponding > 0 au temps t et au temps t+1    
      
      test=0.0d0
      DO I = 1,NNOD
        IF (PONDNOD(i).GT.0.0d0) THEN
C
        VOLSURF(I)=PONDNOD(I)*ARENOD(I)
        VOLSUB(I)=SW(I)*PNODI(I)*VOLNOD(I)
        MASSSURF(I)=VOLSURF(I)*CONCNODUPD(I)
        MASSSUB(I)=VOLSUB(I)*CNNEW(I)
C
        CONCNODUPD_new(I)=(MASSSURF(I)+MASSSUB(I))/
     1                   (MIXPART*VOLSURF(I)+VOLSUB(I))

        source_mixing(i)= ( CONCNODUPD_new(i) - CONCNODUPD(i) )
     1        * PONDNOD(I)*ARENOD(I) !/ deltat
c     
        CNNEW(I)=MIXPART*CONCNODUPD_new(I)     
        CONCNODUPD(I)=CONCNODUPD_new(I)
        END IF
      END DO  
      
500   CONTINUE
       
C     
      RETURN  
      END
