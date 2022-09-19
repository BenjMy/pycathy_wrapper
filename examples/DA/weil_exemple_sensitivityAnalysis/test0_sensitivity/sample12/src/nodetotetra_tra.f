C
C*************************** nodetotetra_tra ***************************
C
C  FOR TRANSPORT: pass the values from the nodes to the centroid of each
C  tretrahedra.
C  This subroutine juste changes the fisrt layer of tetrahedral concentration
C  to avoid numerical dispersion.
C
C***********************************************************************
C
      SUBROUTINE nodetotetra_tra(N,NT,TETRA,nstr,CNEW,CNNEW,PEL,PNODI,
     1           VOLU,VOLNOD,SWNEW,SW,TP,NNOD,COLD,CNOLD)
C
      IMPLICIT  NONE
      INTEGER   i,K,II,INOD,nnod
      INTEGER   N,NT,NSTR
      INTEGER   TETRA(5,*),TP(*),cont(nnod)
      REAL*8    CNEW(*),CNNEW(*),volu(*),volnod(*),pel(*),pnodi(*)
      real*8    sw(*),swnew(*),cold(*),CNOLD(*)
      REAL*8    volneigh(N)
C
        CALL INIT0R(N,volneigh)
c      
      DO K=1,NT
         DO II=1,4
            INOD=TETRA(II,K)
            volneigh(INOD)=volneigh(INOD)+(VOLU(k)*pel(k)*SWNEW(k))
         END DO
      END DO

            do i=1,nnod
            cont(i)=0
            end do

      do i = 1, nt
        DO II=1,4 
        inod=tetra(ii,i)
        if ((inod.LE.NNOD).AND.(cnnew(inod)-cnold(inod).NE.0.0d0)) THEN
        cnew(i) = cnew(i) + volnod(inod)*pnodi(inod)*
     1   sw(inod)/volneigh(inod) * (cnnew(inod)-cnold(inod))
            cont(inod)=cont(inod)+1
        end if
        end do
      end do
      
      do i=1, nnod
      if (cont(i).ne.0) THEN
c      write(*,*) i, cont(i)
      end if
      end do
C
      RETURN
      END
