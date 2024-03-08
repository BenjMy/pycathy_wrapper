C
C**************************  CELL_NOD_TRA ******************************
C
C  transfer cell-values to nodal values taking into account
C  possible coarsening of the triangulation (dostep>1)
C
C***********************************************************************
C
      subroutine cell_nod_tra(ntri,nnod,tp2d,triang,conccel,concnod,
     1                        delta_x,delta_y,arenod,pondcel,pondnod)
C
      implicit none
      integer ntri,nnod
      integer k,ii,inod,i,icell
      integer tp2d(*),triang(4,*)
      real*8  concnod(*),conccel(*)
      real*8  mass_node(nnod)
      real*8  delta_x,delta_y,arenod(*)
      real*8  pondcel(*),pondnod(*)
C
      call init0r(nnod,mass_node)
      call init0r(nnod,concnod)
      do k=1,nnod
        tp2d(k)=0
      end do

      DO K=1,ntri
         icell = k/2 + mod(k,2)
         DO II=1,3
            INOD=TRIANG(II,K)
            mass_node(INOD) = mass_node(INOD) + conccel(icell) *
     1                        delta_x*delta_y*pondcel(icell)
            TP2d(INOD)=TP2d(INOD)+1
         END DO
      END DO
      DO i=1,NNOD
         if (pondnod(i).NE.0.0d0) then
            concnod(i) = mass_node(i)/(TP2d(i)*arenod(i)*pondnod(i))
         end if
      END DO
C
      RETURN
      END
