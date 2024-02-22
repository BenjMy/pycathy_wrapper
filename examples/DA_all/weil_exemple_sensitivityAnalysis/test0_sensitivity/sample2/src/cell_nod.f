C
C**************************  CELL_NOD **********************************
C
C  transfer cell-values to nodal values taking into account
C  possible coarsening of the triangulation (dostep>1)
C
C***********************************************************************
C
      subroutine cell_nod(ncell,nrow,ncol,nnod,ntri,dostep,
     1                    tp2d,triang,indcel,dem_map,
     2                    pondcell,cellcoarse,pondnod)

      implicit none
      include 'CATHY.H'
      
      integer ncell,nrow,ncol,nnod,ntri,dostep
      integer icell,i,j,k,l,ii,inod,nonnull,m,n
      integer tp2d(*),triang(4,*)
      integer indcel(rowmax,*)
      real*8  pondnod(*),cellcoarse(maxcel),pondcell(maxcel)
      real*8  dem_map(rowmax,*)
      do i=1,ncell
         cellcoarse(i) = 0.0d0
      end do
c 
c dalle celle piccole alle grandi
c
      icell = 0
      do i=1,nrow,dostep
         do j=1,ncol,dostep
c       write(6,*) 'i,j=',i,j,'dem=',dem_map(i,j),'indcel=',indcel(i,j)
          nonnull=0
          do m=1,dostep
           do n=1,dostep
            if(dem_map(i+m-1,j+n-1).ne.0 .and. i+dostep-1.le.nrow
     &                      .and. j+dostep-1.le.ncol) then
               nonnull=nonnull+1
            end if
           end do
          end do
cccd            write(6,*) 'i,j=',i,j,'nonnull=',nonnull
            if(nonnull.gt.0) then  
               icell = icell +1
               do k=1,dostep
                  do l=1,dostep
                    if(dem_map(i+k-1,j+l-1) .ne. 0 .and.
     1                 indcel(i+k-1,j+l-1).ne.0)  then
                     cellcoarse(icell) = cellcoarse(icell) + 
     1                                   pondcell(indcel(i+k-1,j+l-1))

                    end if
                  end do 
               end do
               cellcoarse(icell) = cellcoarse(icell)/nonnull
cccd            write(6,*) 'icell=',icell,'cellcoarse=',cellcoarse(icell)
            end if
         end do
      end do
ccc      write(6,*) 'icell=',icell
c
c dalle celle grandi ai nodi
c
      do k=1,nnod
        tp2d(k)=0
      end do
  
      do k=1,nnod
        pondnod(k)=0.0d0
      end do

      DO K=1,ntri
         icell = k/2 + mod(k,2)
         DO II=1,3
            INOD=TRIANG(II,K)
            pondnod(INOD)=pondnod(INOD) + cellcoarse(icell)
            TP2d(INOD)=TP2d(INOD)+1
         END DO
      END DO
      DO K=1,NNOD
         pondnod(K) = pondnod(K)/TP2d(K)
      END DO
C
      RETURN
      END
