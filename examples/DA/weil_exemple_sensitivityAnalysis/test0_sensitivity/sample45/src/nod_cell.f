C
C******************************NOD_CELL*********************************
C
C   transfer nodal values to dem cell values taking into account
C   coarsening of triangulation (dostep>1)
C
C***********************************************************************
C
      subroutine nod_cell(ncell,nrow,ncol,dostep,ncell_coarse,
     1                   nnod,cell,dem_map,indcelwl,cellcoarse,nodvalue,
     2                   cellvalue,arenod,delta_x,delta_y)

      implicit none
      include 'CATHY.H'
      include 'IOUNITS.H'

      integer nrow,ncol
      integer icell,k,l,ncell,nnod,i,j,icell_coarse,m,n
      integer dostep,ncell_coarse,nonnull
      integer cell(5,*),indcelwl(rowmax,*)
      real*8 delta_x,delta_y
      real*8 cellvalue(maxcel),nodvalue(*),arenod(*)
      real*8 cellcoarse(maxcel)
      REAL*8 dem_map(rowmax,colmax) 
 
      do k=1,nnod
         nodvalue(k)=nodvalue(k)/arenod(k)
      end do
      
c
c dai nodi alle celle grandi
c
      do icell_coarse=1,ncell_coarse
         cellcoarse(icell_coarse)=0.0d0
         do k=1,4
            cellcoarse(icell_coarse)=cellcoarse(icell_coarse)+
     1                        nodvalue(cell(k,icell_coarse))
         end do
         cellcoarse(icell_coarse)=cellcoarse(icell_coarse)*0.25d0
      end do

ccccd      write(99,*) 'i   cellcoarse(i)'
ccccd      do icell_coarse=1,ncell_coarse
ccccd         write(99,*) icell_coarse,cellcoarse(icell_coarse)
cccccd      end do
      
c
c dalle celle grandi alle piccole
c      
      icell_coarse = 0
      do i=1,nrow,dostep
         do j=1,ncol,dostep
          nonnull=0
          do m=1,dostep
           do n=1,dostep 
            if(dem_map(i+m-1,j+n-1).ne.0 .and. i+dostep-1.le.nrow
     &                         .and. j+dostep-1.le.ncol) then
               nonnull=nonnull+1
            end if
           end do
          end do
ccccd          write(99,*) 'i=',i,'j=',j,'nonnull=',nonnull
            if (nonnull.gt.0) then   
               icell_coarse = icell_coarse + 1
               do k=1,dostep
                 do l=1,dostep 
                  if (dem_map(i+k-1,j+l-1).ne.0) then
                   icell = indcelwl(i+k-1,j+l-1)
cccccd                   write(99,*) 'cella non nulla=',icell
                   cellvalue(icell) = cellcoarse(icell_coarse)*
     1                                delta_x*delta_y*dostep*dostep/
     2                                nonnull
ccccc                   write(99,*)'i=',icell,'cellvalue=',cellvalue(icell)
                  end if
                 end do
               end do
            end if
         end do
      end do
 
      return

      end
