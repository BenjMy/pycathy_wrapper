!-----------------------------------------------------------------------------
!
! IJ2QFLAG
!
! The program IJ2QFLAG reads the coordinates in cellsID file... 
!
!
!-----------------------------------------------------------------------------

program ij2qflag

use mpar
use mbbio

implicit none

integer(kind=4) i,j,i_basin,ic,cont,iselect,ht
integer(kind=4) xllcorner,yllcorner,nodata
integer(kind=4) north,south,east,west
integer(kind=4),allocatable :: q_flag(:,:)

character(len=24) :: fstr
character(len=24) :: nstr

call rparfile('hap.in')

call load_dtm('basin_b','basin_i')

allocate(q_flag(N,M))

open (15,file='dtm_q_output',access='sequential',status='unknown')
open (16,file='cellsID',access='sequential',status='unknown',iostat=ic)

write(6,*)
write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,'(1x,a)') 'This program creates/overwrites the dtm_q_output file!'
write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*)
write(6,'(1x,a)') 'Select:'
write(6,'(1x,a)') '0) to exit'
write(6,'(1x,a)') '1) to continue'
write(6,'(1x,a)') '(Ctrl C to exit)'
write(6,*)
write(6,'(1x,a,$)') '-> '
read(5,*) iselect
if (iselect.eq.0) stop

write(6,*)
write(6,'(1x,a)') 'Select the header type:'
write(6,'(1x,a)') '0) None'
write(6,'(1x,a)') '1) ESRI ascii file'
write(6,'(1x,a)') '2) GRASS ascii file'
write(6,'(1x,a)') '(Ctrl C to exit)'
write(6,*)
write(6,'(1x,a,$)') '-> '
read(5,*) ht

! select the nodata threshold

write(6,*)
write(6,'(1x,a)') 'Select the nodata value:'
write(6,'(1x,a)') '(Ctrl C to exit)'
write(6,*)
write(6,'(1x,a,$)') '-> '
read(5,*,iostat=ic) nodata
do while (ic /= 0)
   write(6,'(1x,a)') 'Invalid selection: try again...'
   write(6,'(1x,a)') '(Ctrl C to exit)'
   write(6,*)
   write(6,'(1x,a,$)') '-> '
   read(5,*,iostat=ic) nodata
end do

do i=1,N
   do j=1,M
      i_basin=(i-1)*M+j
      if (dtm_index_pr(i_basin) /= 0) then
         q_flag(i,j)=0
      else
         q_flag(i,j)=nodata
      end if
   end do
end do

ic=0

cont=0
do while (ic == 0)
    read(16,*,iostat=ic) i,j
    i_basin=(i-1)*M+j
    if (dtm_index_pr(i_basin) /= 0) then
       q_flag(i,j)=1
       cont=cont+1
    else
       write(6,*)
       write(6,'(1x,a19,i4,i4,a24)') &
       'invalid coordinates',i,j,': cell out of catchment!'
    end if
end do

! ESRI ascii file header writing
if (ht.eq.1) then
   write(15,"(a,8x,i5)") 'ncols',N
   write(15,"(a,9x,i5)") 'nrow',M
   write(15,"(a,4x,i5)") 'xllcorner',1000
   write(15,"(a,4x,i5)") 'yllcorner',1000
   write(15,"(a,7x,f6.2)") 'cellsize',delta_x
   write(15,"(a,1x,i5)") 'NODATA_value',-9999
end if

! GRASS ascii file header writing
if (ht.eq.2) then
   write(15,'(a,1x,i5)') 'north:',0
   write(15,'(a,1x,i5)') 'south:',0
   write(15,'(a,2x,i5)') 'east:',0
   write(15,'(a,2x,i5)') 'west:',0
   write(15,'(a,2x,i5)') 'rows:',M
   write(15,'(a,2x,i5)') 'cols:',N
end if

write(fstr,'(i10)') N 
fstr='('//trim(adjustl(fstr))//'i6'//')'

do j=M,1,-1
   write(15,fstr) (q_flag(i,j),i=1,N)
end do

write(6,*)
write(6,'(1x,i5,1x,a11)') cont-1,'cells found'
write(6,*)

end program ij2qflag
