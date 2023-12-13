!-----------------------------------------------------------------------------
!
! QFLAG2IJ
!
! The program QFLAG2IJ reads the flag in dtm_q_output file and store its 
! coordinates in cellsID file.
!
!-----------------------------------------------------------------------------

program qflag2ij

use mpar

implicit none

integer(kind=4) i,j,i_basin,cont,ht
integer(kind=4),allocatable :: q_flag(:,:)

call rparfile('hap.in')

allocate(q_flag(N,M))

!call load_dtm("basin_b","basin_i")

open (15,file='dtm_q_output',access='sequential',status='unknown')
open (16,file='cellsID',access='sequential',status='unknown')

write(6,*)
write(6,'(1x,a)') 'Select the header type:'
write(6,'(1x,a)') '0) None'
write(6,'(1x,a)') '1) ESRI ascii file'
write(6,'(1x,a)') '2) GRASS ascii file'
write(6,'(1x,a)') '(Ctrl C to exit)'
write(6,*)
write(6,'(1x,a,$)') '-> '
read(5,*) ht

if (ht /= 0) then
   do i=1,6
      read(15,*)
   end do
end if

do j=M,1,-1
   read(15,*) (q_flag(i,j),i=1,N)
end do

cont=0
do i=1,N
   do j=1,M
      if (q_flag(i,j) == 1) then
         cont=cont+1
         write(16,*) i,j 
	 end if
   end do
end do

write(6,*) 
write(6,'(1x,i5,1x,a)') cont,'cells selected'
write(6,*) 

end program qflag2ij
