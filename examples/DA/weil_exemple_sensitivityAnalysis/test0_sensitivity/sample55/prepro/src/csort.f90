!-----------------------------------------------------------------------------
!
! CSORT
!
! The subroutine CSORT (Cells Sort) sorts the cells of a catchment
! in the order of descending elevations. The binary file qoi contains
! the number of processed cells in the first record and the pointers
! i_basin_qo corresponding to the sorted cell elevations qo in the
! subsequent records.
!
!-----------------------------------------------------------------------------

subroutine csort()

use mpar
use mbbio

implicit none

integer(kind=ISP),allocatable :: i_basin_qo(:)
real(kind=REP),allocatable :: qo(:)

integer(kind=ISP) i,j,i_basin,n_rec,i_qo,l

allocate(i_basin_qo(N_celle))
allocate(qo(N_celle))

! read cell elevations from the binary file basin_b

i_qo=0
do i=1,N
   do j=1,M
      i_basin=(i-1)*M+j
      if (dtm_index_pr(i_basin) == 0) cycle
      i_qo=i_qo+1
      i_basin_qo(i_qo)=i_basin
      qo(i_qo)=dtm_quota(i_basin)
   end do
end do

! sort cells in the order of ascending elevations using the algorithm Quicksort

call qsort(N_celle,qo,i_basin_qo)

! open and write the binary file qoi

open(16,file='qoi',access='direct',status='unknown',form='unformatted',recl=4)
n_rec=0
do l=N_celle,1,-1 ! turn on descending order
   n_rec=n_rec+1
   write(16,rec=n_rec+1) i_basin_qo(l)
end do
write(16,rec=1) n_rec
close(16)

deallocate(i_basin_qo)
deallocate(qo)

return
end subroutine csort
