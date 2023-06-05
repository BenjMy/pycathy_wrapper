!-----------------------------------------------------------------------------
!
! DEPIT
!
! The subroutine DEPIT checks that each catchment cell has a lower cell
! among the neighbouring cells. When this is not verified the processed cell
! is raised of the quantity eps above le lowest neighbouring cells.
! The procedure is iterated until all the pits are removed. The number of
! modifications is recorded in the variable n_modifiche.
!
!-----------------------------------------------------------------------------

subroutine depit()

use mpar
use mbbio

implicit none

integer(kind=ISP),allocatable :: i_basin_pit_1(:),i_basin_pit_2(:)
real(kind=REP),allocatable :: q_pit_n(:)
integer(kind=ISP),allocatable :: i_pit_n(:)

integer(kind=ISP) i,j,i_basin,n_rec,i_qo,jr
integer(kind=ISP) ii,jj,ii_basin,ll
integer(kind=ISP) i_basin_l,i_basin_sl
real(kind=REP) eps
real(kind=REP) quota_o_l,quota_o_sl
real(kind=REP) quota_c,quota_c_c,quota_c_c_min
integer(kind=ISP) N_pits,NN_pits,n_pit,nn_pit,n_modifiche,nn_modifiche

allocate(i_basin_pit_1(N_celle),i_basin_pit_2(N_celle))
allocate(q_pit_n(N_celle))
allocate(i_pit_n(N_celle))

! quantity of which cells are raised

eps=pt*delta_x

! open the binary file qoi

open(16,file='qoi',access='direct',status='unknown',form='unformatted',recl=4)

read(16,rec=N_celle+1) i_basin_l
quota_o_l=dtm_quota(i_basin_l)

! verify that the catchment has one outlet cell

read(16,rec=(N_celle-1)+1) i_basin_sl
quota_o_sl=dtm_quota(i_basin_sl)

! write(6,*) "epsilon(quota_o_l)*delta_x = ",epsilon(quota_o_l)*delta_x

if ((quota_o_sl-quota_o_l).lt.epsilon(quota_o_l)*delta_x) then
   write(6,*) 'catchment with more than one outlet cell!'
   write(6,*) 'i_basin_l=',i_basin_l
   write(6,*) 'i_basin_sl=',i_basin_sl
   stop
end if

do i_qo=1,N_celle
   read(16,rec=i_qo+1) i_basin
   quota_c=dtm_quota(i_basin)
   if (quota_c < 0.0E0) cycle
   i_basin_pit_1(N_celle-i_qo+1)=i_basin
end do
N_pits=N_celle

close(16)

n_modifiche=0
nn_modifiche=1
modCycle: do while (nn_modifiche.gt.0)
   nn_modifiche=0
   nn_pit=0
   do ll=1,N_celle
      i_pit_n(ll)=0
   end do
   pitCycle: do n_pit=1,N_pits
      i_basin=i_basin_pit_1(n_pit)
      quota_c=dtm_quota(i_basin)
      if (i_basin.eq.i_basin_l) cycle pitCycle
      jr=mod(i_basin,M)
      if (jr.ne.0) then
         j=jr
         i=(i_basin-j)/M+1
      else
         j=M
         i=i_basin/M
      end if
      quota_c_c_min=huge(quota_c_c_min)
      do ii=i-1,i+1
         do jj=j-1,j+1
            if ((ii.eq.i).and.(jj.eq.j)) cycle
            if ((ii.eq.0).or.(ii.eq.N+1)) cycle
            if ((jj.eq.0).or.(jj.eq.M+1)) cycle
            ii_basin=(ii-1)*M+jj
            quota_c_c=dtm_quota(ii_basin)
            if (quota_c_c < 0.0E0) cycle
            if (quota_c_c.lt.quota_c) cycle pitCycle
            if (quota_c_c.lt.quota_c_c_min) then
               quota_c_c_min=quota_c_c
            end if
         end do
      end do
      if (quota_c.le.quota_c_c_min) then
         quota_c=quota_c_c_min+eps
         n_modifiche=n_modifiche+1
         nn_modifiche=nn_modifiche+1
         call set_quota(i_basin,quota_c)
         do ii=i-1,i+1
            do jj=j-1,j+1
               if ((ii.eq.i).and.(jj.eq.j)) cycle
               if ((ii.eq.0).or.(ii.eq.N+1)) cycle
               if ((jj.eq.0).or.(jj.eq.M+1)) cycle
               ii_basin=(ii-1)*M+jj
               n_rec=dtm_index_pr(ii_basin)
               if (n_rec.eq.0) cycle
               if (i_pit_n(n_rec).eq.0) then
                  i_pit_n(n_rec)=1
                  nn_pit=nn_pit+1
                  i_basin_pit_2(nn_pit)=ii_basin
               end if
            end do
         end do
      end if
   end do pitCycle
   write(6,*) 'dem modifications = ',nn_modifiche
   NN_pits=nn_pit
   if (NN_pits.gt.N_celle) then
      write(6,*) 'NN_pits > N_celle!'
      stop
   end if
   do nn_pit=1,NN_pits
      i_basin=i_basin_pit_2(nn_pit)
      quota_c=dtm_quota(i_basin)
      q_pit_n(nn_pit)=quota_c
   end do
   call qsort(NN_pits,q_pit_n,i_basin_pit_2)
   do nn_pit=1,NN_pits
      i_basin_pit_1(nn_pit)=i_basin_pit_2(nn_pit)
   end do
   N_pits=NN_pits
end do modCycle

! write the number of modifications n_modifiche

write(6,*) 'dem modifications = ',n_modifiche,'(total)'

deallocate(i_basin_pit_1,i_basin_pit_2)
deallocate(q_pit_n)
deallocate(i_pit_n)

return
end subroutine depit
