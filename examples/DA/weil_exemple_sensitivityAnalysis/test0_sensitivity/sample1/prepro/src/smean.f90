!-----------------------------------------------------------------------------
!
! SMEAN
!
! The subroutine SMEAN calculates the mean catchment slope based on 
! maximum (theoretical) values obtained from facet analysis.
!
!-----------------------------------------------------------------------------

subroutine smean()

use mpar
use mbbio

implicit none

integer(kind=ISP) i,j,i_basin,jr,l
integer(kind=ISP) ii,jj,ii_basin
integer(kind=ISP) i_basin_chiusura
integer(kind=ISP) n_o_quota
real(kind=REP) quota_ciijj,quota_max
real(kind=REP) e(9)
real(kind=REP) e0,e1,e2,e1_fmax,e2_fmax,emin,s_max,r_max,s_max_facet,r
real(kind=REP) n_s_max,min_s_max,max_s_max,sum_s_max

! open the binary file qoi

open(16,file='qoi',access='direct',status='unknown',form='unformatted',recl=4)

! identify the outlet cell (the lowest catchment cell)

read(16,rec=N_celle+1) i_basin_chiusura

! write(6,*) 'i_basin_chiusura = ',i_basin_chiusura

n_s_max=+0.0D+00
min_s_max=+1.0D+10
max_s_max=-1.0D+10
sum_s_max=0.0D+00

do n_o_quota=1,N_celle

   read(16,rec=n_o_quota+1) i_basin
   ! if drainage direction are set up, there are not altered
   if (dtm_p_outflow_1(i_basin).ne.0.or.dtm_p_outflow_2(i_basin).ne.0) cycle
   if (i_basin.eq.i_basin_chiusura) cycle ! features are assigned apart
   jr=mod(i_basin,M)
   if (jr.ne.0) then
      j=jr
      i=(i_basin-j)/M+1
   else
      j=M
      i=i_basin/M
   end if

   l=0
   do ii=i-1,i+1
      do jj=j-1,j+1
         l=l+1
         e(l)=0.0_REP
         if (ii.eq.0.or.ii.eq.N+1) cycle
         if (jj.eq.0.or.jj.eq.M+1) cycle
         ii_basin=(ii-1)*M+jj
         quota_ciijj=dtm_quota(ii_basin)
         if (quota_ciijj < 0.0_REP) cycle
         e(l)=quota_ciijj
      end do
   end do

   e0=e(5)
   s_max=0.0_REP
!
! triangle 021
!
   e1=e(2)
   e2=e(1)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) s_max=s_max_facet
   end if
!
! triangle 023
!
   e1=e(2)
   e2=e(3)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) s_max=s_max_facet
   end if
!
! triangle 063
!
   e1=e(6)
   e2=e(3)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) s_max=s_max_facet
   end if
!
! triangle 069
!
   e1=e(6)
   e2=e(9)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) s_max=s_max_facet
   end if
!
! triangle 089
!
   e1=e(8)
   e2=e(9)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) s_max=s_max_facet
   end if
!
! triangle 087
!
   e1=e(8)
   e2=e(7)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) s_max=s_max_facet
   end if
!
! triangle 047
!
   e1=e(4)
   e2=e(7)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) s_max=s_max_facet
   end if
!
! triangle 041
!
   e1=e(4)
   e2=e(1)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) s_max=s_max_facet
   end if


   if (s_max.ge.0.0_REP) then
      n_s_max=n_s_max+1.0
      if (s_max.lt.min_s_max) min_s_max=s_max
      if (s_max.gt.max_s_max) max_s_max=s_max
      sum_s_max=sum_s_max+s_max
   else if (s_max.eq.0.0_REP) then
      write(6,*) 's_max = 0 at i_basin =',i_basin
      n_s_max=n_s_max+1.0
      if (s_max.lt.min_s_max) min_s_max=s_max
      if (s_max.gt.max_s_max) max_s_max=s_max
      sum_s_max=sum_s_max+s_max
   else ! if (s_max.lt.0.0_REP) then
      write(6,*) 's_max < 0 at i_basin =',i_basin
      stop
   end if

end do

mean_s_max=sum_s_max/n_s_max

write(6,'(1x,a,f12.9,a,f12.9,a,f12.9,a)') &
'mean (min,max) facet slope = ',mean_s_max,' (',min_s_max,',',max_s_max,')'

close(16)

return
end subroutine smean
