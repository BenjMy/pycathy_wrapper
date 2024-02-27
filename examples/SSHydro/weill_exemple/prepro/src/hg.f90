!-----------------------------------------------------------------------------
!
! HG (Hydraulic Geometry)
!
! The subroutine HG determines the hillslope an channel cells within the 
! catchment and assign the hydraulic geometry parameters.
!
!-----------------------------------------------------------------------------

subroutine hg()

use mpar
use mbbio

implicit none

integer(kind=ISP) hcID
integer(kind=ISP) i_basin,i_basin_chiusura
integer(kind=ISP) n_o_quota

real(kind=REP) A_inflow,A_outflow
real(kind=RSP) w_1,w_2
real(kind=RSP) Ws1_sf_1,Ws1_sf_2,b1_sf,kSs1_sf_1,kSs1_sf_2,y1_sf,nrc

! open the binary file qoi

open(16,file='qoi',access='direct',status='unknown',form='unformatted',recl=4)
open(17,file='qoi_a',status='unknown')

! identify the outlet cell (the lowest catchment cell)
! read(16,rec=N_celle+1) i_basin_chiusura
! write(6,*) 'i_basin_chiusura = ' , i_basin_chiusura

write(17,*) N_celle

do n_o_quota=1,N_celle
   read(16,rec=n_o_quota+1) i_basin   
   WRITE(17,*) i_basin
   A_inflow=dtm_A_inflow(i_basin)       ! 17
   w_1=dtm_w_1(i_basin)                 ! 18
   w_2=dtm_w_2(i_basin)                 ! 19
   Ws1_sf_1=dtm_Ws1_sf_1(i_basin)       ! 23
   Ws1_sf_2=dtm_Ws1_sf_2(i_basin)       ! 24
   b1_sf=dtm_b1_sf(i_basin)             ! 25
   kSs1_sf_1=dtm_kSs1_sf_1(i_basin)     ! 26
   kSs1_sf_2=dtm_kSs1_sf_2(i_basin)     ! 27
   y1_sf=dtm_y1_sf(i_basin)             ! 28
   hcID=dtm_hcID(i_basin)               ! 31
   nrc=dtm_nrc(i_basin)                 ! 37

   A_outflow=A_inflow+delta_x*delta_y

   call HydraulicGeometry(hcID,A_outflow,w_1,w_2, &
   Ws1_sf_1,Ws1_sf_2,b1_sf,kSs1_sf_1,kSs1_sf_2,y1_sf,nrc)
    
   call set_Ws1_sf_1(i_basin,Ws1_sf_1)     ! 23
   call set_Ws1_sf_2(i_basin,Ws1_sf_2)     ! 24
   call set_b1_sf(i_basin,b1_sf)           ! 25
   call set_KSs1_sf_1(i_basin,KSs1_sf_1)   ! 26
   call set_KSs1_sf_2(i_basin,KSs1_sf_2)   ! 27
   call set_y1_sf(i_basin,y1_sf)           ! 28
   call set_nrc(i_basin,nrc)               ! 37

end do

close(16)
close(17)

return
end subroutine hg


subroutine HydraulicGeometry(hcID,A_outflow,w_1,w_2, &
Ws1_sf_1,Ws1_sf_2,b1_sf,kSs1_sf_1,kSs1_sf_2,y1_sf,nrc)

use mpar
implicit none
integer(kind=ISP),intent(in) :: hcID
real(kind=REP),intent(in) :: A_outflow
real(kind=RSP),intent(in) :: w_1,w_2
real(kind=RSP),intent(inout) :: Ws1_sf_1,Ws1_sf_2,b1_sf
real(kind=RSP),intent(inout) :: kSs1_sf_1,kSs1_sf_2,y1_sf
real(kind=RSP),intent(inout) :: nrc
real(kind=REP) :: RA

if (hcID.eq.0) then ! rill flow
   RA=A_outflow/As_rf
   if (b1_sf.eq.0.0E0) b1_sf=b1_rf
   if (Ws1_sf_1.eq.0.0E0.and.abs(w_1).gt.epsilon(w_1)) &
   Ws1_sf_1=Wsf_rf*Qsf_rf**(-b1_sf)*(RA*w_1)**(w_rf*(b2_rf-b1_sf))
   if (Ws1_sf_2.eq.0.0E0.and.abs(w_2).gt.epsilon(w_2)) &
   Ws1_sf_2=Wsf_rf*Qsf_rf**(-b1_sf)*(RA*w_2)**(w_rf*(b2_rf-b1_sf))
   if (y1_sf.eq.0.0E0) y1_sf=y1_rf
   if (kSs1_sf_1.eq.0.0E0.and.abs(w_1).gt.epsilon(w_1)) &
   kSs1_sf_1=kSsf_rf*Qsf_rf**(-y1_sf)*(RA*w_1)**(w_rf*(y2_rf-y1_sf))
   if (kSs1_sf_2.eq.0.0E0.and.abs(w_2).gt.epsilon(w_2)) &
   kSs1_sf_2=kSsf_rf*Qsf_rf**(-y1_sf)*(RA*w_2)**(w_rf*(y2_rf-y1_sf))
   if (nrc.eq.0.0E0) nrc=delta_x/dr
else ! if (hcID.eq.1) then ! channel flow
   RA=A_outflow/As_cf
   if (b1_sf.eq.0.0E0) b1_sf=b1_cf
   if (Ws1_sf_1.eq.0.0E0.and.abs(w_1).gt.epsilon(w_1)) &
   Ws1_sf_1=Wsf_cf*Qsf_cf**(-b1_sf)*(RA*w_1)**(w_cf*(b2_cf-b1_sf))
   if (Ws1_sf_2.eq.0.0E0.and.abs(w_2).gt.epsilon(w_2)) &
   Ws1_sf_2=Wsf_cf*Qsf_cf**(-b1_sf)*(RA*w_2)**(w_cf*(b2_cf-b1_sf))
   if (y1_sf.eq.0.0E0) y1_sf=y1_cf
   if (kSs1_sf_1.eq.0.0E0.and.abs(w_1).gt.epsilon(w_1)) &
   kSs1_sf_1=kSsf_cf*Qsf_cf**(-y1_sf)*(RA*w_1)**(w_cf*(y2_cf-y1_sf))
   if (kSs1_sf_2.eq.0.0E0.and.abs(w_2).gt.epsilon(w_2)) &
   kSs1_sf_2=kSsf_cf*Qsf_cf**(-y1_sf)*(RA*w_2)**(w_cf*(y2_cf-y1_sf))
   if (nrc.eq.0.0E0) nrc=1.0
end if

end subroutine HydraulicGeometry
