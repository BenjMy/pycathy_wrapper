!-----------------------------------------------------------------------------
!
! MBBIO
!
! The module MBBIO manages the I/O between programs and binary files
! basin_b/basin_i.
!
!-----------------------------------------------------------------------------

module mbbio

use mpar

implicit none

private

public :: load_dtm,save_dtm,close_dtm,init_dtm,dtm_index_pr,set_dtm_index

public :: dtm_vh            !  1, vh: vegetation height
public :: dtm_rl_min        !  2, rl_min: minimum leaf resistance
public :: dtm_LAI           !  3, LAI: leaf area index
public :: dtm_vsc           !  4, vsc: vegetation storage capacity
public :: dtm_theta_s       !  5, theta_s: saturation soil moisture content
public :: dtm_eta           !  6, eta: pore-size distribution index
public :: dtm_psi_s         !  7, psi_s: saturated soil matrix potential
public :: dtm_K_s           !  8, K_s: saturated soil conductivity
public :: dtm_Z_low         !  9, Z_low: lower soil layer depth
public :: dtm_Z_up          ! 10, Z_up: upper soil layer depth
public :: dtm_theta_r       ! 11, theta_r: residual soil moisture content
public :: dtm_rl_max        ! 12, rl_max: maximum leaf resistance
public :: dtm_quota         ! 13, quota: cell elevation
public :: dtm_p_outflow_1   ! 14, pointer to the outflow cell in cardinal direction
public :: dtm_p_outflow_2   ! 15, pointer to the outflow cell in diagonal direction
public :: dtm_dmID          ! 16, drainage method identifier
public :: dtm_A_inflow      ! 17, A_inflow: upstream drainage area
public :: dtm_w_1           ! 18, weighting factor in cardinal direction
public :: dtm_w_2           ! 19, weighting factor in diagonal direction 
public :: dtm_sumdev_num    ! 20, numerator of cumulative deviation 
public :: dtm_local_slope_1 ! 21, local_slope_1: local slope in cardinal direction
public :: dtm_local_slope_2 ! 22, local_slope_2: local slope in diagonal direction
public :: dtm_Ws1_sf_1      ! 23, Ws1_sf_1: surface flow width W_sf for Q=1 
public :: dtm_Ws1_sf_2      ! 24, Ws1_sf_2: surface flow width W_sf for Q=1 
public :: dtm_b1_sf         ! 25, b1_sf: at-a-station exponent for W_sf
public :: dtm_kSs1_sf_1     ! 26, kSs1_sf: GS resistance coefficient kS for Q=1
public :: dtm_kSs1_sf_2     ! 27, kSs1_sf: GS resistance coefficient kS for Q=1
public :: dtm_y1_sf         ! 28, y1_sf: at-a-station exponent for kS
public :: dtm_ASk           ! 29, montgomery and dietrich function
public :: dtm_DN            ! 30, normalized divergence
public :: dtm_hcID          ! 31, hillslope/channel identifier
public :: dtm_epl_1         ! 32, epl_1: elemental path length in cardinal direction
public :: dtm_epl_2         ! 33, epl_2: elemental path length in diagonal direction
public :: dtm_lakes_map     ! 34, CATHY's variables
public :: dtm_zone          ! 35, CATHY's variables
public :: dtm_q_output      ! 36, CATHY's variables
public :: dtm_nrc           ! 37, number of rivulets/channels per cell
public :: dtm_sso_jun       ! 38, Strahler Stream Order memory junction
public :: dtm_sso           ! 39, Strahler Stream Order
public :: dtm_hso           ! 40, Horton Stream Order



public :: set_vh            !  1
public :: set_rl_min        !  2
public :: set_LAI           !  3
public :: set_vsc           !  4
public :: set_theta_s       !  5
public :: set_eta           !  6
public :: set_psi_s         !  7
public :: set_K_s           !  8
public :: set_Z_low         !  9
public :: set_Z_up          ! 10
public :: set_theta_r       ! 11
public :: set_rl_max        ! 12
public :: set_quota         ! 13
public :: set_p_outflow_1   ! 14
public :: set_p_outflow_2   ! 15
public :: set_dmID          ! 16
public :: set_A_inflow      ! 17
public :: set_w_1           ! 18
public :: set_w_2           ! 19
public :: set_sumdev_num    ! 20
public :: set_local_slope_1 ! 21
public :: set_local_slope_2 ! 22
public :: set_Ws1_sf_1      ! 23
public :: set_Ws1_sf_2      ! 24
public :: set_b1_sf         ! 25
public :: set_kSs1_sf_1     ! 26
public :: set_kSs1_sf_2     ! 27
public :: set_y1_sf         ! 28
public :: set_ASk           ! 29
public :: set_DN            ! 30
public :: set_hcID          ! 31
public :: set_epl_1         ! 32
public :: set_epl_2         ! 33
public :: set_lakes_map     ! 34
public :: set_zone          ! 35
public :: set_q_output      ! 36
public :: set_nrc           ! 37
public :: set_sso_jun       ! 38
public :: set_sso           ! 39
public :: set_hso           ! 40


! dtm record composition and length

integer(kind=ISP),parameter :: nbb_rsp=27
integer(kind=ISP),parameter :: nbb_rep=3
integer(kind=ISP),parameter :: nbb_isp=10
integer(kind=ISP),parameter :: nbb_iep=0
!--------------------------------------------------------------------
integer(kind=ISP),parameter :: &
nbrl=nbb_rsp*RSP+nbb_rep*REP+nbb_isp*ISP+nbb_iep*IEP+8
! note that 4 bytes are required to fill some gap when kind=8 is used
!--------------------------------------------------------------------

type cella
   real(kind=RSP) ::    vh            !  1
   real(kind=RSP) ::    rl_min        !  2
   real(kind=RSP) ::    LAI           !  3
   real(kind=RSP) ::    vsc           !  4
   real(kind=RSP) ::    theta_s       !  5
   real(kind=RSP) ::    eta           !  6
   real(kind=RSP) ::    psi_s         !  7
   real(kind=RSP) ::    K_s           !  8
   real(kind=RSP) ::    Z_low         !  9
   real(kind=RSP) ::    Z_up          ! 10
   real(kind=RSP) ::    theta_r       ! 11
   real(kind=RSP) ::    rl_max        ! 12
   real(kind=REP) ::    quota         ! 13
   integer(kind=ISP) :: p_outflow_1   ! 14
   integer(kind=ISP) :: p_outflow_2   ! 15
   integer(kind=ISP) :: dmID          ! 16
   real(kind=REP) ::    A_inflow      ! 17
   real(kind=RSP) ::    w_1           ! 18
   real(kind=RSP) ::    w_2           ! 19
   real(kind=REP) ::    sumdev_num    ! 20
   real(kind=RSP) ::    local_slope_1 ! 21
   real(kind=RSP) ::    local_slope_2 ! 22
   real(kind=RSP) ::    Ws1_sf_1      ! 23
   real(kind=RSP) ::    Ws1_sf_2      ! 24
   real(kind=RSP) ::    b1_sf         ! 25
   real(kind=RSP) ::    kSs1_sf_1     ! 26
   real(kind=RSP) ::    kSs1_sf_2     ! 27
   real(kind=RSP) ::    y1_sf         ! 28
   real(kind=RSP) ::    ASk           ! 29
   real(kind=RSP) ::    DN            ! 30
   integer(kind=ISP) :: hcID          ! 31
   real(kind=RSP) ::    epl_1         ! 32
   real(kind=RSP) ::    epl_2         ! 33
   integer(kind=ISP) :: lakes_map     ! 34
   integer(kind=ISP) :: zone          ! 35
   integer(kind=ISP) :: q_output      ! 36
   real(kind=RSP) ::    nrc           ! 37
   integer(kind=ISP) :: sso_jun       ! 38
   integer(kind=ISP) :: sso           ! 39
   integer(kind=ISP) :: hso           ! 40
end type

type(cella),dimension(:),allocatable :: dtm
integer(kind=ISP),dimension(:),allocatable :: indice

type(cella),parameter :: cella_iniziale=cella(0.0_RSP, & ! 1: vh
                                              0.0_RSP, & ! 2: rl_min
                                              0.0_RSP, & ! 3: LAI
                                              0.0_RSP, & ! 4: vsc
                                              0.0_RSP, & ! 5: theta_s
                                              0.0_RSP, & ! 6: eta
                                              0.0_RSP, & ! 7: psi_s
                                              0.0_RSP, & ! 8: K_s
                                              0.0_RSP, & ! 9: Z_low
                                              0.0_RSP, & ! 10: Z_up
                                              0.0_RSP, & ! 11: theta_r
                                              0.0_RSP, & ! 12: rl_max
                                              0.0_REP, & ! 13: quota
                                                0_ISP, & ! 14: p_outflow_1
                                                0_ISP, & ! 15: p_outflow_2
                                                0_ISP, & ! 16: dmID
                                              0.0_REP, & ! 17: A_inflow
                                              0.0_RSP, & ! 18: w_1
                                              0.0_RSP, & ! 19: w_2
                                              0.0_REP, & ! 20: sumdev_num
                                              0.0_RSP, & ! 21: local_slope_1
                                              0.0_RSP, & ! 22: local_solpe_2
                                              0.0_RSP, & ! 23: Ws1_sf_1
                                              0.0_RSP, & ! 24: Ws1_sf_2
                                              0.0_RSP, & ! 25: b1_sf
                                              0.0_RSP, & ! 26: kSs1_sf_1
                                              0.0_RSP, & ! 27: kSs1_sf_2
                                              0.0_RSP, & ! 28: y1_sf
                                              0.0_RSP, & ! 29: ASk
                                              0.0_RSP, & ! 30: DN
                                                0_ISP, & ! 31: hcID
                                              0.0_RSP, & ! 32: epl_1
                                              0.0_RSP, & ! 33: epl_2
                                                0_ISP, & ! 34: lakes_map
                                                1_ISP, & ! 35: zone
                                                0_ISP, & ! 36: q_output
                                              0.0_RSP, & ! 37: nrc
												0_ISP, & ! 38: sso_jun											    
												0_ISP, & ! 39: sso
												0_ISP)   ! 40: hso
contains

!----------------------------------------------------------------------------- 
!
! INIT_DTM
!
! The subroutine INIT_DTM initializes the dtm data in the variable indice.
!  
!-----------------------------------------------------------------------------

subroutine init_dtm(n_rec)
integer(kind=ISP),intent(in) :: n_rec

if (.not. allocated(indice)) then
   allocate(indice(M*N))
else
   write(6,*) "can't allocate indice"
   stop
end if
indice=0
if (.not. allocated(dtm)) then
   allocate(dtm(n_rec))
else
   write(6,*) "can't allocate dtm"
   stop
end if

dtm=cella_iniziale

end subroutine init_dtm

!----------------------------------------------------------------------------- 
!
! LOAD_DTM
!
! The subroutine LOAD_DTM loads the dtm data from the binary files 
! basin_b/basin_i.
!  
!-----------------------------------------------------------------------------

subroutine load_dtm(NomeDemFile,NomeIndice)
character(len=*),intent(in) :: NomeDemFile,NomeIndice
integer(kind=ISP) :: n_rec,i
open(22,file=NomeDemFile,access='direct',status='old', &
form='unformatted',recl=nbrl)
open(23,file=NomeIndice,access='direct',status='old', &
form='unformatted',recl=4)
read(23,rec=1) n_rec
if (.not. allocated(indice)) then
   allocate(indice(M*N))
else
   write(6,*) "allocated(indice)",ubound(indice,1)
   stop
endif
if (.not. allocated(dtm)) then
   allocate(dtm(n_rec))
else
   write(6,*) "allocated(dtm)",ubound(dtm,1)
   stop
endif
do i=1, M*N
   read(23,rec=i+1) indice(i)
end do
do i=1, n_rec
    read(22,rec=i) dtm(i)
end do
close(22)
close(23)
end subroutine load_dtm

!----------------------------------------------------------------------------- 
!
! SAVE_DTM
!
! The subroutine SAVE_DTM saves the dtm data in the binary files 
! basin_b/basin_i.
!  
!-----------------------------------------------------------------------------

subroutine save_dtm(NomeDemFile,NomeIndice)
character(len=*),intent(in) :: NomeDemFile
character(len=*),intent(in),optional :: NomeIndice
integer(kind=ISP) :: i
open(22,file=NomeDemFile,access='direct',status='unknown', &
form='unformatted',recl=nbrl)
do i=1,ubound(dtm,1)
   write(22,rec=i) dtm(i)
end do
close(22)
if (present(NomeIndice)) then
   open(23,file=NomeIndice,access='direct',status='unknown', &
   form='unformatted',recl=4)
   write(23,rec=1) ubound(dtm,1)
   do i=1,M*N
      write(23,rec=i+1) indice(i)
   end do
   close(23)
end if
deallocate(indice)
deallocate(dtm)
end subroutine save_dtm

!----------------------------------------------------------------------------- 
!
! CLOSE_DTM
!
! The subroutine CLOSE_DTM deallocate the dtm data stored in indice and
! dtm variables
!  
!-----------------------------------------------------------------------------

subroutine close_dtm()
deallocate(indice)
deallocate(dtm)
end subroutine close_dtm

!-----------------------------------------------------------------------------
!
! DTM_INDEX_PR
!
! The function DTM_INDEX_PR provides the index of records in basin_b 
! for parameter read. 
!
!-----------------------------------------------------------------------------

function dtm_index_pr(i)
integer(kind=ISP) :: dtm_index_pr
integer(kind=ISP),intent(in) :: i
if (i <= 0 .or. i > ubound(indice,1)) stop &
"error in dtm_index_pr: input out of range!"
dtm_index_pr=indice(i)
end function dtm_index_pr

!-----------------------------------------------------------------------------
!
! DTM_INDEX_PS
!
! The function DTM_INDEX_PS provides the index of records in basin_b 
! for parameter set. 
!
!-----------------------------------------------------------------------------

function dtm_index_ps(i)
integer(kind=ISP) :: dtm_index_ps
integer(kind=ISP),intent(in) :: i
if (i <= 0 .or. i > ubound(indice,1)) stop &
"error in dtm_index_ps: input out of range"
if (indice(i) <= 0) stop & 
"error in dtm_index_ps: cell write not allowed"
dtm_index_ps=indice(i)
end function dtm_index_ps

!-----------------------------------------------------------------------------
!
! SET_DTM_INDEX
!
! The function SET_DTM_INDEX sets the index of records in basin_b. 
!
!-----------------------------------------------------------------------------

subroutine set_dtm_index(i_basin,i_rec)
integer(kind=ISP),intent(in) :: i_basin,i_rec
if (i_basin <= 0 .or. i_basin > ubound(indice,1)) stop &
"error in set_dtm_index: input out of range!"
indice(i_basin)=i_rec
end subroutine set_dtm_index

!-----------------------------------------------------------------------------
!
! Subroutines for setting parameter values.
!
!-----------------------------------------------------------------------------

subroutine set_vh(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%vh=x
end subroutine

subroutine set_rl_min(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%rl_min=x
end subroutine

subroutine set_LAI(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%LAI=x
end subroutine

subroutine set_vsc(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%vsc=x
end subroutine

subroutine set_theta_s(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%theta_s=x
end subroutine

subroutine set_eta(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%eta=x
end subroutine

subroutine set_psi_s(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%psi_s=x
end subroutine

subroutine set_K_s(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%K_s=x
end subroutine

subroutine set_Z_low(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%Z_low=x
end subroutine

subroutine set_Z_up(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%Z_up=x
end subroutine

subroutine set_theta_r(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%theta_r=x
end subroutine

subroutine set_rl_max(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%rl_max=x
end subroutine

subroutine set_quota(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=REP),intent(in) :: x
dtm(dtm_index_ps(i))%quota=x
end subroutine

subroutine set_p_outflow_1(i,j)
integer(kind=ISP),intent(in) :: i
integer(kind=ISP),intent(in) :: j
dtm(dtm_index_ps(i))%p_outflow_1=j
end subroutine

subroutine set_p_outflow_2(i,j)
integer(kind=ISP),intent(in) :: i
integer(kind=ISP),intent(in) :: j
dtm(dtm_index_ps(i))%p_outflow_2=j
end subroutine

subroutine set_dmID(i,j)
integer(kind=ISP),intent(in) :: i
integer(kind=ISP),intent(in) :: j
dtm(dtm_index_ps(i))%dmID=j
end subroutine

subroutine set_A_inflow(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=REP),intent(in) :: x
dtm(dtm_index_ps(i))%A_inflow=x
end subroutine

subroutine set_w_1(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%w_1=x
end subroutine

subroutine set_w_2(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%w_2=x
end subroutine

subroutine set_sumdev_num(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=REP),intent(in) :: x
dtm(dtm_index_ps(i))%sumdev_num=x
end subroutine

subroutine set_local_slope_1(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%local_slope_1=x
end subroutine

subroutine set_local_slope_2(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%local_slope_2=x
end subroutine

subroutine set_Ws1_sf_1(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%Ws1_sf_1=x
end subroutine

subroutine set_Ws1_sf_2(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%Ws1_sf_2=x
end subroutine

subroutine set_b1_sf(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%b1_sf=x
end subroutine

subroutine set_kSs1_sf_1(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%kSs1_sf_1=x
end subroutine

subroutine set_kSs1_sf_2(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%kSs1_sf_2=x
end subroutine

subroutine set_y1_sf(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%y1_sf=x
end subroutine

subroutine set_ASk(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%ASk=x
end subroutine

subroutine set_DN(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%DN=x
end subroutine

subroutine set_hcID(i,j)
integer(kind=ISP),intent(in) :: i
integer(kind=ISP),intent(in) :: j
dtm(dtm_index_ps(i))%hcID=j
end subroutine

subroutine set_epl_1(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%epl_1=x
end subroutine

subroutine set_epl_2(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%epl_2=x
end subroutine

subroutine set_lakes_map(i,j)
integer(kind=ISP),intent(in) :: i
integer(kind=ISP),intent(in) :: j
dtm(dtm_index_ps(i))%lakes_map=j
end subroutine

subroutine set_zone(i,j)
integer(kind=ISP),intent(in) :: i
integer(kind=ISP),intent(in) :: j
dtm(dtm_index_ps(i))%zone=j
end subroutine

subroutine set_q_output(i,j)
integer(kind=ISP),intent(in) :: i
integer(kind=ISP),intent(in) :: j
dtm(dtm_index_ps(i))%q_output=j
end subroutine

subroutine set_nrc(i,x)
integer(kind=ISP),intent(in) :: i
real(kind=RSP),intent(in) :: x
dtm(dtm_index_ps(i))%nrc=x
end subroutine

subroutine set_sso_jun(i,j)
integer(kind=ISP),intent(in) :: i
integer(kind=ISP),intent(in) :: j
dtm(dtm_index_ps(i))%sso_jun=j
end subroutine

subroutine set_sso(i,j)
integer(kind=ISP),intent(in) :: i
integer(kind=ISP),intent(in) :: j
dtm(dtm_index_ps(i))%sso=j
end subroutine

subroutine set_hso(i,j)
integer(kind=ISP),intent(in) :: i
integer(kind=ISP),intent(in) :: j
dtm(dtm_index_ps(i))%hso=j
end subroutine

!-----------------------------------------------------------------------------
!  
! Functions for reading parameter values.
!
!-----------------------------------------------------------------------------

function dtm_vh(i)
real(kind=RSP) :: dtm_vh
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_vh=dtm(k)%vh
else
   stop "h"
end if
end function

function dtm_rl_min(i)
real(kind=RSP) :: dtm_rl_min
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_rl_min=dtm(k)%rl_min
else
   stop "rl_min"
end if
end function

function dtm_LAI(i)
real(kind=RSP) :: dtm_LAI
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_LAI=dtm(k)%LAI
else
   stop "rl_min"
end if
end function

function dtm_vsc(i)
real(kind=RSP) :: dtm_vsc
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_vsc=dtm(k)%vsc
else
   stop "vsc"
end if
end function

function dtm_theta_s(i)
real(kind=RSP) :: dtm_theta_s
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_theta_s=dtm(k)%theta_s
else
   stop "theta_s"
end if
end function

function dtm_eta(i)
real(kind=RSP) :: dtm_eta
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_eta=dtm(k)%eta
else
   stop "eta"
end if
end function

function dtm_psi_s(i)
real(kind=RSP) :: dtm_psi_s
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_psi_s=dtm(k)%psi_s
else
   stop "psi_s"
end if
end function

function dtm_K_s(i)
real(kind=RSP) :: dtm_K_s
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_K_s=dtm(k)%K_s
else
   stop "K_s"
end if
end function

function dtm_Z_low(i)
real(kind=RSP) :: dtm_Z_low
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_Z_low=dtm(k)%Z_low
else
   stop "Z_low"
end if
end function

function dtm_Z_up(i)
real(kind=RSP) :: dtm_Z_up
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_Z_up=dtm(k)%Z_up
else
   stop "Z_low"
end if
end function

function dtm_theta_r(i)
real(kind=RSP) :: dtm_theta_r
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_theta_r=dtm(k)%theta_r
else
   stop "theta_r"
end if
end function

function dtm_rl_max(i)
real(kind=RSP) :: dtm_rl_max
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_rl_max=dtm(k)%rl_max
else
   stop "rl_max"
end if
end function

function dtm_quota(i)
real(kind=REP) :: dtm_quota
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_quota=dtm(k)%quota
else
   dtm_quota=-1.0
end if
end function

function dtm_p_outflow_1(i)
integer(kind=ISP) :: dtm_p_outflow_1
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_p_outflow_1=dtm(k)%p_outflow_1
else
   stop "dtm_p_outflow_1"
end if
end function

function dtm_p_outflow_2(i)
integer(kind=ISP) :: dtm_p_outflow_2
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_p_outflow_2=dtm(k)%p_outflow_2
else
   stop "dtm_p_outflow_2"
end if
end function

function dtm_dmID(i)
integer(kind=ISP) :: dtm_dmID
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_dmID= dtm(k)%dmID
else
   stop "dtm_dmID "
end if
end function

function dtm_A_inflow(i)
real(kind=REP) :: dtm_A_inflow
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_A_inflow=dtm(k)%A_inflow
else
   stop "dtm_A_inflow"
end if
end function

function dtm_w_1(i)
real(kind=RSP) :: dtm_w_1
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_w_1=dtm(k)%w_1
else
   stop "dtm_w_1"
end if
end function

function dtm_w_2(i)
real(kind=RSP) :: dtm_w_2
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_w_2=dtm(k)%w_2
else
   stop "dtm_w_2"
end if
end function

function dtm_sumdev_num(i)
real(kind=REP) :: dtm_sumdev_num
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_sumdev_num=dtm(k)%sumdev_num
else
   stop "dtm_sumdev_num"
end if
end function

function dtm_local_slope_1(i)
real(kind=RSP) :: dtm_local_slope_1
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_local_slope_1=dtm(k)%local_slope_1
else
   stop "dtm_local_slope_1"
end if
end function

function dtm_local_slope_2(i)
real(kind=RSP) :: dtm_local_slope_2
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then
   dtm_local_slope_2=dtm(k)%local_slope_2
else
   stop "dtm_local_slope_2"
end if
end function

function dtm_Ws1_sf_1(i)
real(kind=RSP) :: dtm_Ws1_sf_1
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_Ws1_sf_1=dtm(k)%Ws1_sf_1
else
   stop "dtm_Ws1_sf_1"
end if
end function

function dtm_Ws1_sf_2(i)
real(kind=RSP) :: dtm_Ws1_sf_2
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_Ws1_sf_2=dtm(k)%Ws1_sf_2
else
   stop "dtm_Ws1_sf_2"
end if
end function

function dtm_b1_sf(i)
real(kind=RSP) :: dtm_b1_sf
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_b1_sf=dtm(k)%b1_sf
else
   stop "dtm_b1_sf"
end if
end function

function dtm_kSs1_sf_1(i)
real(kind=RSP) :: dtm_kSs1_sf_1
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_kSs1_sf_1=dtm(k)%kSs1_sf_1
else
   stop "dtm_kSs1_sf_1"
end if
end function

function dtm_kSs1_sf_2(i)
real(kind=RSP) :: dtm_kSs1_sf_2
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_kSs1_sf_2=dtm(k)%kSs1_sf_2
else
   stop "dtm_kSs1_sf_2"
end if
end function

function dtm_y1_sf(i)
real(kind=RSP) :: dtm_y1_sf
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_y1_sf=dtm(k)%y1_sf
else
   stop "dtm_y1_sf"
end if
end function

function dtm_ASk(i)
real(kind=RSP) :: dtm_ASk
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_ASk= dtm(k)%ASk
else
   stop "dtm_ASk"
end if
end function

function dtm_DN(i)
real(kind=RSP) :: dtm_DN
integer(kind=ISP),intent(in):: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_DN= dtm(k)%DN
else
   stop "dtm_DN"
end if
end function

function dtm_hcID(i)
integer(kind=ISP) :: dtm_hcID
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_hcID= dtm(k)%hcID
else
   stop "dtm_hcID "
end if
end function

function dtm_epl_1(i)
real(kind=RSP) :: dtm_epl_1
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_epl_1=dtm(k)%epl_1
else
   stop "dtm_epl_1"
end if
end function

function dtm_epl_2(i)
real(kind=RSP) :: dtm_epl_2
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_epl_2=dtm(k)%epl_2
else
   stop "dtm_epl_2"
end if
end function

function dtm_lakes_map(i)
integer(kind=ISP) :: dtm_lakes_map
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_lakes_map= dtm(k)%lakes_map
else
   stop "dtm_lakes_map "
end if
end function

function dtm_zone(i)
integer(kind=ISP) :: dtm_zone
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_zone= dtm(k)%zone
else
   stop "dtm_zone "
end if
end function

function dtm_q_output(i)
integer(kind=ISP) :: dtm_q_output
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_q_output= dtm(k)%q_output
else
   stop "dtm_q_output "
end if
end function

function dtm_nrc(i)
real(kind=RSP) :: dtm_nrc
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_nrc=dtm(k)%nrc
else
   stop "dtm_nrc"
end if
end function

function dtm_sso_jun(i)
integer(kind=ISP) :: dtm_sso_jun
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_sso_jun= dtm(k)%sso_jun
else
   stop "dtm_sso_jun"
end if
end function

function dtm_sso(i)
integer(kind=ISP) :: dtm_sso
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_sso= dtm(k)%sso
else
   stop "dtm_sso"
end if
end function

function dtm_hso(i)
integer(kind=ISP) :: dtm_hso
integer(kind=ISP),intent(in) :: i
integer(kind=ISP) :: k
k=dtm_index_pr(i)
if (k /= 0) then  
   dtm_hso= dtm(k)%hso
else
   stop "dtm_hso"
end if
end function

end module mbbio
