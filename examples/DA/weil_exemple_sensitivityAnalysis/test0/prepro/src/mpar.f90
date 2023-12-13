!-----------------------------------------------------------------------------
!
! MPAR
!
! Module containing definitions of parameters assigned in the parameter 
! file <parfile>. This module also contains three subroutines:
!
! 1) RPARFILE
!
! This subroutine reads <parfile> recognizing parameter descriptions 
! (before the "=" mark) and free-format parameters (after the "=" mark).
!
! 2) RROW
!
! This subroutine reads and analyze each row of <parfile>.
!
! 3) WPARFILE
!
! This subroutine rewrites <parfile> with specified changes.
!
!-----------------------------------------------------------------------------
 
module mpar
 
implicit none

integer,parameter :: RSP = selected_real_kind(6,20)  ! Real Single Precision
integer,parameter :: REP = selected_real_kind(14,27) ! Real Extended Precision
integer,parameter :: ISP = selected_int_kind(9)  ! Integer Single Precision
integer,parameter :: IEP = selected_int_kind(9) ! Integer Extended Precision
 
integer(kind=ISP),private,parameter :: uparfile=12

real(kind=REP)    :: delta_x,delta_y,xllcorner,yllcorner
integer(kind=ISP) :: N_celle,N,M,N_celle_read

real(kind=REP)    :: mean_s_max

real(kind=REP)    :: pt
integer(kind=ISP) :: imethod
real(kind=REP)    :: lambda
real(kind=RSP)    :: CC_threshold
integer(kind=ISP) :: ndcf
integer(kind=ISP) :: nchc
real(kind=REP)    :: A_threshold
real(kind=RSP)    :: ASk_threshold
real(kind=RSP)    :: kas
real(kind=RSP)    :: DN_threshold
real(kind=RSP)    :: local_slope_t
integer(kind=ISP) :: p_outflow_vo
real(kind=REP)    :: dr
integer(kind=ISP) :: bcc
real(kind=RSP)    :: cqm
real(kind=RSP)    :: cqg

real(kind=REP)    :: As_rf
real(kind=RSP)    :: Qsf_rf,w_rf
real(kind=RSP)    :: Wsf_rf,b1_rf,b2_rf
real(kind=RSP)    :: kSsf_rf,y1_rf,y2_rf
real(kind=RSP)    :: Qsi_rf
real(kind=REP)    :: As_cf
real(kind=RSP)    :: Qsf_cf,w_cf
real(kind=RSP)    :: Wsf_cf,b1_cf,b2_cf
real(kind=RSP)    :: kSsf_cf,y1_cf,y2_cf
real(kind=RSP)    :: Qsi_cf

real(kind=RSP)    :: Qsf_pf,w_pf
real(kind=RSP)    :: Wsf_pf,b1_pf,b2_pf
real(kind=RSP)    :: Kpsf_pf,d1_pf,d2_pf
real(kind=RSP)    :: Qsi_pf

public :: rparfile,wparfile

contains

!-----------------------------------------------------------------------------
!
! RPARFILE
!
!-----------------------------------------------------------------------------

subroutine rparfile(parfile)

character(len=*),intent(in) :: parfile
character(len=512) :: row
integer(kind=ISP) :: ic

open(uparfile,file=parfile,form="formatted",action="read")

! read all the lines of <parfile>

! structural parameters

read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"
read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"
read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, delta_x"
read(row,*,iostat=ic) delta_x
if (ic/=0) then
   stop "error when reading the parameter file, delta_x"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, delta_y"
read(row,*,iostat=ic) delta_y
if (ic/=0) then
   stop "error when reading the parameter file, delta_y"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, N"
read(row,*,iostat=ic) N
if (ic/=0) then
   stop "error when reading the parameter file, N"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, M"
read(row,*,iostat=ic) M
if (ic/=0) then
   stop "error when reading the parameter file, M"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, N_celle"
read(row,*,iostat=ic) N_celle
if (ic/=0) then
   stop "error when reading the parameter file, N_celle"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, xllcorner"
read(row,*,iostat=ic) xllcorner
if (ic/=0) then
   stop "error when reading the parameter file, xllcorner"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, yllcorner"
read(row,*,iostat=ic) yllcorner
if (ic/=0) then
   stop "error when reading the parameter file, yllcorner"
end if

! terrain analysis parameters

read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"
read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"
read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, pt"
read(row,*,iostat=ic) pt
if (ic/=0) then
   stop "error when reading the parameter file, pt"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, imethod"
read(row,*,iostat=ic) imethod
if (ic/=0) then
   stop "error when reading the parameter file, imethod"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, lambda"
read(row,*,iostat=ic) lambda
if (ic/=0) then
   stop "error when reading the parameter file, lambda"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, CC_threshold"
read(row,*,iostat=ic) CC_threshold
if (ic/=0) then
   stop "error when reading the parameter file, CC_threshold"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, ndcf"
read(row,*,iostat=ic) ndcf
if (ic/=0) then
   stop "error when reading the parameter file, ndcf"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, nchc"
read(row,*,iostat=ic) nchc
if (ic/=0) then
   stop "error when reading the parameter file, nchc"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, A_threshold"
read(row,*,iostat=ic) A_threshold
if (ic/=0) then
   stop "error when reading the parameter file, A_threshold"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, ASk_threshold"
read(row,*,iostat=ic) ASk_threshold
if (ic/=0) then
   stop "error when reading the parameter file, ASk_threshold"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, kas"
read(row,*,iostat=ic) kas
if (ic/=0) then
   stop "error when reading the parameter file, kas"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, DN_threshold"
read(row,*,iostat=ic) DN_threshold
if (ic/=0) then
   stop "error when reading the parameter file, DN_threshold"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, local_slope_t"
read(row,*,iostat=ic) local_slope_t
if (ic/=0) then
   stop "error when reading the parameter file, local_slope_t"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, p_outflow_vo"
read(row,*,iostat=ic) p_outflow_vo
if (ic/=0) then
   stop "error when reading the parameter file, p_outflow_vo"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, bcc"
read(row,*,iostat=ic) bcc
if (ic/=0) then
   stop "error when reading the parameter file, bcc"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, cqm"
read(row,*,iostat=ic) cqm
if (ic/=0) then
   stop "error when reading the parameter file, cqm"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, cqg"
read(row,*,iostat=ic) cqg
if (ic/=0) then
   stop "error when reading the parameter file, cqg"
end if

! rivulet network parameters

read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"
read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"
read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, dr"
read(row,*,iostat=ic) dr
if (ic/=0) then
   stop "error when reading the parameter file, dr"
end if

write(6,*)
if (mod(delta_x,dr).gt.epsilon(delta_x/dr)) then
   write(6,*) 'DEM resolution is not a multiple of the rivulet spacing!'
   write(6,*)
   stop
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, As_rf"
read(row,*,iostat=ic) As_rf
if (ic/=0) then
   stop "error when reading the parameter file, As_rf"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, Qsf_rf,w_rf"
read(row,*,iostat=ic) Qsf_rf,w_rf
if (ic/=0) then
   stop "error when reading the parameter file, Qsf_rf,w_rf"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, Wsf_rf,b1_rf,b2_rf"
read(row,*,iostat=ic) Wsf_rf,b1_rf,b2_rf
if (ic/=0) then
   stop "error when reading the parameter file, Wsf_rf,b1_rf,b2_rf"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, kSsf_rf,y1_rf,y2_rf"
read(row,*,iostat=ic) kSsf_rf,y1_rf,y2_rf
if (ic/=0) then
   stop "error when reading the parameter file, kSsf_rf,y1_rf,y2_rf"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, Qsi_rf"
read(row,*,iostat=ic) Qsi_rf
if (ic/=0) then
   stop "error when reading the parameter file, Qsi_rf"
end if

! channel network parameters

read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"
read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"
read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, As_cf"
read(row,*,iostat=ic) As_cf
if (ic/=0) then
   stop "error when reading the parameter file, As_cf"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, Qsf_cf,w_cf"
read(row,*,iostat=ic) Qsf_cf,w_cf
if (ic/=0) then
   stop "error when reading the parameter file, Qsf_cf,w_cf"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, Wsf_cf,b1_cf,b2_cf"
read(row,*,iostat=ic) Wsf_cf,b1_cf,b2_cf
if (ic/=0) then
   stop "error when reading the parameter file, Wsf_cf,b1_cf,b2_cf"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, kSsf_cf,y1_cf,y2_cf"
read(row,*,iostat=ic) kSsf_cf,y1_cf,y2_cf
if (ic/=0) then
   stop "error when reading the parameter file, kSsf_cf,y1_cf,y2_cf"
end if

call rrow(row,ic)
if (ic/=0) stop "error when reading the parameter file, Qsi_cf"
read(row,*,iostat=ic) Qsi_cf
if (ic/=0) then
   stop "error when reading the parameter file, Qsi_cf"
end if

read(uparfile,"(a)",iostat=ic) row ! read a comment line
if (ic/=0) stop "error when reading the parameter file"

close(uparfile)

end subroutine rparfile

!-----------------------------------------------------------------------------
!
! RROW
!
!-----------------------------------------------------------------------------

subroutine rrow(row,ic)

character(len=*),intent(out) :: row
integer(kind=ISP),intent(out) :: ic
integer(kind=ISP) i

do
   read(uparfile,"(a)",iostat=ic) row
   if (ic/=0) return
   i=index(row,"=")
   if (i>0) exit
end do
row=row(i+1:)

end subroutine rrow

!-----------------------------------------------------------------------------
!
! WPARFILE
!
!-----------------------------------------------------------------------------

subroutine wparfile(parfile)

character(len=*),intent(in) :: parfile

open(uparfile,file=parfile,form="formatted")
rewind(uparfile)

! structural parameters

write(uparfile,"(a)") &
"------------------------------------------------------------------------------"
write(uparfile,"(a)") "STRUCTURAL PARAMETERS"
write(uparfile,"(a)") &
"------------------------------------------------------------------------------"

write(uparfile,"(a,20x,f10.2)") &
"Grid spacing along the x-direction = ",delta_x

write(uparfile,"(a,20x,f10.2)") &
"Grid spacing along the y-direction = ",delta_y

write(uparfile,"(a,14x,i7)") &
"DEM rectangle size along the x-direction = ",N

write(uparfile,"(a,14x,i7)") &
"DEM rectangle size along the y-direction = ",M

write(uparfile,"(a,15x,i10)") &
"Number of cells within the catchment = ",N_celle

write(uparfile,"(a,22x,f20.8)") &
"X low left corner coordinate = ",xllcorner

write(uparfile,"(a,22x,f20.8)") &
"Y low left corner coordinate = ",yllcorner

! terrain analysis parameters

write(uparfile,"(a)") &
"------------------------------------------------------------------------------"
write(uparfile,"(a)") "TERRAIN ANALYSIS PARAMETERS"
write(uparfile,"(a)") &
"------------------------------------------------------------------------------"

write(uparfile,"(a,38x,1e10.3)") &
"Depit threshold slope = ",pt

write(uparfile,"(a,17x,i4)") &
"Drainage directions method (LAD:1,LTD:2) = ",imethod

write(uparfile,"(a,13x,5e10.3)") &
"Upstream deviation memory factor (CBM:0,PBM:1) = ",lambda

write(uparfile,"(a,4x,1e10.3)") &
"Threshold on the contour curvature (NDM:-1E10;DM:+1E10) = ",CC_threshold

write(uparfile,"(a,6x,i1)") &
"Nondispersive channel flow (0:not-required;1:required) = ",ndcf

write(uparfile,"(a,13x,i4)") &
"Channel initiation method (A:1,AS**k:2,ND:3) = ",nchc

write(uparfile,"(a,26x,1e16.9)") &
"Threshold on the support area (A) = ",A_threshold

write(uparfile,"(a,23x,1f10.2)") &
"Threshold on the AS**k function = ",ASk_threshold

write(uparfile,"(a,22x,1f10.2)") &
"Exponent k of the AS**k function = ",kas

write(uparfile,"(a,16x,1e10.3)") &
"Threshold on the normalized divergence (ND) = ",DN_threshold

write(uparfile,"(a,39x,1e10.3)") &
"Path threshold slope = ",local_slope_t

write(uparfile,"(a,4x,i1)") &
"Drainage direction of the outlet cell (if necessary...)  = ",p_outflow_vo

write(uparfile,"(a,19x,i1)") &
"Boundary channel constraction (No:0,Yes:1) =",bcc

write(uparfile,"(a,7x,f5.2)") &
"Coefficient for boundary channel elevation definition =",cqm

write(uparfile,"(a,12x,f5.2)") &
"Coefficient for outlet cell elevation definition =",cqg

! rivulet network parameters

write(uparfile,"(a)") &
"------------------------------------------------------------------------------"
write(uparfile,"(a)") &
"RIVULET NETWORK PARAMETERS (HYDRAULIC GEOMETRY OF THE SINGLE RIVULET)"
write(uparfile,"(a)") &
"------------------------------------------------------------------------------"

write(uparfile,"(a,30x,1f10.3)") &
"Rivulet spacing = ",dr

write(uparfile,"(a,18x,1e19.12)") &
"Reference drainage area (As_rf) = ",As_rf

write(uparfile,"(a,17x,1f10.3,10x,1f10.3)") &
"Flow discharge (Qsf_rf,w_rf) = ",Qsf_rf,w_rf

write(uparfile,"(a,5x,3f10.3)") &
"Water-surface width (Wsf_rf,b1_rf,b2_rf) = ",Wsf_rf,b1_rf,b2_rf

write(uparfile,"(a,1x,3f10.3)") &
"Resistance coefficient (kSsf_rf,y1_rf,y2_rf) = ",kSsf_rf,y1_rf,y2_rf

write(uparfile,"(a,14x,1f10.3)") &
"Initial flow discharge (Qsi_rf) = ",Qsi_rf

! channel network parameters

write(uparfile,"(a)") &
"------------------------------------------------------------------------------"
write(uparfile,"(a)") "CHANNEL NETWORK PARAMETERS"
write(uparfile,"(a)") &
"------------------------------------------------------------------------------"

write(uparfile,"(a,18x,1e19.12)") &
"Reference drainage area (As_cf) = ",As_cf

write(uparfile,"(a,17x,1f10.3,10x,1f10.3)") &
"Flow discharge (Qsf_cf,w_cf) = ",Qsf_cf,w_cf

write(uparfile,"(a,5x,3f10.3)") &
"Water-surface width (Wsf_cf,b1_cf,b2_cf) = ",Wsf_cf,b1_cf,b2_cf

write(uparfile,"(a,1x,3f10.3)") &
"Resistance coefficient (kSsf_cf,y1_cf,y2_cf) = ",kSsf_cf,y1_cf,y2_cf

write(uparfile,"(a,14x,1f10.3)") &
"Initial flow discharge (Qsi_cf) = ",Qsi_cf

write(uparfile,"(a)") & 
"------------------------------------------------------------------------------"

close(uparfile)

end subroutine wparfile

end module mpar
