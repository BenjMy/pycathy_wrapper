!-----------------------------------------------------------------------------
!
! CBB
!
! This program assigns the values read in the ascii DTM file <dtmfile>
! to the ibb-th position of each record of basin_b.
! The following reference system is considered:
!
!              ^
!              |
!              |
!    j=[1,M],y |
!              |
!              |
!             -|--------------->
!                  i=[1,N],x
!
!-----------------------------------------------------------------------------

program cbb

use mpar
use mbbio

implicit none

integer(kind=ISP) i,j,i_basin,k_rec
real(kind=RSP) nodata
integer(kind=ISP) :: ic
integer(kind=ISP) :: ibb
character(len=32) dtmfile

integer(kind=ISP),allocatable,dimension(:) :: iwa
real(kind=RSP),allocatable,dimension(:)    :: rwa_sp
real(kind=REP),allocatable,dimension(:)  :: rwa_ep
integer(kind=ISP),external :: parsel

! select the input DTM file name
! open the ascii catchment DEM file <dtmfile>

write(6,*)
write(6,'(a)') 'Select the input DTM file name:'
write(6,'(a)') '(Ctrl C or Enter to exit)'
write(6,*)
write(6,'(a3,$)') '-> '
read(5,'(a)') dtmfile
if (dtmfile == " ") stop
open(80,file=dtmfile,status='old',iostat=ic)
do while (ic /= 0)
   write(6,'(a)') 'File not found: select another name'
   write(6,'(a)') '(Ctrl C or Enter to exit)'
   write(6,*) 
   write(6,'(a3,$)') '-> '
   read(5,'(a)') dtmfile
   if (dtmfile == " ") stop
   open(80,file=dtmfile,status='old',iostat=ic)
end do

! select the nodata threshold

write(6,*)
write(6,'(a)') 'Select the notada value (only larger values are considered):'
write(6,'(a)') '(Ctrl C to exit)'
write(6,*)
write(6,'(a3,$)') '-> '
read(5,*,iostat=ic) nodata
do while (ic /= 0)
   write(6,'(a)') 'Invalid selection: try again...'
   write(6,'(a)') '(Ctrl C to exit)'
   write(6,*)
   write(6,'(a3,$)') '-> '
   read(5,*,iostat=ic) nodata
end do

! select the parameter

ibb=parsel()

! read <parfile>

call rparfile('hap.in')

select case (ibb)
case(14:16,31,34:36)
   allocate(iwa(N))
case(17)
   allocate(rwa_ep(N))
case default
   allocate(rwa_sp(N))
end select

call load_dtm("basin_b","basin_i")

! write the binary files basin_b

k_rec=0
rewind(80)
do j=M,1,-1
   select case(ibb)
   case (14:16,31,34:36)
      read(80,*,iostat=ic) (iwa(i),i=1,N)
   case (17)
      read(80,*,iostat=ic) (rwa_ep(i),i=1,N)
   case default
      read(80,*,iostat=ic) (rwa_sp(i),i=1,N)
   end select
   select case (ic)
   case (:-1)
      write(6,*) "error when reading the file ",dtmfile
      stop
   case (0)
      do i=1,N
         select case(ibb)
         case (14:16,31,34:36)
            if (iwa(i).gt.int(nodata)) then
               k_rec=k_rec+1
               i_basin=(i-1)*M+j
               select case(ibb)
               case (14)
                  call set_p_outflow_1(i_basin,iwa(i))
               case (15)
                  call set_p_outflow_2(i_basin,iwa(i))
               case (16)
                  call set_dmID(i_basin,iwa(i))
               case (31)
                  call set_hcID(i_basin,iwa(i))
               case (34)
                  call set_lakes_map(i_basin,iwa(i))
               case (35)
                  call set_zone(i_basin,iwa(i))
               case (36)
                  call set_q_output(i_basin,iwa(i))
               case default
                  write(6,*) "unexpected case"
                  stop
               end select
            end if
         case (17)
            if (rwa_ep(i).gt.nodata) then
               k_rec=k_rec+1
               i_basin=(i-1)*M+j
               select case(ibb)
               case (17)
                  call set_A_inflow(i_basin,rwa_ep(i))
               case (20)
                  call set_sumdev_num(i_basin,rwa_ep(i))
               case default
                  write(6,*) "unexpected case"
                  stop
               end select
            end if
         case default
            if (rwa_sp(i).gt.nodata) then
               k_rec=k_rec+1
               i_basin=(i-1)*M+j
               select case(ibb)
               case (1)      
                  call set_vh(i_basin,rwa_sp(i))
               case (2)      
                  call set_rl_min(i_basin,rwa_sp(i))
               case (3)      
                  call set_LAI(i_basin,rwa_sp(i))
               case (4)      
                  call set_vsc(i_basin,rwa_sp(i))
               case (5)      
                  call set_theta_s(i_basin,rwa_sp(i))
               case (6)      
                  call set_eta(i_basin,rwa_sp(i))
               case (7)      
                  call set_psi_s(i_basin,rwa_sp(i))
               case (8)      
                  call set_K_s(i_basin,rwa_sp(i))
               case (9)      
                  call set_Z_low(i_basin,rwa_sp(i))
               case (10)      
                  call set_Z_up(i_basin,rwa_sp(i))
               case (11)      
                  call set_theta_r(i_basin,rwa_sp(i))
               case (12)      
                  call set_rl_max(i_basin,rwa_sp(i))
               case (13)      
                  call set_quota(i_basin,rwa_ep(i))
               case (18)
                  call set_w_1(i_basin,rwa_sp(i))
               case (19)
                  call set_w_2(i_basin,rwa_sp(i))               
               case (21)
                  call set_local_slope_1(i_basin,rwa_sp(i))
               case (22)
                  call set_local_slope_2(i_basin,rwa_sp(i))
               case (23)
                  call set_Ws1_sf_1(i_basin,rwa_sp(i))
               case (24)
                  call set_Ws1_sf_2(i_basin,rwa_sp(i))
               case (25)
                  call set_b1_sf(i_basin,rwa_sp(i))
               case (26)
                  call set_kSs1_sf_1(i_basin,rwa_sp(i))
               case (27)
                  call set_kSs1_sf_2(i_basin,rwa_sp(i))
               case (28)
                  call set_y1_sf(i_basin,rwa_sp(i))
               case (29)
                  call set_ASk(i_basin,rwa_sp(i))
               case (30)
                  call set_DN(i_basin,rwa_sp(i))
               case (32)
                  call set_epl_1(i_basin,rwa_sp(i))
               case (33)
                  call set_epl_2(i_basin,rwa_sp(i))
               case (37)
                  call set_nrc(i_basin,rwa_sp(i))
               case default
                  write(6,*) "unexpected case"
                  stop
               end select
            end if
         end select
      end do
   case (1:)
      write(6,*) "insufficient data in the file ",dtmfile
      stop
   end select
end do
write(6,'(i6,a16)') k_rec,' processed cells'
write(6,*) 

! close the file

close(80)
call save_dtm("basin_b")

end program cbb

!-----------------------------------------------------------------------------
!
! PARSEL
!
! The function PARSEL perform the parameter selction task. 
!
!-----------------------------------------------------------------------------

function parsel()

use mpar
implicit none

integer(kind=ISP) :: parsel
integer(kind=ISP) :: ibb
integer(kind=ISP) :: ic

call maschera()
read(5,*,iostat=ic) ibb
if (ibb <= 0 .or. ibb > 36) ic=1
do while (ic /= 0) 
   write(6,*) 'Invalid selection: try again...'
   call maschera()
   read(5,*,iostat=ic) ibb
   if (ibb <= 0 .or. ibb > 36) ic=1
end do
write(6,*)

parsel=ibb

contains

!-----------------------------------------------------------------------------
!
! MASCHERA
!
! The subroutine MASCHERA provides the mask in the parameter selction task. 
!
!-----------------------------------------------------------------------------

subroutine maschera()
write(6,*) 
write(6,'(a)') 'Select the parameter to be stored:'
write(6,'(a)') '(Ctrl C to exit)'
write(6,*) 
write(6,'(a)') '------------------------------------------------'
write(6,'(a)') ' 1) vh                         21) local_slope_1'
write(6,'(a)') ' 2) rl_min                     22) local_slope_2'
write(6,'(a)') ' 3) LAI                        23) Ws1_sf_1     '
write(6,'(a)') ' 4) vsc                        24) Ws1_sf_2     '
write(6,'(a)') ' 5) theta_s                    25) b1_sf        '
write(6,'(a)') ' 6) eta                        26) kSs1_sf      '
write(6,'(a)') ' 7) psi_s                      27) kSs1_sf      '
write(6,'(a)') ' 8) K_s                        28) y1_sf        '
write(6,'(a)') ' 9) Z_low                      29) ASk          '
write(6,'(a)') '10) Z_up                       30) DN           '
write(6,'(a)') '11) theta_r                    31) hcID         '
write(6,'(a)') '12) rl_max                     32) epl_1        '
write(6,'(a)') '13) quota                      33) epl_2        '
write(6,'(a)') '14) p_outflow_1                34) lakes_map    '
write(6,'(a)') '15) p_outflow_2                35) zone         '
write(6,'(a)') '16) dmID                       36) q_output     '
write(6,'(a)') '17) A_inflow                   37) nrc          '
write(6,'(a)') '18) w_1                                         '
write(6,'(a)') '19) w_2                                         '
write(6,'(a)') '20) sumdev_num                                  '
write(6,'(a)') '------------------------------------------------'
write(6,*)
write(6,'(a3,$)') '-> '
end subroutine

end function
