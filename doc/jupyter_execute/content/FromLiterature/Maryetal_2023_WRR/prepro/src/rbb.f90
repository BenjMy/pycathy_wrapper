!-----------------------------------------------------------------------------
!
! RBB
!
! The program RBB reads the parameters contained in the binary files
! basin_b/basin_i and writes an DTM ascii file as output. This output
! file is called rdtm.dat or idtm.dat depending upon the parameter is
! real or integer, respectively.
!
!-----------------------------------------------------------------------------

program rbb

use mpar
use mbbio

implicit none

integer(kind=ISP) i,j,i_basin
integer(kind=ISP) :: ibb,ips,jd(0:9),ht

integer(kind=IEP) :: idtm_max_i,idtm_min_i,idtm_max,idtm_min
integer(kind=IEP) :: sum_idtm,sum_c
real(kind=REP) :: rdtm_max_i,rdtm_min_i,rdtm_max,rdtm_min
real(kind=REP) :: sum_rdtm,rmed_idtm,rmed_rdtm
real(kind=RSP) :: nodata
integer(kind=ISP) :: ic
character(len=32) :: dtmfile
character(len=24) :: fstr
character(len=24) :: nstr

real(kind=REP),allocatable :: rdtm(:,:)
integer(kind=ISP),allocatable :: idtm(:,:)
integer(kind=ISP),external :: parsel

! read <parfile>

call rparfile('hap.in')

allocate(rdtm(N,M))
allocate(idtm(N,M))

call load_dtm("basin_b","basin_i")

! select the output DTM file name

write(6,*)
write(6,'(1x,a)') 'Select the output DTM file name:'
write(6,'(1x,a)') '(Ctrl C or Enter to exit)'
write(6,*)
write(6,'(1x,a3,$)') '-> '
read(5,'(a)') dtmfile
if (dtmfile == " ") stop
open(80,file=dtmfile,status='unknown',iostat=ic)
do while (ic /= 0)
   write(6,'(1x,a)') 'Invalid selection: try again...'
   write(6,'(1x,a)') '(Ctrl C or Enter to exit)'
   write(6,*) 
   write(6,'(1x,a3,$)') '-> '
   read(5,'(a)') dtmfile
   if (dtmfile == " ") stop
   open(80,file=dtmfile,status='unknown',iostat=ic)
end do

! select the header type

write(6,*)
write(6,'(1x,a)') 'Select the header type:'
write(6,'(1x,a)') '0) None'
write(6,'(1x,a)') '1) ESRI ascii file'
write(6,'(1x,a)') '2) GRASS ascii file'
write(6,'(1x,a)') '(Ctrl C to exit)'
write(6,*)
write(6,'(1x,a3,$)') '-> '
read(5,*,iostat=ic) ht

! select the nodata threshold

write(6,*)
write(6,'(1x,a)') 'Select the nodata value:'
write(6,'(1x,a)') '(Ctrl C to exit)'
write(6,*)
write(6,'(1x,a3,$)') '-> '
read(5,*,iostat=ic) nodata
do while (ic /= 0)
   write(6,'(1x,a)') 'Invalid selection: try again...'
   write(6,'(1x,a)') '(Ctrl C to exit)'
   write(6,*)
   write(6,'(1x,a3,$)') '-> '
   read(5,*,iostat=ic) nodata
end do

! select the parameter

ibb=parsel()

if (ibb == 14 .or. ibb == 15 .or. ibb==50) then
!  two cases are possible:
!  1) HAP pointer system
!  2) Arc/Gis pointer system (Jenson and Domingue, 1988)
   write(6,*) 
   write(6,'(1x,a)') 'Select the pointer system:'
   write(6,'(1x,a)') '1) HAP system'
   write(6,'(1x,a)') '2) Arc/Gis system'
   write(6,'(1x,a)') '(Ctrl C to exit)'
   write(6,*)
   write(6,'(1x,a3,$)') '-> '
   read(5,*,iostat=ic) ips
   do while (ic /= 0)
      write(6,'(1x,a)') 'Invalid selection: try again...'
      write(6,'(1x,a)') '(Ctrl C to exit)'
      write(6,*)
      write(6,'(1x,a3,$)') '-> '
      read(5,*,iostat=ic) ips
   end do
   jd(0)=0 ! outlet cell 
   jd(1)=8
   jd(2)=16
   jd(3)=32
   jd(4)=4
   jd(5)=0
   jd(6)=64
   jd(7)=2
   jd(8)=1
   jd(9)=128
end if

! data creation

sum_c=0
idtm_max_i=-huge(idtm_max)
idtm_min_i=huge(idtm_min)
idtm_max=idtm_max_i
idtm_min=idtm_min_i
sum_idtm=0
rdtm_max_i=-huge(rdtm_max)
rdtm_min_i=huge(rdtm_min)
rdtm_max=rdtm_max_i
rdtm_min=rdtm_min_i
sum_rdtm=0.0_REP
do i=1,N
   do j=1,M
      i_basin=(i-1)*M+j
      if (dtm_index_pr(i_basin) == 0) then
         select case (ibb)
         case (14:16,31,34:36,38:40,50)
            idtm(i,j)=int(nodata)
         case default
            rdtm(i,j)=nodata
         end select
      else
         sum_c=sum_c+1
         select case (ibb)
         case (14:16,31,34:36,38:40,50)
            select case (ibb)
            case (14)
               select case (ips)
               case (1)
                  idtm(i,j)=dtm_p_outflow_1(i_basin)
               case (2)
                  idtm(i,j)=jd(dtm_p_outflow_1(i_basin))
               case default
                  stop "unexpected case"
               end select
            case (15)
               select case (ips)
               case (1)
                  idtm(i,j)=dtm_p_outflow_2(i_basin)
               case (2)
                  idtm(i,j)=jd(dtm_p_outflow_2(i_basin))
               case default
                  stop "unexpected case"
               end select
            case (16)
               idtm(i,j)=dtm_dmID(i_basin)
            case (31)
               idtm(i,j)=dtm_hcID(i_basin)
            case (34)
               idtm(i,j)=dtm_lakes_map(i_basin)
            case (35)
               idtm(i,j)=dtm_zone(i_basin)
            case (36)
               idtm(i,j)=dtm_q_output(i_basin)
            case (38)
               idtm(i,j)=dtm_sso_jun(i_basin)
            case (39)
               idtm(i,j)=dtm_sso(i_basin)
            case (40)
               idtm(i,j)=dtm_hso(i_basin)
            case (50)
               select case (ips)
               case (1)
                  idtm(i,j)=dtm_p_outflow_1(i_basin)*int(dtm_w_1(i_basin))+dtm_p_outflow_2(i_basin)*int(dtm_w_2(i_basin))
               case (2)
                  idtm(i,j)=jd(dtm_p_outflow_1(i_basin)*int(dtm_w_1(i_basin))+dtm_p_outflow_2(i_basin)*int(dtm_w_2(i_basin)))
               case default
                  stop "unexpected case"
               end select
            case default
               stop "unexpected case"
            end select
            sum_idtm=sum_idtm+idtm(i,j)
            if (idtm(i,j).gt.idtm_max) idtm_max=idtm(i,j)
            if (idtm(i,j).lt.idtm_min) idtm_min=idtm(i,j)
         case default
            select case (ibb)
            case (1)
               rdtm(i,j)=dtm_vh(i_basin)
            case (2)
               rdtm(i,j)=dtm_rl_min(i_basin)
            case (3)
               rdtm(i,j)=dtm_LAI(i_basin)
            case (4)
               rdtm(i,j)=dtm_vsc(i_basin)
            case (5)
               rdtm(i,j)=dtm_theta_s(i_basin)
            case (6)
               rdtm(i,j)=dtm_eta(i_basin)
            case (7)
               rdtm(i,j)=dtm_psi_s(i_basin)
            case (8)
               rdtm(i,j)=dtm_K_s(i_basin)
            case (9)
               rdtm(i,j)=dtm_Z_low(i_basin)
            case (10)
               rdtm(i,j)=dtm_Z_up(i_basin)
            case (11)
               rdtm(i,j)=dtm_theta_r(i_basin)
            case (12)
               rdtm(i,j)=dtm_rl_max(i_basin)
            case (13)
               rdtm(i,j)=dtm_quota(i_basin)
            case (17)
               rdtm(i,j)=dtm_A_inflow(i_basin)
            case (18)
               rdtm(i,j)=dtm_w_1(i_basin)
            case (19)
               rdtm(i,j)=dtm_w_2(i_basin)
            case (20)
               rdtm(i,j)=dtm_sumdev_num(i_basin)
            case (21)
               rdtm(i,j)=dtm_local_slope_1(i_basin)
            case (22)
               rdtm(i,j)=dtm_local_slope_2(i_basin)
            case (23)
               rdtm(i,j)=dtm_Ws1_sf_1(i_basin)
            case (24)
               rdtm(i,j)=dtm_Ws1_sf_2(i_basin)
            case (25)
               rdtm(i,j)=dtm_b1_sf(i_basin)
            case (26)
               rdtm(i,j)=dtm_kSs1_sf_1(i_basin)
            case (27)
               rdtm(i,j)=dtm_kSs1_sf_2(i_basin)
            case (28)
               rdtm(i,j)=dtm_y1_sf(i_basin)
            case (29)
               rdtm(i,j)=dtm_ASk(i_basin)
            case (30)
               rdtm(i,j)=dtm_DN(i_basin)
            case (32)
               rdtm(i,j)=dtm_epl_1(i_basin)
            case (33)
               rdtm(i,j)=dtm_epl_2(i_basin)
            case (37)
               rdtm(i,j)=dtm_nrc(i_basin)
            case (51)
               rdtm(i,j)=dtm_w_1(i_basin)+dtm_w_2(i_basin)
            case (52)
               rdtm(i,j)=dtm_local_slope_1(i_basin)*dtm_w_1(i_basin)+dtm_local_slope_2(i_basin)*dtm_w_2(i_basin)
            case (53)
               rdtm(i,j)=dtm_Ws1_sf_1(i_basin)*dtm_w_1(i_basin)+dtm_Ws1_sf_2(i_basin)*dtm_w_2(i_basin)
            case (54)
               rdtm(i,j)=dtm_kSs1_sf_1(i_basin)*dtm_w_1(i_basin)+dtm_kSs1_sf_2(i_basin)*dtm_w_2(i_basin)
            case (55)
               rdtm(i,j)=dtm_epl_1(i_basin)*dtm_w_1(i_basin)+dtm_epl_2(i_basin)*dtm_w_2(i_basin)         
			case default
               stop "unexpected case"
            end select
            sum_rdtm=sum_rdtm+rdtm(i,j)
            if (rdtm(i,j).gt.rdtm_max) rdtm_max=rdtm(i,j)
            if (rdtm(i,j).lt.rdtm_min) rdtm_min=rdtm(i,j)
         end select
      end if
   end do
end do

! integer parameter?

if (idtm_max >= idtm_min) then
   write(6,'(1x,a,i6)') 'min value =',idtm_min
   write(6,'(1x,a,i6)') 'max value =',idtm_max
   if (idtm_max.eq.idtm_max_i) then
      write(6,'(1x,a)') 'initial value unchanged'
      stop
   end if
   if (idtm_min.eq.idtm_min_i) then
      write(6,'(1x,a)') 'initial value unchanged'
      stop
   end if
   rmed_idtm=real(sum_idtm)/sum_c
   write(6,'(1x,a,i6)') 'number of cells =',sum_c
   write(6,'(1x,a,f13.6)') 'mean value =',rmed_idtm
   write(6,*)
   write(6,'(1x,a)') 'writing the output file...'
   write(6,*)

! write an integer output file

   write(fstr,'(i10)') N
   idtm_min=min(idtm_min,int(nodata))  
   idtm_max=max(idtm_max,abs(int(nodata)))  
   if (idtm_max.gt.0) then 
      if (idtm_min.lt.0) then
         write(nstr,'(i10)') int(log10(real(idtm_max)))+3
      else 
         write(nstr,'(i10)') int(log10(real(idtm_max)))+2
      end if
   else if (idtm_max.eq.0) then
      write(nstr,'(i10)') 2
   else
      write(6,*) 'unexpected case!'
      stop
   end if
   fstr='('//trim(adjustl(fstr))//'i'//trim(adjustl(nstr))//')'
   
! ESRI ascii file header writing
   if (ht.eq.1) then
      write(80,'(a,8x,i5)') 'ncols',N
      write(80,'(a,9x,i5)') 'nrow',M
      write(80,'(a,4x,f20.8)') 'xllcorner',xllcorner
      write(80,'(a,4x,f20.8)') 'yllcorner',yllcorner
      write(80,'(a,7x,f6.2)') 'cellsize',delta_x
      write(80,'(a,1x,i5)') 'NODATA_value',-9999
   end if

! GRASS ascii file header writing  
   if (ht.eq.2) then      
      write(80,'(a,1x,i5)') 'north:',0
      write(80,'(a,1x,i5)') 'south:',0
      write(80,'(a,2x,i5)') 'east:',0
      write(80,'(a,2x,i5)') 'west:',0
      write(80,'(a,2x,i5)') 'rows:',M
      write(80,'(a,2x,i5)') 'cols:',N
   end if
   
   do j=M,1,-1
      write(80,fstr) (idtm(i,j),i=1,N)
   end do
end if

! real parameter?

if (rdtm_max >= rdtm_min) then
   select case (ibb)
   case (17)
      write(6,'(1x,a,e19.12)') 'min value =',rdtm_min
      write(6,'(1x,a,e19.12)') 'max value =',rdtm_max
   case default
      write(6,'(1x,a,e13.6)') 'min value =',rdtm_min
      write(6,'(1x,a,e13.6)') 'max value =',rdtm_max
   end select
   if (abs(rdtm_max-rdtm_max_i).le.epsilon(rdtm_max)) then
      write(6,*) 'initial value unchanged'
      stop
   end if
   if (abs(rdtm_min-rdtm_min_i).le.epsilon(rdtm_min)) then
      write(6,*) 'initial value unchanged'
      stop
   end if
   rmed_rdtm=real(sum_rdtm)/sum_c
   select case (ibb)
   case (17)
      write(6,'(1x,a,i6)') 'number of cells =',sum_c
      write(6,'(1x,a,e19.12)') 'mean value =',rmed_rdtm
      write(6,*)
      write(6,'(1x,a)') 'writing the output file...'
   case default
      write(6,'(a,i6)') 'number of cells =',sum_c
      write(6,'(a,e13.6)') 'mean value =',rmed_rdtm
      write(6,*)
      write(6,'(1x,a)') 'writing the output file...'
   end select
   write(6,*)

! write a real output file

   write(fstr,'(i10)') N
   rdtm_min=min(rdtm_min,nodata)  

   select case (ibb)
   case (17)
      if (rdtm_min.lt.0.0) then
!        fstr='('//trim(adjustl(fstr))//'e19.12'//')'       
        fstr='('//trim(adjustl(fstr))//'f15.2'//')'
      else
!        fstr='('//trim(adjustl(fstr))//'e20.12'//')'
        fstr='('//trim(adjustl(fstr))//'f14.2'//')'
      end if
   case default
      if (rdtm_min.lt.0.0) then
!         fstr='('//trim(adjustl(fstr))//'e14.6'//')'
		 fstr='('//trim(adjustl(fstr))//'e20.12'//')'
      else
!         fstr='('//trim(adjustl(fstr))//'e13.6'//')'
		 fstr='('//trim(adjustl(fstr))//'e21.12'//')'
      end if
   end select

! ESRI ascii file header writing
   if (ht.eq.1) then
      write(80,'(a,8x,i5)') 'ncols',N
      write(80,'(a,9x,i5)') 'nrow',M
      write(80,'(a,4x,f20.8)') 'xllcorner',xllcorner
      write(80,'(a,4x,f20.8)') 'yllcorner',yllcorner
      write(80,'(a,7x,f6.2)') 'cellsize',delta_x
      write(80,'(a,1x,i5)') 'NODATA_value',-9999
   end if

! GRASS ascii file header writing  
   if (ht.eq.2) then      
      write(80,*) 'north: 0'
      write(80,*) 'south: 0'
      write(80,*) 'east: 0'
      write(80,*) 'west: 0'
      write(80,*) 'rows:',M
      write(80,*) 'cols:',N
   end if
   
!    
   do j=M,1,-1
      write(80,fstr) (rdtm(i,j),i=1,N)
   end do
end if

! close the file

close(80)

end program rbb

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
if (ibb <= 0) ic=1
if (ibb > 40 .and. ibb < 50) ic=1
if (ibb > 55) ic=1
do while (ic /= 0) 
   write(6,*) 'Invalid selection: try again...'
   call maschera()
   read(5,*,iostat=ic) ibb
   if (ibb <= 0) ic=1
   if (ibb > 40 .and. ibb < 50) ic=1
   if (ibb > 55) ic=1
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
write(6,'(1x,a)') 'Select the parameter to be stored:'
write(6,'(1x,a)') '(Ctrl C to exit)'
write(6,*)
write(6,'(1x,a)') '-------------------------------------------------'
write(6,'(1x,a)') ' 1) vh                         21) local_slope_1 '
write(6,'(1x,a)') ' 2) rl_min                     22) local_slope_2 '
write(6,'(1x,a)') ' 3) LAI                        23) Ws1_sf_1      '
write(6,'(1x,a)') ' 4) vsc                        24) Ws1_sf_2      '
write(6,'(1x,a)') ' 5) theta_s                    25) b1_sf         '
write(6,'(1x,a)') ' 6) eta                        26) kSs1_sf_1     '
write(6,'(1x,a)') ' 7) psi_s                      27) kSs1_sf_2     '
write(6,'(1x,a)') ' 8) K_s                        28) y1_sf         '
write(6,'(1x,a)') ' 9) Z_low                      29) ASk           '
write(6,'(1x,a)') '10) Z_up                       30) DN            '
write(6,'(1x,a)') '11) theta_r                    31) hcID          '
write(6,'(1x,a)') '12) rl_max                     32) epl_1         '
write(6,'(1x,a)') '13) quota                      33) epl_2         '
write(6,'(1x,a)') '14) p_outflow_1                34) lakes_map     '
write(6,'(1x,a)') '15) p_outflow_2                35) zone          '
write(6,'(1x,a)') '16) dmID                       36) q_output      '
write(6,'(1x,a)') '17) A_inflow                   37) nrc           '
write(6,'(1x,a)') '18) w_1                        38) sso_jun       '
write(6,'(1x,a)') '19) w_2                        39) sso           '
write(6,'(1x,a)') '20) sumdev_num                 40) hso           '
write(6,'(1x,a)') '-------------------------------------------------'
write(6,*)

if (ndcf == 1) then
   write(6,'(1x,a)') 'If you want a merged output, you can also select:'
   write(6,'(1x,a)') '-------------------------------------------------'
   write(6,'(1x,a)') '50) p_outflow                                    '
   write(6,'(1x,a)') '51) w                                            '
   write(6,'(1x,a)') '52) local_slope                                  '
   write(6,'(1x,a)') '53) Ws1_sf                                       '
   write(6,'(1x,a)') '54) kSs1_sf                                      '
   write(6,'(1x,a)') '55) epl                                          '
   write(6,'(1x,a)') '-------------------------------------------------'
   write(6,*)
   write(6,'(1x,a3,$)') '-> '
end if

end subroutine

end function
