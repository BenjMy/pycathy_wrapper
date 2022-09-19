!-----------------------------------------------------------------------------
!
! MRBB_SR
!
! The program MRBB (Multiple RBB) reads the parameters contained in the binary 
! files basin_b/basin_i and writes many DTM ascii files as output.
!
!-----------------------------------------------------------------------------

subroutine mrbb_sr

use mpar
use mbbio

implicit none

integer(kind=ISP) :: ips,ht
integer(kind=ISP) :: ic

real(kind=RSP) :: nodata

call rparfile('hap.in')
call load_dtm("basin_b","basin_i")

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

! two cases are possible:
! 1) HAP pointer system
! 2) Arc/Gis pointer system (Jenson and Domingue, 1988)
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

!call rbb(dtmfile,ibb,ht,nodata,ips)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dem file'
write(6,*)
call rbb('dem',13,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'lakes_map file'
write(6,*)
call rbb('lakes_map',34,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'zone file'
write(6,*)
call rbb('zone',35,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_w_1 file'
write(6,*)
call rbb('dtm_w_1',18,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_w_2 file'
write(6,*)
call rbb("dtm_w_2",19,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_p_outflow_1 file'
write(6,*)
call rbb('dtm_p_outflow_1',14,ht,nodata,ips)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_p_outflow_2 file'
write(6,*)
call rbb('dtm_p_outflow_2',15,ht,nodata,ips)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'A_inflow file'
write(6,*)
call rbb('dtm_A_inflow',17,ht,nodata,ips)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_local_slope_1 file'
write(6,*)
call rbb('dtm_local_slope_1',21,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_local_slope_2 file'
write(6,*)
call rbb('dtm_local_slope_2',22,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_epl_1 file'
write(6,*)
call rbb('dtm_epl_1',32,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_epl_2 file'
write(6,*)
call rbb('dtm_epl_2',33,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_kSs1_sf_1 file'
write(6,*)
call rbb('dtm_kSs1_sf_1',26,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_kSs1_sf_2 file'
write(6,*)
call rbb('dtm_kSs1_sf_2',27,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_Ws1_sf file'
write(6,*)
call rbb('dtm_Ws1_sf_1',23,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_Ws1_sf_2 file'
write(6,*)
call rbb('dtm_Ws1_sf_2',24,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_b1_sf file'
write(6,*)
call rbb('dtm_b1_sf',25,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_y1_sf file'
write(6,*)
call rbb('dtm_y1_sf',28,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_hcID file'
write(6,*)
call rbb('dtm_hcID',31,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_q_output file'
write(6,*)
call rbb('dtm_q_output',36,ht,nodata,0)

write(6,'(1x,a)') '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(6,*) 
write(6,'(1x,a)') 'dtm_nrc file'
write(6,*)
call rbb('dtm_nrc',37,ht,nodata,0)

call close_dtm()

end subroutine mrbb_sr

!-----------------------------------------------------------------------------
!
! RBB
!
! The subroutine RBB...
!
!-----------------------------------------------------------------------------

subroutine rbb(dtmfile,ibb,ht,nodata,ips)

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
character(len=*) :: dtmfile
character(len=24) :: fstr
character(len=24) :: nstr

real(kind=REP),allocatable :: rdtm(:,:)
integer(kind=ISP),allocatable :: idtm(:,:)

allocate(rdtm(N,M))
allocate(idtm(N,M))

open(80,file=dtmfile,status='unknown',iostat=ic)

if (ibb == 14.or. ibb == 15 ) then
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
      if (dtm_index_pr(i_basin) == 0)  then
         select case (ibb)
         case (14:16,31,34:36)
            idtm(i,j)=int(nodata)
         case default
            rdtm(i,j)=nodata
         end select
      else
         sum_c=sum_c+1
         select case (ibb)
         case (14:16,31,34:36)
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

! write a integer output file

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
      write(80,"(a,8x,i5)") 'ncols',N
      write(80,"(a,9x,i5)") 'nrow',M
      write(80,"(a,1x,f20.8)") 'xllcorner',xllcorner
      write(80,"(a,1x,f20.8)") 'yllcorner',yllcorner
      write(80,"(a,8x,f6.2)") 'cellsize',delta_x
      write(80,"(a,4x,i5)") 'NODATA_value',-9999  
   end if

! GRASS ascii file header writing  
   if (ht.eq.2) then      
      write(80,'(a,1x,i5)') 'north:',0
      write(80,'(a,1x,f20.8)') 'south:',yllcorner
      write(80,'(a,2x,i5)') 'east:',0
      write(80,'(a,2x,f20.8)') 'west:',xllcorner
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
      write(6,'(1x,a)') 'initial value unchanged'
      stop
   end if
   if (abs(rdtm_min-rdtm_min_i).le.epsilon(rdtm_min)) then
      write(6,'(1x,a)') 'initial value unchanged'
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
      write(6,'(1x,a,i6)') 'number of cells =',sum_c
      write(6,'(1x,a,e13.6)') 'mean value =',rmed_rdtm
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
      write(80,"(a,8x,i5)") 'ncols',N
      write(80,"(a,9x,i5)") 'nrow',M
      write(80,"(a,1x,f20.8)") 'xllcorner',xllcorner
      write(80,"(a,1x,f20.8)") 'yllcorner',yllcorner
      write(80,"(a,8x,f6.2)") 'cellsize',delta_x
      write(80,"(a,4x,i5)") 'NODATA_value',-9999  
   end if

! GRASS ascii file header writing  
   if (ht.eq.2) then      
      write(80,'(a,1x,i5)') 'north:',0
      write(80,'(a,1x,f20.8)') 'south:',yllcorner
      write(80,'(a,2x,i5)') 'east:',0
      write(80,'(a,2x,f20.8)') 'west:',xllcorner
      write(80,'(a,2x,i5)') 'rows:',M
      write(80,'(a,2x,i5)') 'cols:',N
   end if
   
!    
   do j=M,1,-1
      write(80,fstr) (rdtm(i,j),i=1,N)
   end do
end if

! close the file

close(80)

end subroutine rbb

