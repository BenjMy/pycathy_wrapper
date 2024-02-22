!-----------------------------------------------------------------------------
!
! WBB_SR
!
! This program defines the structure of the binary files basin_b/basin_i
! using the ascii DEM file <demfile>, which contains the catchment cell 
! elevations. The following reference system is considered:
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

subroutine wbb_sr

use mpar
use mbbio

implicit none

integer(kind=RSP) i,j,n_rec,i_basin,i_rec
real(kind=REP) nodata
integer(kind=ISP) :: ic
character(len=32) demfile

real(kind=REP),allocatable :: rwa(:)


integer(kind=ISP) ii,jj,ii_basin,last_i_basin,flag_gronda,cnt_celle,cnt_b_cell

real(kind=REP) quota_cij,quota_ciijj,quota_gronda,quota_min,quota_chiusura
real(kind=REP) max_q_depit,cqm_new

! set the DEM file name and nodata threshold
! (catchment cell are those with DEM elevation greater that nodata)

write(6,*)
write(6,*) 'searching the dtm_13.val input file...' 
demfile='dtm_13.val'
nodata=-9999.0_REP
write(6,*) 'assigned nodata value =',nodata 

! read <parfile>

call rparfile('hap.in')

allocate(rwa(N))

! open the ascii catchment DEM file <demfile>

open(73,file=demfile,status='old')

call catchment_cells(nodata,rwa,n_rec,quota_min)

call init_dtm(n_rec)

! write the binary files basin_b/basin_i

i_rec= 0
rewind 73
do j=M,1,-1
   read(73,*,iostat=ic) (rwa(i),i=1,N)
   select case(ic)
   case(:-1)
      write(6,*) "error when reading the file ",demfile
      stop
   case(0)
      do i=1,N
         if (rwa(i).gt.nodata) then
            i_rec=i_rec+1
            i_basin=(i-1)*M+j
            call set_dtm_index(i_basin,i_rec)
            call set_quota(i_basin,rwa(i))
         end if
      end do
   case(1:)
      write(6,*) "insufficient data in the file ",demfile
      stop
   end select
end do
N_celle=n_rec
write(6,*) 'number of processed cells =',n_rec
write(6,*)

if (bcc.eq.0) goto 2000

! procedure for boundary channel constraction

write(6,*) 'boundary channel constraction...'
write(6,*)

quota_gronda=quota_min*cqm
quota_chiusura=quota_gronda*cqg

cnt_b_cell=0
do j=M,1,-1
   do i=1,N
      i_basin=(i-1)*M+j
	  quota_cij=dtm_quota(i_basin)
      if (quota_cij < 0.0_REP) cycle
      flag_gronda=0
      do ii=i-1,i+1
	     if (flag_gronda.eq.1) cycle 
         do jj=j-1,j+1              
		    if (flag_gronda.eq.1) cycle    
            if (ii.eq.0.or.ii.eq.N+1) flag_gronda=1            
            if (jj.eq.0.or.jj.eq.M+1) flag_gronda=1           
            if (flag_gronda.eq.0) then
 		       ii_basin=(ii-1)*M+jj
               quota_ciijj=dtm_quota(ii_basin)
               if (quota_ciijj < 0.0_REP) flag_gronda=1
            end if
			if (flag_gronda.eq.1) then			   
			   cnt_b_cell=cnt_b_cell+1
			   last_i_basin=i_basin
               call set_quota(last_i_basin,quota_gronda)
			end if                  
	     end do
      end do
   end do
end do

call set_quota(last_i_basin,quota_chiusura)

write(6,*) 'boundary channel constraction terminated'
write(6,*)

write(6,*) 'minimum cell elevation (m) =',quota_min
write(6,*) 'boundary channel elevation (m) =',quota_gronda
write(6,*) 'outlet cell elevation (m) =',quota_chiusura
write(6,*) 'number of boundary channel cells =',cnt_b_cell
write(6,*)

max_q_depit=quota_gronda+(cnt_b_cell*delta_x*sqrt(2.0)*pt)
cqm_new=1-(cnt_b_cell*delta_x*sqrt(2.0)*pt)/quota_min

if (max_q_depit.ge.quota_min) then
   write(6,*) '********************************************************************************'
   write(6,*) 'WARNING:'
   write(6,*)
   write(6,*) 'Probably after the depitting, the highest boundary channel cell has elevation'
   write(6,*) 'great then the lowest dem cell.'
   write(6,*) 'A new value of coefficient for boundary channel elevation definition'
   write(6,*) 'is reccomended less then', cqm_new
   write(6,*)
   write(6,*) '********************************************************************************'
   stop
end if

2000 continue

! close the files

close(73)

call save_dtm("basin_b","basin_i")

end subroutine wbb_sr


subroutine catchment_cells(nodata,rwa,n_rec,quota_min)
use mpar
implicit none
integer(kind=ISP) ::  n_rec
real(kind=REP) :: nodata
real(kind=REP) :: quota_min,quota_min_tmp
real(kind=REP),dimension(n) :: rwa
integer(kind=ISP) :: i,ic
rewind 73
n_rec= 0
quota_min=8844.43 !Mount Everest elevation
do
   read(73,*,iostat=ic) (rwa(i),i=1,N)
   if (ic /= 0) exit
   do i=1,N
      if (rwa(i).gt.nodata) then
	     n_rec=n_rec+1
      else
 	     rwa(i)=quota_min  !fittizio per eliminare i nodata per 
 	                       !il calcolo della quota minuma
      end if
   end do
   quota_min_tmp=minval(rwa)
   quota_min=dmin1(quota_min,quota_min_tmp)
end do
rewind 73
end subroutine catchment_cells
