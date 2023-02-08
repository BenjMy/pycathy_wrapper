!-----------------------------------------------------------------------------
!
! CAT_DEL (Catchment Delineation)
!
! The CAT_DEL ensure the delineation of a upsole area defined the i,j of the 
! outlet cell
!
!-----------------------------------------------------------------------------

program cat_del

use mpar
use mbbio

implicit none

integer(kind=ISP) i,j,i_basin,n_rec,nodata,n_o_quota
integer(kind=ISP) i_basin_chiusura,ht,jr
integer(kind=ISP) i_out,j_out,ii,jj,ii_basin,cell_out,cells_out,outlet_cell
integer(kind=ISP) p_inflow,p_outflow_1_ciijj,p_outflow_2_ciijj
integer(kind=ISP),allocatable :: idtm(:,:)

integer(kind=IEP) :: idtm_max,idtm_min
integer(kind=ISP) :: ic
character(len=32) :: dtmfile
character(len=24) :: fstr
character(len=24) :: nstr

call rparfile('hap.in')

call load_dtm("basin_b","basin_i")

allocate(idtm(N,M)) 

open(16,file='qoi',access='direct',status='unknown',form='unformatted',recl=4)

read(16,rec=N_celle+1) i_basin_chiusura

cells_out=-9999
outlet_cell=1

! select outlet cell coordinates

write(6,*)
write(6,'(1x,a)') 'Select outlet cell coordinates (i,j): '
write(6,*)
write(6,'(1x,a3,$)') '-> '
read(5,*) i_out,j_out

! Select new value of the cells out of the catchment

write(6,*)
write(6,'(1x,a)') 'Select the new value of the cells out of the catchment'
write(6,'(1x,a)') '(-9999 is the default)'
write(6,*)
write(6,'(1x,a3,$)') '-> '
read(5,*) cells_out

! Select the value of the outlet cell of the catchment

write(6,*)
write(6,'(1x,a,$)') 'Select the value of the outlet cell of the catchment'
write(6,'(1x,a,$)') '(1 is the default)'
write(6,*)
write(6,'(1x,a3,$)') '-> '
read(5,*) outlet_cell

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

i_basin=(i_out-1)*M+j_out

if (dtm_index_pr(i_basin) == 0) then
   write(6,*) 
   write(6,'(1x,a)') 'cell out of catchment!'
   write(6,*) 
   stop
end if

! Catchment delineation

cell_out=0

do n_o_quota=N_celle,1,-1
   read(16,rec=n_o_quota+1) i_basin   
   jr=mod(i_basin,M)
   if (jr.ne.0) then
      j=jr
      i=(i_basin-j)/M+1
   else
      j=M
      i=i_basin/M
   end if
   if (dtm_index_pr(i_basin) == 0) then      	  
	  idtm(i,j)=nodata
   else
  	  idtm(i,j)=cells_out
      if (i == i_out.and.j == j_out) then
	     cell_out=n_o_quota  
		 idtm(i,j)=outlet_cell
      end if
   end if   
end do

if (cell_out /= 0) then  
  do n_o_quota=cell_out+1,1,-1   
     read(16,rec=n_o_quota) i_basin
	 jr=mod(i_basin,M)
     if (jr.ne.0) then
        j=jr
        i=(i_basin-j)/M+1
     else
        j=M
        i=i_basin/M
     end if
	 if (idtm(i,j) == cells_out) cycle  ! vuol dire che la cella corrente non appartiene al bacino	che si vuole delineare
     do ii=i-1,i+1
        do jj=j-1,j+1
          p_outflow_1_ciijj=0
		  p_outflow_2_ciijj=0
		  if (ii.eq.0.or.ii.eq.N+1) cycle
          if (jj.eq.0.or.jj.eq.M+1) cycle
          ii_basin=(ii-1)*M+jj
          if (dtm_index_pr(ii_basin) == 0) cycle
          p_outflow_1_ciijj=dtm_p_outflow_1(ii_basin)*dtm_w_1(ii_basin)
          p_outflow_2_ciijj=dtm_p_outflow_2(ii_basin)*dtm_w_2(ii_basin)
          p_inflow=3*(ii-i)+(jj-j)+5
          if (p_inflow+p_outflow_1_ciijj.eq.10) then
             idtm(ii,jj)=1
          end if
		  if (p_inflow+p_outflow_2_ciijj.eq.10) then
			 idtm(ii,jj)=1
          end if
        end do
     end do	 
  end do  !do n_o_quota=cell_out+1,1,-1 
end if ! (cell_out /= 0) then

do i=1,N
   do j=1,M
      i_basin=(i-1)*M+j
      if (dtm_index_pr(i_basin) == 0) idtm(i,j)=int(nodata)
   end do
end do   

! write an integer output file

   write(6,*)
   write(6,'(1x,a)') 'writing the output file...'
   write(6,*)

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
      write(80,'(a,1x,i5)') 'NODATA_value',nodata
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

close(16)

end program cat_del