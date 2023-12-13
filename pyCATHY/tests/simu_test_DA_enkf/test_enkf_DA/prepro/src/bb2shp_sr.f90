!-----------------------------------------------------------------------------
!
! BB2SHP_SR
!
! The program reads the parameters contained in the binary files
! basin_b/basin_i and write a *.SHP ArcGIS file.
!
!-----------------------------------------------------------------------------

subroutine bb2shp_sr

use mpar
use mbbio
use ShapeFile
use DBF

implicit none

integer(kind=ISP) i,j,i_basin,iii,jjj,i_cv,j_cv,id
real(kind=REP) w_1,w_2

type(PolyLine)::poly_tmp, polig_tmp
type(DBFfile)::db_main

allocate(poly_tmp%Points(2))
allocate(polig_tmp%Points(4))

call rparfile('hap.in')
call load_dtm("basin_b","basin_i")

! network polyline 

write(6,*) 'writing file river_net.shp'

call IniziaScrittura('river_net',3) ! 3 Polyline

db_main%NomeFile='river_net.dbf'
db_main%uni=40
db_main%nFields=10 
db_main%scrittura=.TRUE.
allocate(db_main%F(db_main%nFields))

db_main%F(1)%FieldName='localslope'
db_main%F(1)%FieldType='F'
db_main%F(1)%FieldLength=11
db_main%F(1)%decimal=6

db_main%F(2)%FieldName='epl'
db_main%F(2)%FieldType='F'
db_main%F(2)%FieldLength=12
db_main%F(2)%decimal=6

db_main%F(3)%FieldName='hcID'
db_main%F(3)%FieldType='N'
db_main%F(3)%FieldLength=8
db_main%F(3)%decimal=0

db_main%F(4)%FieldName='Ws1_sf'
db_main%F(4)%FieldType='F'
db_main%F(4)%FieldLength=12
db_main%F(4)%decimal=6

db_main%F(5)%FieldName='b1_sf'
db_main%F(5)%FieldType='F'
db_main%F(5)%FieldLength=12
db_main%F(5)%decimal=6

db_main%F(6)%FieldName='kSs1_sf'
db_main%F(6)%FieldType='F'
db_main%F(6)%FieldLength=12
db_main%F(6)%decimal=6

db_main%F(7)%FieldName='y1_sf'
db_main%F(7)%FieldType='F'
db_main%F(7)%FieldLength=12
db_main%F(7)%decimal=6

db_main%F(8)%FieldName='p_outflow'
db_main%F(8)%FieldType='N'
db_main%F(8)%FieldLength=4
db_main%F(8)%decimal=0

db_main%F(9)%FieldName='sso'
db_main%F(9)%FieldType='N'
db_main%F(9)%FieldLength=4
db_main%F(9)%decimal=0

db_main%F(10)%FieldName='hso'
db_main%F(10)%FieldType='N'
db_main%F(10)%FieldLength=4
db_main%F(10)%decimal=0

call CreaDB(db_main%NomeFile, db_main%nFields, db_main%F%FieldName, db_main%F%FieldType, &
   & db_main%F%FieldLength, db_main%F%decimal, db_main%uni)

call InitDbf(db_main%NomeFile, db_main, db_main%uni, db_main%scrittura)

do i=1,N
   do j=1,M
      i_basin=(i-1)*M+j
      if (dtm_index_pr(i_basin) == 0) cycle
	  w_1=dtm_w_1(i_basin)
	  w_2=dtm_w_2(i_basin)
      if (w_1 /= 0.0D0) then	     
         iii=nint((dtm_p_outflow_1(i_basin)-5)/3.)
         jjj=dtm_p_outflow_1(i_basin)-5-3*iii
         i_cv=i+iii 
         j_cv=j+jjj
		 poly_tmp%Points(1)%x=xllcorner+delta_x/2+(i-1)*delta_x
         poly_tmp%Points(1)%y=yllcorner+delta_y/2+(j-1)*delta_y		 
		 poly_tmp%Points(2)%x=xllcorner+delta_x/2+(i_cv-1)*delta_x
         poly_tmp%Points(2)%y=yllcorner+delta_y/2+(j_cv-1)*delta_y
		 call ScriviPolyline(poly_tmp)
	     call Add_Record(db_main)
	     call Set_Field_real(db_main, db_main%F(1)%FieldName, dtm_local_slope_1(i_basin))
	     call Set_Field_real(db_main, db_main%F(2)%FieldName, dtm_epl_1(i_basin))
	     call Set_Field_integer(db_main, db_main%F(3)%FieldName, dtm_hcID(i_basin))
	     call Set_Field_real(db_main, db_main%F(4)%FieldName, dtm_Ws1_sf_1(i_basin))
	     call Set_Field_real(db_main, db_main%F(5)%FieldName, dtm_b1_sf(i_basin))
         call Set_Field_real(db_main, db_main%F(6)%FieldName, dtm_KSs1_sf_1(i_basin))
	     call Set_Field_real(db_main, db_main%F(7)%FieldName, dtm_y1_sf(i_basin))
 	     call Set_Field_integer(db_main, db_main%F(8)%FieldName, dtm_p_outflow_1(i_basin))
		 call Set_Field_integer(db_main, db_main%F(9)%FieldName, dtm_sso(i_basin))
		 call Set_Field_integer(db_main, db_main%F(10)%FieldName, dtm_hso(i_basin))
	     call Scrivi_Record(db_main)
	  end if
	  if (w_2 /= 0.0D0) then
		 iii=nint((dtm_p_outflow_2(i_basin)-5)/3.)
         jjj=dtm_p_outflow_2(i_basin)-5-3*iii
         i_cv=i+iii 
         j_cv=j+jjj
		 poly_tmp%Points(1)%x=xllcorner+delta_x/2+(i-1)*delta_x
         poly_tmp%Points(1)%y=yllcorner+delta_y/2+(j-1)*delta_y		 
		 poly_tmp%Points(2)%x=xllcorner+delta_x/2+(i_cv-1)*delta_x
         poly_tmp%Points(2)%y=yllcorner+delta_y/2+(j_cv-1)*delta_y
		 call ScriviPolyline(poly_tmp)
	     call Add_Record(db_main)
	     call Set_Field_real(db_main, db_main%F(1)%FieldName, dtm_local_slope_2(i_basin))
	     call Set_Field_real(db_main, db_main%F(2)%FieldName, dtm_epl_2(i_basin))
	     call Set_Field_integer(db_main, db_main%F(3)%FieldName, dtm_hcID(i_basin))
	     call Set_Field_real(db_main, db_main%F(4)%FieldName, dtm_Ws1_sf_2(i_basin))
	     call Set_Field_real(db_main, db_main%F(5)%FieldName, dtm_b1_sf(i_basin))
         call Set_Field_real(db_main, db_main%F(6)%FieldName, dtm_kSs1_sf_2(i_basin))
	     call Set_Field_real(db_main, db_main%F(7)%FieldName, dtm_y1_sf(i_basin))
 	     call Set_Field_integer(db_main, db_main%F(8)%FieldName, dtm_p_outflow_2(i_basin))
		 call Set_Field_integer(db_main, db_main%F(9)%FieldName, dtm_sso(i_basin))
		 call Set_Field_integer(db_main, db_main%F(10)%FieldName, dtm_hso(i_basin))		 
	     call Scrivi_Record(db_main)
	  end if   
   end do
end do

call ChiudiScrittura()
call CloseDbf(db_main) 

! poligon cells

write(6,*) 'writing file cells.shp'

call IniziaScrittura('cells',5) ! 5 Poligon

db_main%NomeFile='cells.dbf'
db_main%uni=41
db_main%nFields=16 
db_main%scrittura=.TRUE.
allocate(db_main%F(db_main%nFields))

id=1
db_main%F(id)%FieldName='i'
db_main%F(id)%FieldType='N'
db_main%F(id)%FieldLength=10
db_main%F(id)%decimal=0
id=id+1
db_main%F(id)%FieldName='j'
db_main%F(id)%FieldType='N'
db_main%F(id)%FieldLength=10
db_main%F(id)%decimal=0
id=id+1
db_main%F(id)%FieldName='i_basin'
db_main%F(id)%FieldType='N'
db_main%F(id)%FieldLength=10
db_main%F(id)%decimal=0
id=id+1
db_main%F(id)%FieldName='quota'
db_main%F(id)%FieldType='F'
db_main%F(id)%FieldLength=12
db_main%F(id)%decimal=2
id=id+1
db_main%F(id)%FieldName='p_outflow1'
db_main%F(id)%FieldType='N'
db_main%F(id)%FieldLength=10
db_main%F(id)%decimal=0
id=id+1
db_main%F(id)%FieldName='p_outflow2'
db_main%F(id)%FieldType='N'
db_main%F(id)%FieldLength=10
db_main%F(id)%decimal=0
id=id+1
db_main%F(id)%FieldName='dmID'
db_main%F(id)%FieldType='N'
db_main%F(id)%FieldLength=6
db_main%F(id)%decimal=0
id=id+1
db_main%F(id)%FieldName='A_inflow'
db_main%F(id)%FieldType='F'
db_main%F(id)%FieldLength=20
db_main%F(id)%decimal=2
id=id+1
db_main%F(id)%FieldName='AS^k'
db_main%F(id)%FieldType='F'
db_main%F(id)%FieldLength=20
db_main%F(id)%decimal=4
id=id+1
db_main%F(id)%FieldName='DN'
db_main%F(id)%FieldType='F'
db_main%F(id)%FieldLength=20
db_main%F(id)%decimal=10
id=id+1
db_main%F(id)%FieldName='w_1'
db_main%F(id)%FieldType='F'
db_main%F(id)%FieldLength=12
db_main%F(id)%decimal=6
id=id+1
db_main%F(id)%FieldName='w_2'
db_main%F(id)%FieldType='F'
db_main%F(id)%FieldLength=12
db_main%F(id)%decimal=6
id=id+1
db_main%F(id)%FieldName='sumdev_num'
db_main%F(id)%FieldType='F'
db_main%F(id)%FieldLength=12
db_main%F(id)%decimal=6
id=id+1
db_main%F(id)%FieldName='lakes_map'
db_main%F(id)%FieldType='N'
db_main%F(id)%FieldLength=6
db_main%F(id)%decimal=0
id=id+1
db_main%F(id)%FieldName='zone'
db_main%F(id)%FieldType='N'
db_main%F(id)%FieldLength=6
db_main%F(id)%decimal=0
id=id+1
db_main%F(id)%FieldName='q_output'
db_main%F(id)%FieldType='N'
db_main%F(id)%FieldLength=6
db_main%F(id)%decimal=0
id=id+1
call CreaDB(db_main%NomeFile, db_main%nFields, db_main%F%FieldName, db_main%F%FieldType, &
   & db_main%F%FieldLength, db_main%F%decimal, db_main%uni)

call InitDbf(db_main%NomeFile, db_main, db_main%uni, db_main%scrittura)

do i=1,N
   do j=1,M
      i_basin=(i-1)*M+j
	  polig_tmp%Points(1)%x=xllcorner+(i-1)*delta_x
	  polig_tmp%Points(1)%y=yllcorner+(j-1)*delta_y
	  polig_tmp%Points(2)%x=xllcorner+(i)*delta_x
      polig_tmp%Points(2)%y=yllcorner+(j-1)*delta_y
	  polig_tmp%Points(3)%x=xllcorner+(i)*delta_x
	  polig_tmp%Points(3)%y=yllcorner+(j)*delta_y
	  polig_tmp%Points(4)%x=xllcorner+(i-1)*delta_x
      polig_tmp%Points(4)%y=yllcorner+(j)*delta_y
      if (dtm_index_pr(i_basin) == 0) then
		 call ScriviPolyline(polig_tmp)
	     call Add_Record(db_main)
	     call Set_Field_integer(db_main, db_main%F(1)%FieldName, -9999)
	     call Set_Field_integer(db_main, db_main%F(2)%FieldName, -9999)
	     call Set_Field_integer(db_main, db_main%F(3)%FieldName, -9999)
	     call Set_Field_reald(db_main, db_main%F(4)%FieldName, -9999.0_REP)
	     call Set_Field_integer(db_main, db_main%F(5)%FieldName, -9999)
         call Set_Field_integer(db_main, db_main%F(6)%FieldName, -9999)
	     call Set_Field_integer(db_main, db_main%F(7)%FieldName, -9999)
 	     call Set_Field_reald(db_main, db_main%F(8)%FieldName, -9999.0_REP)
 	     call Set_Field_real(db_main, db_main%F(9)%FieldName, -9999.0_RSP)
 	     call Set_Field_real(db_main, db_main%F(10)%FieldName, -9999.0_RSP)
 	     call Set_Field_real(db_main, db_main%F(11)%FieldName, -9999.0_RSP)
 	     call Set_Field_real(db_main, db_main%F(12)%FieldName, -9999.0_RSP)
 	     call Set_Field_reald(db_main, db_main%F(13)%FieldName, -9999.0_REP)
	     call Set_Field_integer(db_main, db_main%F(14)%FieldName, -9999)
	     call Set_Field_integer(db_main, db_main%F(15)%FieldName, -9999)
	     call Set_Field_integer(db_main, db_main%F(16)%FieldName, -9999)
	     call Scrivi_Record(db_main)
	  else
		 call ScriviPolyline(polig_tmp)
	     call Add_Record(db_main)
	     call Set_Field_integer(db_main, db_main%F(1)%FieldName, i)
	     call Set_Field_integer(db_main, db_main%F(2)%FieldName, j)
	     call Set_Field_integer(db_main, db_main%F(3)%FieldName, i_basin)
	     call Set_Field_reald(db_main, db_main%F(4)%FieldName, dtm_quota(i_basin))
	     call Set_Field_integer(db_main, db_main%F(5)%FieldName, dtm_p_outflow_1(i_basin))
         call Set_Field_integer(db_main, db_main%F(6)%FieldName, dtm_p_outflow_2(i_basin))
	     call Set_Field_integer(db_main, db_main%F(7)%FieldName, dtm_dmID(i_basin))
 	     call Set_Field_reald(db_main, db_main%F(8)%FieldName, dtm_A_inflow(i_basin))
  	     call Set_Field_real(db_main, db_main%F(9)%FieldName, dtm_ASk(i_basin))
  	     call Set_Field_real(db_main, db_main%F(10)%FieldName, dtm_DN(i_basin))
 	     call Set_Field_real(db_main, db_main%F(11)%FieldName, dtm_w_1(i_basin))
 	     call Set_Field_real(db_main, db_main%F(12)%FieldName, dtm_w_2(i_basin))
 	     call Set_Field_reald(db_main, db_main%F(13)%FieldName, dtm_sumdev_num(i_basin))
	     call Set_Field_integer(db_main, db_main%F(14)%FieldName, dtm_lakes_map(i_basin))
	     call Set_Field_integer(db_main, db_main%F(15)%FieldName, dtm_zone(i_basin))
	     call Set_Field_integer(db_main, db_main%F(16)%FieldName, dtm_q_output(i_basin))
	     call Scrivi_Record(db_main)
	  end if   
   end do
end do

call ChiudiScrittura()
call CloseDbf(db_main) 


end subroutine bb2shp_sr
