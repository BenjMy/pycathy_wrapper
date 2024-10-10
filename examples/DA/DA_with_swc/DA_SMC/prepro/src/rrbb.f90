!-----------------------------------------------------------------------------
!
! RRBB
!
! The program RRBB reads the record of parameters contained in the binary 
! files basin_b/basin_i for a selected cell (i,j).
!
!-----------------------------------------------------------------------------

program rrbb

use mpar
use mbbio

implicit none

integer(kind=ISP) i,j,i_basin

call rparfile('hap.in')

call load_dtm("basin_b","basin_i")

write(6,*)
write(6,'(1x,a,$)') 'Select cell coordinates (i,j): '
read(5,*) i,j

i_basin=(i-1)*M+j

if (dtm_index_pr(i_basin) == 0) then
   write(6,*) 
   write(6,'(1x,a)') 'cell out of catchment!'
   write(6,*) 
   stop
end if

write(6,*) 
write(6,'(1x,a70)') '----------------------------------------------------------------------'
write(6,'(1x,a6,i4,a1,i4,a1,35x,a11,i7,a1)') 'Cell (',i,',',j,')','(i_basin = ',i_basin,')'
write(6,'(1x,a70)') '----------------------------------------------------------------------'

write(6,10) ' 1) vh = ',dtm_vh(i_basin),'21) local_slope_1 = ',dtm_local_slope_1(i_basin)
10 format(1x,a9,e13.6,15x,a20,e13.6)

write(6,11) ' 2) rl_min = ',dtm_rl_min(i_basin),'22) local_slope_2 = ',dtm_local_slope_2(i_basin)
11 format(1x,a13,e13.6,11x,a19,e13.6)

write(6,12) ' 3) LAI = ',dtm_LAI(i_basin),'23) Ws1_sf_1 = ',dtm_Ws1_sf_1(i_basin)
12 format(1x,a10,e13.6,14x,a15,e13.6)

write(6,13) ' 4) vsc = ',dtm_vsc(i_basin),'24) Ws1_sf_2 = ',dtm_Ws1_sf_2(i_basin)
13 format(1x,a10,e13.6,14x,a15,e13.6)

write(6,14) ' 5) theta_s = ',dtm_theta_s(i_basin),'25) b1_sf = ',dtm_b1_sf(i_basin)
14 format(1x,a14,e13.6,10x,a12,e13.6)

write(6,15) ' 6) eta = ',dtm_eta(i_basin),'26) kSs1_sf_1 = ',dtm_kSs1_sf_1(i_basin)
15 format(1x,a10,e13.6,14x,a16,e13.6)

write(6,16) ' 7) psi_s = ',dtm_psi_s(i_basin),'27) kSs1_sf_2 = ',dtm_kSs1_sf_2(i_basin)
16 format(1x,a12,e13.6,12x,a16,e13.6)

write(6,17) ' 8) K_s = ',dtm_K_s(i_basin),'28) y1_sf = ',dtm_y1_sf(i_basin)
17 format(1x,a10,e13.6,14x,a12,e13.6)

write(6,18) ' 9) Z_low = ',dtm_Z_low(i_basin),'29) ASk = ',dtm_ASk(i_basin)
18 format(1x,a12,e13.6,12x,a10,e13.6)

write(6,19) '10) Z_up = ',dtm_Z_up(i_basin),'30) DN = ',dtm_DN(i_basin)
19 format(1x,a11,e13.6,13x,a9,e13.6)

write(6,20) '11) theta_r = ',dtm_theta_r(i_basin),'31) hcID = ',dtm_hcID(i_basin)
20 format(1x,a14,e13.6,10x,a11,i2)

write(6,21) '12) rl_max = ',dtm_rl_max(i_basin),'32) epl_1 = ',dtm_epl_1(i_basin)
21 format(1x,a13,e13.6,11x,a12,e13.6)

write(6,22) '13) quota = ',dtm_quota(i_basin),'33) epl_2 = ',dtm_epl_2(i_basin)
22 format(1x,a12,f10.3,15x,a12,e13.6)

write(6,23) '14) p_outflow_1 = ',dtm_p_outflow_1(i_basin),'34) lakes_map = ',dtm_lakes_map(i_basin)
23 format(1x,a18,i2,17x,a,i2)

write(6,24) '15) p_outflow_2 = ',dtm_p_outflow_2(i_basin),'35) zone = ',dtm_zone(i_basin)
24 format(1x,a18,i2,17x,a11,i2)

write(6,25) '16) dmID = ',dtm_dmID(i_basin),'36) q_output = ',dtm_q_output(i_basin)
25 format(1x,a11,i2,24x,a15,i2)

write(6,26) '17) A_inflow = ',dtm_A_inflow(i_basin),'37) nrc = ',dtm_nrc(i_basin)
26 format(1x,a15,e13.6,9x,a10,e13.6)

write(6,27) '18) w_1 = ',dtm_w_1(i_basin)
27 format(1x,a10,f9.6)

write(6,28) '19) w_2 = ',dtm_w_2(i_basin)
28 format(1x,a10,f9.6)

write(6,29) '20) sumdev_num = ',dtm_sumdev_num(i_basin)
29 format(1x,a17,e13.6)

write(6,'(1x,a70)') '----------------------------------------------------------------------'
write(6,*) 


end program rrbb
