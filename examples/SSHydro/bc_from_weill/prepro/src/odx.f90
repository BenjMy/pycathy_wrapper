!------------------------------------------------------------------------------
!
! ODX
!
! The program ODX calculates the optimal DELTA_X for the simulation according 
! to the physiographic characteristics obtained from the application of the 
! hydraulic geometry theory.
!
!------------------------------------------------------------------------------

program odx

use mpar
use mbbio

implicit none

integer(kind=ISP) :: i,j,i_basin,nsc,ncc

real(kind=REP) :: A_cell,A_1,A_2,Q_1,Q_2
real(kind=REP) :: G,beta_1,beta_2,u,w,alpha
real(kind=REP) :: ck_1,ck_2,dh_1,dh_2
real(kind=REP) :: wsum_Q,wsum_ck,wsum_dh
real(kind=REP) :: sumlog_ck,sumlog_dh
real(kind=REP) :: avg_ck,avg_dh,optimal_dx
real(kind=REP) :: avglog_ck,avglog_dh,optimallog_dx
real(kind=REP) :: adopted_dx,r1,r2
real(kind=REP) :: min_A,min_kS,min_Ws,min_S0,min_Q,min_ck,min_Dh 
real(kind=REP) :: max_A,max_kS,max_Ws,max_S0,max_Q,max_ck,max_Dh
real(kind=REP) :: sum_A,sum_kS,sum_Ws,sum_S0,sum_Q,sum_ck,sum_Dh

real(kind=RSP) :: sumnrc

call rparfile('hap.in')
call load_dtm("basin_b","basin_i")

A_cell=delta_x*delta_y

write(6,*)
write(6,'(1x,a)') 'Select the Qf = u A**w relationship parameters: '
write(6,*)
write(6,'(1x,a,$)') 'coefficient u -> '
read(5,*) u
write(6,'(1x,a,$)') 'exponent w -> '
read(5,*) w

write(6,*)
write(6,'(1x,a)') 'Select the exponent alpha of the Q**alpha weights: '
write(6,*)
write(6,'(1x,a,$)') 'alpha (0 to weight slope flows, 1 to weight channel flows) -> '
read(5,*) alpha

nsc=0
ncc=0
sumnrc=0.0E+00
min_A=1.0D+10
min_kS=1.0D+10
min_Ws=1.0D+10
min_S0=1.0D+10
min_Q=1.0D+10
min_ck=1.0D+10
min_Dh=1.0D+10
max_A=-1.0D+10
max_kS=-1.0D+10
max_Ws=-1.0D+10
max_S0=-1.0D+10
max_Q=-1.0D+10
max_ck=-1.0D+10
max_Dh=-1.0D+10
sum_A=0.0D+00
sum_kS=0.0D+00
sum_Ws=0.0D+00
sum_S0=0.0D+00
sum_Q=0.0D+00
sum_ck=0.0D+00
sum_Dh=0.0D+00
wsum_Q=0.0D+00
wsum_ck=0.0D+00
wsum_dh=0.0D+00
sumlog_ck=0.0D+00
sumlog_dh=0.0D+00

do i=1,N
   do j=1,M
      i_basin=(i-1)*M+j
      if (dtm_index_pr(i_basin) == 0) cycle
      if (dtm_hcID(i_basin).eq.0) then
         nsc=nsc+1
      else if (dtm_hcID(i_basin).eq.1) then
         ncc=ncc+1
      else
         write(6,*) 'unexpected case!'
         stop
      end if
      A_1=dtm_w_1(i_basin)*(dtm_A_inflow(i_basin)+A_cell)/dtm_nrc(i_basin)
      A_2=dtm_w_2(i_basin)*(dtm_A_inflow(i_basin)+A_cell)/dtm_nrc(i_basin)
      Q_1=u*A_1**w
      Q_2=u*A_2**w
      G=1-dtm_y1_sf(i_basin)+2.0d0/3.0d0*dtm_b1_sf(i_basin)
      beta_1=atan(max(dtm_local_slope_1(i_basin),local_slope_t))
      beta_2=atan(max(dtm_local_slope_2(i_basin),local_slope_t))
      if (dtm_w_1(i_basin) > 0.0_RSP) then
         ! write(6,*) 'i_basin = ',i_basin
         ! write(6,*) 'w1 = ',dtm_w_1(i_basin)
         ! write(6,*) 'A1 = ',A_1
         ! write(6,*) 'Q1 = ',Q_1
         ! write(6,*) 'W1 = ',dtm_Ws1_sf_1(i_basin)
         ck_1=5.0d0/(3.0d0*G)*dtm_kSs1_sf_1(i_basin)**(3.0d0/5.0d0)* &
              dtm_Ws1_sf_1(i_basin)**(-2.0d0/5.0d0)* &
              sin(beta_1)**(3.0d0/10.0d0)*Q_1**(1.0d0-3.0d0/5.0d0*G)
         dh_1=Q_1**(1-dtm_b1_sf(i_basin))/ &
              (2*G*dtm_Ws1_sf_1(i_basin)*tan(beta_1))
         sumnrc=sumnrc+dtm_nrc(i_basin)
         min_A=min(min_A,A_1)
         min_kS=min(min_kS,dtm_kSs1_sf_1(i_basin))
         min_Ws=min(min_Ws,dtm_Ws1_sf_1(i_basin))
         min_S0=min(min_S0,sin(beta_1))
         min_Q=min(min_Q,Q_1)
         min_ck=min(min_ck,ck_1)
         min_Dh=min(min_Dh,dh_1)
         max_A=max(max_A,A_1)
         max_kS=max(max_kS,dtm_kSs1_sf_1(i_basin))
         max_Ws=max(max_Ws,dtm_Ws1_sf_1(i_basin))
         max_S0=max(max_S0,sin(beta_1))
         max_Q=max(max_Q,Q_1)
         max_ck=max(max_ck,ck_1)
         max_Dh=max(max_Dh,dh_1)
         sum_A=sum_A+dtm_nrc(i_basin)*A_1
         sum_kS=sum_kS+dtm_nrc(i_basin)*dtm_kSs1_sf_1(i_basin)
         sum_Ws=sum_Ws+dtm_nrc(i_basin)*dtm_Ws1_sf_1(i_basin)
         sum_S0=sum_S0+dtm_nrc(i_basin)*sin(beta_1)
         sum_Q=sum_Q+dtm_nrc(i_basin)*Q_1
         sum_ck=sum_ck+dtm_nrc(i_basin)*ck_1
         sum_Dh=sum_Dh+dtm_nrc(i_basin)*dh_1
         wsum_Q=wsum_Q+dtm_nrc(i_basin)*Q_1**alpha
         wsum_ck=wsum_ck+dtm_nrc(i_basin)*ck_1*Q_1**alpha
         wsum_dh=wsum_dh+dtm_nrc(i_basin)*dh_1*Q_1**alpha
         sumlog_ck=sumlog_ck+dtm_nrc(i_basin)*log10(ck_1)
         sumlog_dh=sumlog_dh+dtm_nrc(i_basin)*log10(dh_1)
      end if 
      if (dtm_w_2(i_basin) > 0.0_RSP) then
         ! write(6,*) 'i_basin = ',i_basin
         ! write(6,*) 'w2 = ',dtm_w_2(i_basin)
         ! write(6,*) 'A2 = ',A_2
         ! write(6,*) 'Q2 = ',Q_2
         ! write(6,*) 'W2 = ',dtm_Ws1_sf_2(i_basin)
         ck_2=5.0d0/(3.0d0*G)*dtm_kSs1_sf_2(i_basin)**(3.0d0/5.0d0)* &
              dtm_Ws1_sf_2(i_basin)**(-2.0d0/5.0d0)* &
              sin(beta_2)**(3.0d0/10.0d0)*Q_2**(1.0d0-3.0d0/5.0d0*G)
         dh_2=Q_2**(1-dtm_b1_sf(i_basin))/ &
              (2*G*dtm_Ws1_sf_2(i_basin)*tan(beta_2))
         sumnrc=sumnrc+dtm_nrc(i_basin)
         min_A=min(min_A,A_2)
         min_kS=min(min_kS,dtm_kSs1_sf_2(i_basin))
         min_Ws=min(min_Ws,dtm_Ws1_sf_2(i_basin))
         min_S0=min(min_S0,sin(beta_2))
         min_Q=min(min_Q,Q_2)
         min_ck=min(min_ck,ck_2)
         min_Dh=min(min_Dh,dh_2)
         max_A=max(max_A,A_2)
         max_kS=max(max_kS,dtm_kSs1_sf_2(i_basin))
         max_Ws=max(max_Ws,dtm_Ws1_sf_2(i_basin))
         max_S0=max(max_S0,sin(beta_2))
         max_Q=max(max_Q,Q_2)
         max_ck=max(max_ck,ck_2)
         max_Dh=max(max_Dh,dh_2)
         sum_A=sum_A+dtm_nrc(i_basin)*A_2
         sum_kS=sum_kS+dtm_nrc(i_basin)*dtm_kSs1_sf_2(i_basin)
         sum_Ws=sum_Ws+dtm_nrc(i_basin)*dtm_Ws1_sf_2(i_basin)
         sum_S0=sum_S0+dtm_nrc(i_basin)*sin(beta_2)
         sum_Q=sum_Q+dtm_nrc(i_basin)*Q_2
         sum_ck=sum_ck+dtm_nrc(i_basin)*ck_2
         sum_Dh=sum_Dh+dtm_nrc(i_basin)*dh_2
         wsum_Q=wsum_Q+dtm_nrc(i_basin)*Q_2**alpha
         wsum_ck=wsum_ck+dtm_nrc(i_basin)*ck_2*Q_2**alpha
         wsum_dh=wsum_dh+dtm_nrc(i_basin)*dh_2*Q_2**alpha
         sumlog_ck=sumlog_ck+dtm_nrc(i_basin)*log10(ck_2)
         sumlog_dh=sumlog_dh+dtm_nrc(i_basin)*log10(dh_2)
      end if
  end do
end do

avg_ck=wsum_ck/wsum_Q
avg_dh=wsum_dh/wsum_Q

avglog_ck=10**(sumlog_ck/sumnrc)
avglog_dh=10**(sumlog_dh/sumnrc)

optimal_dx=2*avg_dh/avg_ck
optimallog_dx=2*avglog_dh/avglog_ck

write(6,*) 
write(6,'(1x,a,i8,a,i8,a,i8,a)') &
'MEAN (MIN,MAX) OVER',nint(sumnrc),' PATHS,', &
nsc,' SLOPE CELLS,',ncc,' CHANNEL CELLS'
write(6,'(1x,a,f20.6,1x,a1,f20.6,a1,f20.6,a1)') &
'<A> =       ',sum_A/sumnrc,'(',min_A,',',max_A,')'
write(6,'(1x,a,f20.6,1x,a1,f20.6,a1,f20.6,a1)') &
'<kS(A,1)> = ',sum_kS/sumnrc,'(',min_kS,',',max_kS,')'
write(6,'(1x,a,f20.6,1x,a1,f20.6,a1,f20.6,a1)') &
'<W(A,1)> =  ',sum_Ws/sumnrc,'(',min_Ws,',',max_Ws,')'
write(6,'(1x,a,f20.6,1x,a1,f20.6,a1,f20.6,a1)') &
'<S0> =      ',sum_S0/sumnrc,'(',min_S0,',',max_S0,')'
write(6,'(1x,a,f20.6,1x,a1,f20.6,a1,f20.6,a1)') &
'<Q>  =      ',sum_Q/sumnrc,'(',min_Q,',',max_Q,')'
write(6,'(1x,a,f20.6,1x,a1,f20.6,a1,f20.6,a1)') &
'<ck> =      ',sum_ck/sumnrc,'(',min_ck,',',max_ck,')'
write(6,'(1x,a,f20.6,1x,a1,f20.6,a1,f20.6,a1)') &
'<Dh> =      ',sum_Dh/sumnrc,'(',min_Dh,',',max_Dh,')'
write(6,*) 
write(6,'(1x,a)') 'Q**ALPHA-WEIGHTED ARITMETIC AVERAGES'
write(6,'(1x,a,f13.6)') 'average kinematic celerity =    ',avg_ck
write(6,'(1x,a,f13.6)') 'average hydraulic diffusivity = ',avg_dh
write(6,'(1x,a,f13.6)') 'optimal space interval =        ',optimal_dx
write(6,*) 
! write(6,'(1x,a,f13.6)') 'AVERAGES ON LOG VALUES'
! write(6,'(1x,a,f13.6)') 'average kinematic celerity =    ',avglog_ck
! write(6,'(1x,a,f13.6)') 'average hydraulic diffusivity = ',avglog_dh
! write(6,'(1x,a,f13.6)') 'optimal space interval =        ',optimallog_dx
! write(6,*)
write(6,'(1x,a,$)') 'adopted space interval -> '
read(5,*) adopted_dx
r1=adopted_dx/avg_ck
r2=2*avg_dh/avg_ck**2
if (adopted_dx.ge.2*avg_dh/avg_ck) then
   write(6,*)
   write(6,'(1x,a,f13.6,1x,a1,f13.6,a1,f13.6,a1)') &
   'optimal time interval = ',r1,'(',r1-r2,',',r1+r2,')'
   write(6,*)
else
   write(6,*)
   write(6,'(1x,a,f13.6,1x,a1,f13.6,1x,a1,f13.6,a1)') &
   'optimal time interval = ',r2,'(',r2-r1,',',r2+r1,')'
   write(6,*)
end if

end program odx
