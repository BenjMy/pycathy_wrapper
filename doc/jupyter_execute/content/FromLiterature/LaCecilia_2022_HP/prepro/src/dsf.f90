!-----------------------------------------------------------------------------
!
! DSF (Drainage System Features)
!
! The subroutine DSF determines two drainage directions for each DEM cell.
! Depending on parameter dmID, one or two drainage directions are actually
! employed. For each DEM cell, the following physiographic features are 
! assigned: (1) weights associated to drainage directions: w_1 and w_2;
! (2) flows directions: p_outflow_1 and p_outflow_2; (3) elemental lengths
! of channels: epl_1 and epl_2 (it is assumend here that delta_x=delta_y); 
! (4) local slopes (positive downward) of channel elements: local_slope_1 
! and local_slope_2; and (5) the funcion ASk=AS**k.
! For at most the two neighbors spilling cells are also updated: (1) the
! upstream catchment area A_inflow and (2) the upstream error made between 
! selected and theoretical drainage directions sumdev_num.
!
!-----------------------------------------------------------------------------

subroutine dsf()

use mpar
use mbbio

implicit none

integer(kind=ISP) i,j,i_basin,jr,l
integer(kind=ISP) ii,jj,ii_basin,ii_basin_max
integer(kind=ISP) i_basin_chiusura
integer(kind=ISP) n_o_quota
real(kind=REP) quota_ciijj,quota_max
real(kind=REP) e(9)
real(kind=REP) e0,e1,e2,e1_fmax,e2_fmax,emin,s_max,r_max,s_max_facet,r
real(kind=REP) sigma
real(kind=REP) a1,a2
integer(kind=ISP) iii_1,jjj_1,iii_2,jjj_2
integer(kind=ISP) i_basin_cv_1,i_basin_cv_2
real(kind=REP) dev_1,dev_2,sumdev,sumdev_1,sumdev_2
real(kind=RSP) w_1,w_2
real(kind=REP) sumdev_num,sumdev_num_cv_1,sumdev_num_cv_2 
integer(kind=ISP) p_outflow_1,p_outflow_2,p_outflow_outlet,dmID
real(kind=REP) A_inflow,A_cell,A_outflow
real(kind=REP) A_inflow_cv_1,A_inflow_cv_2,A_outflow_n
real(kind=RSP) local_slope_1,local_slope_2,epl_1,epl_2
real(kind=REP) delta_xy
real(kind=REP) pi,rad2
real(kind=RSP) ASk,Kp,DN
real(kind=REP) A_outflow_ciijj,A_inflow_max
real(kind=RSP) local_slope_outlet
integer(kind=ISP) p_inflow,p_outflow_1_ciijj,p_outflow_2_ciijj,p_outflow_L
integer(kind=ISP) ivo,jvo,ivo_basin,nvo
integer(kind=ISP) hcID
integer, external :: ChannelInitiation

integer(kind=ISP) sso,sso_cv,sso_jun_cv

pi=4.0_REP*datan(1.0_REP)
rad2=dsqrt(2.0_REP)
delta_xy=dsqrt(2.0_REP)*delta_x ! assuming delta_x=delta_y
A_cell=delta_x*delta_x

! open the binary file qoi

open(16,file='qoi',access='direct',status='unknown',form='unformatted',recl=4)

! identify the outlet cell (the lowest catchment cell)

read(16,rec=N_celle+1) i_basin_chiusura

! WRITE(6,*) 'i_basin_chiusura = ',i_basin_chiusura

do n_o_quota=1,N_celle

   read(16,rec=n_o_quota+1) i_basin
   ! drainage directions previously set up are not altered
   if (dtm_p_outflow_1(i_basin).ne.0.or.dtm_p_outflow_2(i_basin).ne.0) cycle
   jr=mod(i_basin,M)
   if (jr.ne.0) then
      j=jr
      i=(i_basin-j)/M+1
   else
      j=M
      i=i_basin/M
   end if
   if (i_basin.eq.i_basin_chiusura) cycle ! features are assigned apart

   l=0
   do ii=i-1,i+1
      do jj=j-1,j+1
         l=l+1
         e(l)=0.0_REP
         if (ii.eq.0.or.ii.eq.N+1) cycle
         if (jj.eq.0.or.jj.eq.M+1) cycle
         ii_basin=(ii-1)*M+jj
         quota_ciijj=dtm_quota(ii_basin)
         if (quota_ciijj < 0.0_REP) cycle
         e(l)=quota_ciijj
      end do
   end do

   e0=e(5)
   s_max=0.0_REP
!
! triangle 021
!
   e1=e(2)
   e2=e(1)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) then
         e1_fmax=e1
         e2_fmax=e2
         r_max=r
         s_max=s_max_facet
         p_outflow_1=2
         p_outflow_2=1
         sigma=+1.0_REP
      end if
   end if
!
! triangle 023
!
   e1=e(2)
   e2=e(3)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) then
         e1_fmax=e1
         e2_fmax=e2
         r_max=r
         s_max=s_max_facet
         p_outflow_1=2
         p_outflow_2=3
         sigma=-1.0_REP
      end if
   end if
!
! triangle 063
!
   e1=e(6)
   e2=e(3)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) then
         e1_fmax=e1
         e2_fmax=e2
         r_max=r
         s_max=s_max_facet
         p_outflow_1=6
         p_outflow_2=3
         sigma=+1.0_REP
      end if
   end if
!
! triangle 069
!
   e1=e(6)
   e2=e(9)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) then
         e1_fmax=e1
         e2_fmax=e2
         r_max=r
         s_max=s_max_facet
         p_outflow_1=6
         p_outflow_2=9
         sigma=-1.0_REP
      end if
   end if
!
! triangle 089
!
   e1=e(8)
   e2=e(9)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) then
         e1_fmax=e1
         e2_fmax=e2
         r_max=r
         s_max=s_max_facet
         p_outflow_1=8
         p_outflow_2=9
         sigma=+1.0_REP
      end if
   end if
!
! triangle 087
!
   e1=e(8)
   e2=e(7)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) then
         e1_fmax=e1
         e2_fmax=e2
         r_max=r
         s_max=s_max_facet
         p_outflow_1=8
         p_outflow_2=7
         sigma=-1.0_REP
      end if
   end if
!
! triangle 047
!
   e1=e(4)
   e2=e(7)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) then
         e1_fmax=e1
         e2_fmax=e2
         r_max=r
         s_max=s_max_facet
         p_outflow_1=4
         p_outflow_2=7
         sigma=+1.0_REP
      end if
   end if
!
! triangle 041
!
   e1=e(4)
   e2=e(1)
   if (e1*e2.ne.0.0_REP) then
      call facet(e0,e1,e2,r,s_max_facet)
      if (s_max_facet.gt.s_max) then
         e1_fmax=e1
         e2_fmax=e2
         r_max=r
         s_max=s_max_facet
         p_outflow_1=4
         p_outflow_2=1
         sigma=-1.0_REP
      end if
   end if

   if (s_max.gt.0.0_REP) then

      ! calculate deviations dev_1 and dev_2

      if (imethod.eq.1) then ! D8-LAD method
         dev_1=r_max ! delta_x=delta_y
         dev_2=pi/4-r_max ! delta_x=delta_y
      else if (imethod.eq.2) then ! D8-LTD method
         dev_1=delta_x*dsin(r_max) ! delta_x=delta_y
         dev_2=delta_x*rad2*dsin(pi/4.0_REP-r_max) ! delta_x=delta_y
      else
         write(6,*) 'unespected imethod!'
         stop
      end if
      if (sigma.eq.1.0_REP) then
         dev_2=-dev_2
      else ! if (sigma.eq.-1) then
         dev_1=-dev_1
      end if

      ! load information calculated in previous cycles

      A_inflow=dtm_A_inflow(i_basin)
      A_outflow=A_inflow+A_cell
      sumdev_num=dtm_sumdev_num(i_basin)
      if (A_inflow.eq.0.0_REP) then
         sumdev=0.0_REP   
      else
         sumdev=sumdev_num/A_inflow
      end if 
      if (dabs(dev_1).le.epsilon(dev_1).or.dabs(dev_2).le.epsilon(dev_2)) then
         sumdev=0.0_REP
      end if
      dmID=dtm_dmID(i_basin)
      hcID=dtm_hcID(i_basin)
      Kp=dtm_DN(i_basin)
       
      ! calculate local phisiographic features of the current cell

      epl_1=delta_x
      epl_2=delta_xy
      local_slope_1=(e0-e1_fmax)/epl_1
      local_slope_2=(e0-e2_fmax)/epl_2
      ASk=A_outflow*s_max**kas 
      sumdev_1=lambda*sumdev+dev_1
      sumdev_2=lambda*sumdev+dev_2
      
      DN=Kp/(-mean_s_max)
      if (hcID.eq.0) hcID=ChannelInitiation(A_outflow,ASk,DN)
      if (hcID.eq.1.and.ndcf.eq.1) dmID=2

      a1=dabs(sumdev_1)
      a2=dabs(sumdev_2)

      select case (dmID)

      case (1) ! dispersive method (D_inf)
	      
         if (dabs(dev_1).le.epsilon(dev_1)) then
            w_1=1.0_RSP
            w_2=0.0_RSP
         else if (dabs(dev_2).le.epsilon(dev_2)) then
            w_1=0.0_RSP
            w_2=1.0_RSP
         else 
            w_1=a2/(a1+a2)
            w_2=a1/(a1+a2)
            if (w_1.lt.1.0E-06) then
               write(6,*) &
               'weight w_1 less than 1.0E-06: set w_1 = 0 and w_2 = 1'
               w_1=0.0_RSP
               w_2=1.0_RSP
            end if
            if (w_2.lt.1.0E-06) then
               write(6,*) &
               'weight w_2 less than 1.0E-06: set w_1 = 1 and w_2 = 0'
               w_1=1.0_RSP
               w_2=0.0_RSP
            end if
         end if

      case (2)  ! nondispersive method (D8-LTD)
    
		 if (dabs(a1-a2)/delta_x.lt.10D-14.and.(e0-e1_fmax).gt.0.0_REP) then            
            w_1=1.0_RSP
            w_2=0.0_RSP
            epl_2=0.0
            local_slope_2=0.0
		 else
		    if (a1.lt.a2.and.(e0-e1_fmax).gt.0.0_REP) then 
               w_1=1.0_RSP
               w_2=0.0_RSP
               epl_2=0.0
               local_slope_2=0.0
            else if (a1.gt.a2.or.(e0-e2_fmax).gt.0.0_REP) then 
               w_1=0.0_RSP
               w_2=1.0_RSP
               epl_1=0.0
               local_slope_1=0.0
            else ! non serve
               write(6,*) 's_max < 0, unexpected case!'
               stop
            end if		 
		 end if

      case default

         stop "unexpected case!" 

      end select

      ! calculate features of the two following cells

      ! the function p_outflow is the coordinate z of the plane ax+by+c=z,
      ! with a=3, b=1 and c=5. the function (i_basin_cv-i_basin) is the
      ! coordinate z of the plane ax+by+c=z, with a=M, b=1 e c=0.

      iii_1=nint((p_outflow_1-5)/3.)
      jjj_1=p_outflow_1-5-3*iii_1
      i_basin_cv_1=i_basin+(M*iii_1+jjj_1)
      A_inflow_cv_1=dtm_A_inflow(i_basin_cv_1)
      sumdev_num_cv_1=dtm_sumdev_num(i_basin_cv_1)	  
      iii_2=nint((p_outflow_2-5)/3.)
      jjj_2=p_outflow_2-5-3*iii_2
      i_basin_cv_2=i_basin+(M*iii_2+jjj_2)
      A_inflow_cv_2=dtm_A_inflow(i_basin_cv_2)
      sumdev_num_cv_2=dtm_sumdev_num(i_basin_cv_2)
      A_inflow_cv_1=A_inflow_cv_1+(A_outflow*w_1)
      A_inflow_cv_2=A_inflow_cv_2+(A_outflow*w_2)      
      sumdev_num_cv_1=sumdev_num_cv_1+(A_outflow*w_1*sumdev_1)
      sumdev_num_cv_2=sumdev_num_cv_2+(A_outflow*w_2*sumdev_2)
      call set_w_1(i_basin,w_1)
      call set_w_2(i_basin,w_2)
      call set_p_outflow_1(i_basin,p_outflow_1)
      call set_p_outflow_2(i_basin,p_outflow_2)
      call set_dmID(i_basin,dmID)
      call set_epl_1(i_basin,epl_1)
      call set_epl_2(i_basin,epl_2)
      call set_local_slope_1(i_basin,local_slope_1)
      call set_local_slope_2(i_basin,local_slope_2)
      call set_ASk(i_basin,Ask)
      call set_DN(i_basin,DN)
      call set_hcID(i_basin,hcID)
      call set_A_inflow(i_basin_cv_1,A_inflow_cv_1)
      call set_A_inflow(i_basin_cv_2,A_inflow_cv_2)
      call set_sumdev_num(i_basin_cv_1,sumdev_num_cv_1)
      call set_sumdev_num(i_basin_cv_2,sumdev_num_cv_2)
      if (abs(w_1).gt.epsilon(w_1).and.dtm_hcID(i_basin_cv_1).eq.0) then
         call set_hcID(i_basin_cv_1,hcID)
      end if
      if (abs(w_2).gt.epsilon(w_2).and.dtm_hcID(i_basin_cv_2).eq.0) then
         call set_hcID(i_basin_cv_2,hcID)
      end if

     ! Strahler Stream Order definition

	 if (ndcf == 1) then
	    if (dtm_A_inflow(i_basin) == 0.0E0) then
           call set_sso(i_basin,1)
           call set_hso(i_basin,1)
	    end if
        if (w_1 == 1) then
		   sso=dtm_sso(i_basin)
           sso_cv=dtm_sso(i_basin_cv_1)
		   sso_jun_cv=dtm_sso_jun(i_basin_cv_1)
		   if (sso_cv < sso) then
		      sso_cv=sso
              sso_jun_cv=0
		   else if (sso_cv == sso) then
		      if (sso_jun_cv == 0) then
			     sso_cv=sso_cv+1
                 sso_jun_cv=1
              end if
		   end if		   
		   call set_sso_jun(i_basin_cv_1,sso_jun_cv)
		   call set_sso(i_basin_cv_1,sso_cv)           
		   call set_hso(i_basin_cv_1,sso_cv)
		end if
        if (w_2 == 1) then
		   sso=dtm_sso(i_basin)
           sso_cv=dtm_sso(i_basin_cv_2)
		   sso_jun_cv=dtm_sso_jun(i_basin_cv_2)
		   if (sso_cv < sso) then
		      sso_cv=sso
              sso_jun_cv=0
		   else if (sso_cv == sso) then
		      if (sso_jun_cv == 0) then
			     sso_cv=sso_cv+1
                 sso_jun_cv=1
              end if
		   end if
		   call set_sso_jun(i_basin_cv_2,sso_jun_cv)
		   call set_sso(i_basin_cv_2,sso_cv)
           call set_hso(i_basin_cv_2,sso_cv)
		end if     
	 end if
	   
   else if (abs(s_max).lt.epsilon(s_max)) then

      write(6,*) 's_max = 0 at i_basin =',i_basin

      ! load information calculated in previous cycles

      A_inflow=dtm_A_inflow(i_basin)
      A_outflow=A_inflow+A_cell
      sumdev_num=dtm_sumdev_num(i_basin)
      if (A_inflow.eq.0.0_REP) then
         sumdev=0.0_REP
      else
         sumdev=sumdev_num/A_inflow
      end if 
      Kp=dtm_DN(i_basin)

      ! determination of the output cell
 
      emin=huge(emin)
      p_outflow_L=0
      do l=1,9
         if (e(l).ne.0.0_REP.and.l.ne.5.and.e(l).lt.emin) then
            emin=e(l)
            p_outflow_L=l
         end if
      end do
      if (p_outflow_L.eq.0) then
         write(6,*) 's_max = 0, unexpected case!'
         write(6,*) 'i =',i
         write(6,*) 'j =',j
         stop
      end if

      ! calculate local phisiographic features of the current cell

      if (mod(p_outflow_L,2).eq.0) then
         p_outflow_1=p_outflow_L
         epl_1=delta_x
         local_slope_1=(e0-emin)/epl_1
         ASk=A_outflow*local_slope_1**kas
         DN=Kp/(-mean_s_max)
         if (hcID.eq.0) hcID=ChannelInitiation(A_outflow,ASk,DN)
         sumdev_1=lambda*sumdev
         w_1=1.0_REP
         ! calculate features of the following cell
         ! the function p_outflow is the coordinate z of the plane ax+by+c=z,
         ! with a=3, b=1 and c=5. the function (i_basin_cv-i_basin) is the
         ! coordinate z of the plane ax+by+c=z, with a=M, b=1 e c=0.
         iii_1=nint((p_outflow_1-5)/3.)
         jjj_1=p_outflow_1-5-3*iii_1
         i_basin_cv_1=i_basin+(M*iii_1+jjj_1)
         A_inflow_cv_1=dtm_A_inflow(i_basin_cv_1)
         sumdev_num_cv_1=dtm_sumdev_num(i_basin_cv_1)
         A_inflow_cv_1=A_inflow_cv_1+(A_outflow)*w_1      
         sumdev_num_cv_1=sumdev_num_cv_1+(A_outflow)*w_1*sumdev_1
         call set_w_1(i_basin,w_1)
         call set_p_outflow_1(i_basin,p_outflow_1)
         call set_epl_1(i_basin,epl_1)
         call set_local_slope_1(i_basin,local_slope_1)
         call set_A_inflow(i_basin_cv_1,A_inflow_cv_1)
         if (dtm_hcID(i_basin_cv_1).eq.0) call set_hcID(i_basin_cv_1,hcID)
      else  ! if (mod(p_outflow_L,2).eq.1) then
         p_outflow_2=p_outflow_L
         epl_2=delta_xy
         local_slope_2=(e0-emin)/epl_2
         ASk=A_outflow*local_slope_2**kas
         DN=Kp/(-mean_s_max)
         sumdev_2=lambda*sumdev
         w_2=1.0_REP
         ! calculate features of the following cell
         ! the function p_outflow is the coordinate z of the plane ax+by+c=z,
         ! with a=3, b=1 and c=5. the function (i_basin_cv-i_basin) is the
         ! coordinate z of the plane ax+by+c=z, with a=M, b=1 e c=0.
         iii_2=nint((p_outflow_2-5)/3.)
         jjj_2=p_outflow_2-5-3*iii_2
         i_basin_cv_2=i_basin+(M*iii_2+jjj_2)
         A_inflow_cv_2=dtm_A_inflow(i_basin_cv_2)
         sumdev_num_cv_2=dtm_sumdev_num(i_basin_cv_2)
         A_inflow_cv_2=A_inflow_cv_2+(A_outflow)*w_2      
         sumdev_num_cv_2=sumdev_num_cv_2+(A_outflow)*w_2*sumdev_2
         call set_w_2(i_basin,w_2)
         call set_p_outflow_2(i_basin,p_outflow_2)
         call set_epl_2(i_basin,epl_2)
         call set_local_slope_2(i_basin,local_slope_2)
         call set_A_inflow(i_basin_cv_2,A_inflow_cv_2)
         call set_sumdev_num(i_basin_cv_2,sumdev_num_cv_2)
         if (dtm_hcID(i_basin_cv_2).eq.0) call set_hcID(i_basin_cv_2,hcID)
      end if 
      call set_ASk(i_basin,Ask)
      call set_DN(i_basin,DN)
      call set_hcID(i_basin,hcID)

   else ! if (s_max.lt.0.0_REP) then

      write(6,*) 's_max < 0 at i_basin =',i_basin
      stop

   end if

end do

! for the outlet cell a phantom channel end is set up

nvo=0
A_inflow_max=0.0D+00
do ii=i-1,i+1
   do jj=j-1,j+1
      if (ii.eq.0.or.ii.eq.N+1) cycle
      if (jj.eq.0.or.jj.eq.M+1) cycle
      ii_basin=(ii-1)*M+jj
      if (dtm_index_pr(ii_basin) == 0) cycle
      p_outflow_1_ciijj=dtm_p_outflow_1(ii_basin)
      p_outflow_2_ciijj=dtm_p_outflow_2(ii_basin)
      p_inflow=3*(ii-i)+(jj-j)+5
      if (p_inflow+p_outflow_1_ciijj.eq.10) then
         A_outflow_ciijj=(dtm_A_inflow(ii_basin)+A_cell)*dtm_w_1(ii_basin)
         if (A_outflow_ciijj.gt.A_inflow_max) then
            A_inflow_max=A_outflow_ciijj
            p_outflow_outlet=p_outflow_1_ciijj
            local_slope_outlet=dtm_local_slope_1(ii_basin)
            ivo=i+nint((p_outflow_outlet-5)/3.)
            jvo=j+p_outflow_outlet-5-3*(ivo-i)
            ivo_basin=i_basin+(M*(ivo-i)+(jvo-j))
            if (ivo.eq.0.or.ivo.eq.N+1.or.jvo.eq.0.or.jvo.eq.M+1.or. &
                dtm_index_pr(ivo_basin).eq.0) then
               nvo=1
            else
               nvo=0
            end if
         end if
      end if
      if (p_inflow+p_outflow_2_ciijj.eq.10) then
         A_outflow_ciijj=(dtm_A_inflow(ii_basin)+A_cell)*dtm_w_2(ii_basin)
         if (A_outflow_ciijj.gt.A_inflow_max) then
            A_inflow_max=A_outflow_ciijj
            p_outflow_outlet=p_outflow_2_ciijj
            local_slope_outlet=dtm_local_slope_2(ii_basin)
            ivo=i+nint((p_outflow_outlet-5)/3.)
            jvo=j+p_outflow_outlet-5-3*(ivo-i)
            ivo_basin=i_basin+(M*(ivo-i)+(jvo-j))
            if (ivo.eq.0.or.ivo.eq.N+1.or.jvo.eq.0.or.jvo.eq.M+1.or. &
                dtm_index_pr(ivo_basin).eq.0) then
               nvo=1
            else
               nvo=0
            end if
         end if
      end if
   end do
end do

if (nvo.eq.0) then
   p_outflow_outlet=p_outflow_vo 
   write(6,*) 'the drainage direction of the outlet cell (', &
               p_outflow_outlet,') defined in hap.in is used'
else
   write(6,*) 'the drainage direction of the outlet cell (', &
               p_outflow_outlet,') is used'
end if

if (mod(p_outflow_outlet,2).eq.0) then ! cardinal dainage direction
   p_outflow_1=p_outflow_outlet
   w_1=1.0
   epl_1=delta_x 
   local_slope_1=local_slope_outlet
   call set_p_outflow_1(i_basin,p_outflow_1)
   call set_w_1(i_basin,w_1)   
   call set_epl_1(i_basin,epl_1)
   call set_local_slope_1(i_basin,local_slope_1)
else ! diagonal dainage direction
   p_outflow_2=p_outflow_outlet
   w_2=1.0
   epl_2=delta_xy 
   local_slope_2=local_slope_outlet
   call set_p_outflow_2(i_basin,p_outflow_2)
   call set_w_2(i_basin,w_2)   
   call set_epl_2(i_basin,epl_2)
   call set_local_slope_2(i_basin,local_slope_2)
end if

! Horton Stream Order Definition

if (ndcf == 1) then  
  do n_o_quota=N_celle,1,-1   
     A_inflow_max=0.0D+00
     read(16,rec=n_o_quota+1) i_basin
     jr=mod(i_basin,M)
     if (jr.ne.0) then
        j=jr
        i=(i_basin-j)/M+1
     else
        j=M
        i=i_basin/M
     end if
     do ii=i-1,i+1
        do jj=j-1,j+1
          if (ii.eq.0.or.ii.eq.N+1) cycle
          if (jj.eq.0.or.jj.eq.M+1) cycle
          ii_basin=(ii-1)*M+jj
          if (dtm_index_pr(ii_basin) == 0) cycle
          p_outflow_1_ciijj=dtm_p_outflow_1(ii_basin)*dtm_w_1(ii_basin)
          p_outflow_2_ciijj=dtm_p_outflow_2(ii_basin)*dtm_w_2(ii_basin)
          p_inflow=3*(ii-i)+(jj-j)+5
          if (p_inflow+p_outflow_1_ciijj.eq.10) then
		     A_outflow_ciijj=(dtm_A_inflow(ii_basin)+A_cell)*dtm_w_1(ii_basin)             
             if (A_outflow_ciijj > A_inflow_max) then
                A_inflow_max=A_outflow_ciijj
		        ii_basin_max=ii_basin
		  	 end if
          end if
		  if (p_inflow+p_outflow_2_ciijj.eq.10) then
             A_outflow_ciijj=(dtm_A_inflow(ii_basin)+A_cell)*dtm_w_2(ii_basin)
             if (A_outflow_ciijj > A_inflow_max) then
                A_inflow_max=A_outflow_ciijj
		        ii_basin_max=ii_basin
		  	 end if
          end if
        end do
     end do
     if (A_inflow_max /= 0.0D0) then
	    call set_hso(ii_basin_max,dtm_hso(i_basin))
	 end if
  end do  !do n_o_quota=N_celle,1
end if ! if (ndcf == 1) then

close(16)

return
end subroutine dsf


function ChannelInitiation(A_outflow,ASk,DN)

use mpar
implicit none
integer(kind=ISP) :: ChannelInitiation
real(kind=REP),intent(in) :: A_outflow
real(kind=RSP),intent(in) :: ASk,DN

select case (nchc)
case (1)
   if (A_outflow.le.A_threshold) then
      ChannelInitiation=0
   else
      ChannelInitiation=1
   end if
case (2)
   if (ASk.le.ASk_threshold) then
      ChannelInitiation=0
   else
      ChannelInitiation=1
   end if
case (3)
   if (DN.ge.DN_threshold) then
      ChannelInitiation=0
   else
      ChannelInitiation=1
   end if
case default
   write(6,*) 'nchc out of range!'
   stop
end select

end function ChannelInitiation
