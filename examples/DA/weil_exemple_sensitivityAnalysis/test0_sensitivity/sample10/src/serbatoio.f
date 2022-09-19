C
C**************************  SERBATOIO  ********************************
C
C
C***********************************************************************
C

      subroutine serbatoio(index_rr,i_rr,j_rr,tipo_rr,Q_in_kk,
     1                     Q_in_kkp1,Q_out_kkp1,h_pool_kkp1,
     2                     h_pool_kk,delta_t,surface_water,
     3                     num_rr,n_hA_vec,hres,Ares,h_soglia,
     4                     h_fondo_rr)

      implicit none
      include 'CATHY.H'
      include 'IOUNITS.H'

      integer iin
      integer n_ha,i
      integer index_rr,i_rr,j_rr,tipo_rr,num_rr
      integer n_hA_vec(*)
      real*8  h_soglia(*),hres(maxres,*),Ares(maxres,*)
      real*8  h_fondo_rr
      real*8 t,delta_t,dd,delta_h1,delta_h2,delta_h3
      real*8 delta_h4,hh,Q_in_kk,Q_in_kkp1,t_iniz,t_fin
      real*8 h_pool_kk,h_pool_kkp1
      real*8 hA_r(30),Ah_r(30)
      real*8 Q_out_kkp1,Q_inflow_p,Qh_r_p,Ah_r_p
      real*8 h_sfioro,surface_water

      n_hA= n_hA_vec(num_rr)

      do i=1,n_hA
        hA_r(i)=hres(num_rr,i)
        Ah_r(i)=Ares(num_rr,i)
      end do

      h_sfioro=h_soglia(num_rr)

      hh=h_pool_kk

      t_iniz=0.0d0
      t_fin=delta_t

      t=t_iniz

      call calcoli(hh,dd,delta_t,n_ha,Ah_r,hA_r,Q_in_kk,t,
     1             h_sfioro,Q_in_kkp1,Qh_r_p,Ah_r_p,Q_inflow_p,
     2             surface_water,h_fondo_rr)                                      

      delta_h1=delta_t*dd

      t=t_iniz+0.50d0*delta_t

      hh=h_pool_kk+0.50d0*delta_h1

      if(hh.lt.h_fondo_rr) then
         hh=h_fondo_rr
         write(6,*) 'errore tipo aaaa'
      end if
    
      call calcoli(hh,dd,delta_t,n_ha,Ah_r,hA_r,Q_in_kk,t,
     1             h_sfioro,Q_in_kkp1,Qh_r_p,Ah_r_p,Q_inflow_p,
     2             surface_water,h_fondo_rr)

      delta_h2=delta_t*dd

      t=t_iniz+0.50d0*delta_t
 
      hh=h_pool_kk+0.50d0*delta_h2

      if(hh.lt.h_fondo_rr) then
         hh=h_fondo_rr
         write(6,*) 'errore tipo aaaa'
      end if

      call calcoli(hh,dd,delta_t,n_ha,Ah_r,hA_r,Q_in_kk,t,
     1             h_sfioro,Q_in_kkp1,Qh_r_p,Ah_r_p,Q_inflow_p,
     2             surface_water,h_fondo_rr)

      delta_h3=delta_t*dd

      t=t_fin
 
      hh=h_pool_kk+0.50d0*delta_h3

      if(hh.lt.h_fondo_rr) then
         hh=h_fondo_rr
         write(6,*) 'errore tipo aaaa'
      end if

      call calcoli(hh,dd,delta_t,n_ha,Ah_r,hA_r,Q_in_kk,t,
     1             h_sfioro,Q_in_kkp1,Qh_r_p,Ah_r_p,Q_inflow_p, 
     2             surface_water,h_fondo_rr)

      delta_h4=delta_t*dd

      h_pool_kkp1=h_pool_kk+delta_h1/6+delta_h2/3+delta_h3/3+delta_h4/6
      
      if(h_pool_kkp1.lt.h_fondo_rr) then
         h_pool_kkp1=h_fondo_rr
      end if

      if(h_pool_kkp1.le.h_sfioro) then
         Q_out_kkp1=0.0d0
      else
        if(Q_in_kkp1 + surface_water.ge.0.0d0) then
          Q_out_kkp1=Q_in_kkp1 + surface_water
        else
          Q_out_kkp1=0.0d0
        end if 
      end if


      return

      end


           
