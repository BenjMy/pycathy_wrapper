C
C**************************   CALCOLI   ********************************
C attenta: il contributo diretto di pioggia sul serbatoio e' stato
C inserito nel calcolo di Q_inflow, aggiungendo surface_water
C non moltiplico piu' per l'area di cella!
C***********************************************************************
C
      subroutine calcoli(hh,dd,delta_t,n_ha,Ah_r,hA_r,Q_in_kk,t,
     1                   h_sfioro, Q_in_kkp1,Qh_r_p,Ah_r_p,Q_inflow_p,
     2                   surface_water,h_fondo_rr)

      implicit none
      include 'CATHY.H'
      include 'IOUNITS.H' 

      integer j,n_hA

      real*8  Ah_r_p,hh,Q_inflow_p,Qh_r_p,Q_in_kk
      real*8  Q_in_kkp1,delta_t,t,dd
      real*8  hA_r(n_hA),Ah_r(n_hA)
      real*8 h_sfioro,surface_water,h_fondo_rr

c     interpolo Q entrante

       Q_inflow_p=Q_in_kk+(Q_in_kkp1-Q_in_kk)/delta_t*t

      Q_inflow_p=Q_inflow_p+ surface_water

      if(hh.eq.h_fondo_rr .and. Q_inflow_p.lt.0.0d0) then
         dd=0.0d0
         write(6,*) 'hh=fondo serb=',hh,'Q_in=',Q_inflow_p,'dd=0'
         go to 1000
      end if

c     interpolo A

      do j=1,n_hA-1

         if ((hA_r(j).lt.hh).and.(hA_r(j+1).ge.hh)) then

           Ah_r_p=Ah_r(j)+(Ah_r(j+1)-Ah_r(j))/
     &            (hA_r(j+1)-hA_r(j))*(hh-hA_r(j))
           go to 200
         end if
 
      end do

      if (hh.gt.hA_r(n_hA)) then
          write(6,*) 'fuori tab. in alto!'
          Ah_r_p=Ah_r(n_hA)
          go to 200
      end if

      write(6,*) 'hA fuori campo'
      stop
200   continue

c
c     calcolo Q uscente
c
      if(hh.le.h_sfioro) then
        Qh_r_p=0.0d0
      else
       if(Q_inflow_p.ge.0.0d0) then
          Qh_r_p=Q_inflow_p
       else
cxcx      write(iout41,*) 'attenzione livello serb superiore sfioro'
          Qh_r_p=0.0d0
       end if
      end if

300   continue
      dd=1/Ah_r_p*(Q_inflow_p-Qh_r_p)

1000  continue

      return
      end            
