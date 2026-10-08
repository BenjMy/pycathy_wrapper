      subroutine godunov(ntetra,nface,npfa_tra,side_cnc,
     1                   puntdir,plist,cold_el,prescfa_tra,precl,precr)
      implicit none
      
      integer ntetra,nface,npfa_tra
      integer side_cnc(4,*)
      integer puntdir(*),plist(2,*)

      integer j,itetra 
      integer iface,itetra2
     
      real*8  ubar,vel
      real*8  cold_el(*),prescfa_tra(*)
      real*8  precl(*),precr(nface)

      do itetra = 1,ntetra
         ubar = cold_el(itetra)
        do j=1,4
             iface=side_cnc(j,itetra)
                 if (plist(1,iface).eq.itetra) then
                    precl(iface) = ubar 
c
c on a boundary face we set also the `right' reconstructed value
c
                    if (plist(2,iface).eq.0) then
                        if(puntdir(iface).eq.0) then
c                             precr(iface)= ubar 
c if it is not a Dirichlet boiundary face, the precr value is 0.0
                             precr(iface)= 0.0d0 
                        else
                             precr(iface)= prescfa_tra(puntdir(iface))
                        end if
                    end if
                 else
                       precr(iface)=ubar
                 end if
              end do
      end do
      return
      end

