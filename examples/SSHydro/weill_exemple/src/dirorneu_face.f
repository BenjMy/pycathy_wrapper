c
c************************* DIRNEU_FACE ******************************** 
c
c  individua le facce formate da tre nodi di Dirichlet (o di Neuman,
c  o di altro tipo); tali nodi sono  memorizzati nel vettore  
c  CONTD di dimensione ND , i valori prescritti sono in  PRESCD
c
c ==> TRADUCTION FROM SYLVAIN : convert the node based Dirichlet BC for 
C transport and setthem on faces
c**********************************************************************
      subroutine dirorneu_face(n3d,cont3d,presc3d,nface,plist,iside,
     1                      npq,contpq,prescpq)
      implicit none
      integer n3d,nface,npq,cont3d(*),contpq(*),iv(3),iface,i,j,k
      integer iside(3,*),plist(2,*)
      real*8 presc3d(*),prescpq(*)
      

      call ordn_v(n3d,cont3d,presc3d)
      npq=0

      do iface=1,nface
         if(plist(2,iface).eq.0) then
            iv(1)=iside(1,iface)
            iv(2)=iside(2,iface)
            iv(3)=iside(3,iface)
            call ordn(3,iv)
            do i=1,n3d
               if(cont3d(i).eq.iv(1)) then
                  do j=i+1,n3d
                     if(cont3d(j).eq.iv(2)) then
                        do k=j+1,n3d
                           if(cont3d(k).eq.iv(3)) then
                              npq=npq+1
                              contpq(npq)=iface
c uno dei tre nodi supponendo che il valore prescritto sui tre
c nodi della faccia sia uguale
c                  prescpq(npq)=presc3d(i)
c  average of the values along the three nodes
                              prescpq(npq)=
     1                    (presc3d(i) + presc3d(j) + presc3d(k))/3.0d0
                           endif
                        end do
                     endif
                  end do
               endif
            end do
         endif
      end do
         
      return
      end
