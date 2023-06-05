C ******************** CFLPECNUMBER ****************************************
C  to calculate the CFL number  CFL = Deltat * sup(surface/volume) * sup(dF/du)
C  where F is the function of the advection equation u_t + div(F(u))=0
C  In this case F is a linear function F(u)= v*u
C  and to calculate the Peclet number Pec = CFL/gamma where gamma is the
C  diffusion number gamma=Deltat*sup(|D|/volume)
C  supdiffus is equal to sup (|D|/volume) 
C  **********************************************************************
C    
      subroutine cflpecnumber(NTETRA,DELTATADV,VOLUR,SURFAR,
     1                      VX,VY,VZ,SUPDIFFUS,cflinp,CFLNUMB,PECNUMB,
     2                      CFLTETRA,NADV,nadvfl,DELTAT)

      implicit none
C
      integer ntetra
      integer j, cfltetra, nadv,nadvfl
      real*8 deltat,deltatadv
      real*8 volur(ntetra),surfar(ntetra)
      real*8 vx(ntetra),vy(ntetra),vz(ntetra)
      real*8 cflinp,cflnumb,pecnumb,supdiffus
      real*8 supcfl,supcflt
      
      cfltetra=1
      supcfl= surfar(1)*volur(1)*
     1        dsqrt(vx(1)**2+vy(1)**2+vz(1)**2)
c     write(999,*)surfar(1),volur(1)
      do j=2,ntetra
         supcflt= surfar(j)*volur(j)*
     1            dsqrt(vx(j)**2+vy(j)**2+vz(j)**2)
         if (supcflt.gt.supcfl) then
            supcfl=supcflt
            cfltetra=j
         end if   
c     write(999,*)surfar(j),volur(j)
      end do
      if (nadvfl.eq.1) then
         deltatadv=deltat/float(nadv)
         cflnumb=deltatadv*supcfl
      else if(nadvfl.eq.0) then
         deltatadv= cflinp/supcfl
         nadv=max(int(deltat/deltatadv)+1,2)
         deltatadv= deltat/nadv
         cflnumb=deltatadv*supcfl
      else
         write(6,*) 'wrong value for nadvfl ',nadvfl
         stop
      end if
      if (supdiffus.eq.0.0d0) then
         pecnumb=1.e+20
      else
         pecnumb=cflnumb/(supdiffus*deltatadv)
      end if
         
      return 
      end
      
