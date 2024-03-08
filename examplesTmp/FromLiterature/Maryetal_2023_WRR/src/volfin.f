C
C**************************  VOLFIN  ***********************************
C
C      implements Finite Volume method at each time-step 
C    
C***********************************************************************
C
      subroutine volfin(nnode,ntria,ntetra,nface,npfa_tra,n1,nstr,nzone,
     1            flag_temp,flag_interp,flag_limiter,
     2            plist,puntdir,side_cnc,iside,neigh,neigbc,
     3            tetravert,tvert,tetra,ip4,x,y,z,xc,yc,zc,
     4            deltatadv,volur,volu,vn,swpint_el,swint_el,
     5            poros,cold_el,
     6            cnew_el,cdiff,prescfa_tra,flux,precl,precr,time)
C
      implicit none
      INCLUDE 'CATHY.H'
      integer nnode, ntria,ntetra,nface,nstr,npfa_tra,n1,nzone
      integer flag_interp,flag_limiter,ip4(4,4),flag_temp
      integer plist(2,*),side_cnc(4,*),puntdir(*)
      integer neigh(4,*),neigbc(3,*)
      integer iside(3,*),tetravert(NMAX,*)
      integer tvert(*),tetra(5,*)
      integer istrato,istr,izone,itetra
c
      real*8  tot
      real*8  poros(MAXSTR,*),x(*),y(*),z(*)
      real*8  xc(*),yc(*),zc(*)
      real*8  volur(*),volu(*),deltatadv
      real*8  vn(*),cold_el(*),cnew_el(*),cdiff(*)
      real*8  prescfa_tra(*),flux(*),precl(*)
      real*8  precr(*)
      real*8  swpint_el(*),swint_el(*)
      real*8  massnew,massold
      integer NTLAYER,j,i,iface
      real*8  time,c,a,b
      external istrato
C  flag_temp = 0  --->  Euler scheme in time and Godunov reconstruction
C  flag_temp =1   --->  Euler scheme in time and 2nd order reconstruction
C  flag_temp=2    ---> correction term added in midpoint(Euler) scheme in time
C  flag_temp=3    ---> Mid-point scheme in time, Runge scheme
C  flag_temp=4    --->  Heun scheme in time
C -----------------------------------------------------------------------
C flag_interp .eq. 1 -->  construction of the interpolant as in Durlofsky 
C ------------------      scheme taking into account the 4 possible interpolant
C                        that can be constructed from the refering tetrahedron
C    flag_limiter = 0        Durloflsky limiter
C    flag_limiter = 1        min-limiter
C    flag_limiter = 2        the maximum  gradient interpolant is chosen with 
C                            no limiter
C    flag_limiter =3         the average of the 4 interpolant is applied with 
C                            no limiter
C flag_interp .eq.2  -->  construction of the interpolant considering the 4
C -----------------       neighbours tetrahedra
C    flag_limiter =1         extremum limiter (see paper: compressible large 
C                            eddy simulation using unstructured grid: 
C                            supersonic boundary layer and compression ramps,
C                            Yan, Urbin, Knight, 10th Intern. Conf. on Methods 
C                            of Aerophysical Research,July 9-16 2000, 
C                            Novosibirsk, Russia
C    flag_limiter =2         Limited Central Difference LCD limiter  (see paper 
C                            Multidimensional Slope Limiters for MUSCL-type
C                            finite volume schemes on unstructured grids,
C                             M.E. Hubbard, J. Comp. Phys. 155, pp.54-74, 1999)
C    flag_limiter =3         no limiter is applied
C    flag_limiter=4          Barth-Jespersen limiter (see paper: The design and
C                            application of upwind schemes on unstructerd 
C                            meshes, T. J. Barth, D. C. Jespersen, AIAA-89-0366)
C flag_interp.eq.3  --->  MLG scheme: the interpolant is chosen between 
C ----------------        all the interpolant that can be constructed both
C                        with flag_interp.eq.1. and flag_interp.eq.2
C                        Between them, the 'proper' interpolant is limited
C                        with the LCD limiter (see above)
C flag_interp.eq.4  --->  MLS interpolant: the interpolant is chosen by the
C ----------------       least squares method
C     flag_limiter =1          extremum limiter
C     flag_limiter=2           LCD limiter
C     flag_limiter =3          no limiter
C     flag_limiter=4           Barth-Jespersen limiter 
C  -------------------------------------------------------

      call init0r(ntemax,flux)
      call init0r(nfacemax,precl)
      call init0r(nfacemax,precr)
      if (flag_temp.eq.0) then
C
C      godunov reconstruction and Euler in time
C
C     WITH THIS SUBROUTINE THE VALUE OF CONCENTRATION AT THE RIGHT (precr) 
C     AND LEFT (precl) SIDE OF EACH FACE IS SET  
          call godunov(ntetra,nface,npfa_tra,side_cnc,  
     1                 puntdir,plist,cold_el,prescfa_tra,precl,precr)
C     WITH THIS SUBROUTINE IT IS CALCULATING THE FLUX (flux) GOING OUT/IN
C     FOR EACH ELEMENT       
          call durlo(ntetra,nface,volur,plist,vn,precl,precr,time,flux)
C
C     end of godunov reconstruction
      
C  
      elseif ((flag_temp.eq.1).or.(flag_temp.eq.2)
     1       .or.(flag_temp.eq.10)) then
C
C    Durlofsky reconstruction
C    CARLOTTA NON SONO D'ACCORDO
C    flag_temp.eq.1 ==> Dirichlet BC's calculated at time t^{k+1/2} (mid-point
C                       in time)
C    flag_temp.eq.2 ==> correction term to obtain second order when applying
C                       the time-splitting technique
C    flag_temp.eq.10 ==> Euler in time
       if (flag_interp.eq.1) then
          call   tvd_durlo(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cold_el,cdiff,prescfa_tra,
     5                     precl,precr)
       elseif ((flag_interp.eq.2).or.(flag_interp.eq.4)) then
          call   extr_tvd(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_interp,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cold_el,cdiff,prescfa_tra,
     5                     precl,precr)
       elseif (flag_interp.eq.3) then
          call   mlg_limiter(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cold_el,cdiff,prescfa_tra,
     5                     precl,precr)
       end if
          call durlo(ntetra,nface,volur,plist,vn,precl,precr,time,flux)
      elseif (flag_temp.eq.3) then
C
C  2nd order Runge method:
C 
       if (flag_interp.eq.1) then
          call   tvd_durlo(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cold_el,cdiff,prescfa_tra,
     5                     precl,precr)
       elseif ((flag_interp.eq.2).or.(flag_interp.eq.4)) then
          call extr_tvd(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_interp,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cold_el,cdiff,prescfa_tra,
     5                     precl,precr)
       elseif (flag_interp.eq.3) then
          call mlg_limiter(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cold_el,cdiff,prescfa_tra,
     5                     precl,precr)
       end if 
          call durlo(ntetra,nface,volur,plist,vn,precl,precr,time,flux)

CC    THE NEW CONCENTRATION OF EACH ELEMENT IS PERFORMED BY COMPUTING
CC    A MASS BALANCE OVER IT BETWEEN THE OLD AND NEW ADVECTOVE TIME 
CC    LEVEL

          do itetra=1,ntetra
            istr=istrato(itetra,ntria)
            izone=tetra(5,itetra)
            massold = cold_el(itetra)*swpint_el(itetra)
     1                *poros(istr,izone)
c            massnew = massold + 0.5d0*deltatadv*flux(itetra)
            massnew = massold + deltatadv*flux(itetra)
            cnew_el(itetra)=massnew/(swint_el(itetra)*poros(istr,izone))
         end do
       if (flag_interp.eq.1) then
          call tvd_durlo(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cnew_el,cdiff,prescfa_tra,
     5                     precl,precr)
       elseif ((flag_interp.eq.2).or.(flag_interp.eq.4)) then
          call extr_tvd(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_interp,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cnew_el,cdiff,prescfa_tra,
     5                     precl,precr)
       elseif (flag_interp.eq.3) then
          call   mlg_limiter(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cnew_el,cdiff,prescfa_tra,
     5                     precl,precr)
       end if 
          call durlo(ntetra,nface,volur,plist,vn,precl,precr,time,flux)

      else if (flag_temp.eq.4) then
C 2nd order Runge-Kutta method (Heun method)
       if (flag_interp.eq.1) then
          call   tvd_durlo(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cold_el,cdiff,prescfa_tra,
     5                     precl,precr)
       elseif ((flag_interp.eq.2).or.(flag_interp.eq.4)) then
          call   extr_tvd(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_interp,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cold_el,cdiff,prescfa_tra,
     5                     precl,precr)
       elseif (flag_interp.eq.3) then
          call mlg_limiter(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cold_el,cdiff,prescfa_tra,
     5                     precl,precr)
       end if
          call durlo(ntetra,nface,volur,plist,vn,precl,precr,time,flux)
C CARLOTTA Finish first step Heun scheme: I am calculating c^{ka+1} read
C in a second time as input for the second step of the Heun scheme
C (which updates c^{ka+1})
          do itetra=1,ntetra
            istr=istrato(itetra,ntria)
            izone=tetra(5,itetra)
            massold = cold_el(itetra)*swpint_el(itetra)
     1                *poros(istr,izone)
            massnew = massold + deltatadv*flux(itetra)
            cnew_el(itetra)=massnew/(swint_el(itetra)*poros(istr,izone))
          end do
          tot=0.0d0 
          do itetra=1,ntetra
             istr=istrato(itetra,ntria)
             izone=tetra(5,itetra)
             tot=tot+swint_el(itetra)*poros(istr,izone)*
     1       cnew_el(itetra)*volu(itetra)
          end do 
c          write(*,*)'BILANCIO_PRIMO',time,tot
       if (flag_interp.eq.1) then
          call tvd_durlo(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cnew_el,cdiff,prescfa_tra,
     5                     precl,precr)
       elseif ((flag_interp.eq.2).or.(flag_interp.eq.4)) then
cc CARLOTTA For second step Heun scheme 
          call extr_tvd(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_limiter,flag_interp,flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cnew_el,cdiff,prescfa_tra,
     5                     precl,precr)

       elseif (flag_interp.eq.3) then
          call   mlg_limiter(nnode,ntetra,nface,npfa_tra,n1,
     1                     flag_temp,neigh,
     2                     side_cnc,puntdir,plist,iside,ip4,
     3                     tetravert,tetra,
     4                     x,y,z,xc,yc,zc,cnew_el,cdiff,prescfa_tra,
     5                     precl,precr)
       end if
          call durlo(ntetra,nface,volur,plist,vn,precl,precr,time,flux)

cc CARLOTTA Second step Heun scheme 

        do itetra=1,ntetra
           istr=istrato(itetra,ntria)
           izone=tetra(5,itetra)
           massold = cnew_el(itetra)*swint_el(itetra)*poros(istr,izone)
           massnew = massold + deltatadv*flux(itetra)
           massnew= 0.5d0*(massnew + cold_el(itetra)*swpint_el(itetra)*
     1              poros(istr,izone))
           cnew_el(itetra)=massnew/(swint_el(itetra)*poros(istr,izone))
        end do
          tot=0.0d0 
          do itetra=1,ntetra
             istr=istrato(itetra,ntria)
             izone=tetra(5,itetra)
             tot=tot+swint_el(itetra)*poros(istr,izone)*
     1       cnew_el(itetra)*volu(itetra)
          end do 
c          write(*,*)'BILANCIO_SECONDO',time,tot
      end if
C
C     application of predictor on time integration
C 
      if (flag_temp.ne.4) then
        do itetra=1,ntetra
           istr=istrato(itetra,ntria)
           izone=tetra(5,itetra)
           massold = cold_el(itetra)*swpint_el(itetra)*poros(istr,izone)
           massnew = massold + deltatadv*flux(itetra)
           cnew_el(itetra)=massnew/(swint_el(itetra)*poros(istr,izone))
        end do
       end if
 
      return

      end
