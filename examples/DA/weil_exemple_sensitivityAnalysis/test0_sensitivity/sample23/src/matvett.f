      subroutine matvett(nmax,NN,NTERM,IAT,IAC,COEF,x,r)
      implicit none 
      integer nmax,nterm
      integer i,j,IAC(nmax),IAT(nterm),NN,ind
      real*8 COEF(nterm),x(nmax),r(nmax)
C     prodotto matrice-vettore nel caso di matrice simmetrica scritta
C     in forma compatta  Ax=r
C     COEF coefficienti non nulli della matrice simmetrica A
C     IAT indice di topologia di A
C     IAC indice di colonne di A
C     IAT(i).LE.j.LT.IAT(i+1) -- COEF(j)=A(i,IAC(j)) 
      do i=1,NN
         r(i)=0.d0
      end do
      do i=1,NN
         ind=IAT(i)
         r(i)=r(i)+COEF(ind)*x(IAC(ind))
         do j=ind+1,IAT(i+1)-1
            r(i)=r(i)+COEF(j)*x(IAC(j))
            r(IAC(j))=r(IAC(j))+COEF(j)*x(i)
         end do
      end do
      end
