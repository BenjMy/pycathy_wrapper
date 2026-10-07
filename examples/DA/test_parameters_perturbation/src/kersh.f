C
C************************* KERSH **********************************
C
      subroutine kersh(iout,nequ,nterm,ia,ja,sysmat,prec)
c
c  incomplete Choleski decompostion
c
      implicit none
c
      integer  iout,nequ,nterm
cxcx  integer  ia(nequ+1),JA(nterm)
      integer  ia(*),JA(*)

      integer  i,j,k,kk,k1,i1,j1,k2

cxcx  real*8   prec(nterm),sysmat(nterm)
      real*8   prec(*),sysmat(*)
      real*8   a,zero

      parameter (zero=0.0d0)

      do k=1,nterm
         prec(k) = zero
      end do

      do kk=1,nequ-1

         k = ia(kk)
         a = sysmat(k) - prec(k)
         if(a.le.zero) then
            write(iout,100) kk,a
            write(iout,101) prec(ia(kk-1))
            a = (prec(ia(kk-1)))**2
         end if
         prec(K) = sqrt(a)
 
         i = ia(kk) + 1
         j = ia(kk+1) - 1


         do k1 = i,j
            prec(k1) = (sysmat(k1)-prec(k1))/prec(k)
         end do

         do k2 = i,j-1

            j1 = ia(ja(k2))
            prec(j1) = prec(j1) + prec(k2)**2
            i1 = k2 + 1
            j1 = j1 + 1
            do while(j1.lt.ia(ja(k2)+1).and.i1.le.j)
               if(ja(j1).eq.ja(i1)) then
                  prec(j1) = prec(j1) + prec(k2)*prec(i1)
                  i1 = i1 + 1
                  j1 = j1 + 1
               else if (ja(j1).lt.ja(i1))then
                  j1 = j1 + 1
               else if (ja(j1).gt.ja(i1))then
                   i1 = i1 + 1
               end if
            end do

         end do

         if(j.ge.i)
     1        prec(ia(ja(j))) = prec(ia(ja(j))) + prec(k2)**2

      end do

      k = ia(nequ) 
      a = sysmat(k)-prec(k)
      if(a.le.zero) then
         write(iout,100) nequ ,a
         write(iout,101) prec(ia(nequ-1))
         a = (prec(ia(nequ-1)))**2
      end if
      prec(k) = sqrt(a)

      return
100   format('*** Subroutine Kersh: diagonal',
     1       ' element <= zero at position: ',I5,2X,E16.5)
101   format('***** using previous diagonal value: ',E16.8)
      end
