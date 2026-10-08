c
c********************** ORDN_V ****************************
c
c  ordinamento crescente di un  vettore di interi e in maniera
c  coerente di un vettore di reali
c
c*********************************************************
      subroutine ordn_v(lvec,ivec,vec)

      implicit none

      integer lvec,itemp
      integer ivec(lvec),j,k
      real*8 vec(lvec),rtemp


      do j = 1,lvec-1
         do k = j+1,lvec
            if(ivec(j).gt.ivec(k)) then
               itemp = ivec(j)
               ivec(j) = ivec(k)
               ivec(k) = itemp  
               rtemp = vec(j)
               vec(j) = vec(k)
               vec(k) = rtemp
            end if
         end do
      end do

      return
      end
