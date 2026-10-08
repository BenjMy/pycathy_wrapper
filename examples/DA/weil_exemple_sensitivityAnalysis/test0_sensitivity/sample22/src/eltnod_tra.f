C
C**************************  ELTNOD ************************************
C
C  average the values of a material or solute property assigned by
C  element to obtain the values at each node
C
C***********************************************************************
C
      SUBROUTINE ELTNOD_TRA(N,NT,TETRA,CNEW,CNNEW,PEL,PNODI,VOLU,
     1           VOLNOD,SWNEW,SW,NNOD)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   i,K,II,INOD,NNOD, j1, j2, j3, j4
      INTEGER   N,NT
      INTEGER   TETRA(5,*), cont(nmax)
      REAL*8    CNEW(*),CNNEW(*),volu(*),volnod(*),pel(*),pnodi(*)
      real*8    sw(*),swnew(*),mnnew(n),cnnew_stock(n)
      real*8    mass_test_node,mass_test_tetra
C
      CALL VCOPYR(N,cnnew_stock,CNNEW)
      CALL INIT0I(N,cont)
      CALL INIT0R(N,mnnew)
      CALL INIT0R(N,CNNEW)

      DO K=1,NT
         DO II=1,4
            INOD=TETRA(II,K)
c            if (inod.LE.2*nnod) THEN
            mnnew(INOD)=mnnew(INOD)+
     1                  CNEW(k)*pel(k)*SWNEW(k)*VOLU(k)/4
            cont(inod)=cont(inod)+1
c            else
c            cnnew(INOD)=cnnew_stock(INOD)
c            end if
         END DO
      END DO

c      DO i=1,2*NNOD
      DO i=1,n
         CNNEW(i)=mnnew(i)/(volnod(i)*pnodi(i)*sw(i))
      END DO
      
       do i = 1, nt
          j1=tetra(1,i)
          j2=tetra(2,i)
          j3=tetra(3,i)
          j4=tetra(4,i)
       if ((i.GE.3000) .AND. (i.LE.3010)) THEN
C            write(*,*) 'verif cnew, cold',i, cnew(i),
C     1              j1,cnnew(j1),j2,cnnew(j2),j3,cnnew(j3),j4,cnnew(j4)
        end if
        end do
      
      
c      do i=1,nt
c       mass_test_tetra=mass_test_tetra+
c     1                cnew(i)*volu(i)*pel(i)*swnew(i)
c      end do
c      
c
c      do i=1,n
c       mass_test_node=mass_test_node+
c     1                cnnew(i)*volnod(i)*pnodi(i)*sw(i)
c      end do
c
c      write(*,*) 'mass couche num1 tetra, node',
c     1       mass_test_tetra,mass_test_node



c            INOD=TETRA(II,K)
c            VNOD(INOD)=VNOD(INOD)+disnod(cont(inod),inod)*velt(k)
c            cont(inod)=cont(inod)+1
c         END DO
c      END DO

c      DO i=1,N
c         VNOD(i)=VNOD(i)/weightnod(i)
c      END DO
C
      RETURN
      END
