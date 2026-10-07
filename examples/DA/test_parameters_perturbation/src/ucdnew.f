C
C************************* UCDWRT *************************************
C
C  write in UCD format one scalar field over nodes of tetrahedra
C
C**********************************************************************
C
      SUBROUTINE UCDNEW(IOUT,N,NT,TETRA,X,Y,Z,POT)
      IMPLICIT NONE
      INTEGER  I,J
      INTEGER  N,NT,IOUT
      INTEGER  TETRA(5,NT)
      REAL*8   X(N),Y(N),Z(N),POT(N)
      
      
C
C  write coordinates
C
      WRITE(IOUT,*) N,NT
      DO I=1,N
         WRITE(IOUT,200) I,X(I),Y(I),Z(I)
      END DO
C
C  write tetrahedra info
C
      DO I=1,NT
         WRITE(IOUT,300) I,TETRA(5,I),TETRA(1,I),TETRA(2,I),TETRA(3,I),
     1                TETRA(4,I)
      END DO
C
C  write one nodal scalar component 
C
      WRITE(IOUT,100) 1,0
      WRITE(IOUT,100) 1,1
C
C  write the units of measure
C
      WRITE(IOUT,'(A,A)') '  psi',' m '
C
C write out values of component 'potential'
C
      DO I=1,N
         WRITE(IOUT,400) I,POT(I)
      END DO
C
 100  FORMAT(5I10)
 200  FORMAT(I10,3(1PE15.6))
 300  FORMAT(I10,I3,'  tet ',4I10)
 400  FORMAT(I10,1PE15.6)
C
      RETURN
      END
