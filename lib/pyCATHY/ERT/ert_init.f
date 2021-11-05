C
      SUBROUTINE ERT_INIT(ARCHIE,NERT)

C
      IMPLICIT NONE

      INTEGER NERT

      REAL*8 ARCHIE(4)

      OPEN(18,FILE='input/archie',status='old')
      
      READ(18,*) ARCHIE(1)
      READ(18,*) ARCHIE(2)
      READ(18,*) ARCHIE(3)
      READ(18,*) ARCHIE(4)
      READ(18,*) NERT
      CLOSE(18)
      return
      end
c       INCLUDE 'CATHY.H'
c      INTEGER AUGMX,AUGMY,AUGMZ
c      INTEGER NCOL,NROW,NSTR,NCOL3D,NSTR3D,I,J,NODO
c      INTEGER NROW3D,K,EL
c      INTEGER TRIANG2D(4,2*(NCOL+2*AUGMX)*(NROW+2*AUGMY))
cC attenzione! 
c      INTEGER NNOD2D,NTRI2D,N12D,NZONE2D,NERT
c      INTEGER N3D,NT3D,N13D,NZONE3D,UNO
c      PARAMETER (NNOD2D=(NCOL+2*AUGMX+1)*(NROW+2*AUGMY+1)) 
c      PARAMETER (N3D=NNOD2D*(NSTR+AUGMZ+1))
c      PARAMETER (NT3D= 3*(NSTR+AUGMZ)*(NCOL+2*AUGMX)*(NROW+2*AUGMY)) 
c      INTEGER TETRA(5,NT3D) 
c
c
c      REAL*8 ARCHIE(4),ZRATIO(NSTR),BASE
c      REAL*8 ZRATIO3D(NSTR+AUGMZ)
c      REAL*8 DELTA_X,LENGTH_X,LENGTH_Z,A,ALT
c      REAL*8 DELTA_Y
c      REAL*8 X3D(N3D),Y3D(N3D),Z3D(N3D)
c      REAL*8 DEPTH3D(NNOD2D)
c      REAL*8 XC3(NT3D),YC3(NT3D),ZC3(NT3D)
c      REAL*8 RMAX
c      DATA    RMAX/1.7D+100/
       
C
c the grid for sat3d_ert has 8 rows and 9 columns more than
C the surface grid of the DEM in CatHy. (valid for rectangular grids)
C 
c      NROW3D=NROW+2*AUGMY+1
c      NCOL3D=NCOL+2*AUGMX+1
cc      NSTR3D=NSTR+AUGMZ
cC      
cc at first create the 2d mesh using gridgen2d, then use gen3d for the
cc creation of the 3 dim grid.     
c      OPEN(18,FILE='input_ert/grid',status='unknown')
c      DO I=1,NCOL3D*NROW3D
c         X2D(I)=0.0d0
c         Y2D(I)=0.0d0
c      END DO
c
c      CALL ERT_MESH2D(LENGTH_X,LENGTH_Y,NCOL3D-1,NROW3D-1,TRIANG2D)
c      NROW3D=NROW+8+1
c      NCOL3D=NCOL+8+1
c      LENGTH_Y= 1.0d0
c      LENGTH_X= 1.0d0
c      NNOD2D=NROW3D*NCOL3D
c      NTRI2D=2*(NROW3D-1)*(NCOL3D-1)
cc      CALL ERT_COOR2D(LENGTH_X,LENGTH_Y,NCOL3D-1,NROW3D-1,
cc     1          X2D,Y2D,NNOD2D,NT2D)
cc 
cc the grid is red as there is only one zone. Later in the program for
cC ert measurements assimilation to each element will be assigned 
cc a different zone.
cc
cC 
c      LENGTH_X= DELTA_X*NCOL+2*(2**AUGMX-1)*DELTA_X
c      LENGTH_Y= DELTA_Y*NROW+2*(2**AUGMY-1)*DELTA_Y
c      LENGTH_Z= BASE+ZRATIO(NSTR)*BASE*AUGMZ
c      DO I=1,NSTR3D
c         IF (I.LE.NSTR) THEN
c            ZRATIO3D(I)=ZRATIO(I)*BASE/LENGTH_Z
c         ELSE
c             ZRATIO3D(I)=ZRATIO3D(NSTR)
c         END IF
c      END DO
c      WRITE(*,*) LENGTH_X,LENGTH_Y,LENGTH_Z,'LENGTH_X,LENGTH_Y,LENGTH_Z'
c      WRITE(*,*) (ZRATIO3D(I),I=1,NSTR3D),'ZRATIO3D'
c      LX=0.0d0
c      DO I=1,NCOL3D
c         LY= LENGTH_Y
c         DO J=1,NROW3D
c            NODO= (I-1)*NROW3D+J
c            X3D(NODO)=LX
c            Y3D(NODO)=LY
c            IF (J.LE. AUGMY) THEN
c               LY=LY-2**(AUGMY-J)*DELTA_Y
c            ELSE IF((J.GT.AUGMY).AND.(J.LE.NROW3D-AUGMY) THEN
c               LY=LY-DELTA_Y
c            ELSE
c               LY=LY-2**(J-NROW3D+AUGMY-1)*DELTA_Y
c            END IF
c         END DO
c         IF (I.LE. AUGMX) THEN
c            LX=LX+2**(AUGMX-I)*DELTA_X
c         ELSE IF((I.GT.AUGMX).AND.(I.LE.NCOL3D-AUGMX) THEN
c            LX=LX+DELTA_X
c         ELSE
c            LX=LX-2**(I-NCOL3D+AUGMX-1)*DELTA_X
c         END IF
c      END DO
c      WRITE(*,*) LX,LY,'LX,LZ'
c      DO I=1,NNOD2D
c          DEPTH3D(I)=0.0d0
c      END DO
c      CALL RIORD(NTRI2D,TRIANG2D,3)
c      CALL GEN3D(N3D,NT3D,NNOD2D,NSTR3D,NTRI2D,0,0,
c     1           TRIANG2D,TETRA3D,RMAX,LENGTH_Z,DEPTH3D,ZRATIO3D,
c     2      X3D,Y3D,Z3D,XC3D,YC3D,ZC3D)
c      NZONE3D=1
c      N13D=30
c      UNO=1
c      WRITE(18,*)  N3D,NT3D
c      WRITE(18,*)  NSTR3D,NZONE3D,N13D
c      Do I =1,NSTR3D
c         DO K=1,3*(NCOL2D-1)*(NROW2D-1)
c            EL=(I-1)*3 *(NCOL2D-1)*(NROW2D-1)+K
c            WRITE(18,*) EL,(TETRA(J,EL),J=1,4),I,UNO
c         END DO
c      END DO
c      Do I=1,N3D
c         WRITE(18,*) I, X3D(I),Y3D(I),Z3D(I)
c      END Do
c      CLOSE(18)

c      return
c      end
