C
C**************************  GRDSYS ************************************
C
C  further processing of the 2-D mesh (sorting and area calculation),
C  generation of the 3-D grid, and setup of system matrices
C
C***********************************************************************
C
      SUBROUTINE GRDSYS(IPRT1,IOPT,N,NT,NNOD,NTRI,NTERM,N1,NDZ,
     1                  NSTR,IVERT,LUMP,IMAX,
     2                  TRIANG,TETRA,IA,JA,TOPOL,TETJA,IVOL,
     3                  IP3,IP4,IER,
     4                  RMAX,BASE,DEPTH,ZRATIO,ARENOD,X,Y,Z,XC,YC,ZC,
     5                  AI,BI,CI,DI,VOLNOD,VOLU,VOLUR,LMASS)
C
      IMPLICIT  NONE
      INTEGER   IPRT1,IOPT,N,NT,NNOD,NTRI,NTERM,N1,NDZ,I,J
      INTEGER   NSTR,IVERT,LUMP,IMAX
      INTEGER   TRIANG(4,*),TETRA(5,*),IA(*),JA(*),TOPOL(*),TETJA(4,4,*)
      INTEGER   IVOL(*),IP3(3,3),IP4(4,4),IER(*)
      REAL*8    RMAX,BASE,DEPTH(*)
      REAL*8    ZRATIO(*),ARENOD(*),X(*),Y(*),Z(*),XC(*),YC(*),ZC(*)
      REAL*8    AI(4,*),BI(4,*),CI(4,*),DI(4,*)
      REAL*8    VOLNOD(*),VOLU(*),VOLUR(*),LMASS(4,4)
      INCLUDE   'IOUNITS.H'
C
C  sort into ascending order the nodes connecting to each triangular
C  element
C
      CALL RIORD(NTRI,TRIANG,3)
C     WRITE(666,*) NTRI
C     DO I=1,NTRI
C        WRITE(666,*) (TRIANG(J,I),J=1,4)
C     END DO
C
C  calculate the area assigned to each surface node
C
      
      CALL AREA2D(NNOD,NTRI,TRIANG,IP3,ARENOD,X,Y)
C
C  generate the three-dimensional grid
C
      
      CALL GEN3D(N,NT,NNOD,NSTR,NTRI,IPRT1,IVERT,
     1           TRIANG,TETRA,
     2           RMAX,BASE,DEPTH,ZRATIO,X,Y,Z,XC,YC,ZC)
     
C
C  legge le coordinate della mesh dall'esterno se necessario   
C
      IF (IPRT1 .EQ. -1) THEN
         DO I=1,N
            READ(IIN51,*) X(I),Y(I),Z(I)
         END DO
         CALL NODELT(NT,TETRA,X,XC)
         CALL NODELT(NT,TETRA,Y,YC)
         CALL NODELT(NT,TETRA,Z,ZC)
      END IF
C
C  sort into ascending order the nodes connecting to each tetrahedral
C  element. Sorting is not necessary for the nonsymmetric case.
C
      IF (IOPT .EQ. 1) CALL RIORD(NT,TETRA,4)
C
C  calculate volumes and basis function coefficients
C
      CALL VOLBAS(N,NT,TETRA,IP4,AI,BI,CI,DI,
     1            VOLNOD,VOLU,VOLUR,IVOL,X,Y,Z)
C
C  set up LMASS, the part of the local mass matrix which is constant
C  for all elements
C
      CALL LOCMAS(LMASS,LUMP)
C
C  set up pointers and indices for storage of the system matrices
C
C      IF (IOPT .EQ. 1) THEN
C         CALL STRPIC(N,NTERM,TETRA,JA,TOPOL,NT,N1,IMAX)
C         CALL CHKPIC(N,NTERM,TOPOL,JA,NDZ,IER)
C         CALL TETPIC(NT,TETRA,JA,TOPOL,TETJA)
C      ELSE
C         CALL STRNEW(N,NTERM,TETRA,JA,TOPOL,NT,N1,IMAX)
C         CALL CHKNEW(N,NTERM,TOPOL,JA,NDZ,IER)
C         CALL TETNEW(NT,TETRA,JA,TOPOL,TETJA)
C         CALL TOPIA(N,TOPOL,IA)
C      END IF
C
      RETURN
      END
