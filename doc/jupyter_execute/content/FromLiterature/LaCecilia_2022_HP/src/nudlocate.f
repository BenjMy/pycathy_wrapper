C
C**************************  NUDLOCATE *********************************
C
C  locate the nudging observation points within our grid.
C
C***********************************************************************
C
      SUBROUTINE NUDLOCATE(NT,NUDN,TETRA,NUDTET,
     1                     NUDX,NUDY,NUDZ,X,Y,Z,VOLU)
C
      IMPLICIT NONE
      INTEGER  I,L,TK
      INTEGER  NT,NUDN
      INTEGER  IV(0:3)
      INTEGER  TETRA(5,*),NUDTET(*)
      REAL*8   PHI01_1,PHI01_2,PHI01_3,PHI01_4
      REAL*8   DEN,XX,YY,ZZ,F1,F2,F3,F4
      REAL*8   XXP(0:3),YYP(0:3),ZZP(0:3),XXS(3),YYS(3),ZZS(3)
      REAL*8   A(3,3),XL(3)
      REAL*8   NUDX(*),NUDY(*),NUDZ(*),X(*),Y(*),Z(*),VOLU(*)
      INCLUDE 'IOUNITS.H'
C
C  define the P1 step-wise function on the tetrahedra
C  (0,0,0) (1,0,0) (0,1,0) (0,0,1)
C
      PHI01_1(XX,YY,ZZ) = 1.0D0 - XX - YY - ZZ
      PHI01_2(XX,YY,ZZ) = XX
      PHI01_3(XX,YY,ZZ) = YY
      PHI01_4(XX,YY,ZZ) = ZZ
C
      DO L=1,NUDN
         NUDTET(L)=0
         DO TK=1,NT
            DO I=0,3
               IV(I)=TETRA(I+1,TK)
            END DO
            DO I=0,3
               XXP(I)=X(IV(I))
               YYP(I)=Y(IV(I))
               ZZP(I)=Z(IV(I))
            END DO
C
C  compute coordinates of tetrahedra TK relative to node (0,0,0)
C  and the determinant DEN
C
            DO I=1,3
               XXS(I)=XXP(I) - XXP(0)
               YYS(I)=YYP(I) - YYP(0)
               ZZS(I)=ZZP(I) - ZZP(0)
            END DO
            DEN=6.0D0 * ABS(VOLU(TK))
C
C  matrix of coordinate transformation from local to isoparametric
C
            A(1,1) =  (YYS(2)*ZZS(3) - ZZS(2)*YYS(3)) / DEN
            A(1,2) = -(XXS(2)*ZZS(3) - ZZS(2)*XXS(3)) / DEN
            A(1,3) =  (XXS(2)*YYS(3) - YYS(2)*XXS(3)) / DEN
            A(2,1) = -(YYS(1)*ZZS(3) - ZZS(1)*YYS(3)) / DEN
            A(2,2) =  (XXS(1)*ZZS(3) - ZZS(1)*XXS(3)) / DEN
            A(2,3) = -(XXS(1)*YYS(3) - YYS(1)*XXS(3)) / DEN
            A(3,1) =  (YYS(1)*ZZS(2) - ZZS(1)*YYS(2)) / DEN
            A(3,2) = -(XXS(1)*ZZS(2) - ZZS(1)*XXS(2)) / DEN
            A(3,3) =  (XXS(1)*YYS(2) - YYS(1)*XXS(2)) / DEN
C
C  determine if the observation point L lies within the tetrehedra TK.
C  As soon as we locate the tetrahedra we proceed to the next
C  observation point. If no tetrahedra is found, terminate with error.
C
            XL(1)=NUDX(L) - XXP(0)
            XL(2)=NUDY(L) - YYP(0)
            XL(3)=NUDZ(L) - ZZP(0)
            XX=0.0D0
            YY=0.0D0
            ZZ=0.0D0
            DO I=1,3
               XX=XX + A(1,I)*XL(I)
               YY=YY + A(2,I)*XL(I)
               ZZ=ZZ + A(3,I)*XL(I)
            END DO
            F1=PHI01_1(XX,YY,ZZ)
            F2=PHI01_2(XX,YY,ZZ)
            F3=PHI01_3(XX,YY,ZZ)
            F4=PHI01_4(XX,YY,ZZ)
            IF (F1 .GE. 0.0D0 .AND. F2 .GE. 0.0D0 .AND.
     1          F3 .GE. 0.0D0 .AND. F4 .GE. 0.0D0) THEN
               NUDTET(L)=TK
               GO TO 100
            END IF
         END DO
         WRITE(IOUT2,1000) L,NUDX(L),NUDY(L),NUDZ(L)
         CALL CLOSIO
         STOP
  100    CONTINUE
      END DO
C
      RETURN
 1000 FORMAT(/,' ERROR IN SUBROUTINE NUDLOCATE: OBSERVATION POINT ',I5,
     1       /,' WITH COORDINATES  (',3(1PE12.4),' )  COULD NOT BE',
     2       /,' LOCATED WITHIN THE DOMAIN')
      END
