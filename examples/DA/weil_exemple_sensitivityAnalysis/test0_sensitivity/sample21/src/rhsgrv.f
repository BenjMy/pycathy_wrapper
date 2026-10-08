C
C**************************  RHSGRV ************************************
C
C  add contribution of the unsaturated zone gravitational term 
C  to the RHS vector
C
C***********************************************************************
C
      SUBROUTINE RHSGRV(NT,NTRI,TETRA,TNOTI,DI,PERMZ,IVOL,CKRWE)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  IN,KE,INOD
      INTEGER  IR,ISTR
      INTEGER  NT,NTRI
      INTEGER  TETRA(5,*),IVOL(*)
      REAL*8   T3,T1
      REAL*8   TNOTI(*),PERMZ(MAXSTR,*),DI(4,*),CKRWE(*)
C
      DO KE=1,NT
         ISTR=1+KE/(NTRI*3)
         IR=MOD(KE,NTRI*3)
         IF (IR .EQ. 0) ISTR=ISTR-1
         T3=CKRWE(KE)*PERMZ(ISTR,TETRA(5,KE))
         DO IN=1,4
            INOD=TETRA(IN,KE)
            T1=T3*DI(IN,KE)*IVOL(KE)
            TNOTI(INOD)=TNOTI(INOD)-T1
         END DO
      END DO
C
      RETURN
      END
