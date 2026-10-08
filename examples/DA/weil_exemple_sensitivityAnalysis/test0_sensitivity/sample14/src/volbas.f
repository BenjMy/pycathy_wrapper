C
C**************************  VOLBAS ************************************
C
C  calculate the volume of each element, the volume assigned to each
C  node, and the basis function coefficients. 
C  Basis function coefficients are divided by 6.
C
C***********************************************************************
C
      SUBROUTINE VOLBAS(N,NT,TETRA,IP4,AI,BI,CI,DI,
     1                  VOLNOD,VOLU,VOLUR,IVOL,X,Y,Z)
C
      IMPLICIT  NONE
      INTEGER   J,IEL,IV1,IV2,IV3,IV4
      INTEGER   N,NT
      INTEGER   TETRA(5,*),IP4(4,4),IVOL(*)
      REAL*8    VOL,VOL4,R4
      REAL*8    AI(4,*),BI(4,*),CI(4,*),DI(4,*)
      REAL*8    VOLNOD(*),VOLU(*),VOLUR(*),X(*),Y(*),Z(*)
      INCLUDE  'IOUNITS.H'
C
      R4=0.25D0
      CALL INIT0R(N,VOLNOD)
      DO IEL=1,NT
         CALL VOLUMS(IP4,TETRA(1,IEL),X,Y,Z,VOL)
         CALL BASIS6(IP4,TETRA(1,IEL),X,Y,Z,AI(1,IEL),
     1               BI(1,IEL),CI(1,IEL),DI(1,IEL))
         IF (VOL .EQ. 0.0D0) THEN
            WRITE(IOUT2,1000) IEL,(TETRA(J,IEL),J=1,4)
            CALL CLOSIO
            STOP
         END IF
         VOL4=DABS(VOL)*R4
         IV1=TETRA(1,IEL)
         IV2=TETRA(2,IEL)
         IV3=TETRA(3,IEL)
         IV4=TETRA(4,IEL)
         VOLNOD(IV1)=VOLNOD(IV1) + VOL4
         VOLNOD(IV2)=VOLNOD(IV2) + VOL4
         VOLNOD(IV3)=VOLNOD(IV3) + VOL4
         VOLNOD(IV4)=VOLNOD(IV4) + VOL4
         IVOL(IEL)=1
         IF (VOL .LT. 0.0D0) THEN
            IVOL(IEL)=-1
            VOL=-VOL
         END IF
         VOLU(IEL)=VOL
         VOLUR(IEL)=1.0D0/VOL
      END DO
C
      RETURN
 1000 FORMAT(/,' ERROR IN SUBROUTINE VOLBAS: ZERO AREA CALCULATED',/,
     1         '    AT ELEMENT ',I6,'  NODE NUMBERS:',4I6)
      END
