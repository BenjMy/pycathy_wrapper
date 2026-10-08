C
C**************************  AREA2D ************************************
C
C  calculate the area assigned to each surface node
C
C***********************************************************************
C
      SUBROUTINE AREA2D(NNOD,NTRI,TRIANG,IP3,ARENOD,X,Y)
C
      IMPLICIT  NONE
      INTEGER   J,IEL,IA1,IA2,IA3
      INTEGER   NNOD,NTRI
      INTEGER   TRIANG(4,*),IP3(3,3)
      REAL*8    ARE,ARE3,R3
      REAL*8    ARENOD(*),X(*),Y(*)
      INCLUDE  'IOUNITS.H'
C
      R3=1.0D0/3.0D0
      CALL INIT0R(NNOD,ARENOD)
      DO IEL=1,NTRI
         CALL AREAS(IP3,TRIANG(1,IEL),X,Y,ARE)
         IF (ARE .EQ. 0.0D0) THEN
            WRITE(IOUT2,1000) IEL,(TRIANG(J,IEL),J=1,3)
            CALL CLOSIO
            STOP
         END IF
         ARE3=DABS(ARE)*R3
         IA1=TRIANG(1,IEL)
         IA2=TRIANG(2,IEL)
         IA3=TRIANG(3,IEL)
         ARENOD(IA1)=ARENOD(IA1) + ARE3
         ARENOD(IA2)=ARENOD(IA2) + ARE3
         ARENOD(IA3)=ARENOD(IA3) + ARE3
      END DO
c      do j=1,nnod
c       write(88,*) j,arenod(j)
c      end do

C
      RETURN
 1000 FORMAT(/,' ERROR IN SUBROUTINE AREA2D: ZERO AREA CALCULATED',/,
     1         '    AT ELEMENT ',I6,'  NODE NUMBERS:',3I6)
      END
