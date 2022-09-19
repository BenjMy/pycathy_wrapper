C
C**************************  TIMUPDSUP *********************************
C
C  update TIME and time step size for the next time step if FL3D=FALSE
C
C***********************************************************************
C
      SUBROUTINE TIMUPDSUP(TIME,DELTATS,TIMESTOP)
C
      IMPLICIT  NONE
      REAL*8    TIME,DELTATS
      REAL*8    TIMESTOP
C
C Time step correction
C
      IF ((TIME+DELTATS) .GT. TIMESTOP)THEN
         DELTATS=TIMESTOP-TIME
         TIME=TIMESTOP
      ELSE
         TIME=TIME + DELTATS
      END IF
C
      RETURN
      END
