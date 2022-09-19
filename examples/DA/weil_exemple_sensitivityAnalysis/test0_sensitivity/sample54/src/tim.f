      SUBROUTINE TIM(TIME,ICOD)
      INTEGER ICOD
      REAL*4 TIME,TMP
      real*4 tarray(2)
CM    real*4 etime
c
      IF (ICOD.EQ.1) THEN
c
C  SETTING INIZIALE
c
CM       TMP = etime(tarray)
         CALL ETIME(tarray,TMP)
         TIME = tarray(1)
      ELSE
c
C  RILEVAZIONE PERIODICA
c
CM       TMP = etime(tarray)
         CALL ETIME(tarray,TMP)
         TIME = tarray(1) - TIME
      ENDIF
c
      RETURN
      END
