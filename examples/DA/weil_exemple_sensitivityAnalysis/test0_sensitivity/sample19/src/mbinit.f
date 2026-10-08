C
C**************************  MBINIT ************************************
C
C  initialize mass balance and hydrograph parameters for first
C  time level (time 0.0)
C
C***********************************************************************
C
      SUBROUTINE MBINIT(SFFLWP,AACTP,ADINP,ADOUTP,NDINP,NDOUTP,
     1                  ANINP,ANOUTP,NNINP,NNOUTP,
     2                  NUDINP,NUDOUTP,VNUDTOT,
     3                  NNOD,IFATMP,ATMOLD,NQ,Q,
     4                  VAPOT_T,VAACT_T,VSFTOT,VNDTOT,VNNTOT,VTOT,
     5                  STORE2,CDSTOR,CVIN,CVOUT,CERRAS,CAERAS)
C
      IMPLICIT NONE
      INTEGER  K
      INTEGER  NNOD,NQ
      INTEGER  IFATMP(*)
      REAL*8   SFFLWP,AACTP,ADINP,ADOUTP,NDINP,NDOUTP
      REAL*8   ANINP,ANOUTP,NNINP,NNOUTP
      REAL*8   NUDINP,NUDOUTP,VNUDTOT
      REAL*8   VAPOT_T,VAACT_T,VSFTOT,VNDTOT,VNNTOT,VTOT
      REAL*8   STORE2,CDSTOR,CVIN,CVOUT,CERRAS,CAERAS
      REAL*8   ATMOLD(*),Q(*)
      INCLUDE 'IOUNITS.H'
C
      WRITE(IOUT7,1000)
      WRITE(IOUT8,1010)
      WRITE(IOUT30,1020)
      WRITE(IOUT31,1030)
      WRITE(IOUT32,1040)
      WRITE(IOUT50,1050)
      WRITE(IOUT36,1060)
      SFFLWP=0.0D0
      AACTP=0.0D0
      ADINP=0.0D0
      ADOUTP=0.0D0
      NDINP=0.0D0
      NDOUTP=0.0D0
      ANINP=0.0D0
      ANOUTP=0.0D0
      DO K=1,NNOD
         IF (IFATMP(K) .EQ. 0) THEN
            AACTP=AACTP + ATMOLD(K)
            IF (ATMOLD(K) .GT. 0.0D0) THEN
               ANINP=ANINP + ATMOLD(K)
            ELSE
               ANOUTP=ANOUTP + ATMOLD(K)
            END IF
         END IF
      END DO
      NNINP=0.0D0
      NNOUTP=0.0D0
      DO K=1,NQ
         IF (Q(K) .GT. 0.0D0) THEN
            NNINP=NNINP + Q(K)
         ELSE
            NNOUTP=NNOUTP + Q(K)
         END IF
      END DO
      NUDINP=0.0D0
      NUDOUTP=0.0D0
      VAPOT_T=0.0D0
      VAACT_T=0.0D0
      VSFTOT=0.0D0
      VNDTOT=0.0D0
      VNNTOT=0.0D0
      VNUDTOT=0.0D0
      VTOT=0.0D0
      STORE2=0.0D0
      CDSTOR=0.0D0
      CVIN=0.0D0
      CVOUT=0.0D0
      CERRAS=0.0D0
      CAERAS=0.0D0
C
      RETURN
 1000 FORMAT('#NSTP            DELTAT         TIME    POT. FLUX',
     1       '    ACT. FLUX    OVL. FLUX    RET. FLUX    SEEP FLUX',
     2       '    REC. FLUX     REC.VOL.')
 1010 FORMAT('#NSTEP     DELTAT       TIME   NET NATM,NSF DIR FLUX',
     1       '   NET NATM,NSF NEU FLUX')
 1020 FORMAT('# NSTEP            DELTAT              TIME',
     1       '  NET SEEPFACE VOL  NET SEEPFACE FLX')
 1030 FORMAT('# NSTEP            DELTAT              TIME',
     1       ' NET NANSF DIR VOL NET NANSF DIR FLX')
 1040 FORMAT('# NSTEP            DELTAT              TIME',
     1       ' NET NANSF NEU VOL NET NANSF NEU FLX')
 1050 FORMAT(18X,'#***** Detailed nudging "hydrograph" ***** ',
     1         /,'#nstep    deltat      time   in flux',
     2           '  out flux  in+out f  in volum  out volu',
     3           '  in+out v  invol/dt   outv/dt  iv+ov/dt',
     4           '   cum vol')
 1060 FORMAT(21X,'#***** Cumulative flow volumes ***** ',
     1         /,'# Nstep    Deltat      Time  SeepageF',
     2           '  nansfDir  nansfNeu   Nudging  net (VTOT)')
      END
