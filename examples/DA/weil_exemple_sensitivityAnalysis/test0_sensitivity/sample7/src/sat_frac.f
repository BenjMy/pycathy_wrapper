C
C**************************  SAT_FRAC **********************************
C
C  calculate the fraction of the surface that is Horton or Dunne
C  saturated or ponded (note that these fractions are calculated using
C  PNEW and not IFATM)
C
C***********************************************************************
C
      SUBROUTINE SAT_FRAC(NNOD,NSTR,PONDH_MIN,FHORT,FDUNN,FPOND,FSAT,
     1                    PNEW)
C
      IMPLICIT  NONE
      INTEGER   I,K,KK,HDFLAG,NHORT,NDUNN,NPOND,NSAT
      INTEGER   NNOD,NSTR
      REAL*8    PONDH_MIN,FHORT,FDUNN,FPOND,FSAT
      REAL*8    PNEW(*)
C
      NHORT=0
      NDUNN=0
      NPOND=0
      NSAT=0
      DO I=1,NNOD
         IF (PNEW(I) .GE. 0.0D0) THEN
            NSAT=NSAT+1
            IF (PNEW(I) .GE. PONDH_MIN) NPOND=NPOND+1
	    HDFLAG=0
            DO K=1,NSTR
               KK=K*NNOD + I
               IF (PNEW(KK) .LT. 0.0D0) HDFLAG=1
            END DO
	    IF (HDFLAG .EQ. 0) THEN
	       NDUNN=NDUNN+1
	    ELSE
	       NHORT=NHORT+1
            END IF
         END IF
      END DO
      FHORT=DFLOAT(NHORT)/DFLOAT(NNOD)
      FDUNN=DFLOAT(NDUNN)/DFLOAT(NNOD)
      FPOND=DFLOAT(NPOND)/DFLOAT(NNOD)
      FSAT=DFLOAT(NSAT)/DFLOAT(NNOD)
C
      RETURN
      END
