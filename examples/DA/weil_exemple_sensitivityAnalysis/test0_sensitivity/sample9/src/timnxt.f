C
C**************************  TIMNXT ************************************
C
C  update or re-initialize counters for next time level
C
C***********************************************************************
C
      SUBROUTINE TIMNXT(NSTEP,ITER,NITERT,KBACKT,NSURFT,CPUSUB,CPUSURF)
C
      IMPLICIT  NONE
      INTEGER   NSTEP,ITER,NITERT,KBACKT,NSURFT
      REAL      CPUSUB,CPUSURF
C
      NSTEP=NSTEP + 1
      ITER=1
      NITERT=0
      KBACKT=0
      NSURFT=0
      CPUSUB=0.0
      CPUSURF=0.0
C
      RETURN
      END
