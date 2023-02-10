C
C**************************  HGRAPH ************************************
C
C  calculate atmospheric component of hydrograph and
C  flag possible inconsistent or anomalous situations
C
C  NOTE: this routine has to be completely checked & updated!!!
C
C***********************************************************************
C
      SUBROUTINE HGRAPH(NNOD,TIME,APOT,AACT,
     1                  OVFLOW,REFL,HGFLAG,IFATM,ATMPOT,ATMACT,PNEW,
     2                  evap_eff,inf_tot)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   K
      INTEGER   NNOD
      INTEGER   HGFLAG(9),IFATM(*)
      real*8    evap,infflow,reflow
      REAL*8    TIME,APOT,AACT,OVFLOW,REFL
      REAL*8    evap_eff,inf_tot
      REAL*8    ZERO
      REAL*8    ATMPOT(*),ATMACT(*),PNEW(*)
      INCLUDE  'IOUNITS.H'
      INCLUDE  'SOILCHAR.H'
      PARAMETER (ZERO=0.0D0)
C
      APOT=ZERO
      AACT=ZERO
      refl=ZERO
      infflow=zero
      evap = zero
      OVFLOW=ZERO !ms
C
      DO K=1,NNOD
         APOT = APOT + ATMPOT(K)
         AACT = AACT + ATMACT(K)
C
C  case IFATM = 2 (Dirichlet ponding condition)
C
         IF (IFATM(K) .EQ. 2) THEN
            IF (ATMACT(K) .LT. ZERO) THEN              
C           ................................................ return flow
C                                                       (PNEW>POND_HMIN)
               IF (ATMPOT(K) .GE. ZERO) THEN           
C              ........................................... precipitation
               ELSE                                    
C              ............................................. evaporation
                  IF (ATMACT(K) .LE. ATMPOT(K)) THEN    
C                 ......... enough water from soil alone to satisfy evap
                  evap=evap+atmpot(k)   
                  ELSE 
C                 ........... not enough water from soil to satisfy evap
C                    ...... need to get contribution from surface volume
                  evap=evap+atmpot(k)
                  END IF
               END IF
               refl = refl - ATMACT(K)
               OVFLOW = OVFLOW - ATMACT(K) + ATMPOT(K)
            ELSE
C           ............................................... infiltration
               IF (ATMPOT(K) .GE. ZERO) THEN
C              ........................................... precipitation
                  IF (ATMACT(K) .LE. ATMPOT(K)) THEN
C                 .... enough water from precip alone to satisfy infiltr
                  ELSE 
C                 ...... not enough water from precip to satisfy infiltr
C                    ...... need to get contribution from surface volume
                  END IF
               ELSE
C              ............................................. evaporation
C                    ...... need to get contribution from surface volume
                  evap=evap+atmpot(k)
               END IF
               infflow = infflow + atmact(k)
               OVFLOW = OVFLOW - ATMACT(K) + ATMPOT(K)
            END IF
C
C  case IFATM = 1 (Dirichlet non-ponding condition)
C
         ELSE IF (IFATM(K) .EQ. 1) THEN
            IF (PNEW(K) .GE. ZERO) THEN
C           ------------------------------------------ surface saturated
               IF (ATMACT(K) .LT. ZERO) THEN
C              ............................................ exfiltration
                  IF (ATMPOT(K) .GE. ZERO) THEN
C                 ........................................ precipitation
C                    ...................... only return flow is possible
                     refl = refl - ATMACT(K)
                     OVFLOW = OVFLOW - ATMACT(K) + ATMPOT(K)
                  ELSE
C                 .......................................... evaporation
                     HGFLAG(5)=HGFLAG(5) + 1
                     WRITE(IOUT9,1500) K,TIME,ATMACT(K),ATMPOT(K)
C
                     IF (ATMACT(K) .LE. ATMPOT(K)) THEN
C                    ............................. ATMPOT is evaporation
                        evap=evap+atmpot(k)
                        refl = refl - ATMACT(K) + ATMPOT(K)
                        OVFLOW = OVFLOW - ATMACT(K) + ATMPOT(K)
                     ELSE
C                    ................................ all is evaporation
                        evap=evap+atmact(k)
                     END IF
                  END IF
               ELSE
C              ............................................ infiltration
                  IF (ATMPOT(K) .GE. ZERO) THEN
C                 ........................................ precipitation
                     IF (ATMACT(K) .LE. ATMPOT(K)) THEN
C                    ..................................... overland flow
                        OVFLOW = OVFLOW - ATMACT(K) + ATMPOT(K)
                     ELSE
C                    .................... not enough water to infiltrate
C                       .................................. INCONSISTENCY
                        HGFLAG(1) = HGFLAG(1) + 1
                        WRITE(IOUT9,1100) K,TIME,ATMACT(K),
     1                                    ATMPOT(K)
                     END IF
                  ELSE
C                 .......................................... evaporation
C                    .. evaporation + infiltration in a non ponding case
C                    ..................................... INCONSISTENCY
                     HGFLAG(6)=HGFLAG(6) + 1
                     WRITE(IOUT9,1600) K,TIME,ATMACT(K),ATMPOT(K)
                     evap=evap+atmpot(k)
                  END IF
                  infflow = infflow + atmact(k)
               END IF
            ELSE IF (PNEW(K) .LE. PMIN) THEN
C           ---------------------------------------- surface unsaturated
C                 ...only dry soil (PNEW<PMIN) can exist in this case...
               IF (ATMACT(K) .LT. ZERO) THEN
C              ............................................ exfiltration
                  IF (ATMPOT(K) .GE. ZERO) THEN
C                 ........................................ precipitation
C                    ..................................... INCONSISTENCY
C                    .................... soil can take all precip water
C                    ......................... but there is exfiltration
                     OVFLOW = OVFLOW + ATMPOT(K)
                     HGFLAG(2)=HGFLAG(2) + 1
                     WRITE(IOUT9,1200) K,TIME,ATMACT(K),ATMPOT(K)
                  ELSE
C                 .......................................... evaporation
C                    ............... keep evaporatin as much as possible
C                    ............... with fixed head = PMIN (see switch)
                     evap=evap+atmact(k)
                  END IF
               ELSE
C              ............................................ infiltration
C                    ..................................... INCONSISTENCY
                  HGFLAG(4)=HGFLAG(4) + 1
                  WRITE(IOUT9,1400) K,TIME,ATMACT(K),ATMPOT(K)
                  IF (ATMPOT(K) .GE. ZERO) THEN
C                 ........................................ precipitation
                     IF (ATMACT(K) .LE. ATMPOT(K)) THEN
C                    ........................ INCONSISTENT overland flow
                        OVFLOW = OVFLOW - ATMACT(K) + ATMPOT(K)
                        HGFLAG(9)=HGFLAG(9) + 1
                        WRITE(IOUT9,1900) K,TIME,ATMACT(K),
     1                                    ATMPOT(K)
                     ELSE
C                    ....... INCONSISTENT not enough water to infiltrate
                        HGFLAG(3)=HGFLAG(3) + 1
                        WRITE(IOUT9,1300) K,TIME,ATMACT(K),
     1                                    ATMPOT(K)
                     END IF
                  ELSE
C                 .......................................... evaporation
C                    .......... evaporation + infiltration on a dry soil
C                    ..................................... INCONSISTENCY
                        HGFLAG(8)=HGFLAG(8) + 1
                        WRITE(IOUT9,1800) K,TIME,ATMACT(K),
     1                                    ATMPOT(K)
                  END IF
                  infflow = infflow + atmact(k)
               END IF
            END IF
C
C  case IFATM = 0 (Neumann condition)
C
         ELSE IF (IFATM(K) .EQ. 0) THEN
C
C            in this case we should have always ATMACT=ATMPOT, 
C            and thus no overland flow, except in the case PNEW<PMIN
C            in which case we only have evaporation and again 
C            no overland flow
C
             if(atmpot(k).lt.zero) then
               evap=evap+atmact(k)
             else
               infflow = infflow + atmact(k)
             end if
         END IF
      END DO
cp      evap_eff = evap_eff + evap*deltat
cp      inf_tot = inf_tot + infflow*deltat
cp      reflow = reflow + refl*deltat
C
C
      RETURN
 1100 FORMAT(  ' HGFLAG(1) AT SATURATED SURFACE NODE ',I6,
     1         ' (TIME=',1PE8.2,')',
     2       /,11X,'ACTUAL    INFLOW: INFILTRATION  (',1PE12.5,')',
     3       /,11X,'POTENTIAL INFLOW: PRECIPITATION (',1PE12.5,')',
     4       /,11X,'MAGNITUDE OF ACTUAL EXCEEDS MAGNITUDE OF',
     5         ' POTENTIAL')
 1200 FORMAT(  ' HGFLAG(2) AT AIR DRY   SURFACE NODE ',I6,
     1         ' (TIME=',1PE8.2,' )',
     2       /,11X,'ACTUAL    INFLOW: EXFILTRATION  (',1PE12.5,')',
     3       /,11X,'POTENTIAL INFLOW: PRECIPITATION (',1PE12.5,')',
     4       /,11X,'SURFACE RUNOFF GENERATED AT AN AIR DRY NODE')
 1300 FORMAT(  ' HGFLAG(3) AT AIR DRY   SURFACE NODE ',I6,
     1         ' (TIME=',1PE8.2,' )',
     2       /,11X,'ACTUAL    INFLOW: INFILTRATION  (',1PE12.5,')',
     3       /,11X,'POTENTIAL INFLOW: PRECIPITATION (',1PE12.5,')',
     4       /,11X,'MAGNITUDE OF ACTUAL EXCEEDS MAGNITUDE OF',
     5         ' POTENTIAL')
 1400 FORMAT(  ' HGFLAG(4) AT AIR DRY   SURFACE NODE ',I6,
     1         ' (TIME=',1PE8.2,' )',
     2       /,11X,'ACTUAL    INFLOW: INFILTRATION  (',1PE12.5,')',
     3       /,11X,'POTENTIAL INFLOW: PRECIPITATION (',1PE12.5,')',
     4       /,11X,'INFILTRATION AT AN AIR DRY NODE')
 1500 FORMAT(  ' HGFLAG(5) AT SATURATED SURFACE NODE ',I6,
     1         ' (TIME=',1PE8.2,')',
     2       /,11X,'ACTUAL    INFLOW: EXFILTRATION  (',1PE12.5,')',
     3       /,11X,'POTENTIAL INFLOW: EVAPORATION   (',1PE12.5,')',
     4       /,11X,'EXFILTRATION AT A SATURATED NODE')
 1600 FORMAT(  ' HGFLAG(6) AT SATURATED SURFACE NODE ',I6,
     1         ' (TIME=',1PE8.2,' )',
     2       /,11X,'ACTUAL    INFLOW: INFILTRATION  (',1PE12.5,')',
     3       /,11X,'POTENTIAL INFLOW: EVAPORATION   (',1PE12.5,')',
     4       /,11X,'INFILTRATION AT A SATURATED NODE WITH',
     5         ' EVAPORATIVE POTENTIAL FLUX')
 1700 FORMAT(  ' HGFLAG(7) AT AIR DRY   SURFACE NODE ',I6,
     1         ' (TIME=',1PE8.2,' )',
     2       /,11X,'ACTUAL    INFLOW: EXFILTRATION  (',1PE12.5,')',
     3       /,11X,'POTENTIAL INFLOW: EVAPORATION   (',1PE12.5,')',
     4       /,11X,'MAGNITUDE OF ACTUAL EXCEEDS MAGNITUDE OF',
     5         ' POTENTIAL')
 1800 FORMAT(  ' HGFLAG(8) AT AIR DRY   SURFACE NODE ',I6,
     1         ' (TIME=',1PE8.2,')',
     2       /,11X,'ACTUAL    INFLOW: INFILTRATION  (',1PE12.5,')',
     3       /,11X,'POTENTIAL INFLOW: EVAPORATION   (',1PE12.5,')',
     4       /,11X,'INFILTRATION AT AN AIR DRY NODE WITH',
     5         ' EVAPORATIVE POTENTIAL FLUX')
 1900 FORMAT(  ' HGFLAG(9) AT AIR DRY   SURFACE NODE ',I6,
     1         ' (TIME=',1PE8.2,')',
     2       /,11X,'ACTUAL    INFLOW: INFILTRATION  (',1PE12.5,')',
     3       /,11X,'POTENTIAL INFLOW: PRECIPITATION (',1PE12.5,')',
     4       /,11X,'SURFACE RUNOFF GENERATED AT AN AIR DRY NODE')
      END
