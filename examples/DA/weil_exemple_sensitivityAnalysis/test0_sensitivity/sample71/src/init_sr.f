C
C**************************  INIT_SR ***********************************
C
C  initialization for surface routing module in the case of no coupling
C
C***********************************************************************
C
      SUBROUTINE INIT_SR(NNOD,NCELNL,NUMRES,NSURFT,NSURFT_T,
     1     NSURFT_TB,OVFLNOD,OVFLP)
C
      IMPLICIT NONE
      INCLUDE 'CATHY.H'
      INTEGER  NNOD,NCELNL,NUMRES,NSURFT,NSURFT_T,NSURFT_TB
      REAL*8   OVFLNOD(*),OVFLP(*)

      INCLUDE 'SURFWATER.H'
C
      NSURFT = 0
      NSURFT_T = 0
      NSURFT_TB = 0
      ak_max = 0.0d0
      ak_max_p = 0.0d0
    
      CALL INIT0R(NNOD,OVFLNOD)
      CALL INIT0R(NNOD,OVFLP)

cms   The initialization  of variables in presence of reservoirs needs
cms   to be double-checked
      call init0r(NUMRES,h_pool_kk_vec)
      call init0r(NUMRES,h_pool_kk_vec_p)
 
      CALL INIT0R(MAXCEL,VOLUME_KK_SN)
      CALL INIT0R(MAXCEL,VOLUME_KKP1_SN)
      CALL INIT0R(MAXCEL,VOLUME_KK_SN_P)
      CALL INIT0R(MAXCEL,Q_IN_KK_SN)
      CALL INIT0R(MAXCEL,Q_IN_KK_SN_P)
      CALL INIT0R(MAXCEL,Q_IN_KKP1_SN)
      CALL INIT0R(MAXCEL,Q_OUT_KK_SN_1)
      CALL INIT0R(MAXCEL,Q_OUT_KK_SN_2)
      CALL INIT0R(MAXCEL,Q_OUT_KK_SN_1_P)
      CALL INIT0R(MAXCEL,Q_OUT_KK_SN_2_P)
      CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_1)
      CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_2)
      

      RETURN
      END
