C
C************************** SURF_FLOWTRA  ******************************
C
C surface water routing and solute transport procedure:
C 
C FOR FLOW: takes OVFLNOD as input from FLOW3D and returns `routed'
C           PONDNOD as BC for FLOW3D
C
C FOR TRANSPORT: takes TRAFLNOD as input from subsurface transport
C                and returns routed CONCNOD that is used to define the 
C                Cauchy BC for subsurface transport
C***********************************************************************
C
      subroutine surf_flowtra(ncell,nnod,nrow,ncol,ntri,dostep,numres,
     1                       ncell_coarse,ncelnl,cell,
     2                       indcel,indcelwl,tipo_r,reservr,dem_map,
     3                       lakes_map,cellcol,cellrow,
     4                       celtype,cells_r,tp2d,triang,
     5                       time,deltat,arenod,OVFLNOD,PONDNOD,
     6                       ovflcel,pondcel,TRAFLNOD,TRAFLCEL,CONCNOD,
     7                       CONCCEL,N_HA,NSURF,NSURFT,NSURFT_TB,TRANSP,
     8                       SOURCE_MIXING,MIX_CORRECT,SURFACE_MIX_SN)
C     
      implicit none
      include 'CATHY.H'
      INCLUDE 'SURFWATER.H'
      INCLUDE 'RIVERNETWORK.H'
      INCLUDE 'TRANSPSURF.H'
      LOGICAL  TRANSP
      integer  ncell,nnod,nrow,ncol,ntri,dostep,numres
      integer  ncell_coarse,ncelnl
      INTEGER  NSURF,NSURFT,NSURFT_TB
      integer  i,j,k
      integer  cell(5,maxcel),indcel(rowmax,colmax)
      integer  indcelwl(rowmax,colmax)
      integer  tipo_r(*),reservr(*)
      integer  lakes_map(rowmax,*)
      integer  cellcol(*),cellrow(*)
      integer  celtype(*),cells_r(*)
      integer  tp2d(*),triang(4,*)
      integer  n_hA(*)
      real*8   time,deltat
      real*8   dtsurf
      real*8   arenod(*)
      real*8   OVFLNOD(*),PONDNOD(*)
      real*8   TRAFLNOD(*),CONCNOD(*),SOURCE_MIXING(*)
      real*8   ovflcel(*),pondcel(*),dem_map(rowmax,*)
      real*8   traflcel(*),conccel(*),MIX_CORRECT(*)
      real*8   traflcel_mix(maxcel),surface_mix_sn(maxcel)
C
      call init0r(maxcel,traflcel_mix)
      call init0r(maxcel,surface_mix_sn)
c
c  from node to cell
c    
cm    call nod_cell(ncell,nrow,ncol,dostep,ncell_coarse,
cm   1              nnod,cell,dem_map,indcelwl,cellcoarse,OVFLNOD,
cm   2              ovflcel,arenod,delta_x,delta_y,pondnod)
      call nod_cell(ncell,nrow,ncol,dostep,ncell_coarse,
     1              nnod,cell,dem_map,indcelwl,cellcoarse,OVFLNOD,
     2              ovflcel,arenod,delta_x,delta_y)
c    
c   trasferimento informazioni al surf_route (cambia la numerazione
c   delle celle se ci sono laghi!)
c     
      call transfer_f3d_surf(nrow,ncol,indcel,indcelwl,
     1                       surface_water_sn,ovflcel,
     2                       reservr,lakes_map)
c   
      IF (TRANSP) THEN 
         call nod_cell(ncell,nrow,ncol,dostep,ncell_coarse,
     1        nnod,cell,dem_map,indcelwl,cellcoarse,traflnod,
     2        traflcel,arenod,delta_x,delta_y)

         call nod_cell(ncell,nrow,ncol,dostep,ncell_coarse,
     1        nnod,cell,dem_map,indcelwl,cellcoarse,source_mixing,
     2        traflcel_mix,arenod,delta_x,delta_y)

         call transfer_f3d_surf(nrow,ncol,indcel,indcelwl,
     1        surface_conc_sn,traflcel,
     2        reservr,lakes_map)

         call transfer_f3d_surf(nrow,ncol,indcel,indcelwl,
     1        surface_mix_sn,traflcel_mix,
     2        reservr,lakes_map)
      END IF
c    
c  surface routing
c  
      cu_max=ak_max*deltat
      CUTRGT = 1.0d0
      if(cu_max .gt. CUTRGT) then
         dtsurf = CUTRGT/ak_max
         nsurf = IDINT((deltat/dtsurf))+1
         dtsurf = deltat/nsurf
      else
         dtsurf = deltat
         nsurf=1
      end if
      NSURFT=NSURFT + NSURF
      NSURFT_TB=NSURFT_TB + NSURF
C
      DO I=1,NSURF
C------------------------WATER ROUTING -----------------------------
         CALL ROUTE(NROW,NCELL,TIPO_R,RESERVR,
     1        CELLS_R,N_HA,TIME,DTSURF)
C     CALCULATION OF SURFACE WATER HEIGHTS
         CALL ALTEZZE(NCELNL,CELTYPE,DELTA_X,
     1        DELTA_Y,Q_IN_KK_SN,Q_IN_KKP1_SN,
     2        Q_OUT_KK_SN_1,Q_OUT_KK_SN_2,
     3        Q_OUT_KKP1_SN_1,Q_OUT_KKP1_SN_2,
     4        VOLUME_KK_SN,VOLUME_KKP1_SN,SURFACE_WATER_SN,
     5        H_WATER_KKP1_SN,DTSURF,
     6        H_POOL_KKP1_VEC,RESERVR,ELEVATION)
         
C     
C---------------------------SURFACE TRANSPORT -------------------------
         IF (TRANSP) THEN 
C     CONCENTRATION ROUTING
            CALL ROUTE_TRA(NROW,NCELL,TIPO_R,RESERVR,
     1           CELLS_R,N_HA,TIME,DTSURF)
c            
C     CALCULATION OF CONCENTRATION FOR CELLS
            CALL MASBIL(NCELL,TIPO_R,DELTA_X,DELTA_Y,
     1           QMASS_IN_KK_SN,QMASS_IN_KKP1_SN,
     2           QMASS_OUT_KK_SN_1,QMASS_OUT_KK_SN_2,
     3           QMASS_OUT_KKP1_SN_1,QMASS_OUT_KKP1_SN_2,
     4           MASS_KK_SN,MASS_KKP1_SN,SURFACE_CONC_SN,
     5           CONC_KKP1_SN,VOLUME_KKP1_SN,DTSURF,MIX_CORRECT,
     6           surface_mix_sn)  
c
         END IF
C-----------------------------------------------------------------------



C INTERNAL UPDATES OF SURFACE ROUTING VARIABLES (I.E., FLOW RATE, WATER HEIGHTS)
            IF(NSURF.GT.1 .AND. I.LT.NSURF) THEN
               CALL VCOPYR(MAXCEL,Q_IN_KK_SN,Q_IN_KKP1_SN)
               CALL INIT0R(MAXCEL,Q_IN_KKP1_SN)
               CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_1,Q_OUT_KKP1_SN_1)
               CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_1)
               CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_2,Q_OUT_KKP1_SN_2)
               CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_2)
               CALL VCOPYR(MAXCEL,VOLUME_KK_SN,VOLUME_KKP1_SN)
               CALL INIT0R(MAXCEL,VOLUME_KKP1_SN)
C     
               IF (TRANSP) THEN
                  CALL VCOPYR(MAXCEL,QMASS_IN_KK_SN,QMASS_IN_KKP1_SN)
                  CALL INIT0R(MAXCEL,QMASS_IN_KKP1_SN)
                  CALL VCOPYR(MAXCEL,QMASS_OUT_KK_SN_1,
     1                 QMASS_OUT_KKP1_SN_1)
                  CALL INIT0R(MAXCEL,QMASS_OUT_KKP1_SN_1)
                  CALL VCOPYR(MAXCEL,QMASS_OUT_KK_SN_2,
     1                 QMASS_OUT_KKP1_SN_2)
                  CALL INIT0R(MAXCEL,QMASS_OUT_KKP1_SN_2)
                  CALL VCOPYR(MAXCEL,MASS_KK_SN,MASS_KKP1_SN)
                  CALL INIT0R(MAXCEL,MASS_KKP1_SN)
               END IF
            END IF
         END DO
                        
c
c trasferimento informazioni per il flow3d
c
      call transfer_surf_f3d(nrow,ncol,indcel,indcelwl,
     1                       ncell,h_water_kkp1_sn,
     2                       pondcel,h_pool_kkp1_vec,
     3                       lakes_map,elevation_with_lakes)
c
c  calcola PONDNOD da PONDCEL passando per i triangoli
c
cm    call cell_nod(ncell,nrow,ncol,nnod,ntri,dostep,
cm   1                    cell,tp2d,triang,indcel,dem_map,
cm   2                    PONDCEL,cellcoarse,pondnod,
cm   3                    delta_x,delta_y,arenod)
      call cell_nod(ncell,nrow,ncol,nnod,ntri,dostep,
     1              tp2d,triang,indcelwl,dem_map,
     2              pondcel,cellcoarse,PONDNOD)
c
      IF (TRANSP) THEN 
         call transfer_surf_f3d(nrow,ncol,indcel,indcelwl,
     1        ncell,CONC_KKP1_SN,
     2        CONCCEL,h_pool_kkp1_vec,
     3        lakes_map,elevation_with_lakes)
c           
         call cell_nod_tra(ntri,nnod,tp2d,triang,CONCCEL,concnod,
     3                    delta_x,delta_y,arenod,pondcel,pondnod)
cm       call cell_nod(ncell,nrow,ncol,nnod,ntri,dostep,
cm   1        tp2d,triang,indcelwl,dem_map,
cm   2        CONCCEL,cellcoarse,CONCNOD)
c
      END IF
      return
      end

