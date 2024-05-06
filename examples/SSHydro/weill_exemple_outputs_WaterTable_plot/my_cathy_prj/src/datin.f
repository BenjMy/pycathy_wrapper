C
C**************************  DATIN  ************************************
C
C  read in some of the input data (other data is read in subroutines
C  BCONE, BCNXT, ATMONE, ATMNXT, NUDONE, NUDNXT, EFFONE, EFFNXT, and
C  RAST_INPUT)
C
C***********************************************************************
C
      SUBROUTINE DATIN(ISIMGR,IVERT,ISP,WTPOSITION,BASE,ZRATIO, 
     1                 TRIANG,X,Y,Z,PERMX,PERMY,PERMZ,ELSTOR,
     2                 POROS,VGNCELL,VGRMCCELL,VGPSATCELL,CONTR,NODVP,
     3                 ID_QOUT,TIMPRT,PONDNOD,PTIMEP,TETAF,
     4                 DELTAT,DTMIN,DTMAX,TMAX,DTMAGA,DTMAGM,
     5                 DTREDS,DTREDM,ITUNS,ITUNS1,ITUNS2,
     6                 TOLUNS,TOLSWI,ERNLMX,ITMXCG,TOLCG,
     7                 KSLOPE,LUMP,IPEAT,IVGHU,NLKP,VTKF,
     8                 IOPT,ISOLV,IPRT1,IPRT,IPOND,INDP,NNOD,NTRI,NSTR,
     9                 NZONE,NVEG,N1,NR,NUMVP,NUM_QOUT,NPRT,N,NT,
     A                 ISFONE,ISFCVG,DUPUIT,
     B                 L2NORM,NLRELX,OMEGA,
     C                 PONDH_MIN,
     D                 FL3D,SURF,DEM,GRID,NROW,NCOL,
     E                 NCELNL,NCELL,DOSTEP,NCELL_COARSE,
     F                 DEM_MAP,ZONE,LAKES_MAP,INDEX,INDEX_WITH_LAKES,
     G                 CELL,CELLCOL,CELLROW,TP2D,NODI,
     H                 TIPO_R,RESERVR,N_HA,CELLCOL_WL,CELLROW_WL,
     I                 BASE_MAP,DEPTH,ELTRIA,NCOUT,TRAFLAG,TRANSP,
     K                 VELREC,VEG_TYPE)
C
      IMPLICIT  NONE
      INCLUDE  'CATHY.H'
      INTEGER   IROW,ICOL,INDEX_C,ILAKE
      INTEGER   I,J,K,FLAGASC,PP
      INTEGER   ISIMGR,IVERT,ISP
      INTEGER   ITUNS,ITUNS1,ITUNS2,ITMXCG,KSLOPE,LUMP,IPEAT
      INTEGER   IVGHU,IOPT,ISOLV,IPRT1,IPRT,IPOND,INDP,NZONE
      INTEGER   NNOD,NTRI,NSTR,NVEG
      INTEGER   NLKP,VTKF,VELREC
      INTEGER   N1,NR,NUMVP,NPRT,N,NT,ISFONE,ISFCVG,DUPUIT
      INTEGER   NUM_QOUT,L2NORM,NLRELX
      INTEGER   NROW,NCOL,NCELNL,NCELL,NUM_TOT_R
      INTEGER   DOSTEP,NCELL_COARSE,NCOUT,TRAFLAG
      INTEGER   TRIANG(4,*)
      INTEGER   CONTR(*),NODVP(*),ID_QOUT(*)
      INTEGER   CELLCOL(*),CELLROW(*),TP2D(*)
      INTEGER   CELLCOL_WL(*),CELLROW_WL(*)
      INTEGER   LAKES_MAP(ROWMAX,*)
      INTEGER   ZONE(ROWMAX,*)
      INTEGER   INDEX(ROWMAX,*),INDEX_WITH_LAKES(ROWMAX,*)
      INTEGER   NODI(ROWMAX+1,*),CELL(5,*)
      INTEGER   RESERVR(*)
      INTEGER   TIPO_R(*)
      INTEGER   N_HA(*)
      INTEGER   VEG_TYPE(NODMAX)
      LOGICAL   GRID,DEM,FL3D,SURF,TRANSP
      REAL*8    SUMZ
      REAL*8    WTPOSITION,BASE,TETAF,DELTAT,DTMIN,DTMAX,TMAX
      REAL*8    DTMAGA,DTMAGM,DTREDS,DTREDM
      REAL*8    TOLUNS,TOLSWI,ERNLMX,TOLCG
      REAL*8    OMEGA
      REAL*8    PONDH_MIN
      REAL*8    ZRATIO(*),X(*),Y(*),Z(*),DEPTH(*)
      REAL*8    ELTRIA(NTRMAX)
      REAL*8    PERMX(MAXSTR,*),PERMY(MAXSTR,*),PERMZ(MAXSTR,*)
      REAL*8    ELSTOR(MAXSTR,*),POROS(MAXSTR,*)
      REAL*8    VGNCELL(MAXSTR,*),VGRMCCELL(MAXSTR,*)
      REAL*8    VGPSATCELL(MAXSTR,*)
      REAL*8    TIMPRT(*),PONDNOD(*),PTIMEP(*)
      REAL*8    DEM_MAP(ROWMAX,*),BASE_MAP(ROWMAX,*)
      REAL*8    VEG_MAP(ROWMAX,COLMAX),SCR(NODMAX)

      INCLUDE  'IOUNITS.H'
      INCLUDE  'SOILCHAR.H'
      INCLUDE  'SURFWATER.H'
      INCLUDE  'RIVERNETWORK.H'
C
C  unit IIN1 input 
C
      READ(IIN1,*) IPRT1,NCOUT,TRAFLAG
      READ(IIN1,*) ISIMGR,PONDH_MIN,VELREC
      IF (ISIMGR .EQ. 0) THEN
         FL3D=.TRUE.
         SURF=.FALSE.
         DEM=.FALSE.
         GRID=.TRUE.
      ELSE IF (ISIMGR .EQ. 1) THEN
         FL3D=.TRUE.
         SURF=.FALSE.
         DEM=.TRUE.
         GRID=.FALSE.
      ELSE IF (ISIMGR .EQ. 2) THEN
         FL3D=.TRUE.
         SURF=.TRUE.
         DEM=.TRUE.
         GRID=.FALSE.
      ELSE IF (ISIMGR .EQ. 3) THEN
         FL3D=.FALSE.
         SURF=.TRUE.
         DEM=.TRUE.
         GRID=.FALSE.   
      ELSE
         WRITE(IOUT2,1900) ISIMGR
         CALL CLOSIO
         STOP
      END IF
      READ(IIN1,*) KSLOPE,TOLKSL
      READ(IIN1,*) PKRL,PKRR,PSEL,PSER
      READ(IIN1,*) PDSE1L,PDSE1R,PDSE2L,PDSE2R
      READ(IIN1,*) ISFONE,ISFCVG,DUPUIT
      READ(IIN1,*) TETAF,LUMP,IOPT
      READ(IIN1,*) NLRELX,OMEGA
      READ(IIN1,*) L2NORM,TOLUNS,TOLSWI,ERNLMX
      READ(IIN1,*) ITUNS,ITUNS1,ITUNS2
      READ(IIN1,*) ISOLV,ITMXCG,TOLCG
      READ(IIN1,*) DELTAT,DTMIN,DTMAX,TMAX
      READ(IIN1,*) DTMAGA,DTMAGM,DTREDS,DTREDM
      READ(IIN1,*) IPRT,VTKF,NPRT,(TIMPRT(I),I=1,NPRT)
      READ(IIN1,*) NUMVP,(NODVP(I),I=1,NUMVP)
      READ(IIN1,*) NR
      IF (NR .NE. 0) READ(IIN1,*) (CONTR(I),I=1,NR)
      READ(IIN1,*) NUM_QOUT,(ID_QOUT(I),I=1,NUM_QOUT)
      if (isimgr.le.1) pondh_min=1.0d+10
      WRITE(IOUT2,1000) ISIMGR,PONDH_MIN
      WRITE(IOUT2,1010) IPRT1,IPRT,NPRT,NUMVP,NR
      WRITE(IOUT2,1015) KSLOPE,TOLKSL
      IF (KSLOPE .EQ. 3  .OR.  KSLOPE .EQ. 4) WRITE(IOUT2,1040)
     1                                      PKRL,  PKRR,  PSEL,  PSER,
     2                                      PDSE1L,PDSE1R,PDSE2L,PDSE2R
      WRITE(IOUT2,1045) ISFONE,ISFCVG
      WRITE(IOUT2,1004) TETAF,LUMP,IOPT
      WRITE(IOUT2,1050) NLRELX
      IF (NLRELX .EQ. 1) WRITE(IOUT2,1055) OMEGA
      WRITE(IOUT2,1060) L2NORM,TOLUNS,TOLSWI,ERNLMX
      WRITE(IOUT2,1025) ITUNS,ITUNS1,ITUNS2
      WRITE(IOUT2,1005) ISOLV,ITMXCG,TOLCG
      WRITE(IOUT2,1030) DELTAT,DTMIN,DTMAX,TMAX
      WRITE(IOUT2,1035) DTMAGA,DTMAGM,DTREDS,DTREDM
C
C  unit IIN2 input; if (grid) we just read the grid file; if (dem) we
C                   build up the surface mesh starting from the DEM
C                   of the basin.
C  
      IF (GRID) THEN 
         READ(IIN2,*) NZONE,NVEG,NSTR,N1
         READ(IIN2,*) NNOD,NTRI
         WRITE(IOUT2,1020) NNOD,NTRI,NZONE,NVEG,NSTR,N1
         READ(IIN2,*) IVERT,ISP,BASE
         WRITE(IOUT2,1350) IVERT,ISP,BASE
         READ(IIN2,*) (ZRATIO(I),I=1,NSTR)
         SUMZ=0.0D0
         DO I=1,NSTR
            SUMZ=SUMZ + ZRATIO(I)
         END DO
         WRITE(IOUT2,1340) (I,ZRATIO(I),I=1,NSTR)
         IF (DABS(SUMZ - 1.0D0) .GT. 1.0D-14) THEN
            WRITE(IOUT2,1345) 
            WRITE(IOUT2,*) SUMZ
            CALL CLOSIO
            STOP
         END IF
         IF (ISP .NE. 0) THEN
            READ(IIN2,*) (Z(I),I=1,NNOD)
         ELSE
            READ(IIN2,*) Z(1)
            DO I=2,NNOD
               Z(I)=Z(1)
            END DO
         END IF
         IF (IPRT1 .GE .1) WRITE(IOUT2,1360) (I,Z(I),I=1,NNOD)
         READ(IIN2,*) ((TRIANG(I,K),I=1,4),K=1,NTRI)
         READ(IIN2,*) (X(K),Y(K),K=1,NNOD)
         IF (ISP .GE. 2) THEN
            READ(IIN2,*) (VEG_TYPE(I),I=1,NNOD)
         ELSE
            READ(IIN2,*) VEG_TYPE(1)
            DO I=2,NNOD
               VEG_TYPE(I)=VEG_TYPE(1)
            END DO
         END IF
         IF (NVEG.GT.MAXVEG) THEN
            WRITE(IOUT2,*) 'Error: NVEG is too large=',NVEG
            CALL CLOSIO
            STOP
         END IF
         
      ELSE IF (DEM) THEN
        
         CALL RAST_INPUT_DEM(IIN10,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DEM_MAP)
    
         CALL RAST_INPUT_LZ(IIN21,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        ZONE)
        
         READ(IIN11,*) DELTA_X,DELTA_Y
         WRITE(IOUT40,*) 'DELTA_X=',DELTA_X,'DELTA_Y=',DELTA_Y
         READ(IIN11,*) FACTOR
         WRITE(IOUT40,*)'FACTOR=',FACTOR
         READ(IIN11,*) DOSTEP
         WRITE(IOUT40,*) 'DOSTEP=',DOSTEP                  
        
         IF (FL3D) THEN
            READ(IIN11,*) NZONE,NSTR,N1 
            READ(IIN11,*) IVERT,ISP,BASE
            IF (IVERT.EQ.3) THEN
               CALL RAST_INPUT_DEM(IIN60,NROW,NCOL,NORTH,SOUTH,EAST,
     1              WEST,BASE_MAP)
               CALL TRIANGOLI(NROW,NCOL,DELTA_X,DELTA_Y,WEST,SOUTH,
     1              BASE_MAP,ZONE,FACTOR,NNOD,NTRI,DOSTEP,
     2              NCELL_COARSE,NODI,TRIANG,TP2D,X,Y,DEPTH,ELTRIA,
     3              CELL)
            END IF
         END IF
         
c
c  reading of lakes_map
c              
         CALL RAST_INPUT_LZ(IIN20,NROW,NCOL,NORTH,SOUTH,EAST,
     1        WEST,LAKES_MAP)
c
c  reading of vegetation type (real raster map, then converted to int)
c
         CALL RAST_INPUT_DEM(IIN3,NROW,NCOL,NORTH,SOUTH,EAST,
     1                       WEST,VEG_MAP)
         CALL TRIANGOLI(NROW,NCOL,DELTA_X,DELTA_Y,WEST,SOUTH,
     1                  VEG_MAP,ZONE,FACTOR,NNOD,NTRI,DOSTEP,
     2                  NCELL_COARSE,NODI,TRIANG,TP2D,X,Y,SCR,
     3                  ELTRIA,CELL)
         NVEG=0
         DO I=1,NNOD
            VEG_TYPE(I)=INT(SCR(I))
            NVEG=MAX(NVEG,VEG_TYPE(I))
         END DO
         IF (NVEG.GT.MAXVEG) THEN
            WRITE(IOUT2,*) 'Error: NVEG is too large=',NVEG
            CALL CLOSIO
            STOP
         END IF
c
c  cells numbering with and without lakes
c
        CALL INDEX_DEM(ROWMAX,MAXCEL,NROW,NCOL,NCELNL,NCELL,
     1       ZONE,INDEX,CELLCOL,CELLROW,
     2       LAKES_MAP,INDEX_WITH_LAKES,
     3       CELLCOL_WL,CELLROW_WL)
        WRITE(IOUT40,*) 'NCELNL=',NCELNL
        WRITE(IOUT40,*) 'NCELL=',NCELL
c
c  from cells to triangles: construction of triangles, numbering of 
c  nodes ed elements, assignment of x y z coordinates, assignment of elevation
c  values to triangles
c
         IF (FL3D) THEN
            CALL ASSIGN_DEM(NROW,NCOL,DEM_MAP,INDEX,INDEX_WITH_LAKES,
     1           LAKES_MAP,FACTOR,ELEVATION,ELEVATION_WITH_LAKES)
            
            CALL TRIANGOLI(NROW,NCOL,DELTA_X,DELTA_Y,WEST,SOUTH,
     1           DEM_MAP,ZONE,FACTOR,NNOD,NTRI,DOSTEP,
     2           NCELL_COARSE,NODI,TRIANG,TP2D,X,Y,Z,ELTRIA,
     3           CELL)        
c
c     reading of grid parameters
c
            WRITE(IOUT40,1020) NNOD,NTRI,NZONE,NVEG,NSTR,N1
            WRITE(IOUT40,1350) IVERT,ISP,BASE
            READ(IIN11,*) (ZRATIO(I),I=1,NSTR)
            SUMZ=0.0D0
            DO I=1,NSTR
               SUMZ=SUMZ + ZRATIO(I)
            END DO
            WRITE(IOUT40,1340) (I,ZRATIO(I),I=1,NSTR)
            IF (DABS(SUMZ - 1.0D0) .GT. 1.0D-14) THEN
               WRITE(IOUT2,1345)
               CALL CLOSIO
               STOP
            END IF
c     IF (IPRT1 .GE .1) WRITE(IOUT40,1360) (I,Z(I),I=1,NNOD)
         END IF
      END IF
      N=NNOD*(NSTR + 1)
      NT=3*NTRI*NSTR

      IF (SURF) THEN
         read(IIN17,*) num_tot_r
         write(IOUT40,*) 'numero serbatoi=',num_tot_r
         if (num_tot_r.ne.0) then
            write(IOUT40,*) 'index_c  tipo_r(index_c) reservr(index_c)'
            do j=1,ncelnl
               read(IIN17,*) i,tipo_r(i),reservr(i)
            end do
            write(IOUT40,*) 'livelli iniziali nei serbatoi'
c            do j=1,num_tot_r
c               read(IIN18,*) i,h_pool_kkp1_vec(i)
c               h_pool_kkp1_vec(i)=0
c               write(IOUT40,*) i,h_pool_kkp1_vec(i)
c            end do
         else
            do index_c=1,maxcel
               tipo_r(index_c)=0
               reservr(index_c)=0
            end do
         end if
         if (num_tot_r.ne.0) then
            write(iout40,*) 'caratteristiche',num_tot_r,' serbatoi'
            do i=1,num_tot_r
               read(iin18,*)ilake,n_hA(ilake)
               write(iout40,*) ilake,n_hA(ilake)
               pp=n_hA(ilake)
               write(iout40,*) 'n_ha=',pp
               do k=1,pp
                  read(iin18,*)hres(ilake,k),Ares(ilake,k)
                  write(iout40,*) hres(ilake,k),Ares(ilake,k)   
               end do
               h_fondo(ilake)=hres(ilake,2)
               write(iout40,*) 'h_fondo=',h_fondo(ilake)
               read(iin18,*) h_soglia(ilake)
               write(iout40,*) 'h_soglia=',h_soglia(ilake)
               read(iin18,*) h_pool_kkp1_vec(i)
ccc            h_pool_kkp1_vec(i)=h_fondo(ilake)
               write(iout40,*) 'liv.iniz.=',h_pool_kkp1_vec(i)
            end do
         end if
                 
         
         READ(IIN23,*)
         DO I=1,NCELL
            READ(IIN23,*) QOI_SN(I)
         END DO
    
         CALL RAST_INPUT_REAL(IIN25,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_W_1)

         CALL RAST_INPUT_REAL(IIN26,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_W_2)
        
         CALL RAST_INPUT_INT(IIN27,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_P_OUTFLOW_1)
        
         CALL RAST_INPUT_INT(IIN28,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_P_OUTFLOW_2)

         CALL RAST_INPUT_REAL(IIN29,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_LOCAL_SLOPE_1)
         
         CALL RAST_INPUT_REAL(IIN30,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_LOCAL_SLOPE_2)

         CALL RAST_INPUT_REAL(IIN31,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_EPL_1)
         
         CALL RAST_INPUT_REAL(IIN32,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_EPL_2)
               
         CALL RAST_INPUT_REAL(IIN33,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_KSS1_SF_1)
         
         CALL RAST_INPUT_REAL(IIN34,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_KSS1_SF_2)

         CALL RAST_INPUT_REAL(IIN35,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_WS1_SF_1)
         
         CALL RAST_INPUT_REAL(IIN36,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_WS1_SF_2)

         CALL RAST_INPUT_REAL(IIN37,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_B1_SF)

         CALL RAST_INPUT_REAL(IIN38,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_Y1_SF)  

         CALL RAST_INPUT_REAL(IIN39,NROW,NCOL,NORTH,SOUTH,EAST,WEST,
     1        DTM_NRC)
      END IF
C     
      IF (FL3D) THEN
C
C     unit IIN5 input
C
         READ(IIN5,*) INDP,IPOND
         IF (.NOT. SURF) IPOND = 0
         IF ((INDP .EQ. 3).OR.(INDP .EQ. 4)) THEN
            READ(IIN5,*) WTPOSITION
         END IF
         IF (INDP .EQ. 0) THEN
            READ(IIN5,*) PTIMEP(1)
            WRITE(IOUT2,1070) PTIMEP(1)
            DO K=2,N
               PTIMEP(K)=PTIMEP(1)
            END DO
         ELSE IF (INDP .EQ. 1) THEN
            READ(IIN5,*) (PTIMEP(K),K=1,N)
         END IF
         IF (IPOND .EQ. 1) THEN
            READ(IIN5,*) PONDNOD(1)
            WRITE(IOUT2,1071) PONDNOD(1)   
            IF (PONDNOD(1) .GT. 0.0D0) PTIMEP(1) = PONDNOD(1)
            DO K=2,NNOD
               PONDNOD(K)=PONDNOD(1)   
               IF (PONDNOD(K) .GT. 0.0D0) PTIMEP(K) = PONDNOD(K)
            END DO
         ELSE IF (IPOND .EQ. 2) THEN   
            READ(IIN5,*) (PONDNOD(K),K=1,NNOD)   
            DO K=1,NNOD
               IF (PONDNOD(K) .GT. 0.0D0) PTIMEP(K) = PONDNOD(K)
            END DO
            IF (INDP .EQ. 0 .AND. IPRT1 .GE. 1) THEN
               WRITE(IOUT2,1072) (K,PONDNOD(K),K=1,NNOD)
            END IF
         END IF
         IF (INDP .EQ. 1 .AND. IPRT1 .GE. 1) THEN
            IF (IPOND.EQ.0) THEN
               WRITE(IOUT2,1090) (K,PTIMEP(K),K=1,N)
            ELSE
               WRITE(IOUT2,1091) (K,PTIMEP(K),K=1,N)
            END IF
         END IF
C
C     unit IIN4 input
C
         READ(IIN4,*) PMIN
         READ(IIN4,*) IPEAT,SCF
         READ(IIN4,*) CBETA0,CANG
C    variable vegetation type (NVEG) within the domain
         DO I=1,NVEG
            READ(IIN4,*) PCANA(I),PCREF(I),PCWLT(I),ZROOT(I),
     1                   PZ(I),OMGC(I)
         END DO
c     peat soil deformation is not yet supported for the Newton
c     scheme
         IF (IOPT .NE. 1  .AND.  IPEAT .EQ. 1) THEN
            WRITE(IOUT2,1910)
            CALL CLOSIO
            STOP
         END IF
         READ(IIN4,*) IVGHU
c     peat soil deformation is not yet supported for the chord slope
c     and tangent slope differentiation of moisture curves
         IF (KSLOPE .NE. 0  .AND.  IPEAT .EQ. 1) THEN
            WRITE(IOUT2,1930)
            CALL CLOSIO
            STOP
         END IF
c     peat soil deformation is not yet supported for the extended
c     van Genuchten or Huyakorn moisture curves (unless they are
c     given in lookup table form)
         IF ( (IPEAT .EQ. 1)  .AND.
     1        (IVGHU .EQ. 1 .OR. IVGHU .EQ. 2 .OR. IVGHU .EQ. 3) ) THEN
            WRITE(IOUT2,1940)
            CALL CLOSIO
            STOP
         END IF
c     for IVGHU = -1 the following "VG"/"HU"/"BC" parameters are not
c     needed but are read in anyway. Arbitrary values can be assigned.
         READ(IIN4,*) HUALFA,HUBETA,HUGAMA,HUPSIA,HUSWR
         READ(IIN4,*) HUN
         READ(IIN4,*) HUA,HUB
         READ(IIN4,*) BCBETA,BCRMC,BCPSAT
C
         IF (IVGHU .EQ. -1) THEN
            READ(IIN16,*) NLKP
            IF (NLKP .LT. 3) THEN
               WRITE(IOUT2,1960)
               CALL CLOSIO
               STOP
            END IF
            DO I=1,NSTR
               DO J=1,NZONE
                  DO K=1,NLKP
                     READ(IIN16,*) PCAP(I,J,K),SATC(I,J,K),KRWC(I,J,K)
                  END DO
               END DO
            END DO
            FLAGASC=0
            DO I=1,NSTR
               DO J=1,NZONE
                  DO K=2,NLKP
                     IF (PCAP(I,J,K) .LE. PCAP(I,J,K-1)) FLAGASC=1
                  END DO
               END DO
            END DO
            IF (FLAGASC .EQ. 1) THEN
               WRITE(IOUT2,1970)
               CALL CLOSIO
               STOP
            END IF
         END IF
C
         WRITE(IOUT2,1215) PMIN
         WRITE(IOUT2,1216) IPEAT
         IF (IPEAT .EQ. 1) THEN
            WRITE(IOUT2,1217) CBETA0,CANG
         END IF
         WRITE(IOUT2,1220) IVGHU
         IF (IVGHU .EQ. 0 .OR. IVGHU .EQ. 1) THEN
CM          WRITE(IOUT2,1230) VGN,VGRMC,VGPSAT
            WRITE(IOUT2,*) 'SPATIALLY VARIABLE VAN GENUCHTEN PARAMETERS'
         ELSE IF (IVGHU .EQ. 2 .OR. IVGHU .EQ. 3) THEN
            WRITE(IOUT2,1240) HUALFA,HUBETA,HUGAMA,HUPSIA,HUSWR
            IF (IVGHU .EQ. 2) THEN
               WRITE(IOUT2,1250) HUN
            ELSE
               WRITE(IOUT2,1260) HUA,HUB
            END IF
         ELSE IF (IVGHU .EQ. 4) THEN
            WRITE(IOUT2,1235) BCBETA,BCRMC,BCPSAT
         END IF
C
         WRITE(IOUT2,1100)
         DO I=1,NSTR
            DO J=1,NZONE
              READ(IIN4,*) PERMX(I,J),PERMY(I,J),PERMZ(I,J),ELSTOR(I,J),
     1              POROS(I,J),VGNCELL(I,J),VGRMCCELL(I,J),
     2              VGPSATCELL(I,J)

            IF(IPRT1.GE.2) WRITE(IOUT2,1110) I,J,PERMX(I,J), PERMY(I,J),
     1              PERMZ(I,J),ELSTOR(I,J),POROS(I,J),VGNCELL(I,J),
     2              VGRMCCELL(I,J),VGPSATCELL(I,J)

            END DO
         END DO
      END IF

666   CONTINUE
C
      WRITE(IOUT2,1300) N,NT
C
C SET OF THE TRANSP FLAG, THE READING OF TRANSPORT PARAMETERS WILL 
C  DONE IN DATIN_TRA
      IF (TRAFLAG .EQ. 1) THEN
         TRANSP = .TRUE.
      ELSE
         TRANSP = .FALSE.
      END IF

      RETURN
C
C  format statements
C
 1000 FORMAT(/,5X,'ISIMGR (0 FLOW3D only w/ grid input, ',
     1       /,5X,'        1 FLOW3D only w/ DEM input, ',
     2       /,5X,'        2 FLOW3D and SURF_ROUTE w/ DEM) = ',I6,
     3       /,5X,'PONDH_MIN (MIN. PONDING HEAD)           = ',1PE15.5)
 1004 FORMAT(  5X,'TETAF  (EG: 1 BACKWARD EULER, 0.5 C-N)  = ',1PE15.5,
     1       /,5X,'LUMP   (MASS LUMPING IF NONZERO)        = ',I6,
     2       /,5X,'IOPT   (1 PICARD, 2 NEWTON)             = ',I6)
 1005 FORMAT(/,5X,'ISOLV  (-5 BiCGSTAB w/ diag precond, ',
     1       /,5X,'        -4 BiCGSTAB without precond, ',
     2       /,5X,'        -3 TFQMR   w/ diag precond, ',
     3       /,5X,'        -2 TFQMR   without precond, ',
     4       /,5X,'        -1 TFQMR   w/ K^-1  precond, ',
     5       /,5X,'         0 BiCGSTAB w/ K^-1  precond, ',
     6       /,5X,'         1 GRAMRB (min residual), ',
     7       /,5X,'         2 GCRK(5) (ORTHOMIN), ',
     8       /,5X,'         3 NONSYM (direct solver))      = ',I6,
     9       /,5X,'ITMXCG (MAX ITER FOR CG LINEAR SOLVERS) = ',I6,
     A       /,5X,'TOLCG  (TOLER. FOR CG LINEAR SOLVERS)   = ',1PE15.5)
 1010 FORMAT(  5X,'IPRT1  (FOR OUTPUT OF INPUT DATA)       = ',I6,
     1       /,5X,'IPRT   (FOR ELEM VEL & DET NODAL OUTPUT)= ',I6,
     2       /,5X,'NPRT   (# OF TIME VALUES FOR DET OUTPUT)= ',I6,
     3       /,5X,'NUMVP  (# OF SURF NODES FOR VP OUTPUT)  = ',I6,
     4       /,5X,'NR     (# OF NODES FOR PARTIAL OUTPUT)  = ',I6)
 1015 FORMAT(/,5X,'KSLOPE (0 ANA, 1 KSL/ANA, 2 KSL/C-DIFF,',
     1       /,5X,'        3 LOC KSL/ANA, 4 LOC TAN-SLOPE) = ',I6,
     2       /,5X,'TOLKSL (TOLERANCE FOR CHORD SLOPE)      = ',1PE15.5)
 1020 FORMAT(/,5X,'NNOD   (# OF NODES IN 2-D MESH)         = ',I6,
     1       /,5X,'NTRI   (# OF TRIANGLES IN 2-D MESH)     = ',I6,
     2       /,5X,'NZONE  (NUMERO ZONE (MATERIAL TYPES))   = ',I6,
     3       /,5X,'NVEG   (NUMERO ZONE VEG (VEG TYPES))    = ',I6,
     4       /,5X,'NSTR   (NUMERO STRATI)                  = ',I6,
     5       /,5X,'N1     (NUM. MAX CONTATTI NODALI)       = ',I6)
 1025 FORMAT(  5X,'ITUNS  (MAX NONLINEAR ITER / TIME STEP) = ',I6,
     1       /,5X,'ITUNS1 (DELTAT INCREASE THRESHOLD)      = ',I6,
     2       /,5X,'ITUNS2 (DELTAT DECREASE THRESHOLD)      = ',I6)
 1030 FORMAT(/,5X,'DELTAT (INITIAL TIME STEP SIZE)         = ',1PE15.5,
     1       /,5X,'DTMIN  (MINIMUM TIME STEP SIZE)         = ',1PE15.5,
     2       /,5X,'DTMAX  (MAXIMUM TIME STEP SIZE)         = ',1PE15.5,
     3       /,5X,'TMAX   (TIME AT END OF SIMULATION)      = ',1PE15.5)
 1035 FORMAT(  5X,'DTMAGA (MAG. FACTOR FOR DELTAT, ADD.)   = ',1PE15.5,
     1       /,5X,'DTMAGM (MAG. FACTOR FOR DELTAT, MULT.)  = ',1PE15.5,
     2       /,5X,'DTREDS (RED. FACTOR FOR DELTAT, SUB.)   = ',1PE15.5,
     3       /,5X,'DTREDM (RED. FACTOR FOR DELTAT, MULT.)  = ',1PE15.5)
 1036 FORMAT(  5X,'DTOUT  (TIME STEP FOR OUTPUTS)          = ',1PE15.5)
 1040 FORMAT(  5X,'PKRL   (LEFT  PrHead ENDPT FOR d(Kr)/dP)= ',1PE15.5,
     1       /,5X,'PKRR   (RIGHT PrHead ENDPT FOR d(Kr)/dP)= ',1PE15.5,
     2       /,5X,'PSEL   (LEFT  PrHead ENDPT FOR d(Se)/dP)= ',1PE15.5,
     3       /,5X,'PSER   (RIGHT PrHead ENDPT FOR d(Se)/dP)= ',1PE15.5,
     4       /,5X,'PDSE1L (LEFT  P FOR dd(Se)/dPP, RANGE 1)= ',1PE15.5,
     5       /,5X,'PDSE1R (RIGHT P FOR dd(Se)/dPP, RANGE 1)= ',1PE15.5,
     6       /,5X,'PDSE2L (LEFT  P FOR dd(Se)/dPP, RANGE 2)= ',1PE15.5,
     7       /,5X,'PDSE2R (RIGHT P FOR dd(Se)/dPP, RANGE 2)= ',1PE15.5)
 1045 FORMAT(/,5X,'ISFONE (0 ALL-NODE SFEX UPD, 1 1-NODE)  = ',I6,
     1       /,5X,'ISFCVG (0 NL CONVG W/O SFEX, 1 W/ SFEX) = ',I6)
 1050 FORMAT(  5X,'NLRELX (0 NORELX,1 CONS RELX,2 VAR RELX)= ',I6)
 1055 FORMAT(  5X,'OMEGA  (RELAXATION PARAMETER, CONSTANT) = ',1PE15.5)
 1060 FORMAT(  5X,'L2NORM (0 INFINITY NORM, ELSE L2 NORM)  = ',I6,
     1       /,5X,'TOLUNS (TOLERANCE FOR NONLINEAR ITER)   = ',1PE15.5,
     2       /,5X,'TOLSWI (TOLERANCE FOR BC SWITCHING)     = ',1PE15.5,
     3       /,5X,'ERNLMX (MAX ALLOWABLE CVG OR RESID ERR) = ',1PE15.5)
 1070 FORMAT(/,5X,'CONSTANT INITIAL PRESSURE HEAD          = ',1PE15.5)
 1071 FORMAT(/,5X,'CONSTANT INITIAL PONDING HEAD           = ',1PE15.5)
 1072 FORMAT(/,1X,'INITIAL PONDING HEAD',/,(4(I6,2X,1PE11.3)))
 1090 FORMAT(/,1X,'INITIAL PRESSURE HEAD',/,(4(I6,2X,1PE11.3)))
 1091 FORMAT(/,1X,'INITIAL PONDING+PRESSURE HEAD',/,(4(I6,2X,1PE11.3)))
 1100 FORMAT(//,5X,
     1   'SATURATED HYDRAULIC CONDUCTIVITY, SPECIFIC STORAGE, AND ',
     2   'POROSITY VALUES',/,
     3   1X,' LAYER MAT.TYPE  X-PERM       Y-PERM       Z-PERM',
     4      '       STORAGE      POROSITY          VGN        VGRMC',
     5      '        VGPSAT')
 1110 FORMAT(1X,I4,I8,2X,8(1PE13.5))
C1120 FORMAT(/,5X,'NDIR (# OF NON-ATM, NON-SF DIR NODES 2D)= ',I6)  
C1150 FORMAT(/,5X,'NON-ATM, NON-SF DIR NODES IN 2-D MESH')
C1155 FORMAT(/,5X,'NATM, NSF FIXED DIR NODES IN 3-D MESH')
C1160 FORMAT(1X,10I7)
C1170 FORMAT(/,5X,'NDIRC(# OF FIXED NATM,NSF DIR NODES 3-D)= ',I6)  
C1175 FORMAT(/,5X,'NP   (TOTAL # OF NATM,NSF DIR NODES 3-D)= ',I6)
C1180 FORMAT(/,5X,'NQ   (# OF NON-ATM, NON-SF NEU NODES 3D)= ',I6)  
C1200 FORMAT(/,5X,'NON-ATM, NON-SF NEU NODES IN 3-D MESH')
 1215 FORMAT(/,5X,'PMIN  (AIR DRY PRESSURE HEAD VALUE)     = ',1PE15.5) 
 1216 FORMAT(  5X,'IPEAT (flag for peat deformation)       = ',I6)
 1217 FORMAT(  5X,'CBETA0                                  = ',1PE15.5,
     1       /,5X,'CANG                                    = ',1PE15.5)
 1220 FORMAT(  5X,'IVGHU (-1 moisture curve lookup table, ',
     1       /,5X,'       0 van Genuchten, ',
     2       /,5X,'       1 extended van Genuchten, ',
     3       /,5X,'       2 Huyakorn with Kr=Se**n, ',
     4       /,5X,'       3 Huyakorn with log_10 Kr(Se), ',
     5       /,5X,'       4 Brooks-Corey)                  = ',I6)
 1230 FORMAT(  5X,'VGN                                     = ',1PE15.5,
     1       /,5X,'VGRMC                                   = ',1PE15.5,
     2       /,5X,'VGPSAT                                  = ',1PE15.5)
 1235 FORMAT(  5X,'BCBETA                                  = ',1PE15.5,
     1       /,5X,'BCRMC                                   = ',1PE15.5,
     2       /,5X,'BCPSAT                                  = ',1PE15.5)
 1240 FORMAT(  5X,'HUALFA                                  = ',1PE15.5,
     1       /,5X,'HUBETA                                  = ',1PE15.5,
     2       /,5X,'HUGAMA                                  = ',1PE15.5,
     3       /,5X,'HUPSIA                                  = ',1PE15.5,
     4       /,5X,'HUSWR                                   = ',1PE15.5)
 1250 FORMAT(  5X,'HUN                                     = ',1PE15.5)
 1260 FORMAT(  5X,'HUA                                     = ',1PE15.5,
     1       /,5X,'HUB                                     = ',1PE15.5)
 1300 FORMAT(/,5X,'N     (# OF NODES IN 3-D MESH)          = ',I8,
     1       /,5X,'NT    (# OF TETRAHEDRA IN 3-D MESH)     = ',I8)
 1340 FORMAT(  5X,'LAYER ',I3,5X,'ZRATIO = ',1PE13.5)
 1345 FORMAT(//,' INPUT ERROR : ZRATIO VALUES MUST SUM TO 1.0')
 1350 FORMAT(/,5X,'IVERT  (TYPE OF VERTICAL DISCRETIZATION)= ',I6,
     1       /,5X,'ISP    (0 FLAT SURFACE, ELSE NOT FLAT)  = ',I6,
     2       /,5X,'BASE   (THICKNESS OR BASE OF 3-D MESH)  = ',1PE15.5)
 1360 FORMAT(/,5X,' SURFACE ELEVATION VALUES'/(4(I5,1X,1PE12.4)))
C1400 FORMAT(/,5X,'NSF   (# OF SEEPAGE FACES)              = ',I6)
C1410 FORMAT(  5X,'NUMBER OF NODES ON SEEPAGE FACE ',I6,'  = ',I6)
C1420 FORMAT(  5X,'NODE #''S : ',10I6)
C1500 FORMAT(//,' INPUT ERROR : ELEVATION VALUES NOT IN DESCENDING',
C    1          ' ORDER ON SEEPAGE FACE ',I6,
C    2       /,4X,'NODES',I4,' (NODE #',I6,') AND',I4,' (NODE #',I6,')')
C1600 format(I5,I5,I5/I5,I5,I5/I5,I5,I5) 
 1900 FORMAT(/,5X,'INCORRECT INPUT VALUE FOR ISIMGR: ',I6)
 1910 FORMAT(/,5X,'INVALID IOPT AND IPEAT COMBINATION (peat soil ',
     1            'deformation is not yet ',
     2       /,5X,'supported for the Newton scheme)')
 1930 FORMAT(/,5X,'INVALID KSLOPE AND IPEAT COMBINATION (peat soil ',
     1            'deformation is not yet ',
     2       /,5X,'supported for chord slope and tangent slope ',
     3            'differentiation of moisture curves)')
 1940 FORMAT(/,5X,'INVALID IVGHU AND IPEAT COMBINATION (peat soil ',
     1            'deformation is not yet ',
     2       /,5X,'supported for the extended van Genuchten or ',
     3            'Huyakorn moisture curves ',
     4       /,5X,'(unless they are given in lookup table form))')
 1960 FORMAT(/,5X,'NLKP must be at least 3')
 1970 FORMAT(/,5X,'PCAP values must be in ascending order')
      END
