C
C***************************  PARAMETER INCLUDE FILE CATHY.H ***********
C
C Dimensioning parameters for CATHY
C ---------------------------------
C   ROWMAX - maximum NROW
C            NROW = # of rows in the DEM
C   COLMAX - maximum NCOL
C            NCOL = # of columns in the DEM
C   DEMRES - resolution factor for grid generation from DEM
C   MAXCEL - ROWMAX*COLMAX (maximum NCELL)
C            NCELL = # of cells in the DEM of the catchment, including
C                    "lake" cells
C   MAXRES - maximum NUMRES
C            NUMRES = # of 'reservoirs' defined in the DEM
C   NODMAX - maximum NNOD
C            NNOD  = # of nodes in 2-d mesh
C                  = # of surface nodes in 3-d mesh
C   NTRMAX - maximum NTRI
C            NTRI  = 2*NCELL when SURF_ROUTE is active, otherwise must
C                    be assigned explicitly
C                  = # of triangles in 2-d mesh
C   MAXSTR - maximum NSTR
C            NSTR  = # of vertical layers
C   NMAX   - NODMAX*(MAXSTR + 1)  (maximum N)
C            N     = NNOD*(NSTR + 1)
C                  = # of nodes in 3-d mesh
C   NTEMAX - 3*NTRMAX*MAXSTR  (maximum NT)
C            NT    = 3*NTRI*NSTR
C                  = # of tetrahedra in 3-d mesh
C   MAXTRM - maximum NTERM (value of MAXTRM should be at least 10*N)
C            NTERM = # of nonzero elements in system matrices
C   NP2MAX - maximum NDIR
C            NDIR  = # of non-atmospheric, non-seepage face Dirichlet
C                    nodes in 2-d mesh
C   NPMAX  - maximum NP
C            NP    = NDIR*(NSTR + 1) + NDIRC
C                  = total # of non-atmospheric, non-seepage face
C                    Dirichlet nodes in 3-d mesh
C            NDIRC = # of 'fixed' non-atmospheric, non-seepage face
C                    Dirichlet nodes in 3-d mesh
C   NSFMAX - maximum NSF
C            NSF   = # of seepage faces
C   NNSFMX - maximum # of nodes on a seepage face + 1
C   MAXDIR - NODMAX + NPMAX + NSFMAX*NNSFMX (maximum NUMDIR)
C            NUMDIR= total # of Dirichlet nodes in 3-d mesh
C   NQMAX  - maximum NQ
C            NQ    = # of non-atmospheric, non-seepage face Neumann
C                    nodes in 3-d mesh
C   NPMAX_TRA - MAXIMUM ANP_TRA
C             ANP_TRA= NDIR_TRA*(NSTR+1) + NDIRC_TRA
C                    = total # of transport Dirichlet node in 3d mesh
C             NDIRC_TRA = # of fixed transport Dirichlet nodes in 3d mesh
C   MAXNUDN- maximum NUDN
C            NUDN  = # of observation points for nudging
C                    (NUDN=0 for no nudging)
C   MAXNUDT- maximum NUDT
C            NUDT  = # of observation times for nudging
C   MAXNUDC- maximum NUDC
C            NUDC  = # of concurrent observation datasets for nudging
C                    at any given time
C   MAXZON - maximum NZONE
C            NZONE = # of material types in the porous medium
C   MAXVEG - maximum NVEG
C            NVEG = # of vegetation types
C   MAXIT  - maximum ITUNS
C            ITUNS = maximum nonlinear FLOW3D iterations per time step
C   NRMAX  - maximum NR
C            NR    = # of nodes selected for partial output
C   MAXPRT - maximum NPRT
C            NPRT  = # of time values for detailed output
C   MAXVP  - maximum NUMVP
C            NUMVP = # of surface nodes for vertical profile output
C   MAXQOUT- maximum NUM_QOUT
C            NUM_QOUT = # of surface cells for discharge output
C   N1MAX  - maximum N1
C            N1    = maximum # of element connections to a node
C   NTPMAX - N1MAX*NMAX
C   MAXBOT - maximum IBOT (defined real working storage dimension for
C            NONSYM solver)
C            IBOT  = size of real working storage for NONSYM solver
C   INTBOT - MAXBOT + 6*NMAX + 1 (defined integer working storage
C            dimension for NONSYM solver)
C            Note: the values of MAXBOT and INTBOT should be
C            set to 1 when NONSYM is not used
C   NFACEMAX  - maximum NFACE
C               NFACE  = # of faces
C   NIAUXMAX - ...
C   NRAUXMAX - ...
C   MAXFCONTONODE - ...
C   MAXLKP - maximum NLKP (the value of MAXLKP must be at least 3)
C            NLKP = # of lookup values in the moisture curve tables for
C                   IVGHU = -1
C   MMM_MAX_NRC -   maximum MMM_NRC
C                   MMM_NRC = # of reactions in MMM scheme
C   MMM_MAX_NM  -   maximum MMM_NM
C                   MMM_NM = # of molecules in MMM scheme
C   MMM_MAX_NBG -   maximum MMM_NBG
C                   MMM_NBG = # of bacteria groups in MMM scheme
C   MMM_MAX_NCOMP - maximum MMM_NCOMP, set to MMM_MAX_NM + MMM_MAX_NBG
C                   MMM_NCOMP = # of components (or species) in MMM scheme,
C                               equal to MMM_NM + MMM_NBG
C
C***********************************************************************
C
      INTEGER   ROWMAX,COLMAX,DEMRES
      INTEGER   MAXCEL,MAXRES
      INTEGER   NODMAX,NTRMAX
      INTEGER   MAXSTR,NMAX,NTEMAX,MAXTRM
      INTEGER   NP2MAX,NPMAX,NSFMAX,NNSFMX,MAXDIR,NQMAX
      INTEGER   NPMAX_TRA
      INTEGER   MAXNUDN,MAXNUDT,MAXNUDC
      INTEGER   MAXZON,MAXVEG
      INTEGER   MAXIT
      INTEGER   NRMAX,MAXPRT,MAXVP,MAXQOUT
      INTEGER   N1MAX,NTPMAX
      INTEGER   MAXBOT,INTBOT
      INTEGER   NFACEMAX
      INTEGER   NIAUXMAX,NRAUXMAX
      INTEGER   MAXFCONTONODE
      INTEGER   MAXLKP
      INTEGER   MMM_MAX_NRC,MMM_MAX_NM,MMM_MAX_NBG,MMM_MAX_NCOMP
C
      PARAMETER (ROWMAX=75,COLMAX=10,DEMRES=1)
      PARAMETER (MAXCEL=ROWMAX * COLMAX,MAXRES=1)
      PARAMETER (NODMAX=(ROWMAX / DEMRES + 1)*(COLMAX / DEMRES + 1))
      PARAMETER (NTRMAX=2 * MAXCEL / (DEMRES * DEMRES))
      PARAMETER (MAXSTR=9)
      PARAMETER (NMAX=NODMAX * (MAXSTR + 1),NTEMAX=3 * NTRMAX * MAXSTR)
      PARAMETER (MAXTRM=65000)
      PARAMETER (NP2MAX=1,NPMAX=66,NSFMAX=11,NNSFMX=4)
      PARAMETER (MAXDIR=NODMAX + NPMAX + NSFMAX * NNSFMX)
      PARAMETER (NQMAX=1)
      PARAMETER (NPMAX_TRA=1)
      PARAMETER (MAXNUDN=1,MAXNUDT=3,MAXNUDC=1)
      PARAMETER (MAXZON=750,MAXVEG=1)
      PARAMETER (MAXIT=30)
      PARAMETER (NRMAX=1,MAXPRT=50,MAXVP=40,MAXQOUT=1)
      PARAMETER (N1MAX=25,NTPMAX=N1MAX * NMAX)
      PARAMETER (MAXBOT=1,INTBOT=1)
      PARAMETER (NFACEMAX=100000)
cxcx  PARAMETER (NIAUXMAX=NFACEMAX + MAXTRM + 1)
cxcx  PARAMETER (NRAUXMAX=5 * NFACEMAX + MAXTRM)
      PARAMETER (NIAUXMAX=NMAX + MAXTRM + 1)
      PARAMETER (NRAUXMAX=5 * NMAX + MAXTRM)
      PARAMETER (MAXFCONTONODE=100)
      PARAMETER (MAXLKP=3)
      PARAMETER (MMM_MAX_NRC=5,MMM_MAX_NM=8,MMM_MAX_NBG=3)
      PARAMETER (MMM_MAX_NCOMP=MMM_MAX_NM + MMM_MAX_NBG)
