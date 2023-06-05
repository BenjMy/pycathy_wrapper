C
C***************************  PARAMETER INCLUDE FILE CATHY.H ***********
C
C Dimensioning Parameters for CATHY Code
C --------------------------------------
C   ROWMAX - maximum NROW
C            NROW = # of rows in the DEM
C   COLMAX - maximum NCOL
C            NCOL = # of columns in the DEM
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
C   NP2MAX - maximum NDIR  
C            NDIR  = # of non-atmospheric, non-seepage face Dirichlet
C                    nodes in 2-d mesh
C   MAXSTR - maximum NSTR  
C            NSTR  = # of vertical layers 
C   NMAX   - NODMAX*(MAXSTR + 1)  (maximum N)
C            N     = NNOD*(NSTR + 1)
C                  = # of nodes in 3-d mesh
C   NTEMAX - 3*NTRMAX*MAXSTR  (maximum NT)   
C            NT    = 3*NTRI*NSTR
C                  = # of tetrahedra in 3-d mesh
C   NPMAX  - maximum NP    
C            NP    = NDIR*(NSTR + 1) + NDIRC 
C                  = total # of non-atmospheric, non-seepage face
C                    Dirichlet nodes in 3-d mesh
C            NDIRC = # of 'fixed' non-atmospheric, non-seepage face
C                    Dirichlet nodes in 3-d mesh
C   NPMAX_TRA - MAXIMUM NP-TRA
C             NP_TRA= NDIR_TRA*(NSTR+1) + NDIRC_TRA
C                   = total # of transport Dirichlet node in 3d mesh 
C             NDIRC_TRA = # of fixed transport Dirichlet nodes in 3d mesh 
C   NQMAX  - maximum NQ    
C            NQ    = # of non-atmospheric, non-seepage face Neumann 
C                    nodes in 3-d mesh
C   NSFMAX - maximum NSF
C            NSF   = # of seepage faces
C   NNSFMX - maximum # of nodes on a seepage face + 1
C   MAXDIR - NODMAX + NPMAX + NSFMAX*NNSFMX (maximum NUMDIR)
C            NUMDIR= total # of Dirichlet nodes in 3-d mesh
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
C   MAXTRM - maximum NTERM (value of MAXTRM should be at least 10*N)
C            NTERM = # of nonzero elements in system matrices
C   MAXIT  - maximum ITUNS 
C            ITUNS = maximum nonlinear FLOW3D iterations per time step
C   NRMAX  - maximum NR    
C            NR    = # of nodes selected for partial output
C   MAXPRT - maximum NPRT  
C            NPRT  = # of time values for detailed output
C   MAXVP  - maximum NUMVP
C            NUMVP = # of surface nodes for vertical profile output
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
C   DEMRES - resolution factor for grid generation from DEM
C   MAXQOUT- maximum NUM_QOUT
C            NUM_QOUT = # of surface cells for discharge output 
C   NFACEMAX  - maximum NFACE
C               NFACE  = # of faces
C***********************************************************************
C
      INTEGER   ROWMAX,COLMAX
      INTEGER   MAXCEL,MAXRES
      INTEGER   NODMAX,NTRMAX
      INTEGER   NP2MAX,MAXSTR
      INTEGER   NMAX,NTEMAX
      INTEGER   NPMAX,NQMAX,NSFMAX,NNSFMX
      INTEGER   MAXDIR,NPMAX_TRA
      INTEGER   MAXNUDN,MAXNUDT,MAXNUDC
      INTEGER   MAXZON,MAXTRM,MAXIT,MAXVEG
      INTEGER   NRMAX,MAXPRT,MAXVP
      INTEGER   N1MAX,NTPMAX
      INTEGER   MAXBOT,INTBOT
      INTEGER   DEMRES
      INTEGER   MAXQOUT
      INTEGER   NFACEMAX
      INTEGER   NIAUXMAX,NRAUXMAX,NQMAX_TRA
      INTEGER   MAXVTKPRT
      INTEGER   MAXFCONTONODE,MAXLKP
  
C
      PARAMETER (ROWMAX=20,COLMAX=20,DEMRES=1)
      PARAMETER (MAXCEL=ROWMAX*COLMAX,MAXRES=1)
      PARAMETER (NODMAX=441)
      PARAMETER (NTRMAX=800)
      PARAMETER (NP2MAX=1,MAXSTR=15)
      PARAMETER (NFACEMAX=74000)
      PARAMETER (NMAX=NODMAX*(MAXSTR + 1),NTEMAX=3*NTRMAX*MAXSTR)
      PARAMETER (NPMAX=1,NPMAX_TRA=1,NQMAX=1,NSFMAX=1)
      PARAMETER (NNSFMX=1,MAXDIR=NODMAX+NPMAX+NSFMAX*NNSFMX)
      PARAMETER (MAXNUDN=1,MAXNUDT=1,MAXNUDC=1)
      PARAMETER (MAXZON=1,MAXTRM=180000,MAXIT=30,MAXVEG=1)
      PARAMETER (NRMAX=0,MAXPRT=9,MAXVP=1)
      PARAMETER (N1MAX=25,NTPMAX=N1MAX*NMAX)
      PARAMETER (MAXBOT=1,INTBOT=1,MAXQOUT=1)
cxcx  PARAMETER (NIAUXMAX=NFACEMAX + MAXTRM + 1)
cxcx  PARAMETER (NRAUXMAX=5*NFACEMAX + MAXTRM,NQMAX_TRA=NODMAX)
      PARAMETER (NIAUXMAX=NMAX + MAXTRM + 1)
      PARAMETER (NRAUXMAX=5*NMAX + MAXTRM,NQMAX_TRA=NODMAX)
      PARAMETER (MAXVTKPRT=9)
      PARAMETER (MAXFCONTONODE=100,MAXLKP=3)
