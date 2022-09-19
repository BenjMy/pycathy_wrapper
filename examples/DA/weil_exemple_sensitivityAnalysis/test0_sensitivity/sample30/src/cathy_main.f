C
      PROGRAM CATHY_MAIN
C
C-----------------------------------------------------------------------
C
C CATchment HYdrologic model: coupled subsurface three-dimensional
C variably saturated flow (FLOW3D) and surface routing for hillslopes,
C pools/lakes, and channels (SURF_ROUTE)
C
C   CRS4, Cagliari, Italy
C   DMMMSA, Padova, Italy
C   IMAGE, Padova, Italy
C   INRS-ETE, Quebec, Canada
C   UNIMORE, Reggio Emilia, Italy
C   MIUB, Bonn, Germany
C   LHYGES, Strasbourg, France
C
C   Revision History:
C   ----------------
C   May/2018   : added variable vegetation type (variable PCANA,PCREF,
C                PCWLT,ZROOT,PZ,OMGC within the domain; they are all defined
C                in the SOIL input file). The root_map input file is now
C                VEG_TYPE, an integer raster map, where NVEG is now
C                an INTEGER. MC
C   Jan/18     : new standard version including transport w/ or w/o 
C                dispersion, coupling based on Dirichlet or Cauchy.
C                EnKF and PF removed from the code, as DA is now available
C                as an external shell written in R (anna.botto@unipd.it).
C   Feb/17     : added subroutine ETRAN for computing root water uptake,
C                according to the method implemented in cathy-noah. Feddes
C                parameters must be given in terms of pressure head.
C                WARNING! Not implemented in the EnKF/SIR. MC
C   Apr/15     : added hetereogeneity in VGN parameters (VGN,VGPSAT,VGRMC).
C                These parameters are assigned at the cells and calculated
C                then at the nodes. We read these parameters in the same way
C                of the hydraulic condutivity, porosity, and specific
C                storage. Note that the heterogeneity is not yet supported 
C                (i.e., parameters not perturbed) for the EnKF and SIR DA 
C                schemes. In addition it is implemented just 
C                for the cases IVGHU = 0,1. MS
C   Mar/15     : added hgraph output file in a netCDF format containing 
C                the discharge at the outlet cell and at the internal
C                selected cells. Note that internal grid cells
C                contain the INFLOW values. The file is written in the new
C                subroutine "detoutq_netcdf.f". MS
C   Mar/15     : added a flag in parm (NCOUT) to have a diagnostic
C                output file in a netCDF format containing PSI,
C                SW, PONDNOD variables. The file is written in the
C                new subroutine "detout_netcdf.f". MS
C   Apr/14     : modified Neumann BCs. Re-introduced NNEU for vertically
C                duplicated fluxes. If NNEU<0, free drainage BCs assigned
C                to all NNOD bottom nodes.
C   Jun/12     : added a moisture curve lookup table option for the
C                Newton scheme. Note that the lookup table option is not
C                yet supported for the EnKF and SIR DA schemes. The
C                existing lookup handling option for the Picard scheme
C                has been made more robust (foolproof), consistent with
C                the new implementation for the Newton scheme. Also,
C                mass balance errors are now computed for the lookup
C                table option (ie, the IVGHU=-1 case has been added to
C                subroutines CHMASS and CHMASPT). Finally, for
C                consistency with how heterogeneity is currently handled
C                for saturated conductivity, specific storage, and
C                porosity, the parameter IRETC has been removed, the
C                input of PCAP,SATC,KRWC has been changed from layer-
C                within-zone to zone-within-layer, and the structure of
C                the PCAP,SATC,KRWC arrays has been changed from
C                (NZONE,NSTR,NLKP) to (NSTR,NZONE,NLKP). A parameter
C                similar to IRETC can be re-introduced in a future
C                version to handle heterogeneity by layer only or by
C                zone only for both unsaturated and saturated media
C                parameters. CP
C   May/12     : parameter PMIN removed from common block CHORD in
C                SOILCHAR.H and put into new common block COMPMIN in
C                SOILCHAR.H. CP
C   Mar/12     : bug fix for moisture curve lookup table option.
C                SOILCHAR.H was missing the line
C                "COMMON  /TBLKUP/PCAP,SATC,KRWC". CP
C   May/12     : new subroutine BCREWIND for time dependent Dir and Neu BC 
C                with the EnKF/SIR. DP
C   Apr/12     : new variables ZONE_VERT and STR_ZON for the vertical 
c                stratification of the parameters, to differentiate
c                between parameter zonation and grid layers.
c                If a zone includes more then one layer of the grid,
c                then the parameters in this zone are updated/perturbed in
c                in the same way (used if DAFLAG>1). This new parameters
c                are read at the of the file soil. (Pasetto-Manoli) 
C   Oct/11     : corrected bug of atmospheric BC/variable node BC. 
C                Now the nodes that are either Dir or Neu are not
C                included in the ATM nodes, and they can change with 
C                time following the old BCONE/BCNXT mechanism. MP-MC
C   Aug/11     : introduced new measure for DA schemes: 
C                ENKFNOD(I,2) = 3 -->  measures of ERT (geolettriche)
C                (Pasetto and Manoli)
C   Apr/11     : modify the subroutine updatesir ; new variables NEFFMIN and
C                DSMEASMAX read in the input unit IIN40; each 
C                update has a different DSMEAS
C   Oct/10     : uncertainty of the surface parameters (ks and w) are now
C                taken into account. MC
C   Sep/10     : Implemented Evensen's original routine for EnKF, now the
C                only routine for EnKF (see DAFLAG). Keppenne's algorithm
C                no longer used. New option to include subsurface parameters
C                in the system state, for calibration purposes. MC
C   Dic/09     : New EnKF algorithm DP
C   Jun/09     : New particle filtering algorithm (SIR) MC-DP
C   Nov/08     : Added time correlation of atmospheric input observation
C                errors for the EnKF (see ATMTAU). MC
C   Sep/08     : Fixed a bug in routine atmone. Now if HSPATM=9999, no
C                atmospheric boundary conditions are considered, i.e., 
C                switching procedure is deactivated. MP
C   Jun/08     : Introduced a new flag for the choice of weighting functions 
C                in nudging (see WFLAG). MC
C   May/08     : Implemented Keppenne's algorithm for EnKF (see DAFLAG). MC
C   May/08     : Implemented recharge calculation. MC
C   Dec/07     : Merging of ENKF data assimilation module implemented 
C                by Matteo Camporese and surface routing modifications
C                developed by Mauro Sulis. MC-MS
C   Sep-Oct/07 : Introduction of multiple rivulets per hillslope cell.
C                The number of rivulets per hillslope cell is stored
C                into the dtm_nrc input file (unit IIN39).
C                The Muskingum-Cunge scheme is applied for each rivulet
C                and the total discharge flowing out from the cell is
C                now obtained by summing the values of the all rivulets of
C                the cell. MS-SO  
C   Aug/07     : Improved the EnKF algorithm, by the implementation of the
C                low rank square root by Evensen [Ocean Dynamics, 2004]. MC
C   Nov-Apr/07 (3): Added NUM_QOUT and ID_QOUT variables in parm input file
C                   (the outlet cell ID is automatically considered, 
C                   so only internal cells ID have to specified). MS
C   Nov-Apr/07 (2): Modification of the subroutine route.f. The surface 
C                   runoff propagation is simplified through the definition
C                   of an ordering system based on the descending elevation
C                   value of the cells stored in the qoi_a input file 
C                   (unit IIN23). In order to take into account the case of
C                   dispersion the surface runoff is propagated along the
C                   cardinal and diagonal directions through the definition
C                   of the relative weights assigned into the dtm_w_1 and
C                   dtm_w_2 input files (units IIN25 and IIN26). 
C                   The at-a-station exponents for the water surface width
C                   and ks coefficient are now variable in space and loaded
C                   through the dtm_b1_sf and dtm_y1_sf input files 
C                   (units IIN37 and IIN38). All the subroutines introduced by
C                   Anna-Chiara Bixio and related to the handling of reservoirs
C                   need to be double-checked. MS                    
C   Nov-Apr/07 (1): Definition and characterization of river network
C                   previously made by the network.f subroutine
C                   is now performed externally by the pre-processor
C                   developed by Stefano Orlandini and Giovanni Moretti. 
C                   The results of the pre-processing are loaded in datin.f
C                   subroutine by:
C                   - units IIN27-28 (drainage direction); 
C                   - units IIN29-30 (local slope);
C                   - units IIN31-32 (elemental path length);
C                   - units IIN33-34 (downstream variation of ks coefficient);
C                   - units IIN35-36 (downstream variation of water surface
C                     width)
C                   The calculation of drainage network and of upslope drainage
C                   areas is improved by the implementation of D8-LTD/LAD and 
C                   Dinf methods. The distinction between hillslope and channel
C                   cells is improved by implementing a variation of the 
C                   constant critical support area (Montgomery and Dietrich 
C                   1988,1989),and the normalized divergence technique (Howard,
C                   1994). The river network and watershed characteristics are
C                   now stored in the common block RIVERNETWORK.H  MS-SO-GM
C   Jan/07: Improved ensemble Kalman filter. NUDFLAG is restored, since
C           DAFLAG decides only whether to use EnKF or nudging. MC
C   Oct/06: Implemented new data assimilation scheme, the ensemble Kalman
C           filter (EnKF). Conversion of input file IIN50 (former nudging input
C           file) to more general data assimilation input file. The scheme
C           (nudging or EnKF) is chosen on the basis of the DAFLAG value
C           (former NUDFLAG). MC
C   Ott/06: Added ISIMGR=3 option for surface flow (SURF_ROUTE) 
C           simulation with DEM input. The local contribution rate to
C           surface runoff is directly read from effraininp (unit IIN22)
C           input file by the subroutines EFFONE.f and EFFNXT.f. MS
C   Aug/06: Added variable IETO (read in atmbc): if 0 the atmospheric
C           boundary condition inputs are interpolated linearly between
C           input times, otherwise are considered as a piecewise constant
C           function (ietograph). MC
C   Aug/06: Added IVERT=4 option. The base of the 3D grid is flat,
C           but the first NSTR-1 layers from the top are parallel
C           to the surface. MC
C   Aug/06: Slight modification (to be verified!!) in the subroutine
C           SWITCH, in order to fix a bug which did not allow the correct
C           functioning of the PMIN concept in case of evaporation and
C           exfiltration (case H). MC   
C   May/06: Added INDP=4 and subroutine ICVDWT for partially saturated
C           vertical hydrostatic equilibrium IC's; in this case WTPOSITION
C           is the depth of the water table respect to the surface. MC
C   Apr/06  : added input of independent impermeable basement. The 
C           thickness of the domain is read as a raster map in IIN60, if
C           IVERT=3. MC
C   May/05  : added moisture curve lookup table option (IVGHU=-1 and
C           input file IIN16) to handle variable (by zone and layer)
C           moisture curve parameters. For Picard scheme only! (also for
C           Newton scheme as of May/12). MP-TS-MC
C   May/05  : added input of variable zone, dem containing surface 
C           material zone information and reading of new 3D coordinates
C           after actual mesh generation. MP-TS-MC
C   Mar/05  : IOUT41 file added: hydrograph of the discharge at the outlet
C           cell of the basin. MC
C   Mar/05  : Implemented nudging for pressure head measurements, adding the
C           flag NUDFLAG; in order to be consistent with the Richards 
C           formulation, the difference between measured and modelled
C           pressure heads is multiplied by the overall storage coefficient.
C           MC
C   Feb/05  : Seepage face boundary conditions nodes are now variable in time:
C           added variables sfv, sfvnum, sfvnod, sfvtim, that are used to store
C           the previous values of NSF, NSFNUM, NSFNOD in case of back-stepping.
C           MC
C   Jan/05  : Modified seepage face boundary conditions; DUPUIT flag added:
C           if 0 the seepage face is treated exactly as before, if 1 the nodes
C           below the exit point are Dirichlet nodes with an hydrostatic 
C           pressure profile (the seepage face exit point is hypothesized to
C           be the actual water level of the stream, lake, etc.) MC
C   Jul/04  : Modified non atmopheric non seepage face boundary conditions;
C           now the nodes are time variable, but the interpolation is made by
C           a piecewise constant function as well as the atmospheric boundary
C           conditions. MC
C   Oct-Dec/03: Added swelling/shrinkage model to simulate the reversible
C           displacements of organic soils. New parameters IPEAT, CBETA0,
C           THETA0, CANG, INDE, and PORE.
C           Note that this peat soil deformation option is currently
C           implemented only for: Picard scheme (and not Newton);
C           IVGHU = -1, 0, or 4 soil moisture curve options;
C           KSLOPE = 0 (analytical differentiation of moisture curves);
C           and DAFLAG = 0 (not the EnKF and SIR DA schemes). MC
C   Jul-Sep/02: Added input parameter CUTRGT, ... (still need to input
C           and output this value in DATIN, + make several other changes
C           re coupling); fix BCONE bug for steady state case ... ; ...
C           MP-CP
C   Aug 1/02: Introduced water storage variables STORE0,STORE1,STORE2
C           and subroutine STORCAL to calculate these.
C           MP-CP
C   Jul 30/02: Expanded "mbeconv" (IOUT5) output, introducing several
C           new "cumulative" variables for mass balance and hydrograph
C           diagnostics - CAERAS (used to be ETOT), CERRAS, CDSTOR,
C           CVIN, and CVOUT. CDSTOR for a free drainage simulation, for
C           example, gives us an independent estimate (VTOT is the
C           traditional estimate) of the total volume of water drained.
C           Also added in IOUT2 ("risul") the cumulative (total)
C           mass balance error based on CERRAS in addition to the
C           one already there based on CAERAS. Imo the CERRAS-based
C           global mbe estimate is more correct!
C           CP
C   Feb/02: Added output file IOUT43
C           CP
C   Jul/01: Added INDP=3 and subroutine ICVHWT for partially saturated
C           vertical hydrostatic equilibrium IC's
C           Arno Hilberts & Claudio Paniconi, WUR & CRS4
C   Feb-Sep/01: Nudging algorithm for assimilation of soil moisture
C           observation data. For Picard scheme only!
C           Marino Marrocu, Claudio Paniconi, Mario Putti, CRS4 & DMMMSA
C   May/00: Back-stepping for coupled model
C           CP-MP-AB
C   Apr/00: Completely new version of SWITCH; new routines 
C           PONDUPD (replaces PONDVER) and ADRSTN.
C           CP-MP
C   Jan-Mar/00: added new variables: 
C           IPOND for type of initial conditions.
C           DOSTEP, NCELL_COARSE for coarsening of DEM-based 
C           triangulation.
C           Added different time-stepping for FLOW3D and SURF_ROUTE
C           AB-MP
C   Dec/99: Added setting of ATMACT=ATMPOT in the case IFATM=2
C           in routine SWITCH. 
C           Coded new version of HGRAPH.
C           CP-MP
C   Nov/99: Changes to SURF_ROUTE: this routine now transfers
C           information (OVFLNOD) from nodes to cells, performs
C           routine, and transfer the resulting ponding head 
C           to nodal ponding head PONDNOD.
C           Added comments on SURF_ROUTE variables.
C           AB-MP
C   Oct/99: Major changes to SWITCH, HGRAPH, ATMONE, ATMNXT,
C           and added coupling routine PONDVER.
C           Mario Putti & Claudio Paniconi, DMMMSA & CRS4
C   May-Sep/99: More work on coupling FLOW3D--SURF_ROUTE.
C           Anna Chiara Bixio, Claudio Paniconi, Mario Putti
C   Feb-May/99: Coupling of subsurface flow model (FLOW3D and CMSH/VSAT)
C           with surface routing model of Stefano Orlandini as modified
C           by Anna Chiara Bixio.
C           Mario Putti & Claudio Paniconi, DMMMSA & CRS4
C   Feb/99: Bug fix on declaration of CPUVEC - was dimensioned to 9 
C           but being initialized to 10.
C           Mario Putti & Claudio Paniconi, DMMMSA & CRS4
C   May/96: Added output of cumulative totals of: seepage face flow
C           volume VSFFLW (VSFTOT); non-atmospheric, non-seepage face
C           Dirichlet flow volume VNDIN + VNDOUT (VNDTOT);
C           non-atmospheric, non-seepage face Neumann flow volume 
C           VNNIN + VNNOUT (VNNTOT); and total net flow volume
C           VIN + VOUT (VTOT). VTOT = VSFTOT + VNDTOT + VNNTOT + 
C           (cumulative total of atmospheric Dirichlet and Neumann
C           flow volumes VADIN + VADOUT + VANIN + VANOUT).
C           Unit IOUT36 added for these cumulative flow volume outputs.
C           Changed time step value for signalling a steady state
C           problem from DELTAT >= 1.0e+10 to DELTAT >= 1.0e+15.
C           Claudio Paniconi, CRS4
C   Nov/95: Added Brooks-Corey moisture curves. 
C           Added INDP=2 option and subroutine ICVHE for calculation 
C           of fully saturated vertical hydrostatic equilibrium IC's.
C           Fixed two bugs in the mass balance calculations (see my 
C           email to Mario of 29/11/95).
C           Added output units IOUT30, IOUT31, and IOUT32 for detailed 
C           seepage face (IOUT30), non-atmospheric, non-seepage face 
C           Dirichlet (IOUT31), and non-atmospheric, non-seepage face 
C           Neumann (IOUT32) hydrograph output
C           Claudio Paniconi, CRS4
C   Sep/95: First efforts towards the merging of two existing 3-D 
C           finite element subsurface flow models, FLOW3D and 
C           CMSH/VSAT.
C           Claudio Paniconi, CRS4
C   Aug/95: Better formatting of atmospheric inputs obtained from
C           Luc Debruyckere's meteorological data, and merging of CMSH 
C           and VSAT codes.
C           Jonathan Deckmyn, CRS4
C   Jan/94: Update of FLOW3D to incorporate features from Putti/Paniconi
C           code FLOW2D. See Changes-Dante-Nov.93 in 
C           ~/fortran/cathy/old-versions.flow3d/sep.95.
C           Alberto Dante, DMMMSA
C   Oct/93: Time stepping changes in a local version of VSAT code for
C           a special case simulation. See ~/readme/valentijn-summary.
C           Valentijn Pauwels, CRS4
C   Apr/93: Bug fix in CMSH code. See update-log.apr29.93 in
C           ~/fortran/cathy/old-versions.cmshvsat/sep.95.
C           Claudio Paniconi, CRS4
C   Mar/92: Major revisions to FLOW3D code. See CN-Mar 5/92 notes 
C           and update-log.CN-mar5.92 in 
C           ~/fortran/cathy/old-versions.flow3d/jan.94.
C           Mario Putti & Claudio Paniconi, DMMMSA
C   Mar/92: Major revisions to CMSH/VSAT code. See CN-Mar 5/92 notes 
C           and update-log.CN-mar5.92 in 
C           ~/fortran/cathy/old-versions.cmshvsat/sep.95.
C           Claudio Paniconi, DMMMSA
C   Oct/91: Revisions to FLOW3D code. See files in
C           ~/fortran/cathy/old-versions.flow3d/oct.91.
C           Mario Putti & Claudio Paniconi, DMMMSA
C   Mar/91: First version of FLOW3D code, then called LCM (lean, clean,
C           and mean), based on Gambolati/Pini/Putti codes FLOW3D, 
C           FLOWSUB, and FLOWCH.
C           Mario Putti & Claudio Paniconi, DMMMSA
C   Nov/86-Mar/92: Many revisions, major and minor, to CMSH/VSAT code. 
C           See CN-... notes from Oct 22/90, Jul 16/91, Feb 7/92, 
C           Mar 3/92, etc.
C           Claudio Paniconi, Princeton Univ & DMMMSA
C   Nov/86: First version of CMSH/VSAT code, based on Andrew Binley's 
C           VSAT3D code.
C           Claudio Paniconi, Princeton University
C
C--------------------------------DESCRIPTION OF PARAMETERS AND VARIABLES
C
C Dimensioning Parameters   (defined in parameter include file CATHY.H)
C -----------------------
C   ROWMAX - maximum NROW
C   COLMAX - maximum NCOL
C   MAXCEL - ROWMAX*COLMAX (maximum NCELL)
C   MAXRES - maximum NUMRES
C   NODMAX - maximum NNOD
C   NTRMAX - maximum NTRI
C   NP2MAX - maximum NDIR
C   MAXSTR - maximum NSTR
C   NMAX   - NODMAX*(MAXSTR + 1)  (maximum N)
C   NTEMAX - 3*NTRMAX*MAXSTR  (maximum NT)
C   NPMAX  - maximum NP
C   NQMAX  - maximum NQ
C   NSFMAX - maximum NSF
C   NNSFMX - maximum # of nodes on a seepage face + 1
C   MAXDIR - NODMAX + NPMAX + NSFMAX*NNSFMX (maximum NUMDIR)
C   MAXNUDN- maximum NUDN or NOBS
C   MAXNUDT- maximum NUDT 
C   MAXNUDC- maximum NUDC
C   MAXZON - maximum NZONE
C   MAXLKP - maximum NLKP (value of MAXLKP must be at least 3)
C   MAXTRM - maximum NTERM (value of MAXTRM should be at least 10*N)
C   MAXIT  - maximum ITUNS
C   NRMAX  - maximum NR
C   MAXPRT - maximum NPRT
C   MAXVP  - maximum NUMVP
C   MAXQOUT- maximum NUM_QOUT
C   N1MAX  - maximum N1
C   NTPMAX - N1MAX*NMAX
C   MAXBOT - maximum IBOT (defined real working storage dimension for
C            NONSYM solver)
C   INTBOT - MAXBOT + 6*NMAX + 1 (defined integer working storage 
C            dimension for NONSYM solver)
C            Note: the values of MAXBOT and INTBOT should be 
C            set to 1 when NONSYM is not used
C  
C Actual or Minimum Dimensions
C ----------------------------
C   NROW   - # of rows in the DEM
C   NCOL   - # of columns in the DEM
C   NCELL  - # of cells in the DEM of the catchment, including
C            "lake" cells;
C            set to 0 for the case DEM = FALSE
C   NUMRES - # of 'reservoirs' defined in the DEM
C   NNOD   - # of nodes in 2-d mesh. These are the surface nodes for 
C            the 3-d mesh -- they are all designated as atmospheric
C            boundary condition nodes (rainfall and evaporation inputs),
C            except for those surface nodes which are specifically 
C            designated as non-atmospheric BC's (see description of
C            NDIR, NDIRC, NQ, NSF, and IFATM).
C   NTRI   - 2*NCELL (when SURF_ROUTE is active, otherwise must
C            be assigned explicitly) = # of triangles in 2-d mesh
C   NDIR(3)- # of non-atmospheric, non-seepage face Dirichlet nodes 
C            in 2-d mesh. The BC's assigned to these surface nodes
C            are replicated vertically (compare NDIRC). Same for NNEU,NNEUC
C   NSTR   - # of vertical layers
C   N      - NNOD*(NSTR + 1) = # of nodes in 3-d mesh
C   NT     - 3*NTRI*NSTR = # of tetrahedra in 3-d mesh
C   NDIRC(3) - # of 'fixed' non-atmospheric, non-seepage face Dirichlet
C              nodes in 3-d mesh ('fixed' in the sense that these BC's
C              are not replicated to other nodes - compare NDIR)
C   NP(3)  - NDIR*(NSTR + 1) + NDIRC = total # of non-atmospheric,
C            non-seepage face Dirichlet nodes in 3-d mesh
C   NQ(3)  - NNEU*(NSTR ) + NNEUC = # of non-atmospheric, 
C            non-seepage face Neumann nodes in 3-d mesh
C   NSF    - # of seepage faces (see description of NSFNOD)
C   NUMDIR - total # of Dirichlet nodes in 3-d mesh
C   NUDN   - # of observation points for nudging
C            (NUDN=0 for no nudging) 
C   NUDT   - # of observation times for nudging
C   NUDC   - # of concurrent observation datasets for nudging at any
C            given time. NUDC is updated at the end of every time step.
C            If NUDC=0, nudging is inactive at the current time.
C   NZONE  - # of material types in the porous medium
C ZONE_VERT- # of vertical zones for the soil parameters
C            (used if DAFLAG > 0)
C   NLKP -   # of lookup values in the moisture curve tables for
C            IVGHU = -1. NLKP values of PCAP, SATC, and KRWC will be
C            read in for each zone (inner input loop) and for each layer
C            (outer input loop). NLKP must be at least 3, and for each
C            set of NLKP values, the PCAP values must be in ascending
C            order (e.g., -10.0, -9.9, -9.8, ..., -0.04, -0.02, 0.0).
C   NTERM  - # of nonzero elements in system matrices (symmetric storage
C            used for Picard scheme; nonsymmetric storage for Newton)
C   ITUNS  - maximum nonlinear FLOW3D iterations per time step
C   NR     - # of nodes selected for partial output
C   NPRT   - # of time values for detailed nodal output and element
C            velocity output (see description of IPRT, TIMPRT)
C   NUMVP  - # of surface nodes for vertical profile output
C   NUM_QOUT - # of selected cells for hydrograph output
C   N1     - maximum # of element connections to a node
C   IBOT   - size of real working storage for NONSYM solver (# of
C            nonzero elements in the LU decomposition of system matrix)
C
C General and FLOW3D Integer Parameters, Flags, and Indices
C ---------------------------------------------------------
C   KPRT   - index to current time value for detailed output
C   IPEAT  - flag for peat soil deformation
C            =0 constant porosity (in unsaturated soil)
C            =1 consider porosity variations with water saturation
C   IPRT   - flag for detailed output at all nodes and velocity and
C            water saturation output at all elements (velocity and
C            water saturation output in the case IPRT=4 can be used as
C            input to TRAN3D and DUAL3D codes)
C            =0 don't print nodal pressure, velocity, water saturation,
C               or relative conductivity values
C            =1 print only nodal pressure head values
C            =2 print nodal pressure head and velocity values
C            =3 print nodal pressure, velocity, and relative
C               conductivity values
C            =4 print nodal pressure, velocity, relative conductivity,
C               and overall storage coefficient values, and print
C               element velocity and nodal water saturation values
C   VTKF   - flag for the output in the vtk files at the detailed output
C            times:
C            = 0 do not print the vtk file
C            = 1 prints the grid and nodal pressure head
C            = 2 prints the grid and nodal pressure head and saturation
C            = 3 prints the grid and nodal pressure head and saturation 
C                and the hydraulic conductivity of the elements
C            = 4 prints also the element velocity
C   IPRT1  - flag for output of input and coordinate data 
C            in subroutines DATIN and GEN3D
C            =-1 reads in coordinates from file IIN51 just after   
C                grid generation (in subroutine grdsys)
C            =0 prints parameters only (default)
C            =1 prints parameters + b.c. + geom. char.
C            =2 prints parameters + b.c. + geom. char. + grid info
C            =3 prints parameters + b.c. + geom. char. + grid info,
C               X, Y, Z coordinate values in subroutine GEN3D, and then
C               terminates program execution
C   ISIMGR - flag for type of simulation and type of surface grid
C            =0 subsurface flow only (FLOW3D) with general triangular
C               grid input
C            =1 subsurface flow only (FLOW3D) with DEM input and
C               triangular grid generated from this DEM
C            =2 coupled subsurface flow (FLOW3D) and surface routing
C               (SURF_ROUTE) with DEM input and triangular grid
C               generated from this DEM
C            =3 surface flow only (SURF_ROUTE) with DEM input 
C   IVERT  - =0 each layer will be parallel to the surface, including
C               the base of the 3-d grid. ZRATIO is applied to each
C               vertical cross section.
C            =1 base of the 3-d grid will be flat, and ZRATIO is applied
C               to each vertical cross section
C            =2 base of the 3-d grid will be flat, as will the NSTR-1
C               horizontal cross sections above it. ZRATIO is applied
C               only to the vertical cross section having the lowest
C               elevation.
C            =3 for each cell of the dem a single depth value is read in file
C               input IIN60 (basement). ZRATIO is applied to each vertical
C               cross section.
C            =4 the first NSTR-1 layers from the surface will be parallel
C               to the surface and the base of the 3-d grid will be flat.
C               ZRATIO is applied only to the vertical cross section 
C               having the lowest elevation.
C   ISP    - =0 for flat surface layer (only one Z value is read in, and
C            is replicated to all surface nodes); otherwise surface
C            layer is not flat (Z values read in for each surface node)
C            (for ISP=0, IVERT=0, 1, and 2 yield the same 3-d mesh, 
C            given the same values of BASE and ZRATIO)  - MOREOVER -
C            >=2 read NNOD different values of VEG_TYPE, otherwise one
C                constant value of rooting depth
C   INDP   - flag for pressure head initial conditions (all nodes)
C            =0 for input of uniform initial conditions (one value 
C               read in)
C            =1 for input of nonuniform IC's (one value read in 
C               for each node)
C            =2 for calculation of fully saturated vertical hydrostatic
C               equilibrium IC's (calculated in subroutine ICVHE). In
C               the case of IPOND>0, the fully saturated hydrostatic IC
C               is calculated (in subroutine ICVHEPOND) starting from
C               the ponding head values at the surface nodes, rather
C               than surface pressure heads of 0.
C            =3 for calculation of partially saturated vertical
C               hydrostatic equilibrium IC's (calculated in subroutine
C               ICVHWT) with the water table height (relative to the
C               base of the 3-d grid) given by parameter WTPOSITION
C            =4 for calculation of partially saturated vertical
C               hydrostatic equilibrium IC's (calculated in subroutine
C               ICVDWT) with the water table depth (relative to the
C               surface of the 3-d grid) given by parameter WTPOSITION
C   ANP,   - actual values of NP and NQ that are passed to the 
C   ANQ      FLOW3D module
C   HTIDIR - =0 for temporally variable non-atmospheric, non-seepage
C            face Dirichlet boundary conditions inputs; otherwise
C            non-atmospheric, non-seepage face Dirichlet boundary
C            conditions inputs are homogeneous in time (see also
C            notes following description of QPOLD)
C   HTINEU - =0 for temporally variable non-atmospheric, non-seepage
C            face Neumann boundary conditions inputs; otherwise
C            non-atmospheric, non-seepage face Neumann boundary
C            conditions inputs are homogeneous in time (see also
C            notes following description of QPOLD)
C   HSPATM - =0 for spatially variable atmospheric boundary condition
C            inputs; blank or =9999 if unit IIN6 input is to be ignored;
C            otherwise atmospheric BC's are homogeneous in space
C   HTIATM - =0 for temporally variable atmospheric boundary condition
C            inputs; otherwise atmospheric BC's are homogeneous in time
C            (see also notes following description of ATMINP)
C   IETO   - =0 for linear interpolation of the atmospheric boundary condition
C            inputs between different ATMTIM; otherwise the inputs are assigned
C            as a piecewise constant function (ietograph)
C   LUMP   - =0 for distributed mass matrix; otherwise matrix is lumped
C   LUMPC  - =0 for distributed mass matrix (dispersive part of transport);
C             otherwise matrix is lumped for transport
C   IOPT   - =1 for Picard iteration scheme
C            =2 for Newton iteration scheme
C   NLRELX - flag for nonlinear relaxation
C            =0 no relaxation
C            =1 relaxation with constant relaxation parameter OMEGA
C            =2 relaxation with iteration-dependent relaxation parameter
C               OMEGA, calculated using Huyakorn et al's adaptation
C               (WRR 1986 22(13), pg 1795) of Cooley's empirical scheme
C               (WRR 1983 19(5), pg 1274)
C   L2NORM - =0 to use the infinity norm in the test for convergence of
C            the nonlinear FLOW3D and coupled FLOW3D/SURF_ROUTE
C            iterations; otherwise the L2 norm is used
C   ISOLV  - flag for nonsymmetric linear solver
C            =-5   BiCGSTAB (preconditioned with D^-1)
C            =-4   BiCGSTAB (not preconditioned)
C            =-3   TFQMR (preconditioned with D^-1)
C            =-2   TFQMR (not preconditioned)
C            =-1   TFQMR (preconditioned with K^-1)
C            =0    BiCGSTAB (preconditioned with K^-1)
C            =1    GRAMRB  (minimum residual)
C            =2    GCRK(5) (ORTHOMIN)
C            =3    IBM's NONSYM (direct solver)
C   IVGHU  - =-1 for moisture curve lookup table
C            =0 for van Genuchten moisture curves
C            =1 for extended van Genuchten moisture curves
C            =2 for moisture curves from Huyakorn et al (WRR 20(8) 1984,
C               WRR 22(13) 1986) with Kr=Se**n conductivity relationship
C            =3 for moisture curves from Huyakorn et al (WRR 20(8) 1984,
C               WRR 22(13) 1986) with conductivity relationship from
C               Table 3 of 1984 paper (log_10 Kr(Se) curve)
C            =4 for Brooks-Corey moisture curves
C   KSLOPE - =0 for analytical differentiation of moisture curves
C            =1 for "chord slope" and analytical differentiation
C            =2 for "chord slope" and centered difference formulas
C            =3 for localized "chord slope" and analytical
C               differentiation
C            =4 for localized "tangent slope" differentiation
C            (the "chord slope" formula is the tangent approximation
C            suggested by Huyakorn et al (WRR 20(8) 1984), wherein
C            derivatives are approximated using pressure heads at
C            the current and previous nonlinear iterations; "tangent
C            slope" differentiation is a different tangent approximation
C            wherein derivatives are approximated using pressure heads
C            at the endpoints of a given range (eg: endpoints PKRL, PKRR
C            for the derivative of relative hydraulic conductivity). For
C            KSLOPE=1,2 the chord slope formula is used at every
C            iteration and at all nodes (with some exceptions as
C            dictated by TOLKSL). For KSLOPE=3 or 4 the chord or tangent
C            slope formulas are used only at those nodes whose pressure
C            heads fall within given ranges (see PKRL, PKRR, etc), hence
C            'localized'; for nodes whose pressure heads fall outside
C            these ranges, analytical differentiation is used.)
C   ISFONE - =0 seepage face exit point updating performed by 
C               checking all nodes on a seepage face
C            =1 seepage face exit point updating performed by 
C               checking only the one node above and one node
C               below the current exit point position (deactivated)
C   ISFCVG - =0 convergence of seepage face exit points is not a 
C               condition for convergence of the nonlinear iterative
C               procedure
C            =1 convergence of seepage face exit points is a condition
C               for convergence of the nonlinear iterative procedure
C   DUPUIT - =0 all the nodes below the seepage face exit point are at 
C               atmospheric pressure
C            =1 all the nodes below the seepage face exit point are at
C               hydrostatic pressure. CAUTION: TO USE ONLY IF ISFONE=1!
C   KSFCV  - total number of seepage face exit point convergence failure
C            occurrences (over all nonlinear iterations and all time 
C            steps)
C   KSFCVT - total number of seepage face exit point convergence 
C            failures (over all seepage faces, all nonlinear iterations,
C            and all time steps)
C   ITER   - iteration index for nonlinear FLOW3D iterations for each
C            time step
C   MAXITER- max iteration index for nonlinear FLOW3D iterations for each
C            time step for the ensemble of realizations of EnKF
C   ITUNS1 - if ITER < ITUNS1, time step size is increased
C   ITUNS2 - if ITUNS1 <= ITER < ITUNS2, time step size is not altered
C            if ITUNS2 <= ITER < ITUNS, time step size is decreased
C            if ITER = ITUNS (ie, convergence not achieved in ITUNS
C            iterations), we back-step unless time step size cannot
C            be reduced any further (DELTAT = DTMIN). Back-stepping is
C            also triggered if the linear solver failed (LSFAIL = TRUE)
C            or if the convergence or residual errors become larger
C            than ERNLMX (ERRGMX = TRUE).
C   NSTEP  - time step index 
C   KBACKT - number of back-stepping occurrences at current time level
C   KBACK  - total number of back-stepping occurrences over all time
C            steps
C   ITRTOT - total number of nonlinear FLOW3D iterations over all time
C            steps
C   ITMXCG - maximum # of iterations for conjugate gradient linear
C            system solvers
C   NITER  - number of iterations for the linear solver at each 
C            nonlinear FLOW3D iteration
C   NITERT - number of iterations for the linear solver at each 
C            time step
C   ITLIN  - total number of iterations for the linear solver over all
C            nonlinear FLOW3D iterations and all time steps
C   KLSFAI - total number of linear solver failures
C   IMAX   - largest integer number (machine dependent)
C   MINBOT - minimum IBOT (value returned from the solver)
C   NDZ    - # of zero elements on the diagonal of the system matrices
C            (signals an error condition)
C   NUDFLAG- flag for the choice of the variable to be assimilated by
C            nudging:
C          = 0 soil moisture;
C          = 1 pressure head (in this case NUDG must be properly scaled)
C   WFLAG  - flag for the choice of the weighting functions for nudging:
C          = 0 Cressman-type functions [Paniconi et al., AWR, 2003];
C          = 1 Exponential and Gaussian correlation functions for time
C              and spatial influence, respectively.
C   TRAFLAG- flag for transport modeling
C          = 0 without transport
C          = 1 with transport
C
C  NADV   - integer number of advective time steps
C                   DeltDatadv deltat/nadv
C  NADVFL = 1  : the number of advective time steps is given in input by nadv
C        = 0  :  the number of advective time steps is calculated in the 
C        cflpecnumber subroutine considering the cflnum and the time step size
C  flag_temp = 0    --->  Euler scheme in time and Godunov reconstruction
C  flag_temp =1     --->  Euler scheme in time and 2nd order reconstruction
C  flag_temp=2      ---> correction term added in midpoint(Euler) scheme in time
C  flag_temp=3      ---> Mid-point scheme in time, Runge scheme
C  flag_temp=4     --->  Heun scheme in time
C flag_interp.eq.1 --->  construction of the interpolant as in Durlofsky
C                        scheme taking into account the 4 possible interpolant
C                        that can be constructed from the refering tetrahedron
C    flag_limiter = 0        Durloflsky limiter
C    flag_limiter = 1        min-limiter
C    flag_limiter = 2        the maximum  gradient interpolant is chosen with
C                            no limiter
C    flag_limiter =3         the average of the 4 interpolant is applied with
C                            no limiter
C    flag_interp .eq.2  -->  construction of the interpolant considering the 4
C                         neighbours tetrahedra
C    flag_limiter =1         extremum limiter (see paper: compressible large
C                            eddy simulation using unstructured grid:
C                            supersonic boundary layer and compression ramps,
C                            Yan, Urbin, Knight, 10th Intern. Conf. on Methods
C                            of Aerophysical Research,July 9-16 2000,
C                            Novosibirsk, Russia
C    flag_limiter =2         Limited Central Difference LCD limiter  (see paper
C                            Multidimensional Slope Limiters for MUSCL-type
C                            finite volume schemes on unstructured grids,
C                             M.E. Hubbard, J. Comp. Phys. 155, pp.54-74, 1999)
C    flag_limiter =3         no limiter is applied
C    flag_limiter=4          Barth-Jespersen limiter (see paper: The design and
C                            application of upwind schemes on unstructerd
C                            meshes, T. J. Barth, D. C. Jespersen, AIAA-89-0366)
C   flag_interp.eq.3  --->  MLG scheme: the interpolant is chosen between
C                         all the interpolant that can be constructed both
C                        with flag_interp.eq.1. and flag_interp.eq.2
C                        Between them, the 'proper' interpolant is limited
C                        with the LCD limiter (see above)
C   flag_interp.eq.4  --->  MLS interpolant: the interpolant is chosen by the
C                         least squares method
C     flag_limiter =1          extremum limiter
C     flag_limiter=2           LCD limiter
C     flag_limiter =3          no limiter
C     flag_limiter=4           Barth-Jespersen limiter
C
C   NFACE  - # of faces
C Integer Parameters for SURF_ROUTE
C ---------------------------------
C   NSURF        - # of time steps for the surface routing module
C                  (estimated to satisfy the Courant criterion) at
C                  the current time level of the subsurface module
C   NSURFT       - # of time steps for the surface routing module at
C                  the current subsurface time level including those
C                  from any back-steps at the current time level
C   NSURFT_T     - total # of time steps for the surface routing module
C                  over all subsurface time steps, excluding back-steps
C   NSURFT_TB    - total # of time steps for the surface routing module
C                  over all subsurface time steps, including back-steps
C   NCELNL       - # of cells in the DEM of the catchment excluding the
C                  "lake" cells
C   IPOND        - flag for ponding head initial conditions (surface
C                  nodes)
C                  =0 no input of ponding head initial conditions;
C                     otherwise (IPOND = 1 or 2) ponding head initial
C                     conditions are read into PONDNOD, and, where
C                     PONDNOD > 0, these values are used to update
C                     the surface node values in PTIMEP read in
C                     according to the previous INDP flag
C                  =1 uniform ponding head initial conditions (one value
C                     read in)
C                  =2 nonuniform ponding head initial conditions (one
C                     value read in for each node)
C   DOSTEP       -  step adopted in coarsening the mesh
C   NCELL_COARSE -  # of cells in the coarse mesh
C   HSPEFF       -  =0 for spatially variable effective rainfall inputs,
C                    otherwise effective rainfall is homogeneous in space.
C   HTIEFF       -  =0 for temporally variable effective rainfall inputs; 
C                    otherwise effective rainfall is homogeneous in time
C                    (see also notes following description of ATMINP)
C
C Integer Parameters for Data Assimilation
C ----------------------------------------
C   NUDCTR - counter for most recent nudging observation time. NUDCTR is
C            incremented by one whenever the next value of TIME, which
C            is updated at the end of the current time level, falls
C            within the temporal influence window of the next
C            observation time in NUDTIM.
C   NOBS   - number of observations (measurements) 
C
C Logical Flags
C -------------
C   FL3D   - flag indicating whether the simulation involves the
C            subsurface component or not (FLOW3D)
C             FALSE if no FLOW3D call (SURF_ROUTE simulation only)
C             TRUE if FLOW3D is to be called (coupled or not)
C   SURF   - flag indicating whether the simulation involves the
C            surface routing component or not (SURF_ROUTE)
C             FALSE if no SURF_ROUTE call (FLOW3D simulation only)
C             TRUE  if SURF_ROUTE is to be called (coupled or not)
C   DEM    - flag indicating whether a DEM is input or not
C             FALSE if no DEM input
C             TRUE  if DEM input
C   GRID   - flag indicating whether the triangular surface grid is
C            input or generated
C             FALSE if the grid is generated from DEM input
C             TRUE  if the grid is read in as input
C   PONDING- flag indicating whether there is ponding at the
C            surface and SURF_ROUTE has to be called
C            FALSE no ponding
C            TRUE  ponding -> call SURF_ROUTE
C   PONDP  - PONDING value at previous time level
C   DTGMIN - flag indicating whether the current time step size
C            is greater than the minimum allowed
C             FALSE if not greater
C             TRUE  if greater
C   LSFAIL - flag for linear solver
C             FALSE if linear solver did not fail
C             TRUE  if linear solver failed
C   ERRGMX - flag indicating whether the FLOW3D convergence or residual
C            errors have become greater than the allowed maximum
C             FALSE if not greater
C             TRUE  if greater
C   NORMCV - flag for convergence of the norm of pressure head 
C            differences in the nonlinear FLOW3D iterative procedure
C             FALSE if the norm has not converged
C             TRUE  if the norm has converged
C   ITAGEN - flag indicating whether we can iterate again in the
C            nonlinear FLOW3D iterative procedure
C             FALSE if we cannot iterate again
C             TRUE  if we can iterate again
C   SFCHEK - flag indicating whether it is necessary to check for
C            seepage face exit point convergence as a condition for
C            convergence of the nonlinear FLOW3D iterative procedure
C             FALSE if it is not necessary to check
C             TRUE  if it is necessary to check
C   KSFZER - flag for number of seepage face exit points which did
C            not converge at each nonlinear FLOW3D iteration
C             FALSE if one or more exit points did not converge
C             TRUE  if all exit points converged 
C   NOBACK - flag for back-stepping
C             FALSE if back-stepping still possible
C             TRUE  if no back-stepping possible
C   BCKSTP - flag indicating whether we are back-stepping
C             FALSE we are not in the back-stepping case
C             TRUE  we are doing back-step
C   NCOUT  - flag to dump PSI, SW, PONDNOD output variables in a netCDF file
C            (NCOUT=1 -> "cathy_netcdf_out.nc"). The same flag is used
C            for the hydrograph (NCOUT=1 -> "cathy_netcdf_qout.nc").
C   TRANSP - logical flag for transport modeling   
C             FALSE = no transport
C             TRUE  = transport 
C Real Scalars for Initial Conditions and Vertical Discretization 
C ---------------------------------------------------------------
C   WTPOSITION-for the case INDP=3, specifies the initial water table
C              height relative to the base of the 3-d grid;
C              for the case INDP=4, specifies the initial water table
C              depth relative to the surface
C   BASE     - value which defines the thickness or base of the 3-d
C              mesh. For IVERT=0, BASE is subtracted from each surface
C              elevation value, so that each vertical cross section will
C              be of thickness BASE, and the base of the 3-d mesh will
C              be parallel to the surface. For IVERT=1 or 2, BASE is 
C              subtracted from the lowest surface elevation value, say
C              ZMIN, so that each vertical cross section will be of 
C              thickness (Z - ZMIN) + BASE, where Z is the surface 
C              elevation for that cross section. The base of the 3-d
C              mesh will thus be flat.
C
C Real Scalars for Time Stepping and Linear and Nonlinear Iterations
C ------------------------------------------------------------------
C   TETAF  - weighting parameter for FLOW3D time stepping scheme 
C            (1.0 backward Euler; 0.5 Crank-Nicolson;
C            TETAF is set to 1.0 for steady state problem)
C   DELTAT - initial and current FLOW3D time step size (DELTAT >=
C            1.0e+15 on input indicates steady state problem)
C   DELTAT0- initial FLOW3D time step size, stored in case of EnKF
C            for reinitializing DELTAT when an update occurs
C   DELTATS - initial and current SURF_ROUTE time step size 
C             (this time step is kept constant during the simulation,
C             and the value equals to the DELTAT assigned in the parm
C             input file, IIN1)
C   DTMIN  - minimum FLOW3D time step size allowed
C   DTMAX  - maximum FLOW3D time step size allowed
C   DTAVG  - average FLOW3D time step size used for the simulation
C   DTSMAL - smallest FLOW3D time step size used during the simulation
C   DTBIG  - largest FLOW3D time step size used during the simulation
C   TSMAL  - first time at which DTSMAL is used
C   TBIG   - first time at which DTBIG is used
C   DTMAGA - magnification factor for FLOW3D time step size (additive)
C   DTMAGM - magnification factor for FLOW3D time step size
C            (multiplicative)
C   DTREDS - reduction factor for FLOW3D time step size (subtractive)
C   DTREDM - reduction factor for FLOW3D time step size (multiplicative)
C   TMAX   - time at end of simulation (TMAX is set to 0.0 for 
C            steady state problem)
C   TIME   - time at current time level
C   TIMEP  - time at previous time level
C   OMEGA  - nonlinear relaxation parameter: OMEGA > 1, over-relaxation;
C            OMEGA < 1, under-relaxation. Input value of OMEGA is used
C            only for the case NLRELX=1 (constant relaxation parameter).
C            Input value of OMEGA is ignored otherwise: for NLRELX=0
C            relaxation is not applied; for NLRELX=2 OMEGA is calculated
C            at each nonlinear FLOW3D iteration.
C   OMEGAP - OMEGA value at previous nonlinear FLOW3D iteration
C   TOLUNS - tolerance for convergence of nonlinear FLOW3D iterations
C   TOLSWI - tolerance for boundary condition switching check in FLOW3D
C            iterations (switching check is only performed when PINF or
C            PL2 are smaller than TOLSWI; so for eg if TOLSWI = TOLUNS,
C            switching check is only performed after convergence and not
C            after each iteration)
C   ERNLMX - maximum allowable convergence or residual error in the
C            nonlinear FLOW3D solution. If the convergence or residual 
C            errors become larger than ERNLMX, ERRGMX is set to TRUE and 
C            the code back-steps. This avoids occurrences of overflow or
C            underflow when nonlinear iterations diverge.
C   TOLCG  - tolerance for convergence of conjugate gradient linear
C            system solvers
C 
C Real Scalars for Ponding
C --------------------------------------------------
C   PONDH_MIN - minimum ponding head: if PNEW > PONDH_MIN, then at that
C               node there is ponding; otherwise there is no ponding
C
C Real Scalars for Nudging and EnKF
C --------------------------------------------------
C   NUDG   - nudging factor "G" which determines the relative strength
C            of the nudging (dynamical relaxation) term with respect
C            to the physical forcing term
C   DSRETC,- initial coefficient of variation, in the EnKF/SITR scheme,
C   DSKS,    for the retention curve parameters, the saturated hydraulic
C   DSSTOR,  conductivity (x,y), the elastic storage coefficient, the porosity,
C   DSPOROS, the surface routing parameters, and the measurements,
C   DSSURF,  the saturated hydraulica conductivity (z), respectively
c   DSKSZ
C  STORE_SAT - saturation volume calculated with the porosity in
C              input
C real scalars for transport modelling 
C ------------------------
C DELTATADV  - advective time step
C SUPDIFFUS  - sup of norm of the diffusion matrix to define 
C              the Peclet number 
C CFLNUMB    - CFL number (if CFLNUMB =0 the FV scheme is not used)
C PECNUMB    - Peclet number as calculated by cflpecnumber.f
C CFLINP     - real value, target cfl number when advective time step is
C              calculated automatically (NADVFL=1)
C DIFFUS     - molecular diffusion coefficient      
C GAMMAS     - density of solid material     
C TETAC      - weighting parameter for UT3D time stepping scheme for subsurface 
C              diffusion/dispersion transport model 
C
C Miscellaneous Real Scalars
C --------------------------  
C   AREATOT        - total area of the catchment surface
C   VOLTOT         - total volume of the discretized catchment
C   FHORT          - fraction of surface nodes that are Horton saturated
C   FDUNN          - fraction of surface nodes that are Dunne saturated
C   FPOND          - fraction of surface nodes that are ponded
C   FSAT           - fraction of surface nodes that are saturated or
C                    ponded
C   VOLUME_OUT     - ...
C   VOLUME_SUP     - ...
C   VOLSUPNEW      - ...
C   EVAP_EFF       - ...
C   INF_TOT        - ...
C   RMAX           - largest double precision # (machine dependent)
C   RMIN           - smallest double precision # (machine dependent)
C
C Real*4 Scalars for Cpu Timing
C -----------------------------
C   CPUSUB    - cpu time for subsurface flow module
C   CPUSUB_T  - total cpu time for subsurface flow module
C   CPUSURF   - cpu time for surface routing module
C   CPUSURF_T - total cpu time for surface routing module
C   CPUNL     - total cpu time for nonlinear scheme
C   CPUOVH    - total cpu time for overhead:
C               - data input, initialization, and output of initial
C                 conditions (once)
C               - construction of tetrahedral elements (once)
C               - volume calculations (once)
C               - set up of storage indices and pointers (once)
C               - velocity calculations (every time step for the
C                 case IPRT > 1)
C               - hydrograph calculation (every time step)
C               - input, interpolation, and switching control of
C                 atmospheric boundary conditions for the next time
C                 level (every time step)
C               - update of pressure heads for the next time level
C                 (every time step)
C               - back-stepping procedure (when needed)
C               - final output (once)
C   CPUMN     - total cpu time for the simulation 
C   CPUUPD    - cpu time for all the updates
C   CPUUPD1   - cpu time for a single update
C   CPUUPD2   - cpu time for all the updates over the number of
C                realization
C   CPUNLT    - ...
C   AVGNL     - ...
C   AVGLIN    - ...
C   AVGLNL    - ...
C   ATCTS     - ...
C   ATCNL     - ...
C   ANCTS     - ...
C   ANCNL     - ...
C   PCMN      - ...
C   PCNL      - ...
C   PCOVH     - ...
C   PCNLT     - ...
C   PCVEC1-9  - ...
C
C Real*4 Array for Cpu Timing
C ---------------------------
C   CPUVEC(11) - cpu times for different sections of nonlinear schemes:
C                (1) unsat characteristics
C                (2) initialization of system matrices
C                (3) assembly of local system components into global
C                    matrices
C                (4) calculation of RHS without boundary conditions
C                (5) construction of global LHS system matrix
C                (6) calculation of BC contributions to RHS
C                (7) linear solver and calculation of residual
C                (8) extraction of pressure head solution from the
C                    difference solution and re-setting of solution and
C                    of COEF1 for Dirichlet nodes
C                (9) back-calculation of fluxes at all Dirichlet nodes,
C                    mass balance calculation, application of nonlinear
C                    relaxation scheme (if required), calculation of
C                    nonlinear convergence and residual error norms,
C                    switching control of atmospheric BCs, and
C                    calculation of new position of the exit point along
C                    each seepage face
C               (10) ...timer for SURF_ROUTE...
C               (11) ...timer for nudging...
C  
C Integer Arrays for Grid, BCs, Outputs, and System Matrices
C ----------------------------------------------------------
C   TP    (N)          - # of elements connecting to each node
C   TP2D  (NNOD)       - ...
C   TRIANG(4,NTRI)     - element connectivities in 2-d mesh (TRIANG(4,I)
C                        indicates material type for 2-d element I)
C   CELL  (5,NCELL)    - cell connectivities in DEM to 2-d mesh
C                        (CELL(5,I) indicates material type for cell I)
C   TETRA (5,NT)       - element connectivities in 3-d mesh (TETRA(5,I)
C                        indicates material type for 3-d element I)
C   IVOL  (NT)         - sign of the volume of each element
C   CONTP2(NDIR)       - non-atmospheric, non-seepage face Dirichlet  
C                        node #'s in 2-d mesh
C   CONTP (3,NP)       - non-atmospheric, non-seepage face Dirichlet 
C                        node #'s in 3-d mesh
C   CONTQ (3,NQ)       - non-atmospheric, non-seepage face Neumann 
C                        node #'s in 3-d mesh
C   ACONTP(ANP)        - actual values of CONTP and CONTQ that are 
C   ACONTQ(ANQ)          passed to the FLOW3D module
C   IFATM (NNOD)       - IFATM(I)=0  if surface node I is a Neumann 
C                        atmospheric boundary condition node 
C                        IFATM(I)=1  if surface node I is a Dirichlet
C                        atmospheric boundary condition node 
C                        IFATM(I)=2  if surface node I is "ponded",
C                        treated as a Dirichlet atmospheric boundary 
C                        condition node
C                        IFATM(I)=-1 if surface node I is not an
C                        atmospheric boundary condition node 
C                        Note: surface nodes are numbered 1,...,NNOD in
C                        the 3-d mesh, so there is no need for a pointer
C                        array giving the node #'s for the surface nodes
C   IFATMP(NNOD)       - IFATM values at previous time level
C   SATSUR(NNOD)       - SATSUR(I)=1  if surface node I is unsaturated
C                        SATSUR(I)=2  if surface node I is Horton 
C                        saturated (infiltration excess mechanism)
C                        SATSUR(I)=3  if surface node I is Dunne 
C                        saturated (saturation excess mechanism)
C   NSFNUM(NSF)        - # of nodes on each seepage face
C   NSFNOD(NSF,NNSFMX) - node #'s on each seepage face. The node #'s for
C                        each seepage face must be input in descending
C                        order by elevation. That is, along seepage face
C                        I, Z(NSFNOD(I,J)) .GE. Z(NSFNOD(I,J+1)) must 
C                        hold for J=1,...,NSFNUM(I)-1. Seepage faces can
C                        be defined, for instance, above a well, along a
C                        stream bank, or along a combination of stream
C                        bank and surface nodes. For a configuration 
C                        of seepage face and stream, the stream nodes 
C                        should be designated as non-atmospheric, 
C                        non-seepage face nodes, for instance as 
C                        Dirichlet nodes with a pressure head 
C                        distribution in hydrostatic equilibrium, the 
C                        node at the surface of the stream being 
C                        assigned a pressure head value of zero. 
C                        For output purposes, we set 
C                        NSFNOD(I,NSFNUM(I)+1)=-9999.
C   SFEX  (NSF,NNSFMX) - the exit point on each seepage face. The 
C                        seepage face nodes above the exit point are 
C                        'potential' seepage face nodes, are treated as
C                        zero flux Neumann BC's, and the pressure heads
C                        here should be negative (unsaturated). The 
C                        seepage face nodes below the exit point (and
C                        including the exit point) are 'actual' seepage
C                        face nodes, are treated as zero pressure 
C                        head Dirichlet BC's (saturated), and the 
C                        back-calculated fluxes here should be
C                        negative (outflow). 
C                        The position of the exit point for the 
C                        first time step is calculated from the initial 
C                        conditions. The new position of the exit point 
C                        is calculated after every nonlinear iteration 
C                        of every time step, and the boundary conditions
C                        for the seepage face nodes are adjusted to
C                        reflect changes in the position of the exit 
C                        point. 
C                        For the case where seepage face I is 
C                        completely saturated (all seepage face 
C                        nodes are 'actual'), SFEX(I)=1. For the case
C                        where seepage face I is completely unsaturated
C                        (all seepage face nodes are 'potential' and 
C                        there is no exit point), SFEX(I)=NSFNUM(I)+1.
C                        This convention simplifies the handling of 
C                        seepage face nodes (relying on the fact that 
C                        FORTRAN 77 does not execute a DO loop if the 
C                        iteration count is zero or negative).
C   SFEXIT(NSF,NNSFMX) - SFEX values at previous nonlinear iteration
C   SFEXP (NSF,NNSFMX) - SFEX values at previous time level
C   NODDIR(NUMDIR)     - node #'s for all Dirichlet nodes in 3-d mesh
C   CONTR (NR)         - node #'s for partial output
C   NODVP (NUMVP)      - node #'s for surface nodes selected for 
C                        vertical profile output
C   IA    (NTERM)      - row indices in storage of system matrices for
C                        nonsymmetric case
C   JA    (N1*N)       - column indices (in ascending order) in storage
C                        of system matrices
C   TOPOL (N+1)        - pointer to first nonzero element of each row
C                        which is stored in the system matrices (the
C                        diagonal entry in symmetric storage case)
C   TETJA (4,4,NT)     - gives the index within JA (global position) of
C                        each component of the 4 x 4 local system
C                        matrices (upper triangle of 4 x 4 arrays only
C                        for symmetric case). For the symmetric case,
C                        construction of TETJA requires the reordering
C                        of the nodes of the tetrahedra, so that
C                        assembly of the system matrices can be done
C                        by indexing directly the upper triangular
C                        part of the matrix.
C   IP3   (3,3)        - 3 x 3 permutation matrix
C   IP4   (4,4)        - 4 x 4 permutation matrix
C   SFFLAG(5)          - counter for anomalous, implausible, or 
C                        erroneous occurrences along seepage faces
C                        (see output statements 2100,2200 in subroutine
C                        SFINIT, 2100,2200,2300,2400 in subroutines
C                        EXTONE, EXTALL, and 2500 in subroutine FLUXMB)
C   HGFLAG(9)          - counter for anomalous, implausible, or 
C                        erroneous atmospheric inflow, outflow, and 
C                        runoff occurrences (see subroutine HGRAPH)
C   IER   (7)          - error flags 
C   INSYM (IBOT+6N+1)  - integer scratch vector for NONSYM solver
C 
C Integer Arrays for SURF_ROUTE
C -----------------------------
C   ID_QOUT         (NUM_QOUT)     - I_BASIN index of selected cells for
C                                    hydrograph output
C   ZONE            (NROW,NCOL)    - map of different material zones
C   LAKES_MAP       (NROW,NCOL)    - map of lakes in the DEM. 
C                                    =0 if no lake cell; otherwise lake cell 
C   BASE_MAP        (NROW,NCOL)    - raster map of the catchment impermeable
C                                    basement. Each value represents the 
C                                    thickness of the corresponding cell.
C   DTM_P_OUTFLOW_1 (NCOL,NROW)    - raster map of cardinal flow directions.
C   DTM_P_OUTFLOW_2 (NCOL,NROW)    - raster map of diagonal flow directions.
C                                    3 NW, 6 N, 9 NE
C                                    2 W ,      8 E
C                                    1 SW, 4 S, 7 SE
C                                    |----|----|----|
C                                    | 3  | 6  |  9 |
C                                    |----|----|----|
C                                    | 2  | -  |  8 |
C                                    |----|----|----|
C                                    | 1  | 4  |  7 |
C                                    |----|----|----|
C   NODI            (NROW+1,NCOL+1)- raster containing the numbering of surface
C                                    nodes.
C   INDCEL          (NROW,NCOL)    - cell # for each active cell in the DEM,
C                                    excluding lake cells
C   INDCELWL        (NROW,NCOL)    - cell # for each active cell in the DEM,
C                                    including lakes cells
C   QOI_SN          (NCELL)        - vector containing the I_BASIN index. 
C                                    The cells are ordered in a descending 
C                                    elevation order.
C                                    I_BASIN indexes DEM cells as calculated
C                                    in the pre-processor, in the following way:
C                                    |----|----|----|---|
C                                    | 3  | 6  | 9  |12 |
C                                    |----|----|----|---|
C                                    | 2  | 5  | 8  |11 |
C                                    |----|----|----|---|
C                                    | 1  | 4  | 7  |10 |
C                                    |----|----|----|---|          
C   CELTYPE         (NCELL)        - cell type. =0 if channel cell; otherwise
C                                    reservoir cell of type CELTYPE.          
C   CELLCOL         (NCELL)        - column # for each cell
C   CELLROW         (NCELL)        - row # for each cell
C   CELLCOL_WL      (NCELL)        - column # for each cell counting lake cells
C   CELLROW_WL      (NCELL)        - row # for each cell counting lake cells
C   CELLS_R         (NUMRES)       - ...
C   MODIF           (NCELL)        - ...
C   RESERVR         (NCELL)        - reservoir #. =0 if not a reservoir cell.
C   N_HA            (NUMRES)       - ...
C 
C Integer Arrays for Nudging and EnKF
C -----------------------------------
C   NUDTET(NUDN) - pointer to the element within which the observation
C                  points for nudging are located
C Integer Arrays for transport modelling
C --------------------------
C   SIDE_CNC(4,NT) - element-face connectivities for each tetrale
C   ISIDE(3,NFACE)     - global numbering of nodes that are adjacent
C                        to each face counting from south to north
C   PLIST(2,NFACE)     - west and east element number w.r.t oriented
C                        normal; west and east face local number
C   NEIGH(4,NT)    - list of neighbouring elements to each element
C   TETRAVERT(NNODE,N1)- matrix of elements connected to each node
C   TVERT(NNODE)       - vector of numbers of elements connected to
C                        each node
C   Node_ElementIds(N,N1)- matrix of elements connected to each node
C   Node_ElementCount(N) - vector of numbers of elements connected to
C                        each node
C   GLOB_TO_INT(NFACE) - pointer from global to internal face numbering
C   NODE_FACECOUNT(N)  - number of face connected to a node
C   NODE_FACEIDS(N,MAXFCONTONODE)  - faces # connected to each node
C   IBNDRY(NFACE)      - boundary face flag:
C                        = -1 internal face
C                        =  0 non-ATM non-SF Neumann b.c. face
C                        = 10 ATM Neumann b.c. face
C                        = 20 SF Neumann b.c. face
C                        =  1 non-ATM non-SF Dirichlet b.c. face
C                        = 11 ATM Dirichlet b.c. face
C                        = 21 SF Dirichlet b.c. face
C   IBDEDGE(N)         - boundary node to boundary edge conversion
C                        assuming that given a node I on the boundary,
C                        IBDEDGE(I) is such that I=ISIDE(1,IBDEDGE(I))
C   PTRGEBC(NFACE)     - pointer from global face to boundary condition
C
C CONTP_TRA (3,NP_TRA)  - TRANSPORT Dirichlet  node #'s in 3-d mesh
C CONTPFA_TRA           - transport dirichlet face #'s in 3d mesh
C
C ACONTP_TRA(ANP_TRA)  - actual values of CONTP_TRA that is passed to the
C                        transport model
C
C  CONTQ_TRA(NQMAX_TRA)  - atmospheric cauchy boundary condition node's number
C                        (equals to the surface node's number as we impose a
C                        Cauchy BC for all the surface node)
C  QCAUCHY(NNOD)       - total cauchy bc values (after the switching)
C  Faces_FaceType(NFACEMAX) - Type of Faces For Velocity Reconstruction:
C                          -1  Internal
C                           0  Dirichlet
C                           2  Neumann
C  CTIM  (3)           - most current input time values for
C                        TRANSPORT Dirichlet
C                        BC's, with CTIM(1) < CTIM(2) < CTIM(3)
C                        and CTIM(2) < TIME <= CTIM(3)
C   CINP  (3,NP_TRA)       - TRANSPORT Dirichlet
C                        values corresponding to CTIM times.
C                        PRESC_TRA(I) is obtained from CINP(2,I) and
C                        CINP(3,I) by linear interpolation (not any more, now
C                        by piecewise constant function).
C                        CINP(1,I) values are needed in the event that,
C                        after back-stepping, we have
C                        CTIM(1) < TIME <= CTIM(2)
C
C Real Arrays for Mesh Configuration and Boundary Conditions
C ----------------------------------------------------------
C   X     (N)          - x-coordinates (for 2-d mesh on input)
C   Y     (N)          - y-coordinates (for 2-d mesh on input)
C   Z     (N)          - z-coordinates (surface elevation values on 
C                        input - see description of ISP)
C   DEPTH (NNOD)       - depth value for the 2-d mesh on input, read 
C                        only if IVERT=3
C   XC    (NT)         - x-coord at the centroid of each tetrahedra
C   YC    (NT)         - y-coord at the centroid of each tetrahedra
C   ZC    (NT)         - z-coord at the centroid of each tetrahedra
C   VOLNOD(N)          - absolute value of volume assigned to each node
C   VOLU  (NT)         - absolute value of the volume of each element
C   VOLUR (NT)         - reciprocal of VOLU
C   ZRATIO(NSTR)       - fraction of total grid height that each layer
C                        is to occupy (see also description of IVERT). 
C                        ZRATIO(1) is for the surface-most layer. 
C                        ZRATIO values must sum to 1.
C   PRESC (NP)         - non-atmospheric, non-seepage face Dirichlet 
C                        values at current time level
C   PTIM  (3)          - most current input time values for
C                        non-atmospheric, non-seepage face Dirichlet
C                        BC's, with PTIM(1) < PTIM(2) < PTIM(3)
C                        and PTIM(2) < TIME <= PTIM(3)
C   PINP  (3,NP)       - non-atmospheric, non-seepage face Dirichlet 
C                        values corresponding to PTIM times.
C                        PRESC(I) is obtained from PINP(2,I) and
C                        PINP(3,I) by linear interpolation (not any more, now
C                        by piecewise constant function).
C                        PINP(1,I) values are needed in the event that,
C                        after back-stepping, we have
C                        PTIM(1) < TIME <= PTIM(2)
C   Q     (NQ)         - non-atmospheric, non-seepage face Neumann 
C                        values at current time level
C   QTIM  (3)          - most current input time values for
C                        non-atmospheric, non-seepage face Neumann
C                        BC's, with QTIM(1) < QTIM(2) < QTIM(3)
C                        and QTIM(2) < TIME <= QTIM(3)
C   QINP  (3,NP)       - non-atmospheric, non-seepage face Neumann 
C                        values corresponding to QTIM times.
C                        Q(I) is obtained from QINP(2,I) and
C                        QINP(3,I) by linear interpolation.
C                        QINP(1,I) values are needed in the event that,
C                        after back-stepping, we have
C                        QTIM(1) < TIME <= QTIM(2)
C   QPNEW (NP)         - back-calculated flux values at non-atmospheric,
C                        non-seepage face Dirichlet nodes at current 
C                        time level
C   QPOLD (NP)         - QPNEW values at previous time level
C   QTRANIE(N)         - root-zone water uptake at current time level
C                        always positive, changes sign in BCPIC
C   ZROOT(VEG_TYPE)    - depth of root zone for every vegetation type,
C                        assigned on the basis of VEG_TYPE, read as
C                        raster in IIN3 or in IIN2 after grid info
C   SCF                - soil cover fraction (fraction of soil covered
C                        in vegetation)
c STR_ZON(ZONE_VERT)   - for each vertical zone contain the last layer
C                        (str) of the zone. Examples: if ZONE_VERT=1,
C                        then STR_ZON(1)=NSTR; if ZONE_VERT=NSTR, then
C                        STR_ZON=[2 3 4 ... NSTR]. Used if DAFLAF>0
C   Notes: (a) For a simulation using temporally homogeneous 
C              non-atmospheric, non-seepage face Dirichlet
C              (Neumann) BC's, input data on unit IIN8 (IIN9) should
C              contain a single value of PTIM (QTIM) (0.0)
C              and a single set of PINP (QINP) data.
C              Alternatively, to properly handle the case where the
C              datasets for different simulations are kept in the same
C              file (separated by blank lines), the input data on unit
C              IIN8 (IIN9) for temporally homogeneous
C              non-atmospheric, non-seepage face Dirichlet
C              (Neumann) BC's should contain, as above, a value of PTIM
C              (QTIM) of 0.0 followed by the PINP
C              (QINP) values, and then a value of PTIM
C              (QINP) equal to or larger than TMAX
C              (1.0e+15, say) followed by the same PINP (QINP)
C              values specified at time 0.0.
C          (b) If the first input time value is greater than 0.0, we set
C              the initial (time 0.0) non-atmospheric, non-seepage face
C              Dirichlet (Neumann) BC inputs to 0.0
C          (c) If TIME is larger than the last PTIM (QTIM)
C              value on unit IIN8 (IIN9), HTIDIR (HTINEU)
C              is set to 1 and the last input values are used for the
C              rest of the simulation. To properly handle the case where
C              the datasets for different simulations are kept in the 
C              same file (separated by blank lines), follow the 
C              procedure described in (a).
C   SFQ   (NSF,NNSFMX) - back-calculated flux values at actual seepage
C                        face nodes at current time level
C   SFQP  (NSF,NNSFMX) - SFQ values at previous time level
C   ARENOD(NNOD)       - area assigned to each surface node
C                        (needed for conversion of atmospheric
C                        rainfall/evaporation rates to
C                        volumetric fluxes)
C   ATMPOT(NNOD)       - precipitation (+ve) / evaporation (-ve) fluxes
C                        at current time level for each surface node.
C                        These are potential infiltration/exfiltration
C                        values. 
C   ATMACT(NNOD)       - actual fluxes (infiltration/exfiltration
C                        values) for atmospheric boundary
C                        condition nodes at current time level. 
C                        For IFATM(I)=0,  ATMACT(I) = ATMPOT(I);
C                        For IFATM(I)=1,  ATMACT(I) = back-calculated
C                        flux value;
C                        For IFATM(I)=-1, ATMACT(I) is disregarded.
C   ATMOLD(NNOD)       - ATMACT values at previous time level
C   ATMTIM(3)          - most current input time values for atmospheric
C                        BC's, with ATMTIM(1) < ATMTIM(2) < ATMTIM(3)
C                        and ATMTIM(2) < TIME <= ATMTIM(3)
C   ATMINP(3,NNOD)     - input atmospheric rainfall/evaporation rates
C                        corresponding to ATMTIM times. ATMPOT(I) is
C                        obtained from ATMINP(2,I) and ATMINP(3,I) by
C                        linear interpolation and conversion of rate
C                        to volumetric flux. ATMINP(1,I) values are 
C                        needed in the event that, after back-stepping,
C                        we have ATMTIM(1) < TIME <= ATMTIM(2)
C   EFFTIM(2)          - most current input time values for effective
C                        rainfall inputs, with EFFTIM(1) < EFFTIM(2)
C                        and EFFTIM(1) < TIME <= EFFTIM(2).
C   Notes: (a) For a simulation using temporally homogeneous atmospheric
C              rates, input data on unit IIN6 should contain a single 
C              value of ATMTIM (0.0) and a single set of ATMINP data.
C              Alternatively, to properly handle the case where the
C              datasets for different simulations are kept in the same
C              file (separated by blank lines), the input data on
C              unit IIN6 for temporally homogeneous rates should
C              contain, as above, a value of ATMTIM of 0.0 followed by
C              the ATMINP rates, and then a value of ATMTIM equal to or
C              larger than TMAX (1.0e+15, say) followed by the same
C              ATMINP rates specified at time 0.0.
C          (b) If there is no ATMTIM, ATMINP input, HTIATM is set 
C              to 1 (homogeneous in time) and atmospheric input rates
C              are set to 0.0.
C          (c) If the first input time value is greater than 0.0, we set
C              the initial (time 0.0) atmospheric input rates to 0.0.
C          (d) If TIME is larger than the last ATMTIM value on unit
C              IIN6, HTIATM is set to 1 and the last input atmospheric
C              rates are used for the rest of the simulation. To
C              properly handle the case where the datasets for different
C              simulations are kept in the same file (separated by blank
C              lines), follow the procedure described in (a).
C          (e) If HSPATM is nonzero and not equal to 9999 (spatially
C              homogeneous), each set of ATMINP data should consist of
C              a single value which gets copied to all surface nodes. 
C              If HSPATM is zero (spatially variable), each set of 
C              ATMINP data should consist of NNOD values (note that
C              we read in a value for each surface node, including
C              surface nodes which may be designated as non-atmospheric
C              Dirichlet or Neumann boundary conditions. IFATM controls
C              whether the atmospheric input for a given surface node
C              is actually used). 
C   KD(NSTR,NZONE)          - RETARDATION COEFFICIENT
C   LAMBDA(NSTR,NZONE)      - half-life characterizing first order degradation
C
C   PRESC_TRA (NP_TRA)  - transport Dirichlet (fixed concentration)
C                        values at current time level
C
C    PRESCFA_TRA (NPFA_TRA) -  non-atmospheric, non-seepage face, Dirichlet
C                              values for faces
C
C    ATMCONC(NNOD)  - imposed concentration in the atmospheric forcing - read
C                     in the input file transp_atmbc
C
C   QMASS_IN_KK_SN          (NCELL)    - mass discharge entering the cells a
C                                    previous SURF_FLOWTRA time level
C   QMASS_IN_KK_SN_SAV      (NCELL)   - variable defined to store QMASS_IN_KK_SN
C                                    in case of back-stepping
C   QMASS_IN_KK_SN_P        (NCELL)    - mass discharge entering the cells at
C                                    previous FLOW3D time level
C   QMASS_IN_KKP1_SN        (NCELL)    - mass discharge entering the cells at
C                                    actual time
C   QMASS_OUT_KK_SN_1       (NCELL)    - mass discharge exiting the cells along
C                                    the cardinal direction at
C                                    previous SURF_FLOWTRA time level
C   QMASS_OUT_KK_SN_1_SAV  (NCELL) - variable defined to store QMASS_OUT_KK_SN_1
C                                    in case of back-stepping
C   QMASS_OUT_KK_SN_1_P     (NCELL)    - mass discharge exiting the cells along
C                                    the cardinal direction at
C                                    previous FLOW3D time level
C   QMASS_OUT_KKP1_SN_1     (NCELL)    - mass discharge exiting the cells along
C                                    the cardinal direction at
C                                    actual time
C   QMASS_OUT_KK_SN_2       (NCELL)    - mass discharge exiting the cells along
C                                    the diagonal direction at
C                                    previous SURF_FLOWTRA time level
C   QMASS_OUT_KK_SN_2_SAV  (NCELL) - variable defined to store QMASS_OUT_KK_SN_2
C                                    in case of back-stepping
C   QMASS_OUT_KK_SN_2_P     (NCELL)    - mass discharge exiting the cells along
C                                    the diagonal direction at
C                                    previous FLOW3D time level
C   QMASS_OUT_KKP1_SN_2     (NCELL)    - mass discharge exiting the cells along
C                                    the diagonal direction at
C                                    actual time
C
C Exchange Variables Between subsurface transport and surface transport 
C ------------------------------------------------
C   TRAFLNOD (NNOD)  - node-wise total mass flux from subsurface to surface
C                       sum of advection fluxes and mixing ("dispersion") fluxes
C   TRAFLNOD_FLOW (NNOD)  - node-wise mass flux produced by subsurface transport
C   TRAFLP   (NNOD)  - node-wise mass flux at previous time level
C   CONCNOD (NNOD)   - node-wise surface concentration produced by surface
C                      transport
C   CONCNODUPD (NNOD)- node-wise surface concentration updated to account for
C                      the effect of rainfall or evaopration
C   TRAFLCEL (NCELL) - cell-wise mass flux produced by subsurface transport
C   CONCCEL (NCELL)  - cell-wise surface concentration produced by surface
C                      transport
C   ACTCEL  (NCELL) - cell-wise actual flux from subsurface (ATMACT
C                     transferred to the cells)
C   SOURCE_MIXING(NNOD) - node-wise mass flux produced by mixing
C   MIX_CORRECT(NCELL) - cell-wise correction vector after mixing; filled out
C                        in SURF_FLOWTRA and used in MIXING_CORRECTION to
C                        correct first layer subsurface concentration.
C
C Exchange Variables Between FLOW3D and SURF_ROUTE
C ------------------------------------------------
C   OVFLNOD (NNOD)  - node-wise overland flux produced by FLOW3D
C   OVFLP   (NNOD)  - node-wise overland flux at previous time level
C   PONDNOD (NNOD)  - node-wise ponding pressure head produced by
C                     SURF_ROUTE
C   OVFLCEL (NCELL) - cell-wise overland flux produced by FLOW3D
C   PONDCEL (NCELL) - cell-wise ponding pressure head produced by
C                     SURF_ROUTE
C   ACTCEL  (NCELL) - cell-wise actual flux from subsurface (ATMACT
C                     transferred to the cells)
C   CELLCOARSE(NCELL) - ...can probably be eliminated...
C
C Real Arrays for SURF_ROUTE
C -----------------------
C   DEM_MAP             (NROW,NCOL)    - raster map of the catchment surface. 
C                                        Each value is the elevation of the 
C                                        corresponding cell.
C   DTM_W_1             (NCOL,NROW)    - raster map of the weights. 
C                                        Each value is the weight assigned to 
C                                        the cardinal drainage direction.
C   DTM_W_2             (NCOL,NROW)    - raster map of the weights. 
C                                        Each value is the weight assigned to 
C                                        the diagonal drainage direction.
C   DTM_LOCAL_SLOPE_1   (NCOL,NROW)    - raster map of local slopes. 
C                                        Each value is the local slope assigned
C                                        to the cardinal drainage direction.
C   DTM_LOCAL_SLOPE_2   (NCOL,NROW)    - raster map of local slopes. 
C                                        Each value is the local slope assigned
C                                        to the diagonal drainage direction.
C   DTM_EPL_1           (NCOL,NROW)    - raster map of elemental path length. 
C                                        Each value is the length of the 
C                                        elemenatal path assigned to the 
C                                        cardinal drainage direction.
C   DTM_EPL_2           (NCOL,NROW)    - raster map of elemental path length. 
C                                        Each value is the length of the 
C                                        elemental path assigned to the 
C                                        diagonal drainage direction.
C   DTM_KSS1_SF_1       (NCOL,NROW)    - raster map of roughness coefficient. 
C                                        Each value is the roughness 
C                                        coefficient assigned to the 
C                                        cardinal drainage direction.
C   DTM_KSS1_SF_2       (NCOL,NROW)    - raster map of roughness coefficient. 
C                                        Each value is the roughness 
C                                        coefficient assigned to the 
C                                        diagonal drainage direction.
C   DTM_WS1_SF_1        (NCOL,NROW)    - raster map of surface water width. 
C                                        Each value is the surface water width 
C                                        of the cells scaled in space and 
C                                        assigned to the cardinal 
C                                        drainage direction.
C   DTM_WS1_SF_2        (NCOL,NROW)    - raster map of surface water width. 
C                                        Each value is the surface water width 
C                                        of the cells scaled in space and 
C                                        assigned to the diagonal 
C                                        drainage direction.
C   DTM_B1_SF           (NCOL,NROW)    - raster map of exponent of the 
C                                        at-a-station relationship for 
C                                        the water surface width 
C                                        (i.e., surface water width is scaled 
C                                        in time with a power-law function 
C                                        with that exponent).
C   DTM_Y1_SF           (NCOL,NROW)    - raster map of exponent of the 
C                                        at-a-station relationship for 
C                                        the roughness coefficient 
C                                        (i.e., roughness coefficient is scaled
C                                        in time with a power-law function 
C                                        with that exponent).
C   DTM_NRC             (NCOL,NROW)    - raster map of containing the number 
C                                        of rivulets per hillslope cell.
C-----------------------
C   Real Scalar Variables for Transport MASS BALANCE Calculation
C   (CARLOTTA)
C   MASST - mass of solute in the subsurface at each time step
C            calculated by integrating moisture content*concentration
C            (SW_NODE * CNEW_EL*PNODI)
C            over the entire domain
C   MASS0 - initial (time 0) mass of solute in the subsurface
C   MASSINADV - total mass flown in through advection
C           across the domain boundary at each time step
C   MASSOUTADV - total mass flown out through advection
C            from the domain boundary at each time step
C   MASSINADVTOT - time-cumulated mass total mass flown in through advection
C            across the domain boundary
C   MASSOUTADVTOT - time-cumulated mass total mass flown out through advection
C            from the domain boundary
C   MASSSURF - total mass in the surface at the current time step
C----------------------------------
C   real arrays for transport modeling and velocity reconstruction
C    
C   AREFACE(NFACE) - Surface Area of Each Face
C   Faces_ReferenceVector(3,NFACE) - Normal Vector of each face
C   FaceCentroid(3,NFACEMAX) - X,Y,Z coordinates of the Centroid of each
C                              Face
C   Faces_NeumannFlux (NFACE) - Neumann condition for velocity
C                               reconstruction on Neumann Face
C   Faces_FaceFlux(3,NFACE)   - Three component of the flux
C                               on each face
C
C-----------------------------
C Real Arrays for Nudging
C ------------------------------------
C   NUDTIM (NUDT)      - observation times
C   NUDVAL (NUDC,NUDN) - soil moisture content observation values
C   NUDSMC (NUDN)      - computed soil moisture or pressure head values,
C                        interpolated to the observation points
C   NUDTAU (NUDT,NUDN) - half period of the temporal influence window
C                        for the Cressman-type weighting function W(t)
C                        or temporal integral scale of the correlation
C                        function W(t)
C   NUDRXY (NUDN)      - horizontal radius of influence for the
C                        Cressman-type weighting function W(x,y)
C                        or Gaussian correlation function W(x,y)
C   NUDRZ  (NUDN)      - vertical radius of influence for the
C                        Cressman-type weighting function W(z)
C                        or Gaussian correlation function W(z)
C   NUDEPS (NUDN)      - quality factor "epsilon" of the observation
C                        data
C   NUDX   (NUDN)      - x-coordinates of the observation data points
C   NUDY   (NUDN)      - y-coordinates of the observation data points
C   NUDZ   (NUDN)      - z-coordinates of the observation data points
C   NUDDIF (NUDN)      - dynamical (Newtonian) relaxation or nudging
C                        term without the weighting function
C                        contribution (ie, NUDG * NUDEPS * (NUDVAL -
C                        NUDSMC))
C   NUDNOD (N)         - total nudging term contribution (with
C                        weighting functions) to the RHS system vector
C                        (units L^3/T)
C   NUDCUM (N)         - cumulative (over time) NUDNOD
C
C Real Arrays for Material Properties and Hydraulic Characteristics
C -----------------------------------------------------------------
C   KS(NT)             - saturated hydraulic conductivity-xx in each
C                        element
C   KZ(NT)             - saturated hydraulic conductivity-zz in each
C                        element
C   KZNOD(N)           - saturated hydraulic conductivity-zz in each
C                        node
C   PERMX (NSTR,NZONE) - saturated hydraulic conductivity-xx
C   PERMY (NSTR,NZONE) - saturated hydraulic conductivity-yy
C   PERMZ (NSTR,NZONE) - saturated hydraulic conductivity-zz
C   ELSTOR(NSTR,NZONE) - specific storage
C   POROS (NSTR,NZONE) - porosity (moisture content at saturation)
C   PEL   (NT)         - porosity at each element
C   SEL   (NT)         - elastic storage at each element
C   SNODI (N)          - specific storage at each node
C   PNODI (N)          - porosity at each node
C   PORE  (N)          - variable porosity in case of peat deformation
C   INDE  (N)          - void ratio at each node in case of peat deformation
C   INDE0 (N)          - void ratio at each node at the previous time step
C                        (initial void ratio for the first time step)
C   DEF   (N)          - vertical percentual deformation at each node
C   ETAI  (N)          - overall storage coefficient (general storage
C                        term) at each node
C   CKRW  (N)          - relative hydraulic conductivity at each node
C   CKRWP (N)          - relative hydraulic conductivity at each node at
C                        previous time level
C   SW    (N)          - water saturation (moisture content/porosity)
C                        at each node
C   SENODI(N)          - effective saturation ((moisture content -
C                        residual moisture content) / (porosity -
C                        residual moisture content)) at each node
C   DCKRW (N)          - derivative of CKRW wrt press. head at each node
C   DETAI (N)          - derivative of ETAI wrt press. head at each node
C   ETAE  (NT)         - overall storage coefficient for each element
C   CKRWE (NT)         - relative conductivity for each element
C   SWE   (NT)         - water saturation for each element
C   SEELT (NT)         - effective saturation for each element
C   DCKRWE(NT)         - derivative of CKRW wrt p. head for each element
C   DETAIE(NT)         - derivative of ETAI wrt p. head for each element
C 
C Real Arrays for Pressure Heads, Velocities, and Outputs
C -------------------------------------------------------
C   PNEW  (N)    - pressure heads at current time level, current
C                  nonlinear iteration
C   POLD  (N)    - pressure heads at previous nonlinear iteration
C   PDIFF (N)    - difference in pressure heads between nonlinear
C                  iterations 
C   PTNEW (N)    - weighted pressure heads at current nonlinear
C                  iteration
C   PTOLD (N)    - weighted pressure heads at previous nonlinear
C                  iteration
C   PTIMEP(N)    - pressure heads at previous time level (initial
C                  conditions on input; initial guess for steady state
C                  case)
C   UNOD  (N)    - Darcy velocity-x at each node
C   VNOD  (N)    - Darcy velocity-y at each node
C   WNOD  (N)    - Darcy velocity-z at each node
C   UU    (NT)   - Darcy velocity-x for each element
C   VV    (NT)   - Darcy velocity-y for each element
C   WW    (NT)   - Darcy velocity-z for each element
C   RECNOD(NNOD) - recharge flux computed at nodes immediately above the
C                  water table, reconducted to surface nodes.
C   TIMPRT(NPRT) - time values for detailed output. Detailed output is
C                  produced at initial conditions (TIME=0), at time
C                  values indicated in TIMPRT, and at the end of the
C                  simulation (TIME=TMAX). Detailed output consists of:
C                  values of pressure head, velocity, water saturation,
C                  and relative conductivity (depending on setting of
C                  IPRT) at all nodes; velocity, and water saturation
C                  (depending on setting of IPRT) at all elements;
C                  vertical profiles of pressure head, water saturation,
C                  and relative conductivity for the NODVP surface
C                  nodes; pressure head, water saturation, and SATSUR
C                  values at the surface nodes.
C   
C Real Arrays for System Components
C ---------------------------------
C   AI    (4,NT)       - coefficients 'a-i / 6' of the basis functions
C   BI    (4,NT)       - coefficients 'b-i / 6' of the basis functions
C   CI    (4,NT)       - coefficients 'c-i / 6' of the basis functions
C   DI    (4,NT)       - coefficients 'd-i / 6' of the basis functions
C   LMASS (4,4)        - local mass matrix, without the storage 
C                        coefficient ETAE and without the volume term
C   COEF1 (NTERM)      - global stiffness matrix; also used to store
C                        the LHS system matrix, which is the Jacobian
C                        in the Newton case
C   COEF2 (NTERM)      - global mass matrix
C   COEF3 (NTERM)      - derivative term components of the Jacobian for 
C                        Newton scheme; also used as a scratch vector
C   SCR   (NTERM)      - scratch vector
C   RNSYM (IBOT)       - real scratch vector for NONSYM solver
C   TNOTI (N)          - RHS system vector
C   BKCFLOW_NODE (N)   - Back calculated nodal flow
C   BKCFLOWP_NODE (N)  - Back calculated nodal flow at previous time
C                        level
C   BKCFLOWINT_NODE (N)- Back calculated nodal flow interpolated for
C                        Advective time level
C   XT5   (N)          - TNOTI before imposition of Dirichlet boundary
C                        conditions (needed for back-calculation of 
C                        fluxes used in mass balance calculations)
C   LHSP  (NP)         - values of diagonal elements of LHS system 
C                        matrix corresponding to non-atmospheric,
C                        non-seepage face Dirichlet nodes before 
C                        imposition of Dirichlet BC's (needed for 
C                        back-calculation of fluxes used in mass 
C                        balance calculations)
C   LHSSF (NSF,NNSFMX) - values of diagonal elements of LHS system 
C                        matrix corresponding to seepage face Dirichlet
C                        nodes before imposition of Dirichlet
C                        BC's (needed for back-calculation of fluxes
C                        used in mass balance calculations, and for 
C                        calculation of new position of the exit 
C                        point along each seepage face)
C   LHSATM(NNOD)       - values of diagonal elements of LHS system 
C                        matrix corresponding to atmospheric Dirichlet
C                        nodes before imposition of Dirichlet
C                        BC's (needed for back-calculation of fluxes
C                        used in mass balance calculations, and for 
C                        switching control of atmospheric BC's)
C 
C Real Scalars for Mass Balance   (defined in common block
C and Hydrograph Calculations      include file MB_HGRAPH.H)
C -----------------------------
C   NUDIN  - total "inflow"  flux contribution from the nudging term
C   NUDOUT - total "outflow" flux contribution from the nudging term
C   NUDINP - NUDIN at previous time level
C   NUDOUTP- NUDOUT at previous time level
C   VNUDIN - total "inflow"  volu from nudging term over curr time step
C   VNUDOUT- total "outflow" volu from nudging term over curr time step
C   VNUDTOT- cumulative (over all time steps) nudging term "flow"
C            volume VNUDIN + VNUDOUT
C   OVFLOW - total overland flow (surface runoff) flux produced at 
C            atmospheric surface nodes. Overland flow occurs during
C            rainfall periods when the actual flux is less than the 
C            potential flux, and accounts for both Horton and Dunne
C            saturation mechanisms.
C   REFLOW - total return flow flux produced at atmospheric surface
C            nodes. Return flow occurs during rainfall periods when
C            the actual flux is negative (outflow rather than inflow).
C            In this case all of the potential flux becomes overland 
C            flow, and the magnitude of the actual flux becomes the 
C            return flow component of surface runoff. 
C   RECFLOW- total recharge flux produced at nodes immediately above
C            the water table. We consider recharge only when the
C            Darcy velocity in the z-direction is negative [L^3/T].
C   RECVOL - same as RECFLOW, but in terms of cumulated volume [L^3].
C   SFFLW  - total subsurface flow flux produced at seepage faces
C            at the current time level
C   SFFLWP - total subsurface flow flux produced at seepage faces
C            at the previous time level
C   VSFFLW - total subsurface flow volume produced at seepage faces
C            between current and previous time levels 
C   APOT   - total atmospheric potential flux at the current time level,
C            used for hydrograph output. Note that we disregard 
C            contribution of non-atmospheric, non-seepage face
C            surface nodes in the calculation of APOT.
C   AACT   - total atmospheric actual flux at the current time level,
C            used for hydrograph output. AACT=ADIN+ADOUT+ANIN+ANOUT. 
C   AACTP  - AACT at the previous time level
C   AACTAV - average AACT between the current and previous time levels
C   ADIN   - tot inflow  flux from atmosph Dir nodes at curr time level
C   ADOUT  - tot outflow flux from atmosph Dir nodes at curr time level
C   ADINP  - tot inflow  flux from atmosph Dir nodes at prev time level
C   ADOUTP - tot outflow flux from atmosph Dir nodes at prev time level
C   NDIN   - tot inflow  flux from na, nsf Dir nodes at curr time level
C   NDOUT  - tot outflow flux from na, nsf Dir nodes at curr time level
C   NDINP  - tot inflow  flux from na, nsf Dir nodes at prev time level
C   NDOUTP - tot outflow flux from na, nsf Dir nodes at prev time level
C   ANIN   - tot inflow  flux from atmosph Neu nodes at curr time level
C   ANOUT  - tot outflow flux from atmosph Neu nodes at curr time level
C   ANINP  - tot inflow  flux from atmosph Neu nodes at prev time level
C   ANOUTP - tot outflow flux from atmosph Neu nodes at prev time level
C   NNIN   - tot inflow  flux from na, nsf Neu nodes at curr time level
C   NNOUT  - tot outflow flux from na, nsf Neu nodes at curr time level
C   NNINP  - tot inflow  flux from na, nsf Neu nodes at prev time level
C   NNOUTP - tot outflow flux from na, nsf Neu nodes at prev time level
C   VADIN  - tot inflow  volu from atmosph Dir nodes over curr time step
C   VADOUT - tot outflow volu from atmosph Dir nodes over curr time step
C   VNDIN  - tot inflow  volu from na, nsf Dir nodes over curr time step
C   VNDOUT - tot outflow volu from na, nsf Dir nodes over curr time step
C   VANIN  - tot inflow  volu from atmosph Neu nodes over curr time step
C   VANOUT - tot outflow volu from atmosph Neu nodes over curr time step
C   VNNIN  - tot inflow  volu from na, nsf Neu nodes over curr time step
C   VNNOUT - tot outflow volu from na, nsf Neu nodes over curr time step
C                non-atmospheric---^   ^---non-seepage face
C   VIN    - VADIN  + VNDIN  + VANIN  + VNNIN + VNUDIN  = total inflow
C            volume between current and previous time levels (> 0)
C   VOUT   - VADOUT + VNDOUT + VANOUT + VNNOUT + VSFFLW + VNUDOUT =
C            total outflow volume between current and previous time
C            levels (< 0)
C   VAPOT_T- cumulative (over all time steps) atmospheric potential flux
C            volume
C   VAACT_T- cumulative (over all time steps) atmospheric actual flux
C            volume
C   VSFTOT - cumulative (over all time steps) seepage face flow
C            volume VSFFLW
C   VNDTOT - cumulative (over all time steps) non-atmospheric,
C            non-seepage face Dirichlet flow volume VNDIN + VNDOUT
C   VNNTOT - cumulative (over all time steps) non-atmospheric,
C            non-seepage face Neumann flow volume VNNIN + VNNOUT
C   CVIN   - cumulative (over all time steps) net inflow volume VIN
C   CVOUT  - cumulative (over all time steps) net outflow volume VOUT
C   VTOT   - cumulative (over all time steps) total net flow volume 
C            VIN + VOUT
C   STORE1 - volume of water in the subsurface at each time step
C            calculated by integrating moisture content (SW * PNODI)
C            over the entire domain
C   STORE0 - initial (time 0) volume of water in the subsurface
C   STORE2 - volume of water in the subsurface at each time step
C            calculated by accumulating to STORE0 the current value
C            of DSTORE
C   DSTORE - total volume of storage change between current and 
C            previous time levels (> 0 for net increase in storage)
C   CDSTOR - cumulative (over all time steps) net storage change DSTORE
C            (> 0 for an overall global net increase in storage)
C   ERRAS  - volume ("mass") balance error over the current time step
C            (> 0 for (VIN + VOUT) > DSTORE)
C   CERRAS - cumulative (over all time steps) mass balance error ERRAS
C   CAERAS - cumulative (over all time steps) absolute value of the
C            mass balance error ERRAS
C   ERREL  - relative (percent) mass balance error over the current 
C            time step
C
C Integer and Real Arrays for   (defined in common block include file
C Calculation of Residual and    NORMVL.H used in subroutines FLOW3D
C Difference Norms               and CONVER)
C ---------------------------
C   ITUMAX(ITUNS) - nodes with largest pressure head difference in
C                   absolute value between current and previous
C                   nonlinear iterations for each FLOW3D nonlinear
C                   iteration
C   PIKMXV(ITUNS) - values of pressure head difference at node with
C                   largest pressure head difference in absolute value
C                   between current and previous iterations (for both
C                   FLOW3D and coupled FLOW3D/SURF_ROUTE)
C   PCURRV(ITUNS) - "current iteration" pressure head value used in the
C                   calculation of PIKMXV(ITUNS) values for each
C                   nonlinear iteration
C   PPREVV(ITUNS) - "previous iteration" pressure head value used in the
C                   calculation of PIKMXV(ITUNS) values for each
C                   nonlinear iteration
C   PL2V  (ITUNS) - values of the square root of the sum of squares of
C                   pressure head differences over all nodes (ie, L2
C                   norm of the convergence error), used in comparison
C                   with TOLUNS for convergence test in the case L2NORM
C                   nonzero, for each nonlinear iteration
C   FINFV (ITUNS) - values of residual error in the nonlinear FLOW3D
C                   solution calculated using the infinity norm (for the
C                   nonlinear system f(x)=0, the residual error at
C                   iteration "m" is the norm of f(x^m)) for each
C                   nonlinear iteration
C   FL2V  (ITUNS) - values of residual error in the nonlinear FLOW3D
C                   solution calculated using the L2 norm for each
C                   nonlinear iteration
C
C Real Scalars and Real Arrays for Unsaturated Soil   (defined in common
C Characteristics, for Chord and Tangent Slope         block include
C Formulas, and for Moisture Curve Lookup Table        file SOILCHAR.H)
C -------------------------------------------------
C   PMIN   - 'air dry' pressure head value (for switching control of 
C            atmospheric boundary conditions during evaporation)
C   CBETA0,- parameters  for Camporese adaptation of Pyatt and 
C   THETA0,  John relation for peat soil deformation
C   CANG 
C   PCANA,  Feddes parameters for computation of root water uptake.
C   PCREF,  For details, see Camporese, Daly, Paniconi, WRR 2015.
C   PCWLT,  Note that PCANA, PCREF, and PCWLT are given in terms
C   PZ,OMGC  of pressure head.
C   VGN,   - parameters for van Genuchten and extended van Genuchten 
C   VGM,     moisture curves (other 'VG' parameters - specific storage,
C   VGRMC,   porosity, and VGPNOT - are assigned nodally). VGM is 
C   VGPSAT   derived from VGN. VGRMC is residual moisture content.
C            For IVGHU=0, VGPNOT is (porosity - VGRMC)/porosity,
C            or (1 - residual water saturation).
C            For IVGHU=1, VGPNOT is a continuity parameter, derived by
C            imposing a continuity requirement on the derivative of 
C            moisture content with respect to pressure head. 
C   HUN,   - parameters for moisture curves from
C   HUA,     Huyakorn et al (WRR 20(8) 1984, WRR 22(13) 1986)
C   HUB,     (other 'HU' parameters - specific storage
C   HUALFA,  and porosity - are assigned nodally). HUN is
C   HUBETA,  only used for IVGHU=2; HUA and HUB are only used
C   HUGAMA,  for IVGHU=3. HUSWR is residual water saturation, which
C   HUPSIA,  is equivalent to residual moisture content/porosity.
C   HUSWR
C   BCBETA,- parameters for Brooks-Corey moisture curves (other 'BC' 
C   BCRMC,   parameters - specific storage and porosity - are assigned
C   BCPSAT   nodally). BCRMC is residual moisture content.
C   TOLKSL - tolerance for chord slope formula. Whenever the chord slope
C            formula is to be applied (for KSLOPE=1 or 2 at every iter-
C            ation and at all nodes; for KSLOPE=3 at those nodes whose
C            pressure heads fall within given ranges), it is applied
C            only if the absolute pressure head difference (between the
C            current and previous nonlinear iterations) is larger than
C            TOLKSL. If the difference is smaller than TOLKSL, then
C            differentiation is done either analytically (KSLOPE=1,3) or
C            with a centered difference formula (KSLOPE=2).
C   PKRL,  - left and right endpoints of the pressure head range within
C   PKRR     which the chord slope (case KSLOPE=3) or tangent slope
C            (case KSLOPE=4) formula is used to evaluate the derivative
C            of relative hydraulic conductivity
C   PSEL,  - left and right endpoints of the pressure head range within
C   PSER     which the chord slope (case KSLOPE=3) or tangent slope
C            (case KSLOPE=4) formula is used to evaluate the derivative
C            of effective saturation (moisture content for the case
C            of extended van Genuchten curves, IVGHU=1)
C   PDSE1L,- left and right endpoints of the two pressure head ranges
C   PDSE1R,  within which the chord slope (case KSLOPE=3) or tangent
C   PDSE2L,  slope (case KSLOPE=4) formula is used to evaluate the
C   PDSE2R   second derivative of effective saturation (moisture content
C            for the case of extended van Genuchten curves, IVGHU=1).
C            (Two ranges are specified since in general d(Se)/dP is
C            non-monotonic.)
C   DKRTAN - tangent slope approximation of d(Kr)/dP, the derivative of
C            relative hydraulic conductivity Kr wrt to pressure head P.
C            ie, DKRTAN = (Kr(PKRR) - Kr(PKRL))/(PKRR - PKRL)
C   VGPNOT  (N) - (porosity - VGRMC)/porosity for van Genuchten curves
C                 (IVGHU=0); continuity parameter 'PNOT' for extended
C                 van Genuchten curves (IVGHU=1)
C   BCPORM  (N) - (porosity - BCRMC)/porosity for Brooks-Corey curves
C                 (IVGHU=4)
C   DSETAN  (N) - tangent slope approximation of d(Se)/dP, the
C                 derivative of effective saturation Se (moisture
C                 content for case IVGHU=1) wrt to pressure head P.
C                 ie, DSETAN = (Se(PSER) - Se(PSEL)) /
C                                    (PSER - PSEL)
C   DDSE1T, (N) - tangent slope approximations of dd(Se)/dPP, the second
C   DDSE2T        derivative of effective saturation Se (moisture
C                 content for case IVGHU=1) wrt to pressure head P.
C                 ie, DDSE1T = (DSe(PDSE1R) - DSe(PDSE1L)) /
C                                    (PDSE1R - PDSE1L)
C                     DDSE2T = (DSe(PDSE2R) - DSe(PDSE2L)) /
C                                    (PDSE2R - PDSE2L)
C                 where DSe is the derivative of Se.
C                 (DSETAN, DDSE1T, and DDSE2T contain tangent slope
C                 values at each node only for the case IVGHU=1; for the
C                 other IVGHU cases the tangent slope values are
C                 constant for all nodes and are stored in DSETAN(1),
C                 DDSE1T(1), and DDSE2T(1).)
C   PCAP  (NSTR,NZONE, - pressure head values for the moisture curve
C          NLKP)         lookup table
C   SATC  (NSTR,NZONE, - water saturation values for the moisture curve
C          NLKP)         lookup table
C   KRWC  (NSTR,NZONE, - relative hydraulic conductivity values for the
C          NLKP)         moisture curve lookup table
C
C Real, Scalars and Real Arrays for   (defined in common
C Surface Water Routing and           block include files
C Watershed Characteristics           RIVERNETWORK.H and SURFWATER.H)
C --------------------------------
C   NORTH, SOUTH,    - DEM boundaries as given by GRASS
C   EAST, WEST
C   FACTOR           - multiplicative factor for DEM values (e.g. 
C                      to change the units of the elevation)
C   DELTA_X, DELTA_Y - cell dimensions (assumed equal!)
C   AK_MAX           - maximum wave celerity
C   AK_MAX_SAV       - variable defined to store AK_MAX 
C                      in case of back-stepping
C   AK_MAX_P         - maximum wave celerity at previous FLOW3D time
C                      level
C   CUTRGT           - Courant number "target" (defines the Courant
C                      criterion used in computing NSURF for the surface
C                      routing module; usually = 1.0; a very large
C                      value, eg 1.E+40, suppresses the Courant
C                      criterion so that NSURF=1, ie DTSURF = DELTAT
C                      where DTSURF is the time step size for the
C                      surface routing module)
C   CU_MAX           - maximum Courant number
C   SURFACE_WATER_INP   (2,NCELL)  - input local contribution to surface 
C                                    runoff such as read in the effective 
C                                    rainfall  input file (IIN22)
C   SURFACE_WATER_SN    (NCELL)    - SURFACE_WATER_INP multipled by the 
C                                    cell area (DELTA_X*DELTA_Y)
C   ELTRIA              (NTRI)     - elevation assigned to each triangle 
C                                    of the cell (the two triangles of the same
C                                    cell have the same cell elevation value)
C   H_WATER_KKP1_SN     (NCELL)    - surface water height calculated 
C                                    at the cell
C   Q_IN_KK_SN          (NCELL)    - discharge entering the cells a
C                                    previous SURF_ROUTE time level
C   Q_IN_KK_SN_SAV      (NCELL)    - variable defined to store Q_IN_KK_SN
C                                    in case of back-stepping
C   Q_IN_KK_SN_P        (NCELL)    - discharge entering the cells at
C                                    previous FLOW3D time level
C   Q_IN_KKP1_SN        (NCELL)    - discharge entering the cells at
C                                    actual time
C   Q_OUT_KK_SN_1       (NCELL)    - discharge exiting the cells along
C                                    the cardinal direction at
C                                    previous SURF_ROUTE time level
C   Q_OUT_KK_SN_1_SAV   (NCELL)    - variable defined to store Q_OUT_KK_SN_1
C                                    in case of back-stepping
C   Q_OUT_KK_SN_1_P     (NCELL)    - discharge exiting the cells along
C                                    the cardinal direction at
C                                    previous FLOW3D time level
C   Q_OUT_KKP1_SN_1     (NCELL)    - discharge exiting the cells along 
C                                    the cardinal direction at
C                                    actual time
C   Q_OUT_KK_SN_2       (NCELL)    - discharge exiting the cells along
C                                    the diagonal direction at
C                                    previous SURF_ROUTE time level
C   Q_OUT_KK_SN_2_SAV   (NCELL)    - variable defined to store Q_OUT_KK_SN_2
C                                    in case of back-stepping
C   Q_OUT_KK_SN_2_P     (NCELL)    - discharge exiting the cells along
C                                    the diagonal direction at
C                                    previous FLOW3D time level
C   Q_OUT_KKP1_SN_2     (NCELL)    - discharge exiting the cells along 
C                                    the diagonal direction at
C                                    actual time
C   Q_OVERLAND          (NCELL)    - lateral inflow
C   VOLUME_KK_SN        (NCELL)    - volume of surface water at previous
C                                    SURF_ROUTE time level
C   VOLUME_KK_SN_SAV    (NCELL)    - varible defined to store VOLUME_KK_SN
C                                    in case of back-stepping
C   VOLUME_KK_SN_P      (NCELL)    - volume of surface water at previous
C                                    FLOW3D time level
C   VOLUME_KKP1_SN      (NCELL)    - volume of surface water at actual
C                                    time
C   H_POOL_KK_VEC       (NUMRES)   - water levels in reservoirs at
C                                    previous SURF_ROUTE time level
C   H_POOL_KK_VEC_SAV   (NUMRES)   - variable defined to store H_POOL_KK_VEC
C                                    in case of back-stepping
C   H_POOL_KK_VEC_P     (NUMRES)   - water levels in reservoirs at
C                                    previous FLOW3D time level
C   H_POOL_KKP1_VEC     (NUMRES)   - water levels in reservoirs at
C                                    actual time
C   HRES                (NUMRES,30)- ...
C   ARES                (NUMRES,30)- ...
C   H_SOGLIA            (NUMRES)   - ...
C   H_FONDO             (NUMRES)   - ...
C
C Input and Output   (defined in common block include file IOUNITS.H
C File Units          and in BLOCK DATA subprogram block_data.f)
C ----------------
C   IFN    - I/O file names
C            (see subroutine OPENIO for unit IFN input)
C   IIN1   - control parameters for FLOW3D and SURF_ROUTE
C   IIN2   - grid info
C   IIN3   - root map (ZROOT)
C   IIN4   - soil characteristics
C   IIN5   - initial conditions
C   IIN6   - atmospheric BC's (rainfall/evaporation rates)
C            (see subroutines ATMONE and ATMNXT for unit IIN6 input)
C   IIN7   - seepage face BC's
C   IIN8   - non-atmospheric, non-seepage face Dirichlet BC's
C            (see subroutines BCONE and BCNXT for unit IIN8 input)
C   IIN9   - non-atmospheric, non-seepage face Neumann BC's
C            (see subroutines BCONE and BCNXT for unit IIN9 input)
C   IIN10  - DEM information (dem.dat)
C   IIN11  - DEM parameters and other parameters (if GRID=TRUE they can
C            be found in the grid file)
C   IIN16  - PCAP,SATC,KRWC values for moisture curve lookup table
C            option
C   IIN17  - position of reservoirs and buffer cells
C            (posizione_serb.dat)
C   IIN18  - initial levels in reservoirs (livelli_iniz_serb.dat)
C   IIN19  - depitting parameter epsilon (depit.dat)
C   IIN20  - map of lake cells (to be excluded by the domain)
C            (lakes_map.dat)
C   IIN21  - map of zone cells 
C   IIN22  - effective rainfall input file (this file is read only
C            in the case of surface simulation)
C   IIN23  - vector containing the index of the cells ordered into
C            descending elevation value
C   IIN25  - raster of the weights - cardinal direction 
C   IIN26  - raster of the weights - diagonal direction
C   IIN27  - raster of the drainage directions - cardinal direction 
C   IIN28  - raster of the drainage directions - diagonal direction
C   IIN29  - raster of the local slopes - cardinal direction
C   IIN30  - raster of the local slopes - diagonal direction
C   IIN31  - raster of the elemental path length - cardinal direction 
C   IIN32  - raster of the elemental path length - diagonal direction
C   IIN33  - raster of surface roughness coefficient (ks) - cardinal direction
C   IIN34  - raster of surface roughness coefficient (ks) - diagonal direction
C   IIN35  - raster of surface water width - diagonal direction
C   IIN36  - raster of surface water width - cardinal direction
C   IIN37  - raster of at-a-station scaling coefficients for the 
C            surface water width 
C   IIN38  - raster of at-a-station scaling coefficients for the 
C            surface roughness (ks) coefficient
C   IIN39  - raster containing the number of rivulets per cell
C   IIN40  - EnKF input
C   IIN50  - nudging input
C            (see subroutines NUDONE and NUDNXT for unit IIN50 input)
C   IIN51  - mesh input file
C   IIN60  - raster of basement elevation
C   IIN61  - transp input file, contains transpot numerical and physical
C            parameters
C   IIN62  - transp_ic input file, contains the parameters to define 
C            the initial concentration field  
C   IIN63  - transp_atmbc input file, contains, the atmospheric forcing 
C            for transport
C   IIN64  - transp_dirbc input files, contains the parameter to impose 
C            dirichlet boundary conditions for transport 
C   ITERM  - set ITERM to 6 in BLOCK DATA subprogram for terminal
C            output; otherwise unit ITERM output is to a file
C   IOUT1  - debugging
C   IOUT2  - main output
C   IOUT3  - X, Y, Z coordinate values
C   IOUT4  - convergence behavior and error norms for each iteration of
C            every time step
C   IOUT5  - mass balance and convergence behavior at each time step
C   IOUT6  - vertical profile output
C   IOUT7  - atmospheric and seepage face hydrograph output
C   IOUT8  - non-atmospheric, non-seepage face hydrograph output
C   IOUT9  - detailed HGFLAG output
C   IOUT10 - detailed SFFLAG output
C   IOUT11 - pressure head output at all nodes
C   IOUT12 - velocity output at all nodes
C   IOUT13 - water saturation output at all nodes (SW) for input to
C            TRAN3D and DUAL3D groundwater contaminant transport codes
C   IOUT14 - CKRW output at all nodes
C   IOUT15 - velocity output at all elements, for input to TRAN3D and
C            DUAL3D groundwater contaminant transport codes
C   IOUT16 - pressure head output at surface nodes
C   IOUT17 - SATSUR output at surface nodes
C   IOUT18 - SW output at surface nodes 
C   IOUT19 - non-atmospheric, non-seepage face Dirichlet BCs at each
C            time step
C   IOUT20 - non-atmospheric, non-seepage face Neumann BCs at each
C            time step
C   IOUT30 - detailed seepage face hydrograph output
C   IOUT31 - detailed non-atmospheric, non-seepage face Dirichlet
C            hydrograph output
C   IOUT32 - detailed non-atmospheric, non-seepage face Neumann
C            hydrograph output
C   IOUT36 - output of cumulative flow volumes VSFTOT, VNDTOT, VNNTOT,
C            VNUDTOT, and VTOT
C   IOUT40 - SURF_ROUTE module input data (dem data, geometry data, etc.)
C   IOUT41 - surface runoff hydrograph: plot the computed discharge at 
C            the outlet cell of the basin vs. time (each time value, not
C            only TIMPRT)
C   IOUT42 - SURF_ROUTE ponding head output
C   IOUT43 - CPU, time stepping, iteration and other diagnostics of the
C            surface and subsurface modules at each time step
C   IOUT44 - detailed recharge output
C   IOUT50 - detailed nudging "hydrograph" output
C   IOUT51 - detailed time series output of model results at the
C            nudging observation points. For intercomparison with a
C            model simulation without nudging, run the same simulation
C            but with NUDG=0.0 or NUDEPS=0.0 (don't set NUDN=0 since
C            NUDSMC cannot be calculated without the coordinates of the
C            nudging observation points!), and plot the results in
C            the IOUT51 output file from both runs, together with the
C            NUDTIM and NUDVAL data from the nudging input file.
C            (This output file is designed for NUDN <= 10; for NUDN > 10
C            the output will need re-structuring.)
C   IOUT52 - detailed time series output of the ensemble of PNEW 
C            realizations before the update.
C   IOUT53 - detailed time series output of the ensemble of outlet Q_OUT 
C            realizations after the update.
C   IOUT54 - detailed time series output of the ensemble of PNEW 
C            realizations after the update.
C   IOUT55 - detailed time series output of the ensemble of subsurface
C            water volume.
C   IOUT56 - ensemble parameters, initial conditions, weights and 
c            SIR updates
C   IOUTPT - void ratio output at all nodes in case of deformable peat 
C
C----------------------------------------DATA AND PARAMETER DECLARATIONS
C
C  declare all variables!
C
      IMPLICIT NONE
C
C  dimensioning parameters
C
      INCLUDE 'CATHY.H'

C  common block for surface transport
C
      INCLUDE 'TRANSPSURF.H'
C  
C  actual or minimum dimensions
C
      INTEGER NROW,NCOL,NCELL,NUMRES
      INTEGER NNOD,NTRI,NSTR,N,NT,NSF,NUMDIR
      INTEGER NDIR(3),NDIRC(3),NNEU(3),NNEUC(3),NP(3),NQ(3)
      INTEGER NUDN,NUDT,NUDC
      INTEGER NZONE,NVEG,NLKP,NTERM,ITUNS
      INTEGER NR,NPRT,NUMVP,NUM_QOUT
      INTEGER N1,IBOT
C
C  general and FLOW3D integer parameters, flags, and indices
C
      INTEGER KPRT,IPEAT,IPRT,IPRT1,ISIMGR,IVERT,ISP,INDP,VTKF
      INTEGER ANP,ANQ,HTIDIR,HTINEU,HSPATM,HTIATM,IETO
      INTEGER LUMP,IOPT,NLRELX,L2NORM,ISOLV,IVGHU,KSLOPE
      INTEGER ISFONE,ISFCVG,KSFCV,KSFCVT,DUPUIT
      INTEGER ITER,ITUNS1,ITUNS2,NSTEP,KBACKT,KBACK,ITRTOT
      INTEGER ITMXCG,NITER,NITERT,ITLIN,KLSFAI,IMAX,MAXITER
      INTEGER MINBOT,NDZ,NUDFLAG,WFLAG
      INTEGER NCOUT
      INTEGER VELREC
C
C  integer parameters for SURF_ROUTE 
C
      INTEGER HSPEFF,HTIEFF
      INTEGER NSURF,NSURFT,NSURFT_T,NSURFT_TB
      INTEGER NCELNL,IPOND,DOSTEP,NCELL_COARSE
C
C  integer parameters for nudging and EnKF
C
      INTEGER NUDCTR,NOBS
CM    INTEGER NEFFMIN,NSTEP1,ITRTOT1,NERT,NROUT,CONTDS
C integer parameteres for transport
      INTEGER TRAFLAG,FLAG_TEMP,FLAG_INTERP,FLAG_LIMITER
      INTEGER NADV,NADVFL,LUMPC,CFLTETRA,TADV,CDIFFUS
      INTEGER NTRI3,IPRTCG,IPREC,ISOL,IEXIT
      INTEGER NFACE,NRAUX,NIAUX
      INTEGER NDIR_TRA(3),NDIRC_TRA(3),NP_TRA(3)
      INTEGER NPFA_TRA,HTIDIR_TRA
      INTEGER ANP_TRA,HSPVEL,NQ_TRA
      INTEGER CAUSPATM,CAUTIATM,IETOCAU
      INTEGER IFCONC(NODMAX)
C
C  miscellaneous integer scalars
C
      INTEGER I,J,K,NRE,NRE1
      INTEGER IROW,ICOL
      INTEGER iface
C
C  logical flags
C
      LOGICAL FL3D,SURF,DEM,GRID,PONDING,PONDP
      LOGICAL DTGMIN,LSFAIL,ERRGMX,NORMCV,ITAGEN,SFCHEK,KSFZER
      LOGICAL NOBACK,BCKSTP,ENKF,UPD,RESAMP,TRANSP,REACFLAG
C
C  real scalars for initial conditions and vertical discretization 
C
      REAL*8  WTPOSITION,BASE
C
C  real scalars for time stepping and linear and nonlinear iterations
C
      REAL*8  TETAF,DELTAT,DELTAT0,DTMIN,DTMAX,DTAVG,DTSMAL,DTBIG
      REAL*8  TSMAL,TBIG
      REAL*8  DTMAGA,DTMAGM,DTREDS,DTREDM,TMAX,TIME,TIMEP
      REAL*8  OMEGA,OMEGAP
      REAL*8  TOLUNS,TOLSWI,ERNLMX,TOLCG
C
C real scalars for not-detailed output times
C
CM    REAL*8  DTOUT,TIMEPOUT
C
C  real parameters for SURF_ROUTE 
C 
      REAL*8 DELTATS
C
C real scalars for ponding
C 
      REAL*8  PONDH_MIN
C
C real scalars for nudging and EnKF
C 
      REAL*8  NUDG
CM    REAL*8  DSRETC,DSKS,DSSTOR,DSPOROS,DSIC,DSATM,DSSURF
CM    REAL*8  ENDSSURF_KS,ENDSSURF_WS,DSKSZ
CM    REAL*8  QTIMEP_1,QNEW_1,WSUM,STORE_SAT
CM    REAL*8  QEN_SPREAD,DSMEASMAX
C
C Real scalars for transport modelling      
C
      REAL*8 CFLINP,CFLNUMB,TETAC,DIFFUS,GAMMAS
      REAL*8 DELTATADV,SUPDIFFUS,PECNUMB
      REAL*8 TIMEAUX,TIMWGT,MIXPART   
C
C  miscellaneous real scalars
C
      REAL*8  AREATOT,VOLTOT,SCRATCH
      REAL*8  FHORT,FDUNN,FPOND,FSAT
      REAL*8  VOLUME_OUT,VOLUME_SUP,VOLSUPNEW
      REAL*8  EVAP_EFF,INF_TOT,QTRAN
      REAL*8  RMAX,RMIN
      REAL*8  RECFLOW,RECVOL
C
C  real*4 scalars for cpu timing
C
      REAL    CPUSUB,CPUSUB_T,CPUSURF,CPUSURF_T
      REAL    CPUNL,CPUOVH,CPUMN,CPUNLT,CPUUPD,CPUUPD1,CPUUPD2
      REAL    AVGNL,AVGLIN,AVGLNL,ATCTS,ATCNL,ANCTS,ANCNL,ATCTS1
      REAL    PCMN,PCNL,PCOVH,PCNLT
      REAL    PCVEC1,PCVEC2,PCVEC3,PCVEC4,PCVEC5,PCVEC6,PCVEC7
      REAL    PCVEC8,PCVEC9
C
C  real*4 array for cpu timing
C
      REAL    CPUVEC(11)
C CARLOTTA QUI CACELLA POI                                                                            
C CARLOTTA QUI CACELLA POI 
      REAL*8   PRECL(NFACEMAX),PRECR(NFACEMAX)
      REAL*8   SFMASSOUTTOT,SFMASSOUT
      REAL*8   SFEV,SFEVTOT,SFIN,SFINTOT,DISPTOT
      INTEGER  Node1,Node2,Node3
      INTEGER  FACCIA,ISTR
C
C  integer arrays for grid, BCs, outputs, and system matrices
C
      INTEGER TP(NMAX),TP2D(NODMAX)
      INTEGER TRIANG(4,NTRMAX),CELL(5,MAXCEL)
      INTEGER TETRA(5,NTEMAX),IVOL(NTEMAX)
      INTEGER ZERO(3),ZEROC(3)
      INTEGER CONTP2(NP2MAX),CONTP(3,NPMAX),CONTQ(3,NQMAX)
      INTEGER ACONTP(NPMAX),ACONTQ(NQMAX)
      INTEGER PUNTDIRFLOW_NODE(NMAX),PUNTNEUFLOW_NODE(NMAX)
      INTEGER IFATM(NODMAX),IFATMP(NODMAX),SATSUR(NODMAX)
      INTEGER NSFNUM(NSFMAX),NSFNOD(NSFMAX,NNSFMX)
      INTEGER IFSF(NMAX),IFSFP(NMAX),PUNTDIRSF_NODE(NMAX)
      INTEGER IBNDRY(NFACEMAX),IBNDRYP(NFACEMAX)
      INTEGER SFEX(NSFMAX,NNSFMX),SFEXIT(NSFMAX,NNSFMX)
      INTEGER SFEXP(NSFMAX,NNSFMX)
      INTEGER NODDIR(MAXDIR),CONTR(NRMAX),NODVP(MAXVP)
      INTEGER IA(MAXTRM),JA(NTPMAX),TOPOL(NMAX+1),TETJA(4,4,NTEMAX)
      INTEGER IP3(3,3),IP4(4,4)
      INTEGER SFFLAG(5),HGFLAG(9),IER(7),INSYM(INTBOT)
      INTEGER STR_ZON(MAXSTR),VEG_TYPE(NODMAX)
C
C  integer arrays for SURF_ROUTE
C      
      INTEGER ID_QOUT(MAXQOUT)
      INTEGER ZONE(ROWMAX,COLMAX)
      INTEGER LAKES_MAP(ROWMAX,COLMAX)
      INTEGER NODI(ROWMAX+1,COLMAX+1)
      INTEGER INDCEL(ROWMAX,COLMAX),INDCELWL(ROWMAX,COLMAX)
      INTEGER CELTYPE(MAXCEL)
      INTEGER CELLCOL(MAXCEL),CELLROW(MAXCEL)
      INTEGER CELLCOL_WL(MAXCEL),CELLROW_WL(MAXCEL)
      INTEGER CELLS_R(MAXRES)
c     INTEGER MODIF(MAXCEL)
      INTEGER RESERVR(MAXCEL),N_HA(MAXRES)
      INTEGER TIPO_R(MAXCEL)
C
C  real arrays for SURF_ROUTE
C
      REAL*8 DEM_MAP(ROWMAX,COLMAX),BASE_MAP(ROWMAX,COLMAX)
      REAL*8 ROOT_MAP(ROWMAX,COLMAX)
      REAL*8 EFFTIM(3),SURFACE_WATER_INP(2,MAXCEL)
C
C  integer arrays for subsurface transport
C
      INTEGER   IPARM(10),IAC(MAXTRM)
      INTEGER   SIDE_CNC(4,NTEMAX),PLIST(2,NFACEMAX),ISIDE(3,NFACEMAX)
      INTEGER   NEIGH(4,NTEMAX),TETRAVERT(NMAX,N1MAX),TVERT(NMAX)
      INTEGER   NODE_FACECOUNT(NMAX),NODE_FACEIDS(NMAX,MAXFCONTONODE) 
      INTEGER   Node_ElementIDs(NMAX,N1MAX)
      INTEGER   Node_ElementCount(NMAX)
      INTEGER   IAUX(NIAUXMAX)
      INTEGER   CONTP_TRA(3,NPMAX_TRA)
      INTEGER   CONTPFA_TRA(NFACEMAX),PUNTDIR(NFACEMAX),NEIGBC(3,NPMAX)
      INTEGER   ACONTP_TRA(NPMAX_TRA)
      INTEGER   CONTQ_TRA(NODMAX)
      INTEGER   Faces_FaceType(NFACEMAX)

C  integer arrays for nudging and EnKF
C
      INTEGER NUDTET(MAXNUDN)
C
C  real arrays for mesh configuration and boundary conditions
C
      REAL*8  X(NMAX),Y(NMAX),Z(NMAX)
      REAL*8  XC(NTEMAX),YC(NTEMAX),ZC(NTEMAX)
      REAL*8  VOLNOD(NMAX),VOLU(NTEMAX),VOLUR(NTEMAX),SURFAR(NTEMAX)
      REAL*8  ZRATIO(MAXSTR),DEPTH(NODMAX),ELTRIA(NTRMAX)
      REAL*8  PRESC(NPMAX),PTIM(3),PINP(3,NPMAX)
      REAL*8  Q(NQMAX),QTIM(3),QINP(3,NQMAX)
      REAL*8  QPNEW(NPMAX),QPOLD(NPMAX)
      REAL*8  QTRANIE(NMAX)
      REAL*8  SFQ(NSFMAX,NNSFMX),SFQP(NSFMAX,NNSFMX)
      REAL*8  SFQ_PARTIAL(NSFMAX)
      REAL*8  ARENOD(NODMAX),ATMPOT(NODMAX),ATMACT(NODMAX)
      REAL*8  ATMOLD(NODMAX),ATMTIM(3),ATMINP(3,NODMAX)
      REAL*8  CONNEU_NODE(NMAX),CONDIR_NODE(NMAX) 
C
C  exchange variables between FLOW3D and SURF_ROUTE
C
      REAL*8  OVFLNOD(NODMAX),OVFLP(NODMAX),PONDNOD(NODMAX)
      REAL*8  OVFLCEL(MAXCEL),PONDCEL(MAXCEL),pondnodp(nodmax)
c     REAL*8  ACTCEL(MAXCEL)

C
C  Real Arrays for subsurface transport modeling and velocity reconstruction
C
      REAL*8  AREFACE(NFACEMAX),FACES_REFERENCEVECTOR(3,NFACEMAX)
      REAL*8  FaceCentroid(3,NFACEMAX)
      REAL*8  Faces_NeumannFlux(NFACEMAX)
      REAL*8  Faces_FaceFlux(3,NFACEMAX),Faces_FaceFluxOld(3,NFACEMAX)
      REAL*8  FaceFlux(3,NFACEMAX)
C
C  exchange variables between FLOW3D and SURF_ROUTE
C
      REAL*8  TRAFLNOD(NODMAX),TRAFLP(NODMAX),CONCNOD(NODMAX)
      REAL*8  TRAFLNOD_FLOW(NODMAX),CONCNOD_OLD(NODMAX)
      REAL*8  TRAFLCEL(MAXCEL),CONCCEL(MAXCEL),CONCNODUPD(NODMAX)
      REAL*8  CONCSURFVTK(NMAX),CONCNODBC(NODMAX)
      REAL*8  SOURCE_MIXING(NODMAX),MIX_CORRECT(MAXCEL)
      REAL*8  SURFACE_MIX_SN(MAXCEL)
C
C  real arrays for nudging
C
      REAL*8  NUDTIM(MAXNUDT),NUDVAL(MAXNUDC,MAXNUDN),NUDSMC(MAXNUDN)
      REAL*8  NUDTAU(MAXNUDT,MAXNUDN),NUDRXY(MAXNUDN),NUDRZ(MAXNUDN)
      REAL*8  NUDEPS(MAXNUDN),NUDX(MAXNUDN),NUDY(MAXNUDN),NUDZ(MAXNUDN)
      REAL*8  NUDDIF(MAXNUDN),NUDNOD(NMAX),NUDCUM(NMAX)
C
C  real arrays for material properties and hydraulic characteristics
C
      REAL*8  KS(NTEMAX),KZ(NTEMAX),KZNOD(NMAX)
      REAL*8  PERMX(MAXSTR,MAXZON),PERMY(MAXSTR,MAXZON)
      REAL*8  PERMZ(MAXSTR,MAXZON),ELSTOR(MAXSTR,MAXZON)
      REAL*8  POROS(MAXSTR,MAXZON),VGNCELL(MAXSTR,MAXZON)
      REAL*8  VGRMCCELL(MAXSTR,MAXZON),VGPSATCELL(MAXSTR,MAXZON)
      REAL*8  SNODI(NMAX),PNODI(NMAX)
      REAL*8  PORE(NMAX),INDE(NMAX),INDE0(NMAX),DEF(NMAX)
      REAL*8  ETAI(NMAX),CKRW(NMAX),CKRWP(NMAX),SW(NMAX)
      REAL*8  SENODI(NMAX),DCKRW(NMAX),DETAI(NMAX)
      REAL*8  ETAE(NTEMAX),CKRWE(NTEMAX),SWE(NTEMAX),SEELT(NTEMAX)
      REAL*8  DCKRWE(NTEMAX),DETAIE(NTEMAX)
      REAL*8  ARCHIE(4),ARCHIE_R(2)
      REAL*8  PEL(NTEMAX),SEL(NTEMAX),SWP(NMAX)
      REAL*8  SWINT_EL(NTEMAX),SWPINT_EL(NTEMAX),SW_SAV(NMAX)
      REAL*8  LEL(NTEMAX),KEL(NTEMAX)
      REAL*8  ATMINT(NODMAX)
c CARLOTTA PEL (element porosity) ET1 ET1E(saturation*elastic storage 
c          at nodes and elements) ET2 (dSw/dpnew at nodes)
      REAL*8  ET1(NMAX),ET2(NMAX),ET1E(NTEMAX)
c CARLOTTA
C
C real arrays for subsurface transport
C     
      REAL*8 CNEW(NTEMAX),COLD(NTEMAX),COLDOLD(NTEMAX),CNEW_SAV(NTEMAX)
      REAL*8 CNNEW(NMAX),CNOLD(NMAX),CNEW_VTK(NTEMAX),CSOL(NTEMAX)
      REAL*8 CDIFF(NTEMAX),FLUX(NTEMAX)
      REAL*8 UUOLD(NTEMAX),VVOLD(NTEMAX),WWOLD(NTEMAX)
      REAL*8 UUINT(NTEMAX),VVINT(NTEMAX),WWINT(NTEMAX)
      REAL*8 SWNEW(NTEMAX),SWOLD(NTEMAX),SWINT(NTEMAX)
      REAL*8 DISNOD(24,NMAX),WEIGHTNOD(NMAX)        
      REAL*8 ALFAL(MAXSTR,MAXZON),ALFAT(MAXSTR,MAXZON)
      REAL*8 KD(MAXSTR,MAXZON),LAMBDA(MAXSTR,MAXZON)
      REAL*8 XFAC(3),YFAC(3),ZFAC(3),HERON
      REAL*8 LAMBDANODI(NMAX),KDNODI(NMAX)
      REAL*8 LMASSC(4,4),PRESCFA_TRA(NFACEMAX)
      REAL*8 PRESC_TRA(NPMAX_TRA),CTIM(3),CINP(3,NPMAX_TRA)
      REAL*8 CONCNOD_PERDIR(NODMAX),CONCNOD_PERDIR_TIMEP(NODMAX)
      REAL*8 CONCNOD_PERDIR_INT(NODMAX)
      REAL*8 VN(NFACEMAX),AUX(NRAUXMAX)
      REAL*8 COEF1C(MAXTRM),COEF2C(MAXTRM)
      REAL*8 LHSC(NPMAX_TRA),QCAUCHY(NODMAX)
      REAL*8 XT5C(NMAX),TNOTIC(NMAX)
      REAL*8 ATMCONC(NODMAX),CONCOLD(NODMAX)
      REAL*8 CONCTIM(3),CONCINP(3,NODMAX)
      REAL*8 MASS0,MASST,MASSINADV,MASSOUTADV,MASSINTOT,MASSOUTTOT
      REAL*8 MASSINTOTADV,MASSOUTTOTADV,MASSINTOTCUM,MASSOUTTOTCUM
      REAL*8 MERRAS,CMERRAS,DSMASS,MASSSURF,MASSZONE(MAXZON)
      REAL*8 SFMASS,MASSSOLZONE(MAXZON),VOLZONE(MAXZON)
C  real arrays for pressure heads, velocities, and outputs
C
      REAL*8  PNEW(NMAX),POLD(NMAX),PDIFF(NMAX),PTNEW(NMAX)
      REAL*8  PTOLD(NMAX),PTIMEP(NMAX)
      REAL*8  UNOD(NMAX),VNOD(NMAX),WNOD(NMAX),RECNOD(NODMAX)
      REAL*8  UU(NTEMAX),VV(NTEMAX),WW(NTEMAX),TIMPRT(MAXPRT)
C
C  real arrays for system components
C
      REAL*8  AI(4,NTEMAX),BI(4,NTEMAX),CI(4,NTEMAX),DI(4,NTEMAX)
      REAL*8  LMASS(4,4)
      REAL*8  COEF1(MAXTRM),COEF2(MAXTRM),COEF3(MAXTRM),SCR(MAXTRM)
c CARLOTTA COEF4 
      REAL*8  COEF4(MAXTRM)
c CARLOTTA      
      REAL*8  RNSYM(MAXBOT),RWU(NRMAX)
      REAL*8  TNOTI(NMAX),XT5(NMAX)
      REAL*8  LHSP(NPMAX),LHSSF(NSFMAX,NNSFMX),LHSATM(NODMAX)
      REAL*8  BKCFLOW_NODE(NMAX),BKCFLOWP_NODE(NMAX)
      REAL*8  BKCFLOWINT_NODE(NMAX)
C
C  new variables for time-varying seepage face boundary conditions
C
      INTEGER SFV(2),SFVNUM(2,NSFMAX),SFVNOD(2,NSFMAX,NNSFMX)
      REAL*8  SFVTIM(2)
C
C  real scalars for mass balance and hydrograph calculations
C
      INCLUDE 'MB_HGRAPH.H'
C
C  real scalars and real arrays for unsaturated soil characteristics
C  and for chord and tangent slope formulas
C
      INCLUDE 'SOILCHAR.H'
C
C  real scalars and real arrays for surface water routing 
C
      INCLUDE 'SURFWATER.H'
C
C   integer and real scalars and arrays for surface river network and
C   watershed characteristics
C
      INCLUDE 'RIVERNETWORK.H'
C
C  input and ouput file units
C
      INCLUDE 'IOUNITS.H'
C
C  data declarations
C
      DATA    IMAX/2147483647/
      DATA    RMAX/1.7D+100/
      DATA    RMIN/1.7D-100/
      DATA    IP3/1,2,3,2,3,1,3,1,2/
      DATA    IP4/1,2,3,4,2,3,4,1,3,4,1,2,4,1,2,3/
C
C  open the I/O files
C
      CALL OPENIO
C
C------------------------------------------DATA INPUT AND INITIALIZATION
C
C  initialization call to the cpu timer subroutine.
C  CPU timer subroutine TIM is in program 'tim.f' (for UNIX) and
C  in program 'timcray.f' (for CRAY UNICOS).
C
      CALL TIM(CPUMN,1)
C
C  start reading input data
C
c    
      CALL DATIN(ISIMGR,IVERT,ISP,WTPOSITION,BASE,ZRATIO, 
     1           TRIANG,X,Y,Z,PERMX,PERMY,PERMZ,ELSTOR,
     2           POROS,VGNCELL,VGRMCCELL,VGPSATCELL,CONTR,NODVP,ID_QOUT,
     3           TIMPRT,PONDNOD,PTIMEP,TETAF,
     4           DELTAT,DTMIN,DTMAX,TMAX,DTMAGA,DTMAGM,
     5           DTREDS,DTREDM,ITUNS,ITUNS1,ITUNS2,
     6           TOLUNS,TOLSWI,ERNLMX,ITMXCG,TOLCG,
     7           KSLOPE,LUMP,IPEAT,IVGHU,NLKP,VTKF,
     8           IOPT,ISOLV,IPRT1,IPRT,IPOND,INDP,NNOD,NTRI,NSTR,
     9           NZONE,NVEG,N1,NR,NUMVP,NUM_QOUT,NPRT,N,NT,
     A           ISFONE,ISFCVG,DUPUIT,
     B           L2NORM,NLRELX,OMEGA,
     C           PONDH_MIN,
     D           FL3D,SURF,DEM,GRID,NROW,NCOL,
     E           NCELNL,NCELL,DOSTEP,NCELL_COARSE,
     F           DEM_MAP,ZONE,LAKES_MAP,INDCEL,INDCELWL,
     G           CELL,CELLCOL,CELLROW,TP2D,NODI,
     H           TIPO_R,RESERVR,N_HA,CELLCOL_WL,CELLROW_WL,
     I           BASE_MAP,DEPTH,ELTRIA,NCOUT,TRAFLAG,TRANSP,
     K           VELREC,VEG_TYPE)

C Read transport input files if traflag = 1 
      IF (TRANSP) THEN 
         CALL DATIN_TRA(NSTR,NZONE,NTRI,NTRI3,NADV,NADVFL,LUMPC,
     1        CDIFFUS,FLAG_TEMP,FLAG_INTERP,FLAG_LIMITER,CFLINP,TETAC,
     2        DIFFUS,GAMMAS,IPARM,ALFAL,ALFAT,KD,LAMBDA,MIXPART,
     3        REACFLAG)
C         
         WRITE(IOUT2,*) ' -------------------------'
         WRITE(IOUT2,*) ' SIMULATION WITH TRANSPORT'
         WRITE(IOUT2,*) ' -------------------------'
         WRITE(IOUT2,*) '------ Longitudinal Dispersion-------'
         WRITE(IOUT2,*) (ALFAL(i,1),i=1,nstr)
         write(IOUT2,*) '--------Transverse Dispersion--------'
         WRITE(IOUT2,*) (ALFAT(i,1),i=1,nstr)
         write(IOUT2,*) '--------- Diffusion -----------------'
         WRITE(IOUT2,*) DIFFUS
      END IF      
C
C  Inizialization of NUMRES. This command must be removed when the lake
C  procedure will be active! MC
C
      NUMRES=0
      CALL VCOPYR(NUMRES,H_POOL_KK_VEC,H_POOL_KKP1_VEC)    
      CALL VCOPYR(NUMRES,H_POOL_KK_VEC_SAV,H_POOL_KK_VEC)    
C
C  further processing of the 2-D mesh (sorting and area calculation),
C  generation of the 3-D grid, and setup of system matrices
C
     
      IF (FL3D) THEN     
         CALL GRDSYS(IPRT1,IOPT,N,NT,NNOD,NTRI,NTERM,N1,NDZ,
     1        NSTR,IVERT,LUMP,IMAX,
     2        TRIANG,TETRA,IA,JA,TOPOL,TETJA,IVOL,
     3        IP3,IP4,IER,
     4        RMAX,BASE,DEPTH,ZRATIO,ARENOD,X,Y,Z,XC,YC,ZC,
     5        AI,BI,CI,DI,VOLNOD,VOLU,VOLUR,LMASS)
C     
      END IF
C     
C Connection needed for finite volume - advective transport
      IF (TRANSP) THEN 
         CALL CONNECTION_NEW(N,NT,NFACE,TETRA,SIDE_CNC,PLIST,ISIDE,
     1        NEIGH,AREFACE,FaceCentroid,Node_ElementCount,
     2        Node_ElementIDs,NODE_FACECOUNT,NODE_FACEIDS,
     3        Faces_ReferenceVector,X,Y,Z)
      WRITE(IOUT2,*)'NFACE =',NFACE
C     
c     C     calculation of SURFAR
         DO I=1,NT
            DO J=1,4
               IFACE=SIDE_CNC(J,I)
               DO K=1,3
                  XFAC(K) = X(ISIDE(K,IFACE))
                  YFAC(K) = Y(ISIDE(K,IFACE))
                  ZFAC(K) = Z(ISIDE(K,IFACE))
               END DO
               SURFAR(I)=SURFAR(I) + HERON(XFAC,YFAC,ZFAC)
            END DO
         END DO
C     
C     
         CALL WEIGHTNEIGHNODE(N,NT,X,Y,Z,XC,YC,ZC,TETRA,DISNOD,
     1        WEIGHTNOD)
         CALL LOCMAS(LMASSC,LUMPC)
      END IF
C
C IT IS OUT OF GRDSYS NOW BECAUSE OF THE CONNECTION PROCEDURE !!!!!!
C  set up pointers and indices for storage of the system matrices
C
      IF (IOPT .EQ. 1) THEN
         CALL STRPIC(N,NTERM,TETRA,JA,TOPOL,NT,N1,IMAX)
         CALL CHKPIC(N,NTERM,TOPOL,JA,NDZ,IER)
         CALL TETPIC(NT,TETRA,JA,TOPOL,TETJA)
      ELSE
         CALL STRNEW(N,NTERM,TETRA,JA,TOPOL,NT,N1,IMAX)
         CALL CHKNEW(N,NTERM,TOPOL,JA,NDZ,IER)
         CALL TETNEW(NT,TETRA,JA,TOPOL,TETJA)
         CALL TOPIA(N,TOPOL,IA)
      END IF
C      
      write(iout2,1005) NTERM
cxcx  NRAUX = 5*NFACE + NTERM
cxcx  NIAUX = NTERM + NFACE + 1
      NRAUX = 5*N + NTERM
      NIAUX = NTERM + N + 1
C
C     further input and initialization for ICs, BCs, nudging, parameters,
C  and misc counters, flags, and arrays
C
      IF (.NOT. FL3D) THEN
         CALL INIT_SURF(NCELNL,NUMRES,DELTA_X,DELTA_Y,HSPEFF,
     1        HTIEFF,NSTEP,TIME,DELTAT,DELTATS,QOI_SN,EFFTIM,
     2        SURFACE_WATER_INP)
      ELSE
         CALL INITAL(NNOD,NSTR,N,NT,NTRI,NP,NQ,NSF,NDIR,NDIRC,
     1        NNEU,NNEUC,HTIDIR,HTINEU,HSPATM,HTIATM,IETO,IPRT1,
     2        IPOND,INDP,NITERT,ITLIN,ITRTOT,ITER,NSTEP,
     3        NSURFT,NSURFT_T,NSURFT_TB,
     4        KPRT,ISFCVG,KSFCV,KSFCVT,KBACKT,KBACK,KLSFAI,
     5        IOPT,IPEAT,IVGHU,KSLOPE,HGFLAG,SFFLAG,
     6        TETRA,TP,NODDIR,CONTP,CONTQ,IFATM,IFATMP,
     7        NSFNUM,NSFNOD,SFEX,SFEXP,SFEXIT,DUPUIT,
     8        NUDN,NUDCTR,NUDT,NUDC,NUDG,NUDFLAG,WFLAG,
     9        NUDTIM,NUDRXY,NUDRZ,NUDTAU,
     A        NUDX,NUDY,NUDZ,NUDEPS,NUDTET,NUDCUM,
     B        DTGMIN,SFCHEK,SURF,DEM,PONDING,PONDP,
     C        WTPOSITION,DELTAT,TMAX,TETAF,DTMIN,DTMAX, 
     D        DTSMAL,DTBIG,DTAVG,TIME,TIMEP,PONDH_MIN,
     E        AREATOT,VOLTOT,RMAX,PEL,
     F        X,Y,Z,XC,YC,ZC,VOLU,ARENOD,POROS,ELSTOR,
     H        VGNCELL,VGRMCCELL,VGPSATCELL,SNODI,PNODI,
     I        DEF,PORE,INDE,PNEW,PTIMEP,POLD,PTNEW,PTOLD,
     I        SFQP,Q,QPOLD,QTIM,QINP,PRESC,PTIM,PINP,
     J        ATMPOT,ATMACT,ATMOLD,ATMTIM,ATMINP,
     K        OVFLNOD,OVFLP,PONDNOD,
     L        NCELL,NCELNL,NUMRES,
     M        CPUSUB,CPUSUB_T,CPUSURF,CPUSURF_T,CPUNL,CPUVEC,
     N        ANP,ANQ,ACONTP,ACONTQ,
     O        SFV,SFVNUM,SFVNOD,SFVTIM,QTRANIE,
     P        IFSF,IFSFP,PUNTDIRSF_NODE,SEL,SW,SWP,
     R        PUNTDIRFLOW_NODE,PUNTNEUFLOW_NODE,
     S        CONDIR_NODE,CONNEU_NODE,
     T        LEL,KEL,LAMBDA,LAMBDANODI,KD,KDNODI,
     V        TRIANG)
C  
C  calculate SW and CKRW needed for storage and velocity calculations
C  and detailed output
C
         IF (IPEAT .EQ. 0) THEN
            CALL CHVELO(NLKP,N,NT,NTRI,IVGHU,TETRA,TP,
     1                  PTIMEP,SNODI,PNODI,SW,CKRW,SWE,CKRWE)
            CALL STORCAL(N,STORE1,STORE2,0.0D0,SW,PNODI,VOLNOD)
         ELSE
CM
CM There used to be a CHVELOP0 here for the case of peat soils, with
CM a function CORNR. Need to double-check whether it is not needed
CM
            CALL CHVELOP(NLKP,NNOD,N,NT,NTRI,IVGHU,TETRA,TP,
     1                   PTIMEP,Z,PORE,SW,CKRW,SENODI,CKRWE,SEELT)
            CALL STORCAL(N,STORE1,STORE2,0.0D0,SW,PORE,VOLNOD)
            DO I=1,N
               IF (SW(I).EQ.1.0D0) THEN
                  INDE0(I)=(1+(PNODI(I)/(1-PNODI(I))))*EXP(SNODI(I)
     1                     *PTIMEP(I))-1
               END IF
            END DO
         END IF
         STORE0 = STORE1
         STORE2 = STORE1
         IF (IPRT .GE. 2) THEN
            CALL NODELT(NT,TETRA,CKRW,CKRWE)
            CALL VEL3D(NT,TETRA,NTRI,PTIMEP,PERMX,PERMY,PERMZ,
     1                 UU,VV,WW,BI,CI,DI,CKRWE,VOLUR,IVOL)
            CALL VNOD3D(N,NT,TP,TETRA,UU,VV,WW,UNOD,VNOD,WNOD)
         END IF
C Free drainage
         CALL INIT0R(NTEMAX,KS)
         CALL INIT0R(NTEMAX,KZ)
         CALL INIT0R(NMAX,KZNOD)
         DO I = 1,NSTR
            DO J=(I-1)*NTRI*3+1,I*NTRI*3
               KS(J)=PERMX(I,TETRA(5,J))
               KZ(J)=PERMZ(I,TETRA(5,J))
            END DO
         END DO
         CALL ELTNOD(N,NT,TP,TETRA,KZ,KZNOD)
         IF (NNEU(2).LT.0) THEN
            CALL VCOPYR(N,CKRWP,CKRW)
            CALL NEUMANN(TIME,NNEU,NNEUC,NNOD,NSTR,CKRW,KZNOD,
     1                   ARENOD,ACONTQ,QINP,Q)
            NNINP=0.0D0
            NNOUTP=0.0D0
            DO K=1,ANQ
               IF (Q(K) .GT. 0.0D0) THEN
                  NNINP=NNINP + Q(K)
               ELSE
                  NNOUTP=NNOUTP + Q(K)
               END IF
            END DO
         END IF

         CALL NODETOTETRA(N,NT,TETRA,SWNEW,SW,volu,volnod,
     1        pnodi,pel)
C-------------------------------------------
C  INITIALIZATION OF TRANSPORT
C ------------------------------------------
         IF (TRANSP) THEN            
C
C     SURFACE PARAMETERS
C           
            IF (SURF) CALL INIT_SR_TRA(NNOD,TRAFLNOD,TRAFLP,CONCNOD)
            DO I=1,NODMAX
               CONCNOD_PERDIR(I)=0.0d0               
               CONCNOD_PERDIR_TIMEP(I)=0.0d0               
               CONCNOD_PERDIR_INT(I)=0.0d0               
            END DO
            SUPDIFFUS=0.0d0
C     
C     SUBSURFACE PARAMETERS
C     concentrations
            CALL CONCINI_NODE(N,CNOLD)
            CALL VCOPYR(N,CNNEW,CNOLD)
            CALL NODETOTETRA(N,NT,TETRA,COLD,CNOLD,volu,volnod,
     1                       pnodi,pel)
c concentration first adsorption equilibrium
            IF (REACFLAG) THEN
               DO I=1,NT
                  COLD(I)= COLD(I)*PEL(I)*SWNEW(I)
     2             /(PEL(I)*SWNEW(I)+KEL(I)*(1-PEL(I))*2)
               END DO
            END IF
c
            CALL VCOPYR(NT,COLDOLD,COLD)
            CALL VCOPYR(NT,CNEW,COLD)
            CALL INIT0R(NT,CDIFF)
C     flux and cflnumb
            CALL INIT0R(NT,FLUX)
            CFLNUMB = CFLINP       
c     saturation at tetra           
            CALL NODETOTETRA(N,NT,TETRA,SWOLD,SW,volu,volnod,
     1                       pnodi,pel)
c     velocities at tetra            
            CALL VCOPYR(NT,UUOLD,UU)
            CALL VCOPYR(NT,VVOLD,VV)
            CALL VCOPYR(NT,WWOLD,WW)            
c     
C CALCULATE THE MASS INITIALLY STORED IN THE SYSTEM 
            CALL MASSCAL(NT,MASS0,COLD,
     1                   SWOLD,PEL,VOLU)
c Calcultation Mass0 with adsorbed part and solute part
            MASS0=0.0d0
            MASSOUTTOTADV=0.0d0
            MASSINTOTADV=0.0d0
CM In and out mass due to dispersion is not computed yet
            MASSOUTTOTCUM=0.0d0
            MASSINTOTCUM=0.0d0
            DO I=1,NT
               MASS0=MASS0 + KEL(I)*CNEW(I)*VOLU(I)*(1-PEL(I))*2
     2                 +CNEW(I)*VOLU(I)*PEL(I)
            END DO
            WRITE(IOUT80,1222) MASS0
            WRITE(IOUT81,1223) MASS0
            WRITE(IOUT81,1543)
            WRITE(IOUT81,1541)(I,I=1,NZONE),(I,I=1,NZONE),(I,I=1,NZONE)
C
C Initialisation des variables ajoutées par Laura (mélange couplage)
C
           CALL INIT0R(NNOD,traflnod)
           CALL INIT0R(NNOD,traflnod_flow)
           CALL INIT0R(NNOD,source_mixing)
           CALL INIT0R(NNOD,CONCNODUPD)
           CALL INIT0R(NMAX,CONCSURFVTK)
           CALL INIT0R(NCELL,MIX_CORRECT)
c
c     atmospheric BC for transport
           CALL ATMONE_TRA(NNOD,CAUSPATM,CAUTIATM,IETOCAU,TIME,DELTAT,
     1          ATMCONC,CONCOLD,CONCTIM,CONCINP)                       
c     
c     initialize dirichlet boundary condition for transport            
           CALL BCONE_TRA('   TRANSPORT DIRICHLET',IIN64,NNOD,NPMAX_TRA,
     1           IPRT1,NDIR_TRA,NDIRC_TRA,NP_TRA,NSTR,HTIDIR_TRA,TIME,
     2           DELTAT,CTIM,CINP,CONTP_TRA,PRESC_TRA,ANP_TRA,
     3           ACONTP_TRA)
            IF (ANP_TRA .GT. 0) THEN
               CALL DIRORNEU_FACE(ANP_TRA,ACONTP_TRA,PRESC_TRA,NFACE,
     1              PLIST,ISIDE,NPFA_TRA,CONTPFA_TRA,PRESCFA_TRA)
               CALL BOUNDIR(NFACE,NPFA_TRA,CONTPFA_TRA,PUNTDIR)
            END IF     

            DO I=1,NNOD
               IFCONC(I) = 0
            END DO
C     
            DO I=1,ANP_TRA
               IFCONC(ACONTP_TRA(I)) = -1
            END DO
C               
         END IF
         
C-------------------------------------------
C END OF INITIALIZATION OF TRANSPORT
C ------------------------------------------
C
C  recharge calculation at initial conditions
C
         CALL RECHARGE(NNOD,NSTR,WNOD,ARENOD,PTIMEP,TIME,DELTAT,RECFLOW,
     1                 RECVOL,RECNOD)
C  
C  detailed output at initial conditions
C
         IF (IPRT .GE. 4) THEN
            WRITE(IOUT15,1510)
            WRITE(IOUT13,1520)
         END IF
         CALL DETOUT(0,IPEAT,IPRT,N,NNOD,NUMVP,NSTR,0,0.0D0,NODVP,
     1        SATSUR,PNEW,INDE0,DEF,SW,CKRW,UNOD,VNOD,WNOD,NT,UU,
     2        VV,WW,X,Y,Z,OVFLNOD,ATMACT,ARENOD,PONDNOD,IFATM,PNODI,
     3        RECNOD,QTRANIE,CNOLD)
CM       IF (NCOUT .EQ. 1) THEN
CM      CALL DETOUT_NETCDF(0,NPRT,TIME,NROW,NCOL,NNOD,NSTR,PNEW,SW,
CM   1           PONDNOD)
CM      CALL DIAGN_OUTPUT(0,NPRT,TIME,N,NNOD,NSTR,BASE,ZRATIO,Z,
CM   1          PNODI,ARENOD,VOLNOD,SW,PONDNOD,PNEW)
CM       END IF         
         IF ((IPRT.GE.2).AND.(VTKF.GT.0)) THEN
CM          CALL VTKRIS3DF(NT,N,100,TETRA,0.0d0,PNEW,SW,UU,VV,
CM   1           WW,X,Y,Z)
            CALL VTKRIS3D(NT,N,100,TETRA,PNEW,SW,
     1                    UU,VV,WW,X,Y,Z,KS,0.0d0,VTKF)
            IF (TRANSP) THEN
               CALL VTKRIS3DFC(NT,N,200,TETRA,0.0d0,CNEW,
     1              SW,X,Y,Z)
               CALL VTKRIS3DFCSURF(NT,N,100,TETRA,0.0d0,
     1              CNNEW,X,Y,Z)
            END IF
         END IF
c     IF (NR.GT.0) THEN
c               WRITE(1001,*) TIME,(PNEW(CONTR(I)),I=1,NR)
c              WRITE(1002,*) TIME,(SW(CONTR(I)),I=1,NR)
c               WRITE(1003,*) TIME,(CKRW(CONTR(I)),I=1,NR)
c        END IF
        
C
C  control if arrays are dimensioned correctly
C
         CALL CHECK(N,NMAX,NT,NTEMAX,NNOD,NODMAX,NTRI,NTRMAX,
     1              NZONE,MAXZON,NPRT,MAXPRT,NSTR,MAXSTR,NR,NRMAX,
     2              N1,NTPMAX,NTERM,MAXTRM,ANP,NPMAX,ANQ,NQMAX,
     3              NDIR(2),NP2MAX,ITUNS,MAXIT,NUMVP,MAXVP,
     4              NUM_QOUT,MAXQOUT,NSF,NSFMAX,NSFNUM,NNSFMX,
     5              NROW,ROWMAX,NCOL,COLMAX,
     6              NCELL,MAXCEL,NUMRES,MAXRES,NLKP,MAXLKP,
     7              NUDN,MAXNUDN,NUDT,MAXNUDT,NUDC,MAXNUDC)
C
         WRITE(IOUT4,1060) IOPT,NLRELX,KSLOPE,TOLUNS,TOLSWI,L2NORM
         WRITE(IOUT5,1220) STORE0
         WRITE(IOUT2,1225) STORE0
         IF (NUDN .GT. 0) WRITE(IOUT51,1230)  
         EVAP_EFF=0.0D0     
         INF_TOT=0.0D0     
         VOLUME_OUT=0.0D0
      END IF
      
      WRITE(IOUT57,1232)
      WRITE(IOUT41,1530) NUM_QOUT+1
      IF (SURF) THEN
         WRITE(IOUT41,1540) QOI_SN(NCELNL),(ID_QOUT(I),I=1,NUM_QOUT)
      END IF
C 
C     WRITE(IOUT44,1570)NROW,NCOL 
C     
C---------------------------START NEW TIME STEP (interior loop on time)
C The loop finishes at label INTERIOR 
      DO WHILE (TIME.LE.TMAX)
C      
      IF(FL3D) THEN
C
C  input (if necessary) and interpolate non-atmospheric, non-seepage
C  face Dirichlet boundary conditions for next time level
C
         IF (DELTAT .GE. 1.0E+15) THEN
            TIME=DELTAT
         ELSE
            CALL BCNXT('   NATM, NSF DIRICHLET',IIN8,IOUT2,IOUT19,NNOD,
     1        NPMAX,IPRT1,NDIR,NDIRC,NP,NSTR,HTIDIR,TIME,
     2        PTIM,PINP,CONTP,PRESC,ANP,ACONTP)
C
C  input (if necessary) and interpolate non-atmospheric, non-seepage
C  face Neumann boundary conditions for next time level
C
c           CALL INIT0I(3,ZERO)
c           CALL INIT0I(3,ZEROC)
            CALL BCNXT('     NATM, NSF NEUMANN',IIN9,IOUT2,IOUT20,NNOD,
     1        NQMAX,IPRT1,NNEU,NNEUC,NQ,NSTR,HTINEU,TIME,
     2        QTIM,QINP,CONTQ,Q,ANQ,ACONTQ)
            CALL NEUMANN(TIME,NNEU,NNEUC,NNOD,NSTR,CKRW,KZNOD,
     1                   ARENOD,ACONTQ,QINP,Q)
C  input (if necessary), interpolate, and check for switching of 
C  atmospheric boundary conditions.
C  Note: SWITCH_OLD is the version of SWITCH for the uncoupled,
C  subsurface flow only version of the model. It should be superceded
C  by the new ADRSTN routine, but for the time being we keep the old
C  version of SWITCH active as well.
C         
            IF (TIME.GT.SFVTIM(2)) THEN
               CALL SFVNXT(IIN7,ITERM,IOUT2,IPRT,N,NSFMAX,NNSFMX,
     1           NSF,NSFNUM,NSFNOD,SFV,SFVNUM,SFVNOD,
     2           SFVTIM,Z)
            ELSE
               CALL SFVREC(NSFMAX,NNSFMX,NSF,NSFNUM,NSFNOD,
     1           SFV,SFVNUM,SFVNOD)
            END IF
C
            CALL ATMNXT(NNOD,HSPATM,HTIATM,IETO,TIME,IFATM,ARENOD,
     1               ATMPOT,ATMACT,ATMTIM,ATMINP,DELTAT,
     2               ANP,ANQ,ACONTP,ACONTQ,NSF,NSFNUM,NSFNOD,SCF)
C
C  compute root water uptake for next time step
C  time-variable maximum root depth
C
            CALL ETRAN(N,NNOD,NSTR,ATMPOT,Z,PNEW,PNODI,VEG_TYPE,QTRANIE)
            QTRAN=0.0d0
            DO I=1,N
               QTRAN=QTRAN+QTRANIE(I)
            END DO
            DO I=1,NR
               RWU(I)=0.0d0
               DO J=1,NSTR+1
                  RWU(I)=RWU(I)+QTRANIE((J-1)*NNOD+CONTR(I))
               END DO
            END DO
            WRITE(666,*)TIME,QTRAN,(RWU(I),I=1,NR)
C
            IF (.NOT. SURF) THEN
               CALL SWITCH_OLD(NNOD,IFATM,ATMACT,ATMPOT,PNEW)
            ELSE
               CALL ADRSTN(NNOD,IFATM,ATMPOT,ATMACT,PNEW)
            END IF
         END IF
C
C     input (if necessary) new nudging observation data and update
C     NUDC and NUDCTR
C
         CALL NUDNXT(NUDN,NUDT,NUDC,NUDCTR,WFLAG,
     1               TIME,NUDTIM,NUDTAU,NUDVAL)
C     
         IF (TRANSP) THEN 
            CALL ATMNXT_TRA(NNOD,CAUSPATM,CAUTIATM,IETOCAU,
     1                      TIME,ATMCONC,CONCTIM,CONCINP,DELTAT)
           CALL BCNXT_TRA('   TRANSPORT DIRICHLET',IIN64,NNOD,NPMAX_TRA,
     1           IPRT1,NDIR_TRA,NDIRC_TRA,NP_TRA,NSTR,HTIDIR_TRA,TIME,
     2           CTIM,CINP,CONTP_TRA,PRESC_TRA,ANP_TRA,
     3           ACONTP_TRA)
            CALL DIRORNEU_FACE(ANP_TRA,ACONTP_TRA,PRESC_TRA,NFACE,
     1           PLIST,ISIDE,NPFA_TRA,CONTPFA_TRA,PRESCFA_TRA)
            CALL BOUNDIR(NFACE,NPFA_TRA,CONTPFA_TRA,PUNTDIR)  
            DO I=1,NNOD
               IFCONC(I) = 0
            END DO
C     
            DO I=1,ANP_TRA
               IFCONC(ACONTP_TRA(I)) = -1
            END DO
C     
         END IF
C     update pressure heads for next time level
C
         CALL VCOPYR(N,POLD,PNEW)
         CALL WEIGHT(N,TETAF,PNEW,PTIMEP,PTNEW)
         CALL VCOPYR(N,PTOLD,PTNEW)
      END IF
  100 CONTINUE
      CALL VCOPYI(ANP,NODDIR,ACONTP)
      IF (FL3D) THEN
         WRITE(IOUT2,1010) NSTEP,DELTAT,TIME
         WRITE(ITERM,1010) NSTEP,DELTAT,TIME
         WRITE(IOUT4,1065) NSTEP,DELTAT,TIME
      ELSE IF ( .NOT. FL3D) THEN
         WRITE(ITERM,1010) NSTEP,DELTATS,TIME
      END IF
C
C  call surface routing module coupled or not
C       
      NSURF=0
      IF (SURF .AND. FL3D .AND. PONDING) THEN
         AK_MAX_SAV = AK_MAX
         CALL VCOPYR(MAXCEL,Q_IN_KK_SN_SAV,Q_IN_KK_SN)
         CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_1_SAV,Q_OUT_KK_SN_1)
         CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_2_SAV,Q_OUT_KK_SN_2)
         CALL VCOPYR(NUMRES,H_POOL_KK_VEC_SAV,H_POOL_KK_VEC)
         CALL VCOPYR(MAXCEL,VOLUME_KK_SN_SAV,VOLUME_KK_SN) 
C
C SURFACE TRANSPORT      
         IF (TRANSP) THEN
            AK_MAX_TRA_SAV = AK_MAX_TRA
            CALL VCOPYR(MAXCEL,QMASS_IN_KK_SN_SAV,QMASS_IN_KK_SN)
            CALL VCOPYR(MAXCEL,QMASS_OUT_KK_SN_1_SAV,QMASS_OUT_KK_SN_1)
            CALL VCOPYR(MAXCEL,QMASS_OUT_KK_SN_2_SAV,QMASS_OUT_KK_SN_2)
            CALL VCOPYR(MAXCEL,MASS_KK_SN_SAV,MASS_KK_SN)      
         END IF 
C
         CALL TIM(PCVEC1,1)
C SURFACE ROUTING FOR BOTH SURFACE WATER AND MASS
C
         CALL INIT0R(NNOD,TRAFLNOD)
         DO I=1,NNOD
            TRAFLNOD(I)=(TRAFLNOD_FLOW(I)+SOURCE_mixing(i))/deltat
         END DO
C     
         CALL SURF_FLOWTRA(NCELL,NNOD,NROW,NCOL,NTRI,DOSTEP,NUMRES,
     1        NCELL_COARSE,NCELNL,CELL,
     2        INDCEL,INDCELWL,TIPO_R,RESERVR,DEM_MAP,
     3        LAKES_MAP,CELLCOL,CELLROW,
     4        CELTYPE,CELLS_R,TP2D,TRIANG,
     5        TIME,DELTAT,ARENOD,OVFLNOD,PONDNOD,
     6        OVFLCEL,PONDCEL,TRAFLNOD,TRAFLCEL,CONCNOD,
     7        CONCCEL,N_HA,NSURF,NSURFT,NSURFT_TB,TRANSP,
     8        SOURCE_MIXING,MIX_CORRECT,SURFACE_MIX_SN)
C
C     
         CALL VCOPYR(NT,CNEW_SAV,CNEW)
         CALL VCOPYR(N,SW_SAV,SW)
         IF (PONDING.and.TRANSP) THEN
            CALL MIXING_CORRECTION(NROW,NCOL,indcel,indcelwl,
     1           N,NT,NNOD,TETRA,nstr,NCELL,CELL,MIX_CORRECT,
     2           CNEW,CNNEW,PEL,PNODI,
     3           VOLU,VOLNOD,SWNEW,SW,TP,COLD,CNOLD)
            CALL VCOPYR(NT,COLD,CNEW)
         END IF
         
         CALL TIM(PCVEC1,2)
         CPUSURF=CPUSURF+PCVEC1
         CPUSURF_T=CPUSURF_T+PCVEC1
c     volume_out=volume_out+(Q_outlet_kk+Q_outlet_kkp1)/2*deltat
ccc   write(81,*) time,volume_out
C     
C     update surface BCs for FLOW3D to reflect ponding situation as
C  calculated in SURF_ROUTE, and update PONDING flag
C
         CALL PONDUPD(NNOD,PONDING,IFATM,PONDH_MIN,DELTAT,PONDNOD,
     1        ARENOD,ATMPOT,ATMACT,PNEW)
C
C update of concentration for rainfall
C
         IF (TRANSP) THEN
              CALL CONCUPD(NNOD,IFATM,IFATMP,DELTAT,PONDNOD,
     1        ARENOD,ATMPOT,ATMACT,CONCNOD,ATMCONC,
     2        CONCNODUPD)
         END IF
C     
      ELSE IF (.NOT. FL3D) THEN   
         CALL ROUTE(NROW,NCELL,TIPO_R,RESERVR,
     1        CELLS_R,N_HA,TIME,DELTATS)
         CALL ALTEZZE(NCELNL,CELTYPE,DELTA_X,
     1        DELTA_Y,Q_IN_KK_SN,Q_IN_KKP1_SN,
     2        Q_OUT_KK_SN_1,Q_OUT_KK_SN_2,
     3        Q_OUT_KKP1_SN_1,Q_OUT_KKP1_SN_2,
     4        VOLUME_KK_SN,VOLUME_KKP1_SN,SURFACE_WATER_SN,
     5        H_WATER_KKP1_SN,DELTATS,
     6        H_POOL_KKP1_VEC,RESERVR,ELEVATION)
C ----------------------------------------------
C                 SURFACE TRANSPORT 
C ----------------------------------------------
         IF (TRANSP) THEN 
C     MASS ROUTING AT THE SURFACE
            CALL ROUTE_TRA(NROW,NCELL,TIPO_R,RESERVR,
     1           CELLS_R,N_HA,TIME,DELTATS)
            
C     MASS BALANCE CALCULATION
            CALL MASBIL(NCELL,TIPO_R,DELTA_X,DELTA_Y,
     1           QMASS_IN_KK_SN,QMASS_IN_KKP1_SN,
     2           QMASS_OUT_KK_SN_1,QMASS_OUT_KK_SN_2,
     3           QMASS_OUT_KKP1_SN_1,QMASS_OUT_KKP1_SN_2,
     4           MASS_KK_SN,MASS_KKP1_SN,SURFACE_CONC_SN,
     5           CONC_KKP1_SN,VOLUME_KKP1_SN,DELTATS,MIX_CORRECT,
     6           SURFACE_MIX_SN)     
         END IF
C END OF SURFACE TANSPORT
c     
      END IF
C
      IF (SURF) THEN
         VOLUME_SUP= 0.0d0
         DO I=1,NCELL
            VOLUME_SUP=VOLUME_SUP+PONDCEL(I)*DELTA_X*DELTA_Y
         END DO
      END IF

C  call subsurface flow solver
      IF (FL3D) THEN
         CALL VCOPYR(NNOD,POLD,PNEW)
         CALL TIM(PCVEC1,1)
         CALL FLOW3D(IOPT,IPRT1,NLRELX,L2NORM,N,NT,NTRI,NTERM,NSTR,
     1        NDIR(2),NDIRC(2),ANP,ANQ,HSPATM,HTIATM,HTIDIR,
     2        HTINEU,ITER,ITUNS,NITER,NITERT,ITLIN,ITRTOT,
     3        ISOLV,ITMXCG,IPEAT,IVGHU,KSLOPE,ISFONE,KSFCV,
     4        KSFCVT,IBOT,MINBOT,NNOD,NUMDIR,NSF,KLSFAI,
     5        NLKP,NUDN,NUDC,NUDCTR,NUDFLAG,WFLAG,
     6        TP,TETRA,TETJA,IA,JA,TOPOL,ACONTP,ACONTQ,IVOL,
     7        INSYM,IFATM,IFATMP,NODDIR,NSFNUM,NSFNOD,
     8        SFEX,SFEXIT,SFEXP,SFFLAG,NUDTET,DUPUIT,
     9        SURF,PONDING,
     A        DTGMIN,LSFAIL,ERRGMX,NORMCV,ITAGEN,SFCHEK,
     B        KSFZER,NOBACK,BCKSTP,CPUVEC,CPUNL,TIME,TIMEP,
     C        DELTAT,DTMIN,DTREDM,DTREDS,TETAF,OMEGA,OMEGAP,
     D        RMAX,RMIN,TOLUNS,TOLSWI,TOLCG,ERNLMX,
     E        PONDH_MIN,NUDG,
     F        PONDNOD,OVFLNOD,PERMX,PERMY,PERMZ,INDE0,
     G        SNODI,PNODI,PEL,PORE,INDE,PRESC,PTIM,PINP,DEF,
     H        Q,QTIM,QINP,ATMTIM,ATMINP,AI,BI,CI,DI,LMASS,
     I        ARENOD,VOLNOD,VOLU,VOLUR,
     J        ETAI,ETAE,ET1,ET2,ET1E,DETAI,DETAIE,SW,SWP,SWE,
     K        CKRW,CKRWE,DCKRW,DCKRWE,SENODI,SEELT,
     L        POLD,PDIFF,PTOLD,PNEW,PTNEW,PTIMEP,
     M        TNOTI,XT5,LHSP,COEF1,COEF2,COEF3,COEF4,SCR,
     N        RNSYM,ATMACT,ATMPOT,ATMOLD,LHSATM,LHSSF,
     O        QPNEW,QPOLD,SFQ,SFQP,BKCFLOW_NODE,X,Y,Z,
     P        NUDTIM,NUDRXY,NUDRZ,NUDX,NUDY,NUDZ,NUDCUM,
     Q        NUDEPS,NUDDIF,NUDSMC,NUDNOD,NUDVAL,NUDTAU,
     R        PUNTDIRSF_NODE,QTRANIE)
         CALL TIM(PCVEC1,2)
         CPUSUB=CPUSUB+PCVEC1
         CPUSUB_T=CPUSUB_T+PCVEC1
C ----------------------------------------------------------------------
C  do back-stepping if necessary
C ----------------------------------------------------------------------
         IF (BCKSTP) THEN
            IF (SURF) THEN
               AK_MAX_P=AK_MAX_SAV
               CALL VCOPYR(MAXCEL,Q_IN_KK_SN_P,Q_IN_KK_SN_SAV)
               CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_1_P,Q_OUT_KK_SN_1_SAV)
               CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_2_P,Q_OUT_KK_SN_2_SAV)
               CALL INIT0R(MAXCEL,Q_IN_KKP1_SN)
               CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_1)
               CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_2)           
               CALL VCOPYR(MAXCEL,VOLUME_KK_SN_P,VOLUME_KK_SN_SAV)
               CALL INIT0R(MAXCEL,VOLUME_KKP1_SN)
               CALL VCOPYR(NUMRES,H_POOL_KK_VEC_P,H_POOL_KK_VEC_SAV)
            END IF
            IF (TRANSP) THEN 
               CALL VCOPYR(MAXCEL,QMASS_IN_KK_SN_P,QMASS_IN_KK_SN_SAV)
               CALL VCOPYR(MAXCEL,QMASS_OUT_KK_SN_1_P,
     1              QMASS_OUT_KK_SN_1_SAV)
               CALL VCOPYR(MAXCEL,QMASS_OUT_KK_SN_2_P,
     1              QMASS_OUT_KK_SN_2_SAV)
               CALL INIT0R(MAXCEL,QMASS_IN_KKP1_SN)
               CALL INIT0R(MAXCEL,QMASS_OUT_KKP1_SN_1)
               CALL INIT0R(MAXCEL,QMASS_OUT_KKP1_SN_2)           
               CALL VCOPYR(MAXCEL,MASS_KK_SN_P,MASS_KK_SN_SAV)
               CALL INIT0R(MAXCEL,MASS_KKP1_SN)
               CALL VCOPYR(NT,CNEW,CNEW_SAV)
               CALL VCOPYR(N,SW,SW_SAV)
            END IF  
           CALL BKSTEP(IPRT1,N,NSTR,NDIR,NDIRC,NP,NNEU,NNEUC,NQ,HSPATM,
     1           HTIATM,IETO,HTIDIR,HTINEU,ITER,NITERT,KBACKT,KBACK,
     2           NNOD,NSF,NUMRES,NCELNL,
     3           CONTP,CONTQ,IFATM,IFATMP,NSFNUM,SFEX,SFEXIT,
     4           SFEXP,SURF,PONDING,PONDP,DTGMIN,
     5           TIME,TIMEP,DELTAT,DTMIN,DTREDM,DTREDS,TETAF,
     6           OVFLNOD,OVFLP,PRESC,PTIM,PINP,Q,QTIM,QINP,
     7           ATMTIM,ATMINP,ARENOD,
     8           POLD,PTOLD,PNEW,PTNEW,PTIMEP,
     9           ATMACT,ATMPOT,ATMOLD,
     A           QPNEW,QPOLD,SFQ,SFQP,ANP,ANQ,ACONTP,ACONTQ,
     B           NSFNOD,SFV,SFVNUM,SFVNOD,SFVTIM,TRANSP)
            CALL NEUMANN(TIME,NNEU,NNEUC,NNOD,NSTR,CKRWP,KZNOD,
     1                   ARENOD,ACONTQ,QINP,Q)
            GO TO 100
         END IF
         IF (SURF) THEN
            NSURFT_T=NSURFT_T + NSURF
            DO I=1,NNOD
               IF(PNEW(I) .LE. 0.0d0) PONDNOD(I) = 0.0d0
            END DO
            VOLSUPNEW=0.0d0
            DO I=1,NNOD
               VOLSUPNEW=VOLSUPNEW+PONDNOD(I)*ARENOD(I)
            END DO
         END IF
      END IF
C

C VELOCITY RECONSTRUCTION
      
         CALL UPD_FACETYPE(N,NNOD,NFACE,PLIST,ISIDE,Faces_FaceType,
     1        NODE_FACECOUNT,NODE_FACEIDS,
     2        ATMACT,ARENOD,AREFACE,Faces_NeumannFLux,
     3        PUNTDIRFLOW_NODE,IFSF,PUNTDIRSF_NODE,IFATM,
     5        PUNTNEUFLOW_NODE,CONNEU_NODE)
C              
C        

      IF(VELREC.EQ.1)THEN
c        write(*,*) UU(1),VV(1),WW(1)
c        write(*,*) Node_ElementCount(1),
c    3    Node_ElementIDs(1,1),ISIDE(1,1),PLIST(1,1),FaceCentroid(1,1),
c    4    AREFACE(1),Faces_NeumannFlux(1),Faces_FaceFlux(1,1),
c    5    Node_FaceCount(1),Node_FaceIDs(1,1),Faces_FaceType(1),
c    6    NFACE,Faces_ReferenceVector(1,1),SIDE_CNC(1,1),
c    7    XC(1),YC(1),ZC(1),time
         CALL VELOCITY_MAIN(N,X,Y,Z,PNEW,PTIMEP,SW,SWP,
     1        DELTAT,NT,TETRA,VOLU,IVOL,ET1E,PEL,
     2        UU,VV,WW,Node_ElementCount,
     3        Node_ElementIDs,ISIDE,PLIST,FaceCentroid,
     4        AREFACE,Faces_NeumannFlux,Faces_FaceFlux,
     5        Node_FaceCount,Node_FaceIDs,Faces_FaceType,
     6        NFACE,Faces_ReferenceVector,SIDE_CNC,
     7        XC,YC,ZC,time)
c        write(*,*) UU(1),VV(1),WW(1)
c        write(*,*) Node_ElementCount(1),
c    3    Node_ElementIDs(1,1),ISIDE(1,1),PLIST(1,1),FaceCentroid(1,1),
c    4    AREFACE(1),Faces_NeumannFlux(1),Faces_FaceFlux(1,1),
c    5    Node_FaceCount(1),Node_FaceIDs(1,1),Faces_FaceType(1),
c    6    NFACE,Faces_ReferenceVector(1,1),SIDE_CNC(1,1),
c    7    XC(1),YC(1),ZC(1),time
      END IF
C
C-----------VELOCITY CALCULATION AND HYDROGRAPH AND MASS BALANCE OUTPUT
C
C
C  calculate soil moisture characteristics needed for storage and
C  velocity calculations and for output 
C
      IF (IPEAT .EQ. 0) THEN
         CALL CHVELO(NLKP,N,NT,NTRI,IVGHU,TETRA,TP,
     1               PNEW,SNODI,PNODI,SW,CKRW,SWE,CKRWE)
         CALL STORCAL(N,STORE1,STORE2,DSTORE,SW,PNODI,VOLNOD)
      ELSE
         CALL CHVELOP(NLKP,NNOD,N,NT,NTRI,IVGHU,TETRA,TP,
     1                PNEW,Z,PORE,SW,CKRW,SENODI,CKRWE,SEELT)
         CALL STORCAL(N,STORE1,STORE2,DSTORE,SW,PORE,VOLNOD)
      END IF
      IF (FL3D)THEN
         IF (IPRT .GE. 2) THEN
            CALL NODELT(NT,TETRA,CKRW,CKRWE)
            CALL VEL3D(NT,TETRA,NTRI,PNEW,PERMX,PERMY,PERMZ,
     1           UU,VV,WW,BI,CI,DI,CKRWE,VOLUR,IVOL)
            CALL VNOD3D(N,NT,TP,TETRA,UU,VV,WW,UNOD,VNOD,WNOD)
         END IF 
C
C  output mass balance errors
C
         WRITE(IOUT2,1130) ADIN, ADOUT, NDIN, NDOUT,
     1        ANIN, ANOUT, NNIN, NNOUT,
     2        ADINP,ADOUTP,NDINP,NDOUTP,
     3        ANINP,ANOUTP,NNINP,NNOUTP,
     4        VADIN,VADOUT,VNDIN,VNDOUT,
     5        VANIN,VANOUT,VNNIN,VNNOUT
         IF (NSF .GT. 0) WRITE(IOUT2,1135) SFFLW,SFFLWP,VSFFLW
         WRITE(IOUT2,1140) VIN,VOUT,DSTORE,ERRAS,ERREL
         VSFTOT=VSFTOT + VSFFLW
         VNDTOT=VNDTOT + (VNDIN + VNDOUT)
         VNNTOT=VNNTOT + (VNNIN + VNNOUT)
         VNUDTOT=VNUDTOT + (VNUDIN + VNUDOUT)
         VTOT=VTOT + (VIN + VOUT)
         CVIN=CVIN + VIN
         CVOUT=CVOUT + VOUT
         CDSTOR=CDSTOR + DSTORE
         CERRAS=CERRAS + ERRAS
         CAERAS=CAERAS + DABS(ERRAS)
         WRITE(IOUT5,1240) NSTEP,DELTAT,TIME,ITER,
     1        DFLOAT(NITERT)/DFLOAT(ITER),
     2        STORE1,STORE2,DSTORE,CDSTOR,
     3        VIN,CVIN,VOUT,CVOUT,VIN+VOUT,VTOT,
     4        ERRAS,ERREL,CERRAS,CAERAS
         ADINP =ADIN
         ADOUTP=ADOUT
         NDINP =NDIN
         NDOUTP=NDOUT
         ANINP =ANIN
         ANOUTP=ANOUT
         NNINP =NNIN
         NNOUTP=NNOUT
         SFFLWP=SFFLW
         NUDINP=NUDIN
         NUDOUTP=NUDOUT
         IF (NSF.GT.0) THEN
            DO I=1,NSF
               SFQ_PARTIAL(I)=0.0d0
               DO J=1,NSFNUM(I)
                  SFQ_PARTIAL(I)=SFQ_PARTIAL(I)+SFQ(I,J)
               END DO
            END DO
c            WRITE(111,666) TIME,(SFQ_PARTIAL(I),I=1,NSF)
         END IF
         DO I=1,NSF
            DO J=1,NSFNUM(I)
               SFEXP(I,J)=SFEX(I,J)
               SFEXIT(I,J)=SFEX(I,J)
               SFQP(I,J)=SFQ(I,J)
            END DO
         END DO
      END IF

C--------------------------------------------------     
C               SUBSURFACE TRANSPORT 
C--------------------------------------------------
      IF (TRANSP .AND. FL3D) THEN
C
C     check whether we have diffusion or not
C     
         CALL CONTROLDIFF(NSTR,NZONE,ALFAL,ALFAT,DIFFUS,CDIFFUS)
C     
C     SET BC FOR TRANSPORT
C     
C ADVECTION ONLY - COUPLING BASED ON DIRICHLET
         CALL CONC_BC(NNOD,IFATM,IFATMP,DELTAT,PONDNOD,
     1   PONDNODP,ARENOD,ATMPOT,ATMACT,CONCNOD,CONCNOD_OLD,ATMCONC,
     2   CNNEW,ANP_TRA,ACONTP_TRA,PRESC_TRA,CONCNODBC,CONCNODUPD)
C
CM       IF (CDIFFUS .NE. 0) THEN 
C ADVECTION + DIFFUSION / COUPLING BASED ON CAUCHY - SYLVAIN POST-DOC
CM         CALL CAUCHYUPD(NNOD,NQ_TRA,CONTQ_TRA,QCAUCHY,IFATM,IFCONC,
CM   1     ATMPOT,ATMACT,OVFLNOD,PNEW,CONCNOD,DELTAT,ATMCONC,
CM   2     ARENOD,CNNEW,IFATMP,PONDNOD,CONCNODUPD)
CM         write(111,*)(ARENOD(I),I=1,NNOD)
CM         write(222,*)(IFATM(I),I=1,NNOD)
CM         write(333,*)(QCAUCHY(I),I=1,NNOD)
CM       END IF
C     
C Calculation of SWNEW for each element from the results of FLOW3D SW
C         
         CALL NODETOTETRA(N,NT,TETRA,SWNEW,SW,volu,volnod,
     1                    pnodi,pel)
C     
         CALL CFLPECNUMBER(NT,DELTATADV,VOLUR,SURFAR,UUOLD,VVOLD,WWOLD,
     1        SUPDIFFUS,CFLINP,CFLNUMB,PECNUMB,CFLTETRA,NADV,NADVFL,
     2        DELTAT)
         WRITE(IOUT2,*)'SUPDIFFUS = ',SUPDIFFUS
         WRITE(IOUT2,*)'CFLNUMB = ',CFLNUMB
         WRITE(IOUT2,*)'PECNUMB = ',PECNUMB
         TIMEAUX = TIME - DELTAT
         if (cflnumb.eq.0.d0) go to 350

C--------------------------------------------------------------
C                         ADVECTION                                          
C--------------------------------------------------------------
C           
        
         CALL VCOPYR(NT,SWPINT_EL,SWOLD)
C         
CC CARLOTTA ACCORDING TO THE FLOW CONDITION, I NEED TO SET THE KIND OF
CC       FACE I HAVE (IBNDRY). PER ORA CONSIDERO SF AND ATMOSFERIC
C
         CALL FACETYPE(N,NFACE,PLIST,ISIDE,IFSF,IFSFP,IFATM,IFATMP,
     1        IBNDRYP,IBNDRY,NODE_FACECOUNT,NODE_FACEIDS,AREFACE)
CC
CC INITIALIZE VARIABLES FOR MASSBALANCE CALCULATION 
CC
         MASSINADV=0.0d0
         MASSOUTADV=0.0d0
         MASSINTOT=0.0d0
         MASSOUTTOT=0.0d0
C         
         MASSSURF=0.0D0
         DO I=1,NNOD
            MASSSURF=MASSSURF+(CONCNODUPD(i)*
     1            (arenod(i)*pondnod(i)+(atmpot(i)-atmact(i))*deltat))
         END DO
c
C
         DO TADV=1,NADV
C            
C     depending on flag_temp, interpolation is made at the beginning or in the middle of
C     the advective sub ime step..... 
            IF ((FLAG_TEMP .EQ. 0) .OR. (FLAG_TEMP .EQ. 10)) THEN
               TIMWGT = 0.
            ELSE
               TIMWGT = 1.
            END IF
C     
C     Interpolation of velocity and saturation for the advective sub time step
C     
            CALL TIME_INTERPOLATION(N,NT,NFACE,NNOD,TADV,NADV,
     1           VELREC,SWINT_EL,SWOLD,SWNEW,
     2           UUINT,VVINT,WWINT,UUOLD,VVOLD,WWOLD,UU,VV,WW,
     3           ATMINT,ATMOLD,ATMACT,FACEFLUX,FACES_FACEFLUXOLD,
     4           FACES_FACEFLUX,BKCFLOWINT_NODE,BKCFLOWP_NODE,
     5           BKCFLOW_NODE)
C            
C    BASED ON THE VELOCITIES (UUINT,VVINT,WWINT) or FLUX ON EACH FACE 
C    (FACEFLUX)--VELREC ONLY CASE--CALCULATE THE NORMAL FLUX ON 
C    EACH FACE  
C
            CALL FLUXTRA(N,NFACE,NT,VELREC,
     1           PLIST,ISIDE,SIDE_CNC,
     2           Node_FaceCount,Node_FaceIDS,
     3           X,Y,Z,XC,YC,ZC,Faces_ReferenceVector,FaceCentroid,
     4           AREFACE,UUINT,VVINT,WWINT,FACEFLUX,
     5           Faces_FaceType,Faces_NeumannFlux,
     6           BKCFLOWINT_NODE,ARENOD,VN)
C
           DO I=1,NNOD
             CONCNOD_PERDIR(I)=CONCNODBC(I)
           END DO 
C
CC  CHE FAI FINIRE QUI
c  CARLOTTA UPDATE NODAL CONCENTRATION FOR ATMOSPHERIC DIRICHELT
C  CONDITION. IT IS UPDATED AT THE CURRENT ADVECTIVE TIME

C laura
c         DO I=1,NNOD
c            CONCNOD_PERDIR_INT(I)=CONCNOD_PERDIR_TIMEP(I)+
c     1      (((FLOAT(TADV)/FLOAT(NADV)))
c     2      *(CONCNOD_PERDIR(I)-CONCNOD_PERDIR_TIMEP(I)))
c         END DO
         
          DO I=1,NNOD
            CONCNOD_PERDIR_INT(I)=CONCNOD_PERDIR(I)
          END DO


c  CARLOTTA SUBROUTINE FOR THE IMPOSITION OF DIRICHLET BOUNDARY
c  CONDITIONS. FOR THE MOMENT I AM JUST LOOKING AT ATMOSPHERIC BOUNDARY

           CALL DIRFACE_UPDATE(NNOD,NFACE,ISIDE,ATMINT,
     1        ARENOD,CONCNOD_PERDIR_INT,NPFA_TRA,CONTPFA_TRA,
     2        PRESCFA_TRA,PUNTDIR)
C
C CARLOTTA UPDATE BOUNDARY FLUX. SINCE THERE IS NO VELOCITY
C          RECONSTRUCTION, I MAKE SURE THAT THE FLUX RELATIVE TO
C          THE FACE BOUNDARY EQUALS THE AVERAGE OF THE THREE NODAL 
C          FLUX (ex. at least for surface ATMACT)  
C
           CALL VOLFIN(N,NTRI,NT,NFACE,NPFA_TRA,N1,NSTR,NZONE,
     1          FLAG_TEMP,FLAG_INTERP,FLAG_LIMITER,PLIST,PUNTDIR,
     2          SIDE_CNC,ISIDE,NEIGH,NEIGBC,Node_ElementIDs,
     3          Node_ElementCount,TETRA,
     3          IP4,X,Y,Z,XC,YC,ZC,
     4          DELTATADV,VOLUR,VOLU,VN,SWPINT_EL,SWINT_EL,
     5          POROS,COLD,CNEW,
     6          CDIFF,PRESCFA_TRA,FLUX,precl,precr,timeaux)
C
C
           CALL MASS_ADV(NFACE,DELTATADV,PLIST,VN,PRECL,PRECR,
     1            MASSOUTADV,MASSINADV,MASSOUTTOTADV, 
     2            MASSINTOTADV)
C
            CALL VCOPYR(NT,COLD,CNEW)
            CALL VCOPYR(NT,SWPINT_EL,SWINT_EL)
            CALL VCOPYR(NNOD,CONCNOD_PERDIR_TIMEP,CONCNOD_PERDIR)
c           CALL ELTNOD_TRA(N,NT,TETRA,COLD,CNOLD,PEL,PNODI,VOLU,
c     1          VOLNOD,SWNEW,SW)
C     Update of timeaux/ timwgt updated only at the begining of the new time step
C     
            TIMEAUX = TIMEAUX + DELTATADV
C     
         END DO
C         
c CARLOTTA Calculate total mass of solute stored (MASST) after advective step 
            CALL MASSCAL(NT,MASST,CNEW,SWNEW,PEL,VOLU)
            CALL MASSMB(NT,DSMASS,CNEW,COLDOLD,SWNEW,
     1            SWOLD,PEL,VOLU)
c
C--------------------------------------------------------------
C                       END OF  ADVECTION                                     
C--------------------------------------------------------------
 350     CONTINUE
C         
         IF (CDIFFUS .NE. 0) THEN 
C     
            IF (CFLNUMB .NE. 0.D0) THEN 
               CALL ELTNOD_TRA(N,NT,TETRA,COLD,CNOLD,PEL,PNODI,VOLU,
     1                         VOLNOD,SWNEW,SW,NNOD)
            END IF
c            
            CALL INIT0R(NODMAX,QCAUCHY)
            CALL UT3D(iout1,NNOD,NSTEP,N,NTRI3,NT,NTERM,ANP_TRA,NQ_TRA,
     1           NIAUX,NRAUX,TETRA,IP4,IAUX,IPARM,AUX,
     2           TETJA,JA,ACONTP_TRA,CONTQ_TRA,TOPOL,
     3           DIFFUS,GAMMAS,RMAX,DELTAT,TETAC,TOLCG,
     4           X,Y,Z,UU,VV,WW,SWNEW,VOLU,VOLUR,VOLNOD,CNNEW,
     5           CNOLD,TIME,XT5C,TNOTIC,PRESC_TRA,QCAUCHY,POROS,
     6           KD,ALFAL,ALFAT,LMASSC,LHSC,COEF1C,
     7           COEF2C,CPUVEC,IPRT1,SUPDIFFUS,ATMACT)
C
C     Calculate the concentration in tetrahedra from concentration at nodess
C    
            CALL NODETOTETRA(N,NT,TETRA,CNEW,CNNEW,volu,volnod,
     1                       pnodi,pel)
c
         END IF

c           
C----------------------------------------------------------------------
C                 SUBSURFACE REACTION ...... ON ELEMENTS !!!
C----------------------------------------------------------------------
C LAURA LINEAR ADSORPTION AND FIRST ORDER DECAY
cC        
      IF (REACFLAG) THEN
        DO I=1,NT
           CNEW(I)=(KEL(I)*COLDOLD(I)*(1-PEL(I))*GAMMAS        
     1             +CNEW(I)*PEL(I)*SWNEW(I))
     2             /(PEL(I)*SWNEW(I)+KEL(I)*(1-PEL(I))*GAMMAS)
           CNEW(I)=CNEW(I)*exp(LEL(I)*DELTAT*(-1))
        END DO
c
        CALL VCOPYR(NT,COLD,CNEW)
        CALL VCOPYR(NT,COLDOLD,CNEW)
      END IF
C----------------END OF SUBSURFACE REACTION
C
         DO I=1,NT
            CDIFF(I) = 0.5D0*(CNEW(I) - COLDOLD(I))
         END DO
C
C Compute here the mass under each zone because the subroutine is based 
C on elements 
        CALL MASS_ZONE(NT,TETRA,NZONE,MASSZONE,CNEW,VOLU,PEL,SWNEW,
     1     VOLZONE,MASSSOLZONE,KEL)
             WRITE(IOUT81,1542) TIME,(VOLZONE(I),I=1,NZONE),
     1             (MASSZONE(I),I=1,NZONE),(MASSSOLZONE(I),I=1,NZONE)

C CARLOTTA COPY OLD FACE_FACEFLUX
C
         DO I=1,NFACE
            Faces_FaceFluxOld(1,I)=Faces_FaceFlux(1,I) 
            Faces_FaceFluxOld(2,I)=Faces_FaceFlux(2,I) 
            Faces_FaceFluxOld(3,I)=Faces_FaceFlux(3,I) 
         END DO 
c         
        CALL ELTNOD_TRA(N,NT,TETRA,CNEW,CNNEW,PEL,PNODI,VOLU,
     1           VOLNOD,SWNEW,SW,NNOD)
        CALL VCOPYR(N,CNOLD,CNNEW)
c
C   output : mass through the seepage face
c     
      SFMASS=0.0D0
      DO I=1,NSF
         DO J=1,NSFNUM(I)
            IF (SFEX(I,J) .EQ. 1) THEN
               K=NSFNOD(I,J)
               SFMASS=SFMASS + SFQ(I,J)*CNNEW(K)*DELTAT
               IF ((SFMASS.LT.0.0d0).AND.(SFMASS.GT.-1e-30)) THEN
                  SFMASS=0.0d0
               END IF
            END IF
         END DO
      END DO
C
C TRANSPORT SOURCE TERM DEFINITION  
C
         CALL SOURCE_TRASUR(NNOD,IFATM,IFATMP,IFCONC,DELTAT,ATMPOT,
     1        ATMACT,PONDNOD,CONCNOD,CNNEW,ATMCONC,TRAFLNOD_FLOW,
     2        ARENOD,CONCNODUPD,PONDNODP)
c    
         DO I=1,NNOD
             TRAFLNOD_FLOW(I)=TRAFLNOD_FLOW(I)*deltat
         END DO
c
C IF PONDED, MIXING OF SURFACE WATER WITH FIRST LAYER WATER
C  
        CALL INIT0R(NNOD,source_mixing)
        IF (MIXPART.NE.0.0d0) THEN
           IF (PONDING) THEN
             CALL CONCUPD_MIXING(NNOD,N,NT,TETRA,PONDNOD,VOLNOD,PNODI,
     1       SW,ARENOD,CONCNODUPD,CNNEW,CNEW,MIXPART,TP,PEL,VOLU,SWNEW,
     2       DELTAT,SOURCE_MIXING,ATMACT,TRAFLNOD_FLOW,OVFLNOD,PONDNODP)
C
             CALL nodetotetra_tra(N,NT,TETRA,nstr,CNEW,CNNEW,PEL,
     1            PNODI,VOLU,VOLNOD,SWNEW,SW,TP,NNOD,COLD,CNOLD)
C
           END IF 
        END IF  
c     
        CALL VCOPYR(NT,COLD,CNEW)
        CALL VCOPYR(NT,CNEW_VTK,CNEW)
        CALL VCOPYR(N,CNOLD,CNNEW)
c       CALL VCOPYR(NT,COLDOLD,CNEW)
        CALL VCOPYR(NT,SWOLD,SWNEW)
        CALL VCOPYR(NT,UUOLD,UU)
        CALL VCOPYR(NT,VVOLD,VV)
        CALL VCOPYR(NT,WWOLD,WW)   
C   
        CALL VCOPYR(NNOD,CONCNOD_OLD,CONCNOD)
        CALL VCOPYR(NNOD,CONCNOD,CONCNODUPD)

C CREATION OF CONCSURFVTK FROM CONCNODUPD 
            DO I=1,NNOD
               CONCSURFVTK(I)=CONCNODUPD(I)
            END DO
            DO I=NNOD+1,N
               CONCSURFVTK(I)=0.0d0
            END DO
C
      END IF
C
C
         CALL TRANSPMB(N,NNOD,ANP,ANQ,ACONTP,ACONTQ,QPNEW,Q,ATMACT,
     1                 IFATM,CNNEW,DELTAT,SFMASS,MASSINTOT,MASSOUTTOT,
     2                 MASSINTOTCUM,MASSOUTTOTCUM)
         WRITE(IOUT80,1242) NSTEP,DELTAT,TIME,NADV,
     1        MASST,DSMASS,-1.0d0*MASSINADV,MASSINTOT-DABS(MASSINADV),
     2        -1.0d0*MASSOUTADV,MASSOUTTOT+MASSOUTADV,
     3        MASSINTOTCUM,MASSOUTTOTCUM,MERRAS,CMERRAS,MASSSURF
C
c
C----------------------------------------------------------------------
C                           END OF SUBSURFACE TRANSPORT
C----------------------------------------------------------------------
C FOR TRANSPORT
      CALL VCOPYI(NNOD,IFATMP,IFATM) 

C
C  OUTPUT OF THE DISCHARGE AT THE OUTLET OF THE BASIN AND IN THE OTHER
C  CELLS SUCH AS SELECTED IN PARM INPUT FILE 
C        
      IF (SURF) THEN        
         CALL DETOUTQ(NCELNL,NUM_QOUT,TIME,ID_QOUT,QOI_SN,
     1        Q_IN_KKP1_SN,Q_OUT_KKP1_SN_1,Q_OUT_KKP1_SN_2,TRANSP,
     1        QMASS_IN_KKP1_SN,QMASS_OUT_KKP1_SN_1,QMASS_OUT_KKP1_SN_2)
CM    IF (NCOUT .EQ. 1) THEN
CM       CALL DETOUTQ_NETCDF(NCELNL,NUM_QOUT,NSTEP,TIME,ID_QOUT,QOI_SN,
CM   1        Q_IN_KKP1_SN,Q_OUT_KKP1_SN_1)
CM    END IF
      END IF 
C
C output of the water table depth at the NODVP nodes
C
      IF (NUMVP .GT. 0) CALL WTDEPTH(NUMVP,NODVP,NSTR,NNOD,TIME,Z,PNEW)
C
C update of surface variables for not coupled case
C     
      IF (.NOT. FL3D) THEN
         CALL VCOPYR(MAXCEL,Q_IN_KK_SN,Q_IN_KKP1_SN)
         CALL INIT0R(MAXCEL,Q_IN_KKP1_SN)
         CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_1,Q_OUT_KKP1_SN_1)
         CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_1)
         CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_2,Q_OUT_KKP1_SN_2)
         CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_2)
         IF (TRANSP) THEN
            CALL VCOPYR(MAXCEL,QMASS_IN_KK_SN,QMASS_IN_KKP1_SN)
            CALL INIT0R(MAXCEL,QMASS_IN_KKP1_SN)
            CALL VCOPYR(MAXCEL,QMASS_OUT_KK_SN_1,QMASS_OUT_KKP1_SN_1)
            CALL INIT0R(MAXCEL,QMASS_OUT_KKP1_SN_1)
            CALL VCOPYR(MAXCEL,QMASS_OUT_KK_SN_2,QMASS_OUT_KKP1_SN_2)
            CALL INIT0R(MAXCEL,QMASS_OUT_KKP1_SN_2)
         END IF
      END IF
C
C     hydrograph calculation
C
      IF (FL3D) THEN
         CALL HGRAPH(NNOD,TIME,APOT,AACT,
     1           OVFLOW,REFLOW,HGFLAG,IFATM,ATMPOT,ATMACT,PNEW,
     2           evap_eff,inf_tot)
C
C  recharge calculation
C
         CALL RECHARGE(NNOD,NSTR,WNOD,ARENOD,PNEW,TIME,DELTAT,RECFLOW,
     1              RECVOL,RECNOD)
C
C  hydrograph output
C
         WRITE(IOUT7,1190) NSTEP,DELTAT,TIME,APOT,AACT,OVFLOW,REFLOW,
     1                   SFFLW,RECFLOW,RECVOL/AREATOT
         WRITE(IOUT8,1195) NSTEP,DELTAT,TIME,NDIN+NDOUT,NNIN+NNOUT
         WRITE(IOUT30,1197) NSTEP,DELTAT,TIME,VSFFLW,VSFFLW/DELTAT
         WRITE(IOUT31,1197) NSTEP,DELTAT,TIME,VNDIN+VNDOUT,
     1      (VNDIN+VNDOUT)/DELTAT
         WRITE(IOUT32,1197) NSTEP,DELTAT,TIME,VNNIN+VNNOUT,
     1      (VNNIN+VNNOUT)/DELTAT
      
         IF (NUDN .GT. 0)  WRITE(IOUT50,1198) NSTEP,DELTAT,TIME,
     1      NUDIN,NUDOUT,NUDIN+NUDOUT,
     2      VNUDIN,VNUDOUT,VNUDIN+VNUDOUT,
     3      VNUDIN/DELTAT,VNUDOUT/DELTAT,
     4      (VNUDIN+VNUDOUT)/DELTAT,VNUDTOT
         WRITE(IOUT36,1199) NSTEP,DELTAT,TIME,VSFTOT,VNDTOT,VNNTOT,
     1      VNUDTOT,VTOT

C
C  surface vs subsurface diagnostics output
C
         CALL SAT_FRAC(NNOD,NSTR,PONDH_MIN,FHORT,FDUNN,FPOND,FSAT,
     1           PNEW)
         AACTAV=0.5D0*(AACT + AACTP)
         WRITE(IOUT43,1170) NSTEP,DELTAT,TIME,KBACKT,ITER,
     1           ITER+(ITUNS*KBACKT),NSURF,NSURFT,
     2           APOT,APOT*DELTAT,APOT/AREATOT,
     3           (APOT*DELTAT)/AREATOT,AACTAV,AACTAV*DELTAT,
     4           AACTAV/AREATOT,(AACTAV*DELTAT)/AREATOT,
     5           FHORT,FDUNN,FPOND,FSAT,CPUSUB,CPUSURF

         VAPOT_T=VAPOT_T + (APOT*DELTAT)
         VAACT_T=VAACT_T + (AACTAV*DELTAT)
         AACTP=AACT
C
C  partial output
C           
c        IF (NR .GT. 0) THEN
c             WRITE(IOUT2,1150)
c             WRITE(IOUT2,1160) (CONTR(I),PNEW(CONTR(I)),
c     1              SW(CONTR(I)),CKRW(CONTR(I)),I=1,NR)
C commento per Manoli 
c             WRITE(1001,*) TIME,(PNEW(CONTR(I)),I=1,NR)
c             WRITE(1002,*) TIME,(SW(CONTR(I)),I=1,NR)
c             WRITE(1003,*) TIME,(CKRW(CONTR(I)),I=1,NR)
c             DO I=1,nr
c               WRITE(1001,*) TIME, CONTR(I),PNEW(CONTR(I)),
c     1             SW(CONTR(I))
c             end do
c        END IF
C
C  detailed output
C
         IF (NPRT .GT. 0 .AND. KPRT .LE. NPRT) THEN
            IF (TIME .GE. TIMPRT(KPRT)) THEN
               CALL DETOUT(KPRT,IPEAT,IPRT,N,NNOD,NUMVP,NSTR,NSTEP,TIME,
     1                 NODVP,SATSUR,PNEW,INDE,DEF,SW,CKRW,UNOD,VNOD,
     2                 WNOD,NT,UU,VV,WW,X,Y,Z,OVFLNOD,atmact,arenod,
     3                 PONDNOD,ifatm,PNODI,RECNOD,QTRANIE,CNNEW)
CM          IF (NCOUT .EQ. 1) THEN
CM             CALL DETOUT_NETCDF(KPRT,NPRT,TIME,NROW,NCOL,NNOD,NSTR,
CM   1                 PNEW,SW,PONDNOD)
CM             CALL DIAGN_OUTPUT(KPRT,NPRT,TIME,N,NNOD,NSTR,BASE,ZRATIO,
CM   1                 Z,PNODI,ARENOD,VOLNOD,SW,PONDNOD,PNEW)
CM          END IF
               IF ((IPRT.GE.2).AND.(VTKF.GT.0)) THEN
                  CALL  VTKRIS3D(NT,N,100+kprt,TETRA,PNEW,SW, 
     1                           UU,VV,WW,X,Y,Z,KS,TIME,VTKF)
                  IF (TRANSP) THEN
                    CALL VTKRIS3DFC(NT,N,200+kprt,TETRA,TIME,CNEW,
     1                              SW,X,Y,Z)
                    CALL VTKRIS3DFCSURF(NT,N,100+kprt,TETRA,TIME,
     1                                  CNNEW,X,Y,Z)
                  END IF
               END IF
c              IF (NR.GT.0) THEN
c                     WRITE(1001,*) TIME,(PNEW(CONTR(I)),I=1,NR)
c                    WRITE(1002,*) TIME,(SW(CONTR(I)),I=1,NR)
c                     WRITE(1003,*) TIME,(CKRW(CONTR(I)),I=1,NR)
c              END IF
               KPRT=KPRT + 1
            END IF
         END IF
C     
C     detailed time series output of model results at the nudging
C     observation points
C     
         IF (NUDN .GT. 0) THEN
            IF (IPEAT. EQ. 0) THEN
              CALL NUDCPT(NUDFLAG,NUDN,TETRA,NUDTET,IVOL,VOLUR,SW,PNODI,
     1              PTNEW,AI,BI,CI,DI,NUDX,NUDY,NUDZ,NUDSMC)
            ELSE
               CALL NUDCPT(NUDFLAG,NUDN,TETRA,NUDTET,IVOL,VOLUR,SW,PORE,
     1              PTNEW,AI,BI,CI,DI,NUDX,NUDY,NUDZ,NUDSMC)
            END IF
            WRITE(IOUT51,1245) NSTEP,DELTAT,TIME,(NUDSMC(I),I=1,NUDN)
         END IF
      END IF
C     
      IF(FL3D) THEN
         CALL VCOPYR(ANP,QPOLD,QPNEW)
         CALL VCOPYI(NNOD,IFATMP,IFATM)
         CALL VCOPYR(NNOD,ATMOLD,ATMACT)
         CALL VCOPYR(N,PTIMEP,PNEW)
         CALL VCOPYR(N,CKRWP,CKRW)
C
C update of surface variables for coupled case
C     
         IF (SURF) THEN
            PONDP = PONDING
            CALL VCOPYR(NNOD,OVFLP,OVFLNOD)
            CALL VCOPYR(MAXCEL,Q_IN_KK_SN,Q_IN_KKP1_SN)
            CALL INIT0R(MAXCEL,Q_IN_KKP1_SN)
            CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_1,Q_OUT_KKP1_SN_1)
            CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_1)
            CALL VCOPYR(MAXCEL,Q_OUT_KK_SN_2,Q_OUT_KKP1_SN_2)
            CALL INIT0R(MAXCEL,Q_OUT_KKP1_SN_2)
            CALL VCOPYR(MAXCEL,VOLUME_KK_SN,VOLUME_KKP1_SN)
            CALL INIT0R(MAXCEL,VOLUME_KKP1_SN)
            IF (TRANSP) THEN
               CALL VCOPYR(NNOD,TRAFLP,TRAFLNOD)
               CALL VCOPYR(MAXCEL,QMASS_IN_KK_SN,QMASS_IN_KKP1_SN)
               CALL INIT0R(MAXCEL,QMASS_IN_KKP1_SN)
               CALL VCOPYR(MAXCEL,QMASS_OUT_KK_SN_1,QMASS_OUT_KKP1_SN_1)
               CALL INIT0R(MAXCEL,QMASS_OUT_KKP1_SN_1)
               CALL VCOPYR(MAXCEL,QMASS_OUT_KK_SN_2,QMASS_OUT_KKP1_SN_2)
               CALL INIT0R(MAXCEL,QMASS_OUT_KKP1_SN_2)
               CALL VCOPYR(MAXCEL,MASS_KK_SN,MASS_KKP1_SN)
               CALL INIT0R(MAXCEL,MASS_KKP1_SN)       
            END IF
            IF (NSTEP.GT.1) THEN
               CALL VCOPYR(NUMRES,H_POOL_KK_VEC,H_POOL_KKP1_VEC)
            END IF
         END IF
C 
C  determine largest, smallest, and average time step sizes
C
         CALL DTSTAT(TIME,DELTAT,DTBIG,TBIG,DTSMAL,TSMAL,DTAVG)
C
C  if convergence failed and no back-stepping is 
C  possible (case (iv)) we end the simulation. Otherwise we 
C  have case (ii) (convergence achieved) and we prepare for
C  the next time step (unless we've reached TMAX or we're running
C  a steady state problem). 
C
         IF (NOBACK) GO TO 300
C
C  update TIME and time step size for the next time step.
C  If TMAX reached, end the simulation.       
c
c update time
c
         IF (DELTAT .GE. 1.0D+15) THEN
             TIME=TMAX+1
             GO TO 200
         END IF
         IF (DABS(TIME-TMAX) .LE. 0.001D0*DELTAT) THEN
             TIME=TMAX+1
             GO TO 200
         END IF    
         CALL TIMUPD(TIMEP,TIME,DELTAT,DTMIN,DTMAX,TMAX,
     1               DTMAGA,DTMAGM,DTREDS,DTREDM,DTGMIN,
     2               ITER,ITUNS1,ITUNS2)
 
C
C  update or re-initialize counters for next time level
C       
         CALL TIMNXT(NSTEP,ITER,NITERT,KBACKT,NSURFT,CPUSUB,CPUSURF)
      ELSE IF (.NOT. FL3D) THEN
         IF (DABS(TIME-TMAX) .LE. 0.001D0*DELTATS) GO TO 300
         CALL TIMUPDSUP(TIME,DELTATS,TMAX)
         IF ((TIME+DELTATS) .GT. TMAX) THEN
           DELTATS=TMAX - TIME
           TIME=TMAX
         ELSE
            TIME=TIME + DELTATS
         END IF
         CALL EFFNXT(NCELNL,HSPEFF,HTIEFF,DELTA_X,DELTA_Y,TIME,EFFTIM,
     1        QOI_SN,SURFACE_WATER_INP,SURFACE_WATER_SN)
         NSTEP= NSTEP + 1
      END IF
  200 CONTINUE 
      END DO
C 
C END of THE INTERIOR TIME LOOP
C
C-------------------------------END OF SIMULATION, OUTPUT FINAL SOLUTION
C
  300 CONTINUE
      TIME=TMAX
C
C  output final results
C
      IF (FL3D) THEN
         CALL DETOUT(KPRT,IPEAT,IPRT,N,NNOD,NUMVP,NSTR,NSTEP,TIME,NODVP,
     1        SATSUR,PNEW,INDE,DEF,SW,CKRW,UNOD,VNOD,WNOD,NT,
     2        UU,VV,WW,X,Y,Z,OVFLNOD,atmact,arenod,PONDNOD,ifatm,
     3        PNODI,RECNOD,QTRANIE,CNNEW)
         IF ((IPRT.GE.2).AND.(VTKF.GT.0)) THEN
            CALL VTKRIS3D(NT,N,100+kprt, TETRA,PNEW,SW, 
     1                    UU,VV,WW,X,Y,Z,KS,TIME,VTKF)
            IF (TRANSP) THEN
               CALL VTKRIS3DFC(NT,N,200+kprt,TETRA,TIME,CNEW,
     1                         SW,X,Y,Z)
               CALL VTKRIS3DFCSURF(NT,N,100+kprt,TETRA,TIME,
     1                             CNNEW,X,Y,Z)
            END IF
         END IF
c        IF (NR.GT.0) THEN
c               WRITE(1001,*) TIME,(PNEW(CONTR(I)),I=1,NR)
c              WRITE(1002,*) TIME,(SW(CONTR(I)),I=1,NR)
c               WRITE(1003,*) TIME,(CKRW(CONTR(I)),I=1,NR)
c        END IF

         WRITE(IOUT43,1175) KBACK,ITRTOT,ITRTOT+(ITUNS*KBACK),
     1        NSURFT_T,NSURFT_TB,VAPOT_T,VAPOT_T/AREATOT,
     2        VAACT_T,VAACT_T/AREATOT,CPUSUB_T,CPUSURF_T
C
C  output cumulative flow volumes and mass balance errors, 
C  and flag and failure summaries
C
         WRITE(IOUT2,1200) VSFTOT,VNDTOT,VNNTOT,VNUDTOT,VTOT
         IF (DELTAT .LT. 1.0D+15) THEN
            IF (VTOT .NE. 0.0D0) WRITE(IOUT2,1210) CAERAS,
     1           100.0D0*CAERAS/VTOT
            IF (VTOT .NE. 0.0D0) WRITE(IOUT2,1215) CERRAS,
     1           100.0D0*CERRAS/VTOT
         END IF
         WRITE(IOUT2,1260) (HGFLAG(I),I=1,9)
         WRITE(IOUT9,1260) (HGFLAG(I),I=1,9)
         IF (NSF .GT. 0) THEN
            WRITE(IOUT2,1270) (SFFLAG(I),I=1,5)
c            WRITE(IOUT10,1270) (SFFLAG(I),I=1,5)
         END IF
         WRITE(IOUT2,1250) KBACK,KLSFAI
         IF (NSF .GT. 0) WRITE(IOUT2,1255) KSFCV,KSFCVT
C
C  output cpu times for code sections and total cpu time 
C  for the simulation
C
         IF (IOPT .NE. 1 .AND. ISOLV .EQ. 3) WRITE(IOUT2,1310) MINBOT
         CALL TIM(CPUMN,2)
         AVGNL=FLOAT(ITRTOT)/NSTEP
         AVGLIN=FLOAT(ITLIN)/NSTEP
         AVGLNL=FLOAT(ITLIN)/ITRTOT
         ATCTS=CPUMN/NSTEP
         ATCNL=CPUMN/ITRTOT
         ANCTS=CPUNL/NSTEP
         ANCNL=CPUNL/ITRTOT
         PCMN=100.0
         PCNL=CPUNL/CPUMN*100.0
         PCVEC1=CPUVEC(1)/CPUMN*100.0
         PCVEC2=CPUVEC(2)/CPUMN*100.0
         PCVEC3=CPUVEC(3)/CPUMN*100.0
         PCVEC4=CPUVEC(4)/CPUMN*100.0
         PCVEC5=CPUVEC(5)/CPUMN*100.0
         PCVEC6=CPUVEC(6)/CPUMN*100.0
         PCVEC7=CPUVEC(7)/CPUMN*100.0
         PCVEC8=CPUVEC(8)/CPUMN*100.0
         PCVEC9=CPUVEC(9)/CPUMN*100.0
         CPUOVH=CPUMN-CPUNL
         PCOVH=CPUOVH/CPUMN*100.0
         WRITE(IOUT2,1400) CPUMN,PCMN,CPUNL,PCNL,CPUOVH,PCOVH
         WRITE(IOUT2,1420) CPUVEC(1),PCVEC1,CPUVEC(2),PCVEC2,
     1        CPUVEC(3),PCVEC3,CPUVEC(4),PCVEC4,
     2        CPUVEC(5),PCVEC5,CPUVEC(6),PCVEC6,
     3        CPUVEC(7),PCVEC7,CPUVEC(8),PCVEC8,
     4        CPUVEC(9),PCVEC9
         CPUNLT=0.0
         DO I=1,11
            CPUNLT=CPUNLT+CPUVEC(I)
         END DO
         PCNLT=CPUNLT/CPUMN*100.0
         WRITE(IOUT2,1440) CPUNLT,PCNLT
         WRITE(IOUT2,1450) NSTEP,DTSMAL,TSMAL,DTBIG,TBIG,
     1                     DTAVG/DFLOAT(NSTEP),ITRTOT,ITLIN,
     2                     AVGNL,AVGLIN,AVGLNL,ATCTS,ATCNL,ANCTS,ANCNL
C
C     CALL FLUSH(IOUT2)
C  
      END IF
      CALL CLOSIO
      STOP
C
C------------------------------------------------------FORMAT STATEMENTS
C
  666 FORMAT(1000(1PE15.6))
  667 FORMAT(10(1x,I5))
  668 FORMAT(5E15.7)
c 1000 FORMAT(I7,1PE16.8,'     NSTEP   TIME')
c 2000 FORMAT(' SURFACE NODE              X              Y',
c     1     '  PRESSURE HEAD')
c 2060 FORMAT(7X,I8,3(1PE15.6),i6,1pe15.6)
 1005 FORMAT(/,5X,'NTERM (# OF NONZERO TERMS IN SYS. MAT.) = ',I8)
 1010 FORMAT(//,'      TIME STEP: ',I8,'    DELTAT: ',1PE12.4,
     1              '    TIME: ',1PE12.4,
     2        /,'     *********************************************',
     3          '*********************',/)
 1060 FORMAT(  5X,'IOPT   (1 PICARD, 2 NEWTON)             = ',I6,
     1       /,5X,'NLRELX (0 NORELX,1 CONS RELX,2 VAR RELX)= ',I6,
     2       /,5X,'KSLOPE (0 ANA, 1 KSL/ANA, 2 KSL/C-DIFF,',
     3       /,5X,'        3 LOC KSL/ANA, 4 LOC TAN-SLOPE) = ',I6,
     4       /,5X,'TOLUNS (TOLERANCE FOR NONLINEAR ITER)   = ',1PE15.5,
     5       /,5X,'TOLSWI (TOLERANCE FOR BC SWITCHING)     = ',1PE15.5,
     6       /,5X,'L2NORM (0 INFINITY NORM, ELSE L2 NORM)  = ',I6,
     7      //,' nlinr  linr converg error norms  node   PNEW at',
     8         '   POLD at   resid error norms',
     9       /,'  iter  iter       PL2      PINF IKMAX     IKMAX',
     A         '     IKMAX       FL2      FINF',
     B       /,' ===============================================',
     C         '==============================')
 1065 FORMAT(23X,' (NSTEP: ',I5,'  DELTAT: ',1PE12.4,'  TIME: ',
     1         1PE12.4,')')
 1130 FORMAT(/,'  INFLOW (I) AND OUTFLOW (O) FROM',
     1         ' ATM (A) AND NON-ATM, NON-SEEP FACE (N) BC''S;',
     2       /,'  ''C F'' CURRENT FLUX; ''P F'' PREVIOUS FLUX;',
     3         ' ''VOL'' VOLUME',
     4       /,'      IA DIRIC OA DIRIC IN DIRIC ON DIRIC',
     5         ' IA NEUMN OA NEUMN IN NEUMN ON NEUMN',
     6       /,1X,'C F',1X,8(1PE9.1),
     7       /,1X,'P F',1X,8(1PE9.1),
     8       /,1X,'VOL',1X,8(1PE9.1))
 1135 FORMAT(/,' SUBSURFACE OUTFLOW AT SEEP. FACE BC''S:',
     1         ' C F=',1PE8.1,' P F=',1PE8.1,
     2         ' VOL=',1PE8.1)
 1140 FORMAT(/,' TOTAL I VOL  TOTAL O VOL  STOR CHNG VOL',
     1         ' | ABS MASS BAL ERR   REL MASS BAL ERR,%',
     2       /,1X,1PE11.3,2X,1PE11.3,2X,1PE13.5,' | ',
     3         1PE16.8,3X,1PE18.10)
 1150 FORMAT(//,1X,'SOLUZIONE SUI NODI DI OUTPUT',/,
     1          1X,2(' NODE  PRESSIONE  SW         CKRW      '))
 1160 FORMAT(2(1X,I5,3(1PE11.3)))
 1170 FORMAT(I7,2(1PE10.3),5(I6),8(1PE11.3),4(0PF7.3),2(1PE11.3))
 1175 FORMAT('#',192('='),
     1       /,'#Total:',20X,5(I6),4(11X,1PE11.3),28X,2(1PE11.3))
 1190 FORMAT(I10,2(1PE13.5),7(1PE13.5))
 1195 FORMAT(I6,2(1PE11.3),2(11X,1PE13.5))
 1197 FORMAT(I7,4(4X,1PE17.9))
 1198 FORMAT(I6,12(1PE10.2))
 1199 FORMAT(I7,6(1PE10.2),1PE12.4)
 1200 FORMAT(/////,
     1          '  CUMULATIVE SEEPAGE FACE FLOW VOLUME     ',
     2          '       : ',1PE12.5,
     3        /,'  CUMULATIVE NON-ATM, NON-SF DIRICH FLOW V',
     4          'OLUME  : ',1PE12.5,
     5        /,'  CUMULATIVE NON-ATM, NON-SF NEUMANN FLOW ',
     6          'VOLUME : ',1PE12.5,
     7        /,'  CUMULATIVE NUDGING TERM "FLOW" VOLUME   ',
     8          '       : ',1PE12.5,
     9        /,'  CUMULATIVE TOTAL NET FLOW VOLUME (VIN + ',
     A          'VOUT)  : ',1PE12.5)
 1210 FORMAT(//,'  CUMULATIVE (TOTAL) MBE BASED ON CAERAS  ',
     1          '       : ',1PE12.5,
     2        /,'  CUMULATIVE (TOTAL) RELATIVE MBE BASED ON',
     3          ' CAERAS: ',1PE12.5,' %')
 1215 FORMAT( /,'  CUMULATIVE (TOTAL) MBE BASED ON CERRAS  ',
     1          '       : ',1PE12.5,
     2        /,'  CUMULATIVE (TOTAL) RELATIVE MBE BASED ON',
     3          ' CERRAS: ',1PE12.5,' %')
 1220 FORMAT('#INITIAL VOLUME OF WATER IN THE SUBSURFACE = ',1PE13.5,
     1     /,'#NSTEP       DELTAT         TIME  NLIN   AVG.LIN',
     2       '       STORE1       STORE2       DSTORE',
     3       '   CUM.DSTORE       VIN  CUM.VIN       VOUT',
     4       '    CUM.VOUT   VIN+VOUT   CUM. VIN  M.BAL.ERR',
     5       '   REL. MBE   CUM. MBE      CUM.',
     6     /,'#                                 ITER.    ITER.',
     7       '                             (>0 x inc)',
     8       '                                            ',
     9       '                          + VOUT ',
     A       ' (Vi+Vo-Ds)        (%)             ABS(MBE)')
 1222 FORMAT('#INITIAL MASS OF SOLUTE IN THE SUBSURFACE = ',1PE13.5,
     1     /,'#NSTEP       DELTAT         TIME  NADV   ',
     2       '      MASST        DSMASS     MASSINADV  ',
     3       '  MASSINDISP    MASSOUTADV   MASSOUTDISP    ',
     4       '    MINTOT    MASSOUTTOT        MERRAS       CMERRAS',
     5       '     MASSSURF')
 1223 FORMAT('#INITIAL MASS OF SOLUTE IN THE SUBSURFACE = ',1PE13.5)
 1225 FORMAT(/,5X,'INITIAL VOLUME OF WATER IN SUBSURFACE   = ',1PE13.5)
 1230 FORMAT('#NSTEP    DELTAT      TIME       ',
     1       ' <--- NUDSMC(I), I=1,2,...,NUDN --->')
 1232 FORMAT('#TIME      <--- WTDEPTH(NODVP(I)), I=1,2,...,NUMVP --->')
 1240 FORMAT(I6,2(1PE13.6),I6,1PE10.3,2(1PE13.6),2(1PE13.5),2(1PE10.3),
     1       7(1PE11.3),1PE10.3)
 1242 FORMAT(I6,2(1PE13.6),I6,12(1PE14.5E3))
 1245 FORMAT(I6,12(1PE13.5))
 1250 FORMAT( /,'  TOTAL NUMBER OF BACK-STEPPING OCCURRENCES',
     1          '      : ',I6,
     2        /,'  TOTAL NUMBER OF LINEAR SOLVER FAILURES   ',
     3          '      : ',I6)
 1255 FORMAT( /,'  TOTAL NUMBER OF S-F EX-PT CONV FAIL OCCUR',
     1          'RENCES: ',I6,
     2        /,'  TOTAL NUMBER OF S-F EX-PT CONVERGENCE FAI',
     3          'LURES : ',I6)
 1260 FORMAT(/,' HGFLAG:    (1)    (2)    (3)    (4)    (5)',  
     +                 '    (6)    (7)    (8)    (9)',
     +       /,8X,9(I6,1X))
 1270 FORMAT(/,' SFFLAG:    (1)    (2)    (3)    (4)    (5)',  
     +       /,8X,5(I6,1X))
 1310 FORMAT(//,1X,' MINIMUM IBOT VALUE FOR NONSYM SOLVER = ',I10)
 1400 FORMAT(//,1X,' TOT CPU FOR THE SIMULATION           = ',0PF15.2,
     1             ' SECONDS  ',0PF6.2,' %',/,
     2          1X,' TOT CPU FOR NONLINEAR SCHEME         = ',0PF15.2,
     3             ' SECONDS  ',0PF6.2,' %',/,
     4          1X,' TOT CPU FOR OVERHEAD (SEE CODE DESC.)= ',0PF15.2,
     5             ' SECONDS  ',0PF6.2,' %')
 1420 FORMAT( /,1X,' TOT CPU FOR UNSAT CHARACTERISTICS    = ',0PF15.2,
     1             ' SECONDS  ',0PF6.2,' %',/,
     2          1X,' TOT CPU FOR INIT. OF SYSTEM MATRICES = ',0PF15.2,
     3             ' SECONDS  ',0PF6.2,' %',/,
     4          1X,' TOT CPU FOR LOCAL SYSTEM ASSEMBLY    = ',0PF15.2,
     5             ' SECONDS  ',0PF6.2,' %',/,
     6          1X,' TOT CPU FOR RHS CALCULATION W/O BC''S = ',0PF15.2,
     7             ' SECONDS  ',0PF6.2,' %',/,
     8          1X,' TOT CPU FOR GLOBAL LHS SYSTEM MATRIX = ',0PF15.2,
     9             ' SECONDS  ',0PF6.2,' %',/,
     A          1X,' TOT CPU FOR BC CONTRIBUTIONS TO RHS  = ',0PF15.2,
     B             ' SECONDS  ',0PF6.2,' %',/,
     C          1X,' TOT CPU FOR LINEAR SOLVER & RESIDUAL = ',0PF15.2,
     D             ' SECONDS  ',0PF6.2,' %',/,
     E          1X,' TOT CPU FOR DIR. RESETTING AND DIFMAX= ',0PF15.2,
     F             ' SECONDS  ',0PF6.2,' %',/,
     G          1X,' TOT CPU FOR B-FLX, M-BAL, SW, & EX-PT= ',0PF15.2,
     H             ' SECONDS  ',0PF6.2,' %')
 1440 FORMAT(41X,'-----------',9X,'-------',/,41X,0PF11.2,
     1             ' SECONDS  ',0PF6.2,' %')
 1450 FORMAT( /,1X,' TOTAL NUMBER OF TIME STEPS           = ',I8,/,
     1          1X,' SMALLEST TIME STEP SIZE              = ',1PE11.3,
     2             '  AT TIME ',1PE11.3,/,
     3          1X,' LARGEST  TIME STEP SIZE              = ',1PE11.3,
     4             '  AT TIME ',1PE11.3,/,
     5          1X,' AVERAGE  TIME STEP SIZE              = ',1PE11.3,/,
     6          1X,' TOTAL NUMBER OF NONLINEAR ITERATIONS = ',I8,/,
     7          1X,' TOTAL NUMBER OF LINEAR ITERATIONS    = ',I8,
     8             '  (= 0 FOR NONSYM)',/,
     9          1X,' AVG NL ITERATIONS PER TIME STEP      = ',0PF8.2,/,
     A          1X,' AVG LINEAR ITERATIONS PER TIME STEP  = ',0PF8.2,/,
     B          1X,' AVG LIN ITERATIONS PER NL ITERATION  = ',0PF8.2,/,
     C          1X,' AVG TOT CPU PER TIME STEP            = ',0PF8.2,
     D             ' SECONDS',/,
     E          1X,' AVG TOT CPU PER NL ITERATION         = ',0PF8.2,
     F             ' SECONDS',/,
     G          1X,' AVG TOT NL CPU PER TIME STEP         = ',0PF8.2,
     H             ' SECONDS',/,
     I          1X,' AVG TOT NL CPU PER NL ITERATION      = ',0PF8.2,
     J             ' SECONDS')
 1460 FORMAT( /,1X,'NENS  TOTAL NUMBER OF TIME STEPS   = ',I8,/,
     1          1X,'NENS0 TOTAL NUMBER OF TIME STEPS   = ',I8,/,
     1          1X,' SMALLEST NUMBER OF TIME STEP      = ',I8,
     2             '  OF THE REALIZATION ',I8,/,
     3          1X,' LARGEST  NUMBER OF TIME STEP      = ',I8,
     4             '  OF THE REALIZATION ',I8,/,
     5          1X,' NENS AVERAGE  TIME STEP SIZE      = ',1PE11.3,/,
     6          1X,' NENS TOT NUMBER OF NONLINEAR ITERATIONS = ',I8,/,
     6          1X,'NENS0 TOT NUMBER OF NONLINEAR ITERATIONS = ',I8,/,
     7          1X,' TOT NUMBER OF LINEAR ITERATIONS   = ',I8,
     8             '  (= 0 FOR NONSYM)',/,
     9          1X,' AVG NL ITERATIONS PER TIME STEP   = ',0PF8.2,/,
     A          1X,' AVG LINEAR ITERATIONS PER TIME STEP  = ',0PF8.2,/,
     B          1X,' AVG LIN ITERATIONS PER NL ITERATION  = ',0PF8.2,/,
     C          1X,' AVG TOT CPU PER TIME STEP (NENS)     = ',0PF8.2,
     D             ' SECONDS',/,
     C          1X,' AVG TOT CPU PER TIME STEP (NENS0)    = ',0PF8.2,
     D             ' SECONDS',/,
     E          1X,' AVG TOT CPU PER NL ITERATION         = ',0PF8.2,
     J             ' SECONDS',/,
     K          1X,' TOT CPU UPDATES                      = ',0PF8.2,
     L             ' SECONDS',/,
     M          1X,' AVG CPU UPDATE PER REALIZATION       = ',0PF8.2,
     N             ' SECONDS',/,
     P          1X,' FINAL ENSEMBLE SIZE                  = ',I8,/,
     Q          1X,' TOT CPU TIME                         = ',0PF10.2)
 1510 FORMAT('  0   HSPVEL')
 1520 FORMAT('  0   HSPSW')
 1530 FORMAT('#',I8)
 1540 FORMAT('#          TIME',1X,5000I16)
 1541 FORMAT('# TIME  ',1X,5000I16)
 1542 FORMAT(1PE16.8,5000E16.7)
 1543 FORMAT('#VOLUME, SOLUTE MASS, ADSORBED MASS')
c 1560 FORMAT(1X,1PE15.6)
c 1570 FORMAT(1X,'NROW = ',I8,1X,'NCOL = ',I8)
 1600 FORMAT(//,5X,
     1   'SATURATED HYDRAULIC CONDUCTIVITY, SPECIFIC STORAGE, AND ',
     2   'POROSITY VALUES',/,
     3   1X,' LAYER MAT.TYPE  X-PERM       Y-PERM       Z-PERM',
     4      '       STORAGE      POROSITY')
 1610 FORMAT(1X,I4,I8,2X,5(1PE13.5))
 1630 FORMAT(  5X,'VGN                                     = ',1PE15.5,
     1       /,5X,'VGRMC                                   = ',1PE15.5,
     2       /,5X,'VGPSAT                                  = ',1PE15.5)
      END
