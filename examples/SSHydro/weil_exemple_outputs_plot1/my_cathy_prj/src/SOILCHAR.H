C
C**************************  COMMON BLOCK SOILCHAR *********************
C            
C  real scalars and real arrays for unsaturated soil characteristics
C  and for chord and tangent slope formulas
C            
C***********************************************************************
C
      REAL*8  PMIN
      REAL*8  VGN(NMAX),VGM(NMAX)
      REAL*8  VGRMC(NMAX),VGPSAT(NMAX)
      REAL*8  VGN1(NMAX),VGM52(NMAX)
      REAL*8  VGNR(NMAX),VGMM1(NMAX)
      REAL*8  VGMR(NMAX),VGPSN(NMAX)
      REAL*8  CBETA0,CANG
      REAL*8  HUN,HUA,HUB,HUALFA,HUBETA,HUGAMA,HUPSIA,HUSWR
      REAL*8  HUSWR1,HUALB,HUGAM1,HUGB,HUGB1,HU1BET,HUGAM2
      REAL*8  HUN1,HU2A,HUB2A,HUAB,HULN10
      REAL*8  BCBETA,BCRMC,BCPSAT
      REAL*8  BCB1,BCB2,BCBB1P,BCBPS,BC23B,BC23BP,BC33B
      REAL*8  TOLKSL,PKRL,PKRR,PSEL,PSER
      REAL*8  PDSE1L,PDSE1R,PDSE2L,PDSE2R,DKRTAN
      REAL*8  VGPNOT(NMAX),BCPORM(NMAX)
      REAL*8  DSETAN(NMAX),DDSE1T(NMAX),DDSE2T(NMAX)
      REAL*8  PCAP(MAXZON,MAXSTR,MAXLKP),SATC(MAXZON,MAXSTR,MAXLKP)
      REAL*8  KRWC(MAXZON,MAXSTR,MAXLKP)
      REAL*8  SCF,PCREF(MAXVEG),PCWLT(MAXVEG),PCANA(MAXVEG)
      REAL*8  PZ(MAXVEG),OMGC(MAXVEG),ZROOT(MAXVEG)
C
C  common blocks for van Genuchten, extended van Genuchten,
C  Huyakorn, and Brooks-Corey moisture curve parameters
C
      COMMON  /CAMPO/CBETA0,CANG
      COMMON  /VG/VGPNOT,VGN,VGM,VGRMC,VGPSAT,
     1            VGN1,VGM52,VGNR,VGMM1,VGMR,VGPSN
      COMMON  /HU/HUN,HUA,HUB,HUALFA,HUBETA,HUGAMA,HUPSIA,HUSWR,
     1            HUSWR1,HUALB,HUGAM1,HUGB,HUGB1,HU1BET,HUGAM2,
     2            HUN1,HU2A,HUB2A,HUAB,HULN10
      COMMON  /BC/BCPORM,BCBETA,BCRMC,BCPSAT,
     1            BCB1,BCB2,BCBB1P,BCBPS,BC23B,BC23BP,BC33B
      COMMON  /TBLKUP/PCAP,SATC,KRWC
C
C  common block for PMIN and for chord and tangent slope formulas
C
      COMMON  /CHORD/PMIN,
     1               TOLKSL,PKRL,PKRR,PSEL,PSER,
     2               PDSE1L,PDSE1R,PDSE2L,PDSE2R,DKRTAN,
     3               DSETAN,DDSE1T,DDSE2T
C
C  common block for root water uptake parameters
C
      COMMON  /FEDDES/SCF,PCREF,PCWLT,PCANA,PZ,OMGC,ZROOT
