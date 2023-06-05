C
C**************************  COMMON BLOCK IOUNIT ***********************
C
C  input and output file units.
C  Values are assigned in BLOCK DATA subprogram block_data.f.
C  Files are opened and closed in subroutines openio.f and closio.f. 
C
C***********************************************************************
C
      INTEGER    IFN
      INTEGER    IIN1,  IIN2,  IIN3, IIN4,  IIN5,  IIN6
      INTEGER    IIN7,  IIN8,  IIN9
      INTEGER    IIN10, IIN11
      INTEGER    IIN16, IIN17, IIN18, IIN20, IIN21
      INTEGER    IIN22, IIN23, IIN25, IIN26
      INTEGER    IIN27, IIN28, IIN29, IIN30, IIN31
      INTEGER    IIN32, IIN33, IIN34, IIN35, IIN36
      INTEGER    IIN37, IIN38, IIN39
      INTEGER    IIN50, IIN51, IIN60
      INTEGER    IIN61, IIN62, IIN63, IIN64 
      INTEGER    ITERM
      INTEGER    IOUT1, IOUT2, IOUT3, IOUT4, IOUT5, IOUT6
      INTEGER    IOUT7, IOUT8, IOUT9, IOUT10,IOUT11,IOUT12
      INTEGER    IOUT13,IOUT14,IOUT15,IOUT16,IOUT17,IOUT18
      INTEGER    IOUT19,IOUT20,IOUT30,IOUT31,IOUT32,IOUT36
      INTEGER    IOUT40,IOUT41,IOUT42,IOUT43,IOUT44,IOUT50
      INTEGER    IOUT51,IOUT57,IOUT80,IOUT81,IOUTPT
c
      COMMON    /IOUNIT/ IFN,
     1                   IIN1,  IIN2,  IIN3, IIN4,  IIN5,  IIN6,
     2                   IIN7,  IIN8,  IIN9,
     3                   IIN10, IIN11, 
     4                   IIN16, IIN17, IIN18, IIN20, 
     5                   IIN21, IIN22, IIN23, IIN25,
     6                   IIN26, IIN27, IIN28, IIN29, IIN30,
     7                   IIN31, IIN32, IIN33, IIN34, IIN35,
     8                   IIN36, IIN37, IIN38, IIN39,
     9                   IIN50, IIN51, IIN60, IIN61, IIN62, IIN63, 
     A                   IIN64, ITERM,
     B                   IOUT1, IOUT2, IOUT3, IOUT4, IOUT5, IOUT6,
     C                   IOUT7, IOUT8, IOUT9, IOUT10,IOUT11,IOUT12,
     D                   IOUT13,IOUT14,IOUT15,IOUT16,IOUT17,IOUT18,
     E                   IOUT19,IOUT20,IOUT30,IOUT31,IOUT32,IOUT36,
     F                   IOUT40,IOUT41,IOUT42,IOUT43,IOUT44,IOUT50,
     G                   IOUT51,IOUT57,IOUT80,IOUT81,IOUTPT
			     
