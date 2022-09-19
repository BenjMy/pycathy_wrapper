@echo off
   
   :RICHIESTE
   echo.
   echo HAP version 11.5.1
   echo.
   echo (1) wpar      write parameter file
   echo (2) wbb       write dtm_13.val in basin_b 
   echo (3) cbb       write dtm_*.val in basin_b
   echo (4) rn        river network
   echo (5) rrbb      read record basin_b
   echo (6) rbb       read basin_b
   echo (7) mrbb      multiple rbb
   echo (8) bb2shp    basin_b to shape
   echo (9) cat_del   catchment delineation
   echo.
   echo For CATHY simulations is enough only:
   echo.
   echo (10) cppp      CATHY preprocessing program
   echo.
   set choice=
   set /p choice=Select the program to compile 
   echo.
   
   if '%choice%'=='1' goto UNO
   if '%choice%'=='2' goto DUE
   if '%choice%'=='3' goto TRE
   if '%choice%'=='4' goto QUATTRO
   if '%choice%'=='5' goto CINQUE
   if '%choice%'=='6' goto SEI
   if '%choice%'=='7' goto SETTE
   if '%choice%'=='8' goto OTTO
   if '%choice%'=='9' goto NOVE
   if '%choice%'=='10' goto DIECI
   if '%choice%'=='0' goto END
   
   :UNO
   gfortran wpar.f90 mpar.f90
   goto END
   
   :DUE
   gfortran wbb.f90 wbb_sr.f90 mpar.f90 mbbio.f90
   goto END
   
   :TRE
   gfortran cbb.f90 mpar.f90 mbbio.f90
   goto END
   
   :QUATTRO
   gfortran rn.f90 mpar.f90 mbbio.f90 csort.f90 qsort.f90 depit.f90 cca.f90 smean.f90 dsf.f90 facet.f90 hg.f90
   goto END
      
   :CINQUE
   gfortran rrbb.f90 mpar.f90 mbbio.f90
   goto END
   
   :SEI
   gfortran rbb.f90 mpar.f90 mbbio.f90
   goto END
   
   :SETTE
   gfortran mrbb.f90 mrbb_sr.f90 mpar.f90 mbbio.f90
   goto END
   
   :OTTO
   gfortran bb2shp.f90 bb2shp_sr.f90 mpar.f90 mbbio.f90 shape.f90 dbase.f90 streamer.f90
   goto END

   :NOVE
   gfortran cat_del.f90 mpar.f90 mbbio.f90
   goto END
      
   :DIECI
   gfortran cppp.f90 mpar.f90 mbbio.f90 wbb_sr.f90 csort.f90 qsort.f90 depit.f90 cca.f90 smean.f90 dsf.f90 facet.f90 hg.f90 mrbb_sr.f90 bb2shp_sr.f90 shape.f90 dbase.f90 streamer.f90
   goto END
  
     
   :END
   echo.
@echo on
