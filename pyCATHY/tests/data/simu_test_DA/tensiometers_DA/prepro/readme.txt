
============
NRW_DTM (North Rhine-Westphalia domain) (example for delineating a DEM)
1st pass: Boundary channel construction (No:0,Yes:1) = 1
          Catchment delineation construction (No:0,Yes:1) = 0
ncview CATHY_prepro.nc
i,j = 90,281
(use DEM, MAGX, bilin, etc buttons)
2nd pass: Catchment delineation construction (No:0,Yes:1) = 1
          Outlet cell i-coordinate =                                     90
          Outlet cell j-coordinate =                                     281
Select the header type: 2
Select the nodata value: 0
Select the pointer system: 1
ncview CATHY_prepro_cut.nc

DA_720 (des Anglais) (already cut because 720 m too coarse to do so)

============

Note that the units for the pre-processor are pre-set to meters and seconds
(m, s)! Therefore when using DEM input for the simulator we must by default
use m, s!

From codes/preprocessor:
  - ./compgfortran
    (This is our "Makefile" for now. An actual Makefile is still to be
    completed. With the current Makefile, "make clean" for example is
    already set up - it will remove all executables, etc except "cppp"
    and "cat_del".)
    Note that it is usually necessary to run ./compgfortran twice.
    ./compgfortran produces "cppp" (Cathy PreProcessing Program executable).

From datasets/desAnglais_360/preprocessor, and starting with a "full" DEM:
(or from datasets/simple-test/preprocessor, or any other dataset directory)
  - Input files are "hap.in" and "dtm_13.val".
  - "cp dtm_13.val dtm_original" (or any other name) to save the original
    "full" DEM.
  - In "hap.in" set "Boundary channel construction" to 1.
    (For "Coefficient for boundary channel elevation definition" and
    "Coefficient for outlet cell elevation definition", the values 0.50
    should be OK. These parameters are used to ensure that the boundary channel
    is able to drain the entire DEM, and that a single cell within this
    boundary channel ends up being the outlet cell for the "full" DEM.)
    Note: spelling correction is needed in subroutines "mpar" and wbb_sr"
    ("constraction" -> "construction").
  - ../../../codes/preprocessor/cppp
    At the prompts from "cppp", select "GRASS" as the "header type", "0" as
    the "nodata value", and "HAP" as the "pointer system".
    This first step determines the flow directions for all cells of the DEM.
    To do so it will add a "boundary channel" around the DEM.
    Amongst the many output files produced by "cppp" there is also the
    netcdf-formatted file "CATHY_prepro.nc" that Mauro added in May/12.
    This file can be visualized with "ncview CATHY_prepro.nc" and it can be
    viewed as an ascii or text file using "ncdump CATHY_prepro.nc | more".
  - "cp CATHY_prepro.nc CATHY_prepro_original.nc" (or any other name) to
    save the original "full" DEM CDF-formatted file.
  - Determine the coordinates (i,j) of the catchment outlet cell. Can do
    this by visualizing variable "HCID" in "ncview" and comparing this map
    for example with an actual contour map showing stream channels, etc.
  - ../../../codes/preprocessor/cat_del
    At the prompt from "cat_del", specify the (i,j) coordinates of the
    catchment outlet cell.
    This second step performs the catchment delineation.
  - ../../../codes/preprocessor/mask_dem
    This third step produces the correct DEM input file for the final
    "cppp" step.
    The output from this step is a *new* "dtm_13.val" DEM file that contains
    elevation values only for the cells that belong to the catchment (and
    -9999 values for all cells outside the catchment). I.e., "dtm_13.val" is
    now the "catchment" DEM file (as opposed to the "full" DEM file).
  - In "hap.in" set "Boundary channel construction" to 0.
    [in Mauro's new version this parameter should always be 1!]
  - ../../../codes/preprocessor/cppp
    It is now necessary to run "cppp" again, this time on the "catchment" DEM
    file. 
    This final step determines the flow directions for the catchment cells,
    as well as all the other files that are needed as input for the CATHY
    simulator. There is also a new "CATHY_prepro.nc" CDF-formatted file for
    the "catchment" DEM variables.

If we start with a DEM file ("dtm_13.val") for an already delineated
catchment (i.e., a "catchment" DEM file instead of a "full" DEM file), then
only the last step in the above procedure (i.e., run "cppp" just once) is
needed (make sure that "Boundary channel construction" is set to 0 in
"hap.in").

Parameter settings in hap.in:
  - For the sample DEM  2  4  8  12
                        1  3  6  10
                        9  5  4   6
    we have the following "structural parameters":
    DEM rectangle size along the x-direction = 4 (this is the coordinate "i")
    DEM rectangle size along the y-direction = 3 (this is the coordinate "j")
    Number of cells within the catchment = (# rows x # columns for a "full" DEM,
      arbitrary for a "catchment" DEM - the "cppp" code will determine this
      value)
    X low left corner coordinate = (x coordinate of cell with elevation 9 above)
    Y low left corner coordinate = (y coordinate of cell with elevation 9 above)
    Note that (X low left corner, Y low left corner) is the (i=1, j=1) cell.
    For the des Anglais dataset, the "X low left corner" and "Y low left corner"
    values specified in hap.in are in the UTM coordinate system.
  - "terrain analysis parameters":
    Drainage directions method (LAD:1,LTD:2) = 2
    Upstream deviation memory factor (CBM:0,PBM:1) = 0
    Threshold on the contour curvature (NDM:-1E10;DM:+1E10) = -0.100E+11
    Nondispersive channel flow (0:not-required;1:required) = 1
    (The above 4 values correspond to a "classic", nondispersive D8 algorithm
    for drainage direction determination.)
    Channel initiation method (A:1,AS**k:2,ND:3) = 1 (this is the Montgomery
      and Foufoula-Georgiou (1993) algorithm)
    Threshold on the support area (A) = (for "1" chosen above, this threshold
      area is best determined by trial-and-error, comparing the resulting
      stream vs hillslope drainage network with a map showing the actual
      stream network, for example)
  - "rivulet network parameters":
    Rivulet spacing = (e.g., "10" x a 360 m DEM implies 36 rivulets / cell)
    Reference drainage area (As_rf) = (can use for e.g. the value used for
      the "Threshold on the support area" above)
    Flow discharge (Qsf_rf,w_rf) = (Q is the reference discharge for *each*
      rivulet on a hillslope cell; see below; w_rf = 1 is OK)
    Water-surface width (Wsf_rf,b1_rf,b2_rf) = (see below)
    Resistance coefficient (kSsf_rf,y1_rf,y2_rf) = (see below)
    Initial flow discharge (Qsi_rf) = (always 0.0)
  - "channel network parameters":
    Reference drainage area (As_cf) = (can use the area of the entire catchment)
    Initial flow discharge (Qsi_cf) = (alwoys 0.0)
    The Qsf_rf, Wsf_rf, kSsf_rf, Qsf_cf, Wsf_cf, and kSsf_cf parameter values
    should be representative of bankfull discharge flow conditions (return
    period of ~1-2 years).
    For spatially uniform rivulet or channel geometry settings, see the
    Leopold and Maddock scaling equations for the appropriate values to
    assign to parameters b1_rf, b2_rf, y1_rf, y2_rf, b1_cf, b2_cf, y1_cf,
    and y2_cf.
    To simulate "sheet flow" conditions, see Sulis et al, AWR, 2010.

