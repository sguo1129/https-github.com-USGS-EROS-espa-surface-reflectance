Installation Guide 
(tested in Linux bash, August 2012)

1. Install dependent libraries - GCTP, HDF4, HDF-EOS2, TIFF, GeoTIFF

2. Set up environment.  Can look at and modify ledapsSrc/src/env.sh or add
the following to your bash shell.  For C shell, use 'setenv VAR "directory"'.
export GCTP_INC="path_to_GCTP_include_files"
export GCTP_LIB="path_to_GCTP_libraries"
export TIFF_INC="path_to_TIFF_include_files"
export TIFF_LIB="path_to_TIFF_libraries"
export GEOTIFF_INC="path_to_GEOTIFF_include_files"
export GEOTIFF_LIB="path_to_GEOTIFF_libraries"
export HDF_INC="path_to_HDF_include_files"
export HDF_LIB="path_to_HDF_libraries"
export HDFEOS_INC="path_to_HDFEOS_include_files"
export HDFEOS_LIB="path_to_HDFEOS_libraries"
export BIN="path_to_directory_for_LEDAPS_binaries"

3. Install ancillary files
tar -xvzf ledapsAnc.1980-2012.tar.gz

4. Checkout (from Google ledaps project) and install source files
cd ledapsSrc/src
make
make install
make clean

This will create a list of executable files under $BIN:
lndapp lndcal  lndpm  lndsr  lndsrbm.ksh sixsV1.0B
cmrbv1.0  compadjn  comptemp  geo2xy SDSreader3.0 xy2geo
(tested in gcc and gfortran compiler)

* Note that new cloud mask program (lndsrbm) calls command 
"ncdump" which is a standard HDF command. This release includes 
a linux version ncdump command in your path. You need replace 
ncdump program if it runs from different OS systems. If you don't 
want to run optional lndsr-based cloud mask (lndsrbm), you will 
not need the ncdump program.

3. Setup environment
export ANC_PATH="directory_saved_ancillary_files"
(or in c shell use 
setenv ANC_PATH "directory_saved_ancillary_files")

4. Test 
download Landsat files and then run following commands separately
lndpm <Landsat_meta_file>
lndcal <lndcal_input_text>
lndsr <lndsr_input_text>
lndsrbm.ksh <lndsr_input_text>

OR simply run the combined commands to generate all (includes lndsrbm)
do_ledaps.csh <Landsat_meta_file>

Create a file called do_ledaps.csh in your $BIN directory.
+++++++++++++++
#! /bin/csh -f

if $#argv != 1 then
    echo "Usage: do_ledaps.csh <Landsat_MTL_file>"
    exit
else
    set meta_file = $argv[1]
    set meta = `echo $meta_file | sed -e 's/.txt//' -e 's/_MTL//' -e 's/.met//'`
endif

# run LEDAPS modules
$BIN/lndpm $meta_file
$BIN/lndcal lndcal.$meta.txt
$BIN/lndsr lndsr.$meta.txt
$BIN/lndsrbm.ksh lndsr.$meta.txt
+++++++++++++++

5. Check output
lndcal.*.hdf: top-of-atmosphere (TOA) reflectance in HDF format
lndth.*.hdf: thermal brightness temperature in HDF format
lndsr.*.hdf: surface reflectance in HDF format
(each file is also associated with a ENVI header file that contains geographic
information)

REFERENCES
Jeffrey G. Masek, Eric F. Vermote, Nazmi E. Saleous, Robert Wolfe,
Forrest G. Hall, Karl F. Huemmrich, Feng Gao, Jonathan Kutler, and 
Teng-Kui Lim, 2006, A Landsat Surface Reflectance Dataset for North
America, 1990-2000, IEEE Geoscience and Remote Sensing Letters, 
Vol. 3, No. 1, pp. 68-72.

