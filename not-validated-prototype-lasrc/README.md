## LaSRC Version 0.7.0 Release Notes
Release Date: July 2016

### Downloads
LaSRC (Landsat Surface Reflectance Code) source code

    git clone https://github.com/USGS-EROS/espa-surface-reflectance.git

LaSRC auxiliary files

    http://espa.cr.usgs.gov/downloads/auxiliaries/l8sr_auxiliary/l8sr_auxiliary.tar.gz

See git tag [lasrc-version_0.7.0]

### Installation
  * Install dependent libraries - ESPA product formatter (https://github.com/USGS-EROS/espa-product-formatter)
  * Set up environment variables.  Can create an environment shell file or add the following to your bash shell.  For C shell, use 'setenv VAR "directory"'.
```
    export PREFIX="path_to_directory_for_lasrc_build_data"
```

  * Install baseline auxiliary files and set up the environment variables.
```
    tar -xvzf l8sr_auxiliary.tar.gz
    export L8_AUX_DIR="directory_saved_auxiliary_files"
    (or in c shell use 
    setenv L8_AUX_DIR "directory_saved_auxiliary_files")
```

  * Download (from Github USGS-EROS surface-reflectance project) and install source files. The following build will create a list of executable files under $PREFIX/bin (tested in Linux with the gcc compiler). It will also copy various scripts from the scripts directory up the the $PREFIX/bin directory.
```
    cd lasrc\c_version\src
    make
    make install
    make clean

    cd lasrc\landsat_aux\src
    make
    make install
    make clean
```

  * Test - Download Landsat Level 1 files.  Run the do\_lasrc Python script in the PREFIX bin directory to run the applications.  Use do\_lasrc.py --help for the usage information.  This script requires that your LaSRC binaries are in your $PATH or that you have a $BIN environment variable set up to point to the PREFIX bin directory.
```
    convert_lpgs_to_espa --mtl <Landsat_MTL_file> --xml <Landsat_ESPA_XML_file>
    do_lasrc.py --xml <Landsat_ESPA_XML_file>
```

  * Check output
```
    {scene_name}_toa_*: top-of-atmosphere (TOA) reflectance (or brightness temperature for thermal bands) in internal ESPA file format
    {scene_name}_sr_*: surface reflectance in internal ESPA file format
```

### Dependencies
  * ESPA raw binary and ESPA common libraries from ESPA product formatter and associated dependencies
  * XML2 library
  * Auxiliary data products
    1. LAADS Terra and Aqua CMG and CMA data
    2. CMGDEM HDF file
    3. Various input files and model information provided with the LaSRC auxiliary .tar.gz file

### Auxiliary Data Updates
The baseline auxiliary files provided don't include the daily climate data.  In order to generate or update the auxiliary files to the most recent day of year (actually the most current auxiliary files available will be 2-3 days prior to the current day of year do to the latency of the underlying LAADS products) the user will want to run the updatelads.py script available in $PREFIX/bin.  This script can be run with the "--help" argument to print the usage information.  In general the --quarterly argument will reprocess/update all the LAADS data back to 2013.  This is good to do every once in a while to make sure any updates to the LAADS data products are captured.  The --today command-line argument will process the LAADS data for the most recent year.  In general, it is suggested to run the script with --quarterly once a quarter.  Then run the script with --today on a nightly basis.  This should provide an up-to-date version of the auxiliary input data for LaSRC.  The easiest way to accomplish this is to set up a nightly and quarterly cron job.

The updatelads script requires a username/password to access the ladssci.nascom.nasa.gov FTP site.  The user will need to contact USGS EROS Customer Services to obtain a username/password for the LAADS FTP site.  In your email explain that you will be using this ftp access to obtain LAADS data for processing Landsat 8 products using the LaSRC application provided by the USGS EROS.  For questions regarding this information, please contact the Landsat Contact Us page and specify USGS CDR/ECV in the "Regarding" section. https://landsat.usgs.gov/contactus.php

The provided username and password should be used in the --username and --password command-line arguments for the updatelads.py script.  If not specified the source code will try to use the ESPA_LAADS_CONFIG http service to automatically determine the username/password, which is only available to the USGS LSRD systems.

### Data Preprocessing
This version of the LaSRC application requires the input Landsat products to be in the ESPA internal file format.  After compiling the product formatter raw\_binary libraries and tools, the convert\_lpgs\_to\_espa command-line tool can be used to create the ESPA internal file format for input to the LaSRC application.

### Data Postprocessing
After compiling the product-formatter raw\_binary libraries and tools, the convert\_espa\_to\_gtif and convert\_espa\_to\_hdf command-line tools can be used to convert the ESPA internal file format to HDF or GeoTIFF.  Otherwise the data will remain in the ESPA internal file format, which includes each band in the ENVI file format (i.e. raw binary file with associated ENVI header file) and an overall XML metadata file.

### Verification Data

### User Manual

### Product Guide

## Release Notes
  1. Updated the FORTRAN code to be the latest version (v3.0) of software
     received from NASA GSFC for LaSRC.
     a. Band ratios are interpolated at the pixel level versus the CMG level,
        which helps solve the blockiness results from the previous algorithm.
     b. Aerosols are not retrieved over cirrus pixels, however they are
        retrieved for all other non-fill pixels (including water).  The results
        of the aerosol retrieval are tested and flagged if the retrieval does
        not meet residual and NDVI criteria.  These flagged pixels are
        attempted to be interpolated via aerosol interpolation.  Cirrus, cloud,
        and water pixels are not used as part of the interpolation.  The
        aerosol interpolation process has changed and allows the results of the
        interpolation to be at the pixel level versus a block level.  The final
        step is to perform the atmospheric correction based on the calculated or
        interpolated aerosols.  This level of correction is not applied to
        cirrus or cloud pixels.
  2. Updated the C version of the LaSRC code to include the modifications
     delivered as part of v3.0 of the FORTRAN source code.
  3. Merged in changes from version 0.6.1 and 0.6.2 for the updatelads.py
     script.
