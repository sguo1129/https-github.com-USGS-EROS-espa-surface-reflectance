#! /bin/csh -f
##################################################################
# Modified on September 11, 2012 by Gail Schmidt, USGS/EROS
# Remove the call to lndcsm, as the ACCA cloud cover will no longer
# be used in the processing.  The internal cloud mask will be the
# only source of QA information.
##################################################################

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
