#! /bin/csh -f
##################################################################
# Modified on September 11, 2012 by Gail Schmidt, USGS/EROS
# Remove the call to lndcsm, as the ACCA cloud cover will no longer
# be used in the processing.  The internal cloud mask will be the
# only source of QA information.
#
# Modified on March 7, 2014 by Gail Schmidt, USGS/EROS
# Updated to use the input XML file as part of the switch to the
# ESPA internal file format.
##################################################################

if $#argv != 1 then
    echo "Usage: do_ledaps.csh <Landsat_XML_file>"
    exit
else
    set xml_file = $argv[1]
    set base = `echo $xml_file | sed -e 's/.xml//'`
endif

# run LEDAPS modules
$BIN/lndpm $xml_file
$BIN/lndcal lndcal.$base.txt
$BIN/lndsr lndsr.$base.txt
$BIN/lndsrbm.ksh lndsr.$base.txt
