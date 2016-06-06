#!/bin/ksh
#build the input file for processing ldcm data
#processldcm.ksh $rootdir $scenedir
# Example - processldcm.ksh ./l8sr_test/FORTRAN_version/ LC81440402015045LGN00
#   The last line support subsetting.  Using all zeros will indicate full res.
rootdir=$1
case=$2
mtlfile=`ls $rootdir$case/*.txt`
echo 0 >ldcm.inp
ls $rootdir$case/*_B?.TIF >>ldcm.inp
ls $rootdir$case/*_B1?.TIF >>ldcm.inp
ls $rootdir$case/*_BQA.TIF >>ldcm.inp
year=`echo $case | cut -c 10-13`
yearday=`echo $case | cut -c 10-16`
echo $yearday
fileanc=`ls /usr/local/ledaps/L8ANC/LADS/$year/L8ANC$yearday.hdf_fused`
echo $fileanc >>ldcm.inp
mtlfile=`ls $rootdir$case/*.txt`
rnl=`grep REFLECTIVE_LINES $mtlfile | awk '{print $3}'`
rnc=`grep REFLECTIVE_SAMPLES $mtlfile | awk '{print $3}'`
pnl=`grep  PANCHROMATIC_LINES $mtlfile | awk '{print $3}'`
pnc=`grep PANCHROMATIC_SAMPLES $mtlfile | awk '{print $3}'`
echo $rnl $rnc $pnl $pnc >>ldcm.inp
ts=`grep SUN_ELEVATION $mtlfile | awk '{print 90.-$3}'`
fs=`grep SUN_AZIMUTH $mtlfile | awk '{print $3}'`
echo $ts $fs >>ldcm.inp
utmzone=`grep  "UTM_ZONE" $mtlfile | awk '{print $3}'`
x0=`grep  "CORNER_UL_PROJECTION_X_PRODUCT" $mtlfile | awk '{print $3}'`
y0=`grep  "CORNER_UL_PROJECTION_Y_PRODUCT" $mtlfile | awk '{print $3}'`
echo $utmzone 1 1 $y0 $x0 >>ldcm.inp
echo correcteddata.hdf >>ldcm.inp
echo 0 0 0 0 >>ldcm.inp
#exit

#LDCMSR-v3.0 <ldcm.inp >$case.log
/home/gschmidt/l8sr_version3/buildme/LDCMSR-v3.0 <ldcm.inp 
mv correcteddata.hdf sr$case.hdf


