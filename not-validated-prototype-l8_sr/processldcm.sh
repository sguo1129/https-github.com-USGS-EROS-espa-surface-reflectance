#!/bin/sh
#build the input file for processing ldcm data
rootdir=$1
case=$2
mtlfile=`ls $rootdir$case/*.txt`
echo 0 >ldcm.inp
ls $rootdir$case/*_B?.TIF >>ldcm.inp
ls $rootdir$case/*_B1?.TIF >>ldcm.inp
ls $rootdir$case/*_BQA.TIF >>ldcm.inp
yearday=`echo $case | cut -c 10-16`
echo $yearday
fileanc=`ls LANDSATANC/L8ANC$yearday.hdf_fused`
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
#exit

#LDCMSR-v1.0 <ldcm.inp >$case.log
#LDCMSR-v1.0 <ldcm.inp 
#mv correcteddata.hdf sr$case.hdf


