#!/bin/sh
#build the input file for processing ldcm data
rootdir=$1
case=$2
scene_dir=$rootdir/$case
mtlfile=`ls $scene_dir/*.txt`
inpfile=$case\.input
echo 0 >$inpfile
ls $scene_dir/*_B?.TIF >>$inpfile
ls $scene_dir/*_B1?.TIF >>$inpfile
ls $scene_dir/*_BQA.TIF >>$inpfile
yearday=`echo $case | cut -c 10-16`
echo $yearday
fileanc=`ls LANDSATANC/L8ANC$yearday.hdf_fused`
echo $fileanc >>$inpfile

rnl=`grep REFLECTIVE_LINES $mtlfile | awk '{print $3}'`
rnc=`grep REFLECTIVE_SAMPLES $mtlfile | awk '{print $3}'`
pnl=`grep  PANCHROMATIC_LINES $mtlfile | awk '{print $3}'`
pnc=`grep PANCHROMATIC_SAMPLES $mtlfile | awk '{print $3}'`
echo $rnl $rnc $pnl $pnc >>$inpfile
ts=`grep SUN_ELEVATION $mtlfile | awk '{print 90.-$3}'`
fs=`grep SUN_AZIMUTH $mtlfile | awk '{print $3}'`
echo $ts $fs >>$inpfile
utmzone=`grep  "UTM_ZONE" $mtlfile | awk '{print $3}'`
x0=`grep  "CORNER_UL_PROJECTION_X_PRODUCT" $mtlfile | awk '{print $3}'`
y0=`grep  "CORNER_UL_PROJECTION_Y_PRODUCT" $mtlfile | awk '{print $3}'`
echo $utmzone 1 1 $y0 $x0 >>$inpfile
#exit

#LDCMSR-v1.0 <$inpfile >$case.log
#LDCMSR-v1.0 <$inpfile
#mv correcteddata.hdf sr$case.hdf


