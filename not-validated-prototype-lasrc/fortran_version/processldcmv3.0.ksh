#!/bin/ksh
#build the input file for processing ldcm data
rootdir=$1
for case in `ls $rootdir`
do
mtlfile=`ls $rootdir$case/*.txt`
fileinp=ldcm$case.inp
echo 0 >$fileinp
ls $rootdir$case/*_B?.TIF >>$fileinp
ls $rootdir$case/*_B1?.TIF >>$fileinp
ls $rootdir$case/*_BQA.TIF >>$fileinp
yearday=`echo $case | cut -c 10-16`
echo $yearday
year=`echo $yearday | cut -c 1-4`
fileanc=`ls /data2/Project-2013/LANDSATANC_COMB/$year/L8ANC$yearday.hdf_fused`
echo $fileanc >>$fileinp
mtlfile=`ls $rootdir$case/*.txt`
rnl=`grep REFLECTIVE_LINES $mtlfile | awk '{print $3}'`
rnc=`grep REFLECTIVE_SAMPLES $mtlfile | awk '{print $3}'`
pnl=`grep  PANCHROMATIC_LINES $mtlfile | awk '{print $3}'`
pnc=`grep PANCHROMATIC_SAMPLES $mtlfile | awk '{print $3}'`
echo $rnl $rnc $pnl $pnc >>$fileinp
ts=`grep SUN_ELEVATION $mtlfile | awk '{print 90.-$3}'`
fs=`grep SUN_AZIMUTH $mtlfile | awk '{print $3}'`
echo $ts $fs >>$fileinp
utmzone=`grep  "UTM_ZONE" $mtlfile | awk '{print $3}'`
x0=`grep  "CORNER_UL_PROJECTION_X_PRODUCT" $mtlfile | awk '{print $3}'`
y0=`grep  "CORNER_UL_PROJECTION_Y_PRODUCT" $mtlfile | awk '{print $3}'`
echo $utmzone 1 1 $y0 $x0 >>$fileinp
echo srv3.0$case.hdf  >>$fileinp
#echo 2900 4200 400 2400 >>$fileinp
echo 1 $rnl 1 $rnc >>$fileinp
#echo srv3.0$case.hdf  >>$fileinp
LDCMSR-v3.0 <$fileinp >log$case &
done

