#!/bin/ksh
#build the input file for processing ldcm data
rootdir=$1
ver=$2
for case in `ls  $rootdir | awk -F / '{print $NF}'`
do
mtlfile=`ls $rootdir$case/*.txt`
fileinp=ldcm$ver$case.inp
echo 0 >$fileinp
ls $rootdir$case/*_B?.TIF >>$fileinp
ls $rootdir$case/*_B1?.TIF >>$fileinp
ls $rootdir$case/*_BQA.TIF >>$fileinp
yearday=`echo $case | cut -c 10-16`
echo $yearday
year=`echo $yearday | cut -c 1-4`
day=`echo $yearday | cut -c 5-7`
fileanc=`ls /LANDSRx/f008/jroger/MODIS-ANC/L8ANC$yearday.hdf_fused`
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
time=`grep SCENE_CENTER_TIME $mtlfile |  awk '{print $3}' | tr -d '"' | tr -d Z | tr -s : " "`
echo $day $time >>$fileinp
echo sr$ver$case.hdf  >>$fileinp
#echo 2603 3870 2759 3771 >>$fileinp
#./LDCMSR-$ver <$fileinp >log$case$ver
./LDCMSR-$ver <$fileinp 
done

