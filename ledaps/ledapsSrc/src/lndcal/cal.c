#include "cal.h"
#include "const.h"
#include "error.h"
#define nint(A)(A<0?(int)(A-0.5):(int)(A+0.5))

/* Functions */
/* !Revision:
 *
 * revision 1.1.0 9/13/2012  Gail Schmidt, USGS
 * - modified cal6 application to flag the saturated thermal pixels to be
 *   consistent with the processing of the reflective bands
 * revision 1.2.0 2/20/2013  Gail Schmidt, USGS
 * - added computeBounds routine for computing the bounding coordinates
 * revision 1.2.1 3/22/2013  Gail Schmidt, USGS
 * - writing UL and LR corners to the output metadata to be able to detect
 *   ascending scenes or scenes where the image is flipped North to South
 * revision 1.2.2 8/1/2013  Gail Schmidt, USGS
 * - validated the TOA reflectance values to make sure they were within
 *   the valid range of values
 * - validated the thermal values to make sure they were within the valid
 *   range of values
 * revision 2.0.0 1/30/2014  Gail Schmidt, USGS
 * - modified the brightness temp values to be written in Kelvin vs. degrees
 *   Celsius
 *
 * NOTES:
 * 1. TOA radiance and reflectance equations are available in
 *    http://landsathandbook.gsfc.nasa.gov/data_prod/prog_sect11_3.html
 * 2. The TOA reflectance gain/bias values from the MTL file (stored in the
 *    XML file after converting from LPGS to ESPA) do not account for the
 *    solar angle.  Thus the gain and bias need to be applied and then we
 *    still need to account for the solar angle.
 */

bool Cal(Lut_t *lut, int iband, Input_t *input, unsigned char *line_in, 
         int16 *line_out, unsigned char *line_out_qa, Cal_stats_t *cal_stats,
         int iy) {
  int is,val;
  float rad_gain, rad_bias;           /* TOA radiance gain/bias */
  float refl_gain = 0.0,
        refl_bias = 0.0;              /* TOA reflectance gain/bias */
  float rad;                          /* TOA radiance value */
  float ref_conv = 0.0;               /* TOA reflectance conversion value */
  float ref;                          /* TOA reflectance value */
  float fval;                         /* temporary float value */
  int nsamp= input->size.s;
  int ifill= (int)lut->in_fill;

  /* Get the TOA radiance gain/bias */
  rad_gain = lut->meta.rad_gain[iband];
  rad_bias = lut->meta.rad_bias[iband];

  /* Get the TOA reflectance gain/bias if they are available, otherwise use
     the TOA reflectance equation from the Landsat handbook. */
  if (input->meta.use_toa_refl_consts) {
    refl_gain = lut->meta.refl_gain[iband];
    refl_bias = lut->meta.refl_bias[iband];

    if ( iy==0 ) {
      printf("*** band=%1d refl gain=%f refl bias=%f cos_sun_zen=%f\n", iband+1,
             refl_gain, refl_bias, lut->cos_sun_zen);
      fflush(stdout);
    }
  }
  else {
    ref_conv = (PI * lut->dsun2) / (lut->esun[iband] * lut->cos_sun_zen);
  
    if ( iy==0 ) {
      printf("*** band=%1d rad gain=%f rad bias=%f dsun2=%f\n"
             "    ref_conv=%f=(PI*%f)/(%f*%f) ***\n", iband+1,
             rad_gain, rad_bias, lut->dsun2, ref_conv, lut->dsun2,
             lut->esun[iband], lut->cos_sun_zen);
      fflush(stdout);
    }
  }

  /* Loop through the samples in the line */
  for (is = 0; is < nsamp; is++) {
    val= getValue((unsigned char *)line_in, is);
    if (val == ifill || line_out_qa[is]==lut->qa_fill ) {
      line_out[is] = lut->out_fill;
      cal_stats->nfill[iband]++;
      continue;
    }

    /* flag saturated pixels, added by Feng (3/23/09) */
    if (val == SATU_VAL[iband]) {
      line_out[is] = lut->out_satu;
      continue;
    }

    cal_stats->nvalid[iband]++;
    fval= (float)val;

    /* If the TOA reflectance gain/bias values are available, then use them.
       Otherwise compute the TOA radiance then reflectance, per the Landsat
       handbook equations. */
    rad = (rad_gain * fval) + rad_bias;
    if (input->meta.use_toa_refl_consts) {
      ref = ((refl_gain * fval) + refl_bias) / lut->cos_sun_zen;
    }
    else {
      ref = rad * ref_conv;
    }

    /* Apply a scaling of 10000 (tied to the lut->scale_factor). Valid ranges
       are set up in lut.c as well. */
    line_out[is] = (int16)(ref * 10000.0 + 0.5);

    /* Cap the output using the min/max values.  Then reset the toa reflectance
       value so that it's correctly reported in the stats and the min/max
       range matches that of the image data. */
    if (line_out[is] < lut->valid_range_ref[0]) {
      line_out[is] = lut->valid_range_ref[0];
      ref = line_out[is] * 0.0001;
    }
    else if (line_out[is] > lut->valid_range_ref[1]) {
      line_out[is] = lut->valid_range_ref[1];
      ref = line_out[is] * 0.0001;
    }

    if (cal_stats->first[iband]) {
      cal_stats->idn_min[iband] = val;
      cal_stats->idn_max[iband] = val;

      cal_stats->rad_min[iband] = rad;
      cal_stats->rad_max[iband] = rad;

      cal_stats->ref_min[iband] = ref;
      cal_stats->ref_max[iband] = ref;

      cal_stats->iref_min[iband] = line_out[is];
      cal_stats->iref_max[iband] = line_out[is];

      cal_stats->first[iband] = false;
    } else {
      if (val < cal_stats->idn_min[iband]) 
        cal_stats->idn_min[iband] = val;
      if (val > cal_stats->idn_max[iband]) 
        cal_stats->idn_max[iband] = val;

      if (rad < cal_stats->rad_min[iband]) cal_stats->rad_min[iband] = rad;
      if (rad > cal_stats->rad_max[iband]) cal_stats->rad_max[iband] = rad;

      if (ref < cal_stats->ref_min[iband]) cal_stats->ref_min[iband] = ref;
      if (ref > cal_stats->ref_max[iband]) cal_stats->ref_max[iband] = ref;

      if (line_out[is] < cal_stats->iref_min[iband]) 
        cal_stats->iref_min[iband] = line_out[is];
      if (line_out[is] > cal_stats->iref_max[iband]) 
        cal_stats->iref_max[iband] = line_out[is];
    }
  }  /* end for is */

  return true;
}

bool Cal6(Lut_t *lut, Input_t *input, unsigned char *line_in, int16 *line_out, 
         unsigned char *line_out_qa, Cal_stats6_t *cal_stats, int iy) {
  int is, val;
  float rad_gain, rad_bias, rad, temp;
  int nsamp= input->size_th.s;
  int ifill= (int)lut->in_fill;

  rad_gain = lut->meta.rad_gain_th;
  rad_bias = lut->meta.rad_bias_th;
  
  if ( iy==0 ) {
    printf("*** band=%1d gain=%f bias=%f ***\n", 6, rad_gain, rad_bias);
  }

  for (is = 0; is < nsamp; is++) {
    val= getValue((unsigned char *)line_in, is);
    if (val == ifill || line_out_qa[is]==lut->qa_fill ) {
      line_out[is] = lut->out_fill;
      cal_stats->nfill++;
      continue;
    }

    /* for saturated pixels */
    if (val >= SATU_VAL6) {
      line_out[is] = lut->out_satu;
      continue;
    }

    cal_stats->nvalid++;
 
    /* compute the brightness temperature in Kelvin and apply scaling of
       10.0 (tied to lut->scale_factor_th). valid ranges are set up in lut.c
       as well. */
    rad = (rad_gain * (float)val) + rad_bias;
    temp = lut->K2 / log(1.0 + (lut->K1/rad));
    line_out[is] = (int16)(temp * 10.0 + 0.5);

    /* Cap the output using the min/max values.  Then reset the temperature
       value so that it's correctly reported in the stats and the min/max
       range matches that of the image data. */
    if (line_out[is] < lut->valid_range_th[0]) {
      line_out[is] = lut->valid_range_th[0];
      temp = line_out[is] * 0.1;
    }
    else if (line_out[is] > lut->valid_range_th[1]) {
      line_out[is] = lut->valid_range_th[1];
      temp = line_out[is] * 0.1;
    }

    if (cal_stats->first) {
      cal_stats->idn_min = val;
      cal_stats->idn_max = val;

      cal_stats->rad_min = rad;
      cal_stats->rad_max = rad;

      cal_stats->temp_min = temp;
      cal_stats->temp_max = temp;

      cal_stats->itemp_min = line_out[is];
      cal_stats->itemp_max = line_out[is];

      cal_stats->first = false;
    } else {
      if (val < (int)cal_stats->idn_min) 
        cal_stats->idn_min = val;
      if (val > cal_stats->idn_max) 
        cal_stats->idn_max = val;

      if (rad < cal_stats->rad_min) cal_stats->rad_min = rad;
      if (rad > cal_stats->rad_max) cal_stats->rad_max = rad;

      if (temp < cal_stats->temp_min) cal_stats->temp_min = temp;
      if (temp > cal_stats->temp_max) cal_stats->temp_max = temp;

      if (line_out[is] < cal_stats->itemp_min) 
        cal_stats->itemp_min = line_out[is];
      if (line_out[is] > cal_stats->itemp_max) 
        cal_stats->itemp_max = line_out[is];
    }
  }  /* end for is */

  return true;
}

/*************************************************************************
 *** this program returns the correct value (as an int)                ***
 *************************************************************************/
int getValue(unsigned char* line_in, int ind)
{
  return (int) line_in[ind];
}
