#include "sr.h"
#include "ar.h"
#include "const.h"
#include "sixs_runs.h"

/* !Revision:
 *
 * revision 1.2.1 3/22/2013  Gail Schmidt, USGS
 * - writing UL and LR corners to the output metadata to be able to detect
 *   ascending scenes or scenes where the image is flipped North to South
 * revision 8/11/2015  Gail Schmidt, USGS
 * - input saturated pixels are flagged as such and output as saturated
 */

extern atmos_t atmos_coef;
int SrInterpAtmCoef(Lut_t *lut, Img_coord_int_t *input_loc, atmos_t *atmos_coef,atmos_t *interpol_atmos_coef); 

bool Sr
(
    Lut_t *lut,
    int nsamp,
    int il,
    int16 **line_in,
    int16 **line_out,
    Sr_stats_t *sr_stats
)
{
    int is;                 /* current sample in the line */
    int ib;                 /* current band for this pixel */
    Img_coord_int_t loc;    /* line/sample location */
    float rho;              /* surface reflectance value */
    float tmpflt;           /* temporary float value for corrections */
    atmos_t interpol_atmos_coef;

    allocate_mem_atmos_coeff(1,&interpol_atmos_coef);
    loc.l = il;

    /* loop through the samples in this line */
    for (is = 0; is < nsamp; is++) {
        loc.s = is;

        /* NAZMI 6/2/04 : correct even cloudy pixels */
        SrInterpAtmCoef(lut, &loc, &atmos_coef, &interpol_atmos_coef);
    
        /* Loop through each band, correcting the pixel.  Fill and saturated
           pixels are skipped and flagged. */
        for (ib = 0; ib < lut->nband; ib++) {
            if (line_in[ib][is] == lut->in_fill) {
                /* fill pixel */
                line_out[ib][is] = lut->output_fill;
                sr_stats->nfill[ib]++;
                continue;
            }
            else if (line_in[ib][is] == lut->in_satu) {
                /* saturated pixel */
                line_out[ib][is] = lut->output_satu;
                sr_stats->nsatu[ib]++;
                continue;
            }
            else {
                rho=(float)line_in[ib][is]/10000.;
                rho=(rho/interpol_atmos_coef.tgOG[ib][0]-
                    interpol_atmos_coef.rho_ra[ib][0]);
                tmpflt=(interpol_atmos_coef.tgH2O[ib][0]*
                    interpol_atmos_coef.td_ra[ib][0]*
                    interpol_atmos_coef.tu_ra[ib][0]);
                rho /= tmpflt;
                rho /= (1.+interpol_atmos_coef.S_ra[ib][0]*rho);
    
                /* Scale the reflectance value and store it as an int16 */
                line_out[ib][is] = (short)(rho*10000.);  /* scale for output */
    
                /* Verify the reflectance value is within the valid range */
                if (line_out[ib][is] < lut->min_valid_sr) {
                    sr_stats->nout_range[ib]++;
                    line_out[ib][is] = lut->min_valid_sr;
                }
                else if (line_out[ib][is] > lut->max_valid_sr) {
                    sr_stats->nout_range[ib]++;
                    line_out[ib][is] = lut->max_valid_sr;
                }
            }
    
            /* Keep track of the min/max value for the stats */
            if (sr_stats->first[ib]) {
                sr_stats->sr_min[ib] = sr_stats->sr_max[ib] = line_out[ib][is];
                sr_stats->first[ib] = false;
            }
            else {
                if (line_out[ib][is] < sr_stats->sr_min[ib])
                    sr_stats->sr_min[ib] = line_out[ib][is];
    
                else if (line_out[ib][is] > sr_stats->sr_max[ib])
                    sr_stats->sr_max[ib] = line_out[ib][is];
            } 
        }  /* end for ib */
    }  /* end for is */

    free_mem_atmos_coeff(&interpol_atmos_coef);
    return true;
}


int SrInterpAtmCoef
(
    Lut_t *lut,                    /* I: lookup table info */
    Img_coord_int_t *input_loc,    /* I: input line/sample location */
    atmos_t *atmos_coef,           /* I: actual atmospheric coefficients */
    atmos_t *interpol_atmos_coef   /* O: interpolated atmospheric
                                         coefficients */
)
/* 
  Point order:

    0 ---- 1    +--> sample
    |      |    |
    |      |    v
    2 ---- 3   line

 */
{
    Img_coord_int_t p[4];
    int i, n,ipt, ib;
    double dl, ds, w;
    double sum[7][13];  /* sum for each coefficient and each band */
    double sum_w;       /* sum of the weights */
    Img_coord_int_t ar_region_half;

    ar_region_half.l = (lut->ar_region_size.l + 1) / 2;
    ar_region_half.s = (lut->ar_region_size.s + 1) / 2;

    p[0].l = (input_loc->l - ar_region_half.l) / lut->ar_region_size.l;

    p[2].l = p[0].l + 1;
    if (p[2].l >= lut->ar_size.l) {
        p[2].l = lut->ar_size.l - 1;
        if (p[0].l > 0)
            p[0].l--;
    }
      
    p[1].l = p[0].l;
    p[3].l = p[2].l;

    p[0].s = (input_loc->s - ar_region_half.s) / lut->ar_region_size.s;
    p[1].s = p[0].s + 1;

    if (p[1].s >= lut->ar_size.s) {
        p[1].s = lut->ar_size.s - 1;
        if (p[0].s > 0)
            p[0].s--;
    }

    p[2].s = p[0].s;
    p[3].s = p[1].s;

    /* Initialize the variables to 0 */
    n = 0;
    sum_w = 0.0;
/* TODO: GAIL -- This ipt < 9 should really be ipt < 13.  Otherwise some of
   the sums are uninitialized */
    for (ib = 0; ib < 7; ib++)
        for (ipt = 0; ipt < 9; ipt++)
            sum[ib][ipt] = 0.;

    /* Loop through the four points to be used in the interpolation */
    for (i = 0; i < 4; i++) {
        /* If the points are valid */
        if (p[i].l != -1 && p[i].s != -1) {
            ipt = p[i].l * lut->ar_size.s + p[i].s;
            if (!(atmos_coef->computed[ipt]))
                continue; 

            dl = (input_loc->l - ar_region_half.l) -
                 (p[i].l * lut->ar_region_size.l);
            dl = fabs(dl) / lut->ar_region_size.l;
            ds = (input_loc->s - ar_region_half.s) -
                 (p[i].s * lut->ar_region_size.s);
            ds = fabs(ds) / lut->ar_region_size.s;
            w = (1.0 - dl) * (1.0 - ds);

            /* Increment the count of valid points and add the current weight
               to the sum of weights */
            n++;
            sum_w += w;

            /* Loop through each band and add in the coefficient * weight */
            for (ib = 0; ib < 6; ib++) {
                sum[ib][0] += (atmos_coef->tgOG[ib][ipt] * w);
                sum[ib][1] += (atmos_coef->tgH2O[ib][ipt] * w);
                sum[ib][2] += (atmos_coef->td_ra[ib][ipt] * w);
                sum[ib][3] += (atmos_coef->tu_ra[ib][ipt] * w);
                sum[ib][4] += (atmos_coef->rho_mol[ib][ipt] * w);
                sum[ib][5] += (atmos_coef->rho_ra[ib][ipt] * w);
                sum[ib][6] += (atmos_coef->td_da[ib][ipt] * w);
                sum[ib][7] += (atmos_coef->tu_da[ib][ipt] * w);
                sum[ib][8] += (atmos_coef->S_ra[ib][ipt] * w);
                sum[ib][9] += (atmos_coef->td_r[ib][ipt] * w);
                sum[ib][10] += (atmos_coef->tu_r[ib][ipt] * w);
                sum[ib][11] += (atmos_coef->S_r[ib][ipt] * w);
                sum[ib][12] += (atmos_coef->rho_r[ib][ipt] * w);
            }
        }  /* end if */
    } /* end for */

    /* If there were valid points */
    if (n > 0) {
        /* Loop through each band and compute the coefficient */
        for (ib = 0; ib < 6; ib++) {
            interpol_atmos_coef->tgOG[ib][0] = sum[ib][0] / sum_w;
            interpol_atmos_coef->tgH2O[ib][0] = sum[ib][1] / sum_w;
            interpol_atmos_coef->td_ra[ib][0] = sum[ib][2] / sum_w;
            interpol_atmos_coef->tu_ra[ib][0] = sum[ib][3] / sum_w;
            interpol_atmos_coef->rho_mol[ib][0] = sum[ib][4] / sum_w;
            interpol_atmos_coef->rho_ra[ib][0] = sum[ib][5] / sum_w;
            interpol_atmos_coef->td_da[ib][0] = sum[ib][6] / sum_w;
            interpol_atmos_coef->tu_da[ib][0] = sum[ib][7] / sum_w;
            interpol_atmos_coef->S_ra[ib][0] = sum[ib][8] / sum_w;
            interpol_atmos_coef->td_r[ib][0] = sum[ib][9] / sum_w;
            interpol_atmos_coef->tu_r[ib][0] = sum[ib][10] / sum_w;
            interpol_atmos_coef->S_r[ib][0] = sum[ib][11] / sum_w;
            interpol_atmos_coef->rho_r[ib][0] = sum[ib][12] / sum_w;
        }
    }

    return 0;
}

