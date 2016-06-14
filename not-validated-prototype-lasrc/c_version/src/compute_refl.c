#include "lasrc.h"
#include "time.h"

/******************************************************************************
MODULE:  compute_toa_refl

PURPOSE:  Computes the TOA reflectance and at-sensor brightness temps for all
the bands except the pan band.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error computing the reflectance
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
  1. These TOA and BT algorithms match those as published by the USGS Landsat
     team in http://landsat.usgs.gov/Landsat8_Using_Product.php
******************************************************************************/
int compute_toa_refl
(
    Input_t *input,     /* I: input structure for the Landsat product */
    uint16 *qaband,     /* I: QA band for the input image, nlines x nsamps */
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    float xmus,         /* I: cosine of solar zenith angle */
    char *instrument,   /* I: instrument to be processed (OLI, TIRS) */
    int16 **sband       /* O: output TOA reflectance and brightness temp
                              values (scaled) */
)
{
    char errmsg[STR_SIZE];                   /* error message */
    char FUNC_NAME[] = "compute_toa_refl";   /* function name */
    int i;               /* looping variable for pixels */
    int ib;              /* looping variable for input bands */
    int sband_ib;        /* looping variable for output bands */
    int iband;           /* current band */
    float rotoa;         /* top of atmosphere reflectance */
    float tmpf;          /* temporary floating point value */
    float refl_mult;     /* reflectance multiplier for bands 1-9 */
    float refl_add;      /* reflectance additive for bands 1-9 */
    float xcals;         /* radiance multiplier for bands 10 and 11 */
    float xcalo;         /* radiance additive for bands 10 and 11 */
    float k1b10;         /* K1 temperature constant for band 10 */
    float k1b11;         /* K1 temperature constant for band 11 */
    float k2b10;         /* K2 temperature constant for band 10 */
    float k2b11;         /* K2 temperature constant for band 11 */
    uint16 *uband = NULL;  /* array for input image data for a single band,
                              nlines x nsamps */

    /* Allocate space for band data */
    uband = calloc (nlines*nsamps, sizeof (uint16));
    if (uband == NULL)
    {
        sprintf (errmsg, "Error allocating memory for uband");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through all the bands (except the pan band) and compute the TOA
       reflectance and at-sensor brightness temp */
    for (ib = DN_BAND1; ib <= DN_BAND11; ib++)
    {
        /* Don't process the pan band */
        if (ib == DN_BAND8)
            continue;
        printf ("%d ... ", ib+1);

        /* Read the current band and calibrate bands 1-9 (except pan) to
           obtain TOA reflectance. Bands are corrected for the sun angle at
           the center of the scene. */
        if (ib <= DN_BAND9)
        {
            if (ib <= DN_BAND7)
            {
                iband = ib;
                sband_ib = ib;
            }
            else
            {  /* don't count the pan band */
                iband = ib - 1;
                sband_ib = ib - 1;
            }

            if (get_input_refl_lines (input, iband, 0, nlines, uband) !=
                SUCCESS)
            {
                sprintf (errmsg, "Reading band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Get TOA reflectance coefficients for this band from XML file */
            refl_mult = input->meta.gain[iband];
            refl_add = input->meta.bias[iband];

#ifdef _OPENMP
            #pragma omp parallel for private (i, rotoa)
#endif
            for (i = 0; i < nlines*nsamps; i++)
            {
                /* If this pixel is not fill */
                if (qaband[i] != 1)
                {
                    /* Compute the TOA reflectance based on the scene center sun
                       angle.  Scale the value for output. */
                    rotoa = (uband[i] * refl_mult) + refl_add;
                    rotoa = rotoa * MULT_FACTOR / xmus;

                    /* Save the scaled TOA reflectance value, but make
                       sure it falls within the defined valid range. */
                    if (rotoa < MIN_VALID)
                        sband[sband_ib][i] = MIN_VALID;
                    else if (rotoa > MAX_VALID)
                        sband[sband_ib][i] = MAX_VALID;
                    else
                        /* TODO FORTRAN code doesn't round here, but it's
                           probably a good idea */
//                        sband[sband_ib][i] = (int) (roundf (rotoa));
                        sband[sband_ib][i] = (int) rotoa;
                }
                else
                    sband[sband_ib][i] = FILL_VALUE;
            }
        }  /* end if band <= band 9 */

        /* Read the current band and calibrate thermal bands.  Not available
           for OLI-only scenes. */
        else if (ib == DN_BAND10 && strcmp (instrument, "OLI"))
        {
            if (get_input_th_lines (input, 0, 0, nlines, uband) != SUCCESS)
            {
                sprintf (errmsg, "Reading band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Get brightness temp coefficients for this band from XML file */
            xcals = input->meta.gain_th[0];
            xcalo = input->meta.bias_th[0];
            k1b10 = input->meta.k1_const[0];
            k2b10 = input->meta.k2_const[0];

            /* Compute brightness temp for band 10.  Make sure it falls
               within the min/max range for the thermal bands. */
#ifdef _OPENMP
            #pragma omp parallel for private (i, tmpf)
#endif
            for (i = 0; i < nlines*nsamps; i++)
            {
                /* If this pixel is not fill */
                if (qaband[i] != 1)
                {
                    /* Compute the TOA spectral radiance */
                    tmpf = xcals * uband[i] + xcalo;

                    /* Compute the at-satellite brightness temp (K) and
                       scale for output */
                    tmpf = k2b10 / log (k1b10 / tmpf + 1.0);
                    tmpf = tmpf * MULT_FACTOR_TH;  /* scale the value */

                    /* Make sure the brightness temp falls within the specified
                       range */
                    if (tmpf < MIN_VALID_TH)
                        sband[SR_BAND10][i] = MIN_VALID_TH;
                    else if (tmpf > MAX_VALID_TH)
                        sband[SR_BAND10][i] = MAX_VALID_TH;
                    else
                        /* TODO FORTRAN code doesn't round here, but it's
                           probably a good idea */
//                        sband[SR_BAND10][i] = (int) (roundf (tmpf));
                        sband[SR_BAND10][i] = (int) (tmpf);
                }
                else
                    sband[SR_BAND10][i] = FILL_VALUE;
            }
        }  /* end if band 10 */

        else if (ib == DN_BAND11 && strcmp (instrument, "OLI"))
        {
            if (get_input_th_lines (input, 1, 0, nlines, uband) != SUCCESS)
            {
                sprintf (errmsg, "Reading band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Get brightness temp coefficients for this band from XML file */
            xcals = input->meta.gain_th[1];
            xcalo = input->meta.bias_th[1];
            k1b11 = input->meta.k1_const[1];
            k2b11 = input->meta.k2_const[1];

            /* Compute brightness temp for band 11.  Make sure it falls
               within the min/max range for the thermal bands. */
#ifdef _OPENMP
            #pragma omp parallel for private (i, tmpf)
#endif
            for (i = 0; i < nlines*nsamps; i++)
            {
                /* If this pixel is not fill */
                if (qaband[i] != 1)
                {
                    /* Compute the TOA spectral radiance */
                    tmpf = xcals * uband[i] + xcalo;

                    /* Compute the at-satellite brightness temp (K) and
                       scale for output */
                    tmpf = k2b11 / log (k1b11 / tmpf + 1.0);
                    tmpf = tmpf * MULT_FACTOR_TH;  /* scale the value */

                    /* Make sure the brightness temp falls within the specified
                       range */
                    if (tmpf < MIN_VALID_TH)
                        sband[SR_BAND11][i] = MIN_VALID_TH;
                    else if (tmpf > MAX_VALID_TH)
                        sband[SR_BAND11][i] = MAX_VALID_TH;
                    else
                        /* TODO FORTRAN code doesn't round here, but it's
                           probably a good idea */
//                        sband[SR_BAND11][i] = (int) (roundf (tmpf));
                        sband[SR_BAND11][i] = (int) (tmpf);
                }
                else
                    sband[SR_BAND11][i] = FILL_VALUE;
            }
        }  /* end if band 11 */
    }  /* end for ib */
    printf ("\n");

    /* The input data has been read and calibrated. The memory can be freed. */
    free (uband);

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  compute_sr_refl

PURPOSE:  Computes the surfance reflectance for all the reflectance bands.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error computing the reflectance
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
1. Initializes the variables and data arrays from the lookup table and
   auxiliary files.
2. The tauray array was originally read in from a static ASCII file, but it is
   now hardcoded to save time from reading the file each time.  This file was
   generated (like many of the other auxiliary input tables) by running 6S and
   storing the coefficients.
3. Snow pixels are commonly flagged as cloud in this algorithm.
4. Aerosols are not retrieved for cirrus pixels.  They are retrieved for all
   other non-fill pixels, including water.  If the model residual using the
   retrieved aerosols is too high, then that pixel is flagged for potential
   aerosol interpolation.  The pixel is also flagged for potential aerosol
   interpolation if the NDVI test fails.  However, the aerosols and residuals
   are retained for those pixels, even though they were flagged.
   After flagging clouds, shadows, adjacent cloud, a water test is done.  If
   the current pixel is not cirrus, cloud, shadow, or adjacent cloud, then it's
   flagged as water if the NDVI < 0.01.
   Next we loop through all the pixels and do the aerosol interpolation.
   Aerosol interpolation is attempted on pixels that were flagged due to the
   residual or NDVI.  Pixels that are cirrus, cloud, or water are not used for
   the interpolation.
   The last step is to perform the atmospheric correction based on the
   aerosols.  This level of correction is not applied to cirrus or cloud pixels.
******************************************************************************/
int compute_sr_refl
(
    Input_t *input,     /* I: input structure for the Landsat product */
    Espa_internal_meta_t *xml_metadata,
                        /* I: XML metadata structure */
    char *xml_infile,   /* I: input XML filename */
    uint16 *qaband,     /* I: QA band for the input image, nlines x nsamps */
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    float pixsize,      /* I: pixel size for the reflectance bands */
    int16 **sband,      /* I/O: input TOA and output surface reflectance */
    float xts,          /* I: solar zenith angle (deg) */
    float xfs,          /* I: solar azimuth angle (deg) */
    float xmus,         /* I: cosine of solar zenith angle */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm         /* I: auxiliary filename for ozone and water vapor */
)
{
    char errmsg[STR_SIZE];                   /* error message */
    char FUNC_NAME[] = "compute_sr_refl";   /* function name */
    int retval;          /* return status */
    int i, j, k, l;      /* looping variable for pixels */
    int ib;              /* looping variable for input bands */
    int iband;           /* current band */
    int curr_pix;        /* current pixel in 1D arrays of nlines * nsamps */
    int win_pix;         /* current pixel in the line,sample window */
    int inf_pix;         /* location of inferior pixel in line,sample window */
    int sup_pix;         /* location of superior pixel in line,sample window */
    bool isuccess;       /* was the refined corrected value successful? */
    float tmpf;          /* temporary floating point value */
    float rotoa;         /* top of atmosphere reflectance */
    float roslamb;       /* lambertian surface reflectance */
    float ros2b1;        /* temp storage of the roslamb value for band 1 */
    float tgo;           /* other gaseous transmittance (tgog * tgoz) */
    float roatm;         /* atmospheric intrinsic reflectance */
    float ttatmg;        /* total atmospheric transmission */
    float satm;          /* atmosphere sphereical albedo */
    float xrorayp;       /* reflectance of the atmosphere due to molecular
                            (Rayleigh) scattering */
    float next;
    float erelc[NSR_BANDS];    /* band ratio variable for bands 1-7 */
    float troatm[NSR_BANDS];   /* atmospheric reflectance table for bands 1-7 */
    float btgo[NSR_BANDS];     /* other gaseous transmittance for bands 1-7 */
    float broatm[NSR_BANDS];   /* atmospheric reflectance for bands 1-7 */
    float bttatmg[NSR_BANDS];  /* ttatmg for bands 1-7 */
    float bsatm[NSR_BANDS];    /* atmosphere spherical albedo for bands 1-7 */

    int iband1, iband3; /* band indices (zero-based) */
    float raot;         /* AOT reflectance */
    float residual;     /* model residual */
    float rsurf;        /* surface reflectance */
    float corf;         /* aerosol impact (higher values represent high
                           aerosol) */
    long nbclear;       /* count of the clear (non-cloud) pixels */
    long nbval;         /* count of the non-fill pixels */
    double anom;        /* band 3 and 5 combination */
    double mall;        /* average/mean temp of all the pixels */
    double mclear;      /* average/mean temp of the clear pixels */
    float fack, facl;   /* cloud height factor in the k,l dim */
    int cldhmin;        /* minimum bound of the cloud height */
    int cldhmax;        /* maximum bound of the cloud height */
    float cldh;         /* cloud height */
    int icldh;          /* looping variable for cloud height */
    int mband5, mband5k, mband5l;    /* band 6 value and k,l locations */
    float tcloud;       /* temperature of the current pixel */

    float cfac = 6.0;     /* cloud factor */
    float fndvi;          /* NDVI value */
    int inf_index;        /* inferior index for aerosol interpolation */
    int sup_index;        /* superior index for aerosol interpolation */
    int int_start;        /* index for the start of interpolation */
    int int_end;          /* index for the end of interpolation */
    float ros4, ros5;     /* surface reflectance for band 4 and band 5 */
    int tmp_percent;      /* current percentage for printing status */
#ifndef _OPENMP
    int curr_tmp_percent; /* percentage for current line */
#endif

    float lat, lon;       /* pixel lat, long location */
    int lcmg, scmg;       /* line/sample index for the CMG */
    int lcmg1, scmg1;     /* line+1/sample+1 index for the CMG */
    float u, v;           /* line/sample index for the CMG */
    float th1, th2;       /* values for NDWI calculations */
    float xcmg, ycmg;     /* x/y location for CMG */
    float xndwi;          /* calculated NDWI value */
    int uoz11, uoz21, uoz12, uoz22;  /* ozone at line,samp; line, samp+1;
                           line+1, samp; and line+1, samp+1 */
    float pres11, pres12, pres21, pres22;  /* pressure at line,samp;
                           line, samp+1; line+1, samp; and line+1, samp+1 */
    uint8 *cloud = NULL;  /* bit-packed value that represent clouds,
                             nlines x nsamps */
    uint8 *ipflag = NULL; /* QA flag to assist with aerosol interpolation,
                              nlines x nsamps */
    float *twvi = NULL;   /* interpolated water vapor value,
                             nlines x nsamps */
    float *tozi = NULL;   /* interpolated ozone value, nlines x nsamps */
    float *tp = NULL;     /* interpolated pressure value, nlines x nsamps */
    float *tresi = NULL;  /* residuals for each pixel, nlines x nsamps;
                             tresi < 0.0 flags water pixels and pixels with
                             high residuals */
    float *taero = NULL;  /* aerosol values for each pixel, nlines x nsamps */
    int16 *aerob1 = NULL; /* atmospherically corrected band 1 data
                             (TOA refl), nlines x nsamps */
    int16 *aerob2 = NULL; /* atmospherically corrected band 2 data
                             (TOA refl), nlines x nsamps */
    int16 *aerob4 = NULL; /* atmospherically corrected band 4 data
                             (TOA refl), nlines x nsamps */
    int16 *aerob5 = NULL; /* atmospherically corrected band 5 data
                             (TOA refl), nlines x nsamps */
    int16 *aerob7 = NULL; /* atmospherically corrected band 7 data
                             (TOA refl), nlines x nsamps */

    /* Vars for forward/inverse mapping space */
    Geoloc_t *space = NULL;       /* structure for geolocation information */
    Space_def_t space_def;        /* structure to define the space mapping */
    Img_coord_float_t img;        /* coordinate in line/sample space */
    Geo_coord_t geo;              /* coordinate in lat/long space */

    /* Lookup table variables */
    float xtv;           /* observation zenith angle (deg) -- NOTE: set to 0.0
                            and never changes */
    float xmuv;          /* cosine of observation zenith angle */
    float xfi;           /* azimuthal difference between the sun and
                            observation angle (deg) */
    float cosxfi;        /* cosine of azimuthal difference */
    float xtsstep;       /* solar zenith step value */
    float xtsmin;        /* minimum solar zenith value */
    float xtvstep;       /* observation step value */
    float xtvmin;        /* minimum observation value */
    float ****rolutt = NULL;    /* intrinsic reflectance table
                                   [NSR_BANDS][7][22][8000] */
    float ****transt = NULL;    /* transmission table
                                   [NSR_BANDS][7][22][22] */
    float ***sphalbt = NULL;    /* spherical albedo table [NSR_BANDS][7][22] */
    float ***normext = NULL;    /* aerosol extinction coefficient at the
                                   current wavelength (normalized at 550nm)
                                   [NSR_BANDS][7][22] */
    float **tsmax = NULL;       /* maximum scattering angle table [20][22] */
    float **tsmin = NULL;       /* minimum scattering angle table [20][22] */
    float **nbfi = NULL;        /* number of azimuth angles [20][22] */
    float **nbfic = NULL;       /* communitive number of azimuth angles
                                   [20][22] */
    float **ttv = NULL;         /* view angle table [20][22] */
    float tts[22];              /* sun angle table */
    int32 indts[22];

    /* Auxiliary file variables */
    int16 *dem = NULL;        /* CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 *andwi = NULL;      /* avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *sndwi = NULL;      /* standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob1 = NULL;    /* mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob2 = NULL;    /* mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob7 = NULL;    /* mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob1 = NULL;   /* intercept band1 ratio,
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *intratiob2 = NULL;   /* intercept band2 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *intratiob7 = NULL;   /* intercept band7 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob1 = NULL;   /* slope band1 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob2 = NULL;   /* slope band2 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob7 = NULL;   /* slope band7 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    uint16 *wv = NULL;       /* water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 *oz = NULL;        /* ozone values [CMG_NBLAT x CMG_NBLON] */
    float raot550nm;    /* nearest input value of AOT */
    float uoz;          /* total column ozone */
    float uwv;          /* total column water vapor (precipital water vapor) */
    float pres;         /* surface pressure */
    float rb1;          /* band ratio 1 (unscaled) */
    float rb2;          /* band ratio 2 (unscaled) */
    float slpr11, slpr12, slpr21, slpr22;  /* band ratio slope at line,samp;
                           line, samp+1; line+1, samp; and line+1, samp+1 */
    float intr11, intr12, intr21, intr22;  /* band ratio intercept at line,samp;
                           line, samp+1; line+1, samp; and line+1, samp+1 */
    float slprb1, slprb2, slprb7;  /* interpolated band ratio slope values for
                                      band ratios 1, 2, 7 */
    float intrb1, intrb2, intrb7;  /* interpolated band ratio intercept values
                                      for band ratios 1, 2, 7 */
    int ratio_pix11;  /* pixel location for ratio products [lcmg][scmg] */
    int ratio_pix12;  /* pixel location for ratio products [lcmg][scmg+1] */
    int ratio_pix21;  /* pixel location for ratio products [lcmg+1][scmg] */
    int ratio_pix22;  /* pixel location for ratio products [lcmg+1][scmg+1] */
    int cmg_pix11;    /* pixel location for CMG/DEM products [lcmg][scmg] */
    int cmg_pix12;    /* pixel location for CMG/DEM products [lcmg][scmg+1] */
    int cmg_pix21;    /* pixel location for CMG/DEM products [lcmg+1][scmg] */
    int cmg_pix22;    /* pixel location for CMG/DEM products [lcmg+1][scmg+1] */

    /* Output file info */
    Output_t *sr_output = NULL;  /* output structure and metadata for the SR
                                    product */
    Envi_header_t envi_hdr;      /* output ENVI header information */
    char envi_file[STR_SIZE];    /* ENVI filename */
    char *cptr = NULL;           /* pointer to the file extension */

    /* Table constants */
    float aot550nm[22] =  /* AOT look-up table */
        {0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 1.20,
         1.40, 1.60, 1.80, 2.00, 2.30, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00};
    float tpres[7] =      /* surface pressure table */
        {1050.0, 1013.0, 900.0, 800.0, 700.0, 600.0, 500.0};

    /* Atmospheric correction variables */
    /* Look up table for atmospheric and geometric quantities */
    float tauray[NSR_BANDS] =  /* molecular optical thickness coefficients --
        produced by running 6S */
        {0.23638, 0.16933, 0.09070, 0.04827, 0.01563, 0.00129, 0.00037,
         0.07984};
    double oztransa[NSR_BANDS] =   /* ozone transmission coeff */
        {-0.00255649, -0.0177861, -0.0969872, -0.0611428, 0.0001, 0.0001,
          0.0001, -0.0834061};
    double wvtransa[NSR_BANDS] =   /* water vapor transmission coeff */
        {2.29849e-27, 2.29849e-27, 0.00194772, 0.00404159, 0.000729136,
         0.00067324, 0.0177533, 0.00279738};
    double wvtransb[NSR_BANDS] =   /* water vapor transmission coeff */
        {0.999742, 0.999742, 0.775024, 0.774482, 0.893085, 0.939669, 0.65094,
         0.759952};
    double ogtransa1[NSR_BANDS] =  /* other gases transmission coeff */
        {4.91586e-20, 4.91586e-20, 4.91586e-20, 1.04801e-05, 1.35216e-05,
         0.0205425, 0.0256526, 0.000214329};
    double ogtransb0[NSR_BANDS] =  /* other gases transmission coeff */
        {0.000197019, 0.000197019, 0.000197019, 0.640215, -0.195998, 0.326577,
         0.243961, 0.396322};
    double ogtransb1[NSR_BANDS] =  /* other gases transmission coeff */
        {9.57011e-16, 9.57011e-16, 9.57011e-16, -0.348785, 0.275239, 0.0117192,
         0.0616101, 0.04728};

    time_t rawtime;
    struct tm *timeinfo;

    /* Allocate memory for the many arrays needed to do the surface reflectance
       computations */
    retval = memory_allocation_sr (nlines, nsamps, &aerob1, &aerob2, &aerob4,
        &aerob5, &aerob7, &cloud, &ipflag, &twvi, &tozi, &tp, &tresi, &taero,
        &dem, &andwi, &sndwi, &ratiob1, &ratiob2, &ratiob7, &intratiob1,
        &intratiob2, &intratiob7, &slpratiob1, &slpratiob2, &slpratiob7, &wv,
        &oz, &rolutt, &transt, &sphalbt, &normext, &tsmax, &tsmin, &nbfic,
        &nbfi, &ttv);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error allocating memory for the data arrays needed "
            "for surface reflectance calculations.");
        error_handler (false, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Initialize the geolocation space applications */
    if (!get_geoloc_info (xml_metadata, &space_def))
    {
        sprintf (errmsg, "Getting the space definition from the XML file");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    space = setup_mapping (&space_def);
    if (space == NULL)
    {
        sprintf (errmsg, "Setting up the geolocation mapping");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Initialize the look up tables and atmospheric correction variables */
    retval = init_sr_refl (nlines, nsamps, input, space, anglehdf, intrefnm,
        transmnm, spheranm, cmgdemnm, rationm, auxnm, &xtv, &xmuv, &xfi,
        &cosxfi, &raot550nm, &pres, &uoz, &uwv, &xtsstep, &xtsmin, &xtvstep,
        &xtvmin, tsmax, tsmin, tts, ttv, indts, rolutt, transt, sphalbt,
        normext, nbfic, nbfi, dem, andwi, sndwi, ratiob1, ratiob2, ratiob7,
        intratiob1, intratiob2, intratiob7, slpratiob1, slpratiob2, slpratiob7,
        wv, oz);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error initializing the lookup tables and "
            "atmospheric correction variables.");
        error_handler (false, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through all the reflectance bands and perform atmospheric
       corrections based on climatology */
    printf ("Performing atmospheric corrections for each reflectance "
        "band ...");
    for (ib = 0; ib <= SR_BAND7; ib++)
    {
        printf (" %d ...", ib+1);

        /* Get the parameters for the atmospheric correction */
        /* rotoa is not defined for this call, which is ok, but the
           roslamb value is not valid upon output. Just set it to 0.0 to
           be consistent. */
        rotoa = 0.0;
        retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
            raot550nm, ib, pres, tpres, aot550nm, rolutt, transt, xtsstep,
            xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic,
            nbfi, tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
            ogtransb1, wvtransa, wvtransb, oztransa, rotoa, &roslamb,
            &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Performing lambertian atmospheric correction "
                "type 2.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Save these band-related parameters for later */
        btgo[ib] = tgo;
        broatm[ib] = roatm;
        bttatmg[ib] = ttatmg;
        bsatm[ib] = satm;

        /* Perform atmospheric corrections for bands 1-7 */
#ifdef _OPENMP
        #pragma omp parallel for private (i, j, curr_pix, rotoa, roslamb)
#endif
        for (i = 0; i < nlines; i++)
        {
            curr_pix = i * nsamps;
            for (j = 0; j < nsamps; j++, curr_pix++)
            {
                /* If this pixel is not fill.  Otherwise fill pixels have
                   already been marked in the TOA calculations. */
                if (qaband[curr_pix] != 1)
                {
                    /* Store the TOA scaled TOA reflectance values for later
                       use before completing atmospheric corrections */
                    if (ib == DN_BAND1)
                        aerob1[curr_pix] = sband[ib][curr_pix];
                    else if (ib == DN_BAND2)
                        aerob2[curr_pix] = sband[ib][curr_pix];
                    else if (ib == DN_BAND4)
                        aerob4[curr_pix] = sband[ib][curr_pix];
                    else if (ib == DN_BAND5)
                        aerob5[curr_pix] = sband[ib][curr_pix];
                    else if (ib == DN_BAND7)
                        aerob7[curr_pix] = sband[ib][curr_pix];
    
                    /* Apply the atmospheric corrections (ignoring the Rayleigh
                       scattering component and water vapor), and store the
                       scaled value for further corrections.  (NOTE: the full
                       computations are in atmcorlamb2) */
                    rotoa = sband[ib][curr_pix] * SCALE_FACTOR;
                    roslamb = rotoa / tgo;
                    roslamb = roslamb - roatm;
                    roslamb = roslamb / ttatmg;
                    roslamb = roslamb / (1.0 + satm * roslamb);
                    sband[ib][curr_pix] = (int) (roslamb * MULT_FACTOR);
                }
            }  /* end for j */
        }  /* end for i */
    }  /* for ib */
    printf ("\n");

    /* Initialize the band ratios */
    for (ib = 0; ib < NSR_BANDS; ib++)
    {
        erelc[ib] = -1.0;
        troatm[ib] = 0.0;
    }

    /* Interpolate the auxiliary data for each pixel location */
    printf ("Interpolating the auxiliary data ...\n");
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    printf ("DEBUG: Start of aux interpolation: %s", asctime (timeinfo));
    tmp_percent = 0;
#ifdef _OPENMP
    #pragma omp parallel for private (i, j, curr_pix, img, geo, lat, lon, xcmg, ycmg, lcmg, scmg, lcmg1, scmg1, u, v, cmg_pix11, cmg_pix12, cmg_pix21, cmg_pix22, ratio_pix11, ratio_pix12, ratio_pix21, ratio_pix22, uoz11, uoz12, uoz21, uoz22, pres11, pres12, pres21, pres22, rb1, rb2, slpr11, slpr12, slpr21, slpr22, intr11, intr12, intr21, intr22, slprb1, slprb2, slprb7, intrb1, intrb2, intrb7, xndwi, th1, th2, fndvi, iband, iband1, iband3, retval, corf, raot, residual, next, rotoa, raot550nm, roslamb, tgo, roatm, ttatmg, satm, xrorayp, ros5, ros4) firstprivate(erelc, troatm)
#endif
    for (i = 0; i < nlines; i++)
    {
#ifndef _OPENMP
        /* update status, but not if multi-threaded */
        curr_tmp_percent = 100 * i / nlines;
        if (curr_tmp_percent > tmp_percent)
        {
            tmp_percent = curr_tmp_percent;
            if (tmp_percent % 10 == 0)
            {
                printf ("%d%% ", tmp_percent);
                fflush (stdout);
            }
        }
#endif

        curr_pix = i * nsamps;
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            /* If this pixel is fill, then don't process */
            if (qaband[curr_pix] == 1)
            {
                ipflag[curr_pix] = IPFLAG_FILL;
                continue;
            }

            /* Get the lat/long for the current pixel, for the center of
               the pixel */
            /* TODO the line/sample conversion should use the center of the
               pixel, however that's not being done in the FORTRAN code */
            img.l = i - 0.5;
            img.s = j + 0.5;
            img.is_fill = false;
            if (!from_space (space, &img, &geo))
            {
                sprintf (errmsg, "Mapping line/sample (%d, %d) to "
                    "geolocation coords", i, j);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
            lat = geo.lat * RAD2DEG;
            lon = geo.lon * RAD2DEG;

            /* Use that lat/long to determine the line/sample in the
               CMG-related lookup tables, using the center of the UL
               pixel. Note, we are basically making sure the line/sample
               combination falls within -90, 90 and -180, 180 global climate
               data boundaries.  However, the source code below uses lcmg+1
               and scmg+1, which for some scenes may wrap around the dateline
               or the poles.  Thus we need to wrap the CMG data around to the
               beginning of the array. */
            /* Each CMG pixel is 0.05 x 0.05 degrees.  Use the center of the
               pixel for each calculation. */
            /* TODO the line/sample calculation from the x/ycmg values should
               be rounded.  The FORTRAN code does this, then overwrites it.
               To be consistent, we will leave it out for now. */
            ycmg = (89.975 - lat) * 20.0;   /* vs / 0.05 */
            xcmg = (179.975 + lon) * 20.0;  /* vs / 0.05 */
            lcmg = (int) roundf (ycmg);
            scmg = (int) roundf (xcmg);
            if ((lcmg < 0 || lcmg >= CMG_NBLAT) ||
                (scmg < 0 || scmg >= CMG_NBLON))
            {
                sprintf (errmsg, "Invalid line/sample combination for the "
                    "CMG-related lookup tables - line %d, sample %d "
                    "(0-based). CMG-based tables are %d lines x %d "
                    "samples.", lcmg, scmg, CMG_NBLAT, CMG_NBLON);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            /* If the current CMG pixel is at the edge of the CMG array,
               then allow the next pixel for interpolation to wrap around
               the array */
            if (scmg >= CMG_NBLON-1)  /* 180 degrees so wrap around */
                scmg1 = 0;
            else
                scmg1 = scmg + 1;

            if (lcmg >= CMG_NBLAT-1)  /* -90 degrees so wrap around */
                lcmg1 = 0;
            else
                lcmg1 = lcmg + 1;

            /* Determine the fractional difference between the integer location
               and floating point pixel location */
            u = (ycmg - lcmg);
            v = (xcmg - scmg);
            cmg_pix11 = lcmg * CMG_NBLON + scmg;
            cmg_pix12 = lcmg * CMG_NBLON + scmg1;
            cmg_pix21 = lcmg1 * CMG_NBLON + scmg;
            cmg_pix22 = lcmg1 * CMG_NBLON + scmg1;

            /* Interpolate water vapor.  If the water vapor value is fill (=0),
               then use it as-is. */
            twvi[curr_pix] = wv[cmg_pix11] * (1.0 - u) * (1.0 - v) +
                             wv[cmg_pix12] * (1.0 - u) * v +
                             wv[cmg_pix21] * u * (1.0 - v) +
                             wv[cmg_pix22] * u * v;
            twvi[curr_pix] = twvi[curr_pix] * 0.01;   /* vs / 100 */

            /* Interpolate ozone.  If the ozone value is fill (=0), then use a
               default value of 120. */
            uoz11 = oz[cmg_pix11];
            if (uoz11 == 0)
                uoz11 = 120;

            uoz12 = oz[cmg_pix12];
            if (uoz12 == 0)
                uoz12 = 120;

            uoz21 = oz[cmg_pix21];
            if (uoz21 == 0)
                uoz21 = 120;

            uoz22 = oz[cmg_pix22];
            if (uoz22 == 0)
                uoz22 = 120;

            tozi[curr_pix] = uoz11 * (1.0 - u) * (1.0 - v) +
                             uoz12 * (1.0 - u) * v +
                             uoz21 * u * (1.0 - v) +
                             uoz22 * u * v;
            tozi[curr_pix] = tozi[curr_pix] * 0.0025;   /* vs / 400 */

            /* Get the surface pressure from the global DEM.  Set to 1013.0
               (sea level) if the DEM is fill (= -9999), which is likely
               ocean. Also flag the deep water pixels.  The dimensions on the
               DEM array is the same as that of the CMG arrays. Use the current
               pixel locations already calculated. */
            if (dem[cmg_pix11] != -9999)
                pres11 = 1013.0 * exp (-dem[cmg_pix11] * ONE_DIV_8500);
            else
            {
                pres11 = 1013.0;
                cloud[curr_pix] = -128; /* set the water bit in the cloud QA */
            }

            if (dem[cmg_pix12] != -9999)
                pres12 = 1013.0 * exp (-dem[cmg_pix12] * ONE_DIV_8500);
            else
                pres12 = 1013.0;

            if (dem[cmg_pix21] != -9999)
                pres21 = 1013.0 * exp (-dem[cmg_pix21] * ONE_DIV_8500);
            else
                pres21 = 1013.0;

            if (dem[cmg_pix22] != -9999)
                pres22 = 1013.0 * exp (-dem[cmg_pix22] * ONE_DIV_8500);
            else
                pres22 = 1013.0;

            tp[curr_pix] = pres11 * (1.0 - u) * (1.0 - v) +
                           pres12 * (1.0 - u) * v +
                           pres21 * u * (1.0 - v) +
                           pres22 * u * v;

            /* Inverting aerosols -- not retrieved for cirrus pixels */
            if (sband[SR_BAND9][curr_pix] >
                (100.0 / (tp[curr_pix] * ONE_DIV_1013)))
            {  /* Set cirrus bit */
                cloud[curr_pix]++;
            }
            else
            {
                /* Determine the band ratios */
                ratio_pix11 = lcmg * RATIO_NBLON + scmg;
                ratio_pix12 = ratio_pix11 + 1;
                ratio_pix21 = lcmg1 * RATIO_NBLON + scmg;
                ratio_pix22 = ratio_pix21 + 1;

                rb1 = ratiob1[ratio_pix11] * 0.001;  /* vs. / 1000. */
                rb2 = ratiob2[ratio_pix11] * 0.001;  /* vs. / 1000. */
                if ((fabs((rb1 - 0.454878*rb2*rb2 - 0.459559*rb2)) > 0.15) ||
                    (ratiob7[ratio_pix11] < 1000))
                {
                    slpratiob1[ratio_pix11] = 0;
                    slpratiob2[ratio_pix11] = 0;
                    slpratiob7[ratio_pix11] = 0;
                    intratiob1[ratio_pix11] = 327;
                    intratiob2[ratio_pix11] = 482;
                    intratiob7[ratio_pix11] = 2000;
                }
                else if (sndwi[ratio_pix11] < 200)
                {
                    slpratiob1[ratio_pix11] = 0;
                    slpratiob2[ratio_pix11] = 0;
                    slpratiob7[ratio_pix11] = 0;
                    intratiob1[ratio_pix11] = ratiob1[ratio_pix11];
                    intratiob2[ratio_pix11] = ratiob2[ratio_pix11];
                    intratiob7[ratio_pix11] = ratiob7[ratio_pix11];
                }

                rb1 = ratiob1[ratio_pix12] * 0.001;  /* vs. / 1000. */
                rb2 = ratiob2[ratio_pix12] * 0.001;  /* vs. / 1000. */
                if ((fabs((rb1 - 0.454878*rb2*rb2 - 0.459559*rb2)) > 0.15) ||
                    (ratiob7[ratio_pix12] < 1000))
                {
                    slpratiob1[ratio_pix12] = 0;
                    slpratiob2[ratio_pix12] = 0;
                    slpratiob7[ratio_pix12] = 0;
                    intratiob1[ratio_pix12] = 327;
                    intratiob2[ratio_pix12] = 482;
                    intratiob7[ratio_pix12] = 2000;
                }
                else if (sndwi[ratio_pix12] < 200)
                {
                    slpratiob1[ratio_pix12] = 0;
                    slpratiob2[ratio_pix12] = 0;
                    slpratiob7[ratio_pix12] = 0;
                    intratiob1[ratio_pix12] = ratiob1[ratio_pix12];
                    intratiob2[ratio_pix12] = ratiob2[ratio_pix12];
                    intratiob7[ratio_pix12] = ratiob7[ratio_pix12];
                }

                rb1 = ratiob1[ratio_pix21] * 0.001;  /* vs. / 1000. */
                rb2 = ratiob2[ratio_pix21] * 0.001;  /* vs. / 1000. */
                if ((fabs((rb1 - 0.454878*rb2*rb2 - 0.459559*rb2)) > 0.15) ||
                    (ratiob7[ratio_pix21] < 1000))
                {
                    slpratiob1[ratio_pix21] = 0;
                    slpratiob2[ratio_pix21] = 0;
                    slpratiob7[ratio_pix21] = 0;
                    intratiob1[ratio_pix21] = 327;
                    intratiob2[ratio_pix21] = 482;
                    intratiob7[ratio_pix21] = 2000;
                }
                else if (sndwi[ratio_pix21] < 200)
                {
                    slpratiob1[ratio_pix21] = 0;
                    slpratiob2[ratio_pix21] = 0;
                    slpratiob7[ratio_pix21] = 0;
                    intratiob1[ratio_pix21] = ratiob1[ratio_pix21];
                    intratiob2[ratio_pix21] = ratiob2[ratio_pix21];
                    intratiob7[ratio_pix21] = ratiob7[ratio_pix21];
                }

                rb1 = ratiob1[ratio_pix22] * 0.001;  /* vs. / 1000. */
                rb2 = ratiob2[ratio_pix22] * 0.001;  /* vs. / 1000. */
                if ((fabs((rb1 - 0.454878*rb2*rb2 - 0.459559*rb2)) > 0.15) ||
                    (ratiob7[ratio_pix22] < 1000))
                {
                    slpratiob1[ratio_pix22] = 0;
                    slpratiob2[ratio_pix22] = 0;
                    slpratiob7[ratio_pix22] = 0;
                    intratiob1[ratio_pix22] = 327;
                    intratiob2[ratio_pix22] = 482;
                    intratiob7[ratio_pix22] = 2000;
                }
                else if (sndwi[ratio_pix22] < 200)
                {
                    slpratiob1[ratio_pix22] = 0;
                    slpratiob2[ratio_pix22] = 0;
                    slpratiob7[ratio_pix22] = 0;
                    intratiob1[ratio_pix22] = ratiob1[ratio_pix22];
                    intratiob2[ratio_pix22] = ratiob2[ratio_pix22];
                    intratiob7[ratio_pix22] = ratiob7[ratio_pix22];
                }

                slpr11 = slpratiob1[ratio_pix11] * 0.001;  /* vs / 1000 */
                intr11 = intratiob1[ratio_pix11] * 0.001;  /* vs / 1000 */
                slpr12 = slpratiob1[ratio_pix12] * 0.001;  /* vs / 1000 */
                intr12 = intratiob1[ratio_pix12] * 0.001;  /* vs / 1000 */
                slpr21 = slpratiob1[ratio_pix21] * 0.001;  /* vs / 1000 */
                intr21 = intratiob1[ratio_pix21] * 0.001;  /* vs / 1000 */
                slpr22 = slpratiob1[ratio_pix22] * 0.001;  /* vs / 1000 */
                intr22 = intratiob1[ratio_pix22] * 0.001;  /* vs / 1000 */
                slprb1 = slpr11 * (1.0 - u) * (1.0 - v) +
                         slpr12 * (1.0 - u) * v +
                         slpr21 * u * (1.0 - v) +
                         slpr22 * u * v;
                intrb1 = intr11 * (1.0 - u) * (1.0 - v) +
                         intr12 * (1.0 - u) * v +
                         intr21 * u * (1.0 - v) +
                         intr22 * u * v;

                slpr11 = slpratiob2[ratio_pix11] * 0.001;  /* vs / 1000 */
                intr11 = intratiob2[ratio_pix11] * 0.001;  /* vs / 1000 */
                slpr12 = slpratiob2[ratio_pix12] * 0.001;  /* vs / 1000 */
                intr12 = intratiob2[ratio_pix12] * 0.001;  /* vs / 1000 */
                slpr21 = slpratiob2[ratio_pix21] * 0.001;  /* vs / 1000 */
                intr21 = intratiob2[ratio_pix21] * 0.001;  /* vs / 1000 */
                slpr22 = slpratiob2[ratio_pix22] * 0.001;  /* vs / 1000 */
                intr22 = intratiob2[ratio_pix22] * 0.001;  /* vs / 1000 */
                slprb2 = slpr11 * (1.0 - u) * (1.0 - v) +
                         slpr12 * (1.0 - u) * v +
                         slpr21 * u * (1.0 - v) +
                         slpr22 * u * v;
                intrb2 = intr11 * (1.0 - u) * (1.0 - v) +
                         intr12 * (1.0 - u) * v +
                         intr21 * u * (1.0 - v) +
                         intr22 * u * v;

                slpr11 = slpratiob7[ratio_pix11] * 0.001;  /* vs / 1000 */
                intr11 = intratiob7[ratio_pix11] * 0.001;  /* vs / 1000 */
                slpr12 = slpratiob7[ratio_pix12] * 0.001;  /* vs / 1000 */
                intr12 = intratiob7[ratio_pix12] * 0.001;  /* vs / 1000 */
                slpr21 = slpratiob7[ratio_pix21] * 0.001;  /* vs / 1000 */
                intr21 = intratiob7[ratio_pix21] * 0.001;  /* vs / 1000 */
                slpr22 = slpratiob7[ratio_pix22] * 0.001;  /* vs / 1000 */
                intr22 = intratiob7[ratio_pix22] * 0.001;  /* vs / 1000 */
                slprb7 = slpr11 * (1.0 - u) * (1.0 - v) +
                         slpr12 * (1.0 - u) * v +
                         slpr21 * u * (1.0 - v) +
                         slpr22 * u * v;
                intrb7 = intr11 * (1.0 - u) * (1.0 - v) +
                         intr12 * (1.0 - u) * v +
                         intr21 * u * (1.0 - v) +
                         intr22 * u * v;

                /* Use a version of NDWI to calculate the band ratio */
                xndwi = ((double) sband[SR_BAND5][curr_pix] -
                         (double) (sband[SR_BAND7][curr_pix] * 0.5)) /
                        ((double) sband[SR_BAND5][curr_pix] +
                         (double) (sband[SR_BAND7][curr_pix] * 0.5));

                th1 = (andwi[ratio_pix11] + 2.0 * sndwi[ratio_pix11]) * 0.001;
                th2 = (andwi[ratio_pix11] - 2.0 * sndwi[ratio_pix11]) * 0.001;
                if (xndwi > th1)
                    xndwi = th1;
                if (xndwi < th2)
                    xndwi = th2;

                erelc[DN_BAND1] = (xndwi * slprb1 + intrb1);
                erelc[DN_BAND2] = (xndwi * slprb2 + intrb2);
                erelc[DN_BAND4] = 1.0;
                erelc[DN_BAND7] = (xndwi * slprb7 + intrb7);

                /* Retrieve the TOA reflectance values for the current pixel */
                troatm[DN_BAND1] = aerob1[curr_pix] * SCALE_FACTOR;
                troatm[DN_BAND2] = aerob2[curr_pix] * SCALE_FACTOR;
                troatm[DN_BAND4] = aerob4[curr_pix] * SCALE_FACTOR;
                troatm[DN_BAND7] = aerob7[curr_pix] * SCALE_FACTOR;

                /* Retrieve the aerosol information */
                iband1 = DN_BAND4;
                iband3 = DN_BAND1;
                retval = subaeroret (iband1, iband3, xts, xtv, xmus, xmuv,
                    xfi, cosxfi, pres, uoz, uwv, erelc, troatm, tpres,
                    aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
                    xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi,
                    tts, indts, ttv, tauray, ogtransa1, ogtransb0,
                    ogtransb1, wvtransa, wvtransb, oztransa, &raot,
                    &residual, &next);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing atmospheric correction.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
                corf = raot / xmus;

                /* Check the model residual.  Corf represents aerosol impact.
                   Test the quality of the aerosol inversion. */
                if (residual < (0.010 + 0.005 * corf))
                {
                    /* Test if band 5 makes sense */
                    iband = DN_BAND5;
                    rotoa = aerob5[curr_pix] * SCALE_FACTOR;
                    raot550nm = raot;
                    retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                        raot550nm, iband, pres, tpres, aot550nm, rolutt,
                        transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                        normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
                        ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                        ogtransb1, wvtransa, wvtransb, oztransa, rotoa,
                        &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp,
                        &next);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing lambertian "
                            "atmospheric correction type 2.");
                        error_handler (true, FUNC_NAME, errmsg);
                        exit (ERROR);
                    }
                    ros5 = roslamb;

                    /* Test if band 4 makes sense */
                    iband = DN_BAND4;
                    rotoa = aerob4[curr_pix] * SCALE_FACTOR;
                    raot550nm = raot;
                    retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                        raot550nm, iband, pres, tpres, aot550nm, rolutt,
                        transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                        normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
                        ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                        ogtransb1, wvtransa, wvtransb, oztransa, rotoa,
                        &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp,
                        &next);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing lambertian "
                            "atmospheric correction type 2.");
                        error_handler (true, FUNC_NAME, errmsg);
                        exit (ERROR);
                    }
                    ros4 = roslamb;

                    /* Use the NDVI to validate the reflectance values */
                    if ((ros5 > 0.1) && ((ros5 - ros4) / (ros5 + ros4) > 0))
                    {
                        taero[curr_pix] = raot;
                        tresi[curr_pix] = residual;
                    }
                    else
                    {
                        taero[curr_pix] = raot;
                        tresi[curr_pix] = residual;
                        ipflag[curr_pix] = IPFLAG_NDVI_FAIL;
                    }
                }
                else
                {
                    taero[curr_pix] = raot;
                    tresi[curr_pix] = residual;
                    ipflag[curr_pix] = IPFLAG_RESIDUAL_FAIL;
                }
            }  /* end if cirrus */
        }  /* end for j */
    }  /* end for i */

#ifndef _OPENMP
    /* update status */
    printf ("100%%\n");
    fflush (stdout);
#endif

    /* Done with the aerob* arrays */
    free (aerob1);  aerob1 = NULL;
    free (aerob2);  aerob2 = NULL;
    free (aerob4);  aerob4 = NULL;
    free (aerob5);  aerob5 = NULL;
    free (aerob7);  aerob7 = NULL;

    /* Done with the ratiob* arrays */
    free (andwi);  andwi = NULL;
    free (sndwi);  sndwi = NULL;
    free (ratiob1);  ratiob1 = NULL;
    free (ratiob2);  ratiob2 = NULL;
    free (ratiob7);  ratiob7 = NULL;
    free (intratiob1);  intratiob1 = NULL;
    free (intratiob2);  intratiob2 = NULL;
    free (intratiob7);  intratiob7 = NULL;
    free (slpratiob1);  slpratiob1 = NULL;
    free (slpratiob2);  slpratiob2 = NULL;
    free (slpratiob7);  slpratiob7 = NULL;

    /* Done with the DEM, water vapor, and ozone arrays */
    free (dem);  dem = NULL;
    free (wv);  wv = NULL;
    free (oz);  oz = NULL;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    printf ("DEBUG: End of aux interpolation: %s", asctime (timeinfo));

    /* Refine the cloud mask */
    /* Compute the average temperature of the clear, non-water, non-filled
       pixels */
    printf ("Refining the cloud mask ...\n");
    nbval = 0;
    nbclear = 0;
    mclear = 0.0;
    mall = 0.0;
    for (i = 0; i < nlines*nsamps; i++)
    {
        /* If this pixel is fill, then don't process */
        if (qaband[i] != 1)
        {
            /* Keep track of the number of total (non-fill) pixels in addition
               to the sum of the unscaled thermal values */
            nbval++;
            mall += sband[SR_BAND10][i] * SCALE_FACTOR_TH;

            /* Check for clear pixels */
            if ((!btest (cloud[i], CIR_QA)) &&
                (sband[SR_BAND5][i] > 300))
            {
                /* Check to see if this is a clear pixel */
                anom = sband[SR_BAND2][i] - (sband[SR_BAND4][i] * 0.5);
                if (anom < 300)
                {
                    /* Keep track of the number of clear pixels in addition to
                       the sum of the unscaled thermal values */
                    nbclear++;
                    mclear += sband[SR_BAND10][i] * SCALE_FACTOR_TH;
                }
            }
        }
    }  /* end for i */

    /* Compute the average/mean temperature of the clear pixels, otherwise set
       to 275 Kelvin */
    if (nbclear > 0)
        mclear = mclear / nbclear;
    else
        mclear = 275.0;

    /* Compute the average/mean temperature of the clear pixels */
    if (nbval > 0)
        mall = mall / nbval;

    printf ("Average clear temperature %%clear %f Kelvin %f%%\n", mclear,
        nbclear * 100.0 / (nlines * nsamps));
    printf ("Average temperature %f Kelvin %ld total pixels\n", mall, nbval);

    /* Determine the cloud mask */
    for (i = 0; i < nlines*nsamps; i++)
    {
        /* Test all bad retrievals (except for fill value) */
        /* TODO -- FORTRAN code does not skip fill here */
        if (ipflag[i] > IPFLAG_CLEAR && ipflag[i] != IPFLAG_FILL)
        {
            if (((sband[SR_BAND2][i] - sband[SR_BAND4][i] * 0.5) > 500) &&
                ((sband[SR_BAND10][i] * SCALE_FACTOR_TH) < (mclear - 2.0)))
            {  /* Snow or cloud */
                cloud[i] += 2;
            }
        }
    }

    /* Set up the adjacent to something bad (snow or cloud) bit */
    printf ("Setting up the adjacent to something bit ...\n");
    for (i = 0; i < nlines; i++)
    {
        curr_pix = i * nsamps;
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            /* If this pixel is cloud or cirrus, then look at the 11x11
               surrounding window as set adjacent to cloud pixels */
            if (btest (cloud[curr_pix], CLD_QA) ||
                btest (cloud[curr_pix], CIR_QA))
            {
                /* Check the 11x11 window around the current pixel */
                for (k = i-5; k <= i+5; k++)
                {
                    /* Make sure the line is valid */
                    if (k < 0 || k >= nlines)
                        continue;

                    win_pix = k * nsamps + j-5;
                    for (l = j-5; l <= j+5; l++, win_pix++)
                    {
                        /* Make sure the sample is valid */
                        if (l < 0 || l >= nsamps)
                            continue;

                        /* If it's not already set as cloud or cirrus and it's
                           not fill, then set as adjacent cloud */
                        /* TODO FORTRAN code does not skip fill pixels */
                        if (!btest (cloud[win_pix], CLD_QA) &&
                            !btest (cloud[win_pix], CIR_QA) &&
                            !btest (cloud[win_pix], CLDA_QA) &&
                            qaband[win_pix] != 1)
                        {  /* Set the adjacent cloud bit */
                            cloud[win_pix] += 4;
                        }
                    }  /* for l */
                }  /* for k */
            }  /* if btest */
        }  /* for j */
    }  /* for i */

    /* Compute the cloud shadow using the temperature to determine the height
       of the cloud */
    printf ("Determining cloud shadow ...\n");
    facl = cosf(xfs * DEG2RAD) * tanf(xts * DEG2RAD) / pixsize;  /* lines */
    fack = sinf(xfs * DEG2RAD) * tanf(xts * DEG2RAD) / pixsize;  /* samps */
    for (i = 0; i < nlines; i++)
    {
        curr_pix = i * nsamps;
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            if (btest (cloud[curr_pix], CLD_QA) ||
                btest (cloud[curr_pix], CIR_QA))
            {
                tcloud = sband[SR_BAND10][curr_pix] * SCALE_FACTOR_TH;
                cldh = (mclear - tcloud) * 1000.0 / cfac;
                if (cldh < 0.0)
                    cldh = 0.0;
                cldhmin = cldh - 1000.0;
                cldhmax = cldh + 1000.0;
                mband5 = 9999;
                mband5k = -9999;
                mband5l = -9999;
                if (cldhmin < 0)
                    cldhmin = 0.0;
                for (icldh = cldhmin * 0.1; icldh <= cldhmax * 0.1; icldh++)
                {
                    cldh = icldh * 10.0;
                    k = i + facl * cldh;  /* lines */
                    l = j - fack * cldh;  /* samps */
                    /* Make sure the line and sample is valid */
                    if (k < 0 || k >= nlines || l < 0 || l >= nsamps)
                        continue;

                    /* Test the pixel for cloud shadow */
                    win_pix = k * nsamps + l;
                    if ((sband[SR_BAND6][win_pix] < 800) &&
                        ((sband[SR_BAND3][win_pix] -
                          sband[SR_BAND4][win_pix]) < 100))
                    {
                        /* If it's cloud, cirrus, shadow, or fill then skip */
                        /* TODO FORTRAN code does not skip fill pixels */
                        if (btest (cloud[win_pix], CLD_QA) ||
                            btest (cloud[win_pix], CIR_QA) ||
                            btest (cloud[win_pix], CLDS_QA) ||
                            qaband[win_pix] == 1)
                        {
                            continue;
                        }
                        else
                        { /* store the value of band6 as well as the
                             l and k value */
                            if (sband[SR_BAND6][win_pix] < mband5)
                            {
                                 mband5 = sband[SR_BAND6][win_pix];
                                 mband5k = k;
                                 mband5l = l;
                            }
                        }
                    }
                }  /* for icldh */

                /* Set the cloud shadow bit */
                if (mband5 < 9999)
                    cloud[mband5k*nsamps + mband5l] += 8;
            }  /* end if btest */
        }  /* end for j */
    }  /* end for i */

    /* Expand the cloud shadow using the residual */
    printf ("Expanding cloud shadow ...\n");
    for (i = 0; i < nlines; i++)
    {
        curr_pix = i * nsamps;
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            /* If this is a cloud shadow pixel */
            if (btest (cloud[curr_pix], CLDS_QA))
            {
                /* Check the 13x13 window around the current pixel */
                for (k = i-6; k <= i+6; k++)
                {
                    /* Make sure the line is valid */
                    if (k < 0 || k >= nlines)
                        continue;

                    win_pix = k * nsamps + j-6;
                    for (l = j-6; l <= j+6; l++, win_pix++)
                    {
                        /* Make sure the sample is valid */
                        if (l < 0 || l >= nsamps)
                            continue;

                        /* If this is already cloud, cloud shadow, or fill then
                           skip */
                        /* TODO FORTRAN code does not skip fill pixels */
                        if (btest (cloud[win_pix], CLD_QA) ||
                            btest (cloud[win_pix], CLDS_QA) ||
                            qaband[win_pix] == 1)
                            continue;
                        else
                        {
                            if (btest (cloud[win_pix], CLDT_QA))
                                continue;
                            else
                            {
                                /* Set the temporary bit */
                                if (tresi[win_pix] < 0)
                                    cloud[win_pix] += 16;
                            }
                        }
                    }  /* end for l */
                }  /* end for k */
            }  /* end if btest */
        }  /* end for j */
    }  /* end for i */

    /* Update the cloud shadow */
    printf ("Updating cloud shadow ...\n");
    for (i = 0; i < nlines*nsamps; i++)
    {
        /* If the temporary bit was set in the above loop */
        if (btest (cloud[i], CLDT_QA))
        {
            /* Remove the temporary bit and set the cloud shadow bit */
            /* ==> cloud[i] += 8; cloud[i] -= 16; */
            cloud[i] -= 8;
        }
    }  /* end for i */

    /* Detect water */
    printf ("Detecting water ...\n");
    for (i = 0; i < nlines*nsamps; i++)
    {
        /* If not cloudy, cirrus, shadow, adjacent cloud, or fill then look
           for water */
        /* TODO FORTRAN code does not skip fill pixels */
        if (!btest (cloud[i], CIR_QA) &&
            !btest (cloud[i], CLD_QA) &&
            !btest (cloud[i], CLDS_QA) &&
            !btest (cloud[i], CLDA_QA) &&
            qaband[i] != 1)
        {
            /* Compute the NDVI */
            if (sband[SR_BAND5][i] < 100)
                fndvi = -0.01;
            else
                fndvi = ((double) sband[SR_BAND5][i] -
                         (double) sband[SR_BAND4][i]) /
                        ((double) sband[SR_BAND5][i] +
                         (double) sband[SR_BAND4][i]);

            /* Flag water pixels */
            if (fndvi < 0.01)
                ipflag[i] = IPFLAG_WATER;  /* water */
        }
    }  /* for i */

    /* Aerosol interpolation -- first interpolate across the line in the sample
       direction. Does not use water, cloud, or cirrus pixels. */
    printf ("Performing aerosol interpolation across each line ...\n");
    for (i = 0; i < nlines; i++)
    {
        curr_pix = i * nsamps;
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            /* If the current pixel failed aerosol retrieval and is not cloud,
               cirrus, or water then interpolate */
            if ((ipflag[curr_pix] == IPFLAG_NDVI_FAIL ||
                 ipflag[curr_pix] == IPFLAG_RESIDUAL_FAIL) &&
                !btest (cloud[curr_pix], CIR_QA) &&
                !btest (cloud[curr_pix], CLD_QA) &&
                !btest (cloud[curr_pix], WAT_QA))
            {
                /* Look at the samples in the current line and find the inferior
                   and superior boundaries where aerosol retrieval failed, but
                   the pixel hasn't yet been interpolated */
                k = j;
                win_pix = curr_pix;
                while (ipflag[win_pix] > IPFLAG_INTERP && k > 0)
                {
                    win_pix--;
                    k--;
                }
                inf_index = k;  /* inferior index */

                k = j;
                win_pix = curr_pix;
                while (ipflag[win_pix] > IPFLAG_INTERP && k < nsamps-1)
                {
                    win_pix++;
                    k++;
                }
                sup_index = k;  /* superior index */

                /* Interpolation starts one pixel after the inferior index and
                   continues until one pixel before the superior index */
                int_start = inf_index + 1;
                int_end = sup_index - 1;

                /* Make sure IPFLAG values at both of the indices are clear
                   (aerosol retrieval was successful) and not already
                   interpolated */
                inf_pix = i * nsamps + inf_index;
                sup_pix = i * nsamps + sup_index;
                if (ipflag[inf_pix] > IPFLAG_INTERP ||
                    ipflag[sup_pix] > IPFLAG_INTERP)
                    continue;

                /* Interpolate the aerosols across the line */
                win_pix = i * nsamps + int_start;
                for (k = int_start; k <= int_end; k++, win_pix++)
                {
                    /* If this pixel is not fill or water and therefore is a
                       pixel where the aerosol retrieval failed then
                       interpolate the aerosol value */
                    if (ipflag[win_pix] != IPFLAG_FILL &&
                        ipflag[win_pix] != IPFLAG_WATER)
                    {
                        taero[win_pix] = taero[inf_pix] +
                            (taero[sup_pix] - taero[inf_pix]) *
                            (k - inf_index) / (sup_index - inf_index);
                        ipflag[win_pix] = IPFLAG_INTERP;
                    }
                }
            }  /* end if */
        }  /* end for j < nsamps */
    }  /* end for i < nlines */

    /* Aerosol interpolation -- next interpolate down the samples in the line
       direction. Does not use water, cloud, or cirrus pixels. */
    printf ("Performing aerosol interpolation down each sample ...\n");
    for (j = 0; j < nsamps; j++)
    {
        for (i = 0; i < nlines; i++)
        {
            curr_pix = i * nsamps + j;
            /* If the current pixel failed aerosol retrieval and is not cloud,
               cirrus, or water then interpolate */
            if ((ipflag[curr_pix] == IPFLAG_NDVI_FAIL ||
                 ipflag[curr_pix] == IPFLAG_RESIDUAL_FAIL) &&
                !btest (cloud[curr_pix], CIR_QA) &&
                !btest (cloud[curr_pix], CLD_QA) &&
                !btest (cloud[curr_pix], WAT_QA))
            {
                /* Look at the line values for the current sample and find the
                   inferior and superior boundaries where aerosol retrieval
                   failed, but the pixel hasn't yet been interpolated */
                l = i;
                win_pix = curr_pix;
                while (ipflag[win_pix] > IPFLAG_INTERP && l > 0)
                {
                    l--;
                    win_pix = l * nsamps + j;
                }
                inf_index = l;  /* inferior index */

                l = i;
                win_pix = curr_pix;
                while (ipflag[win_pix] > IPFLAG_INTERP && l < nlines-1)
                {
                    l++;
                    win_pix = l * nsamps + j;
                }
                sup_index = l;  /* superior index */

                /* Interpolation starts one pixel after the inferior index and
                   continues until one pixel before the superior index */
                int_start = inf_index + 1;
                int_end = sup_index - 1;

                /* Make sure IPFLAG values at both of the indices are clear
                   (aerosol retrieval was successful) and not already
                   interpolated */
                inf_pix = inf_index * nsamps + j;
                sup_pix = sup_index * nsamps + j;
                if (ipflag[inf_pix] > IPFLAG_INTERP ||
                    ipflag[sup_pix] > IPFLAG_INTERP)
                    continue;

                /* Interpolate the aerosols down the samples */
                for (l = int_start; l <= int_end; l++)
                {
                    /* If this pixel is not fill or water and therefore is a
                       pixel where the aerosol retrieval failed then
                       interpolate the aerosol value */
                    win_pix = l * nsamps + j;
                    if (ipflag[win_pix] != IPFLAG_FILL &&
                        ipflag[win_pix] != IPFLAG_WATER)
                    {
                        taero[win_pix] = taero[inf_pix] +
                            (taero[sup_pix] - taero[inf_pix]) *
                            (l - inf_index) / (sup_index - inf_index);
                        ipflag[win_pix] = IPFLAG_INTERP;
                    }
                }
            }  /* end if */
        }  /* end for i < nlines */
    }  /* end for j < nsamps */

    /* Perform the second level of atmospheric correction for the aerosols.
       This is not applied to cirrus or cloud pixels. */
    printf ("Performing atmospheric correction ...\n");
    /* 0 .. DN_BAND7 is the same as 0 .. SR_BAND7 here, since the pan band
       isn't spanned */
    for (ib = 0; ib <= DN_BAND7; ib++)
    {
        printf ("  Band %d\n", ib+1);
#ifdef _OPENMP
        #pragma omp parallel for private (i, rsurf, rotoa, raot550nm, pres, uwv, uoz, retval, roslamb, tgo, roatm, ttatmg, satm, xrorayp, next, isuccess, ros2b1)
#endif
        for (i = 0; i < nlines*nsamps; i++)
        {
            /* If this pixel is fill, then don't process. Otherwise the
               fill pixels have already been marked in the TOA process. */
            if (qaband[i] != 1)
            {
                /* If this isn't a cirrus or cloud pixel */
                if (!btest (cloud[i], CIR_QA) && !btest (cloud[i], CLD_QA))
                {
                    rsurf = sband[ib][i] * SCALE_FACTOR;
                    rotoa = (rsurf * bttatmg[ib] / (1.0 - bsatm[ib] * rsurf)
                        + broatm[ib]) * btgo[ib];
                    raot550nm = taero[i];
                    pres = tp[i];
                    uwv = twvi[i];
                    uoz = tozi[i];
                    retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                        raot550nm, ib, pres, tpres, aot550nm, rolutt,
                        transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                        normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
                        ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                        ogtransb1, wvtransa, wvtransb, oztransa, rotoa,
                        &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp,
                        &next);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing lambertian "
                            "atmospheric correction type 2.");
                        error_handler (true, FUNC_NAME, errmsg);
                        exit (ERROR);
                    }

                    /* If this is the coastal aerosol band then set the
                       aerosol bits in the QA band */
                    if (ib == DN_BAND1)
                    {
                        /* Recompute based on predefined taero value */
                        ros2b1 = roslamb;
                        if (roslamb < 0.01 && raot550nm > 0.06)
                        {
                            /* This should be refined to get better values */
                            isuccess = false;
                            while (roslamb < 0.0)
                            {
                                raot550nm = taero[i] - 0.05;
                                ros2b1 = roslamb;
                                pres = tp[i];
                                uwv = twvi[i];
                                uoz = tozi[i];
                                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi,
                                    cosxfi, raot550nm, ib, pres, tpres,
                                    aot550nm, rolutt, transt, xtsstep, xtsmin,
                                    xtvstep, xtvmin, sphalbt, normext, tsmax,
                                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
                                    uwv, tauray, ogtransa1, ogtransb0,
                                    ogtransb1, wvtransa, wvtransb, oztransa,
                                    rotoa, &roslamb, &tgo, &roatm, &ttatmg,
                                    &satm, &xrorayp, &next);
                                if (retval != SUCCESS)
                                {
                                    sprintf (errmsg, "Performing lambertian "
                                        "atmospheric correction type 2.");
                                    error_handler (true, FUNC_NAME, errmsg);
                                    exit (ERROR);
                                }

                                /* Save the results and recheck */
                                taero[i] = raot550nm;
                                if (roslamb >= 0.0)
                                    isuccess = true;
                            }  /* end while roslamb < 0.0 */

                            /* If not successful, run one more time */
                            if (!isuccess)
                            {
                                raot550nm = taero[i] - 0.05;
                                ros2b1 = roslamb;
                                pres = tp[i];
                                uwv = twvi[i];
                                uoz = tozi[i];
                                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi,
                                    cosxfi, raot550nm, ib, pres, tpres,
                                    aot550nm, rolutt, transt, xtsstep, xtsmin,
                                    xtvstep, xtvmin, sphalbt, normext, tsmax,
                                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
                                    uwv, tauray, ogtransa1, ogtransb0,
                                    ogtransb1, wvtransa, wvtransb, oztransa,
                                    rotoa, &roslamb, &tgo, &roatm, &ttatmg,
                                    &satm, &xrorayp, &next);
                                if (retval != SUCCESS)
                                {
                                    sprintf (errmsg, "Performing lambertian "
                                        "atmospheric correction type 2.");
                                    error_handler (true, FUNC_NAME, errmsg);
                                    exit (ERROR);
                                }

                                /* Save the results */
                                taero[i] = raot550nm;
                            }

                            /* Refined optical depth */
                            raot550nm = taero[i] - (roslamb - 0.01) *
                                 (-0.05 / (roslamb - ros2b1));
                            taero[i] = raot550nm;
                            retval = atmcorlamb2 (xts, xtv, xmus, xmuv,
                                xfi, cosxfi, raot550nm, ib, pres, tpres,
                                aot550nm, rolutt, transt, xtsstep, xtsmin,
                                xtvstep, xtvmin, sphalbt, normext, tsmax,
                                tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
                                uwv, tauray, ogtransa1, ogtransb0,
                                ogtransb1, wvtransa, wvtransb, oztransa,
                                rotoa, &roslamb, &tgo, &roatm, &ttatmg,
                                &satm, &xrorayp, &next);
                            if (retval != SUCCESS)
                            {
                                sprintf (errmsg, "Performing lambertian "
                                    "atmospheric correction type 2.");
                                error_handler (true, FUNC_NAME, errmsg);
                                exit (ERROR);
                            }
                        }  /* if roslamb && raot550nm */

                        /* Set up aerosol QA bits */
                        tmpf = fabs (rsurf - roslamb);
                        if (tmpf <= 0.015)
                        {  /* Set the first aerosol bit (low aerosols) */
                            cloud[i] += 16;
                        }
                        else
                        {
                            if (tmpf < 0.03)
                            {  /* Set the second aerosol bit (average
                                  aerosols) */
                                cloud[i] += 32;
                            }
                            else
                            {  /* Set both aerosol bits (high aerosols) */
                                cloud[i] += 48;
                            }
                        }
                    }  /* end if this is the coastal aerosol band */

                    /* Save the scaled surface reflectance value, but make
                       sure it falls within the defined valid range. */
                    roslamb = roslamb * MULT_FACTOR;  /* scale the value */
                    if (roslamb < MIN_VALID)
                        sband[ib][i] = MIN_VALID;
                    else if (roslamb > MAX_VALID)
                        sband[ib][i] = MAX_VALID;
                    else
                          /* TODO FORTRAN code doesn't round here, but it's
                             probably a good idea */
//                        sband[ib][i] = (int) (roundf (roslamb));
                        sband[ib][i] = (int) (roslamb);
                }  /* end if not cirrus and not cloud */
            }  /* end if not fill */
        }  /* end for i */
    }  /* end for ib */

    /* Free memory for band data */
    free (twvi);
    free (tozi);
    free (tp);
    free (tresi);
    free (taero);
 
    /* Write the data to the output file */
    printf ("Writing surface reflectance corrected data to the output "
        "files ...\n");

    /* Open the output file */
    sr_output = open_output (xml_metadata, input, false /*surf refl*/);
    if (sr_output == NULL)
    {   /* error message already printed */
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Loop through the reflectance bands and write the data */
    for (ib = 0; ib <= DN_BAND7; ib++)
    {
        printf ("  Band %d: %s\n", ib+1,
            sr_output->metadata.band[ib].file_name);
        if (put_output_lines (sr_output, sband[ib], ib, 0, nlines,
            sizeof (int16)) != SUCCESS)
        {
            sprintf (errmsg, "Writing output data for band %d", ib);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Create the ENVI header file this band */
        if (create_envi_struct (&sr_output->metadata.band[ib],
            &xml_metadata->global, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Creating ENVI header structure.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Write the ENVI header */
        strcpy (envi_file, sr_output->metadata.band[ib].file_name);
        cptr = strchr (envi_file, '.');
        strcpy (cptr, ".hdr");
        if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Writing ENVI header file.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Append the surface reflectance bands (1-7) to the XML file */
    if (append_metadata (7, sr_output->metadata.band, xml_infile) !=
        SUCCESS)
    {
        sprintf (errmsg, "Appending surface reflectance bands to the "
            "XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Write the cloud mask band */
    printf ("  Band %d: %s\n", SR_CLOUD+1,
            sr_output->metadata.band[SR_CLOUD].file_name);
    if (put_output_lines (sr_output, cloud, SR_CLOUD, 0, nlines,
        sizeof (uint8)) != SUCCESS)
    {
        sprintf (errmsg, "Writing cloud mask output data");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Free memory for cloud data */
    free (cloud);

    /* Create the ENVI header for the cloud mask band */
    if (create_envi_struct (&sr_output->metadata.band[SR_CLOUD],
        &xml_metadata->global, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Creating ENVI header structure.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Write the ENVI header */
    strcpy (envi_file, sr_output->metadata.band[SR_CLOUD].file_name);
    cptr = strchr (envi_file, '.');
    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Writing ENVI header file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Append the cloud mask band to the XML file */
    if (append_metadata (1, &sr_output->metadata.band[SR_CLOUD],
        xml_infile) != SUCCESS)
    {
        sprintf (errmsg, "Appending cloud mask band to XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Write the ipflag mask band */
    printf ("  Band %d: %s\n", SR_IPFLAG+1,
            sr_output->metadata.band[SR_IPFLAG].file_name);
    if (put_output_lines (sr_output, ipflag, SR_IPFLAG, 0, nlines,
        sizeof (uint8)) != SUCCESS)
    {
        sprintf (errmsg, "Writing ipflag mask output data");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Free memory for ipflag data */
    free (ipflag);

    /* Create the ENVI header for the ipflag mask band */
    if (create_envi_struct (&sr_output->metadata.band[SR_IPFLAG],
        &xml_metadata->global, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Creating ENVI header structure.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Write the ENVI header */
    strcpy (envi_file, sr_output->metadata.band[SR_IPFLAG].file_name);
    cptr = strchr (envi_file, '.');
    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Writing ENVI header file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Append the ipflag mask band to the XML file */
    if (append_metadata (1, &sr_output->metadata.band[SR_IPFLAG],
        xml_infile) != SUCCESS)
    {
        sprintf (errmsg, "Appending ipflag mask band to XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Close the output surface reflectance products */
    close_output (sr_output, false /*sr products*/);
    free_output (sr_output);

    /* Free the spatial mapping pointer */
    free (space);

    /* Free the data arrays */
    for (i = 0; i < NSR_BANDS; i++)
    {
        for (j = 0; j < 7; j++)
        {
            for (k = 0; k < 22; k++)
            {
                free (rolutt[i][j][k]);
                free (transt[i][j][k]);
            }
            free (rolutt[i][j]);
            free (transt[i][j]);
            free (sphalbt[i][j]);
        }
        free (rolutt[i]);
        free (transt[i]);
        free (sphalbt[i]);
    }
    free (rolutt);
    free (transt);
    free (sphalbt);

    /* tsmax[20][22] and float tsmin[20][22] and float nbfic[20][22] and
       nbfi[20][22] and float ttv[20][22] */
    for (i = 0; i < 20; i++)
    {
        free (tsmax[i]);
        free (tsmin[i]);
        free (nbfic[i]);
        free (nbfi[i]);
        free (ttv[i]);
    }
    free (tsmax);
    free (tsmin);
    free (nbfic);
    free (nbfi);
    free (ttv);

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  init_sr_refl

PURPOSE:  Initialization for the atmospheric corrections.  Initialization for
look up tables, auxiliary data, mapping, and geolocation information is used
for the surface reflectance correction.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error initializing the atmospheric parameters
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
1. The view angle is set to 0.0 and this never changes.
2. The DEM is used to calculate the surface pressure.
******************************************************************************/
int init_sr_refl
(
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    Input_t *input,     /* I: input structure for the Landsat product */
    Geoloc_t *space,    /* I: structure for geolocation information */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm,        /* I: auxiliary filename for ozone and water vapor */
    float *xtv,         /* O: observation zenith angle (deg) */
    float *xmuv,        /* O: cosine of observation zenith angle */
    float *xfi,         /* O: azimuthal difference between sun and
                              observation (deg) */
    float *cosxfi,      /* O: cosine of azimuthal difference */
    float *raot550nm,   /* O: nearest value of AOT */
    float *pres,        /* O: surface pressure */
    float *uoz,         /* O: total column ozone */
    float *uwv,         /* O: total column water vapor (precipital water
                              vapor) */
    float *xtsstep,     /* O: solar zenith step value */
    float *xtsmin,      /* O: minimum solar zenith value */
    float *xtvstep,     /* O: observation step value */
    float *xtvmin,      /* O: minimum observation value */
    float **tsmax,      /* O: maximum scattering angle table [20][22] */
    float **tsmin,      /* O: minimum scattering angle table [20][22] */
    float tts[22],      /* O: sun angle table */
    float **ttv,        /* O: view angle table [20][22] */
    int32 indts[22],    /* O: index for the sun angle table */
    float ****rolutt,   /* O: intrinsic reflectance table
                              [NSR_BANDS][7][22][8000] */
    float ****transt,   /* O: transmission table [NSR_BANDS][7][22][22] */
    float ***sphalbt,   /* O: spherical albedo table [NSR_BANDS][7][22] */
    float ***normext,   /* O: aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm)
                              [NSR_BANDS][7][22] */
    float **nbfic,      /* O: communitive number of azimuth angles [20][22] */
    float **nbfi,       /* O: number of azimuth angles [20][22] */
    int16 *dem,         /* O: CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 *andwi,       /* O: avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *sndwi,       /* O: standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob1,     /* O: mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob2,     /* O: mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob7,     /* O: mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob1,  /* O: integer band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob2,  /* O: integer band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob7,  /* O: integer band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob1,  /* O: slope band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob2,  /* O: slope band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob7,  /* O: slope band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    uint16 *wv,         /* O: water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 *oz           /* O: ozone values [CMG_NBLAT x CMG_NBLON] */
)
{
    char errmsg[STR_SIZE];                   /* error message */
    char FUNC_NAME[] = "init_sr_refl";       /* function name */
    int retval;          /* return status */
    int lcmg, scmg;      /* line/sample index for the CMG */
    int cmg_pix;         /* pixel location in the CMG array for [lcmg][scmg] */
    int dem_pix;         /* pixel location in the DEM array for [lcmg][scmg] */
    float xcmg, ycmg;    /* x/y location for CMG */

    /* Vars for forward/inverse mapping space */
    Img_coord_float_t img;        /* coordinate in line/sample space */
    Geo_coord_t geo;              /* coordinate in lat/long space */
    float center_lat, center_lon; /* lat/long for scene center */

    /* Initialize the look up tables */
    *xtv = 0.0;
    *xmuv = cos (*xtv * DEG2RAD);
    *xfi = 0.0;
    *cosxfi = cos (*xfi * DEG2RAD);
    *xtsmin = 0;
    *xtsstep = 4.0;
    *xtvmin = 2.84090;
    *xtvstep = 6.52107 - *xtvmin;
    retval = readluts (tsmax, tsmin, ttv, tts, nbfi, nbfic, indts, rolutt,
        transt, sphalbt, normext, *xtsstep, *xtsmin, anglehdf, intrefnm,
        transmnm, spheranm);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Reading the LUTs");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    printf ("The LUTs for urban clean case v2.0 have been read.  We can "
        "now perform atmospheric correction.\n");

    /* Read the auxiliary data files used as input to the reflectance
       calculations */
    retval = read_auxiliary_files (anglehdf, intrefnm, transmnm, spheranm,
        cmgdemnm, rationm, auxnm, dem, andwi, sndwi, ratiob1, ratiob2,
        ratiob7, intratiob1, intratiob2, intratiob7, slpratiob1,
        slpratiob2, slpratiob7, wv, oz);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Reading the auxiliary files");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Getting parameters for atmospheric correction */
    /* Update to get the parameter of the scene center */
    *raot550nm = 0.12;
    *pres = 1013.0;
    *uoz = 0.30;
    *uwv = 0.5;

    /* Use scene center (and center of the pixel) to compute atmospheric
       parameters */
    /* TODO -- FORTRAN code adds one to the calculation.  Is this correct?? */
    img.l = (int) (1 + nlines * 0.5) - 0.5;
    img.s = (int) (1 + nsamps * 0.5) + 0.5;
    img.is_fill = false;
    if (!from_space (space, &img, &geo))
    {
        sprintf (errmsg, "Mapping scene center to geolocation coords");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    center_lat = geo.lat * RAD2DEG;
    center_lon = geo.lon * RAD2DEG;
    printf ("Scene center line/sample: %f, %f\n", img.l, img.s);
    printf ("Scene center lat/long: %f, %f\n", center_lat, center_lon);

    /* Use the scene center lat/long to determine the line/sample in the
       CMG-related lookup tables, using the center of the UL pixel */
    ycmg = (89.975 - center_lat) * 20.0;    /* vs / 0.05 */
    xcmg = (179.975 + center_lon) * 20.0;   /* vs / 0.05 */
    lcmg = (int) roundf (ycmg);
    scmg = (int) roundf (xcmg);
    if ((lcmg < 0 || lcmg >= CMG_NBLAT) || (scmg < 0 || scmg >= CMG_NBLON))
    {
        sprintf (errmsg, "Invalid line/sample combination for the "
            "CMG-related lookup tables - line %d, sample %d (0-based).  "
            "CMG-based tables are %d lines x %d samples.", lcmg, scmg,
            CMG_NBLAT, CMG_NBLON);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    cmg_pix = lcmg * CMG_NBLON + scmg;
    if (wv[cmg_pix] != 0)
        *uwv = wv[cmg_pix] / 200.0;
    else
        *uwv = 0.5;

    if (oz[cmg_pix] != 0)
        *uoz = oz[cmg_pix] / 400.0;
    else
        *uoz = 0.3;

    dem_pix = lcmg * DEM_NBLON + scmg;
    if (dem[dem_pix] != -9999)
        *pres = 1013.0 * exp (-dem[dem_pix] * ONE_DIV_8500);
    else
        *pres = 1013.0;
    *raot550nm = 0.05;

    /* Successful completion */
    return (SUCCESS);
}
