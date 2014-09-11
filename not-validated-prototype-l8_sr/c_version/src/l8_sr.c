#include <sys/stat.h>
#include <unistd.h>
#include "l8_sr.h"

/******************************************************************************
MODULE:  l8_sr

PURPOSE:  Computes the surface reflectance values for the Landsat 8 products.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing of the surface reflectance
SUCCESS         Processing was successful

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date          Programmer       Reason
----------    ---------------  -------------------------------------
6/30/2014     Gail Schmidt     Conversion of the original L8 SR code delivered
                               by Eric Vermote to adhere to ESPA software
                               guidelines and file formats
7/15/2014     Gail Schmidt     Cleaned up some of the divide by constants so
                               that we could speed up the processing
7/21/2014     Gail Schmidt     Set all int16 output bands to use the same fill
                               value, which is a change for bands 10 and 11.
                               The fill value used is -9999 vs. -1000.
7/22/2014     Gail Schmidt     Made the tauray array a hard-coded array vs.
                               reading it from a static ASCII file.
7/22/2014     Gail Schmidt     Cleaned up unused ogtransa0, ogtransc0,
                               ogtransc1, wvtransc arrays.  Made the rest of
                               these transmission arrays doubles and hard-coded
                               their static values in this code vs. reading
                               from a static ASCII file.
7/22/2014     Gail Schmidt     Changed the 2D arrays to 1D arrays for the
                               image data to speed up processing.
7/29/2014     Gail Schmidt     Defined a static NSR_BANDS variable for the
                               variables that refer to the surface reflectance
                               band-related bands (ogtrans, wvtrans, tauray,
                               erelc, etc.)  These previously were of size 16.
8/1/2014      Gail Schmidt     Add check on the solar zenith to make sure the
                               scene can be processed for surface reflectance.
                               If solar zenith is too large, then only process
                               TOA reflectance.
8/1/2014      Gail Schmidt     Added flag to allow user to specify only TOA
                               reflectance corrections to be completed and
                               written.  Also added flag to allow the user to
                               specify TOA reflectance bands (bands 1-7) should
                               be written in addition to SR bands.
8/14/2014    Gail Schmidt      Updated for v1.3 delivered by Eric Vermote
8/25/2014    Gail Schmidt      Split the main application into smaller modules
                               for allocating memory and reading the auxiliary
                               data files.

NOTES:
1. Bands 1-7 are corrected to surface reflectance.  Band 8 (pand band) is not
   processed.  Band 9 (cirrus band) is corrected to TOA reflectance.  Bands
   10 and 11 are corrected to brightness temperature.
2. The TOA reflectance corrections are made with a correction for the sun angle.
   The sun angle correction is currently only made based on the angle at the
   center of the scene, not the per pixel angle.
3. SDstart and SDreaddata have minor memory leaks.  Ultimately both call
   HAregister_atom which makes a malloc call and the memory is never freed.
4. Conversion algorithms for TOA reflectance and at-sensor brightness
   temperature are available from
   http://landsat.usgs.gov/Landsat8_Using_Product.php
******************************************************************************/
int main (int argc, char *argv[])
{
    bool verbose;            /* verbose flag for printing messages */
    char FUNC_NAME[] = "main"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char envi_file[STR_SIZE];/* ENVI filename */
    char *aux_path = NULL;   /* path for Landsat auxiliary data */
    char *xml_infile = NULL; /* input XML filename */
    char *aux_infile = NULL; /* input auxiliary filename for water vapor
                                and ozone*/
    char *cptr = NULL;       /* pointer to the file extension */
    char aux_year[5];        /* string to contain the year of auxiliary file */

    int retval;              /* return status */
    int ib;                  /* looping variable for input bands */
    int sband_ib;            /* looping variable for output bands */
    int curr_pix;            /* current pixel in the 1D arrays of
                                nlines * nsamps */      
    int win_pix;             /* current pixel in the line,sample window */
    Input_t *input = NULL;       /* input structure for the Landsat product */
    Output_t *sr_output = NULL;  /* output structure and metadata for the SR
                                    product */
    Output_t *toa_output = NULL; /* output structure and metadata for the TOA
                                    product */
    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
    Espa_global_meta_t *gmeta = NULL;   /* pointer to global meta */
    Envi_header_t envi_hdr;      /* output ENVI header information */

    /* Vars for forward/inverse mapping space */
    Geoloc_t *space = NULL;       /* structure for geolocation information */
    Space_def_t space_def;        /* structure to define the space mapping */
    Img_coord_float_t img;        /* coordinate in line/sample space */
    Geo_coord_t geo;              /* coordinate in lat/long space */
    float center_lat, center_lon; /* lat/long for scene center */
    float lat, lon;               /* pixel lat, long location */

    struct stat statbuf;      /* buffer for the file stat function */
    uint16 *uband = NULL;     /* array of input image data for a current band,
                                 nlines x nsamps */
    uint16 *qaband = NULL;    /* QA band for the input image, nlines x nsamps */
    int16 *aerob1 = NULL;     /* atmospherically corrected band 1 data
                                 (TOA refl), nlines x nsamps */
    int16 *aerob2 = NULL;     /* atmospherically corrected band 2 data
                                 (TOA refl), nlines x nsamps */
    int16 *aerob4 = NULL;     /* atmospherically corrected band 4 data
                                 (TOA refl), nlines x nsamps */
    int16 *aerob5 = NULL;     /* atmospherically corrected band 5 data
                                 (TOA refl), nlines x nsamps */
    int16 *aerob7 = NULL;     /* atmospherically corrected band 7 data
                                 (TOA refl), nlines x nsamps */
    int16 **sband = NULL;     /* output surface reflectance and brightness
                                 temp bands, qa band is separate as a uint16 */
    int16 **dem = NULL;       /* CMG DEM data array [DEM_NBLAT][DEM_NBLON] */
    int16 **andwi = NULL;     /* avg NDWI [RATIO_NBLAT][RATIO_NBLON] */
    int16 **sndwi = NULL;     /* standard NDWI [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob1 = NULL;   /* mean band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob2 = NULL;   /* mean band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob7 = NULL;   /* mean band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob1 = NULL;   /* ??? band1 ratio
                                    [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob2 = NULL;   /* ??? band2 ratio
                                    [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob7 = NULL;   /* ??? band7 ratio
                                    [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob1 = NULL;   /* slope band1 ratio
                                    [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob2 = NULL;   /* slope band2 ratio
                                    [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob7 = NULL;   /* slope band7 ratio
                                    [RATIO_NBLAT][RATIO_NBLON] */
    uint16 **wv = NULL;       /* water vapor values [CMG_NBLAT][CMG_NBLON] */
    uint8 **oz = NULL;        /* ozone values [CMG_NBLAT][CMG_NBLON] */
    uint8 *cloud = NULL;      /* bit-packed value that represent clouds,
                                 nlines x nsamps */
    float *twvi = NULL;       /* interpolated water vapor value,
                                 nlines x nsamps */
    float *tozi = NULL;       /* interpolated ozone value, nlines x nsamps */
    float *tp = NULL;         /* interpolated pressure value, nlines x nsamps */
    float *tresi = NULL;      /* residuals for each pixel, nlines x nsamps */
    float *taero = NULL;      /* aerosol values for each pixel,
                                 nlines x nsamps */

    int i, j, k, l;      /* looping variables */
    int lcmg, scmg;      /* line/sample index for the CMG */
    int iband;           /* current band */
    float u, v;
    float th1, th2;      /* values for NDWI calculations */
    float xcmg, ycmg;    /* x/y location for CMG */
    float xndwi;         /* calculated NDWI value */
    float xts;           /* solar zenith angle (deg) */
    float xfs;           /* solar azimuth angle (deg) */
    float xtv;           /* observation zenith angle (deg) */
    float xfi;           /* azimuthal difference between the sun and
                            observation angle (deg) */
    float xtsstep;       /* solar zenith step value */
    float xtsmin;        /* minimum solar zenith value */
    float xtvstep;       /* observation step value */
    float xtvmin;        /* minimum observation value */
    bool process_sr = true;  /* this is set to false if the solar zenith
                                is too large and the surface reflectance
                                cannot be calculated or if the user specifies
                                that surface reflectance processing will not
                                be completed and only TOA processing will be
                                done */
    bool write_toa = false;  /* this is set to true if the user specifies
                                TOA products should be output for delivery */

    /* The following arguments are all names of the LUTs */
    char anglehdf[STR_SIZE];  /* angle HDF filename */
    char intrefnm[STR_SIZE];  /* intrinsic reflectance filename */
    char transmnm[STR_SIZE];  /* transmission filename */
    char spheranm[STR_SIZE];  /* spherical albedo filename */
    char cmgdemnm[STR_SIZE];  /* climate modeling grid DEM filename */
    char rationm[STR_SIZE];   /* ratio averages filename */
    char auxnm[STR_SIZE];     /* auxiliary filename for ozone and water vapor*/

    /* Atmospheric correction variables */
    /* Look up table for atmospheric and geometric quantities */
    float tauray[NSR_BANDS] =  /* molecular optical thickness coeff */
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
    float **nbfic = NULL;       /* communitive number of azimuth angles
                                   [20][22] */
    float **nbfi = NULL;        /* number of azimuth angles [20][22] */
    float **ttv = NULL;         /* view angle table [20][22] */
    float tts[22];              /* sun angle table */
    int32 indts[22];
    float aot550nm[22] =  /* AOT look-up table */
        {0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 1.20,
         1.40, 1.60, 1.80, 2.00, 2.30, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00};
    float tpres[7] =      /* surface pressure table */
        {1050.0, 1013.0, 900.0, 800.0, 700.0, 600.0, 500.0};

    float rotoa;        /* top of atmosphere reflectance */
    float raot550nm;    /* nearest input value of AOT */
    float uoz;          /* total column ozone */
    float uwv;          /* total column water vapor (precipital water vapor) */
    float pres;         /* surface pressure */
    float roslamb;      /* lambertian surface reflectance */
    float tgo;          /* other gaseous transmittance */
    float roatm;        /* atmospheric reflectance */
    float ttatmg;
    float satm;         /* spherical albedo */
    float xrorayp;      /* molecular reflectance */
    float next;         /* ???? */

    float pixsize;      /* pixel size for the reflectance files */
    int nlines, nsamps; /* number of lines and samples in the reflectance and
                           thermal bands */
    int row, col;       /* pixel row, column (line, sample) location */
    int colp, rowp;     /* row/col for true north adjustment */
    int uoz11, uoz21, uoz12, uoz22;  /* ozone at line,samp; line, samp+1;
                           line+1, samp; and line+1, samp+1 */
    float pres11, pres12, pres21, pres22;  /* pressure at line,samp;
                           line, samp+1; line+1, samp; and line+1, samp+1 */
    float erelc[NSR_BANDS];    /* band ratio variable for bands 1-7 */
    float troatm[NSR_BANDS];   /* atmospheric reflectance table for bands 1-7 */
    float btgo[NSR_BANDS];     /* other gaseous transmittance for bands 1-7 */
    float broatm[NSR_BANDS];   /* atmospheric reflectance for bands 1-7 */
    float bttatmg[NSR_BANDS];  /* ttatmg for bands 1-7 */
    float bsatm[NSR_BANDS];    /* spherical albedo for bands 1-7 */
    int iband1, iband3; /* band indices (zero-based) */
    float raot;
    float residual;     /* model residual */
    float rsurf;
    float xmus;         /* cosine of solar zenith */
    float corf;
    float tmpf;         /* temporary floating point value */

    /* K[1|2]b1[0|1] constants, additive, and multiplier are also found in the
       MTL file */
    const float xcals = 3.3420E-04;  /* radiance multiplier for bands
                                        10 and 11 */
    const float xcalo = 0.10000;     /* radiance additive for bands
                                        10 and 11 */
    const float refl_mult = 2.0E-05; /* reflectance multiplier for bands 1-9 */
    const float refl_add = -0.1;     /* reflectance additive for bands 1-9 */
    const float k1b10 = 774.89;      /* temperature constant for band 10 */
    const float k1b11 = 480.89;      /* temperature constant for band 11 */
    const float k2b10 = 1321.08;     /* temperature constant for band 10 */
    const float k2b11 = 1201.14;     /* temperature constant for band 11 */

    float dy, dx;                    /* delta x/y for true north adjustment */
    float ang;                       /* angle for true north adjustment */
    long nbclear;                    /* count of the clear (non-cloud) pixels */
    long nbval;                      /* count of the non-fill pixels */
    double anom;                     /* band 3 and 5 combination */
    double mall;                     /* average/mean temp of all the pixels */
    double mclear;                   /* average/mean temp of the clear pixels */
    float fack, facl;                /* cloud height factor in the k,l dim */
    int cldhmin;                     /* minimum bound of the cloud height */
    int cldhmax;                     /* maximum bound of the cloud height */
    float cldh;                      /* cloud height */
    int icldh;                       /* looping variable for cloud height */
    int mband5, mband5k, mband5l;    /* band 6 value and k,l locations */
    float tcloud;                    /* temperature of the current pixel */

    float cfac = 6.0;  /* cloud factor */
    double aaot;
    double sresi;      /* sum of 1 / residuals */
    float fndvi;       /* NDVI value */
    int nbaot;
    int step;
    int hole;
    float ros4, ros5;    /* surface reflectance for band 4 and band 5 */
    int tmp_percent;     /* current percentage for printing status */

    printf ("Starting TOA and surface reflectance processing ...\n");

    /* Read the command-line arguments */
    retval = get_args (argc, argv, &xml_infile, &aux_infile, &process_sr,
        &write_toa, &verbose);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }

    /* Provide user information if verbose is turned on */
    if (verbose)
    {
        printf ("  XML input file: %s\n", xml_infile);
        printf ("  AUX input file: %s\n", aux_infile);
        if (!process_sr)
        {
            printf ("    **Surface reflectance corrections will not be "
                "completed.  Only top of atmosphere corrections will be "
                "completed.\n");
        }
    }

    /* Validate the input metadata file */
    if (validate_xml_file (xml_infile) != SUCCESS)
    {  /* Error messages already written */
        exit (ERROR);
    }

    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (xml_infile, &xml_metadata) != SUCCESS)
    {  /* Error messages already written */
        exit (ERROR);
    }

    /* Open the reflectance product, set up the input data structure, and
       allocate memory for the data buffers */
    input = open_input (&xml_metadata);
    if (input == (Input_t *) NULL)
    {
        sprintf (errmsg, "Error opening/reading the input DN data: %s",
            xml_infile);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    gmeta = &xml_metadata.global;

    /* Output some information from the input files if verbose */
    if (verbose)
    {
        printf ("  Number of lines/samples: %d/%d\n", input->size.nlines,
            input->size.nsamps);
        printf ("  Nband: %d\n", input->nband);
        printf ("  Pixsize: %f,%f\n", input->size.pixsize[0],
            input->size.pixsize[1]);
        printf ("  Number of thermal lines/samples: %d/%d\n",
            input->size_th.nlines, input->size_th.nsamps);
        printf ("  Nband thermal: %d\n", input->nband_th);
        printf ("  Pixsize: %f,%f\n", input->size_th.pixsize[0],
            input->size_th.pixsize[1]);
        printf ("  Number of qa lines/samples: %d/%d\n",
            input->size_qa.nlines, input->size_qa.nsamps);
        printf ("  Nband QA: %d\n", input->nband_qa);
        printf ("  Pixsize: %f,%f\n", input->size_qa.pixsize[0],
            input->size_qa.pixsize[1]);
        printf ("  Fill value: %d\n", input->meta.fill);
        printf ("  Solar zenith: %f\n", xml_metadata.global.solar_zenith);
        printf ("  Solar azimuth: %f\n", xml_metadata.global.solar_azimuth);
    }

    /* Pull the needed metadata from the XML file and input structure */
    xts = gmeta->solar_zenith;
    xfs = gmeta->solar_azimuth;
    pixsize = (float) input->size.pixsize[0];
    raot550nm = 0.06;
    xmus = cos (xts * DEG2RAD);
    nlines = input->size.nlines;
    nsamps = input->size.nsamps;

    /* The surface reflectance algorithm cannot be implemented for solar
       zenith angles greater than 76 degrees.  Need to flag if the current
       scene falls into that category. */
    if (xts > 76.0)
    {
        process_sr = false;
        sprintf (errmsg, "Solar zenith angle is too large to allow for surface "
            "reflectance processing.  Corrections will be limited to top of "
            "atmosphere and at-sensor brightness temperature corrections.");
        error_handler (false, FUNC_NAME, errmsg);
    }

    /* Allocate memory for all the data arrays */
    if (verbose)
        printf ("Allocating memory for the data arrays ...\n");
    retval = memory_allocation (nlines, nsamps, &dem, &andwi, &sndwi, &ratiob1,
        &ratiob2, &ratiob7, &intratiob1, &intratiob2, &intratiob7, &slpratiob1,
        &slpratiob2, &slpratiob7, &wv, &oz, &rolutt, &transt, &sphalbt,
        &normext, &tsmax, &tsmin, &nbfic, &nbfi, &ttv, &uband, &qaband,
        &aerob1, &aerob2, &aerob4, &aerob5, &aerob7, &sband, &cloud, &twvi,
        &tozi, &tp, &tresi, &taero);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        sprintf (errmsg, "Error allocating memory for the data arrays.");
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Get the path for the auxiliary products from the L8_AUX_DIR environment
       variable.  If it isn't defined, then assume the products are in the
       local directory. */
    aux_path = getenv ("L8_AUX_DIR");
    if (aux_path == NULL)
    {
        aux_path = ".";
        sprintf (errmsg, "L8_AUX_DIR environment variable isn't defined. It is "
            "assumed the auxiliary products will be available from the local "
            "directory.");
        error_handler (false, FUNC_NAME, errmsg);
    }

    /* Grab the year of the auxiliary input file to be used for the correct
       location of the auxiliary file in the auxliary directory */
    strncpy (aux_year, &aux_infile[5], 4);
    aux_year[4] = '\0';

    /* Set up the look-up table files and make sure they exist */
    sprintf (anglehdf, "%s/LDCMLUT/ANGLE_NEW.hdf", aux_path);
    sprintf (intrefnm, "%s/LDCMLUT/RES_LUT_V3.0-URBANCLEAN-V2.0.hdf", aux_path);
    sprintf (transmnm, "%s/LDCMLUT/TRANS_LUT_V3.0-URBANCLEAN-V2.0.ASCII",
        aux_path);
    sprintf (spheranm, "%s/LDCMLUT/AERO_LUT_V3.0-URBANCLEAN-V2.0.ASCII",
        aux_path);
    sprintf (cmgdemnm, "%s/CMGDEM.hdf", aux_path);
    sprintf (rationm, "%s/ratiomapndwiexp.hdf", aux_path);
    sprintf (auxnm, "%s/LADS/%s/%s", aux_path, aux_year, aux_infile);

    if (stat (anglehdf, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find anglehdf data file: %s\n  Check "
            "L8_AUX_DIR environment variable.", anglehdf);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (intrefnm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find intrefnm data file: %s\n  Check "
            "L8_AUX_DIR environment variable.", intrefnm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (transmnm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find transmnm data file: %s\n  Check "
            "L8_AUX_DIR environment variable.", transmnm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (spheranm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find spheranm data file: %s\n  Check "
            "L8_AUX_DIR environment variable.", spheranm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (cmgdemnm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find cmgdemnm data file: %s\n  Check "
            "L8_AUX_DIR environment variable.", cmgdemnm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (rationm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find rationm data file: %s\n  Check "
            "L8_AUX_DIR environment variable.", rationm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (auxnm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find auxnm data file: %s\n  Check "
            "L8_AUX_DIR environment variable.", auxnm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Initialization for look up tables */
    if (verbose)
        printf ("Initializing the look-up tables ...\n");
    xtv = 0.0;
    xfi = 0.0;
    xtsmin = 0;
    xtsstep = 4.0;
    xtvmin = 2.84090;
    xtvstep = 6.52107 - 2.84090;
    retval = readluts (tsmax, tsmin, ttv, tts, nbfi, nbfic, indts, rolutt,
        transt, sphalbt, normext, xtsstep, xtsmin, anglehdf, intrefnm, transmnm,
        spheranm);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Reading the LUTs");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    printf ("The LUTs for urban clean case v2.0 have been read.  We can now "
        "perform atmospheric correction.\n");

    /* Read the auxiliary data files used as input to the reflectance
       calculations */
    if (verbose)
        printf ("Reading the auxiliary files ...\n");
    retval = read_auxiliary_files (anglehdf, intrefnm, transmnm, spheranm,
        cmgdemnm, rationm, auxnm, dem, andwi, sndwi, ratiob1, ratiob2, ratiob7,
        intratiob1, intratiob2, intratiob7, slpratiob1, slpratiob2, slpratiob7,
        wv, oz);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Reading the auxiliary files");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Initialize the geolocation space applications */
    if (!get_geoloc_info (&xml_metadata, &space_def))
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

    /* Getting parameters for atmospheric correction */
    /* Update to get the parameter of the scene center */
    raot550nm = 0.12;
    pres = 1013.0;
    uoz = 0.30;
    uwv = 0.5;
    xtv = 0.0;
    xfi = 0.0;

    /* Read the QA band */
    if (get_input_qa_lines (input, 0, 0, nlines, qaband) != SUCCESS)
    {
        sprintf (errmsg, "Reading QA band");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Use scene center (and center of the pixel) to compute atmospheric
       parameters */
    img.l = (int) (nlines * 0.5) - 0.5;
    img.s = (int) (nsamps * 0.5) + 0.5;
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

    /* Use that lat/long to determine the line/sample in the
       CMG-related lookup tables, using the center of the UL pixel */
    ycmg = (89.975 - center_lat) * 20.0;    /* vs / 0.05 */
    xcmg = (179.975 + center_lon) * 20.0;   /* vs / 0.05 */
    lcmg = (int) (ycmg + 0.5);
    scmg = (int) (xcmg + 0.5);
    if ((lcmg < 0 || lcmg >= CMG_NBLAT) || (scmg < 0 || scmg >= CMG_NBLON))
    {
        sprintf (errmsg, "Invalid line/sample combination for the CMG-related "
            "lookup tables - line %d, sample %d (0-based).  CMG-based tables "
            "are %d lines x %d samples.", lcmg, scmg, CMG_NBLAT, CMG_NBLON);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (wv[lcmg][scmg] != 0)
        uwv = wv[lcmg][scmg] / 200.0;
    else
        uwv = 0.5;

    if (oz[lcmg][scmg] != 0)
        uoz = oz[lcmg][scmg] / 400.0;
    else
        uoz = 0.3;

    if (dem[lcmg][scmg] != -9999)
        pres = 1013.0 * exp (-dem[lcmg][scmg] * ONE_DIV_8500);
    else
        pres = 1013.0;

    raot550nm = 0.05;

    /* Loop through all the bands (except the QA band) and compute the TOA
       reflectance and at-sensor brightness temp */
    printf ("Calculating TOA reflectance and at-sensor brightness "
            "temperatures. Band ");
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
                exit (ERROR);
            }

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
                        sband[sband_ib][i] = (int) rotoa;
                }
                else
                    sband[sband_ib][i] = FILL_VALUE;
            }
        }  /* end if band <= band 9 */

        /* Read the current band and calibrate thermal bands */
        else if (ib == DN_BAND10)
        {
            if (get_input_th_lines (input, 0, 0, nlines, uband) != SUCCESS)
            {
                sprintf (errmsg, "Reading band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            /* Compute brightness temp for band 10.  Make sure it falls
               within the min/max range for the thermal bands. */
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
                        sband[SR_BAND10][i] = (int) (tmpf + 0.5);
                }
                else
                    sband[SR_BAND10][i] = FILL_VALUE;
            }
        }  /* end if band 10 */

        else if (ib == DN_BAND11)
        {
            if (get_input_th_lines (input, 1, 0, nlines, uband) != SUCCESS)
            {
                sprintf (errmsg, "Reading band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            /* Compute brightness temp for band 11.  Make sure it falls
               within the min/max range for the thermal bands. */
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
                        sband[SR_BAND11][i] = (int) (tmpf + 0.5);
                }
                else
                    sband[SR_BAND11][i] = FILL_VALUE;
            }
        }  /* end if band 11 */
    }  /* end for ib */
    printf ("\n");

    /* The input data has been read and calibrated. The memory can be freed. */
    free (uband);
    uband = NULL;

    /* Open the TOA output file, and set up the bands according to whether
       the TOA reflectance bands will be written. */
    toa_output = open_output (&xml_metadata, input, true /*toa*/);
    if (toa_output == NULL)
    {   /* error message already printed */
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    printf ("Writing TOA reflectance corrected data to the output files ...\n");

    /* If we are writing the TOA data, do so now for bands 1-7.  This will
       occur if the user specified TOA to be written or if the surface
       reflectance processing will not be completed. */
    if (write_toa || !process_sr)
    {
        for (ib = SR_BAND1; ib <= SR_BAND7; ib++)
        {
            printf ("  Band %d: %s\n", ib+1,
                toa_output->metadata.band[ib].file_name);
            if (put_output_lines (toa_output, sband[ib], ib, 0, nlines,
                sizeof (int16)) != SUCCESS)
            {
                sprintf (errmsg, "Writing output TOA data for band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            /* Create the ENVI header file this band */
            if (create_envi_struct (&toa_output->metadata.band[ib],
                &xml_metadata.global, &envi_hdr) != SUCCESS)
            {
                sprintf (errmsg, "Creating ENVI header structure.");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
      
            /* Write the ENVI header */
            strcpy (envi_file, toa_output->metadata.band[ib].file_name);
            cptr = strchr (envi_file, '.');
            strcpy (cptr, ".hdr");
            if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
            {
                sprintf (errmsg, "Writing ENVI header file.");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
        }

        /* Append the TOA reflectance bands, bands 1-7, to the XML file */
        if (append_metadata (7, toa_output->metadata.band, xml_infile) !=
            SUCCESS)
        {
            sprintf (errmsg, "Appending TOA reflectance bands to XML file.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Write bands 9-11 (cirrus and thermals), which don't get any further
       processing. */
    for (ib = SR_BAND9; ib <= SR_BAND11; ib++)
    {
        printf ("  Band %d: %s\n", ib+2,
            toa_output->metadata.band[ib].file_name);
        if (put_output_lines (toa_output, sband[ib], ib, 0, nlines,
            sizeof (int16)) != SUCCESS)
        {
            sprintf (errmsg, "Writing output TOA data for band %d", ib+2);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Create the ENVI header file this band */
        if (create_envi_struct (&toa_output->metadata.band[ib],
            &xml_metadata.global, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Creating ENVI header structure.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
  
        /* Write the ENVI header */
        strcpy (envi_file, toa_output->metadata.band[ib].file_name);
        cptr = strchr (envi_file, '.');
        strcpy (cptr, ".hdr");
        if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Writing ENVI header file.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Append the TOA cirrus/thermal band to the XML file */
        if (append_metadata (1, &toa_output->metadata.band[ib], xml_infile) !=
            SUCCESS)
        {
            sprintf (errmsg, "Appending TOA cirrus/thermal band to XML file.");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Close the output TOA products, cleanup bands, and free the memory */
    close_output (toa_output, true /*toa products*/);
    if (process_sr && !write_toa)
    {
        /* Remove the TOA bands 1-7 that were created by the open routine,
           since they aren't actually used */
        for (ib = SR_BAND1; ib <= SR_BAND7; ib++)
            unlink (toa_output->metadata.band[ib].file_name);
    }
    free_output (toa_output);

    /* Only continue with the surface reflectance corrections if SR processing
       has been requested and is possible due to the solar zenith angle */
    if (process_sr)
    {
        /* Loop through all the reflectance bands and perform atmospheric
           corrections */
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
            retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, ib, pres, tpres,
                aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
                sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv,
                uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                wvtransb, oztransa, rotoa, &roslamb, &tgo, &roatm, &ttatmg,
                &satm, &xrorayp, &next);
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
            for (i = 0; i < nlines*nsamps; i++)
            {
                /* If this pixel is not fill.  Otherwise fill pixels have
                   already been marked in the TOA calculations. */
                if (qaband[i] != 1)
                {
                    /* Store the TOA reflectance values, unscaled, for later
                       use before completing atmospheric corrections */
                    rotoa = sband[ib][i] * SCALE_FACTOR;
                    if (ib == DN_BAND1)
                        aerob1[i] = sband[ib][i];
                    else if (ib == DN_BAND2)
                        aerob2[i] = sband[ib][i];
                    else if (ib == DN_BAND4)
                        aerob4[i] = sband[ib][i];
                    else if (ib == DN_BAND5)
                        aerob5[i] = sband[ib][i];
                    else if (ib == DN_BAND7)
                        aerob7[i] = sband[ib][i];

                    /* Apply the atmospheric corrections, and store the scaled
                       value for later corrections */
                    roslamb = rotoa / tgo;
                    roslamb = roslamb - roatm;
                    roslamb = roslamb / ttatmg;
                    roslamb = roslamb / (1.0 + satm * roslamb);
                    sband[ib][i] = (int) (roslamb * MULT_FACTOR);
                }
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
        tmp_percent = 0;
        for (i = 0; i < nlines; i++)
        {
            /* update status? */
            if (100 * i / nlines > tmp_percent)
            {
                tmp_percent = 100 * i / nlines;
                if (tmp_percent % 10 == 0)
                {
                    printf ("%d%% ", tmp_percent);
                    fflush (stdout);
                }
            }

            curr_pix = i * nsamps;
            for (j = 0; j < nsamps; j++, curr_pix++)
            {
                /* If this pixel is fill, then don't process */
                if (qaband[curr_pix] == 1)
                    continue;

                /* Get the lat/long for the current pixel, for the center of
                   the pixel */
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
                   pixel */
                ycmg = (89.975 - lat) * 20.0;   /* vs / 0.05 */
                xcmg = (179.975 + lon) * 20.0;  /* vs / 0.05 */
                lcmg = (int) (ycmg);
                scmg = (int) (xcmg);
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

                u = (ycmg - lcmg);
                v = (xcmg - scmg);
                twvi[curr_pix] = wv[lcmg][scmg] * (1.0 - u) * (1.0 - v) +
                                 wv[lcmg][scmg+1] * (1.0 - u) * v +
                                 wv[lcmg+1][scmg] * u * (1.0 - v) +
                                 wv[lcmg+1][scmg+1] * u * v;
                twvi[curr_pix] = twvi[curr_pix] * 0.01;   /* vs / 100 */

                uoz11 = oz[lcmg][scmg];
                if (uoz11 == 0)
                    uoz11 = 120;

                uoz12 = oz[lcmg][scmg+1];
                if (uoz12 == 0)
                    uoz12 = 120;

                uoz21 = oz[lcmg+1][scmg];
                if (uoz21 == 0)
                    uoz21 = 120;

                uoz22 = oz[lcmg+1][scmg+1];
                if (uoz22 == 0)
                    uoz22 = 120;

                tozi[curr_pix] = uoz11 * (1.0 - u) * (1.0 - v) +
                                 uoz12 * (1.0 - u) * v +
                                 uoz21 * u * (1.0 - v) +
                                 uoz22 * u * v;
                tozi[curr_pix] = tozi[curr_pix] * 0.0025;   /* vs / 400 */

                if (dem[lcmg][scmg] != -9999)
                    pres11 = 1013.0 * exp (-dem[lcmg][scmg] * ONE_DIV_8500);
                else
                {
                    pres11 = 1013.0;
                    cloud[curr_pix] = 128;    /* set water bit */
                    tresi[curr_pix] = -1.0;
                }

                if (dem[lcmg][scmg+1] != -9999)
                    pres12 = 1013.0 * exp (-dem[lcmg][scmg+1] * ONE_DIV_8500);
                else
                    pres12 = 1013.0;

                if (dem[lcmg+1][scmg] != -9999)
                    pres21 = 1013.0 * exp (-dem[lcmg+1][scmg] * ONE_DIV_8500);
                else
                    pres21 = 1013.0;

                if (dem[lcmg+1][scmg+1] != -9999)
                    pres22 = 1013.0 * exp (-dem[lcmg+1][scmg+1] * ONE_DIV_8500);
                else
                    pres22 = 1013.0;

                tp[curr_pix] = pres11 * (1.0 - u) * (1.0 - v) +
                               pres12 * (1.0 - u) * v +
                               pres21 * u * (1.0 - v) +
                               pres22 * u * v;

                /* Inverting aerosols */
                /* Filter cirrus pixels */
                if (sband[SR_BAND9][curr_pix] >
                    (100.0 / (tp[curr_pix] * ONE_DIV_1013)))
                {  /* Set cirrus bit */
                    cloud[curr_pix]++;
                }
                else
                {  /* Inverting aerosol */
                    if (ratiob1[lcmg][scmg] == 0)
                    {
                        /* Average the valid ratio around the location */
                        erelc[DN_BAND1] = 0.4817;
                        erelc[DN_BAND2] = erelc[DN_BAND1] / 0.844239;
                        erelc[DN_BAND4] = 1.0;
                        erelc[DN_BAND7] = 1.79;
                    }
                    else
                    {
                        /* Use the NDWI to calculate the band ratio */
                        xndwi = ((double) sband[SR_BAND5][curr_pix] -
                                 (double) (sband[SR_BAND7][curr_pix] * 0.5)) /
                                ((double) sband[SR_BAND5][curr_pix] +
                                 (double) (sband[SR_BAND7][curr_pix] * 0.5));

                        th1 = (andwi[lcmg][scmg] + 2.0 * sndwi[lcmg][scmg]) *
                            0.001;
                        th2 = (andwi[lcmg][scmg] - 2.0 * sndwi[lcmg][scmg]) *
                            0.001;
                        if (xndwi > th1)
                            xndwi = th1;
                        if (xndwi < th2)
                            xndwi = th2;

                        erelc[DN_BAND1] = (xndwi * slpratiob1[lcmg][scmg] +
                            intratiob1[lcmg][scmg]) * 0.001;
                        erelc[DN_BAND2] = (xndwi * slpratiob2[lcmg][scmg] +
                            intratiob2[lcmg][scmg]) * 0.001;
                        erelc[DN_BAND4] = 1.0;
                        erelc[DN_BAND7] = (xndwi * slpratiob7[lcmg][scmg] +
                            intratiob7[lcmg][scmg]) * 0.001;
                    }

                    troatm[DN_BAND1] = aerob1[curr_pix] * SCALE_FACTOR;
                    troatm[DN_BAND2] = aerob2[curr_pix] * SCALE_FACTOR;
                    troatm[DN_BAND4] = aerob4[curr_pix] * SCALE_FACTOR;
                    troatm[DN_BAND7] = aerob7[curr_pix] * SCALE_FACTOR;

                    /* If this is water ... */
                    if (btest (cloud[curr_pix], WAT_QA))
                    {
                        /* Check the NDVI */
                        fndvi = ((double) sband[SR_BAND5][curr_pix] -
                                 (double) sband[SR_BAND4][curr_pix]) /
                                ((double) sband[SR_BAND5][curr_pix] +
                                 (double) sband[SR_BAND4][curr_pix]);
                        if (fndvi < 0.1)
                        {  /* skip the rest of the processing */
                            taero[curr_pix] = 0.0;
                            tresi[curr_pix] = -0.01;
                            continue;
                        }
                    }
           
                    iband1 = DN_BAND4;
                    iband3 = DN_BAND1;
                    retval = subaeroret (iband1, iband3, xts, xtv, xfi, pres,
                        uoz, uwv, erelc, troatm, tpres, aot550nm, rolutt,
                        transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                        normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv,
                        tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                        wvtransb, oztransa, &raot, &residual, &next);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing atmospheric correction.");
                        error_handler (true, FUNC_NAME, errmsg);
                        exit (ERROR);
                    }
                    corf = raot / xmus;

                    if (residual < (0.015 + 0.005 * corf))
                    {  /* test if band 5 makes sense */
                        iband = DN_BAND5;
                        rotoa = aerob5[curr_pix] * SCALE_FACTOR;
                        raot550nm = raot;
                        retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband,
                            pres, tpres, aot550nm, rolutt, transt, xtsstep,
                            xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                            tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv,
                            tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                            wvtransb, oztransa, rotoa, &roslamb, &tgo, &roatm,
                            &ttatmg, &satm, &xrorayp, &next);
                        if (retval != SUCCESS)
                        {
                            sprintf (errmsg, "Performing lambertian "
                                "atmospheric correction type 2.");
                            error_handler (true, FUNC_NAME, errmsg);
                            exit (ERROR);
                        }
                        ros5 = roslamb;

                        iband = DN_BAND4;
                        rotoa = aerob4[curr_pix] * SCALE_FACTOR;
                        raot550nm = raot;
                        retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband,
                            pres, tpres, aot550nm, rolutt, transt, xtsstep,
                            xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                            tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv,
                            tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                            wvtransb, oztransa, rotoa, &roslamb, &tgo, &roatm,
                            &ttatmg, &satm, &xrorayp, &next);
                        if (retval != SUCCESS)
                        {
                            sprintf (errmsg, "Performing lambertian "
                                "atmospheric correction type 2.");
                            error_handler (true, FUNC_NAME, errmsg);
                            exit (ERROR);
                        }
                        ros4 = roslamb;

                        if ((ros5 > 0.1) && ((ros5 - ros4) / (ros5 + ros4) > 0))
                        {
                            taero[curr_pix] = raot;
                            tresi[curr_pix] = residual;
                        }
                        else
                        {
                            taero[curr_pix] = 0.0;
                            tresi[curr_pix] = -0.01;
                        }
                    }
                    else
                    {
                        taero[curr_pix] = 0.0;
                        tresi[curr_pix] = -0.01;
                    }
                }  /* end if cirrus */
            }  /* end for i */
        }  /* end for j */

        /* update status */
        printf ("100%%\n");
        fflush (stdout);

        /* Done with the ratiob* arrays */
        for (i = 0; i < RATIO_NBLAT; i++)
        {
            free (andwi[i]);
            free (sndwi[i]);
            free (ratiob1[i]);
            free (ratiob2[i]);
            free (ratiob7[i]);
            free (intratiob1[i]);
            free (intratiob2[i]);
            free (intratiob7[i]);
            free (slpratiob1[i]);
            free (slpratiob2[i]);
            free (slpratiob7[i]);
        }
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

        /* Done with the aerob* arrays */
        free (aerob1);  aerob1 = NULL;
        free (aerob2);  aerob2 = NULL;
        free (aerob4);  aerob4 = NULL;
        free (aerob5);  aerob5 = NULL;
        free (aerob7);  aerob7 = NULL;

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
                nbval++;
                mall += sband[SR_BAND10][i] * SCALE_FACTOR_TH;
                if ((!btest (cloud[i], CIR_QA)) &&
                    (sband[SR_BAND5][i] > 300))
                {
                    anom = sband[SR_BAND2][i] - sband[SR_BAND4][i] * 0.5;
                    if (anom < 300)
                    {
                        nbclear++;
                        mclear += sband[SR_BAND10][i] * SCALE_FACTOR_TH;
                    }
                }
            }
        }  /* end for i */

        if (nbclear > 0)
            mclear = mclear / nbclear;
        else
            mclear = 275.0;

        if (nbval > 0)
            mall = mall / nbval;

        printf ("Average clear temperature %%clear %f %f %f %ld\n", mclear,
            nbclear * 100.0 / (nlines * nsamps), mall, nbval);

        /* Determine the cloud mask */
        for (i = 0; i < nlines*nsamps; i++)
        {
            if (tresi[i] < 0.0)
            {
                if (((sband[SR_BAND2][i] - sband[SR_BAND4][i] * 0.5) > 500) &&
                    ((sband[SR_BAND10][i] * SCALE_FACTOR_TH) < (mclear - 2.0)))
                {  /* Snow or cloud for now */
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
                if (btest (cloud[curr_pix], CLD_QA) ||
                    btest (cloud[curr_pix], CIR_QA))
                {
                    /* Check the 5x5 window around the current pixel */
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

                            if (!btest (cloud[win_pix], CLD_QA) &&
                                !btest (cloud[win_pix], CIR_QA) &&
                                !btest (cloud[win_pix], CLDA_QA))
                            {  /* Set the adjacent cloud bit */
                                cloud[win_pix] += 4;
                            }
                        }  /* for l */
                    }  /* for k */
                }  /* if btest */
            }  /* for j */
        }  /* for i */

#ifdef NOT_USED
        /* Compute adjustment to true North */
        /* Use scene center */
        img.l = (int) (nlines * 0.5);
        img.s = (int) (nsamps * 0.5);
        img.is_fill = false;
        row = img.l;
        col = img.s;
        printf ("Scene center line, sample: %d, %d\n", row, col);
        if (!from_space (space, &img, &geo))
        {
            sprintf (errmsg, "Mapping scene center to geolocation coords");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
        center_lat = geo.lat * RAD2DEG;
        center_lon = geo.lon * RAD2DEG;
        printf ("Scene center lat/long: %f, %f\n", center_lat, center_lon);

        /* Move 100 pixels to the north */
        img.l -= 100;
        if (!from_space (space, &img, &geo))
        {
            sprintf (errmsg, "Mapping 100 lines north of scene center to "
                "geolocation coords");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
        lat = geo.lat * RAD2DEG;
        lon = geo.lon * RAD2DEG;
        printf ("100 lines north of scene center lat/long: %f, %f\n", lat, lon);

        /* Use the longitude from the scene center and the latitude from the
           point 100 lines north to compute the line, sample */
        geo.lon = center_lon * DEG2RAD;
        geo.lat = lat * DEG2RAD;
        if (!to_space (space, &geo, &img))
        {
            sprintf (errmsg, "Mapping geolocation coords to line/sample space");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
        rowp = (int) img.l;
        colp = (int) img.s;
        printf ("Line, sample true north adj: %d, %d\n", rowp, colp);
        dy = row - rowp;
        dx = colp - col;
        ang = atan (dx / dy) * RAD2DEG;
        printf ("Adjustment to true North: %f\n", ang);
#endif

        /* Compute the cloud shadow */
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

                        win_pix = k * nsamps + l;
                        if ((sband[SR_BAND6][win_pix] < 800) &&
                            ((sband[SR_BAND3][win_pix] -
                              sband[SR_BAND4][win_pix]) < 100))
                        {
                            if (btest (cloud[win_pix], CLD_QA) ||
                                btest (cloud[win_pix], CIR_QA) ||
                                btest (cloud[win_pix], CLDS_QA))
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
                    /* Check the 6x6 window around the current pixel */
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

                            if (btest (cloud[win_pix], CLD_QA) ||
                                btest (cloud[win_pix], CLDS_QA))
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
                cloud[i] += 8;
                cloud[i] -= 16;
            }
        }  /* end for i */

        /* Aerosol interpolation */
        printf ("Performing aerosol interpolation ...\n");
        hole = 1;
        step = 10;
        while ((hole != 0) && (step < 1000))
        {
            hole = 0;
            for (i = 0; i < nlines; i += step)
            {
                for (j = 0; j < nsamps; j += step)
                {
                    nbaot = 0;
                    aaot = 0.0;
                    sresi = 0.0;

                    /* Check the window around the current pixel */
                    for (k = i; k <= i+step-1; k++)
                    {
                        /* Make sure the line is valid */
                        if (k < 0 || k >= nlines)
                            continue;

                        win_pix = k * nsamps + j;
                        for (l = j; l <= j+step-1; l++, win_pix++)
                        {
                            /* Make sure the sample is valid */
                            if (l < 0 || l >= nsamps)
                                continue;

                            if ((tresi[win_pix] > 0) && (cloud[win_pix] == 0))
                            {
                                nbaot++;
                                aaot += taero[win_pix] / tresi[win_pix];
                                sresi += 1.0 / tresi[win_pix];
                            }
                        }
                    }

                    /* If pixels were found */
                    if (nbaot != 0)
                    {
                        aaot /= sresi;

                        /* Check the window around the current pixel */
                        for (k = i; k <= i+step-1; k++)
                        {
                            /* Make sure the line is valid */
                            if (k < 0 || k >= nlines)
                                continue;

                            win_pix = k * nsamps + j;
                            for (l = j; l <= j+step-1; l++, win_pix++)
                            {
                                /* Make sure the sample is valid */
                                if (l < 0 || l >= nsamps)
                                    continue;

                                if ((tresi[win_pix] < 0) &&
                                    (!btest (cloud[win_pix], CIR_QA)) &&
                                    (!btest (cloud[win_pix], CLD_QA)) &&
                                    (!btest (cloud[win_pix], WAT_QA)))
                                {
                                    taero[win_pix] = aaot;
                                    tresi[win_pix] = 1.0;
                                }
                            }  /* for l */
                        }  /* for k */
                    }
                    else
                    {  /* this is a hole */
                        hole++;
                    }
                }  /* end for j */
            }  /* end for i */

            /* Modify the step value */
            step *= 2;
        }  /* end while */

        /* Perform the atmospheric correction */
        printf ("Performing atmospheric correction ...\n");
        /* 0 .. DN_BAND7 is the same as 0 .. SR_BAND7 here, since the pan band
           isn't spanned */
        for (ib = 0; ib <= DN_BAND7; ib++)
        {
            printf ("  Band %d\n", ib+1);
            for (i = 0; i < nlines * nsamps; i++)
            {
                /* If this pixel is fill, then don't process. Otherwise the
                   fill pixels have already been marked in the TOA process. */
                if (qaband[i] != 1)
                {
                    if (tresi[i] > 0.0 &&
                        !btest (cloud[i], CIR_QA) &&
                        !btest (cloud[i], CLD_QA))
                    {
                        rsurf = sband[ib][i] * SCALE_FACTOR;
                        rotoa = (rsurf * bttatmg[ib] / (1.0 - bsatm[ib] * rsurf)
                            + broatm[ib]) * btgo[ib];
                        raot550nm = taero[i];
                        pres = tp[i];
                        uwv = twvi[i];
                        uoz = tozi[i];
                        retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, ib,
                            pres, tpres, aot550nm, rolutt, transt, xtsstep,
                            xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                            tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv,
                            tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                            wvtransb, oztransa, rotoa, &roslamb, &tgo, &roatm,
                            &ttatmg, &satm, &xrorayp, &next);
                        if (retval != SUCCESS)
                        {
                            sprintf (errmsg, "Performing lambertian "
                                "atmospheric correction type 2.");
                            error_handler (true, FUNC_NAME, errmsg);
                            exit (ERROR);
                        }

                        /* Handle the aerosol computation in the cloud mask if
                           this is the cirrus band */
                        if (ib == DN_BAND1)
                        {
                            if (roslamb < -0.005)
                            {
                                taero[i] = 0.05;
                                raot550nm = 0.05;
                                pres = tp[i];
                                uwv = twvi[i];
                                uoz = tozi[i];
                                retval = atmcorlamb2 (xts, xtv, xfi, raot550nm,
                                    ib, pres, tpres, aot550nm, rolutt, transt,
                                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                                    normext, tsmax, tsmin, nbfic, nbfi, tts,
                                    indts, ttv, uoz, uwv, tauray, ogtransa1,
                                    ogtransb0, ogtransb1, wvtransa, wvtransb,
                                    oztransa, rotoa, &roslamb, &tgo, &roatm,
                                    &ttatmg, &satm, &xrorayp, &next);
                                if (retval != SUCCESS)
                                {
                                    sprintf (errmsg, "Performing lambertian "
                                        "atmospheric correction type 2.");
                                    error_handler (true, FUNC_NAME, errmsg);
                                    exit (ERROR);
                                }
                            }
                            else
                            {  /* Set up aerosol QA bits */
                                if (fabs (rsurf - roslamb) <= 0.015)
                                {  /* Set the first aerosol bit */
                                    cloud[i] += 16;
                                }
                                else
                                {
                                    if (fabs (rsurf - roslamb) < 0.03)
                                    {  /* Set the second aerosol bit */
                                        cloud[i] += 32;
                                    }
                                    else
                                    {  /* Set both aerosol bits */
                                        cloud[i] += 48;
                                    }
                                }
                            }  /* end if/else roslamb */
                        }  /* end if ib */

                        /* Save the scaled surface reflectance value, but make
                           sure it falls within the defined valid range. */
                        roslamb = roslamb * MULT_FACTOR;  /* scale the value */
                        if (roslamb < MIN_VALID)
                            sband[ib][i] = MIN_VALID;
                        else if (roslamb > MAX_VALID)
                            sband[ib][i] = MAX_VALID;
                        else
                            sband[ib][i] = (int) roslamb;
                    }  /* end if */
                }  /* end if qaband */
            }  /* end for i */
        }  /* end for ib */

        /* Write the data to the output file */
        printf ("Writing surface reflectance corrected data to the output "
            "files ...\n");

        /* Open the output file */
        sr_output = open_output (&xml_metadata, input, false /*surf refl*/);
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
                &xml_metadata.global, &envi_hdr) != SUCCESS)
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

        /* Create the ENVI header for the cloud mask band */
        if (create_envi_struct (&sr_output->metadata.band[SR_CLOUD],
            &xml_metadata.global, &envi_hdr) != SUCCESS)
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

        /* Close the output surface reflectance products */
        close_output (sr_output, false /*sr products*/);
        free_output (sr_output);
    }  /* end if process_sr */
  
    /* Free the spatial mapping pointer */
    free (space);

    /* Free the metadata structure */
    free_metadata (&xml_metadata);

    /* Close the input product */
    printf ("Closing input/output and freeing pointers ...\n");
    close_input (input);
    free_input (input);

    /* Free the filename pointers */
    free (xml_infile);
    free (aux_infile);

    /* Free the data arrays */
    for (i = 0; i < DEM_NBLAT; i++)
        free (dem[i]);
    free (dem);

    for (i = 0; i < CMG_NBLAT; i++)
    {
        free (wv[i]);
        free (oz[i]);
    }
    free (wv);
    free (oz);

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

    /* Free memory for band data */
    free (qaband);
    free (twvi);
    free (tozi);
    free (tp);
    free (tresi);
    free (taero);
    free (cloud);

    for (i = 0; i < NBAND_TTL_OUT-1; i++)
        free (sband[i]);
    free (sband);

    /* Indicate successful completion of processing */
    printf ("Surface reflectance processing complete!\n");
    exit (SUCCESS);
}


/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date        Programmer       Reason
---------   ---------------  -------------------------------------
7/6/2014    Gail Schmidt     Original Development
7/31/2014   Gail Schmidt     Added flag to write the TOA and process option
                             for surface reflectance

NOTES:
******************************************************************************/
void usage ()
{
    printf ("l8_sr computes the surface reflectance values for the input "
            "Landsat 8 DN products.  Surface reflectance correction and/or "
            "top of atmosphere correction is applied and written for bands "
            "1-7.  Top of atmosphere and at-sensor corrections are applied "
            "and written for bands 9 (cirrus), 10 (thermal), and 11 "
            "(thermal).\n\n");
    printf ("usage: l8_sr "
            "--xml=input_xml_filename "
            "--aux=input_auxiliary_filename "
            "--process_sr=true:false --write_toa [--verbose]\n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -xml: name of the input XML file to be processed\n");
    printf ("    -aux: name of the input auxiliary file containing ozone "
            "and water vapor for the scene date.  The file is expected to "
            "live in the $L8_AUX_DIR/LADS directory or in the local "
            "directory.\n");

    printf ("\nwhere the following parameters are optional:\n");
    printf ("    -process_sr: the default is to process surface reflectance, "
            "however if this flag is set to false then only the TOA "
            "reflectance processing and brightness temperature will be "
            "done.\n");
    printf ("    -write_toa: the intermediate TOA reflectance products "
            "for bands 1-7 are written to the output file\n");
    printf ("    -verbose: should intermediate messages be printed? (default "
            "is false)\n");

    printf ("\nl8_sr --help will print the usage statement\n");
    printf ("\nExample: l8_sr --xml=LC80410272013181LGN00.xml "
            "--aux=L8ANC2013181.hdf_fused --verbose\n");
    printf ("   ==> Writes bands 9-11 as TOA reflectance and brightness "
            "temperature.  Writes bands 1-7 as surface reflectance.\n\n");

    printf ("\nExample: l8_sr --xml=LC80410272013181LGN00.xml "
            "--aux=L8ANC2013181.hdf_fused --write_toa --verbose\n");
    printf ("   ==> Writes bands 1-11 as TOA reflectance and brightness "
            "temperature.  Writes bands 1-7 as surface reflectance.\n");

    printf ("\nExample: l8_sr --xml=LC80410272013181LGN00.xml "
            "--aux=L8ANC2013181.hdf_fused --process_sr=false --verbose\n");
    printf ("   ==> Writes bands 1-11 as TOA reflectance and brightness "
            "temperature.  Surface reflectance corrections are not applied.\n");
}


/******************************************************************************
MODULE:  btest

PURPOSE:  Tests to see if bit n is set in the byte_val variable.

RETURN VALUE:
Type = bool
Value      Description
-----      -----------
false      bit n is not set in byte_val
true       bit n is set in byte_val

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
7/10/2014    Gail Schmidt     Original Development

NOTES:
******************************************************************************/
bool btest
(
    uint8 byte_val,   /* I: byte value to be tested with the bit n */
    byte n            /* I: bit number to be tested (0 is rightmost bit) */
)
{
    /* Take 2 ** n, then AND that result with the byte value */
    return (byte_val & (1 << n));
}

