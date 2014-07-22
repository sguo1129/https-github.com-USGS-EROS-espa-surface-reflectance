#include <sys/stat.h>
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

NOTES:
******************************************************************************/
int main (int argc, char *argv[])
{
    bool verbose;            /* verbose flag for printing messages */
    char FUNC_NAME[] = "main"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char envi_file[STR_SIZE];/* ENVI filename */
    char sds_name[STR_SIZE]; /* name of the SDS being read */
    char *aux_path = NULL;   /* path for Landsat auxiliary data */
    char *xml_infile = NULL; /* input XML filename */
    char *aux_infile = NULL; /* input auxiliary filename for water vapor
                                and ozone*/
    char *cptr = NULL;       /* pointer to the file extension */

    int retval;              /* return status */
    int ib;                  /* looping variable for input bands */
    int sband_ib;            /* looping variable for output bands */
    Input_t *input = NULL;       /* input structure for the Landsat product */
    Output_t *sr_output = NULL;  /* output structure and metadata for the SR
                                    product */
    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
    Espa_global_meta_t *gmeta = NULL;   /* pointer to global meta */
    Envi_header_t envi_hdr;   /* output ENVI header information */

    /* Vars for forward/inverse mapping space */
    Geoloc_t *space = NULL;       /* structure for geolocation information */
    Space_def_t space_def;        /* structure to define the space mapping */
    Img_coord_float_t img;        /* coordinate in line/sample space */
    Geo_coord_t geo;              /* coordinate in lat/long space */
    float center_lat, center_lon; /* lat/long for scene center */
    float lat, lon;               /* pixel lat, long location */

    struct stat statbuf;      /* buffer for the file stat function */
    int nit;                  /* number of iterations */
    uint16 **uband = NULL;    /* array of input image data for a current band */
    uint16 **qaband = NULL;   /* QA band for the input image */
    int16 **aerob1 = NULL;    /* atmospherically corrected band 1 data */
    int16 **aerob2 = NULL;    /* atmospherically corrected band 2 data */
    int16 **aerob4 = NULL;    /* atmospherically corrected band 4 data */
    int16 **aerob5 = NULL;    /* atmospherically corrected band 5 data */
    int16 **aerob7 = NULL;    /* atmospherically corrected band 7 data */
    int16 ***sband = NULL;    /* output surface reflectance and brightness
                                 temp bands */
    int16 **dem = NULL;       /* CMG DEM data array [DEM_NBLAT][DEM_NBLON] */
    int16 **ratiob1 = NULL;   /* mean band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob2 = NULL;   /* mean band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob7 = NULL;   /* mean band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    uint16 **wv = NULL;       /* water vapor values [CMG_NBLAT][CMG_NBLON] */
    uint8 **oz = NULL;        /* ozone values [CMG_NBLAT][CMG_NBLON] */
    uint8 **cloud = NULL;     /* signed value to represent clouds */
    float **twvi = NULL;      /* interpolated water vapor value */
    float **tozi = NULL;      /* interpolated ozone value */
    float **tp = NULL;        /* interpolated pressure value */
    float **tresi = NULL;     /* residuals for each pixel */
    float **taero = NULL;     /* aerosol values for each pixel */

    int i, j, k, l;      /* looping variables */
    int lcmg, scmg;      /* line/sample index for the CMG */
    int iband;           /* current band */
    float u, v;
    float xcmg, ycmg;    /* x/y location for CMG */
    float xts;           /* solar zenith angle (deg) */
    float xfs;           /* solar azimuth angle (deg) */
    float xtv;           /* observation zenith angle (deg) */
    float xfi;           /* azimuthal difference between the sun and
                            observation angle (deg) */
    float xtsstep;       /* solar zenith step value */
    float xtsmin;        /* minimum solar zenith value */
    float xtvstep;       /* observation step value */
    float xtvmin;        /* minimum observation value */

    /* Variables for reading the HDF files */
    int status;                 /* return status of the HDF function */
    int start[5];               /* starting point to read SDS data */
    int edges[5];               /* number of values to read in SDS data */
    int sd_id;                  /* file ID for the HDF file */
    int sds_id;                 /* ID for the current SDS */
    int sds_index;              /* index for the current SDS */

    /* The following arguments are all names of the LUTs */
    char tauraynm[STR_SIZE];    /* molecular optical thickness filename */
    char gscoefnm[STR_SIZE];    /* gaseous transmission coef filename */
    char anglehdf[STR_SIZE];    /* angle HDF filename */
    char intrefnm[STR_SIZE];    /* intrinsic reflectance filename */
    char transmnm[STR_SIZE];    /* transmission filename */
    char spheranm[STR_SIZE];    /* spherical albedo filename */
    char cmgdemnm[STR_SIZE];    /* climate modeling grid DEM filename */
    char rationm[STR_SIZE];     /* ratio averages filename */
    char auxnm[STR_SIZE];     /* auxiliary filename for ozone and water vapor*/
    char sbandname[16][STR_SIZE] =  /* "band" names for input files */
       {"ldcmb1", "ldcmb2", "ldcmb3", "ldcmb4", "ldcmb5", "ldcmb6", "ldcmb7",
        "ldcmb8", "", "", "", "", "", "", "", ""}; /* only 8 bands are used */

    /* Atmospheric correction variables */
    /* Look up table for atmospheric and geometric quantities */
    float tauray[16];           /* molecular optical thickness coeff */
    float oztransa[16];         /* ozone transmission coeff */
    float wvtransa[16];         /* water vapor transmission coeff */
    float wvtransb[16];         /* water vapor transmission coeff */
    float wvtransc[16];         /* water vapor transmission coeff */
    float ogtransa0[16];        /* other gases transmission coeff */
    float ogtransa1[16];        /* other gases transmission coeff */
    float ogtransb0[16];        /* other gases transmission coeff */
    float ogtransb1[16];        /* other gases transmission coeff */
    float ogtransc0[16];        /* other gases transmission coeff */
    float ogtransc1[16];        /* other gases transmission coeff */

    float ****rolutt = NULL;    /*** I: intrinsic reflectance table
                                        [16][7][22][8000] */
    float ****transt = NULL;    /*** I: transmission table [16][7][22][22] */
    float ***sphalbt = NULL;    /*** I: spherical albedo table [16][7][22] */
    float **tsmax = NULL;       /* [20][22] */
    float **tsmin = NULL;       /* [20][22] */
    float **nbfic = NULL;       /* [20][22] */
    float **nbfi = NULL;        /* [20][22] */
    float **ttv = NULL;         /* [20][22] */
    float tts[22];
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

    float pixsize;      /* pixel size for the reflectance files */
    int nlines, nsamps; /* number of lines and samples in the reflectance and
                           thermal bands */
    int row, col;       /* pixel row, column (line, sample) location */
    int colp, rowp;     /* row/col for true north adjustment */
    int uoz11, uoz21, uoz12, uoz22;  /* ozone at line,samp; line, samp+1;
                           line+1, samp; and line+1, samp+1 */
    float pres11, pres12, pres21, pres22;  /* pressure at line,samp;
                           line, samp+1; line+1, samp; and line+1, samp+1 */
    float erelc[16];    /* band ratio variable GAIL - 8?? or 16?? */
    float troatm[16];   /* atmospheric reflectance table */
    float btgo[8];      /* other gaseous transmittance for bands 1-8 */
    float broatm[8];    /* atmospheric reflectance for bands 1-8 */
    float bttatmg[8];   /* ttatmg for bands 1-8 */
    float bsatm[8];     /* spherical albedo for bands 1-8 */
    int iband1, iband3; /* band indices (zero-based) */
    float raot;
    float residual;     /* model residual */
    float rsurf;
    float xmus;         /* cosine of solar zenith */
    float corf;
    float tmpf;         /* temporary floating point value */
    float xcals = 3.3420E-04;
    float xcalo = 0.10000;
    float k1b10 = 774.89;            /* temperature constant for band 10 */
    float k1b11 = 480.89;            /* temperature constant for band 11 */
    float k2b10 = 1321.08;           /* temperature constant for band 10 */
    float k2b11 = 1201.14;           /* temperature constant for band 11 */
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
    float aaot;
    float sresi;       /* sum of 1 / residuals */
    int nbaot;
    int step;
    int hole;
    float ros4, ros5;    /* surface reflectance for band 4 and band 5 */
    int tmp_percent;     /* current percentage for printing status */

    printf ("Starting surface reflectance processing ...\n");

    /* Read the command-line arguments */
    retval = get_args (argc, argv, &xml_infile, &aux_infile, &verbose);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }

    /* Provide user information if verbose is turned on */
    if (verbose)
    {
        printf ("  XML input file: %s\n", xml_infile);
        printf ("  AUX input file: %s\n", aux_infile);
    }

    /* Validate the input metadata file */
    if (validate_xml_file (xml_infile, ESPA_SCHEMA) != SUCCESS)
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

    /* Get the path for the auxiliary products from the ANC_PATH environment
       variable.  If it isn't defined, then assume the products are in the
       local directory. */
    aux_path = getenv ("ANC_PATH");
    if (aux_path == NULL)
    {
        aux_path = ".";
        sprintf (errmsg, "ANC_PATH environment variable isn't defined. It is "
            "assumed the auxiliary products will be available from the local "
            "directory.");
        error_handler (false, FUNC_NAME, errmsg);
    }

    /* Allocate memory for all the climate modeling grid files */
    dem = calloc (DEM_NBLAT, sizeof (int16*));
    if (dem == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the DEM");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    for (i = 0; i < DEM_NBLAT; i++)
    {
        dem[i] = calloc (DEM_NBLON, sizeof (int16));
        if (dem[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the DEM");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    ratiob1 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (ratiob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob1");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    ratiob2 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (ratiob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob2");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    ratiob7 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (ratiob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob7");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    for (i = 0; i < RATIO_NBLAT; i++)
    {
        ratiob1[i] = calloc (RATIO_NBLON, sizeof (int16));
        if (ratiob1[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the ratiob1");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        ratiob2[i] = calloc (RATIO_NBLON, sizeof (int16));
        if (ratiob2[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the ratiob2");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        ratiob7[i] = calloc (RATIO_NBLON, sizeof (int16));
        if (ratiob7[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the ratiob7");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    wv = calloc (CMG_NBLAT, sizeof (int16*));
    if (wv == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the wv");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    oz = calloc (CMG_NBLAT, sizeof (uint8*));
    if (oz == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the oz");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    for (i = 0; i < CMG_NBLAT; i++)
    {
        wv[i] = calloc (CMG_NBLON, sizeof (int16));
        if (wv[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the wv");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        oz[i] = calloc (CMG_NBLON, sizeof (uint8));
        if (oz[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the oz");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* rolutt[16][7][22][8000] and transt[16][7][22][22] and
       sphalbt[16][7][22] */
    rolutt = calloc (16, sizeof (float***));
    if (rolutt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for rolutt");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    transt = calloc (16, sizeof (float***));
    if (transt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for transt");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    sphalbt = calloc (16, sizeof (float**));
    if (sphalbt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for sphalbt");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    for (i = 0; i < 16; i++)
    {
        rolutt[i] = calloc (7, sizeof (float**));
        if (rolutt[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for rolutt");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        transt[i] = calloc (7, sizeof (float**));
        if (transt[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for transt");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        sphalbt[i] = calloc (7, sizeof (float*));
        if (sphalbt[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for sphalbt");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        for (j = 0; j < 7; j++)
        {
            rolutt[i][j] = calloc (22, sizeof (float*));
            if (rolutt[i][j] == NULL)
            {
                sprintf (errmsg, "Error allocating memory for rolutt");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            transt[i][j] = calloc (22, sizeof (float*));
            if (transt[i][j] == NULL)
            {
                sprintf (errmsg, "Error allocating memory for transt");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            sphalbt[i][j] = calloc (22, sizeof (float));
            if (sphalbt[i][j] == NULL)
            {
                sprintf (errmsg, "Error allocating memory for sphalbt");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            for (k = 0; k < 22; k++)
            {
                rolutt[i][j][k] = calloc (8000, sizeof (float));
                if (rolutt[i][j][k] == NULL)
                {
                    sprintf (errmsg, "Error allocating memory for rolutt");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }

                transt[i][j][k] = calloc (22, sizeof (float));
                if (transt[i][j][k] == NULL)
                {
                    sprintf (errmsg, "Error allocating memory for transt");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
            }  /* for k */
        }  /* for j */
    }  /* for i */

    /* tsmax[20][22] and float tsmin[20][22] and float nbfic[20][22] and
       nbfi[20][22] and float ttv[20][22] */
    tsmax = calloc (20, sizeof (float*));
    if (tsmax == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tsmax");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    tsmin = calloc (20, sizeof (float*));
    if (tsmin == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tsmin");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    nbfic = calloc (20, sizeof (float*));
    if (nbfic == NULL)
    {
        sprintf (errmsg, "Error allocating memory for nbfic");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    nbfi = calloc (20, sizeof (float*));
    if (nbfi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for nbfi");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    ttv = calloc (20, sizeof (float*));
    if (ttv == NULL)
    {
        sprintf (errmsg, "Error allocating memory for ttv");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    for (i = 0; i < 20; i++)
    {
        tsmax[i] = calloc (22, sizeof (float));
        if (tsmax[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for tsmax");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        tsmin[i] = calloc (22, sizeof (float));
        if (tsmin[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for tsmin");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        nbfic[i] = calloc (22, sizeof (float));
        if (nbfic[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for nbfic");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        nbfi[i] = calloc (22, sizeof (float));
        if (nbfi[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for nbfi");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        ttv[i] = calloc (22, sizeof (float));
        if (ttv[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for ttv");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Set up the look-up table files and make sure they exist */
    sprintf (tauraynm, "%s/LDCMLUT/tauray-ldcm.ASC", aux_path);
    sprintf (gscoefnm, "%s/LDCMLUT/gascoef-ldcm.ASC", aux_path);
    sprintf (anglehdf, "%s/LDCMLUT/ANGLE_NEW.hdf", aux_path);
    sprintf (intrefnm, "%s/LDCMLUT/RES_LUT_V3.0-URBANCLEAN-V2.0.hdf", aux_path);
    sprintf (transmnm, "%s/LDCMLUT/TRANS_LUT_V3.0-URBANCLEAN-V2.0.ASCII",
        aux_path);
    sprintf (spheranm, "%s/LDCMLUT/AERO_LUT_V3.0-URBANCLEAN-V2.0.ASCII",
        aux_path);
    sprintf (cmgdemnm, "%s/CMGDEM.hdf", aux_path);
    sprintf (rationm, "%s/newratio_averagesSD.hdf", aux_path);
    sprintf (auxnm, "%s/LANDSATANC/%s", aux_path, aux_infile);

    if (stat (tauraynm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find tauraynm data file: %s\n  Check "
            "ANC_PATH environment variable.", tauraynm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (gscoefnm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find gscoefnm data file: %s\n  Check "
            "ANC_PATH environment variable.", gscoefnm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (anglehdf, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find anglehdf data file: %s\n  Check "
            "ANC_PATH environment variable.", anglehdf);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (intrefnm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find intrefnm data file: %s\n  Check "
            "ANC_PATH environment variable.", intrefnm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (transmnm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find transmnm data file: %s\n  Check "
            "ANC_PATH environment variable.", transmnm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (spheranm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find spheranm data file: %s\n  Check "
            "ANC_PATH environment variable.", spheranm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (cmgdemnm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find cmgdemnm data file: %s\n  Check "
            "ANC_PATH environment variable.", cmgdemnm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (rationm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find rationm data file: %s\n  Check "
            "ANC_PATH environment variable.", rationm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    if (stat (auxnm, &statbuf) == -1)
    {
        sprintf (errmsg, "Could not find auxnm data file: %s\n  Check "
            "ANC_PATH environment variable.", auxnm);
        error_handler (false, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Initialization for look up tables */
    printf ("Initializing the look-up tables ...\n");
    xtv = 0.0;
    xfi = 0.0;
    xtsmin = 0;
    xtsstep = 4.0;
    xtvmin = 2.84090;
    xtvstep = 6.52107-2.84090;
    retval = readluts (tauray, oztransa, wvtransa, wvtransb, wvtransc,
        ogtransa0, ogtransa1, ogtransb0, ogtransb1, ogtransc0, ogtransc1,
        tsmax, tsmin, ttv, tts, nbfi, nbfic, indts, rolutt, transt, sphalbt,
        xtsstep, xtsmin, sbandname, tauraynm, gscoefnm, anglehdf, intrefnm,
        transmnm, spheranm);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Reading the LUTs");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    printf ("The LUTs for urban clean case v2.0 have been read.  We can now "
        "perform atmospheric correction.\n");

    /* Read the DEM */
    sd_id = SDstart (cmgdemnm, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", cmgdemnm);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "averaged elevation");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the DEM file %s", sds_name,
            cmgdemnm);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < DEM_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = DEM_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, dem[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Close the DEM file */
    status = SDend (sd_id);
    if (status != 0)
    {
        sprintf (errmsg, "Closing DEM file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the RATIO file */
    sd_id = SDstart (rationm, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", rationm);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "ratiob9 Mean");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, ratiob1[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "ratiob3 Mean");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, ratiob2[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "ratiob7 Mean");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, ratiob7[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Close the RATIO file */
    status = SDend (sd_id);
    if (status != 0)
    {
        sprintf (errmsg, "Closing RATIO file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read ozone and water vapor from the user-specified auxiliary file */
    sd_id = SDstart (auxnm, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", auxnm);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "Coarse Resolution Ozone");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the AUX file %s", sds_name,
            auxnm);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < CMG_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = CMG_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, oz[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "Coarse Resolution Water Vapor");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the AUX file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < CMG_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = CMG_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, wv[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Close the AUX file */
    status = SDend (sd_id);
    if (status != 0)
    {
        sprintf (errmsg, "Closing AUX file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Allocate space for band data */
    uband = calloc (nlines, sizeof (uint16*));
    if (uband == NULL)
    {
        sprintf (errmsg, "Error allocating memory for uband");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    qaband = calloc (nlines, sizeof (uint16*));
    if (qaband == NULL)
    {
        sprintf (errmsg, "Error allocating memory for qaband");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    aerob1 = calloc (nlines, sizeof (int16*));
    if (aerob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob1");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    aerob2 = calloc (nlines, sizeof (int16*));
    if (aerob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob2");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    aerob4 = calloc (nlines, sizeof (int16*));
    if (aerob4 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob4");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    aerob5 = calloc (nlines, sizeof (int16*));
    if (aerob5 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob5");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    aerob7 = calloc (nlines, sizeof (int16*));
    if (aerob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob7");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    twvi = calloc (nlines, sizeof (float*));
    if (twvi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for twvi");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    tozi = calloc (nlines, sizeof (float*));
    if (tozi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tozi");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    tp = calloc (nlines, sizeof (float*));
    if (tp == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tp");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    tresi = calloc (nlines, sizeof (float*));
    if (tresi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tresi");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    taero = calloc (nlines, sizeof (float*));
    if (taero == NULL)
    {
        sprintf (errmsg, "Error allocating memory for taero");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    cloud = calloc (nlines, sizeof (uint8*));
    if (cloud == NULL)
    {
        sprintf (errmsg, "Error allocating memory for cloud");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    for (i = 0; i < nlines; i++)
    {
        uband[i] = calloc (nsamps, sizeof (uint16));
        if (uband[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for uband");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        qaband[i] = calloc (nsamps, sizeof (uint16));
        if (qaband[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for qaband");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        aerob1[i] = calloc (nsamps, sizeof (int16));
        if (aerob1[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for aerob1");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        aerob2[i] = calloc (nsamps, sizeof (int16));
        if (aerob2[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for aerob2");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        aerob4[i] = calloc (nsamps, sizeof (int16));
        if (aerob4[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for aerob4");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        aerob5[i] = calloc (nsamps, sizeof (int16));
        if (aerob5[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for aerob5");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        aerob7[i] = calloc (nsamps, sizeof (int16));
        if (aerob7[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for aerob7");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        twvi[i] = calloc (nlines, sizeof (float));
        if (twvi[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for twvi");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        tozi[i] = calloc (nlines, sizeof (float));
        if (tozi[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for tozi");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        tp[i] = calloc (nlines, sizeof (float));
        if (tp[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for tp");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        tresi[i] = calloc (nlines, sizeof (float));
        if (tresi[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for tresi");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        taero[i] = calloc (nlines, sizeof (float));
        if (taero[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for taero");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        cloud[i] = calloc (nlines, sizeof (uint8));
        if (cloud[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for cloud");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    sband = calloc (NBAND_TTL_OUT, sizeof (int16**));
    if (sband == NULL)
    {
        sprintf (errmsg, "Error allocating memory for sband");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    for (i = 0; i < NBAND_TTL_OUT; i++)
    {
        sband[i] = calloc (nlines, sizeof (int16*));
        if (sband[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for sband");
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        for (j = 0; j < nlines; j++)
        {
            sband[i][j] = calloc (nsamps, sizeof (int16));
            if (sband[i][j] == NULL)
            {
                sprintf (errmsg, "Error allocating memory for sband");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
        }
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
    for (i = 0; i < nlines; i++)
    {
        if (get_input_qa_lines (input, 0, i, 1, qaband[i]) != SUCCESS)
        {
            sprintf (errmsg, "Reading line %d for QA band", i);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
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
    ycmg = (89.975 - center_lat) * 20.0;   /* vs / 0.05 */
    xcmg = (179.975 + center_lon) * 20.0;  /* vs / 0.05 */
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

    /* Loop through all the bands (except the QA band) and perform atmospheric
       corrections */
    printf ("Calibrating reflectance and thermal bands.\n");
    for (ib = 0; ib <= NBAND_TTL_MAX-1; ib++)
    {
        /* Don't process the pan band */
        if (ib == DN_BAND8)
            continue;

        /* Get the parameters for the atmospheric correction */
        if (ib < DN_BAND9)
        {
            /* rotoa is not defined for this call, which is ok, but the
               roslamb value is not valid upon output. Just set it to 0.0 to
               be consistent. */
            rotoa = 0.0;
            retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, ib, pres, tpres,
                aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
                sphalbt, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv,
                tauray, ogtransa0, ogtransa1, ogtransb0, ogtransb1, ogtransc0,
                ogtransc1, wvtransa, wvtransb, wvtransc, oztransa, rotoa, 
                &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp);
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
        }

        /* Read the current band and calibrate bands 1-9 (except pan) */
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

            for (i = 0; i < nlines; i++)
            {
                if (get_input_refl_lines (input, iband, i, 1, uband[i]) !=
                    SUCCESS)
                {
                    sprintf (errmsg, "Reading line %d for band %d", i, ib+1);
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }

                for (j = 0; j < nsamps; j++)
                {
                    rotoa = (uband[i][j] * 2.0000E-05) - 0.1;
                    sband[sband_ib][i][j] = (int) (rotoa * 10000.0 / xmus);
                }
            }
        }

        /* Perform atmospheric corrections for bands 1-7 */
        if (ib <= DN_BAND7)
        {
            for (i = 0; i < nlines; i++)
            {
                for (j = 0; j < nsamps; j++)
                {
                    /* If this pixel is not fill -- GAIL -- should we use btest with bit 0 here?? */
                    if (qaband[i][j] != 1)
                    {
                        rotoa = sband[sband_ib][i][j] * 0.0001;  /* div 10000 */
                        if (ib == DN_BAND1)
                            aerob1[i][j] = sband[sband_ib][i][j];
                        if (ib == DN_BAND2)
                            aerob2[i][j] = sband[sband_ib][i][j];
                        if (ib == DN_BAND4)
                            aerob4[i][j] = sband[sband_ib][i][j];
                        if (ib == DN_BAND5)
                            aerob5[i][j] = sband[sband_ib][i][j];
                        if (ib == DN_BAND7)
                            aerob7[i][j] = sband[sband_ib][i][j];

                        roslamb = rotoa / tgo;
 	                    roslamb = roslamb - roatm;
 	                    roslamb = roslamb / ttatmg;
 	                    roslamb = roslamb / (1.0 + satm * roslamb);
	                    sband[sband_ib][i][j] = (int) (roslamb * 10000.0);
                    }
                    else
                        sband[sband_ib][i][j] = FILL_VALUE;
                }  /* end for i */
            }  /* end for j */
        }  /* if ib */

        /* Read the current band and calibrate thermal bands */
        if (ib == DN_BAND10)
        {
            for (i = 0; i < nlines; i++)
            {
                if (get_input_th_lines (input, 0, i, 1, uband[i]) != SUCCESS)
                {
                    sprintf (errmsg, "Reading line %d for band %d", i, ib+1);
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }

                for (j = 0; j < nsamps; j++)
                {
                    /* If this pixel is not fill */
                    if (qaband[i][j] != 1)
                    {
                        tmpf = xcals * uband[i][j] + xcalo;
                        tmpf = k2b10 / log (1.0 + k1b10 / tmpf);
                        sband[SR_BAND10][i][j] = (int) (tmpf * 10.0);
                    }
                    else
                        sband[SR_BAND10][i][j] = FILL_VALUE;
                }  /* end for i */
            }  /* end for j */
        }  /* end if ib */

        if (ib == DN_BAND11)
        {
            for (i = 0; i < nlines; i++)
            {
                if (get_input_th_lines (input, 1, i, 1, uband[i]) != SUCCESS)
                {
                    sprintf (errmsg, "Reading line %d for band %d", i, ib+1);
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }

                for (j = 0; j < nsamps; j++)
                {
                    /* If this pixel is not fill */
                    if (qaband[i][j] != 1)
                    {
                        tmpf = xcals * uband[i][j] + xcalo;
                        tmpf = k2b11 / log (1.0 + k1b11 / tmpf);
                        sband[SR_BAND11][i][j] = (int) (tmpf * 10.0);
                    }
                    else
                        sband[SR_BAND11][i][j] = FILL_VALUE;
                }  /* end for i */
            }  /* end for j */
        }  /* end if ib */
    }  /* end for ib */

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

        for (j = 0; j < nsamps; j++)
        {
            /* If this pixel is fill, then don't process */
            if (qaband[i][j] == 1)
                continue;

            /* Get the lat/long for the current pixel, for the center of the
               pixel */
            img.l = i - 0.5;  /* GAIL */
            img.s = j + 0.5;
            img.is_fill = false;
            if (!from_space (space, &img, &geo))
            {
                sprintf (errmsg, "Mapping line/sample (%d, %d) to geolocation "
                    "coords", i, j);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
            lat = geo.lat * RAD2DEG;
            lon = geo.lon * RAD2DEG;

            /* Use that lat/long to determine the line/sample in the
               CMG-related lookup tables, using the center of the UL pixel */
            ycmg = (89.975 - lat) * 20.0;   /* vs / 0.05 */
            xcmg = (179.975 + lon) * 20.0;  /* vs / 0.05 */
            lcmg = (int) (ycmg);
            scmg = (int) (xcmg);
            if ((lcmg < 0 || lcmg >= CMG_NBLAT) ||
                (scmg < 0 || scmg >= CMG_NBLON))
            {
                sprintf (errmsg, "Invalid line/sample combination for the "
                    "CMG-related lookup tables - line %d, sample %d (0-based). "
                    "CMG-based tables are %d lines x %d samples.", lcmg, scmg,
                    CMG_NBLAT, CMG_NBLON);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            u = (ycmg - lcmg);
            v = (xcmg - scmg);
            twvi[i][j] = wv[lcmg][scmg] * (1.0 - u) * (1.0 - v) +
                         wv[lcmg][scmg+1] * (1.0 - u) * v +
                         wv[lcmg+1][scmg] * u * (1.0 - v) +
                         wv[lcmg+1][scmg+1] * u * v;
            twvi[i][j] = twvi[i][j] * 0.01;   /* vs / 100 */

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

            tozi[i][j] = uoz11 * (1.0 - u) * (1.0 - v) +
                         uoz12 * (1.0 - u) * v +
                         uoz21 * u * (1.0 - v) +
                         uoz22 * u * v;
            tozi[i][j] = tozi[i][j] * 0.0025;   /* vs / 400 */

            if (dem[lcmg][scmg] != -9999)
                pres11 = 1013.0 * exp (-dem[lcmg][scmg] * ONE_DIV_8500);
            else
            {
                pres11 = 1013.0;
                cloud[i][j] = 128;    /* set water bit */
                tresi[i][j] = -1.0;
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

            tp[i][j] = pres11 * (1.0 - u) * (1.0 - v) +
                       pres12 * (1.0 - u) * v +
                       pres21 * u * (1.0 - v) +
                       pres22 * u * v;

            /* Inverting aerosols */
            /* Filter cirrus pixels */
            if (sband[SR_BAND9][i][j] > (100.0 / (tp[i][j] * ONE_DIV_1013)))
            {  /* Set cirrus bit */
                cloud[i][j]++;
            }
            else
            {  /* Inverting aerosol */
                for (ib = 0; ib < 8; ib++)
                    erelc[ib] = -1.0;

                if (ratiob1[lcmg][scmg] == 0)
                {
                    erelc[DN_BAND4] = 1.0;
                    erelc[DN_BAND1] = 0.417;
                    erelc[DN_BAND2] = 0.476;
                    erelc[DN_BAND7] = 1.79;
                }
                else
                {
                    erelc[DN_BAND4] = 1.0;
                    erelc[DN_BAND1] = ratiob1[lcmg][scmg] * 0.001;/* vs /1000 */
                    erelc[DN_BAND2] = ratiob2[lcmg][scmg] * 0.001;/* vs /1000 */
                    erelc[DN_BAND7] = ratiob7[lcmg][scmg] * 0.001;/* vs /1000 */
                }

                troatm[0] = aerob1[i][j] * 0.0001;  /* vs / 10000 */
                troatm[1] = aerob2[i][j] * 0.0001;  /* vs / 10000 */
                troatm[3] = aerob4[i][j] * 0.0001;  /* vs / 10000 */
                troatm[6] = aerob7[i][j] * 0.0001;  /* vs / 10000 */

                /* If this is water ... */
                if (btest (cloud[i][j], WAT_QA))
                {
                    /* Check the NDVI */
                    if (((sband[SR_BAND5][i][j] - sband[SR_BAND4][i][j]) /
                         (sband[SR_BAND5][i][j] + sband[SR_BAND4][i][j])) < 0.1)
                    {  /* skip the rest of the processing */
                        taero[i][j] = 0.0;
                        tresi[i][j] = -0.01;
                        continue;
                    }
                }
       
                iband1 = DN_BAND4;
                iband3 = DN_BAND1;
                retval = subaeroret (iband1, iband3, xts, xtv, xfi, pres, uoz,
                    uwv, erelc, troatm, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, tsmax, tsmin,
                    nbfic, nbfi, tts, indts, ttv, tauray, ogtransa0, ogtransa1,
                    ogtransb0, ogtransb1, ogtransc0, ogtransc1, wvtransa,
                    wvtransb, wvtransc, oztransa, &raot, &residual, &nit);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing atmospheric correction.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
                corf = raot / xmus;

                if (residual < (0.01 + 0.005 * corf))
                {  /* test if band 5 makes sense */
                    iband = DN_BAND5;
                    rotoa = aerob5[i][j] * 0.0001;  /* vs / 10000 */
                    raot550nm = raot;
                    retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres,
                        tpres, aot550nm, rolutt, transt, xtsstep, xtsmin,
                        xtvstep, xtvmin, sphalbt, tsmax, tsmin, nbfic, nbfi,
                        tts, indts, ttv, uoz, uwv, tauray, ogtransa0,
                        ogtransa1, ogtransb0, ogtransb1, ogtransc0, ogtransc1,
                        wvtransa, wvtransb, wvtransc, oztransa, rotoa, &roslamb,
                        &tgo, &roatm, &ttatmg, &satm, &xrorayp);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing lambertian atmospheric "
                            "correction type 2.");
                        error_handler (true, FUNC_NAME, errmsg);
                        exit (ERROR);
                    }
                    ros5 = roslamb;

                    iband = DN_BAND4;
                    rotoa = aerob4[i][j] * 0.0001;  /* vs / 10000 */
                    raot550nm = raot;
                    retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres,
                        tpres, aot550nm, rolutt, transt, xtsstep, xtsmin,
                        xtvstep, xtvmin, sphalbt, tsmax, tsmin, nbfic, nbfi,
                        tts, indts, ttv, uoz, uwv, tauray, ogtransa0,
                        ogtransa1, ogtransb0, ogtransb1, ogtransc0, ogtransc1,
                        wvtransa, wvtransb, wvtransc, oztransa, rotoa, &roslamb,
                        &tgo, &roatm, &ttatmg, &satm, &xrorayp);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Performing lambertian atmospheric "
                            "correction type 2.");
                        error_handler (true, FUNC_NAME, errmsg);
                        exit (ERROR);
                    }
                    ros4 = roslamb;

                    if ((ros5 > 0.1) && ((ros5 - ros4) / (ros5 + ros4) > 0))
                    {
                        taero[i][j] = raot;
                        tresi[i][j] = residual;
                    }
                    else
                    {
                        taero[i][j] = 0.0;
                        tresi[i][j] = -0.01;
                    }
                }
                else
                {
                    taero[i][j] = 0.0;
                    tresi[i][j] = -0.01;
                }
            }  /* end if cirrus */
        }  /* end for i */
    }  /* end for j */

    /* update status */
    printf ("100%%\n");
    fflush (stdout);

    /* Refine the cloud mask */
    /* Compute the average temperature of the clear, non-water, non-filled
       pixels */
    printf ("Refining the cloud mask ...\n");
    nbval = 0;
    nbclear = 0;
    mclear = 0.0;
    mall = 0.0;
    for (i = 0; i < nlines; i++)
    {
        for (j = 0; j < nsamps; j++)
        {
            /* If this pixel is fill, then don't process */
            if (qaband[i][j] != 1)
            {
                nbval++;
                mall += sband[SR_BAND10][i][j] * 0.1;  /* vs / 10 */
                if ((!btest (cloud[i][j], CIR_QA)) &&
                    (sband[SR_BAND5][i][j] > 300))
                {
                    anom = sband[SR_BAND2][i][j] -
                           sband[SR_BAND4][i][j] * 0.5;  /* vs / 2 */
                    if (anom < 300)
                    {
                        nbclear++;
                        mclear += sband[SR_BAND10][i][j] * 0.1;  /* vs / 10 */
                    }
                }
            }
        }  /* end for i */
    }  /* end for j */

    if (nbclear > 0)
        mclear = mclear / nbclear;
    else
        mclear = 275.0;

    if (nbval > 0)
        mall = mall / nbval;

    printf ("Average clear temperature %%clear %f %f %f %ld\n", mclear,
        nbclear * 100.0 / (nlines * nsamps), mall, nbval);

    /* Determine the cloud mask */
    for (i = 0; i < nlines; i++)
    {
        for (j = 0; j < nsamps; j++)
        {
            if (tresi[i][j] < 0.0)
            {
                if (((sband[SR_BAND2][i][j] - sband[SR_BAND4][i][j] * 0.5) >
                    500) && ((sband[SR_BAND10][i][j] * 0.1) < (mclear - 2.0)))
                {  /* Snow or cloud for now */
                    cloud[i][j] += 2;
                }
            }
        }
    }

    /* Set up the adjacent to something bad (snow or cloud) bit */
    printf ("Setting up the adjacent to something bit ...\n");
    for (i = 0; i < nlines; i++)
    {
        for (j = 0; j < nsamps; j++)
        {
            if (btest (cloud[i][j], CLD_QA) || (btest (cloud[i][j], CIR_QA)))
            {
                /* Check the 5x5 window around the current pixel */
                for (k = i-5; k <= i+5; k++)
                {
                    for (l = j-5; l <= j+5; l++)
                    {
                        if ((k >= 0) && (k < nlines) &&
                            (l >= 0) && (l < nsamps))
                        {
                            if ((!btest (cloud[k][l], CLD_QA)) &&
                                (!btest (cloud[k][l], CIR_QA)) &&
                                (!btest (cloud[k][l], CLDA_QA)))
                            {  /* Set the adjacent cloud bit */
                                cloud[k][l] += 4;
                            }
                        }
                    }  /* for l */
                }  /* for k */
            }
        }  /* for j */
    }  /* for i */

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
        sprintf (errmsg, "Mapping geolocation coords to line, sample space");
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

    /* Compute the cloud shadow */
    printf ("Determining cloud shadow ...\n");
    facl = cosf (xfs * DEG2RAD) * tanf (xts * DEG2RAD) / pixsize;  /* lines */
    fack = sinf (xfs * DEG2RAD) * tanf (xts * DEG2RAD) / pixsize;  /* samps */
    for (i = 0; i < nlines; i++)
    {
        for (j = 0; j < nsamps; j++)
        {
            if (btest (cloud[i][j], CLD_QA) || btest (cloud[i][j], CIR_QA))
            {
                tcloud = sband[SR_BAND10][i][j] * 0.1;  /* vs / 10 */
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
                    k = j - fack * cldh;
                    l = i + facl * cldh;
                    if ((l >= 0) && (l < nlines) && (k >= 0) && (k < nsamps))
                    {
                        if ((sband[SR_BAND6][l][k] < 800) &&
                            ((sband[SR_BAND3][l][k] - sband[SR_BAND4][l][k]) <
                            100))
                        {
                            if (btest (cloud[l][k], CLD_QA) ||
                                btest (cloud[l][k], CIR_QA) ||
                                btest (cloud[l][k], CLDS_QA))
                            {
                                continue;
                            }
                            else
                            { /* store the value of band6 as well as the
                                 l and k value */
                                if (sband[SR_BAND6][l][k] < mband5)
                                {
                                     mband5 = sband[SR_BAND6][l][k];
                                     mband5k = k;
                                     mband5l = l;
                                }
                            }
                        }
                    }
                }  /* for icldh */

                /* Set the cloud shadow bit */
                if (mband5 < 9999)
                    cloud[mband5l][mband5k] += 8;
            }  /* end if btest */
        }  /* end for j */
    }  /* end for i */
	 
    /* Expand the cloud shadow using the residual */
    printf ("Expanding cloud shadow ...\n");
    for (i = 0; i < nlines; i++)
    {
        for (j = 0; j < nsamps; j++)
        {
            /* If this is a cloud shadow pixel */
            if (btest (cloud[i][j], CLDS_QA))
            {
                /* Check the 6x6 window around the current pixel */
                for (k = i-6; k <= i+6; k++)
                {
                    for (l = j-6; l <= j+6; l++)
                    {
                        if ((k >= 0) && (k < nlines) &&
                            (l >= 0) && (l < nsamps))
                        {
                            if (btest (cloud[k][l], CLD_QA) ||
                                btest (cloud[k][l], CLDS_QA))
                                continue;
                            else
                            {
                                if (btest (cloud[k][l], CLDT_QA))
                                    continue;
                                else
                                {
                                    /* Set the temporary bit */
                                    if (tresi[k][l] < 0)
                                        cloud[k][l] += 16;
                                }
                            }
                        }
                    }  /* end for l */
                }  /* end for k */
            }  /* end if btest */
        }  /* end for j */
    }  /* end for i */

    /* Update the cloud shadow */
    printf ("Updating cloud shadow ...\n");
    for (i = 0; i < nlines; i++)
    {
        for (j = 0; j < nsamps; j++)
        {
            /* If the temporary bit was set in the above loop */
            if (btest (cloud[i][j], CLDT_QA))
            {
                /* Remove the temporary bit and set the cloud shadow bit */
                cloud[i][j] += 8;
                cloud[i][j] -= 16;
            }
        }  /* end for j */
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
                    for (l = j; l <= j+step-1; l++)
                    {
                        if ((k >= 0) && (k < nlines) &&
                            (l >= 0) && (l < nsamps))
                        {
                            if ((tresi[k][l] > 0) && (cloud[k][l] == 0))
                            {
                                nbaot++;
                                aaot += taero[k][l] / tresi[k][l];
                                sresi += 1.0 / tresi[k][l];
                            }
                        }
                    }
                }

                /* If pixels were found */
                if (nbaot != 0)
                {
                    aaot = aaot / sresi;

                    /* Check the window around the current pixel */
                    for (k = i; k <= i+step-1; k++)
                    {
                        for (l = j; l <= j+step-1; l++)
                        {
                            if ((k >= 0) && (k < nlines) &&
                                (l >= 0) && (l < nsamps))
                            {
                                if ((tresi[k][l] < 0) &&
                                    (!btest (cloud[k][l], CIR_QA)) &&
                                    (!btest (cloud[k][l], CLD_QA)) &&
                                    (!btest (cloud[k][l], WAT_QA)))
                                {
                                    taero[k][l] = aaot;
                                    tresi[k][l] = 1.0;
                                }
                            }
                        }
                    }
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
        printf ("  Band %d: ", ib+1);
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

            for (j = 0; j < nsamps; j++)
            {
                /* If this pixel is fill, then don't process */
                if (qaband[i][j] != 1)
                {
                    if ((tresi[i][j] > 0.0) &&
                        (!btest (cloud[i][j], CIR_QA)) &&
                        (!btest (cloud[i][j], CLD_QA)))
                    {
                        rsurf = sband[ib][i][j] * 0.0001;  /* vs / 10000 */
                        rotoa = (rsurf * bttatmg[ib] / (1.0 - bsatm[ib] * rsurf)
                            + broatm[ib]) * btgo[ib];
                        raot550nm = taero[i][j];
                        pres = tp[i][j];
                        uwv = twvi[i][j];
                        uoz = tozi[i][j];
                        retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, ib,
                            pres, tpres, aot550nm, rolutt, transt, xtsstep,
                            xtsmin, xtvstep, xtvmin, sphalbt, tsmax, tsmin,
                            nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                            ogtransa0, ogtransa1, ogtransb0, ogtransb1,
                            ogtransc0, ogtransc1, wvtransa, wvtransb, wvtransc,
                            oztransa, rotoa, &roslamb, &tgo, &roatm, &ttatmg,
                            &satm, &xrorayp);
                        if (retval != SUCCESS)
                        {
                            sprintf (errmsg, "Performing lambertian "
                                "atmospheric correction type 2.");
                            error_handler (true, FUNC_NAME, errmsg);
                            exit (ERROR);
                        }

                        if (ib == 0)
                        {
                            if (roslamb < -0.005)
                            {
                                taero[i][j] = 0.05;
                                raot550nm = 0.05;
                                pres = tp[i][j];
                                uwv = twvi[i][j];
                                uoz = tozi[i][j];
                                retval = atmcorlamb2 (xts, xtv, xfi, raot550nm,
                                    ib, pres, tpres, aot550nm, rolutt, transt,
                                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                                    tsmax, tsmin, nbfic, nbfi, tts, indts, ttv,
                                    uoz, uwv, tauray, ogtransa0, ogtransa1,
                                    ogtransb0, ogtransb1, ogtransc0, ogtransc1,
                                    wvtransa, wvtransb, wvtransc, oztransa,
                                    rotoa, &roslamb, &tgo, &roatm, &ttatmg,
                                    &satm, &xrorayp);
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
                                if (raot550nm < 0.2)
                                {  /* Set the first aerosol bit */
                                    cloud[i][j] += 16;
                                }
                                else
                                {
                                    if (raot550nm < 0.5)
                                    {  /* Set the second aerosol bit */
                                        cloud[i][j] += 32;
                                    }
                                    else
                                    {  /* Set both aerosol bits */
                                        cloud[i][j] += 48;
                                    }
                                }
                            }  /* end if/else roslamb */
                        }  /* end if ib */

                        /* Save the surface reflectance value */
                        sband[ib][i][j] = (int) (roslamb * 10000.0);
                    }  /* end if */
                }  /* end if qaband */
            }  /* end for j */
        }  /* end for i */

        /* update status */
        printf ("100%%\n");
        fflush (stdout);
    }  /* end for ib */

    /* Write the data to the output file */
    printf ("Writing corrected data to the output files ...\n");

    /* Open the output file */
    sr_output = open_output (&xml_metadata, input);
    if (sr_output == NULL)
    {   /* error message already printed */
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Loop through the output bands */
    for (ib = 0; ib < NBAND_TTL_OUT; ib++)
    {
        /* Determine the number of bytes to be written per pixel */
        if (ib == SR_CLOUD)
        {
            for (i = 0; i < nlines; i++)
            {
                if (put_output_lines (sr_output, cloud[i], ib, i, 1,
                    sizeof (uint8)) != SUCCESS)
                {
                    sprintf (errmsg, "Writing output data for line %d", i);
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
            }
        }
        else
        {
            for (i = 0; i < nlines; i++)
            {
                if (put_output_lines (sr_output, sband[ib][i], ib, i, 1,
                    sizeof (int16)) != SUCCESS)
                {
                    sprintf (errmsg, "Writing output data for line %d", i);
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
            }
        }

    }

    /* Write the ENVI header for spectral indices files */
    printf ("Writing ENVI headers ...\n");
    for (ib = 0; ib < sr_output->nband; ib++)
    {
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
  
    /* Append the spectral index bands to the XML file */
    printf ("Appending metadata ...\n");
    if (append_metadata (sr_output->nband, sr_output->metadata.band,
        xml_infile) != SUCCESS)
    {
        sprintf (errmsg, "Appending spectral index bands to XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
  
    /* Close the reflectance product */
    printf ("Closing input/output and freeing pointers ...\n");
    close_input (input);
    free_input (input);

    /* Free the metadata structure */
    free_metadata (&xml_metadata);

    /* Close the output spectral indices products */
    close_output (sr_output);
    free_output (sr_output);

    /* Free the filename pointers */
    free (xml_infile);
    free (aux_infile);

    /* Free the data arrays */
    for (i = 0; i < DEM_NBLAT; i++)
        free (dem[i]);
    free (dem);

    for (i = 0; i < RATIO_NBLAT; i++)
    {
        free (ratiob1[i]);
        free (ratiob2[i]);
        free (ratiob7[i]);
    }
    free (ratiob1);
    free (ratiob2);
    free (ratiob7);

    for (i = 0; i < CMG_NBLAT; i++)
    {
        free (wv[i]);
        free (oz[i]);
    }
    free (wv);
    free (oz);

    for (i = 0; i < 16; i++)
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

    /* Allocate space for band data */
    for (i = 0; i < nlines; i++)
    {
        free (uband[i]);
        free (qaband[i]);
        free (aerob1[i]);
        free (aerob2[i]);
        free (aerob4[i]);
        free (aerob5[i]);
        free (aerob7[i]);
        free (twvi[i]);
        free (tozi[i]);
        free (tp[i]);
        free (tresi[i]);
        free (taero[i]);
        free (cloud[i]);
    }
    free (uband);
    free (qaband);
    free (aerob1);
    free (aerob2);
    free (aerob4);
    free (aerob5);
    free (aerob7);
    free (twvi);
    free (tozi);
    free (tp);
    free (tresi);
    free (taero);
    free (cloud);

    for (i = 0; i < NBAND_TTL_OUT; i++)
    {
        for (j = 0; j < nlines; j++)
            free (sband[i][j]);
        free (sband[i]);
    }
    free (sband);

    /* Free the spatial mapping pointer */
    free (space);

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
--------    ---------------  -------------------------------------
7/6/2014    Gail Schmidt     Original Development

NOTES:
******************************************************************************/
void usage ()
{
    printf ("l8_sr computes the surface reflectance values for the input "
            "Landsat 8 DN products.\n\n");
    printf ("usage: l8_sr "
            "--xml=input_xml_filename "
            "--aux=input_auxiliary_filename [--verbose]\n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -xml: name of the input XML file to be processed\n");
    printf ("    -aux: name of the input auxiliary file containing ozone "
        "and water vapor for the scene date.\n");

    printf ("\nwhere the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed? (default "
            "is false)\n");
    printf ("\nl8_sr --help will print the usage statement\n");
    printf ("\nExample: l8_sr --xml=LC80410272013181LGN00.xml "
            "--aux=L8ANC2013181.hdf_fused --verbose\n");
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

