/*****************************************************************************
FILE: lndpm.c
  
PURPOSE: Contains functions for reading the input XML metadata file and
creating the parameter and metadata files needed for downstream processing
by LEDAPS applications.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/16/2014    Gail Schmidt     Modified to use ESPA internal file format.
                              Removed support for non-LPGS file formats.
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
  1. The XML metadata format written via this library follows the ESPA internal
     metadata format found in ESPA Raw Binary Format v1.0.doc.  The schema for
     the ESPA internal metadata format is available at
     http://espa.cr.usgs.gov/static/schema/espa_internal_metadata_v1_0.xsd.
*****************************************************************************/
#include "lndpm.h"

FILE *LOG_FP = NULL;     /* pointer to the log file */
int open_log (char name[]);
void LOG (char *fmt, ...);
int conv_date (int *mm, int *dd, int yyyy);
int get_meta_from_lpgs (char *mtl_file, char *scene_name, char *acquisition_date, int *inst_type);
int find_file(char *path, char *name);
char *cm_name_for_id(int id);

int main (int argc, char *argv[])
{
    char FUNC_NAME[] = "lndpm";    /* function name */
    char errmsg[STR_SIZE];         /* error message */
    char input_mtl[STR_SIZE];      /* name of the input MTL file */
    char tmpstr[STR_SIZE];         /* temporary storage of the scene name */
    char tmpfile[STR_SIZE];        /* temporary storage of the MTL file */
    char scene_name[STR_SIZE];     /* name of the scene */
    char lndcal_name[STR_SIZE];    /* name of the lndcal input file */
    char lndsr_name[STR_SIZE];     /* name of the lndsr input file */
    char dem[STR_SIZE];            /* name of DEM file */
    char ozone[STR_SIZE];          /* name of ozone file */
    char reanalysis[STR_SIZE];     /* name of NCEP file */
    char path_buf[DIR_BUF_SIZE];   /* path to the auxillary/cal file */
    char cal_file[STR_SIZE];       /* name of the calibration file */
    char acquisition_date[MAX_STRING_LENGTH] = "\0";  /* acquisition date */
    char *anc_path = NULL;         /* path for LEDAPS ancillary data */
    char *token_ptr = NULL;        /* pointer used for obtaining scene name */
    char *file_ptr = NULL;         /* pointer used for obtaining file name */
    int year, month, day;          /* year, month, day of acquisition date */
    int inst = -1;                 /* instrument type */
    bool anc_missing = false;      /* is the ancillary data missing? */
    FILE *out = NULL;              /* pointer to the output header file */

    /* Open the log file */
    printf ("\nRunning lndpm ...\n");
    if (open_log ("lndpm : Landsat Metadata Parser") == ERROR)
    {
        sprintf (errmsg, "Create log report file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Check the command-line arguments and get the name of the XML file */
    if (argc != 2)
    {
        sprintf (errmsg, "Usage: lndpm <input_mtl_file>");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    strcpy (input_mtl, argv[1]);

    /* Get the scene name.  Strip off the path and file extension as it's
       assumed the XML filename is the scene name followed by the extension. */
    strcpy (tmpfile, input_mtl);
    file_ptr = strrchr (tmpfile, '/');
    if (file_ptr != NULL)
        file_ptr++;
    else
        file_ptr = tmpfile;
    token_ptr = strtok (file_ptr, ".");
    sprintf (scene_name, "%s", token_ptr);
    if (strstr (scene_name, "_MTL"))
    {
        /* remove last 4 chars "_MTL" from scene name */
        strncpy (tmpstr, scene_name, strlen(scene_name)-4);
        strcpy (scene_name, tmpstr);
    }
  
    /* Set up the names of the input files for downstream processing */
    sprintf (lndcal_name, "lndcal.%s.txt", scene_name);
    sprintf (lndsr_name, "lndsr.%s.txt", scene_name);

    /* Generate parameters for LEDAPS downstream apps from MTL metadata */
    if (get_meta_from_lpgs (input_mtl, scene_name, acquisition_date, &inst) !=
        SUCCESS)
    {
        sprintf (errmsg, "Reading metadata from input MTL: %s", input_mtl);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the header file for lndcal for writing */
    out = fopen (lndcal_name, "w");
    if (out == NULL)
    {
        sprintf (errmsg, "Opening lndcal header file for writing: %s",
            lndcal_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Get the path for the LEDAPS auxillary products (NCEP, TOMS, DEM, etc.)
       from the ANC_PATH environment variable.  If it isn't defined, then
       assume the products are in the local directory. */
    anc_path = getenv ("ANC_PATH");
    if (anc_path == NULL)
    {
        anc_path = ".";
        sprintf (errmsg, "ANC_PATH environment variable isn't defined. It is "
            "assumed the LEDAPS ancillary products will be available from "
            "the local directory.");
        error_handler (false, FUNC_NAME, errmsg);
    }

    /* Write the header data to the lndcal header file */
    fprintf (out, "PARAMETER_FILE\n");
    fprintf (out, "HEADER_FILE = %s.metadata.txt\n", scene_name);
    fprintf (out, "REF_FILE = lndcal.%s.hdf\n", scene_name);
    fprintf (out, "THERM_FILE = lndth.%s.hdf\n", scene_name);
    if (inst == LEDAPS_TM)
    {
        /* Find and prepare calibration files */
        strcpy (cal_file, "gold.dat");
        strcpy (path_buf, anc_path);
        if (find_file (path_buf, cal_file))
            strcpy (cal_file, path_buf);
        else
        {
            sprintf (errmsg, "Could not find gold.dat");
            error_handler (false, FUNC_NAME, errmsg);
        }
        fprintf (out, "GOLD_FILE = %s\n", cal_file);

        strcpy (cal_file, "gnew.dat");
        strcpy (path_buf, anc_path);
        if (find_file (path_buf, cal_file))
            strcpy (cal_file, path_buf);
        else
        {
            sprintf (errmsg, "Could not find gnew.dat");
            error_handler (false, FUNC_NAME, errmsg);
        }
        fprintf (out, "GNEW_FILE = %s\n", cal_file);

        strcpy (cal_file, "gold_2003.dat");
        strcpy (path_buf, anc_path);
        if (find_file (path_buf, cal_file))
            strcpy (cal_file, path_buf);
        else
        {
            sprintf (errmsg, "Could not find gold_2003.dat");
            error_handler (false, FUNC_NAME, errmsg);
        }
        fprintf (out, "GOLD_2003 = %s\n", cal_file);
    }
    fprintf (out, "LEDAPSVersion = %s\n", LEDAPS_VERSION);
    fprintf (out, "DNOUT = \"FALSE\"\n");
    fprintf (out, "RE_CAL = \"NO\"\n");
    fprintf (out, "END\n");
    fclose (out);

    /* Get year, month, and day from acquisition date */
    token_ptr = strtok (acquisition_date, "-");
    sscanf (token_ptr, "%d", &year);
    token_ptr = strtok (NULL, "-");
    sscanf (token_ptr, "%d", &month);
    token_ptr = strtok (NULL, "-");
    sscanf (token_ptr, "%d", &day);
  
    /* Convert to day of year */
    if (conv_date (&month, &day, year) != SUCCESS)
    {
        sprintf (errmsg, "Not able to convert the month, day, year from the "
            "acquisition date to Julian DOY.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find and prepare auxillary files */
    /* DEM file */
    strcpy (dem, "CMGDEM.hdf");
    strcpy (path_buf, anc_path);
    if (find_file (path_buf, dem))
    {
        strcpy (dem, path_buf);
        printf ("using DEM : %s\n", dem);
    }
    else
    {
        sprintf (errmsg, "Could not find DEM auxillary data: %s\n  Check "
            "ANC_PATH environment variable.", dem);
        error_handler (false, FUNC_NAME, errmsg);
        anc_missing = true;
    }

    /* TOMS ozone file */
    sprintf (ozone, "TOMS_%d%03d.hdf", year, day);
    strcpy (path_buf, anc_path);
    if (find_file (path_buf, ozone))
    {
        strcpy (ozone, path_buf);
        printf ("using TOMS : %s\n", ozone);
    }
    else
    {
        sprintf (errmsg, "Could not find TOMS auxillary data: %s\n  Check "
            "ANC_PATH environment variable.", ozone);
        error_handler (false, FUNC_NAME, errmsg);
        anc_missing = true;
    }
    
    /* NCEP file */
    sprintf (reanalysis, "REANALYSIS_%d%03d.hdf", year, day);
    strcpy (path_buf, anc_path);
    if (find_file (path_buf, reanalysis))
    {
        strcpy (reanalysis, path_buf);
        printf ("using REANALYSIS : %s\n", reanalysis);
    }
    else
    {
        sprintf (errmsg, "Could not find NCEP REANALYSIS auxillary data: %s\n"
            "  Check ANC_PATH environment variable.", reanalysis);
        error_handler (false, FUNC_NAME, errmsg);
        anc_missing = true;
    }

    /* Check to see if missing ancillary data */
    if (anc_missing)
    {
        sprintf (errmsg, "Verify the missing auxillary data products, then "
            "try reprocessing.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the header file for lndsr for writing */
    out = fopen (lndsr_name, "w");
    if (out == NULL)
    {
        sprintf (errmsg, "Opening lndsr header file for writing: %s",
            lndsr_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    fprintf (out, "PARAMETER_FILE\n");
    fprintf (out, "DEM_FILE = %s\n", dem);
    if (fopen (ozone, "r") != NULL)
    {
        /* if ozone file doesn't exist then don't write it to the header file
           and instead use climatology estimation */
        fprintf (out, "OZON_FIL = %s\n", ozone);
    }
    fprintf (out, "PRWV_FIL = %s\n", reanalysis);
    fprintf (out, "REF_FILE = lndcal.%s.hdf\n", scene_name);
    fprintf (out, "TEMP_FILE = lndth.%s.hdf\n", scene_name);
    fprintf (out, "SR_FILE = lndsr.%s.hdf\n", scene_name);
    fprintf (out, "META_FILE = %s\n", input_mtl);
    fprintf (out, "LEDAPSVersion = %s\n", LEDAPS_VERSION);
    fprintf (out, "END\n");
    fclose (out);
    fclose (LOG_FP);

    /* Successful completion */
    printf ("lndpm complete.\n");
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_meta_from_lpgs

PURPOSE: Obtains the necessary fields from the LPGS MTL metadata and generates
the metadata.txt file needed for downstream LEDAPS applications.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error with the LPGS MTL file
SUCCESS         Successfully generated the metadata.txt file

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/1/2013     Gail Schmidt     Made changes to the LPGS metadata tags to support
                              the new LPGS metadata format and tags.  This
                              software will also continue to support the old
                              metadata tags in case the user wants to use older
                              data in their personal archive.

NOTES:
******************************************************************************/
int get_meta_from_lpgs
(
    char *mtl_file,              /* I: MTL filename */
    char *scene_name,            /* I: scene name */
    char *acquisition_date,      /* O: acquisition date */
    int *inst_type               /* O: instrument type */
)
{
    char FUNC_NAME[] = "get_meta_from_lpgs";    /* function name */
    char errmsg[STR_SIZE];         /* error message */
    char buffer[STR_SIZE] = "\0";
    char tmpstr[STR_SIZE] = "\0";
    char tmpstr2[STR_SIZE] = "\0";
    char production_date[STR_SIZE] = "\0";
    char sensor[STR_SIZE] = "\0";
    char met_name[STR_SIZE] = "\0";
    char *label = NULL;
    char *token_ptr = NULL;
    char *seperator = "=\" \t";
    int i, year, pyear, zone;
    int bno;                    /* band number */
    float  lmax[8], lmin[8], qmin[8], qmax[8];
    double res;                    /* resolution of the image */
    float vert_long_from_pole, lat_true_scale, false_easting, false_northing;
    float fnum;
    float gain[8], bias[8];        /* gain and bias values for the bands */
    FILE *mtl_ptr=NULL, *meta_ptr=NULL;
    bool param_error;              /* was there an error with the params */

    /* metadata */
    int   gctp_proj_num=-99;          /* GCTP projection number */
    int   chkFlg[NMETA];
    char  meta_list[NMETA][STR_SIZE];
    char  band_fname[8][STR_SIZE];
    char  band_gain[8][STR_SIZE];
    char  ul_xy[2][STR_SIZE];

    /* NOTE: Make sure to keep the #defines in lndpm.h up-to-date to match the
       metadata fields in this list */
    char  *meta_lut[NMETA][3] =
    {
        {"HEADER_FILE",       "", ""},
        {"FILE_TYPE",         "OUTPUT_FORMAT", ""},
        {"DATA_PROVIDER",     "STATION_ID", "GROUND_STATION"},
        {"SATELLITE",         "SPACECRAFT_ID", ""},
        {"INSTRUMENT",        "SENSOR_ID", ""},
        {"ACQUISITION_DATE",  "DATE_ACQUIRED", "ACQUISITION_DATE"},
        {"ACQUISITION_TIME",  "SCENE_CENTER_TIME", "SCENE_CENTER_SCAN_TIME"},
        {"PRODUCTION_DATE",   "FILE_DATE", "PRODUCT_CREATION_TIME"},
        {"SOLAR_ZENITH",      "SUN_ELEVATION", ""},
        {"SOLAR_AziMUTH",     "SUN_AZIMUTH", ""},
        {"WRS_SYSTEM",        "", ""},
        {"WRS_PATH",          "WRS_PATH", ""},
        {"WRS_ROW",           "WRS_ROW", "STARTING_ROW"},
        {"NBAND",             "", ""},
        {"BANDS",             "", ""},
        {"GAIN_SETTINGS",     "", ""},
        {"NSAMPLE",           "REFLECTIVE_SAMPLES", "PRODUCT_SAMPLES_REF"},
        {"NLINE",             "REFLECTIVE_LINES", "PRODUCT_LINES_REF"},
        {"FILE_NAMES",        "", ""},
        {"NBAND_TH",          "", ""},
        {"BANDS_TH",          "", ""},
        {"GAIN_SETTINGS_TH",  "GAIN_BAND_6_VCID_1", "BAND6_GAIN1"},
        {"NSAMPLE_TH",        "THERMAL_SAMPLES", "PRODUCT_SAMPLES_THM"},
        {"NLINE_TH",          "THERMAL_LINES", "PRODUCT_LINES_THM"},
        {"FILE_NAMES_TH",     "", ""},
        {"PROJECTION_NUMBER", "MAP_PROJECTION", ""},
        {"PIXEL_SIZE",        "GRID_CELL_SIZE_REFLECTIVE",
         "GRID_CELL_SIZE_REF"},
        {"UPPER_LEFT_CORNER", "", ""},
        {"PROJECTION_ZONE",   "UTM_ZONE", "ZONE_NUMBER"},
        {"PROJECTION_SPHERE", "", ""},
        {"GAIN",              "", ""},
        {"BIAS",              "", ""},
        {"GAIN_TH",           "", ""},
        {"BIAS_TH",           "", ""},
        {"VERTICAL_LON_FROM_POLE", "VERTICAL_LON_FROM_POLE",
         "VERTICAL_LONGITUDE_FROM_POLE"},
        {"TRUE_SCALE_LAT",    "TRUE_SCALE_LAT", "LATITUDE_OF_TRUE_SCALE"},
        {"FALSE_EASTING",     "FALSE_EASTING", ""},
        {"FALSE_NORTHING",    "FALSE_NORTHING", ""},
        {"PROJECTION_PARAMETERS", "", ""},
        {"END",               "", ""}
    };

    /* Open the MTL file for reading data needed by the output metadata file */
    mtl_ptr = fopen (mtl_file, "r");
    if (mtl_ptr == NULL)
    {
        sprintf (errmsg, "Can't open the input MTL file: %s", mtl_file);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Open the metadata.txt file */
    sprintf (met_name, "%s.metadata.txt", scene_name);
    meta_ptr = fopen (met_name, "w");
    if (meta_ptr == NULL)
    {
        sprintf (errmsg, "Can't open output header file for write: %s",
            met_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Initialize the metadata flags and fields */
    for (i = 0; i < NMETA; i++)
        chkFlg[i] = 0;
    for (i = 0; i < MAX_N_BANDS; i++)
    {
        qmin[i] = 0.0;
        qmax[i] = 255.0;
    }

    /* Process the MTL file line by line */
    while (fgets (buffer, STR_SIZE, mtl_ptr) != NULL)
    {
        /* Get string token */
        token_ptr = strtok (buffer, seperator);
        label = token_ptr;
 
        while (token_ptr != NULL)
        {
            token_ptr = strtok (NULL, seperator);

            /* Information provided just after label, for some tags we are
               supporting both the old and the new metadata tags */
            for (i = 1; i < NMETA-1; i++)  /* don't look for start/end labels */
            {
                if (strcmp (label, meta_lut[i][1]) == 0 ||
                    strcmp (label, meta_lut[i][2]) == 0)
                {
                    chkFlg[i] = 1;
                    switch (i)
                    { /* some metadata need more work to satisfy input */
                        case FILE_TYPE:
                            /* Regardless of the file type for the LPGS product
                               LEDAPS will utilize the ESPA internal format */
                            sprintf (meta_list[i], "%s = BINARY\n",
                                meta_lut[i][0]);
                            break;
                        case SATELLITE:
                            if (strcmp (token_ptr, "LANDSAT_7") == 0 ||
                                strcmp (token_ptr, "Landsat7") == 0)
                                sprintf (meta_list[i], "%s = LANDSAT_7\n",
                                    meta_lut[i][0]);
                            else if (strcmp (token_ptr, "LANDSAT_5") == 0 ||
                                     strcmp (token_ptr, "Landsat5") == 0)
                                sprintf (meta_list[i], "%s = LANDSAT_5\n",
                                    meta_lut[i][0]);
                            else if (strcmp (token_ptr,"LANDSAT_4") == 0 ||
                                     strcmp (token_ptr, "Landsat4") == 0)
                                sprintf (meta_list[i], "%s = LANDSAT_4\n",
                                    meta_lut[i][0]);
                            else
                                chkFlg[i] = 0;
                            break;
                        case INSTRUMENT:
                            strcpy (sensor, token_ptr);
                            if (strstr (token_ptr, "ETM") != NULL ||
                                strstr (token_ptr, "ETM+") != NULL)
                            {
                                sprintf (meta_list[i], "%s = ETM\n",
                                    meta_lut[i][0]);
                                *inst_type = LEDAPS_ETM;
                            }
                            else if (strstr (token_ptr, "TM") != NULL)
                            {
                                sprintf (meta_list[i], "%s = TM\n",
                                    meta_lut[i][0]);
                                *inst_type = LEDAPS_TM;
                            }
                            else
                                chkFlg[i] = 0;
                            break;
                        case PRODUCTION_DATE:
                            sscanf (token_ptr, "%10s", production_date);
                            break;
                        case SOLAR_ZENITH:
                            sscanf (token_ptr, "%f", &fnum);
                            sprintf (meta_list[i], "%s = %f\n", meta_lut[i][0],
                                90-fnum);
                            break;
                        case GAIN_SETTINGS_TH:
                            if (strcmp (token_ptr, "H") == 0)
                                sprintf (meta_list[i], "%s = HIGH\n",
                                    meta_lut[i][0]);
                            else
                                sprintf (meta_list[i], "%s = LOW\n",
                                    meta_lut[i][0]);
                            break;
                        case PROJECTION_NUMBER:
                            if (strstr (token_ptr, "UTM"))
                            {
                                gctp_proj_num = GCTP_UTM_PROJ;
                                sprintf (meta_list[i], "%s = %d\n",
                                    meta_lut[i][0], gctp_proj_num);
                            }
                            else if (strstr (token_ptr, "PS"))
                            {
                                gctp_proj_num = GCTP_PS_PROJ;
                                sprintf (meta_list[i], "%s = %d\n",
                                    meta_lut[i][0], gctp_proj_num);
                            }
                            else
                            {
                                sprintf (errmsg, "Please check projection.  "
                                    "Only UTM and PS are supported.");
                                error_handler (true, FUNC_NAME, errmsg);
                                return (ERROR);
                            }
                            break; 
                        case PIXEL_SIZE:
                            sscanf (token_ptr, "%lf", &res);
                            sprintf (meta_list[i], "%s = %f\n", meta_lut[i][0],
                                res);
                            break;
                        case PROJECTION_ZONE:
                            sscanf (token_ptr, "%d", &zone);
                            sprintf (meta_list[i], "%s = %d\n", meta_lut[i][0],
                                zone);
                            break;
                        case VERT_LON_FROM_POLE:
                            sscanf (token_ptr, "%f", &vert_long_from_pole);
                            sprintf (meta_list[i], "%s = %f\n", meta_lut[i][0],
                                vert_long_from_pole);
                            break;
                        case TRUE_SCALE_LAT:
                            sscanf (token_ptr, "%f", &lat_true_scale);
                            sprintf (meta_list[i],"%s = %f\n", meta_lut[i][0],
                                lat_true_scale);
                            break;
                        case FALSE_EASTING:
                            sscanf (token_ptr, "%f", &false_easting);
                            sprintf (meta_list[i],"%s = %f\n", meta_lut[i][0],
                                false_easting);
                            break;
                        case FALSE_NORTHING:
                            sscanf (token_ptr, "%f", &false_northing);
                            sprintf (meta_list[i],"%s = %f\n", meta_lut[i][0],
                                false_northing);
                            break;
                        default:  /* other parameters can be used directly */
                            sscanf (token_ptr, "%s", tmpstr);
                            sprintf (meta_list[i],"%s = %s\n", meta_lut[i][0],
                                tmpstr);
                    }  /* end switch */

                    /* save acquisition date for later use */
                    if (i == ACQUISITION_DATE)
                        strcpy (acquisition_date, tmpstr);
                }  /* end if strcmp */
            }  /* end for */
 
            /* Extract coordinates of upleft corner */
            if (strcmp (label, "CORNER_UL_PROJECTION_X_PRODUCT") == 0 ||
                strcmp (label, "PRODUCT_UL_CORNER_MAPX") == 0 ||
                strcmp (label, "SCENE_UL_CORNER_MAPX") == 0)
            {
                chkFlg[UPPER_LEFT_CORNER] = 1; 
                sscanf (token_ptr, "%s", ul_xy[0]);
            }

            if (strcmp (label, "CORNER_UL_PROJECTION_Y_PRODUCT") == 0 ||
               strcmp (label, "PRODUCT_UL_CORNER_MAPY") == 0 ||
               strcmp (label, "SCENE_UL_CORNER_MAPY") == 0)
                sscanf(token_ptr, "%s", ul_xy[1]);

            /* Get thermal band filename */
            if (*inst_type == LEDAPS_ETM)
            {
                if (strcmp (label, "FILE_NAME_BAND_6_VCID_1") == 0 ||
                    strcmp (label, "BAND61_FILE_NAME") == 0 ||
                    strcmp (label, "BAND6L_FILE_NAME") == 0)
                {
                    sscanf (token_ptr, "%s", tmpstr);
                    sprintf(meta_list[FILE_NAMES_TH],"%s = %s\n",
                        meta_lut[FILE_NAMES_TH][0], tmpstr);
                    chkFlg[FILE_NAMES_TH] = 1;
                }
            }
            else
            {
                if (strcmp (label, "FILE_NAME_BAND_6") == 0 ||
                    strcmp (label, "BAND6_FILE_NAME") == 0)
                {
                    sscanf (token_ptr, "%s", tmpstr);
                    sprintf (meta_list[FILE_NAMES_TH],"%s = %s\n",
                        meta_lut[FILE_NAMES_TH][0], tmpstr);
                    chkFlg[FILE_NAMES_TH] = 1;
                }
            }

            /* Extract gain settings and file names for ETM+ band 1-5 and 7 */
            for (i = 1; i <= 7; i++)
            {
                if (*inst_type == LEDAPS_ETM)
                {
                    if (i != 6)
                    {
                        sprintf (tmpstr, "FILE_NAME_BAND_%d", i);
                        sprintf (tmpstr2, "BAND%d_FILE_NAME", i);
                        if (strcmp (label, tmpstr) == 0 ||
                            strcmp (label, tmpstr2) == 0)
                        {
                            sscanf (token_ptr, "%s", band_fname[i]);
                            chkFlg[FILE_NAMES] = 1;
                        }

                        sprintf (tmpstr,"GAIN_BAND_%d",i);
                        sprintf (tmpstr2,"BAND%d_GAIN",i);
                        if (strcmp (label, tmpstr) == 0 ||
                            strcmp (label, tmpstr2) == 0)
                        {
                            if (strcmp(token_ptr,"H") == 0)
                                strcpy (band_gain[i], "HIGH");
                            else
                                strcpy (band_gain[i], "LOW");
                            chkFlg[GAIN_SETTINGS] = 1;
                        }      
                    }
                }
            }

            /* Extract filename and LMIN, LMAX, QCALMIN and QCALMAX */
            for (i = 1; i <= 7; i++)
            {
                if (*inst_type == LEDAPS_ETM && i == 6)
                {
                    /* use low gain thermal band for ETM+ */
                    bno = 61;
                    sprintf (tmpstr, "FILE_NAME_BAND_6_VCID_1");
                    sprintf (tmpstr2, "BAND%d_FILE_NAME", bno);
                }
                else
                {
                    bno = i;
                    sprintf (tmpstr, "FILE_NAME_BAND_%d", bno);
                    sprintf (tmpstr2, "BAND%d_FILE_NAME", bno);
                }

                /* get file names */
                if (strcmp (label, tmpstr) == 0 ||
                    strcmp (label, tmpstr2) == 0)
                {
                    sscanf(token_ptr, "%s", band_fname[i]);
                    if (bno != 61)
                        chkFlg[FILE_NAMES] = 1;
                    else
                        chkFlg[FILE_NAMES_TH] = 1;
                }

                if (bno == 61)
                    sprintf (tmpstr, "RADIANCE_MAXIMUM_BAND_6_VCID_1");
                else
                    sprintf (tmpstr, "RADIANCE_MAXIMUM_BAND_%d", bno);
                sprintf (tmpstr2, "LMAX_BAND%d", bno);
                if (strcmp (label, tmpstr) == 0 || strcmp (label, tmpstr2) == 0)
                {
                    sscanf (token_ptr, "%f", &(lmax[i]));
                    chkFlg[GAIN] = chkFlg[BIAS] = 1;
                    chkFlg[GAIN_TH] = chkFlg[BIAS_TH] = 1;
                    /* disable gain settings */
                    chkFlg[GAIN_SETTINGS] = chkFlg[GAIN_SETTINGS_TH] = -1;
                }

                if (bno == 61)
                    sprintf (tmpstr, "RADIANCE_MINIMUM_BAND_6_VCID_1");
                else
                    sprintf (tmpstr, "RADIANCE_MINIMUM_BAND_%d", bno);
                sprintf (tmpstr2, "LMIN_BAND%d", bno);
                if (strcmp (label, tmpstr) == 0 ||
                    strcmp(label, tmpstr2) == 0)
                    sscanf (token_ptr, "%f", &(lmin[i]));

                if (bno == 61)
                    sprintf (tmpstr, "QUANTIZE_CAL_MIN_BAND_6_VCID_1");
                else
                    sprintf (tmpstr, "QUANTIZE_CAL_MIN_BAND_%d", bno);
                sprintf (tmpstr2, "QCALMIN_BAND%d", bno);
                if (strcmp (label, tmpstr) == 0 || strcmp (label, tmpstr2) == 0)
                    sscanf (token_ptr, "%f", &(qmin[i]));

                if (bno == 61)
                    sprintf (tmpstr, "QUANTIZE_CAL_MAX_BAND_6_VCID_1");
                else
                    sprintf (tmpstr, "QUANTIZE_CAL_MAX_BAND_%d", bno);
                sprintf (tmpstr2, "QCALMAX_BAND%d", bno);
                if (strcmp(label, tmpstr) == 0 || strcmp (label, tmpstr2) == 0)
                    sscanf (token_ptr, "%f", &(qmax[i]));
            }

            /* In case label (key words) is not the first word in a line */
            label = token_ptr;
        }  /* end while token_ptr */
    }  /* end while fgets */

    /* Start and end tags */
    sprintf (meta_list[0], "%s\n", meta_lut[0][0]);
    sprintf (meta_list[NMETA-1], "%s\n", meta_lut[NMETA-1][0]);
    chkFlg[0] = 1;
    chkFlg[NMETA-1] = 1;

    /* Some parameters aren't available in the metadata file */
    sprintf (meta_list[DATA_PROVIDER], "%s = USGS/EROS\n",
      meta_lut[DATA_PROVIDER][0]);
    sprintf (meta_list[WRS_SYSTEM], "%s = 2\n", meta_lut[WRS_SYSTEM][0]);
    sprintf (meta_list[NBAND], "%s = 6\n", meta_lut[NBAND][0]);
    sprintf (meta_list[BANDS], "%s = 1, 2, 3, 4, 5, 7\n", meta_lut[BANDS][0]);
    sprintf (meta_list[NBAND_TH], "%s = 1\n", meta_lut[NBAND_TH][0]);
    sprintf (meta_list[BANDS_TH], "%s = 6\n", meta_lut[BANDS_TH][0]);
    sprintf (meta_list[PROJECTION_SPHERE], "%s = 12\n",
        meta_lut[PROJECTION_SPHERE][0]);  /* WGS84=12 in GCTPC */
    chkFlg[DATA_PROVIDER] = 1;
    chkFlg[WRS_SYSTEM] = 1;
    chkFlg[NBAND] = 1;
    chkFlg[BANDS] = 1;
    chkFlg[NBAND_TH] = 1;
    chkFlg[BANDS_TH] = 1;
    chkFlg[PROJECTION_SPHERE] = 1;

    /* Write band specific parameters */
    /* If defined lmin/lmax/qcalmin/qcalmax then use it instead of gain
       setting */
    if (chkFlg[GAIN] == 1)
    {
        /* Disable gain settings */
        chkFlg[GAIN_SETTINGS] = -1;
        chkFlg[GAIN_SETTINGS_TH] = -1;

        /* Use actual defined values */
        for (i = 1; i <= 7; i++)
        {
            gain[i] = (lmax[i] - lmin[i]) / (qmax[i] - qmin[i]);
            bias[i] = lmin[i] -gain[i] * qmin[i];
        }
        sprintf (meta_list[GAIN], "%s = %f, %f, %f, %f, %f, %f\n",
            meta_lut[GAIN][0], gain[1], gain[2], gain[3], gain[4], gain[5],
            gain[7]);
        sprintf (meta_list[BIAS], "%s = %f, %f, %f, %f, %f, %f\n",
            meta_lut[BIAS][0], bias[1], bias[2], bias[3], bias[4], bias[5],
            bias[7]);
        sprintf (meta_list[GAIN_TH], "%s = %f\n", meta_lut[GAIN_TH][0],
            gain[6]);
        sprintf (meta_list[BIAS_TH], "%s = %f\n", meta_lut[BIAS_TH][0],
            bias[6]);
    }
    else
    {
        sprintf (meta_list[GAIN_SETTINGS], "%s = %s, %s, %s, %s, %s, %s\n",
            meta_lut[GAIN_SETTINGS][0], band_gain[1],band_gain[2],
            band_gain[3],band_gain[4],band_gain[5],band_gain[7]);

        /* Use HIGH or LOW gain setting and set gain and bias flag invalid */
        chkFlg[GAIN] = -1;
        chkFlg[BIAS] = -1;
        chkFlg[GAIN_TH] = -1;
        chkFlg[BIAS_TH] = -1;
    } 
    sprintf (meta_list[FILE_NAMES],"%s = %s, %s, %s, %s, %s, %s\n",
        meta_lut[FILE_NAMES][0], band_fname[1],band_fname[2],
        band_fname[3],band_fname[4],band_fname[5],band_fname[7]); 

    /* Write coordinates of upleft corner.  Always convert to the UL coordinate
       system since the corner points listed in the MTL file represent the
       center of the pixel. */
    sprintf(meta_list[UPPER_LEFT_CORNER],"%s = %f, %f\n", 
        meta_lut[UPPER_LEFT_CORNER][0], atof (ul_xy[0]) - 0.5*res,
        atof (ul_xy[1]) + 0.5*res); 

    /* Some inputs have production year prior to acquisition date and need
       special care here. Get year from acquisition date and production date. */
    sscanf (acquisition_date, "%d", &year);
    sscanf (production_date, "%d", &pyear);
    if (pyear < year)
        strcpy (production_date, acquisition_date);
    sprintf (meta_list[PRODUCTION_DATE], "%s = %s\n",
        meta_lut[PRODUCTION_DATE][0], production_date);

    /* If this is the UTM projection then disable the PS projection fields */
    if (gctp_proj_num == GCTP_UTM_PROJ)
    {
        chkFlg[VERT_LON_FROM_POLE] = -1;
        chkFlg[TRUE_SCALE_LAT] = -1;
        chkFlg[FALSE_EASTING] = -1;
        chkFlg[FALSE_NORTHING] = -1;
        chkFlg[PROJECTION_PARAMS] = -1;
    }

    /* If this is the PS projection then disable the UTM zone field and set up
       the projection parameters */
    if (gctp_proj_num == GCTP_PS_PROJ)
    {
        chkFlg[PROJECTION_ZONE] = -1;
        sprintf (meta_list[PROJECTION_PARAMS], "%s = 0.0, 0.0, 0.0, 0.0, %f, "
            "%f, %f, %f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0\n",
            meta_lut[PROJECTION_PARAMS][0], vert_long_from_pole, lat_true_scale,
            false_easting, false_northing);
        chkFlg[PROJECTION_PARAMS] = 1;
    }

    /* Write all parameters to input file */
    param_error = false;
    for (i = 0; i < NMETA; i++)
    {
        if (chkFlg[i] == 1)
        {
            fprintf (meta_ptr, "%s", meta_list[i]);
            LOG ("FOUND: %s", meta_list[i]);
        }
        else if (chkFlg[i] == 0)
        {
            fprintf (meta_ptr, "%s = \n", meta_lut[i][0]);
            sprintf (errmsg, "Metadata field %s was not found in metadata.",
                meta_lut[i][0]);
            error_handler (false, FUNC_NAME, errmsg);
            param_error = true;
        }
    }
    if (param_error)
        return (ERROR);

    /* Close the file pointers */
    if (mtl_ptr)
        fclose (mtl_ptr);
    if (meta_ptr)
        fclose (meta_ptr);

    /* Successful completion */
    return (SUCCESS);
}  /* end get_meta_from_lpgs */


/******************************************************************************
MODULE:  cm_name_for_id

PURPOSE: Given the ID for a metadata.txt param name, return a pointer to the
name of the parameter.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
NULL            ID not found
not-NULL        Name of parameter for the specified ID

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
******************************************************************************/
char *cm_name_for_id
(
    int id           /* I: ID for the parameter name to be returned; follow
                           the #defines from lndpm.h for the parameter IDs */
)
{
    int ii = 0;                           /* looping variable */
    char FUNC_NAME[] = "cm_name_for_id";  /* function name */
    char errmsg[STR_SIZE];                /* error message */

    /* Make sure these match the #defines for the header file provided in
       lndpm.h */
    static struct cm_symbol
    {
      int  id;
      char *name;
    } symbols[] = {
      {FILE_TYPE,          "FILE_TYPE"},
      {DATA_PROVIDER,      "DATA_PROVIDER"},
      {SATELLITE,          "SATELLITE"},
      {INSTRUMENT,         "INSTRUMENT"},
      {ACQUISITION_DATE,   "ACQUISITION_DATE"},
      {ACQUISITION_TIME,   "ACQUISITION_TIME"},
      {PRODUCTION_DATE,    "PRODUCTION_DATE"},
      {SOLAR_ZENITH,       "SOLAR_ZENITH"},
      {SOLAR_AZIMUTH,      "SOLAR_AziMUTH"},      /* yes, case is correct */
      {WRS_SYSTEM,         "WRS_SYSTEM"},
      {WRS_PATH,           "WRS_PATH"},
      {WRS_ROW,            "WRS_ROW"},
      {NBAND,              "NBAND"},
      {BANDS,              "BANDS"},
      {GAIN_SETTINGS,      "GAIN_SETTINGS"},
      {NSAMPLE,            "NSAMPLE"},
      {NLINE,              "NLINE"},
      {FILE_NAMES,         "FILE_NAMES"},
      {NBAND_TH,           "NBAND_TH"},
      {BANDS_TH,           "BANDS_TH"},
      {GAIN_SETTINGS_TH,   "GAIN_SETTINGS_TH"},
      {NSAMPLE_TH,         "NSAMPLE_TH"},
      {NLINE_TH,           "NLINE_TH"},
      {FILE_NAMES_TH,      "FILE_NAMES_TH"},
      {PROJECTION_NUMBER,  "PROJECTION_NUMBER"},
      {PIXEL_SIZE,         "PIXEL_SIZE"},
      {UPPER_LEFT_CORNER,  "UPPER_LEFT_CORNER"},
      {PROJECTION_ZONE,    "PROJECTION_ZONE"},
      {PROJECTION_SPHERE,  "PROJECTION_SPHERE"},
      {GAIN,               "GAIN"},
      {BIAS,               "BIAS"},
      {GAIN_TH,            "GAIN_TH"},
      {BIAS_TH,            "BIAS_TH"},
      {VERT_LON_FROM_POLE, "VERTICAL_LON_FROM_POLE"},
      {TRUE_SCALE_LAT,     "TRUE_SCALE_LAT",},
      {FALSE_EASTING,      "FALSE_EASTING"},
      {FALSE_NORTHING,     "FALSE_NORTHING"},
      {PROJECTION_PARAMS,  "PROJECTION_PARAMETERS"},
    };
    
    /* Loop through the available parameters and find the specified parameter
       ID */
    for (ii = 0; ii < MAX_CM_PARAMS; ++ii)
    {
      if (symbols[ii].id != id)
        continue;
      return (symbols[ii].name);
    }
    
    /* Error in processing */
    sprintf (errmsg, "No parameter name for ID %d", id);
    error_handler (true, FUNC_NAME, errmsg);
    return NULL;
}


/******************************************************************************
MODULE:  cm_val_for_id

PURPOSE: Given the ID and metadata structure, return the value from the
metadata for this parameter.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
NULL            ID not found
not-NULL        Metadata value for the specified ID

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
******************************************************************************/
char *cm_val_for_id
(
    int id,             /* I: ID for the parameter name to be returned; follow
                              the #defines from lndpm.h for the parameter IDs */
    METADATA *metadata  /* I: metadata values for each of the parameters */
)
{
    char FUNC_NAME[] = "cm_val_for_id";    /* function name */
    char errmsg[STR_SIZE];                 /* error message */

    /* Validate the ID is within bounds */
    if (id < 1 || id >= MAX_CM_PARAMS)
    {
        sprintf (errmsg, "Invalid ID %d", id);
        error_handler (true, FUNC_NAME, errmsg);
        return NULL;
    }

    /* Return the value of the metadata for that ID */
    return metadata[id].val;
}


/******************************************************************************
MODULE:  conv_date

PURPOSE: Convert year-mm-dd to Julian day of year (DOY) or vice versa.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error with the conversion
SUCCESS         Successfully converted the date

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
             Feng Gao         Revised from program "jdoy.c" by MODIS LDOPE QA
                              team
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
  1. *mm and *dd are input/output varialbes
     Input julian day number:
        input *mm = 0;
        input *dd = julian day
        output in mm-dd-yyyy format
     Input in mm-dd-yyyy format
        output *mm = 0;
        output *dd = julian day number
******************************************************************************/
int conv_date
(
    int *mm,         /* I/O: month value */
    int *dd,         /* I/O: day of month or Julian day value */
    int yyyy         /* I: year value */
)
{
    char FUNC_NAME[] = "conv_date";    /* function name */
    char errmsg[STR_SIZE];             /* error message */
    int nm, im;        /* looping vars */
    int ndays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
                       /* number of days in each month */

    /* Check for leap year */
    if ((yyyy%400 == 0) || ((yyyy%4 == 0) && (yyyy%100 != 0))) ndays[1] = 29; 

    /* Do the conversion */
    if (*mm == 0)
    { /* Convert from Julian DOY to mm-dd-yyyy */
        if (*dd <= 0) {
            sprintf (errmsg, "Invalid input date: %d %d", *dd, yyyy);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        else
        {
            for (im = 0; ((im < 12) && (*dd > 0)); im++)
                *dd -= ndays[im];
            if ((im > 12) || ((im == 12) && (*dd > 0)))
            {
                sprintf (errmsg, "Invalid input date: %d %d", *dd, yyyy);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            else
            {
                *mm = im;
                *dd += ndays[*mm - 1];
            }
        }
    }
    else
    { /* Convert from mm-dd-yyyy to Julian DOY */
        if ((*mm <= 0) || (*dd <= 0))
        {
            sprintf (errmsg, "Invalid input date: %d-%d-%d", *mm, *dd, yyyy);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        else
        {
            nm = *mm - 1;
            for (im = 0; im < nm; im++)
                *dd += ndays[im];
            *mm = 0;
        }
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  open_log

PURPOSE: Open the log file as a text file for appending and add a time stamp.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error opening log file
SUCCESS         Successfully opened log file

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
1. LOG_FP and logFile are global variables for all routines in this file to
   log to.
******************************************************************************/
int open_log
(
    char app_name[]     /* I: application name/info string */
)
{
    struct tm *currtime;    /* current time in GMT */
    time_t t;               /* current time variable */
    char str[STR_SIZE];     /* current time in HMS */
    char FUNC_NAME[] = "open_log";    /* function name */
    char errmsg[STR_SIZE];            /* error message */

    /* Open the log file for appending */
    LOG_FP = fopen (logFile, "a");
    if (LOG_FP == NULL)
    {
        sprintf (errmsg, "Cannot open log file for appending: %s", logFile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Log the time and application name or string */
    t = time(NULL);
    currtime = (struct tm *) gmtime (&t);
    strftime (str, 100, "%FT%H:%M:%SZ", currtime);
    LOG ("\n\n##################################\n");
    LOG ("Starting <%s> on %s\n", app_name, str);
    return (SUCCESS);
}


/******************************************************************************
MODULE:  LOG

PURPOSE: Log a message to the log file

RETURN VALUE: None

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
1. LOG_FP is a global variable for all routines in this file to log to.
******************************************************************************/
void LOG
(
    char *fmt, ...   /* format of log information */
)
{
    va_list ap;
    
    va_start (ap, fmt);
    vfprintf (LOG_FP, fmt, ap);
    va_end (ap);
}


/******************************************************************************
MODULE:  scan_dir

PURPOSE: Search the specified directory, recursively, for the specified
filename.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
non-zero        File is found, and path points to full path
zero            File is not found, and path is unchanged

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
******************************************************************************/
int scan_dir
(
    char *path,   /* I: directory/path to search; upon successfully finding
                        the file, path contains the location of the file */
    char *name    /* I: filename for which to search */
)
{
    char FUNC_NAME[] = "scan_dir";    /* function name */
    char errmsg[STR_SIZE];             /* error message */
    DIR *fd;                  /* pointer to directory */
    struct dirent *dirent_p;  /* pointer to directory */
    char *nbp;                /* pointer to the end of the path */
    int found = 0;            /* was the file found? */
    
    /* Add directory separator to end of directory name if not already there */
    nbp = path + strlen (path);
    if (*nbp != '/')
        *nbp++ = '/';
      
    if (nbp+MAXNAMLEN+2 >= path + DIR_BUF_SIZE)
    {
        sprintf (errmsg, "Path name too long -- cannot search: %s", path);
        error_handler (true, FUNC_NAME, errmsg);
        return (found);
    }

    /* Check the directory to see if it's accessible */
    if ((fd = opendir (path)) == NULL)
    {
        sprintf (errmsg, "Could not read directory: %s", path);
        error_handler (true, FUNC_NAME, errmsg);
        return (found);
    }

    /* Recursively search for the file while traversing through the directory
       structure */
    while ((dirent_p = readdir(fd)) != NULL)   /* search directory */
    {
        if (dirent_p->d_ino == 0)                /* slot not in use */
            continue;
        if (strcmp (dirent_p->d_name, ".") == 0   /* ignore current ... */
        || strcmp (dirent_p->d_name, "..") == 0)  /* and parent directory */
            continue;
        
        strcat (path, dirent_p->d_name);          /* check this path */
        if ((found = find_file (path, name)) != 0)
            break;                                /* found it */
        else
            *nbp = '\0';                          /* restore directory name */
    }

    /* Close the directory pointer and successful completion */
    closedir (fd);
    return (found);
}


/******************************************************************************
MODULE:  find_file

PURPOSE: Look in the current directory for the specified file.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
non-zero        File is found, and path points to full path
zero            File is not found, and path is unchanged

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
******************************************************************************/
int find_file
(
    char *path,   /* I: directory/path to search; upon successfully finding
                        the file, path contains the location of the file */
    char *name    /* I: filename for which to search */
)
{
    char FUNC_NAME[] = "find_file";    /* function name */
    char errmsg[STR_SIZE];             /* error message */
    struct stat stbuf;                 /* buffer for file/directory stat */
    char pbuf[DIR_BUF_SIZE] = {0};     /* path buffer */
    int found = 0;                     /* was the file found? */
    
    /* This is the path we are checking */
    strcpy (pbuf, path);

    /* Make sure the path exists */
    if (stat (pbuf, &stbuf) != 0)
    {
        sprintf (errmsg, "Can't stat directory: %s", pbuf);
        error_handler (true, FUNC_NAME, errmsg);
        return (found);
    }

    /* If this is a directory, then search it.  Otherwise it's a file so
       check it for our filename. */
    if ((stbuf.st_mode & S_IFMT) == S_IFDIR)
        found = scan_dir (pbuf, name);
    else
        found = (strcmp(pbuf + strlen(pbuf) - strlen(name), name) == 0);

    /* If the file was found remember the location of that file */
    if (found)
        strcpy (path, pbuf);

    return (found);
}

