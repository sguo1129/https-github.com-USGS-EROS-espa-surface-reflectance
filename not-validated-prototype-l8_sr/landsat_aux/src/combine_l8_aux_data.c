/******************************************************************************
FILE: combine_l8_aux_data.c
  
PURPOSE: Contains functions for reading the daily, global Aqua and Terra CMG
and CMA files and "fusing" them into a single output HDF file.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
Unknown      Jim Ray          Original development by Jim Ray, Sigma Space Corp.
8/26/2014    Gail Schmidt     Update of the original code delivered by Eric
                              Vermote, NASA GSFC, to improve the documentation
                              and error handling for use within the ESPA
                              system

NOTES:
  1. MODIS CMG files are daily global surface reflectance products.
  2. MODIS CMA files are daily global aerosol optical thickness products.
******************************************************************************/
#include "combine_l8_aux_data.h"

/* Program will look for these SDSs in the CMG/CMA inputs */
#define N_SDS 5
char list_of_sds[N_SDS][50] = {
    "Coarse Resolution Ozone",
    "Coarse Resolution Granule Time", 
    "Coarse Resolution Water Vapor",
    "Coarse Resolution Air Temperature (2m)", 
    "Coarse Resolution AOT at 550 nm"};
#define OZONE 0
#define WV 2
#define AIR_TEMP_2M 3
   
/* Global variables */
bool global_yearday_is_set = false;
char global_yearday[10]; 

/******************************************************************************
MODULE:  combine_l8_anc_data

PURPOSE:  Reads the daily, global Aqua and Terra CMG and CMA files and "fuses"
them into a single output HDF file.  The application fills in the holes of the
Terra data with the Aqua data.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the inputs or writing the fused output
SUCCESS        Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
8/26/2014    Gail Schmidt     Conversion of the original code delivered by
                              Eric Vermote, NASA GSFC, for use within ESPA

NOTES:
******************************************************************************/
int main (int argc, char **argv)
{    
    bool verbose;              /* verbose flag for printing messages */
    char FUNC_NAME[] = "main"; /* function name */
    char errmsg[STR_SIZE];     /* error message */
    char dim0name[] = "YDim_MOD09CMG";   /* y dimension name */
    char dim1name[] = "XDim_MOD09CMG";   /* x dimension name */
    char *terra_cmg_file = NULL;  /* input Terra CMG file */
    char *aqua_cmg_file = NULL;   /* input Aqua CMG file */
    char *terra_cma_file = NULL;  /* input Terra CMA file */
    char *aqua_cma_file = NULL;   /* input Aqua CMA file */
    char *output_dir = NULL;      /* output directory for the auxiliary file */
    char tmpstr[STR_SIZE];        /* temporary string for creating file
                                     attributes */
    char outfilename[STR_SIZE];       /* name of the output HDF file */
    io_param terra_params[N_SDS];     /* array of Terra SDS parameters */
    io_param aqua_params[N_SDS];      /* array of Aqua SDS parameters */
    long pix;                /* current pixel location in the 1D array */
    int i, j;                /* looping variables */
    int nbits;               /* number of bits per pixel for this data array */
    int line;                /* current line in the CMG data array */
    int samp;                /* current sample in the line */
    int left, right;         /* pixel locations for interpolation */ 
    int n_pixels;            /* number of pixels in this 1D array */
    int n_bad;               /* number of bad/mismatches SDSs */
    int retval;              /* return status */
    int32 dims[2] = {IFILL, IFILL}; /* dimensions of desired CMG/CMA SDSs */
    int32 sd_out;            /* SD ID for the output file */
    int32 sds_id[N_SDS+1];   /* SDS IDs for the output file */
    int32 dimid;             /* dimension ID */
    int32 start[2];          /* starting location in each dimension */
    int32 where[N_SDS];      /* location of any missing SDSs */
    int8 *wherefrom = NULL;  /* array to identify where the pixel value was
                                pulled from - AQUA or TERRA */
    uint8 terra_pix;         /* terra pixel */
    uint8 aqua_pix;          /* aqua pixel */
    uint8 *tmask = NULL;     /* mask for the Terra pixel values */
    uint8 *amask = NULL;     /* mask for the Aqua pixel values */

    /* Read the command-line arguments */
    retval = get_args (argc, argv, &terra_cmg_file, &aqua_cmg_file,
        &terra_cma_file, &aqua_cma_file, &output_dir, &verbose);
    if (retval != SUCCESS)
    {   /* get_args already printed the error message */
        exit (ERROR);
    }

    /* Initialize the SDS information for the input files */
    for (i = 0; i < N_SDS; i++)
    {
        strcpy (terra_params[i].sdsname, "(missing SDS)");
        terra_params[i].sd_id = -1;
        terra_params[i].sds_id = -1;
        terra_params[i].data_type = -1;
        terra_params[i].data = NULL;

        strcpy (aqua_params[i].sdsname, "(missing SDS)");
        aqua_params[i].sd_id = -1;
        aqua_params[i].sds_id = -1;
        aqua_params[i].data_type = -1;
        aqua_params[i].data = NULL;
    }

    /* Read the input files */
    retval = parse_sds_info (terra_cmg_file, terra_params, aqua_params);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error parsing file: %s", terra_cmg_file);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
       
    retval = parse_sds_info (aqua_cmg_file, terra_params, aqua_params);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error parsing file: %s", aqua_cmg_file);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
       
    retval = parse_sds_info (terra_cma_file, terra_params, aqua_params);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error parsing file: %s", terra_cma_file);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
       
    retval = parse_sds_info (aqua_cma_file, terra_params, aqua_params);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error parsing file: %s", aqua_cma_file);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Check for missing SDSs, flag the missing SDSs in the 'where' array.
       Loop through each of the SDSs in the Terra and Aqua input structures.
       Keep track of where the Aqua params line up with the same Terra
       params. */
    for (i = 0; i < N_SDS; i++)
    {
        where[i] = -1;
        for (j = 0; j < N_SDS; j++)
        {
            if (terra_params[i].data_type == aqua_params[j].data_type &&
                !strcmp (terra_params[i].sdsname, aqua_params[j].sdsname))
                where[i] = j;
        }
    }

    /* Do we have any missing attributes between Terra and Aqua? */
    n_bad = 0;   
    for (i = 0; i < N_SDS; i++)
    {
        if (where[i] == -1)
            n_bad++;
    }
    if (n_bad > 0)
    {
        sprintf (errmsg, "Different sets of SDSs have been staged.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);

#ifdef DEBUG
        printf ("\nTerra:\n");
        for (i = 0; i < N_SDS; i++)
            printf (" NAME %s, TYPE %d\n", terra_params[i].sdsname,
                terra_params[i].data_type);
        printf ("\nAqua:\n");
        for (i = 0; i < N_SDS; i++)
            printf (" NAME %s, TYPE %d\n", aqua_params[i].sdsname,
                aqua_params[i].data_type);
#endif
    }

    /* Set the dimensions for the output product.  Use one of the input
       files. */
    dims[0] = terra_params[0].sds_dims[0];
    dims[1] = terra_params[0].sds_dims[1];

    /* Allocate memory for the data arrays, separate memory for each of the
       SDSs we are going to read and output */
    nbits = 0;
    n_pixels = dims[1] * dims[0];
    for (i = 0; i < N_SDS; i++)
    {
        if (terra_params[i].data_type == DFNT_INT16)
            nbits = sizeof (int16);
        else if (terra_params[i].data_type == DFNT_UINT16)
            nbits = sizeof (uint16);
        else if (terra_params[i].data_type == DFNT_INT8)
            nbits = sizeof (int8);
        else if (terra_params[i].data_type == DFNT_UINT8)
            nbits = sizeof (uint8);
        else
        {
            sprintf (errmsg, "Unsupported data type for SDS %s.  Only int16 "
                "uint16, int8, and uint8 are supported.",
                terra_params[i].sdsname);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        terra_params[i].data = calloc (n_pixels, nbits);
        if (terra_params[i].data == NULL)
        {
            sprintf (errmsg, "Allocating memory (%d bits) for Terra SDS: %s",
                nbits, terra_params[i].sdsname);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        aqua_params[i].data = calloc (n_pixels, nbits);
        if (aqua_params[i].data == NULL)
        {
            sprintf (errmsg, "Allocating memory (%d bits) for Aqua SDS: %s",
                nbits, aqua_params[i].sdsname);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }  /* end for i */

    /* Allocate memory for the wherefrom array, with bits of size int8 */
    wherefrom = calloc (n_pixels, sizeof (int8));
    if (wherefrom == NULL)
    {
        sprintf (errmsg, "Allocating memory for the wherefrom array");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Create the output file */
    make_outfile_name (global_yearday, output_dir, outfilename);
    if (verbose)
        printf ("Creating output auxiliary file: %s\n", outfilename);
    sd_out = SDstart (outfilename, DFACC_CREATE);
    if (sd_out == -1)
    {
        sprintf (errmsg, "Unable to create the output file %s", outfilename);
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Loop through the SDSs that we intend to read/write, and create an SDS
       in the output file for that SDS */
    for (i = 0; i < N_SDS; i++)
    {
        /* Create the SDS using information from the Terra file for this SDS */
        sds_id[i] = SDcreate (sd_out, terra_params[i].sdsname,
            terra_params[i].data_type, 2, dims);
        if (sds_id[i] == -1)
        {
            sprintf (errmsg, "Creating SDS %s in the output file",
                terra_params[i].sdsname);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Set the dimension names to the dimension ID */
        dimid = SDgetdimid (sds_id[i], 0);
        if (dimid != -1)
            SDsetdimname (dimid, dim0name);
        dimid = SDgetdimid (sds_id[i], 1);
        if (dimid != -1)
            SDsetdimname (dimid, dim1name); 
    }

    /* Create the wherefrom SDS to keep track of where each pixel came from */
    sds_id[i] = SDcreate (sd_out, "wherefrom", DFNT_INT8, 2, dims);
    if (sds_id[i] == -1)
    {
        sprintf (errmsg, "Unable to create the 'wherefrom' SDS in the output "
            "file");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

    /* Set the dimension names to the dimension ID */
    dimid = SDgetdimid (sds_id[i], 0);
    if (dimid != -1)
        SDsetdimname (dimid, dim0name);
    dimid = SDgetdimid (sds_id[i], 1);
    if (dimid != -1)
        SDsetdimname (dimid, dim1name); 

    /* Set the output file attributes */
    tmpstr[0] = '\0';
    for (i = 0; i < argc; i++)
        sprintf (tmpstr + strlen (tmpstr), " %s", argv[i]);
    SDsetattr (sd_out, "command", DFNT_CHAR, strlen (tmpstr), tmpstr);

    /* Start of processing the inputs .... */
    if (verbose)
        printf ("Reading each SDS ...\n");

    /* Read each SDS */
    start[0] = 0;
    start[1] = 0;
    for (i = 0; i < N_SDS; i++)
    {
        if (verbose)
            printf ("    %s\n", terra_params[i].sdsname);

        /* Read the Terra data for this SDS */
        retval = SDreaddata (terra_params[i].sds_id, start, NULL, dims,
            terra_params[i].data);
        if (retval == -1)
        {
            sprintf (errmsg, "Unable to read the SDS %s from the Terra file.",
                terra_params[i].sdsname);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }

        /* Read the Aqua data for this SDS */
        retval = SDreaddata (aqua_params[i].sds_id, start, NULL, dims,
            aqua_params[i].data);
        if (retval == -1)
        {
            sprintf (errmsg, "Unable to read the SDS %s from the Aqua file.",
                aqua_params[i].sdsname);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Use the Coarse Resolution Ozone SDS to determine if the pixel will come
       from Terra or Aqua.  This SDS is a uint8 data array. */
    if (verbose)
        printf ("Combining Aqua and Terra products for each SDS ...\n");
    tmask = (uint8 *) terra_params[OZONE].data;
    amask = (uint8 *) aqua_params[OZONE].data;
    for (i = 0; i < n_pixels; i++)
    {
        /* Initialize the masks */
        wherefrom[i] = UNSET;
        terra_pix = tmask[i];
        aqua_pix = amask[i];

        /* If the Terra pixel is not fill, then use Terra.  Otherwise if the
           Terra pixel is fill and the Aqua pixel isn't, then use Aqua.  In
           the latter case, the Aqua pixels for each SDS are copied over to the
           Terra array so that at the end the Terra array has all the output
           info. */
        if (terra_pix != 0)
        {  /* do nothing but set wherefrom */
            wherefrom[i] = TERRA;
        }
        else if (terra_pix == 0 && aqua_pix != 0)
        {  /* copy Aqua pixels over to Terra pixels for each SDS and set
              wherefrom */
            for (j = 0; j < N_SDS; j++)
            {
                copy_param (terra_params[j].data, aqua_params[j].data,
                    terra_params[j].data_type, i);
            }
            wherefrom[i] = AQUA;
        }
    }

    /* Interpolate water vapor, ozone, and temperature at 2m data.  But, only
       for lines 1000 to 2600, assuming CMGs (exclude the poles). */
    if (verbose)
        printf ("Interpolating combined products for WV, OZ, and "
            "AIR_TEMP_2M ...\n");
    if (dims[0] == 3600)
    {  /* then, yeah, we have a CMG */
        for (line = 1000; line < 2600; line++)
        {
            /* Get the current pixel location in the 1D array for this line */
            pix = (long) line * dims[1];

            /* Loop through the pixels in this line */
            left = right = -1;
            for (samp = 0; samp < dims[1]; samp++)
            {
                /* If the pixel is not fill then continue.  Recall that the
                   tmask now contains the final output data array, combined
                   from Terra and Aqua. */
                if (tmask[pix+samp] != 0)
                    continue;

                /* Find the left and right pixels in the line to use for
                   interpolation.  Basically need the non-fill pixels
                   surrounding the current pixel. */
                left = right = samp;
                while (tmask[pix+right] == 0) right++;
                samp = right;
                left--;

                /* Interpolate all the fill pixels between the left and right
                   non-fill pixels for the ozone, water vapor, and air temp
                   data */
                interpolate (DFNT_UINT8, terra_params[OZONE].data, pix, left,
                    right);
                interpolate (DFNT_UINT16, terra_params[WV].data, pix, left,
                    right);
                interpolate (DFNT_UINT16, terra_params[AIR_TEMP_2M].data, pix,
                    left, right);
            }
        }
    }  /* if dims[0] */

    /* Write each SDS to the output file */
    start[0] = 0;
    start[1] = 0;
    for (i = 0; i < N_SDS; i++)
    {
        retval = SDwritedata (sds_id[i], start, NULL, dims,
            terra_params[i].data); 
        if (retval == -1)
        {
            sprintf (errmsg, "Unable to write the %s SDS to the output file.",
                terra_params[i].sdsname);
            error_handler (true, FUNC_NAME, errmsg);
            exit (ERROR);
        }
    }

    /* Write the wherefrom SDS to the output file, and set the key attribute
       to provide information on these pixel values */
    retval = SDwritedata(sds_id[i], start, NULL, dims, wherefrom); 
    if (retval == -1)
    {
        sprintf (errmsg, "Unable to write the wherefrom SDS to the output "
            "file.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }
    strcpy (tmpstr, "0=none, 1=Terra, 2=Aqua"); 
    SDsetattr (sds_id[i], "key", DFNT_CHAR, strlen (tmpstr), tmpstr);

    /* Close and clean up */
    for (i = 0; i < N_SDS; i++)
    {
        SDendaccess (terra_params[i].sds_id);
        free (terra_params[i].data);
        SDend (terra_params[i].sd_id);

        SDendaccess (aqua_params[i].sds_id);
        free (aqua_params[i].data);
        SDend (aqua_params[i].sd_id);

        SDendaccess (sds_id[i]);
    }   
    SDendaccess (sds_id[i]);
    SDend (sd_out);
    free (terra_cmg_file);
    free (aqua_cmg_file);
    free (terra_cma_file);
    free (aqua_cma_file);
    free (output_dir);
    free (wherefrom);

    /* Successful completion */
    exit (SUCCESS);
}


/******************************************************************************
MODULE:  interpolate

PURPOSE:  Interpolates all fill pixels between the left and right pixels.

RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
8/28/2014    Gail Schmidt     Conversion of the original code delivered by
                              Eric Vermote, NASA GSFC, for use within ESPA

NOTES:
  1. Only supports uint8 and uint16.
******************************************************************************/
void interpolate
(
    int32 data_type,     /* I: data type of the data array */
    void *data,          /* I: data array */
    long lineoffset,     /* I: pixel location for the start of this line */
    int left,            /* I: location in the line of the left pixel */
    int right            /* I: location in the line of the right pixel */
)
{
    uint8 *ui8x = NULL;     /* uint8 pointer */
    uint16 *ui16x = NULL;   /* uint16 pointer */
    int i;                  /* looping variable */
    int diff;               /* distance between the left and right pixels */
    float slope;            /* slope for this pixel */

    /* Determine the distance between the left and right pixels */
    diff = right - left;

    /* Handle the interpolation between the pixels based on the data type */
    if (data_type == DFNT_UINT8)
    {
        ui8x = (uint8 *)data;
        if (ui8x[lineoffset+right] > ui8x[lineoffset+left])
        {
            slope = ((float) ui8x[lineoffset+right] -
                     (float) ui8x[lineoffset+left]) / (float) (diff);
            for (i = 0; i < diff; i++)
            {
                ui8x[lineoffset+i+left] = (uint8)
                    ((float) ui8x[lineoffset+left] + (slope * i));
            }
        }
        else
        {
            slope = ((float) ui8x[lineoffset+left] - 
                     (float) ui8x[lineoffset+right]) / (float) (diff);
            for (i = 0; i < diff; i++)
            {
                ui8x[lineoffset+i+left] = (uint8)
                    ((float) ui8x[lineoffset+left] - (slope * i));
            }
        }
    }
    else if (data_type == DFNT_UINT16)
    {
        ui16x = (uint16 *)data;
        if (ui16x[lineoffset+right] > ui16x[lineoffset+left])
        {
            slope = ((float) ui16x[lineoffset+right] - 
                     (float) ui16x[lineoffset+left]) / (float) (diff);
            for (i = 0; i < diff; i++)
            {
                ui16x[lineoffset+i+left] = (uint16)
                    ((float) ui16x[lineoffset+left] + (slope * i));
            }
        }
        else
        {
            slope = ((float) ui16x[lineoffset+left] -
                     (float) ui16x[lineoffset+right]) / (float)(diff);
            for (i = 0; i < diff; i++)
            {
                ui16x[lineoffset+i+left] = (uint16)
                    ((float) ui16x[lineoffset+left] - (slope * i));
            }
        }
    }

    return;
}


/******************************************************************************
MODULE:  make_outfile_name

PURPOSE:  Creates the output filename for the auxiliary products, using the
  input year-DOY string and the source directory.

RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
8/28/2014    Gail Schmidt     Conversion of the original code delivered by
                              Eric Vermote, NASA GSFC, for use within ESPA

NOTES:
******************************************************************************/
void make_outfile_name
(
    char *yearday_str,      /* I: string containing the year and DOY */
    char *output_dir,       /* I: output directory for the auxiliary prods */
    char outfile[STR_SIZE]  /* O: output filename for the auxiliary products */
)
{
    sprintf (outfile, "%s/L8ANC%s.hdf_fused", output_dir, yearday_str);
    return;
}


/******************************************************************************
MODULE:  copy_param

PURPOSE:  Creates the output filename for the auxiliary products, using the
  input year-DOY string and the source directory.

RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
8/28/2014    Gail Schmidt     Conversion of the original code delivered by
                              Eric Vermote, NASA GSFC, for use within ESPA

NOTES:
******************************************************************************/
void copy_param
(
    void *dest,        /* O: destination data array */
    void *source,      /* I: source data array */
    int32 data_type,   /* I: data type */
    int32 offset       /* I: pixel in source/dest data arrays to be copied */
)
{
    int8 *i8dest = NULL;
    int8 *i8source = NULL;
    uint8 *ui8dest = NULL;
    uint8 *ui8source = NULL;
    int16 *i16dest = NULL;
    int16 *i16source = NULL;
    uint16 *ui16dest = NULL;
    uint16 *ui16source = NULL;

    /* Copy the source pixel to the destination pixel, based on the data
       type */
    switch (data_type)
    {
        case DFNT_INT8:
            i8dest = (int8 *)dest;
            i8source = (int8 *)source;
            i8dest[offset] = i8source[offset];
            break;

        case DFNT_UINT8:
            ui8dest = (uint8 *)dest;
            ui8source = (uint8 *)source;
            ui8dest[offset] = ui8source[offset];
            break;

        case DFNT_INT16:
            i16dest = (int16 *)dest;
            i16source = (int16 *)source;
            i16dest[offset] = i16source[offset];
            break;

        case DFNT_UINT16:
            ui16dest = (uint16 *)dest;
            ui16source = (uint16 *)source;
            ui16dest[offset] = ui16source[offset];

        default:
            break;
    }

    return;
}


/******************************************************************************
MODULE:  parse_sds_info

PURPOSE:  Reads the daily, global Aqua and Terra CMG and CMA files.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the inputs or writing the fused output
SUCCESS        Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
8/26/2014    Gail Schmidt     Conversion of the original code delivered by
                              Eric Vermote, NASA GSFC, for use within ESPA

NOTES:
  1. If memory for the Terra_parameters or Aqua_parameters structure has not
     been allocated, this function will allocate memory for an array of
     N_SDS io_param structures.  Otherwise the passed in array will be
     utilized.
******************************************************************************/
int parse_sds_info
(
    char *filename,            /* I: Aqua/Terra file to be read */
    io_param Terra_params[],   /* O: array of structs for Terra params */
    io_param Aqua_params[]     /* O: array of structs for Aqua params */
)
{
    char FUNC_NAME[] = "parse_sds_info"; /* function name */
    char errmsg[STR_SIZE];  /* error message */
    int retval;             /* return status */
    int i, j;               /* looping variables */
    int sd_id;              /* file ID for the HDF file */
    int sds_id;             /* ID for the current SDS */
    int nsds;               /* number of SDSs in the file */
    int nattr;              /* number of attributes for this file */
    int rank;               /* rank of the dimensions in this SDS */
    int sat;                /* satellite type - TERRA, AQUA */
    int localdims[2];       /* stored version of the dimensions */
    int sds_dims[2];        /* SDS dimension sizes */
    int data_type;          /* value representing the data type of this SDS */
    int n_val = 0;          /* number of values read in the string */
    char sds_name[STR_SIZE]; /* name of the SDS at the specified index */
    char lgid[STR_SIZE];    /* local granule ID */
    char product_type[20];  /* MODIS product type */
    char yearday[10];       /* year/day string */

    /* Open the input file for reading the metadata and attributes */
    sd_id = SDstart (filename, DFACC_RDONLY);
    if (sd_id == -1)
    {
        sprintf (errmsg, "Error reading file: %s", filename);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    
    retval = SDfileinfo (sd_id, &nsds, &nattr);
    if (retval == -1)
    {
        sprintf (errmsg, "SDfileinfo error reading file: %s", filename);
        error_handler (true, FUNC_NAME, errmsg);
        SDend (sd_id);
        return (ERROR);
    }

    if (nsds < 1)
    {
        sprintf (errmsg, "File %s contains no SDSs", filename);
        error_handler (true, FUNC_NAME, errmsg);
        SDend (sd_id);
        return (ERROR);
    }

    /* Get the Local Granule ID of the file, and parse it into what we need */
    for (i = 0; i < STR_SIZE; i++)
        lgid[i] = '\0';
    if (metareader (sd_id, "COREMETADATA", "LOCALGRANULEID", &n_val, lgid) !=
        SUCCESS)
    {
        sprintf (errmsg, "Error parsing the core metadata for LOCALGRANULEID");
        error_handler (true, FUNC_NAME, errmsg);
        SDend (sd_id);
        return (ERROR);
    }

    /* Validate this Local Granule ID is what we expect.  It should only be
       40 characters long (41 with the end of string).  If it's valid, then
       parse the granule ID for date and day of year.
       Example - MOD09CMA.A2014133.006.2014135103800.hdf */
    if (strlen (lgid) != 41)
    {
        sprintf (errmsg, "File %s has unrecognizable LocalGranuleID: %s",
            filename, lgid);
        error_handler (true, FUNC_NAME, errmsg);
        SDend (sd_id);
        return (ERROR);
    }
    parse_lgid (lgid, product_type, yearday);

    /* Set this global variable */
    if (!global_yearday_is_set)
    {
        global_yearday_is_set = true;
        strcpy(global_yearday, yearday); 
    }
    else
    {
        if (strcmp (global_yearday, yearday))
        {
            sprintf (errmsg, "File %s has date of %s, which differs from the "
                "other files which have been parsed.", filename, yearday);
            error_handler (true, FUNC_NAME, errmsg);
            SDend (sd_id);
            return (ERROR);
        }
    }

    /* Validate the product type in the Local Granule ID */
    if (!strcmp (product_type, "MOD09CMA") ||
        !strcmp (product_type, "MOD09CMG"))
        sat = TERRA;
    else if (!strcmp (product_type, "MYD09CMG") ||
             !strcmp (product_type, "MYD09CMA"))
        sat = AQUA;
    else
    {
        sprintf (errmsg, "File %s has unrecognizable product type: %s.  Expect "
            "M[OY]D09CM[GA].", filename, product_type);
        error_handler (true, FUNC_NAME, errmsg);
        SDend (sd_id);
        return (ERROR);
    }

    /* Now check out the SDSs and their associated information */
    localdims[0] = IFILL;
    localdims[1] = IFILL;
    for (i = 0; i < nsds; i++)
    {
        sds_id = SDselect (sd_id, i);
        if (sds_id == -1)
        {
            sprintf (errmsg, "SDselect error for SDS %d", i);
            error_handler (true, FUNC_NAME, errmsg);
            SDend (sd_id);
            return (ERROR);
        }
      
        /* CMG and CMA files should all have SDSs of the same rank (2) and
           dimensions (3600 by 7200).  If some file has as SDS with a different
           dimension, report it as a warning. */
        retval = SDgetinfo (sds_id, sds_name, &rank, sds_dims, &data_type,
            &nattr);   
        if (retval == -1)
        {
            sprintf (errmsg, "SDgetinfo error for SDS %d", i);
            error_handler (true, FUNC_NAME, errmsg);
            SDend (sd_id);
            return (ERROR);
        }

        if (rank != 2)
        {
            sprintf (errmsg, "SDS %d has unanticipated rank of %d, "
                "skipping ...", i, rank);
            error_handler (false, FUNC_NAME, errmsg);
            continue;
        }

        if (localdims[0] == IFILL && localdims[1] == IFILL)
        {
            localdims[0] = sds_dims[0];
            localdims[1] = sds_dims[1];
        }
        else
        {
            if (localdims[0] != sds_dims[0])
            {
                 sprintf (errmsg, "SDS has unanticipated x-dimension size of "
                     "%d, skipping ...", sds_dims[0]);
                 error_handler (false, FUNC_NAME, errmsg);
                 continue;
            }

            if (localdims[1] != sds_dims[1])
            {
                 sprintf (errmsg, "SDS has unanticipated y-dimension size of "
                     "%d, skipping ...", sds_dims[1]);
                 error_handler (false, FUNC_NAME, errmsg);
                 continue;
            }
        }
    
        /* Check against names of SDSs we need */
        for (j = 0; j < N_SDS; j++)
        {
            if (!strcmp (sds_name, list_of_sds[j]))
            {  /* keep this SDS info */
                if (sat == TERRA)
                {
                    Terra_params[j].sd_id = sd_id;
                    Terra_params[j].sds_id = sds_id;
                    Terra_params[j].data_type = data_type;
                    Terra_params[j].data = (void *)NULL;    
                    strcpy (Terra_params[j].sdsname, sds_name);
                    Terra_params[j].sds_dims[0] = sds_dims[0];
                    Terra_params[j].sds_dims[1] = sds_dims[1];
                }
                else if (sat == AQUA)
                {
                    Aqua_params[j].sd_id = sd_id;
                    Aqua_params[j].sds_id = sds_id;
                    Aqua_params[j].data_type = data_type;
                    Aqua_params[j].data = (void *)NULL;    
                    strcpy (Aqua_params[j].sdsname, sds_name);
                    Aqua_params[j].sds_dims[0] = sds_dims[0];
                    Aqua_params[j].sds_dims[1] = sds_dims[1];
                }
            }
        }  /* for j */
    }  /* for i */

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  parse_lgid

PURPOSE:  Parses the Local Granule ID to obtain key metadata information.

RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
8/26/2014    Gail Schmidt     Conversion of the original code delivered by
                              Eric Vermote, NASA GSFC, for use within ESPA

NOTES:
******************************************************************************/
void parse_lgid
(
    char lgid[],          /* I: local granule ID */
    char product_type[],  /* O: MODIS product type */
    char yearday[]        /* O: year/day string */
)
{
    int32 i;             /* looping variable */
    int32 n;             /* number of characters in the lgid string */
    int32 x = 0;         /* current character location in the output string */
    int32 n_token = 0;   /* number of "tokens" in the string */
    char c;              /* value of the current character in the string */
    char mystr[100];     /* current string */

    /* Parse the local granule ID string.  Start at char 1 (0-based) to skip
       the "\"" character that metareader() always returns.  Split the granule
       ID into "tokens", split by the '.'s in the string.
       Example - MOD09CMA.A2014133.006.2014135103800.hdf */
    n = strlen (lgid);
    for (i = 1; i < n; i++)
    {
        c = lgid[i];
        mystr[x] = c;

        /* Is this the start of a new "token"? */
        if (c == '.')
        {
            /* NULL terminate the string, then copy it to the desired variable
               depending on whether or not it's the first or second token in
               the local granule ID */
            mystr[x] = '\0';
            if (n_token == 0)
                strcpy (product_type, mystr);
            else if (n_token == 1)
            {
                /* Skip the first character in the year/day token */
                strncpy (yearday, &mystr[1], 7);
                yearday[7] = '\0';
            }

            /* Increment the token count and reset the location of the
               next character in the current string */
            n_token++;
            x = 0;
        }
        else x++;
    }

    return;
}


/******************************************************************************
MODULE:  metareader

PURPOSE:  Parses the metadata to obtain the desired metadata attribute.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the metadata
SUCCESS        Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
8/26/2014    Gail Schmidt     Conversion of the original code delivered by
                              Eric Vermote, NASA GSFC, for use within ESPA

NOTES:
******************************************************************************/
int metareader
(
    int sd_id,               /* I: file ID */
    char *type_of_meta,      /* I: which metadata will be read */
    char *metastring,        /* I: metadata variable to be found */
    int32 *count,            /* O: count of characters in output string */
    char *data               /* O: string returned for metastring */
)
{
    char errmsg[STR_SIZE];       /* error message */
    char FUNC_NAME[] = "metareader"; /* function name */
    char attr_name[STR_SIZE];    /* name of the attribute */
    char *charattr = NULL;       /* character string for the attribute */
    char line[MAXLENGTH2];       /* line of metadata */
    char objs[MAXLENGTH];        /* current metadata string */
    char lhs[MAXLENGTH];         /* left-hand side of the metadata line */
    char rhs[MAXLENGTH2];        /* right-hand side of the metadata line */
    int retval;             /* return status */
    int nsds;               /* number of SDSs in the file */
    int nattr;              /* number of attributes for this file */
    int i, ii;              /* looping variables for parsing the metadata
                               line */
    int j;                  /* looping variable for the attributes */
    int data_type;          /* data type of the attribute */
    int nvals;              /* number of values in the attribute */
    int n_obj;              /* number of OBJECT strings we've counted */
    int n_val;              /* number of values for the current attribute */
    int start;              /* location of the character string to start
                               reading lines */
    int obj_offset[10];     /* offset of location of OBJECT token */
    bool within;            /* are we within the attribute string */
    bool wasjustobject;     /* was the left-hand side of the string the
                               OBJECT token? */

    /* Get the information for the current SDS */
    retval = SDfileinfo (sd_id, &nsds, &nattr);
    if (retval == -1)
    {
        sprintf (errmsg, "SDgetinfo error for SDS %d", sd_id);
        error_handler (true, FUNC_NAME, errmsg);
        SDend (sd_id);
        return (ERROR);
    }

    /* Loop through the attributes */
    for (j = 0; j < nattr; j++)
    {
        /* Get the information for the current attribute */
        retval = SDattrinfo (sd_id, j, attr_name, &data_type, &nvals);
        if (retval == -1)
        {
            sprintf (errmsg, "SDattrinfo error for SDS %d", sd_id);
            error_handler (true, FUNC_NAME, errmsg);
            SDend (sd_id);
            return (ERROR);
        }
 
        /* Make the attribute name all uppercase, then check to see if it's
           the attribute we are looking for */
        start = 0;
        for (i = 0; i < strlen (attr_name); i++)
            attr_name[i] = toupper (attr_name[i]);  
        attr_name[i] = '\0';
        if (strstr (attr_name, type_of_meta))
        {
            /* Allocate memory for reading the attribute; assume it's a
               character type since those are the attributes we are reading. */
            charattr = malloc ((nvals+1) * sizeof (char));
            if (charattr == NULL)
            {
                sprintf (errmsg, "Unable to allocate memory for reading "
                    "attribute: %s", attr_name);
                error_handler (true, FUNC_NAME, errmsg);
                SDend (sd_id);
                return (ERROR);
            }

            /* Read the attribute */
            retval = SDreadattr (sd_id, j, charattr);
            if (retval == -1)
            {
                sprintf (errmsg, "SDreadattr error for attribute %s", charattr);
                error_handler (true, FUNC_NAME, errmsg);
                SDend (sd_id);
                return (ERROR);
            }

            /* Parse the attribute */
            n_obj = 0;
            wasjustobject = false;
            objs[0] = '\0';
            do
            { 
                line[0] = '\0';
                get_a_line (charattr, &start, line);

                /* Get rid of whitespace, if not inside "..." characters */
                within = false;
                for (i = 0, ii = 0; i < strlen (line); i++)
                {
                    /* Is this the start of the attribute line (starts and
                       ends with quotes) */
                    if (line[i] == '"')
                    {
                        if (!within)
                            within = true;
                        else if (within)
                            within = false;
                    }

                    /* Copy the string, but skip the white spaces if we are
                       outside the quotes */
                    if (within)
                        line[ii++] = line[i];
                    else
                    {
                        if (line[i] != ' ')
                            line[ii++] = line[i];
                    }
                }
                line[ii] = '\0';

                /* Get left-hand-side and right-hand-side of "equation", which
                   is split by the '=' sign. */
                lhs[0] = '\0';
                rhs[0] = '\0';

                /* If there is an equal sign in the string, then copy the
                   characters up to the equal sign into the left-hand side */
                if (strchr (line, '='))
                {
                    for (i = 0, ii = 0; i < (strchr (line, '=') - line); i++)
                        lhs[ii++] = line[i];
                    lhs[ii] = '\0';

                    /* Now grab the right-hand side, but skip over the equal
                       sign first */
                    i++;
                    for (ii = 0; i < strlen (line) - 1; i++, ii++)
                        rhs[ii] = line[i];
                    rhs[ii] = '\0';
                }

                /* Is this the OBJECT token */
                if (!strcmp (lhs, "OBJECT"))
                {
                    wasjustobject = true;
                    obj_offset[n_obj++] = strlen (objs);
                    strcat (objs, rhs);
                }

                /* Is this the GROUP token */
                if (!strcmp (lhs, "GROUP"))
                    wasjustobject = false;

                /* If this is the CLASS token within the OBJECT token then
                   copy the right-hand side of the string to our overall
                   string */
                if (!strcmp (lhs, "CLASS") && wasjustobject)
                {
                    for (i = 0; i < strlen(rhs) - 1; i++)
                        rhs[i] = rhs[i+1];
                    rhs[i-1] = '\0';
                    strcat (objs, rhs);
                }

                /* If this is the END_OBJECT token, then decrement the count
                   of open objects and clear the object string */
                if (!strcmp (lhs, "END_OBJECT"))
                    objs[obj_offset[--n_obj]] = '\0';

                /* If this is the NUM_VAL object, then grab the number of
                   values */
                if (!strcmp (lhs, "NUM_VAL"))
                    n_val = atoi (rhs);
          
                /* If this is the VALUE token, and it's the metadata token
                   the user has requested. */
                if (!strcmp (lhs, "VALUE"))
                {
                    if (!strcmp (objs, metastring))
                    {
                        strcpy (data, rhs);
                        *count = n_val;
                    }
                }
            } while (line[0] != '\0');

            free (charattr);  
        }  /* If this is the desired metadata attribute */
    }  /* for j */

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_a_line

PURPOSE:  Reads a line from the input text string

RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
8/28/2014    Gail Schmidt     Conversion of the original code delivered by
                              Eric Vermote, NASA GSFC, for use within ESPA

NOTES:
******************************************************************************/
void get_a_line
(
    char text[MAXLENGTH2],  /* I: text string to be read from */
    int *start,             /* I/O: location where to start reading the line;
                                    updated for where to start reading the next
                                    line after the current line is read */
    char *line              /* O: line that was read from the text string */
)
{
    int i;               /* current location in the output line */
    int where = 0;       /* where to start reading in the input text string */
    bool getout;         /* are we done reading - at end of the string or at
                            the end of the line */

    /* Initialize the line to an empty line */
    where = *start;
    line[0] = '\0';

    /* If this is an empty line, then just return */
    if (text[where] == '\0')
        return;

    /* Loop through the input text string, reading each character, and stopping
       when and end of string char, end of line char, or the end of the
       actual string is found */
    i = 0;
    getout = false;
    while (!getout)
    {
        /* Are we at the end of the string or have we found an end of string
           or end of line character? */
        if (text[where] == '\0' || text[where] == '\n' || i >= MAXLENGTH2)
            getout = true;

        /* Copy the current character to the output line and keep track of
           where we are in the input and output strings */
        line[i++] = text[where++];
    }; 
    line[i] = '\0';

    /* Return the location where the next line should start in the input
       text string */
    *start = where;
    return;
}


/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date        Programmer       Reason
---------   ---------------  -------------------------------------
8/26/2014   Gail Schmidt     Original Development

NOTES:
******************************************************************************/
void usage ()
{
    printf ("combine_l8_aux_data derives auxiliary data from the Aqua/Terra "
            "CMG and CMA files (for the same day).  The output is a \"fused\" "
            "HDF file containing the desired SDSs and data for that day. "
            "The application fills in the holes of the Terra data with the "
            "Aqua data.\n\n");
    printf ("usage: combine_l8_aux_data "
            "--terra_cmg=input_terra_cmg_filename "
            "--aqua_cmg=input_aqua_cmg_filename "
            "--terra_cma=input_terra_cma_filename "
            "--aqua_cma=input_aqua_cma_filename "
            "--output_dir=output_directory "
            "[--verbose]\n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -terra_cmg: name of the input Terra CMG file to be "
            "processed\n");
    printf ("    -aqua_cmg: name of the input Terra CMG file to be "
            "processed\n");
    printf ("    -terra_cma: name of the input Terra CMA file to be "
            "processed\n");
    printf ("    -aqua_cma: name of the input Terra CMA file to be "
            "processed\n");
    printf ("    -output_dir: name of the output directory for the combined "
            "auxiliary file to be written\n");

    printf ("\nwhere the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed? (default "
            "is false)\n");

    printf ("\ncombine_l8_aux_data --help will print the usage statement\n");
}
