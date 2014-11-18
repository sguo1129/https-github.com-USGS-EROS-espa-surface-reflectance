/*****************************************************************************
FILE: input.c
  
PURPOSE: Contains functions for handling of the input data files for this
application.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
6/23/2014    Gail Schmidt     Original development

NOTES:
*****************************************************************************/

#include "input.h"

/******************************************************************************
MODULE:  open_input

PURPOSE:  Sets up the input data structure, opens the input reflectance file
for read access, allocates space, and stores some of the metadata for later
reference.

RETURN VALUE:
Type = Input_t*
Value      Description
-----      -----------
NULL       Error occurred opening or reading the file
non-NULL   Successful completion

HISTORY:
Date         Programmer       Reason
----------   ---------------  -------------------------------------
6/20/2014    Gail Schmidt     Original Development
11/17/2014   Gail Schmidt     Modified to support OLI-only scenes

NOTES:
  1. This routine opens the input L8 files.  It also allocates memory for
     pointers in the input structure.  It is up to the caller to use
     close_input and free_input to close the files and free up the memory when
     done using the input data structure.
******************************************************************************/
Input_t *open_input
(
    Espa_internal_meta_t *metadata      /* I: input metadata */
)
{
    char FUNC_NAME[] = "open_input";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    Input_t *this = NULL;     /* input data structure to be initialized,
                                 populated, and returned to the caller */
    int ib;                   /* loop counter for bands */

    /* Create the Input data structure */
    this = malloc (sizeof (Input_t));
    if (this == NULL) 
    {
        strcpy (errmsg, "Error allocating memory for Input data structure");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Initialize and get input from metadata file */
    if (get_xml_input (metadata, this) != SUCCESS)
    {
        free (this);
        this = NULL;
        strcpy (errmsg, "Error getting input information from the metadata "
            "file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (NULL);
    }

    /* Open files for access */
    for (ib = 0; ib < this->nband; ib++)
    {
        this->fp_bin[ib] = open_raw_binary (this->file_name[ib], "rb");
        if (this->fp_bin[ib] == NULL)
        {
            free_input (this);
            sprintf (errmsg, "Opening reflectance raw binary file: %s",
                this->file_name[ib]);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }
        this->open[ib] = true;
    }

    for (ib = 0; ib < this->nband_th; ib++)
    {  /* NOTE: nband_th will be 0 for OLI-only scenes */
        this->fp_bin_th[ib] = open_raw_binary (this->file_name_th[ib], "rb");
        if (this->fp_bin_th[ib] == NULL)
        {
            free_input (this);
            sprintf (errmsg, "Opening thermal raw binary file: %s",
                this->file_name_th[ib]);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }
        this->open_th[ib] = true;
    }

    for (ib = 0; ib < this->nband_pan; ib++)
    {
        this->fp_bin_pan[ib] = open_raw_binary (this->file_name_pan[ib], "rb");
        if (this->fp_bin_pan[ib] == NULL)
        {
            free_input (this);
            sprintf (errmsg, "Opening pan raw binary file: %s",
                this->file_name_pan[ib]);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }
        this->open_pan[ib] = true;
    }

    for (ib = 0; ib < this->nband_qa; ib++)
    {
        this->fp_bin_qa[ib] = open_raw_binary (this->file_name_qa[ib], "rb");
        if (this->fp_bin_qa[ib] == NULL)
        {
            free_input (this);
            sprintf (errmsg, "Opening QA raw binary file: %s",
                this->file_name_qa[ib]);
            error_handler (true, FUNC_NAME, errmsg);
            return (NULL);
        }
        this->open_qa[ib] = true;
    }

    return this;
}


/******************************************************************************
MODULE:  close_input

PURPOSE:  Ends SDS access and closes the input file.

RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
5/19/2014    Gail Schmidt     Original Development (based on input routines
                              from the spectral indices application)

NOTES:
******************************************************************************/
void close_input
(
    Input_t *this    /* I: pointer to input data structure */
)
{
    int ib;      /* loop counter for bands */
  
    /* Close the reflectance files */
    for (ib = 0; ib < this->nband; ib++)
    {
        if (this->open[ib])
        {
            close_raw_binary (this->fp_bin[ib]);
            this->open[ib] = false;
        }
    }

    /* Close the thermal files */
    for (ib = 0; ib < this->nband_th; ib++)
    {
        if (this->open_th[ib])
        {
            close_raw_binary (this->fp_bin_th[ib]);
            this->open_th[ib] = false;
        }
    }

    /* Close the pan files */
    for (ib = 0; ib < this->nband_pan; ib++)
    {
        if (this->open_pan[ib])
        {
            close_raw_binary (this->fp_bin_pan[ib]);
            this->open_pan[ib] = false;
        }
    }

    /* Close the QA files */
    for (ib = 0; ib < this->nband_qa; ib++)
    {
        if (this->open_qa[ib])
        {
            close_raw_binary (this->fp_bin_qa[ib]);
            this->open_qa[ib] = false;
        }
    }
}


/******************************************************************************
MODULE:  free_input

PURPOSE:  Frees memory in the input data structure.

RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/24/2014    Gail Schmidt     Original Development

NOTES:
******************************************************************************/
void free_input
(
    Input_t *this    /* I: pointer to input data structure */
)
{
    char FUNC_NAME[] = "free_input";   /* function name */
    char errmsg[STR_SIZE];             /* error message */
    int ib;                            /* loop counter for bands */
   
    if (this != NULL)
    {
        if (this->open[0] || this->open_th[0] || this->open_pan[0]) 
        {
            strcpy (errmsg, "Freeing input data structure, but reflectance, "
                "thermal, and/or pan file(s) is/are still open. Use "
                "close_input to close the file");
            error_handler (false, FUNC_NAME, errmsg);
        }
  
        /* Free image band pointers */
        for (ib = 0; ib < this->nband; ib++)
            free (this->file_name[ib]);
        for (ib = 0; ib < this->nband_th; ib++)
            free (this->file_name_th[ib]);
        for (ib = 0; ib < this->nband_pan; ib++)
            free (this->file_name_pan[ib]);
        for (ib = 0; ib < this->nband_qa; ib++)
            free (this->file_name_qa[ib]);

        /* Free the data structure */
        free (this);
    } /* end if */
}


/******************************************************************************
MODULE:  get_input_refl_lines

PURPOSE:  Reads the reflectance data for the current refl band and lines, and
populates the output buffer.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading data for this band
SUCCESS    Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/24/2014    Gail Schmidt     Original Development

NOTES:
  1. The Input_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_input to do that.
******************************************************************************/
int get_input_refl_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current refl band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint16 *out_arr  /* O: output array to populate */
)
{
    char FUNC_NAME[] = "get_input_refl_line";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->open[iband])
    {
        strcpy (errmsg, "Reflectance band has not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iband < 0 || iband >= this->nband)
    {
        strcpy (errmsg, "Invalid reflectance band number for the input date");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->size.nlines)
    {
        strcpy (errmsg, "Invalid line number for reflectance band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Read the data, but first seek to the correct line */
    loc = (long) iline * this->size.nsamps * sizeof (uint16);
    if (fseek (this->fp_bin[iband], loc, SEEK_SET))
    {
        strcpy (errmsg, "Seeking to the current line in the input file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_bin[iband], nlines, this->size.nsamps,
        sizeof (uint16), out_arr) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from reflectance band %d starting "
            "at line %d", nlines, iband, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_input_th_lines

PURPOSE:  Reads the thermal data for the current thermal band and lines, and
populates the output buffer.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading data for this band
SUCCESS    Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/24/2014    Gail Schmidt     Original Development

NOTES:
  1. The Input_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_input to do that.
******************************************************************************/
int get_input_th_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current thermal band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint16 *out_arr  /* O: output array to populate */
)
{
    char FUNC_NAME[] = "get_input_th_line";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->open_th[iband])
    {
        strcpy (errmsg, "Thermal band has not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iband < 0 || iband >= this->nband_th)
    {
        strcpy (errmsg, "Invalid thermal band number for the input data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->size_th.nlines)
    {
        strcpy (errmsg, "Invalid line number for thermal band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Read the data, but first seek to the correct line */
    loc = (long) iline * this->size_th.nsamps * sizeof (uint16);
    if (fseek (this->fp_bin_th[iband], loc, SEEK_SET))
    {
        strcpy (errmsg, "Seeking to the current line in the input file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_bin_th[iband], nlines, this->size_th.nsamps,
        sizeof (uint16), out_arr) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from thermal band %d starting at "
            "line %d", nlines, iband, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_input_pan_lines

PURPOSE:  Reads the pan data for the current pan band and lines, and populates
the output buffer.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading data for this band
SUCCESS    Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/24/2014    Gail Schmidt     Original Development

NOTES:
  1. The Input_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_input to do that.
******************************************************************************/
int get_input_pan_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current pan band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint16 *out_arr  /* O: output array to populate */
)
{
    char FUNC_NAME[] = "get_input_pan_line";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->open_pan[iband])
    {
        strcpy (errmsg, "Pan band has not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iband < 0 || iband >= this->nband_pan)
    {
        strcpy (errmsg, "Invalid pan band number for the input data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->size_pan.nlines)
    {
        strcpy (errmsg, "Invalid line number for pan band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Read the data, but first seek to the correct line */
    loc = (long) iline * this->size_pan.nsamps * sizeof (uint16);
    if (fseek (this->fp_bin_pan[iband], loc, SEEK_SET))
    {
        strcpy (errmsg, "Seeking to the current line in the input file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_bin_pan[iband], nlines, this->size_pan.nsamps,
        sizeof (uint16), out_arr) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from pan band %d starting at "
            "line %d", nlines, iband, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


/******************************************************************************
MODULE:  get_input_qa_lines

PURPOSE:  Reads the QA data for the current QA band and lines, and populates
the output buffer.

RETURN VALUE:
Type = int
Value      Description
-----      -----------
ERROR      Error occurred reading data for this band
SUCCESS    Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/24/2014    Gail Schmidt     Original Development

NOTES:
  1. The Input_t data structure needs to be populated and memory allocated
     before calling this routine.  Use open_input to do that.
******************************************************************************/
int get_input_qa_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current QA band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint16 *out_arr  /* O: output array to populate */
)
{
    char FUNC_NAME[] = "get_input_qa_line";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    long loc;                 /* current location in the input file */
  
    /* Check the parameters */
    if (this == NULL) 
    {
        strcpy (errmsg, "Input structure has not been opened/initialized");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (!this->open_qa[iband])
    {
        strcpy (errmsg, "QA band has not been opened");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iband < 0 || iband >= this->nband_qa)
    {
        strcpy (errmsg, "Invalid QA band number for the input data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    if (iline < 0 || iline >= this->size_qa.nlines)
    {
        strcpy (errmsg, "Invalid line number for QA band");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    /* Read the data, but first seek to the correct line */
    loc = (long) iline * this->size_qa.nsamps * sizeof (uint16);
    if (fseek (this->fp_bin_qa[iband], loc, SEEK_SET))
    {
        strcpy (errmsg, "Seeking to the current line in the input file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (read_raw_binary (this->fp_bin_qa[iband], nlines, this->size_qa.nsamps,
        sizeof (uint16), out_arr) != SUCCESS)
    {
        sprintf (errmsg, "Reading %d lines from QA band %d starting at "
            "line %d", nlines, iband, iline);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
  
    return (SUCCESS);
}


#define DATE_STRING_LEN (50)
#define TIME_STRING_LEN (50)

/******************************************************************************
MODULE:  get_xml_input

PURPOSE:  Pulls input values from the XML structure.

RETURN VALUE:
Type = int
Value        Description
-------      -----------
ERROR        Error occurred opening or reading the file
SUCCESS      Successful completion

HISTORY:
Date         Programmer       Reason
----------   ---------------  -------------------------------------
6/20/2014    Gail Schmidt     Original Development
11/17/2014   Gail Schmidt     Modified to support OLI-only scenes

NOTES:
******************************************************************************/
int get_xml_input
(
    Espa_internal_meta_t *metadata,  /* I: XML metadata */
    Input_t *this                    /* O: data structure for the input file */
)
{
    char FUNC_NAME[] = "get_xml_input";   /* function name */
    char errmsg[STR_SIZE];    /* error message */
    int ib;                   /* looping variable for bands */
    char acq_date[DATE_STRING_LEN + 1];    /* acquisition date */
    char prod_date[DATE_STRING_LEN + 1];   /* production date */
    char acq_time[TIME_STRING_LEN + 1];    /* acquisition time */
    char temp[STR_SIZE]; /* temporary string */
    int i;               /* looping variable */
    int refl_indx=0;     /* band index in XML file for the reflectance band */
    int th_indx=9;       /* band index in XML file for the thermal band */
    int pan_indx=7;      /* band index in XML file for the pan band */
    int qa_indx=11;      /* band index in XML file for the QA band */
    Espa_global_meta_t *gmeta = &metadata->global; /* pointer to global meta */

    /* Initialize the input fields */
    this->meta.sat = SAT_NULL;
    this->meta.inst = INST_NULL;
    this->meta.acq_date.fill = true;
    this->meta.time_fill = true;
    this->meta.prod_date.fill = true;
    this->meta.sun_zen = ANGLE_FILL;
    this->meta.sun_az = ANGLE_FILL;
    this->meta.wrs_sys = (Wrs_t)WRS_NULL;
    this->meta.ipath = -1;
    this->meta.irow = -1;
    this->meta.fill = INPUT_FILL;
    this->size.nsamps = this->size.nlines = -1;
    this->meta.gain_set = false;

    this->nband = 0;
    for (ib = 0; ib < NBAND_REFL_MAX; ib++)
    {
        this->meta.iband[ib] = -1;
        this->meta.gain[ib] = GAIN_BIAS_FILL;
        this->meta.bias[ib] = GAIN_BIAS_FILL;
        this->file_name[ib] = NULL;
        this->open[ib] = false;
        this->fp_bin[ib] = NULL;
    }

    this->nband_th = 0;
    for (ib = 0; ib < NBAND_THM_MAX; ib++)
    {
        this->meta.iband_th[ib] = -1;
        this->meta.gain_th[ib] = GAIN_BIAS_FILL;
        this->meta.bias_th[ib] = GAIN_BIAS_FILL;
        this->file_name_th[ib] = NULL;
        this->open_th[ib] = false;
        this->fp_bin_th[ib] = NULL;
    }

    this->nband_pan = 0;
    for (ib = 0; ib < NBAND_PAN_MAX; ib++)
    {
        this->meta.iband_pan[ib] = -1;
        this->meta.gain_pan[ib] = GAIN_BIAS_FILL;
        this->meta.bias_pan[ib] = GAIN_BIAS_FILL;
        this->file_name_pan[ib] = NULL;
        this->open_pan[ib] = false;
        this->fp_bin_pan[ib] = NULL;
    }

    this->nband_qa = 0;
    for (ib = 0; ib < NBAND_QA_MAX; ib++)
    {
        this->meta.iband_qa[ib] = -1;
        this->file_name_qa[ib] = NULL;
        this->open_qa[ib] = false;
        this->fp_bin_qa[ib] = NULL;
    }

    /* Pull the appropriate data from the XML file */
    acq_date[0] = acq_time[0] = '\0';
    prod_date[0] = '\0';
    if (!strcmp (gmeta->satellite, "LANDSAT_8"))
        this->meta.sat = SAT_LANDSAT_8;
    else
    {
        sprintf (errmsg, "Unsupported satellite: %s", gmeta->satellite);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (!strcmp (gmeta->instrument, "OLI_TIRS"))
        this->meta.inst = INST_OLI_TIRS;
    else if (!strcmp (gmeta->instrument, "OLI"))
        this->meta.inst = INST_OLI;
    else
    {
        sprintf (errmsg, "Unsupported instrument: %s", gmeta->instrument);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    strcpy (acq_date, gmeta->acquisition_date);
    strcpy (acq_time, gmeta->scene_center_time);
    this->meta.time_fill = false;

    /* Make sure the acquisition time is not too long (i.e. contains too
       many decimal points for the date/time routines).  The time should be
       hh:mm:ss.ssssssZ (see DATE_FORMAT_DATEA_TIME in date.h) which is 16
       characters long.  If the time is longer than that, just chop it off. */
    if (strlen (acq_time) > 16)
        sprintf (&acq_time[15], "Z");

    this->meta.sun_zen = gmeta->solar_zenith;
    if (this->meta.sun_zen < -90.0 || this->meta.sun_zen > 90.0)
    {
        sprintf (errmsg, "Solar zenith angle is out of range: %f",
            this->meta.sun_zen);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    this->meta.sun_az = gmeta->solar_azimuth;
    if (this->meta.sun_az < -360.0 || this->meta.sun_az > 360.0)
    {
        sprintf (errmsg, "Solar azimuth angle is out of range: %f",
            this->meta.sun_az);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    switch (gmeta->wrs_system)
    {
        case 1: this->meta.wrs_sys = WRS_1; break;
        case 2: this->meta.wrs_sys = WRS_2; break;
        default:
            sprintf (errmsg, "Invalid WRS system: %d", gmeta->wrs_system);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
    }
    this->meta.ipath = gmeta->wrs_path;
    this->meta.irow = gmeta->wrs_row;

    if (this->meta.inst == INST_OLI_TIRS)
    {
        this->nband = 8;        /* number of reflectance bands */
        for (ib = 0; ib < this->nband-1; ib++)
            this->meta.iband[ib] = ib+1;
        this->meta.iband[7] = 9;

        this->nband_th = 2;     /* number of thermal bands */
        this->meta.iband_th[0] = 10;
        this->meta.iband_th[1] = 11;

        this->nband_pan = 1;    /* number of pan bands */
        this->meta.iband_pan[0] = 8;

        this->nband_qa = 1;     /* number of QA bands */
        this->meta.iband_pan[0] = 12;
    }
    else if (this->meta.inst == INST_OLI)
    {
        this->nband = 8;        /* number of reflectance bands */
        for (ib = 0; ib < this->nband-1; ib++)
            this->meta.iband[ib] = ib+1;
        this->meta.iband[7] = 9;

        this->nband_th = 0;     /* number of thermal bands */

        this->nband_pan = 1;    /* number of pan bands */
        this->meta.iband_pan[0] = 8;

        this->nband_qa = 1;     /* number of QA bands */
        this->meta.iband_pan[0] = 12;
    }

    /* Find band 1, band 10, and band 8 in the input XML file to obtain
       band-related information for the reflectance, thermal, and pan bands */
    for (i = 0; i < metadata->nbands; i++)
    {
        if (!strcmp (metadata->band[i].name, "band1"))
        {
            /* this is the index we'll use for reflectance band info */
            refl_indx = i;

            /* get the band1 info */
            this->meta.gain[0] = metadata->band[i].toa_gain;
            this->meta.bias[0] = metadata->band[i].toa_bias;
            this->file_name[0] = strdup (metadata->band[i].file_name);

            /* get the production date but only the date portion
               (yyyy-mm-dd) */
            strncpy (prod_date, metadata->band[i].production_date, 10);
            prod_date[10] = '\0';
        }
        else if (!strcmp (metadata->band[i].name, "band2"))
        {
            /* get the band2 info */
            this->meta.gain[1] = metadata->band[i].toa_gain;
            this->meta.bias[1] = metadata->band[i].toa_bias;
            this->file_name[1] = strdup (metadata->band[i].file_name);
        }
        else if (!strcmp (metadata->band[i].name, "band3"))
        {
            /* get the band3 info */
            this->meta.gain[2] = metadata->band[i].toa_gain;
            this->meta.bias[2] = metadata->band[i].toa_bias;
            this->file_name[2] = strdup (metadata->band[i].file_name);
        }
        else if (!strcmp (metadata->band[i].name, "band4"))
        {
            /* get the band4 info */
            this->meta.gain[3] = metadata->band[i].toa_gain;
            this->meta.bias[3] = metadata->band[i].toa_bias;
            this->file_name[3] = strdup (metadata->band[i].file_name);
        }
        else if (!strcmp (metadata->band[i].name, "band5"))
        {
            /* get the band5 info */
            this->meta.gain[4] = metadata->band[i].toa_gain;
            this->meta.bias[4] = metadata->band[i].toa_bias;
            this->file_name[4] = strdup (metadata->band[i].file_name);
        }
        else if (!strcmp (metadata->band[i].name, "band6"))
        {
            /* get the band6 info */
            this->meta.gain[5] = metadata->band[i].toa_gain;
            this->meta.bias[5] = metadata->band[i].toa_bias;
            this->file_name[5] = strdup (metadata->band[i].file_name);
        }
        else if (!strcmp (metadata->band[i].name, "band7"))
        {
            /* get the band7 info */
            this->meta.gain[6] = metadata->band[i].toa_gain;
            this->meta.bias[6] = metadata->band[i].toa_bias;
            this->file_name[6] = strdup (metadata->band[i].file_name);
        }
        else if (!strcmp (metadata->band[i].name, "band9"))
        {
            /* get the band9 info */
            this->meta.gain[7] = metadata->band[i].toa_gain;
            this->meta.bias[7] = metadata->band[i].toa_bias;
            this->file_name[7] = strdup (metadata->band[i].file_name);
        }

        else if (!strcmp (metadata->band[i].name, "band8"))
        {
            /* this is the index we'll use for pan band info */
            pan_indx = i;

            /* get the band8 info */
            this->meta.gain_pan[0] = metadata->band[i].toa_gain;
            this->meta.bias_pan[0] = metadata->band[i].toa_bias;
            this->file_name_pan[0] = strdup (metadata->band[i].file_name);
        }

        /* NOTE: band10 and band11 won't exist in the input XML file */
        else if (!strcmp (metadata->band[i].name, "band10"))
        {
            /* this is the index we'll use for thermal band info */
            th_indx = i;

            /* get the band10 info */
            this->meta.gain_th[0] = metadata->band[i].toa_gain;
            this->meta.bias_th[0] = metadata->band[i].toa_bias;
            this->file_name_th[0] = strdup (metadata->band[i].file_name);
        }
        else if (!strcmp (metadata->band[i].name, "band11"))
        {
            /* get the band11 info */
            this->meta.gain_th[1] = metadata->band[i].toa_gain;
            this->meta.bias_th[1] = metadata->band[i].toa_bias;
            this->file_name_th[1] = strdup (metadata->band[i].file_name);
        }

        else if (!strcmp (metadata->band[i].name, "qa"))
        {
            /* this is the index we'll use for qa band info */
            qa_indx = i;

            /* get the QA band info */
            this->file_name_qa[0] = strdup (metadata->band[i].file_name);
        }
    }  /* for i */

    /* Get the size of the reflectance, thermal, pan, etc. bands by using
       the representative band in the XML file */
    this->size.nsamps = metadata->band[refl_indx].nsamps;
    this->size.nlines = metadata->band[refl_indx].nlines;
    this->size.pixsize[0] = metadata->band[refl_indx].pixel_size[0];
    this->size.pixsize[1] = metadata->band[refl_indx].pixel_size[1];
    this->scale_factor = metadata->band[refl_indx].scale_factor;

    if (this->meta.inst == INST_OLI_TIRS)
    {  /* skip for OLI */
        this->size_th.nsamps = metadata->band[th_indx].nsamps;
        this->size_th.nlines = metadata->band[th_indx].nlines;
        this->size_th.pixsize[0] = metadata->band[th_indx].pixel_size[0];
        this->size_th.pixsize[1] = metadata->band[th_indx].pixel_size[1];
        this->scale_factor_th = metadata->band[th_indx].scale_factor;
    }

    this->size_pan.nsamps = metadata->band[pan_indx].nsamps;
    this->size_pan.nlines = metadata->band[pan_indx].nlines;
    this->size_pan.pixsize[0] = metadata->band[pan_indx].pixel_size[0];
    this->size_pan.pixsize[1] = metadata->band[pan_indx].pixel_size[1];
    this->scale_factor_pan = metadata->band[pan_indx].scale_factor;

    this->size_qa.nsamps = metadata->band[qa_indx].nsamps;
    this->size_qa.nlines = metadata->band[qa_indx].nlines;
    this->size_qa.pixsize[0] = metadata->band[qa_indx].pixel_size[0];
    this->size_qa.pixsize[1] = metadata->band[qa_indx].pixel_size[1];

    /* Check WRS path/rows */
    if (this->meta.wrs_sys == WRS_1)
    {
        if (this->meta.ipath > 251)
        {
            sprintf (errmsg, "WRS path number out of range: %d",
                this->meta.ipath);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        else if (this->meta.irow > 248)
        {
            sprintf (errmsg, "WRS row number out of range: %d",
                this->meta.irow);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }
    else if (this->meta.wrs_sys == WRS_2)
    {
        if (this->meta.ipath > 233)
        {
            sprintf (errmsg, "WRS path number out of range: %d",
                this->meta.ipath);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        else if (this->meta.irow > 248)
        {
            sprintf (errmsg, "WRS row number out of range: %d",
                this->meta.irow);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Check satellite/instrument combination */
    if (this->meta.inst == INST_OLI_TIRS ||
        this->meta.inst == INST_OLI)
    {
        if (this->meta.sat != SAT_LANDSAT_8)
        {
            sprintf (errmsg, "Invalid instrument/satellite combination");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Convert the acquisition date/time values from string to date struct */
    sprintf (temp, "%sT%s", acq_date, acq_time);
    if (!date_init (&this->meta.acq_date, temp, DATE_FORMAT_DATEA_TIME))
    {
        sprintf (errmsg, "Converting the acquisition date and time");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Convert the production date value from string to date struct */
    if (!date_init (&this->meta.prod_date, prod_date, DATE_FORMAT_DATEA))
    {
        sprintf (errmsg, "Converting the production date and time");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    return (SUCCESS);
}

