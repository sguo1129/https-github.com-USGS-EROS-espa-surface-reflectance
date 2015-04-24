/*
!C****************************************************************************

!File: input.c
  
!Description: Functions reading data from the input data file.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

 Revision 1.1 2002/05/02
 Robert Wolfe
 Added handling for SDS's with ranks greater than 2.

 Revision 1.2 2013/03/12
 Gail Schmidt, USGS EROS
 Modified the read header to process the acquisition date, even if it's too
 long by parsing it down to the appropriate number of allowed characters.

 Revision 2.0 2014/01/24
 Gail Schmidt, USGS EROS
 Modified to use ESPA internal raw binary format

 Revision 2015/03/31
 Gail Schmidt, USGS EROS
 Added support for pulling the TOA reflectance parameters, K1/K2 consts, and
    earth-sun distance from the XML, if they exist.

!Team Unique Header:
  This software was developed by the MODIS Land Science Team Support 
  Group for the Laboratory for Terrestrial Physics (Code 922) at the 
  National Aeronautics and Space Administration, Goddard Space Flight 
  Center, under NASA Task 92-012-00.

 ! References and Credits:
  ! MODIS Science Team Member:
      Christopher O. Justice
      MODIS Land Science Team           University of Maryland
      justice@hermes.geog.umd.edu       Dept. of Geography
      phone: 301-405-1600               1113 LeFrak Hall
                                        College Park, MD, 20742

  ! Developers:
      Robert E. Wolfe (Code 922)
      MODIS Land Team Support Group     Raytheon ITSS
      robert.e.wolfe.1@gsfc.nasa.gov    4400 Forbes Blvd.
      phone: 301-614-5508               Lanham, MD 20770  
  
 ! Design Notes:
   1. The following public functions handle the input data:

	OpenInput - Setup 'input' data structure and open file for access.
	CloseInput - Close the input file.
	FreeInput - Free the 'input' data structure memory.

   2. 'OpenInput' must be called before any of the other routines.  
   3. 'FreeInput' should be used to free the 'input' data structure.

!END****************************************************************************
*/

#include <stdlib.h>
#include "input.h"
#include "util.h"
#include "error.h"
#include "mystring.h"
#include "const.h"
#include "date.h"
#define INPUT_FILL (0)

/* Functions */
Input_t *OpenInput(Espa_internal_meta_t *metadata)
/* 
!C******************************************************************************

!Description: 'OpenInput' sets up the 'input' data structure, opens the
 input raw binary files for read access.
 
!Input Parameters:
 metadata     'Espa_internal_meta_t' data structure with XML info

!Output Parameters:
 (returns)      'input' data structure or NULL when an error occurs

!Team Unique Header:

!END****************************************************************************
*/
{
  Input_t *this = NULL;
  char *error_string = NULL;
  int ib;

  /* Create the Input data structure */
  this = (Input_t *)malloc(sizeof(Input_t));
  if (this == NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", NULL);

  /* Initialize and get input from header file */
  if (!GetXMLInput (this, metadata)) {
    free(this);
    this = NULL;
    RETURN_ERROR("getting input from header file", "OpenInput", NULL);
  }

  /* Open files for access */
  if (this->file_type == INPUT_TYPE_BINARY) {
    for (ib = 0; ib < this->nband; ib++) {
      this->fp_bin[ib] = fopen(this->file_name[ib], "r");
      if (this->fp_bin[ib] == NULL) {
        error_string = "opening binary file";
        break;
      }
      this->open[ib] = true;
    }
    if ( this->nband_th == 1 ) {
      this->fp_bin_th = fopen(this->file_name_th, "r");
      if (this->fp_bin_th == NULL) 
        error_string = "opening thermal binary file";
      else
        this->open_th = true;
    }
  } else 
    error_string = "invalid file type";

  if (error_string != NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      free(this->file_name[ib]);
      this->file_name[ib] = NULL;

      if (this->open[ib]) {
        if ( this->file_type == INPUT_TYPE_BINARY )
          fclose(this->fp_bin[ib]);
        this->open[ib] = false;
      }
    }
    free(this->file_name_th);
    this->file_name_th = NULL;
    if ( this->file_type == INPUT_TYPE_BINARY )
      fclose(this->fp_bin_th);  
    this->open_th = false;
    free(this);
    this = NULL;
    RETURN_ERROR(error_string, "OpenInput", NULL);
  }

  return this;
}

bool GetInputLine(Input_t *this, int iband, int iline, unsigned char *line) 
{
  long loc;
  void *buf_void = NULL;

  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "GetInputLine", false);
  if (iband < 0  ||  iband >= this->nband) 
    RETURN_ERROR("band index out of range", "GetInputLine", false);
  if (iline < 0  ||  iline >= this->size.l) 
    RETURN_ERROR("line index out of range", "GetInputLine", false);
  if (!this->open[iband])
    RETURN_ERROR("band not open", "GetInputLine", false);

  buf_void = (void *)line;
  if (this->file_type == INPUT_TYPE_BINARY) {
    loc = (long) (iline * this->size.s * sizeof(uint8));
    if (fseek(this->fp_bin[iband], loc, SEEK_SET))
      RETURN_ERROR("error seeking line (binary)", "GetInputLine", false);
    if (fread(buf_void, sizeof(uint8), (size_t)this->size.s, 
              this->fp_bin[iband]) != (size_t)this->size.s)
      RETURN_ERROR("error reading line (binary)", "GetInputLine", false);
  }

  return true;
}

bool GetInputLineTh(Input_t *this, int iline, unsigned char *line) 
{
  long loc;
  void *buf_void = NULL;

  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "GetInputLine", false);
  if ( this->nband_th < 1 ) 
    RETURN_ERROR("no thermal input band", "GetInputLine", false);
  if (iline < 0  ||  iline >= this->size_th.l) 
    RETURN_ERROR("line index out of range", "GetInputLine", false);
  if (!this->open_th)
    RETURN_ERROR("band not open", "GetInputLine", false);

  buf_void = (void *)line;
  if (this->file_type == INPUT_TYPE_BINARY) {
    loc = (long) (iline * this->size_th.s * sizeof(uint8));
    if (fseek(this->fp_bin_th, loc, SEEK_SET))
      RETURN_ERROR("error seeking line (binary)", "GetInputLine", false);
    if (fread(buf_void, sizeof(uint8), (size_t)this->size_th.s, 
              this->fp_bin_th) != (size_t)this->size_th.s)
      RETURN_ERROR("error reading line (binary)", "GetInputLine", false);
  }

  return true;
}


bool CloseInput(Input_t *this)
/* 
!C******************************************************************************

!Description: 'CloseInput' ends SDS access and closes the input file.
 
!Input Parameters:
 this           'input' data structure

!Output Parameters:
 this           'input' data structure; the following fields are modified:
                   open
 (returns)      status:
                  'true' = okay
                  'false' = error return

!Team Unique Header:

!END****************************************************************************
*/
{
  int ib;
  bool none_open;

  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "CloseInput", false);

  none_open = true;
  for (ib = 0; ib < this->nband; ib++) {
    if (this->open[ib]) {
      none_open = false;
      if (this->file_type == INPUT_TYPE_BINARY)
        fclose(this->fp_bin[ib]);
      this->open[ib] = false;
    }
  }

  /*** now close the thermal file ***/
  if (this->open_th) 
  {
    if (this->file_type == INPUT_TYPE_BINARY)
      fclose(this->fp_bin_th);
    this->open_th = false;
  }

  if (none_open)
    RETURN_ERROR("no files open", "CloseInput", false);

  return true;
}


bool FreeInput(Input_t *this)
/* 
!C******************************************************************************

!Description: 'FreeInput' frees the 'input' data structure memory.
 
!Input Parameters:
 this           'input' data structure

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)

!Team Unique Header:

!END****************************************************************************
*/
{
  int ib;

  if (this != NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      free(this->file_name[ib]);
      this->file_name[ib] = NULL;
    }
    free(this->file_name_th);
    this->file_name_th = NULL;

    free(this);
    this = NULL;
  }

  return true;
}

bool InputMetaCopy(Input_meta_t *this, int nband, Input_meta_t *copy) 
{
  int ib;

  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "InputMetaCopy", false);

  copy->sat = this->sat;
  copy->inst = this->inst;
  if (!DateCopy(&this->acq_date, &copy->acq_date)) 
    RETURN_ERROR("copying acquisition date/time", "InputMetaCopy", false);
  if (!DateCopy(&this->prod_date, &copy->prod_date)) 
    RETURN_ERROR("copying production date/time", "InputMetaCopy", false);
  copy->time_fill = this->time_fill;
  copy->sun_zen = this->sun_zen;
  copy->sun_az = this->sun_az;
  copy->earth_sun_dist = this->earth_sun_dist;
  copy->wrs_sys = this->wrs_sys;
  copy->ipath = this->ipath;
  copy->irow = this->irow;
  copy->fill = this->fill;

  copy->iband_th = this->iband_th;
  for (ib = 0; ib < nband; ib++) {
    copy->iband[ib] = this->iband[ib];
    copy->rad_gain[ib] = this->rad_gain[ib];
    copy->rad_bias[ib] = this->rad_bias[ib];
    copy->refl_gain[ib] = this->refl_gain[ib];
    copy->refl_bias[ib] = this->refl_bias[ib];
  }
  copy->rad_gain_th = this->rad_gain_th;
  copy->rad_bias_th = this->rad_bias_th;
  copy->k1_const = this->k1_const;
  copy->k2_const = this->k2_const;
  copy->use_toa_refl_consts = this->use_toa_refl_consts;
  return true;
}


#define DATE_STRING_LEN (50)
#define TIME_STRING_LEN (50)

bool GetXMLInput(Input_t *this, Espa_internal_meta_t *metadata)
/* 
!C******************************************************************************

!Description: 'GetXMLInput' pulls input values from the XML structure.
 
!Input Parameters:
 this         'Input_t' data structure to be populated
 metadata     'Espa_internal_meta_t' data structure with XML info

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)
                  'false' = error determining if the gains/biases were provided

!Team Unique Header:

! Design Notes:
  1. This replaces the previous GetHeaderInput so the input values are pulled
     from the XML file instead of the header file (*.metadata.txt).
  2. Given that LPGS writes the gain values, the gain settings (HIGH, LOW) are
     no longer needed.

!END****************************************************************************
*/
{
    char *error_string = NULL;
    int ib;
    char acq_date[DATE_STRING_LEN + 1];
    char prod_date[DATE_STRING_LEN + 1];
    char acq_time[TIME_STRING_LEN + 1];
    char temp[MAX_STR_LEN + 1];
    int i;               /* looping variable */
    int refl_indx=0;     /* band index in XML file for the reflectance band */
    int th_indx=5;       /* band index in XML file for the thermal band */
    Espa_global_meta_t *gmeta = &metadata->global; /* pointer to global meta */

    /* Initialize the input fields.  Set file type to binary, since that is
       the ESPA internal format for the input L1G/T products. */
    this->file_type = INPUT_TYPE_BINARY;
    this->meta.sat = SAT_NULL;
    this->meta.inst = INST_NULL;
    this->meta.acq_date.fill = true;
    this->meta.time_fill = true;
    this->meta.prod_date.fill = true;
    this->meta.sun_zen = ANGLE_FILL;
    this->meta.sun_az = ANGLE_FILL;
    this->meta.wrs_sys = (Wrs_t)WRS_FILL;
    this->meta.ipath = -1;
    this->meta.irow = -1;
    this->meta.fill = INPUT_FILL;
    this->nband = 0;
    this->nband_th = 0;
    this->size.s = this->size.l = -1;
    for (ib = 0; ib < NBAND_REFL_MAX; ib++)
    {
        this->meta.iband[ib] = -1;
        this->meta.rad_gain[ib] = GAIN_BIAS_FILL;
        this->meta.rad_bias[ib] = GAIN_BIAS_FILL;
        this->meta.refl_gain[ib] = GAIN_BIAS_FILL;
        this->meta.refl_bias[ib] = GAIN_BIAS_FILL;
        this->file_name[ib] = NULL;
        this->open[ib] = false;
        this->fp_bin[ib] = NULL;
    }
    this->nband_th = 0;
    this->open_th = false;
    this->meta.rad_gain_th = GAIN_BIAS_FILL;
    this->meta.rad_bias_th = GAIN_BIAS_FILL;
    this->meta.k1_const = GAIN_BIAS_FILL;
    this->meta.k2_const = GAIN_BIAS_FILL;
    this->file_name_th = NULL;
    this->fp_bin_th = NULL;

    /* Pull the appropriate data from the XML file */
    acq_date[0] = acq_time[0] = '\0';
    prod_date[0] = '\0';
    if (!strcmp (gmeta->satellite, "LANDSAT_1"))
        this->meta.sat = SAT_LANDSAT_1;
    else if (!strcmp (gmeta->satellite, "LANDSAT_2"))
        this->meta.sat = SAT_LANDSAT_2;
    else if (!strcmp (gmeta->satellite, "LANDSAT_3"))
        this->meta.sat = SAT_LANDSAT_3;
    else if (!strcmp (gmeta->satellite, "LANDSAT_4"))
        this->meta.sat = SAT_LANDSAT_4;
    else if (!strcmp (gmeta->satellite, "LANDSAT_5"))
        this->meta.sat = SAT_LANDSAT_5;
    else if (!strcmp (gmeta->satellite, "LANDSAT_7"))
        this->meta.sat = SAT_LANDSAT_7;
    else
    {
        sprintf (temp, "invalid satellite; value = %s", gmeta->satellite);
        RETURN_ERROR (temp, "GetXMLInput", true);
    }

    if (!strcmp (gmeta->instrument, "TM"))
        this->meta.inst = INST_TM;
    else if (!strncmp (gmeta->instrument, "ETM", 3))
        this->meta.inst = INST_ETM;
    else
    {
        sprintf (temp, "invalid instrument; value = %s", gmeta->instrument);
        RETURN_ERROR (temp, "GetXMLInput", true);
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
        error_string = "solar zenith angle out of range";
        RETURN_ERROR (error_string, "GetXMLInput", true);
    }
    this->meta.sun_zen *= RAD;   /* convert to radians */

    this->meta.sun_az = gmeta->solar_azimuth;
    if (this->meta.sun_az < -360.0 || this->meta.sun_az > 360.0)
    {
        error_string = "solar azimuth angle out of range";
        RETURN_ERROR (error_string, "GetXMLInput", true);
    }
    this->meta.sun_az *= RAD;    /* convert to radians */

    this->meta.earth_sun_dist = gmeta->earth_sun_dist;

    switch (gmeta->wrs_system)
    {
        case 1: this->meta.wrs_sys = WRS_1; break;
        case 2: this->meta.wrs_sys = WRS_2; break;
        default:
            sprintf (temp, "invalid WRS system; value = %d",
                gmeta->wrs_system);
            RETURN_ERROR (temp, "GetXMLInput", true);
    }
    this->meta.ipath = gmeta->wrs_path;
    this->meta.irow = gmeta->wrs_row;

    if (this->meta.inst == INST_TM || this->meta.inst == INST_ETM)
    {
        this->nband = 6;     /* number of reflectance bands */
        this->meta.iband[0] = 1;
        this->meta.iband[1] = 2;
        this->meta.iband[2] = 3;
        this->meta.iband[3] = 4;
        this->meta.iband[4] = 5;
        this->meta.iband[5] = 7;

        this->nband_th = 1;  /* number of thermal bands; only use 6L for ETM */
        this->meta.iband_th = 6;
    }

    /* Find band 1 and band 6/61 in the input XML file to obtain band-related
       information */
    this->meta.use_toa_refl_consts = false;
    for (i = 0; i < metadata->nbands; i++)
    {
        if (!strcmp (metadata->band[i].name, "band1") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* L1G or L1T */
        {
            /* this is the index we'll use for reflectance band info */
            refl_indx = i;

            /* get the band1 info */
            this->meta.rad_gain[0] = metadata->band[i].rad_gain;
            this->meta.rad_bias[0] = metadata->band[i].rad_bias;
            this->file_name[0] = strdup (metadata->band[i].file_name);

            /* get the production date but only the date portion (yyyy-mm-dd) */
            strncpy (prod_date, metadata->band[i].production_date, 10);
            prod_date[10] = '\0';

            /* are the TOA reflectance coefficients available in XML file? */
            if (existReflGB (metadata))
            {
                this->meta.use_toa_refl_consts = true;
                this->meta.refl_gain[0] = metadata->band[i].refl_gain;
                this->meta.refl_bias[0] = metadata->band[i].refl_bias;
            }
        }
        else if (!strcmp (metadata->band[i].name, "band2") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* L1G or L1T */
        {
            /* get the band2 info */
            this->meta.rad_gain[1] = metadata->band[i].rad_gain;
            this->meta.rad_bias[1] = metadata->band[i].rad_bias;
            this->file_name[1] = strdup (metadata->band[i].file_name);

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.refl_gain[1] = metadata->band[i].refl_gain;
                this->meta.refl_bias[1] = metadata->band[i].refl_bias;
            }
        }
        else if (!strcmp (metadata->band[i].name, "band3") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* L1G or L1T */
        {
            /* get the band3 info */
            this->meta.rad_gain[2] = metadata->band[i].rad_gain;
            this->meta.rad_bias[2] = metadata->band[i].rad_bias;
            this->file_name[2] = strdup (metadata->band[i].file_name);

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.refl_gain[2] = metadata->band[i].refl_gain;
                this->meta.refl_bias[2] = metadata->band[i].refl_bias;
            }
        }
        else if (!strcmp (metadata->band[i].name, "band4") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* L1G or L1T */
        {
            /* get the band4 info */
            this->meta.rad_gain[3] = metadata->band[i].rad_gain;
            this->meta.rad_bias[3] = metadata->band[i].rad_bias;
            this->file_name[3] = strdup (metadata->band[i].file_name);

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.refl_gain[3] = metadata->band[i].refl_gain;
                this->meta.refl_bias[3] = metadata->band[i].refl_bias;
            }
        }
        else if (!strcmp (metadata->band[i].name, "band5") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* L1G or L1T */
        {
            /* get the band5 info */
            this->meta.rad_gain[4] = metadata->band[i].rad_gain;
            this->meta.rad_bias[4] = metadata->band[i].rad_bias;
            this->file_name[4] = strdup (metadata->band[i].file_name);

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.refl_gain[4] = metadata->band[i].refl_gain;
                this->meta.refl_bias[4] = metadata->band[i].refl_bias;
            }
        }
        else if (!strcmp (metadata->band[i].name, "band7") &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* L1G or L1T */
        {
            /* get the band7 info */
            this->meta.rad_gain[5] = metadata->band[i].rad_gain;
            this->meta.rad_bias[5] = metadata->band[i].rad_bias;
            this->file_name[5] = strdup (metadata->band[i].file_name);

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.refl_gain[5] = metadata->band[i].refl_gain;
                this->meta.refl_bias[5] = metadata->band[i].refl_bias;
            }
        }

        if (!strcmp (metadata->band[i].name, "band6") &&
            this->meta.inst == INST_TM &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* L1G or L1T */
        {
            /* this is the index we'll use for thermal band info */
            th_indx = i;

            /* get the band6 info */
            this->meta.rad_gain_th = metadata->band[i].rad_gain;
            this->meta.rad_bias_th = metadata->band[i].rad_bias;
            this->file_name_th = strdup (metadata->band[i].file_name);

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.k1_const = metadata->band[i].k1_const;
                this->meta.k2_const = metadata->band[i].k2_const;
            }
        }
        else if (!strcmp (metadata->band[i].name, "band61") &&
            this->meta.inst == INST_ETM &&
            !strncmp (metadata->band[i].product, "L1", 2))  /* L1G or L1T */
        {
            /* this is the index we'll use for thermal band info */
            th_indx = i;

            /* get the band6 info */
            this->meta.rad_gain_th = metadata->band[i].rad_gain;
            this->meta.rad_bias_th = metadata->band[i].rad_bias;
            this->file_name_th = strdup (metadata->band[i].file_name);

            if (this->meta.use_toa_refl_consts)
            {
                this->meta.k1_const = metadata->band[i].k1_const;
                this->meta.k2_const = metadata->band[i].k2_const;
            }
        }
    }  /* for i */

    /* Let the user know if the XML params are being used or the LUT params */
    if (this->meta.use_toa_refl_consts)
    {
        printf ("Using the TOA reflectance coefficients, K1/K2 thermal "
            "constants, and earth-sun distance from the XML file.\n");
    }
    else
    {
        printf ("Using the hard-coded TOA reflectance coefficients, K1/K2 "
            "thermal constants, and earth-sun distance table.\n");
    }

    /* Pull the reflectance info from band1 in the XML file */
    this->size.s = metadata->band[refl_indx].nsamps;
    this->size.l = metadata->band[refl_indx].nlines;
    this->size_th.s = metadata->band[th_indx].nsamps;
    this->size_th.l = metadata->band[th_indx].nlines;

    /* Check WRS path/rows */
    if (this->meta.wrs_sys == WRS_1)
    {
        if (this->meta.ipath > 251)
            error_string = "WRS path number out of range";
        else if (this->meta.irow > 248)
            error_string = "WRS row number out of range";
    }
    else if (this->meta.wrs_sys == WRS_2)
    {
        if (this->meta.ipath > 233)
            error_string = "WRS path number out of range";
        else if (this->meta.irow > 248)
            error_string = "WRS row number out of range";
    }
    else
        error_string = "invalid WRS system";

    if (error_string != NULL)
    {
        RETURN_ERROR (error_string, "GetHeaderInput", true);
    }

    /* Check satellite/instrument combination */
    if (this->meta.inst == INST_MSS)
    {
        if (this->meta.sat != SAT_LANDSAT_1 &&
            this->meta.sat != SAT_LANDSAT_2 &&
            this->meta.sat != SAT_LANDSAT_3 &&
            this->meta.sat != SAT_LANDSAT_4 &&
            this->meta.sat != SAT_LANDSAT_5)
            error_string = "invalid insturment/satellite combination";
    }
    else if (this->meta.inst == INST_TM)
    {
        if (this->meta.sat != SAT_LANDSAT_4 &&
            this->meta.sat != SAT_LANDSAT_5)
            error_string = "invalid insturment/satellite combination";
    }
    else if (this->meta.inst == INST_ETM)
    {
        if (this->meta.sat != SAT_LANDSAT_7)
            error_string = "invalid insturment/satellite combination";
    }
    else
        error_string = "invalid instrument type";

    if (error_string != NULL)
    {
        RETURN_ERROR (error_string, "GetHeaderInput", true);
    }

    /* Convert the acquisition date/time values */
    sprintf (temp, "%sT%s", acq_date, acq_time);
    if (!DateInit (&this->meta.acq_date, temp, DATE_FORMAT_DATEA_TIME))
    {
        error_string = "converting acquisition date/time";
        RETURN_ERROR (error_string, "GetHeaderInput", false);
    }

    /* Convert the production date value */
    if (!DateInit (&this->meta.prod_date, prod_date, DATE_FORMAT_DATEA))
    {
        error_string = "converting production date";
        RETURN_ERROR (error_string, "GetHeaderInput", false);
    }

    return true;
}

