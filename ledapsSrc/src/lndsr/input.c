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

 Modified on 8/17/2012 by Gail Schmidt, USGS EROS
  When freeing the thermal band in FreeInput, the QA SDS is not valid
  and should not be attempted to be freed.  OpenInput does not read the
  QA for the thermal band.
  
 Revision 2.0 2014/01/24
 Gail Schmidt, USGS EROS
 Modified to use ESPA internal raw binary format

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
	FreeOutput - Free the 'input' data structure memory.

   2. 'OpenInput' must be called before any of the other routines.  
   3. 'FreeInput' should be used to free the 'input' data structure.

!END****************************************************************************
*/

#include <stdlib.h>
#include <math.h>
#include "input.h"
#include "error.h"
#include "mystring.h"
#include "myhdf.h"
#include "const.h"

#define INPUT_FILL (-9999)

/* Functions */
Input_t *OpenInput(Espa_internal_meta_t *metadata, bool thermal)
/* 
!C******************************************************************************

!Description: 'OpenInput' sets up the 'input' data structure, opens the
 input raw binary files for read access.
 
!Input Parameters:
 metadata     'Espa_internal_meta_t' data structure with XML info
 thermal      boolean to indicate if thermal data is being processed

!Output Parameters:
 (returns)      'input' data structure or NULL when an error occurs

!Team Unique Header:

!END****************************************************************************
*/
{
  Input_t *this = NULL;
  char *error_string = (char *)NULL;
  int ib;

  /* Create the Input data structure */
  this = (Input_t *)malloc(sizeof(Input_t));
  if (this == (Input_t *)NULL) 
    RETURN_ERROR("allocating Input data structure", "OpenInput", NULL);

  /* Initialize and get input from header file */
  if (!GetXMLInput (this, metadata, thermal)) {
    free(this);
    this = NULL;
    RETURN_ERROR("getting input from header file", "OpenInput", NULL);
  }

  /* Open TOA reflectance files for access */
  for (ib = 0; ib < this->nband; ib++) {
    this->fp_bin[ib] = fopen(this->file_name[ib], "r");
    if (this->fp_bin[ib] == NULL) {
      error_string = "opening input TOA binary file";
      break;
    }
    this->open[ib] = true;
  }

  /* Open QA file for access */
  this->fp_bin_qa = fopen(this->file_name_qa, "r");
  if (this->fp_bin_qa == NULL) 
    error_string = "opening QA binary file";
  else
    this->open_qa = true;

  if (error_string != NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      free(this->file_name[ib]);
      this->file_name[ib] = NULL;

      if (this->open[ib]) {
        fclose(this->fp_bin[ib]);
        this->open[ib] = false;
      }
    }
    free(this->file_name_qa);
    this->file_name_qa = NULL;
    fclose(this->fp_bin_qa);  
    this->open_qa = false;
    free(this);
    this = NULL;
    RETURN_ERROR(error_string, "OpenInput", NULL);
  }

  return this;
}


bool CloseInput(Input_t *this)
/* 
!C******************************************************************************

!Description: 'CloseInput' closes the input file.
 
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
      fclose(this->fp_bin[ib]);
      this->open[ib] = false;
    }
  }

  /*** now close the QA file ***/
  if (this->open_qa) 
  {
    fclose(this->fp_bin_qa);
    this->open_qa = false;
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

  if (this != (Input_t *)NULL) {
    for (ib = 0; ib < this->nband; ib++) {
      free(this->file_name[ib]);
      this->file_name[ib] = NULL;
    }
    free(this->file_name_qa);
    this->file_name_qa = NULL;

    free(this);
    this = NULL;
  }

  return true;
}


bool GetInputLine(Input_t *this, int iband, int iline, int16 *line)
{
  long loc;
  void *buf_void = NULL;

  /* Check the parameters */
  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "GetIntputLine", false);
  if (iband < 0 || iband >= this->nband)
    RETURN_ERROR("band number out of range", "GetInputLine", false);
  if (iline < 0 || iline >= this->size.l)
    RETURN_ERROR("line number out of range", "GetInputLine", false);
  if (!this->open[iband])
    RETURN_ERROR("band not open", "GetInputLine", false);

  /* Read the data */
  buf_void = (void *)line;
  loc = (long) (iline * this->size.s * sizeof(int16));
  if (fseek(this->fp_bin[iband], loc, SEEK_SET))
    RETURN_ERROR("error seeking line (binary)", "GetInputLine", false);
  if (fread(buf_void, sizeof(int16), (size_t)this->size.s, 
            this->fp_bin[iband]) != (size_t)this->size.s)
    RETURN_ERROR("error reading line (binary)", "GetInputLine", false);

  return true;
}


bool GetInputQALine(Input_t *this, int iline, uint8 *line) 
{
  long loc;
  void *buf_void = NULL;

  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "GetInputQALine", false);
  if (iline < 0  ||  iline >= this->size.l) 
    RETURN_ERROR("line index out of range", "GetInputQALine", false);
  if (!this->open_qa)
    RETURN_ERROR("QA band not open", "GetInputQALine", false);

  buf_void = (void *)line;
  loc = (long) (iline * this->size.s * sizeof(uint8));
  if (fseek(this->fp_bin_qa, loc, SEEK_SET))
    RETURN_ERROR("error seeking line (binary)", "GetInputQALine", false);
  if (fread(buf_void, sizeof(uint8), (size_t)this->size.s, 
            this->fp_bin_qa) != (size_t)this->size.s)
    RETURN_ERROR("error reading line (binary)", "GetInputQALine", false);

  return true;
}


bool InputMetaCopy(Input_meta_t *this, int nband, Input_meta_t *copy) 
{
  int ib;

  if (this == (Input_meta_t *)NULL) 
    RETURN_ERROR("invalid input structure (original)", "InputMetaCopy", false);
  if (copy == (Input_meta_t *)NULL) 
    RETURN_ERROR("invalid input structure (copy)", "InputMetaCopy", false);
  if (nband < 1  ||  nband > NBAND_REFL_MAX)
    RETURN_ERROR("invalid number of bands", "InputMetaCopy", false);

  copy->sat = this->sat;
  copy->inst = this->inst;
  if (!DateCopy(&this->acq_date, &copy->acq_date)) 
    RETURN_ERROR("copying acquisition date/time", "InputMetaCopy", false);
  copy->sun_zen = this->sun_zen;
  copy->sun_az = this->sun_az;
  copy->wrs_sys = this->wrs_sys;
  copy->ipath = this->ipath;
  copy->irow = this->irow;
  copy->fill = this->fill;

  for (ib = 0; ib < nband; ib++) {
    copy->iband[ib] = this->iband[ib];
  }
  copy->iband_qa = this->iband_qa;

  return true;
}


#define DATE_STRING_LEN (50)
#define TIME_STRING_LEN (50)

bool GetXMLInput(Input_t *this, Espa_internal_meta_t *metadata, bool thermal)
/* 
!C******************************************************************************

!Description: 'GetXMLInput' pulls input values from the XML structure.
 
!Input Parameters:
 this         'Input_t' data structure to be populated
 metadata     'Espa_internal_meta_t' data structure with XML info
 thermal      boolean to indicate if thermal data is being processed

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)
                  'false' = error getting metadata from the XML file

!Team Unique Header:

! Design Notes:
  1. This replaces the previous GetInputMeta so the input values are pulled
     from the XML file instead of the HDF and MTL files.

!END****************************************************************************
*/
{
    char *error_string = NULL;
    int ib;
    char acq_date[DATE_STRING_LEN + 1];
    char acq_time[TIME_STRING_LEN + 1];
    char temp[MAX_STR_LEN + 1];
    int i;               /* looping variable */
    int indx=-1;         /* band index in XML file for band1 or band6 */
    Espa_global_meta_t *gmeta = &metadata->global; /* pointer to global meta */

    /* Initialize the input fields.  Set file type to binary, since that is
       the ESPA internal format for the input TOA products. */
    acq_date[0] = acq_time[0] = '\0';
    this->meta.sat = SAT_NULL;
    this->meta.inst = INST_NULL;
    this->meta.acq_date.fill = true;
    this->meta.sun_zen = ANGLE_FILL;
    this->meta.sun_az = ANGLE_FILL;
    this->meta.wrs_sys = (Wrs_t)WRS_FILL;
    this->meta.ipath = -1;
    this->meta.irow = -1;
    this->meta.fill = INPUT_FILL;
    this->nband = 0;
    this->size.s = this->size.l = -1;
    for (ib = 0; ib < NBAND_REFL_MAX; ib++)
    {
        this->meta.iband[ib] = -1;
        this->file_name[ib] = NULL;
        this->open[ib] = false;
        this->fp_bin[ib] = NULL;
    }
    this->open_qa = false;
    this->file_name_qa = NULL;
    this->fp_bin_qa = NULL;

    /* Pull the appropriate data from the XML file */
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
        if (!thermal)
        {  /* reflectance bands */
            this->nband = 6;     /* number of reflectance bands */
            this->meta.iband[0] = 1;
            this->meta.iband[1] = 2;
            this->meta.iband[2] = 3;
            this->meta.iband[3] = 4;
            this->meta.iband[4] = 5;
            this->meta.iband[5] = 7;
        }
        else
        {  /* thermal band */
            this->nband = 1;     /* number of thermal bands */
            this->meta.iband[0] = 6;
        }
    }

    /* Find TOA band 1 in the input XML file to obtain band-related
       information */
    if (!thermal)
    {  /* reflectance bands */
        for (i = 0; i < metadata->nbands; i++)
        {
            if (!strcmp (metadata->band[i].name, "toa_band1") &&
                !strcmp (metadata->band[i].product, "toa_refl"))
            {
                /* this is the index we'll use for reflectance band info */
                indx = i;

                /* get the band1 info */
                this->file_name[0] = strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "toa_band2") &&
                !strcmp (metadata->band[i].product, "toa_refl"))
            {
                this->file_name[1] = strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "toa_band3") &&
                !strcmp (metadata->band[i].product, "toa_refl"))
            {
                this->file_name[2] = strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "toa_band4") &&
                !strcmp (metadata->band[i].product, "toa_refl"))
            {
                this->file_name[3] = strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "toa_band5") &&
                !strcmp (metadata->band[i].product, "toa_refl"))
            {
                this->file_name[4] = strdup (metadata->band[i].file_name);
            }
            else if (!strcmp (metadata->band[i].name, "toa_band7") &&
                !strcmp (metadata->band[i].product, "toa_refl"))
            {
                this->file_name[5] = strdup (metadata->band[i].file_name);
            }

            if (!strcmp (metadata->band[i].name, "toa_qa") &&
                !strcmp (metadata->band[i].product, "toa_refl"))
            {
                this->file_name_qa = strdup (metadata->band[i].file_name);
            }
        }  /* for i */
    }  /* for i */
    else
    {  /* thermal band */
        for (i = 0; i < metadata->nbands; i++)
        {
            if (!strcmp (metadata->band[i].name, "toa_band6") &&
                !strcmp (metadata->band[i].product, "toa_bt"))
            {
                /* this is the index we'll use for reflectance band info */
                indx = i;

                /* get the band1 info */
                this->file_name[0] = strdup (metadata->band[i].file_name);
            }

            if (!strcmp (metadata->band[i].name, "toa_bt_qa") &&
                !strcmp (metadata->band[i].product, "toa_bt"))
            {
                this->file_name_qa = strdup (metadata->band[i].file_name);
            }
        }  /* for i */
    }

    if (indx == -1)
    {
        error_string = "not able to find the reflectance/thermal index band";
        RETURN_ERROR (error_string, "GetXMLInput", true);
    }

    /* Pull the reflectance info from band1 in the XML file */
    this->size.s = metadata->band[indx].nsamps;
    this->size.l = metadata->band[indx].nlines;

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

    return true;
}

