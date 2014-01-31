/*
!C****************************************************************************

!File: output.c
  
!Description: Functions creating and writing data to the product output file.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

 Revision 1.1 2012/10/18
 Gail Schmidt, USGS EROS
 Modified PutMetadata* to mark the bits that are actually set and not set
   (i.e. for lndcal, bit 6 actually isn't set and for lndth, bit 6 and fill
    are the only pixels that are set).

 Revision 1.1 2012/08/01
 Gail Schmidt, USGS EROS
 Modified OpenOutput to make sure the Output_t->buf array doesn't go past
   NBAND_REFL_MAX.  Also, initialized Output_t->qabuf at the same time.

 Revision 1.2 2013/01/22
 Gail Schmidt, USGS EROS
 Modified applications to use only one version and that is the
   LEDAPSVersion tag which will get updated with each release of LEDAPS

 Revision 1.3 2013/02/20
 Gail Schmidt, USGS EROS
 Modified the put metadata routine to write out the bounding coordinates.

 Revision 1.4 2013/03/22
 Gail Schmidt, USGS EROS
 Modified to output the UL and LR corner lat/longs.  We are already writing
   the bounding coords, however for ascending scenes and scenes in the polar
   regions, the scenes are flipped upside down.  The bounding coords will be
   correct in North represents the northernmost latitude and South represents
   the southernmost latitude.  However, the UL corner in this case would be
   more south than the LR corner.  Comparing the UL and LR corners will allow
   the user to determine if the scene is flipped.

 Revision 1.5 2013/08/05
 Gail Schmidt, USGS EROS
 Modified to output the gain and bias values for the reflectance and
   brightness temperature bands.

 Revision 2.0 2014/01/27
 Gail Schmidt, USGS EROS
 Modified application to utilize the ESPA internal raw binary format.

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
  
!END****************************************************************************
*/

#include <stdlib.h>
#include "output.h"
#include "const.h"
#include "names.h"
#include "error.h"
#include "cal.h"

Output_t *OpenOutput(Espa_internal_meta_t *in_meta, Input_t *input,
  Param_t *param, Lut_t *lut, bool thermal, int mss_flag)
/* 
!C******************************************************************************

!Description: 'OutputFile' sets up the 'output' data structure, opens the
 output file for write access, and creates the output Science Data Set (SDS).
 
!Input Parameters:
 in_meta        input XML metadata structure (band-related info)
 input          input structure with input image metadata (nband, iband, size)
 param          input paramter information (LEDAPS version)
 lut            lookup table (fill and saturation values)
 thermal        is this the thermal data?
 mss_flag       is this MSS data?

!Output Parameters:
 (returns)      'output' data structure or NULL when an error occurs

!Revision:
revision 1.0.0 9/12/2012
 - modified the application to write the thermal band QA bits to the output
   thermal product and use lndth_qa for the band name vs. lndcal_qa
revision 2.0.0 1/27/2014
 - modified to utilize the ESPA internal raw binary structure

!END****************************************************************************
*/
{
  Output_t *this = NULL;       /* pointer to output structure */
  char *mychar = NULL;         /* pointer to '_' */
  char scene_name[STR_SIZE];   /* scene name for the current scene */
  char rep_band[STR_SIZE];     /* representative band in the XML file */
  int ib;             /* looping variables */
  int nband;          /* number of bands for this dataset */
  int nband_tot;      /* number of total bands with refl/thermal and QA */
  int rep_indx=0;     /* band index in XML file for the current product */
  char production_date[MAX_DATE_LEN+1]; /* current date/time for production */
  time_t tp;          /* time structure */
  struct tm *tm;      /* time structure for UTC time */
  Espa_band_meta_t *bmeta = NULL;  /* pointer to the band metadata array
                         within the output structure */

  /* Determine the number of output bands */
  if (!thermal)
    nband = input->nband;
  else
    nband = input->nband_th;

  if (mss_flag == 1) 
    nband_tot = nband;
  else
    nband_tot = nband + NBAND_QA;

  /* Check parameters */
  if (input->size.l < 1)
    RETURN_ERROR("invalid number of output lines", "OpenOutput", NULL);

  if (input->size.s < 1)
    RETURN_ERROR("invalid number of samples per output line", "OpenOutput",
      NULL);

  if (nband < 1 || nband > NBAND_REFL_MAX)
    RETURN_ERROR("invalid number of bands", "OpenOutput", NULL);

  /* Create the Output data structure */
  this = (Output_t *) malloc (sizeof(Output_t));
  if (this == NULL) 
    RETURN_ERROR("allocating Output data structure", "OpenOutput", NULL);

  /* Find the representative band for metadata information. Reflective - band1.
     Thermal - band2. */
  if (!thermal)
    strcpy (rep_band, "band1");
  else
    strcpy (rep_band, "band6");
  for (ib = 0; ib < in_meta->nbands; ib++)
  {
    if (!strcmp (in_meta->band[ib].name, rep_band) &&
        !strncmp (in_meta->band[ib].product, "L1", 2))  /* L1G or L1T */
    {
      /* this is the index we'll use for band info from the XML strcuture */
      rep_indx = ib;
      break;
    }
  }

  /* Initialize the internal metadata for the output product. The global
     metadata won't be updated, however the band metadata will be updated
     and used later for appending to the original XML file. */
  init_metadata_struct (&this->metadata);

  /* Allocate memory for the total bands */
  if (allocate_band_metadata (&this->metadata, nband_tot) != SUCCESS)
    RETURN_ERROR("allocating band metadata", "OpenOutput", NULL);
  bmeta = this->metadata.band;

  /* Determine the scene name */
  strcpy (scene_name, in_meta->band[rep_indx].file_name);
  mychar = strrchr (scene_name, '_');
  if (mychar != NULL)
    *mychar = '\0';

  /* Get the current date/time (UTC) for the production date of each band */
  if (time (&tp) == -1)
    RETURN_ERROR("getting time", "OpenOutput", NULL);

  tm = gmtime (&tp);
  if (tm == NULL)
    RETURN_ERROR("converting time to UTC", "OpenOutput", NULL);

  if (strftime (production_date, MAX_DATE_LEN, "%Y-%m-%dT%H:%M:%SZ", tm) == 0)
    RETURN_ERROR("formating production date/time", "OpenOutput", NULL);

  /* Populate the data structure */
  this->open = false;
  this->nband = nband_tot;
  this->size.l = input->size.l;
  this->size.s = input->size.s;
  for (ib = 0; ib < nband_tot; ib++) {
    strncpy (bmeta[ib].short_name, in_meta->band[rep_indx].short_name,
      3);
    bmeta[ib].short_name[3] = '\0';
    if (!thermal)
    {
      strcpy (bmeta[ib].product, "toa_refl");
      strcat (bmeta[ib].short_name, "REF");
    }
    else
    {
      strcpy (bmeta[ib].product, "toa_bt");
      strcat (bmeta[ib].short_name, "BT");
    }
    bmeta[ib].nlines = this->size.l;
    bmeta[ib].nsamps = this->size.s;
    bmeta[ib].pixel_size[0] = in_meta->band[rep_indx].pixel_size[0];
    bmeta[ib].pixel_size[1] = in_meta->band[rep_indx].pixel_size[1];
    strcpy (bmeta[ib].pixel_units, "meters");
    sprintf (bmeta[ib].app_version, "LEDAPS_%s", param->LEDAPSVersion);
    strcpy (bmeta[ib].production_date, production_date);

    if (ib < nband)  /* image band */
    {
      bmeta[ib].data_type = ESPA_INT16;
      bmeta[ib].fill_value = lut->out_fill;
      bmeta[ib].saturate_value = lut->out_satu;
      if (!thermal)
      {
        sprintf (bmeta[ib].name, "toa_band%d", input->meta.iband[ib]);
        bmeta[ib].scale_factor = lut->scale_factor_ref;
        bmeta[ib].add_offset = lut->add_offset_ref;
        sprintf (bmeta[ib].long_name, "band %d TOA reflectance",
          input->meta.iband[ib]);
        strcpy (bmeta[ib].data_units, lut->units_ref);
        bmeta[ib].valid_range[0] = lut->valid_range_ref[0];
        bmeta[ib].valid_range[1] = lut->valid_range_ref[1];
        bmeta[ib].calibrated_nt = lut->calibrated_nt_ref;
      }
      else
      {
        sprintf (bmeta[ib].name, "toa_band%d", input->meta.iband_th);
        bmeta[ib].scale_factor = lut->scale_factor_ref;
        bmeta[ib].scale_factor = lut->scale_factor_th;
        bmeta[ib].add_offset = lut->add_offset_th;
        sprintf (bmeta[ib].long_name, "band %d brightness temperature",
          input->meta.iband_th);
        strcpy (bmeta[ib].data_units, lut->units_th);
        bmeta[ib].valid_range[0] = lut->valid_range_th[0];
        bmeta[ib].valid_range[1] = lut->valid_range_th[1];
        bmeta[ib].calibrated_nt = lut->calibrated_nt_th;
      }
    }
    else  /* QA band */
    {
      if (!thermal)
      {
        /* Set up QA bitmap information */
        strcpy (bmeta[ib].name, "toa_qa");
        if (allocate_bitmap_metadata (&bmeta[ib], 8) != SUCCESS)
          RETURN_ERROR("allocating 8 bits for the bitmap", "OpenOutput", NULL); 

        bmeta[ib].bitmap_description[0] =
          "Data Fill Flag (0 = valid data, 1 = invalid data)";
	    strcpy (bmeta[ib].bitmap_description[1],
          "Band 1 Data Saturation Flag (0 = valid data, 1 = saturated data)");
	    strcpy (bmeta[ib].bitmap_description[2],
          "Band 2 Data Saturation Flag (0 = valid data, 1 = saturated data)");
	    strcpy (bmeta[ib].bitmap_description[3],
          "Band 3 Data Saturation Flag (0 = valid data, 1 = saturated data)");
	    strcpy (bmeta[ib].bitmap_description[4],
          "Band 4 Data Saturation Flag (0 = valid data, 1 = saturated data)");
	    strcpy (bmeta[ib].bitmap_description[5],
          "Band 5 Data Saturation Flag (0 = valid data, 1 = saturated data)");
	    strcpy (bmeta[ib].bitmap_description[6],
          "Band 6 Data Saturation Flag (not set)");
	    strcpy (bmeta[ib].bitmap_description[7],
          "Band 7 Data Saturation Flag (0 = valid data, 1 = saturated data)");
      }
      else
      {
        strcpy (bmeta[ib].name, "toa_bt_qa");
        if (allocate_bitmap_metadata (&bmeta[ib], 8) != SUCCESS)
          RETURN_ERROR("allocating 8 bits for the bitmap", "OpenOutput", NULL); 
        bmeta[ib].bitmap_description[0] =
          "Data Fill Flag (0 = valid data, 1 = invalid data)";
	    strcpy (bmeta[ib].bitmap_description[1],
          "Band 1 Data Saturation Flag (not set)");
	    strcpy (bmeta[ib].bitmap_description[2],
          "Band 2 Data Saturation Flag (not set)");
	    strcpy (bmeta[ib].bitmap_description[3],
          "Band 3 Data Saturation Flag (not set)");
	    strcpy (bmeta[ib].bitmap_description[4],
          "Band 4 Data Saturation Flag (not set)");
	    strcpy (bmeta[ib].bitmap_description[5],
          "Band 5 Data Saturation Flag (not set)");
	    strcpy (bmeta[ib].bitmap_description[6],
          "Band 6 Data Saturation Flag (0 = valid data, 1 = saturated data)");
	    strcpy (bmeta[ib].bitmap_description[7],
          "Band 7 Data Saturation Flag (not set)");
      }
      bmeta[ib].data_type = ESPA_UINT8;
      bmeta[ib].fill_value = lut->qa_fill;
      strcpy (bmeta[ib].long_name, "QA band");
      strcpy (bmeta[ib].data_units, "bitmap");
      bmeta[ib].valid_range[0] = 0;
      bmeta[ib].valid_range[1] = 255;
    }

    /* Set up the filename with the scene name and band name and open the
       file for write access */
    sprintf (bmeta[ib].file_name, "%s_%s.img", scene_name,
      bmeta[ib].name);
    this->fp_bin[ib] = open_raw_binary (bmeta[ib].file_name, "w");
    if (this->fp_bin[ib] == NULL)
      RETURN_ERROR("unable to open output band file", "OpenOutput", NULL);
  }  /* for ib */
  this->open = true;

  /* Successful completion */
  return this;
}


bool CloseOutput(Output_t *this)
/* 
!C******************************************************************************

!Description: 'CloseOutput' ends SDS access and closes the output file.
 
!Input Parameters:
 this           'output' data structure

!Output Parameters:
 this           'output' data structure; the following fields are modified:
                   open
 (returns)      status:
                  'true' = okay
                  'false' = error return

!END****************************************************************************
*/
{
  int ib;

  if (!this->open)
    RETURN_ERROR("image files not open", "CloseOutput", false);

  for (ib = 0; ib < this->nband; ib++)
    close_raw_binary (this->fp_bin[ib]);

  this->open = false;
  return true;
}


bool FreeOutput(Output_t *this)
/* 
!C******************************************************************************

!Description: 'FreeOutput' frees the 'output' data structure memory.
 
!Input Parameters:
 this           'output' data structure for which the fields are freed

!Output Parameters:
 this           'output' data structure
 (returns)      status:
                  'true' = okay
                  'false' = error occurred

!END****************************************************************************
*/
{
  if (this->open) 
    RETURN_ERROR("file still open", "FreeOutput", false);

  free(this);
  this = NULL;

  return true;
}


bool PutOutputLine(Output_t *this, int iband, int iline, void *line)
/* 
!C******************************************************************************

!Description: 'WriteOutput' writes a line of data to the output file.
 
!Input Parameters:
 this           'output' data structure
 iband          index (within Output_t struct) of output band to be written
 iline          output line number (used for validation only)
 line           buffer of data to be written

!Output Parameters:
 (returns)      status:
                  'true' = okay
                  'false' = error return

!END****************************************************************************
*/
{
  int nbytes = 0;      /* number of bytes in each pixel */
  Espa_band_meta_t *bmeta = this->metadata.band;  /* pointer to band metadata */

  /* Check the parameters */
  if (this == NULL) 
    RETURN_ERROR("invalid input structure", "PutOutputLine", false);
  if (!this->open)
    RETURN_ERROR("file not open", "PutOutputLine", false);
  if (iband < 0 || iband >= this->nband)
    RETURN_ERROR("invalid band number", "PutOutputLine", false);
  if (iline < 0 || iline >= this->size.l)
    RETURN_ERROR("invalid line number", "PutOutputLine", false);

  /* Write the data, only the current line (i.e. one line at a time) */
  if (bmeta[iband].data_type == ESPA_INT16)
    nbytes = sizeof (int16);
  else
    nbytes = sizeof (unsigned char);
  if (write_raw_binary (this->fp_bin[iband], 1, this->size.s, nbytes, line)
      != SUCCESS)
    RETURN_ERROR("writing output line", "PutOutputLine", false);

  return true;
}
