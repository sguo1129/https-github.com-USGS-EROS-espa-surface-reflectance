#include <stdio.h> 
#include <stdlib.h> 
#include <ctype.h>
#include "espa_metadata.h"
#include "parse_metadata.h"

#define D2R     1.745329251994328e-2

int main(int argc, char **argv)
{
  char projection[256];
  float coordinates[8];
  double parm[13];
  double radius, lat, lon, dl, ds;
  double corner[2];
  int ret;
  int zonecode, sphercode, rows, cols;
  float orientationangle, pixelsize, upperleftx, upperlefty;
  double arg2, arg3;
  char *error_file = "geo_xy.ERROR";
  FILE *error_ptr=NULL;
  
  int LSsphdz(char *projection, float coordinates[8], double *parm, double *radius, double corner[2]);
  int LSutminv(double s, double l, double *lon, double *lat);
  int LSutmfor(double *s, double *l, double lon, double lat);
  int LSpsinv(double s, double l, double *lon, double *lat);
  int LSpsfor(double *s, double *l, double lon, double lat);
  int get_data(char *filename, char *projection, int *zonecode, int *sphercode,
    float *orientationangle, float *pixelsize, float *upperleftx,
    float *upperlefty, int *rows, int *cols, double *projparms);
  
  if (argc < 4) {
  #ifdef INV
     printf("usage: %s <XML file> <sample> <line>\n", argv[0]);
  #else
     printf("usage: %s <XML file> <longitude> <latitude>\n", argv[0]);
  #endif
     printf("Jim Ray, SSAI, %s\n\n", __DATE__);
     exit(0);
  }
  
  if ( (ret = get_data(argv[1], projection, &zonecode, &sphercode,
      &orientationangle, &pixelsize, &upperleftx, &upperlefty, &rows, &cols,
      parm)) != 0) {
     printf("Error reading file %s, cannot continue\n", argv[1]);
     error_ptr = fopen (error_file, "w");
     fprintf(error_ptr, "Error reading file %s, cannot continue\n", argv[1]);
     fclose (error_ptr);
     exit(1);
  } 
  
  /* if processing PS projection, then convert the angular projection params
     to radians */
  if (!strcmp (projection, "GCTP_PS")) {
      parm[4] *= D2R;
      parm[5] *= D2R;
  }
  
  arg2 = atof(argv[2]);
  arg3 = atof(argv[3]);
  
  coordinates[4] = (double)zonecode;
  coordinates[5] = (double)sphercode;
  coordinates[6] = (double)orientationangle;
  coordinates[7] = (double)pixelsize;
  corner[0] = (double)upperleftx;
  corner[1] = (double)upperlefty;
  
  LSsphdz(projection, coordinates, parm, &radius, corner);
  
  #ifdef INV
  ds = arg2;
  dl = arg3;
  
  if (ds > (double)cols) {
     printf("Sample argument (%s) exceeds number of columns in file (%d): will "
         "use %d\n", argv[2], cols, cols);
     ds = (double)cols;
  }
  
  if (dl > (double)rows) {
     printf("Sample argument (%s) exceeds number of rows in file (%d): will "
         "use %d\n", argv[3], rows, rows);
     dl = (double)rows;   
  }
  
  if (!strcmp (projection, "GCTP_UTM"))
     ret = LSutminv(ds, dl, &lon, &lat);
  else if (!strcmp (projection, "GCTP_PS"))
     ret = LSpsinv(ds, dl, &lon, &lat);
  printf("line   %5.1f  samp   %5.1f  => long %f lat %f\n", dl, ds, lon, lat);
  #else
  lon = arg2;
  lat = arg3;
  
  /* We need sanity checks on these as well */
  
  if (!strcmp (projection, "GCTP_UTM"))
      ret = LSutmfor(&ds, &dl, lon, lat);
  else if (!strcmp (projection, "GCTP_PS"))
      ret = LSpsfor(&ds, &dl, lon, lat);
  printf("long %f lat %f => line   %f  samp   %f  \n",  lon, lat, dl, ds);
  #endif
  
  exit (1);
}


/******************************************************************************
Module: get_data

Description: Reads the metadata from the input XML file

Inputs:
  filename:  name of XML file

Outputs:
  projection:     projection name (GCTP_UTM, GCTP_UTM, etc.)
  zonecode:       UTM zone number
  spherecode:     spheroid number
  orientationangle:  orientation of the scene (degrees)
  pixelsize:       size of each pixel (assumed square)
  upperleft[x/y]: UL x,y corner point (meters)
  rows:           number of lines in the scene
  cols:           number of samples in the scene
  projparms:      array of 13 projection parameters for the projection

History:
  2/6/2014  Gail Schmidt, USGS/EROS
  Modified to use the ESPA internal file format
******************************************************************************/
int get_data(char *filename, char *projection, int *zonecode, int *sphercode,
  float *orientationangle, float *pixelsize, float *upperleftx,
  float *upperlefty, int *rows, int *cols, double *projparms)
{
  int i;              /* looping variable */
  int ib;             /* band looping variable */
  int rep_indx=-1;    /* band index in XML file for the current product */
  Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
  Espa_global_meta_t *gmeta = NULL;   /* pointer to global metadata */
  Espa_band_meta_t *bmeta = NULL;     /* pointer to the band metadata array
                                         within the output structure */

  /* Initialize the outputs */
  *zonecode = *sphercode = *rows = *cols = -1;
  *orientationangle = *pixelsize = -999.0;
  for (i = 0; i < 13; i++)
    projparms[i] = 0.0;

  /* Validate the input metadata file */
  if (validate_xml_file (filename, ESPA_SCHEMA) != SUCCESS)
  {  /* Error messages already written */
    printf("Error validating XML file: %s", filename);
    return (-2);
  }

  /* Initialize the metadata structure */
  init_metadata_struct (&xml_metadata);

  /* Parse the metadata file into our internal metadata structure; also
     allocates space as needed for various pointers in the global and band
     metadata */
  if (parse_metadata (filename, &xml_metadata) != SUCCESS)
  {  /* Error messages already written */
    printf("Error parsing XML file: %s", filename);
    return (-4);
  }
  gmeta = &xml_metadata.global;

  /* Look for band1 in the TOA product and use for our representative band */
  for (ib = 0; ib < xml_metadata.nbands; ib++)
  {
    if (!strcmp (xml_metadata.band[ib].name, "toa_band1") &&
        !strcmp (xml_metadata.band[ib].product, "toa_refl"))
    {
      /* this is the index we'll use for band info from the XML strcuture */
      rep_indx = ib;
      break;
    }
  }
  if (rep_indx == -1)
  {
    printf("Error finding toa_band1 band in the XML file");
    return (-5);
  }
  bmeta = &xml_metadata.band[rep_indx];

  /* Pull the projection and key metadata information from the XML file. For
     the UL corner make sure to addjust the center of the pixel appropriately
     to produce coords for the UL of the pixel. */
  *sphercode = gmeta->proj_info.sphere_code;
  if (gmeta->proj_info.proj_type == GCTP_UTM_PROJ)
  {
    strcpy (projection, "GCTP_UTM");
    *zonecode = gmeta->proj_info.utm_zone;
  }
  else if (gmeta->proj_info.proj_type == GCTP_PS_PROJ)
  {
    strcpy (projection, "GCTP_PS");
    projparms[4] = gmeta->proj_info.longitude_pole;
    projparms[5] = gmeta->proj_info.latitude_true_scale;
    projparms[6] = gmeta->proj_info.false_easting;
    projparms[7] = gmeta->proj_info.false_northing;

  }
  else
  {
    printf ("Error in projection code. Only GCTP_UTM and GCTP_PS are currently "
      "supported.\n");
    return (-5);
  } 

  *orientationangle = gmeta->orientation_angle;
  *upperleftx = gmeta->proj_info.ul_corner[0];
  *upperlefty = gmeta->proj_info.ul_corner[1];
  *pixelsize = bmeta->pixel_size[0];
  if (!strcmp (gmeta->proj_info.grid_origin, "center"))
  { /* adjust by pixel size */
    *upperleftx -= bmeta->pixel_size[0];
    *upperlefty += bmeta->pixel_size[1];
  }
  *rows = bmeta->nlines;
  *cols = bmeta->nsamps;

  if (!strcmp (projection, "GCTP_UTM") && (*zonecode == -1))
  {
      printf("ERROR reading UTM zone code, cannot continue...\n");
      return (-5);
  } 

  return (0);
}

