/*
!C****************************************************************************

!File: param.c
  
!Description: Functions for accepting parameters from the command line or 
 a file.

!Revision History:
 Revision 1.0 2000/11/07
 Robert Wolfe
 Original Version.

 Revision 1.1 2000/12/13
 Sadashiva Devadiga
 Modified to accept parameters from command line or file.

 Revision 1.2 2001/05/08
 Sadashiva Devadiga
 Added checks for required parameters.

 Revision 1.3 2002/03/02
 Robert Wolfe
 Added special handling for input ISINUS case.

 Revision 1.4 2002/05/10
 Robert Wolfe
 Added separate output SDS name.

 Revision 1.5 01/22/2013 Gail Schmidt, USGS EROS
 Modified applications to use only one version and that is the
 LEDAPSVersion tag which will get updated with each release of LEDAPS

 Revision 2.0 01/21/2014 Gail Schmidt, USGS EROS
 Modified applications to use the ESPA internal raw binary file format.

 Revision 03/31/2015 Gail Schmidt, USGS EROS
 Added an existance check for the TOA reflectance (and K1/K2 consts).

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

      Sadashiva Devadiga (Code 922)
      MODIS Land Team Support Group     SSAI
      devadiga@ltpmail.gsfc.nasa.gov    5900 Princess Garden Pkway, #300
      phone: 301-614-5549               Lanham, MD 20706
  
 ! Design Notes:
   1. The following public functions handle the input data:

	GetParam - Setup 'param' data structure and populate with user
	           parameters.
	FreeParam - Free the 'param' data structure memory.

   2. 'GetParam' must be called before 'FreeParam'.  
   3. 'FreeParam' should be used to free the 'param' data structure.

!END****************************************************************************
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "lndcal.h"
#include "param.h"
#include "mystring.h"
#include "util.h"
#include "error.h"

typedef enum {
  PARAM_NULL = -1,
  PARAM_START = 0,
  PARAM_XML_FILE,
  PARAM_LEDAPSVERSION,
  PARAM_END,
  PARAM_MAX
} Param_key_t;

Key_string_t Param_string[PARAM_MAX] = {
  {(int)PARAM_START,       "PARAMETER_FILE"},
  {(int)PARAM_XML_FILE,    "XML_FILE"},
  {(int)PARAM_LEDAPSVERSION,  "LEDAPSVersion"},
  {(int)PARAM_END,         "END"}
};

/* Functions */

Param_t *GetParam(int argc, const char **argv)
/* 
!C******************************************************************************

!Description: 'GetParam' sets up the 'param' data structure and populate with user
 parameters, either from the command line or from a parameter file.
 
!Input Parameters:
 argc           number of command line arguments
 argv           command line argument list

!Output Parameters:
 (returns)      'param' data structure or NULL when an error occurs

!Team Unique Header:

!END****************************************************************************
*/
{
  Param_t *this;
  char *error_string = NULL;
  FILE *fp;
  Key_t key;
  int len;
  char line[MAX_STR_LEN + 1];
  char temp[MAX_STR_LEN + 1];
  Param_key_t param_key;
  char *param_file_name;
  bool got_start, got_end;

  if (argc < 2) 
    RETURN_ERROR("no command line parameter", "GetParam", NULL);
  if (argc > 3) 
    RETURN_ERROR("too many command line parameters", "GetParam", NULL);
  if (strlen(argv[1]) < 1)
    RETURN_ERROR("no parameter file name", "GetParam", NULL);
  param_file_name = (char *)argv[1];

  /* Open the parameter file */
  if ((fp = fopen(param_file_name, "r")) == NULL)
    RETURN_ERROR("unable to open parameter file", "GetParam", NULL);

  /* Create the Param data structure */
  this = (Param_t *)malloc(sizeof(Param_t));
  if (this == NULL) {
    fclose(fp);
    RETURN_ERROR("allocating Input structure", "GetParam", NULL);
  }

  /* set default parameters */
  this->param_file_name         = NULL;
  this->input_xml_file_name     = NULL;
  this->LEDAPSVersion           = NULL;

  /* Populate the data structure */
  this->param_file_name = DupString(param_file_name);
  if (this->param_file_name == NULL)
    error_string = "duplicating parameter file name";

  if (error_string != NULL) {
    free(this->param_file_name);
    this->param_file_name = NULL;
    FreeParam(this);
    RETURN_ERROR(error_string, "GetParam", NULL);
  }

  /* Parse the header file */
  got_start = got_end = false;

  while((len = GetLine(fp, line)) > 0) {

    if (!StringParse(line, &key)) {
      sprintf(temp, "parsing header file; line = %s", line);
      error_string = temp;
      break;
    }
    if (key.len_key <= 0) continue;
    if (key.key[0] == '#') continue;

    param_key = (Param_key_t) KeyString(key.key, key.len_key, Param_string, 
       (int)PARAM_NULL, (int)PARAM_MAX);
    if (param_key == PARAM_NULL) {
      key.key[key.len_key] = '\0';
      sprintf(temp, "invalid key; key = %s", key.key);
      error_string = temp;
      break;
    }
    if (!got_start) {
      if (param_key == PARAM_START) {
        if (key.nval != 0) {
          error_string = "no value expected (start key)";
          break;
        }
        got_start = true;
        continue;
      } else {
        error_string  = "no start key in parameter file";
        break;
      }
    }

    /* Get the value for each keyword */
    switch (param_key) {

      case PARAM_XML_FILE:
        if (key.nval <= 0) {
          error_string = "no input XML metadata file name";
          break; 
        } else if (key.nval > 1) {
          error_string = "too many input XML metadata file names";
          break; 
        }
        if (key.len_value[0] < 1) {
          error_string = "no input XML metadata file name";
          break;
        }
        key.value[0][key.len_value[0]] = '\0';
        this->input_xml_file_name = DupString(key.value[0]);
        if (this->input_xml_file_name == NULL) {
          error_string = "duplicating input XML metadata file name";
          break;
        }
        break;

      case PARAM_LEDAPSVERSION:
        if (key.nval <= 0) {
          error_string = "no LEDAPSVersion number";
          break;
        } else if (key.nval > 1) {
          error_string = "too many LEDAPSVersion numbers";
          break;
        }
        if (key.len_value[0] < 1) {
          error_string = "no LEDAPSVersion number";
          break;
        }
        key.value[0][key.len_value[0]] = '\0';
        this->LEDAPSVersion = DupString(key.value[0]);
        if (this->LEDAPSVersion == NULL) {
          error_string = "duplicating LEDAPSVersion number";
          break;
        }
        break;

      case PARAM_END:
        if (key.nval != 0) {
          error_string = "no value expected (end key)";
          break; 
        }
        got_end = true;
        break;

      default:
        error_string = "key not implmented";
    }
    if (error_string != NULL) break;
    if (got_end) break;
  }

  /* Close the parameter file */
  fclose(fp);

  if (error_string == NULL) {
    if (!got_start) 
      error_string = "no start key in header";
    else if (!got_end)
      error_string = "no end key in header";
  }

  /* Handle null values */
  if (error_string == NULL) {
    if (this->input_xml_file_name == NULL) 
      error_string = "no input XML metadata file name given";
    if (this->LEDAPSVersion == NULL)
      error_string = "no LEDAPS Version given";
  }

  /* Handle errors */
  if (error_string != NULL) {
    free(this->param_file_name);
    free(this->input_xml_file_name);
    free(this->LEDAPSVersion);
    free(this);
    RETURN_ERROR(error_string, "GetParam", NULL);
  }
  
  return this;
}


bool FreeParam(Param_t *this)
/* 
!C******************************************************************************

!Description: 'FreeParam' frees the 'param' data structure memory.
 
!Input Parameters:
 this           'param' data structure

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)

!Team Unique Header:

 ! Design Notes:
   1. 'GetParam' must be called before this routine is called.

!END****************************************************************************
*/
{
  if (this != NULL) {
    free(this->param_file_name);
    free(this->input_xml_file_name);
    free(this);
    this = NULL;
  }
  return true;
} 
  
bool existRadGB(Espa_internal_meta_t *metadata)
/* 
!C******************************************************************************

!Description: 'existRadGB' determines if the gains and biases for TOA radiance
conversion exist and were set from the input XML file.
 
!Input Parameters:
 metadata     'Espa_internal_meta_t' data structure

!Output Parameters:
 (returns)      N/A

!Team Unique Header:

!END****************************************************************************
*/
{
  int refl_indx = 0;  /* index of band1 in the input file */
  int i;              /* looping variable */

  /* Find band1 in the input XML file */
  for (i = 0; i < metadata->nbands; i++)
  {
    if (!strcmp (metadata->band[i].name, "band1") &&
        !strncmp (metadata->band[i].product, "L1", 2))  /* L1G or L1T */
    {
      /* this is the index we'll use for reflectance band info */
      refl_indx = i;
    }
  }

  /* If the gain or bias for band1 in the input file is not set, then assume
     none are set and therefore need to be manually set before continuing. */
  if (fabs (metadata->band[refl_indx].rad_gain - ESPA_FLOAT_META_FILL) <
        ESPA_EPSILON ||
      fabs (metadata->band[refl_indx].rad_bias - ESPA_FLOAT_META_FILL) <
        ESPA_EPSILON)
  {
    return false;
  }

  return true;
}
  
bool existReflGB(Espa_internal_meta_t *metadata)
/* 
!C******************************************************************************

!Description: 'existReflGB' determines if the gains and biases for TOA
reflectance conversion exist and were set from the input XML file.  If these
conversion parameters exist, then it is also assumed the K1 and K2 thermal
constants were set as well.
 
!Input Parameters:
 metadata     'Espa_internal_meta_t' data structure

!Output Parameters:
 (returns)      N/A

!Team Unique Header:

!END****************************************************************************
*/
{
  int refl_indx = 0;  /* index of band1 in the input file */
  int i;              /* looping variable */

  /* Find band1 in the input XML file */
  for (i = 0; i < metadata->nbands; i++)
  {
    if (!strcmp (metadata->band[i].name, "band1") &&
        !strncmp (metadata->band[i].product, "L1", 2))  /* L1G or L1T */
    {
      /* this is the index we'll use for reflectance band info */
      refl_indx = i;
    }
  }

  /* If the gain or bias for band1 in the input file is not set, then assume
     none are set and therefore need to be manually set before continuing. */
  if (fabs (metadata->band[refl_indx].refl_gain - ESPA_FLOAT_META_FILL) <
        ESPA_EPSILON ||
      fabs (metadata->band[refl_indx].refl_bias - ESPA_FLOAT_META_FILL) <
        ESPA_EPSILON)
  {
    return false;
  }

  return true;
}
