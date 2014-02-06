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

 Revision 1.5 2012/12/27
 Gail Schmidt, USGS/EROS
 Added support for reading the metadata filename from the parameter file

 Revision 1.6 2013/01/22
 Gail Schmidt, USGS/EROS
 Modified to utilize only one Version value - LEDAPSVersion

 Revision 2.0 02/03/2014 Gail Schmidt, USGS EROS
 Modified applications to use the ESPA internal raw binary file format.

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

#include "param.h"
#include "mystring.h"
#include "error.h"

typedef enum {
  PARAM_NULL = -1,
  PARAM_START = 0,
  PARAM_XML_FILE,
  PARAM_NCEP_FILE,
  PARAM_PRWV_FILE,
  PARAM_OZON_FILE,
  PARAM_DEM_FILE,
  PARAM_LEDAPSVERSION,
  PARAM_END,
  PARAM_MAX
} Param_key_t;

Key_string_t Param_string[PARAM_MAX] = {
  {(int)PARAM_START,     "PARAMETER_FILE"},
  {(int)PARAM_XML_FILE,  "XML_FILE"},
  {(int)PARAM_NCEP_FILE, "NCEP_FIL"},
  {(int)PARAM_PRWV_FILE, "PRWV_FIL"},
  {(int)PARAM_OZON_FILE, "OZON_FIL"},
  {(int)PARAM_DEM_FILE,  "DEM_FILE"},
  {(int)PARAM_LEDAPSVERSION,  "LEDAPSVersion"},
  {(int)PARAM_END,       "END"}
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
  char *error_string = (char *)NULL;
  FILE *fp;
  Key_t key;
  int i,len;
  char line[MAX_STR_LEN + 1];
  char temp[MAX_STR_LEN + 1];
  Param_key_t param_key;
  char *param_file_name;
  bool got_start, got_end;

  if (argc < 2) 
    RETURN_ERROR("no command line parameter", "GetParam", NULL);
  if (argc > 2) 
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
  this->param_file_name  = NULL;
  this->input_xml_file_name  = NULL;
  this->LEDAPSVersion   = NULL;
  this->num_ncep_files   = 0;            /* number of NCEP files     */
  this->num_prwv_files   = 0;            /* number of PRWV hdf files */
  this->num_ozon_files   = 0;            /* number of OZONe hdf files */
  this->dem_file = NULL;
  this->dem_flag = false;
  this->thermal_band=false;

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

      case PARAM_NCEP_FILE:
        this->num_ncep_files   = key.nval;            
        if (key.nval > 4) {
          error_string = "too many NCEP file names";
          break; 
        }
        for (i=0;i<key.nval;i++) {
          if (key.len_value[i] < 1) {
           error_string = "no NCEP file name";
           break;
          }
          key.value[i][key.len_value[i]] = '\0';
          this->ncep_file_name[i] = DupString(key.value[i]);
          if (this->ncep_file_name[i] == NULL) {
          error_string = "duplicating NCEP file name";
            break;
          }
        }
        break;

      case PARAM_PRWV_FILE:
        this->num_prwv_files   = key.nval;       
        if (key.nval > 1) {
          error_string = "too many PRWV file names";
          break; 
        }
        if (key.nval > 0) {
          if (key.len_value[0] < 1) {
            error_string = "no PRWV hdf file name";
            break;
          }
          key.value[0][key.len_value[0]] = '\0';
          this->prwv_file_name = DupString(key.value[0]);
          if (this->prwv_file_name == NULL) {
            error_string = "duplicating PRWV hdf file name";
            break;
          }
        }
        break;

      case PARAM_OZON_FILE:
        this->num_ozon_files   = key.nval;       
        if (key.nval > 1) {
          error_string = "too many OZON file names";
          break; 
        }
        if (key.nval > 0) {
          if (key.len_value[0] < 1) {
            error_string = "no OZON hdf file name";
            break;
          }
          key.value[0][key.len_value[0]] = '\0';
          this->ozon_file_name = DupString(key.value[0]);
          if (this->ozon_file_name == NULL) {
            error_string = "duplicating OZON hdf file name";
            break;
          }
        }
        break;

      case PARAM_DEM_FILE:
        this->dem_flag = true;
        if (key.nval <= 0) {
          this->dem_flag = false;
          break; 
        } else if (key.nval > 1) {
          error_string = "too many DEM file names";
          break; 
        }
        if (key.len_value[0] < 1) {
          this->dem_flag = false;
          break;
        }
        key.value[0][key.len_value[0]] = '\0';
        this->dem_file  = DupString(key.value[0]);
        if (this->dem_file == NULL) {
          error_string = "duplicating dem file name";
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
  }
  return true;
}
