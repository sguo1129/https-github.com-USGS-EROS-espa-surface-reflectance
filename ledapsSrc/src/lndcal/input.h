/*
!C****************************************************************************

!File: input.h

!Description: Header file for 'input.c' - see 'input.c' for more information.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

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
   1. Structure is declared for the 'input' data type.
  
!END****************************************************************************
*/

#ifndef INPUT_H
#define INPUT_H

#include <stdlib.h>
#include <stdio.h>
#include "lndcal.h"
#include "bool.h"
#include "const.h"
#include "date.h"
#include "param.h"

#define ANGLE_FILL (-999.0)
#define WRS_FILL (-1)
#define GAIN_BIAS_FILL (-999.0)

/* Input file type definition */

typedef enum {
  INPUT_TYPE_NULL = -1,
  INPUT_TYPE_BINARY = 0, 
  INPUT_TYPE_MAX
} Input_type_t;

/* Structure for the 'input metadata' data type */

typedef struct {
  Sat_t sat;               /* Satellite */
  Inst_t inst;             /* Instrument */
  Date_t acq_date;         /* Acqsition date/time (scene center) */
  bool time_fill;          /* Acqsition time fill; true = fill value (0h) */
  Date_t prod_date;        /* Production date (must be available for ETM) */
  float sun_zen;           /* Solar zenith angle (radians; scene center) */
  float sun_az;            /* Solar azimuth angle (radians; scene center) */
  float earth_sun_dist;    /* Earth-sun distance */
  Wrs_t wrs_sys;           /* WRS system */
  int ipath;               /* WRS path number */
  int irow;                /* WRS row number */
  unsigned char fill;      /* Fill value */
  int iband[NBAND_REFL_MAX]; /* Band numbers */
  int iband_th;            /* Thermal Band number= (6) */
  float rad_gain[NBAND_REFL_MAX]; /* TOA radiance band gain */
  float rad_bias[NBAND_REFL_MAX]; /* TOA radiance band bias */
  float rad_gain_th;           /* Thermal TOA radiance band gain */
  float rad_bias_th;           /* Thermal TOA radiance band bias */
  bool use_toa_refl_consts;    /* Are the TOA reflectance gain/bias and K1/K2
                                  constants available? Same with earth-sun
                                  distance */
  float refl_gain[NBAND_REFL_MAX]; /* TOA reflectance band gain */
  float refl_bias[NBAND_REFL_MAX]; /* TOA reflectance band bias */
  float k1_const;          /* K1 thermal constant */
  float k2_const;          /* K2 thermal constant */
} Input_meta_t;

/* Structure for the 'input' data type */

typedef struct {
  Input_type_t file_type;  /* Type of the input image files */
  Input_meta_t meta;       /* Input metadata */
  int nband;               /* Number of input image files (bands) */
  int nband_th;            /* Number of thermal input image files (0 or 1) */
  Img_coord_int_t size;    /* Input file size */
  Img_coord_int_t size_th; /* Input (thermal) file size */
  char *file_name[NBAND_REFL_MAX];  
                           /* Name of the input image files */
  char *file_name_th;      /* Name of the thermal input image files */
  bool open[NBAND_REFL_MAX]; 
                           /* Flag to indicate whether the specific input file 
			      is open for access; 'true' = open, 
			     'false' = not open */
  bool open_th;            /* thermal open flag */
  FILE *fp_bin[NBAND_REFL_MAX];  /* File pointer for binary files */
  FILE *fp_bin_th;         /* File pointer for thermal binary file */
} Input_t;

/* Prototypes */

Input_t *OpenInput(Espa_internal_meta_t *metadata);
bool GetInputLine(Input_t *this, int iband, int iline, unsigned char *line);
bool GetInputLineTh(Input_t *this, int iline, unsigned char *line);
bool CloseInput(Input_t *this);
bool FreeInput(Input_t *this);
bool InputMetaCopy(Input_meta_t *this, int nband, Input_meta_t *copy);
bool GetXMLInput(Input_t *this, Espa_internal_meta_t *metadata);

#endif
