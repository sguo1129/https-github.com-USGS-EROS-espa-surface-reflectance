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
#include "lndsr.h"
#include "const.h"
#include "date.h"

#define ANGLE_FILL -999.0
#define WRS_FILL -1

typedef struct {
  Sat_t sat;               /* Satellite */
  Inst_t inst;             /* Instrument */
  Date_t acq_date;         /* Acqsition date/time (scene center) */
  float sun_zen;           /* Solar zenith angle (radians; scene center) */
  float sun_az;            /* Solar azimuth angle (radians; scene center) */
  Wrs_t wrs_sys;           /* WRS system */
  int ipath;               /* WRS path number */
  int irow;                /* WRS row number */
  int fill;                /* Fill value */
  int iband[NBAND_REFL_MAX]; /* Reflectance band numbers */
  int iband_qa;            /* QA band number */
} Input_meta_t;

/* Structure for the 'input' data type */

typedef struct {
  Input_meta_t meta;       /* Input metadata */
  int nband;               /* Number of input image bands */
  Img_coord_int_t size;    /* Input file size */
  char *file_name[NBAND_REFL_MAX];  /* Name of the input image files */
  char *file_name_qa;      /* Name of the input QA file */
  FILE *fp_bin[NBAND_REFL_MAX];  /* File pointer for input binary files */
  bool open[NBAND_REFL_MAX]; /* Flag to indicate whether the specific input
                                file is open for access; 'true' = open, 
                                'false' = not open */
  FILE *fp_bin_qa;         /* File pointer for QA binary file */
  bool open_qa;            /* Flag to indicate whether the specific input
                              file is open for access; 'true' = open, 
                              'false' = not open */
} Input_t;

/* Prototypes */

Input_t *OpenInput(Espa_internal_meta_t *metadata, bool thermal);
bool GetInputLine(Input_t *this, int iband, int iline, int16 *line);
bool CloseInput(Input_t *this);
bool FreeInput(Input_t *this);
bool InputMetaCopy(Input_meta_t *this, int nband, Input_meta_t *copy);
bool GetXMLInput(Input_t *this, Espa_internal_meta_t *metadata, bool thermal);
bool GetInputQALine(Input_t *this, int iline, uint8 *line);

#endif
