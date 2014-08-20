#ifndef _L8_SR_H_
#define _L8_SR_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "common.h"
#include "input.h"
#include "output.h"
#include "lut_subr.h"
#include "espa_metadata.h"
#include "espa_geoloc.h"
#include "parse_metadata.h"
#include "write_metadata.h"
#include "envi_header.h"
#include "error_handler.h"

/* Prototypes */
void usage ();

int get_args
(
    int argc,             /* I: number of cmd-line args */
    char *argv[],         /* I: string of cmd-line args */
    char **xml_infile,    /* O: address of input XML file */
    char **aux_infile,    /* O: address of input auxiliary file containing
                                water vapor and ozone */
    bool *process_sr,     /* O: process the surface reflectance products */
    bool *write_toa,      /* O: write intermediate TOA products flag */
    bool *verbose         /* O: verbose flag */
);

void usage ();

bool btest
(
    uint8 byte_val,   /* I: byte value to be tested with the bit n */
    byte n            /* I: bit number to be tested (0 is rightmost bit) */
);


#endif
