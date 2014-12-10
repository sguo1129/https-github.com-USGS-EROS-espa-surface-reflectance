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

int compute_toa_refl
(
    Input_t *input,     /* I: input structure for the Landsat product */
    uint16 *qaband,     /* I: QA band for the input image, nlines x nsamps */
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    float xmus,         /* I: cosine of solar zenith angle */
    char *instrument,   /* I: instrument to be processed (OLI, TIRS) */
    int16 **sband       /* O: output surface reflectance and brightness
                              temp bands */
);

int compute_sr_refl
(
    Input_t *input,     /* I: input structure for the Landsat product */
    Espa_internal_meta_t *xml_metadata,
                        /* I: XML metadata structure */
    char *xml_infile,   /* I: input XML filename */
    uint16 *qaband,     /* I: QA band for the input image, nlines x nsamps */
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    float pixsize,      /* I: pixel size for the reflectance bands */
    int16 **sband,      /* I/O: input TOA and output surface reflectance */
    Geoloc_t *space,    /* I: structure for geolocation information */
    Space_def_t *space_def, /* I: structure to define the space mapping */
    float xts,          /* I: solar zenith angle (deg) */
    float xfs,          /* I: solar azimuth angle (deg) */
    float xtv,          /* I: observation zenith angle (deg) */
    float xmus,         /* I: cosine of solar zenith angle */
    float xmuv,         /* I: cosine of observation zenith angle */
    float xfi,          /* I: azimuthal difference between sun and
                              observation (deg) */
    float cosxfi,       /* I: cosine of azimuthal difference */
    float raot550nm,    /* I: nearest value of AOT */
    float pres,         /* I: surface pressure */
    float uoz,          /* I: total column ozone */
    float uwv,          /* I: total column water vapor (precipital water
                              vapor) */
    float **tsmax,      /* I: maximum scattering angle table [20][22] */
    float **tsmin,      /* I: minimum scattering angle table [20][22] */
    float xtsstep,      /* I: solar zenith step value */
    float xtsmin,       /* I: minimum solar zenith value */
    float xtvstep,      /* I: observation step value */
    float xtvmin,       /* I: minimum observation value */
    float tts[22],      /* I: sun angle table */
    float **ttv,        /* I: view angle table [20][22] */
    int32 indts[22],    /* I: index for the sun angle table */
    float ****rolutt,   /* I: intrinsic reflectance table
                              [NSR_BANDS][7][22][8000] */
    float ****transt,   /* I: transmission table [NSR_BANDS][7][22][22] */
    float ***sphalbt,   /* I: spherical albedo table [NSR_BANDS][7][22] */
    float ***normext,   /* I: aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm)
                              [NSR_BANDS][7][22] */
    float **nbfic,      /* I: communitive number of azimuth angles [20][22] */
    float **nbfi,       /* I: number of azimuth angles [20][22] */
    int16 **dem,        /* I: CMG DEM data array [DEM_NBLAT][DEM_NBLON] */
    int16 **andwi,      /* I: avg NDWI [RATIO_NBLAT][RATIO_NBLON] */
    int16 **sndwi,      /* I: standard NDWI [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob1,    /* I: mean band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob2,    /* I: mean band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob7,    /* I: mean band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob1, /* I: integer band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob2, /* I: integer band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob7, /* I: integer band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob1, /* I: slope band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob2, /* I: slope band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob7, /* I: slope band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    uint16 **wv,        /* I: water vapor values [CMG_NBLAT][CMG_NBLON] */
    uint8 **oz          /* I: ozone values [CMG_NBLAT][CMG_NBLON] */
);

#endif
