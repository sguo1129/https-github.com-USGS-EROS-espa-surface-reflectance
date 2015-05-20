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
#include "l8_angles.h"

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
    short *xts_arr,     /* I: scaled solar zenith angle (deg) from a 
                              representative band or band average */
    char *instrument,   /* I: instrument to be processed (OLI, TIRS) */
    int16 **sband       /* O: output TOA reflectance and brightness temp
                              values (scaled) */
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
    short *xts_arr,     /* I: scaled solar zenith angle (deg) from a 
                              representative band or band average */
    short *xfs_arr,     /* I: scaled solar azimuth angle (deg) from a
                              representative band or band average */
    short *xtv_arr,     /* I: scaled observation zenith angle (deg) from a
                              representative band or band average */
    short *xfv_arr,     /* I: scaled observation azimuth angle (deg) from a
                              representative band or band average */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm         /* I: auxiliary filename for ozone and water vapor */
);

int init_sr_refl
(
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    Input_t *input,     /* I: input structure for the Landsat product */
    Geoloc_t *space,    /* I: structure for geolocation information */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm,        /* I: auxiliary filename for ozone and water vapor */
    float *raot550nm,   /* O: nearest value of AOT */
    float *pres,        /* O: surface pressure */
    float *uoz,         /* O: total column ozone */
    float *uwv,         /* O: total column water vapor (precipital water
                              vapor) */
    float *xtsstep,     /* O: solar zenith step value */
    float *xtsmin,      /* O: minimum solar zenith value */
    float *xtvstep,     /* O: observation step value */
    float *xtvmin,      /* O: minimum observation value */
    float **tsmax,      /* O: maximum scattering angle table [20][22] */
    float **tsmin,      /* O: minimum scattering angle table [20][22] */
    float tts[22],      /* O: sun angle table */
    float **ttv,        /* O: view angle table [20][22] */
    int32 indts[22],    /* O: index for the sun angle table */
    float ****rolutt,   /* O: intrinsic reflectance table
                              [NSR_BANDS][7][22][8000] */
    float ****transt,   /* O: transmission table [NSR_BANDS][7][22][22] */
    float ***sphalbt,   /* O: spherical albedo table [NSR_BANDS][7][22] */
    float ***normext,   /* O: aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm)
                              [NSR_BANDS][7][22] */
    float **nbfic,      /* O: communitive number of azimuth angles [20][22] */
    float **nbfi,       /* O: number of azimuth angles [20][22] */
    int16 **dem,        /* O: CMG DEM data array [DEM_NBLAT][DEM_NBLON] */
    int16 **andwi,      /* O: avg NDWI [RATIO_NBLAT][RATIO_NBLON] */
    int16 **sndwi,      /* O: standard NDWI [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob1,    /* O: mean band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob2,    /* O: mean band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob7,    /* O: mean band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob1, /* O: integer band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob2, /* O: integer band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob7, /* O: integer band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob1, /* O: slope band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob2, /* O: slope band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob7, /* O: slope band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    uint16 **wv,        /* O: water vapor values [CMG_NBLAT][CMG_NBLON] */
    uint8 **oz          /* O: ozone values [CMG_NBLAT][CMG_NBLON] */
);

#endif
