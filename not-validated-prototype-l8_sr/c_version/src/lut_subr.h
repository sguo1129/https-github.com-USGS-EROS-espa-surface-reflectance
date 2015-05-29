#ifndef _LUT_SUBR_H_
#define _LUT_SUBR_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "common.h"
#include "espa_metadata.h"
#include "error_handler.h"

/* Prototypes */
int atmcorlamb2
(
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xmus,                      /* I: cosine of solar zenith angle */
    float xmuv,                      /* I: cosine of observation zenith angle */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float cosxfi,                    /* I: cosine of azimuthal difference */
    float raot550nm,                 /* I: nearest value of AOT */
    int iband,                       /* I: band index (0-based) */
    float pres,                      /* I: surface pressure */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /* I: intrinsic reflectance table
                                           [NSR_BANDS][7][22][8000] */
    float ****transt,                /* I: transmission table
                                           [NSR_BANDS][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float ***sphalbt,                /* I: spherical albedo table
                                           [NSR_BANDS][7][22] */
    float ***normext,                /* I: ?????
                                           [NSR_BANDS][7][22] */
    float **tsmax,                   /* I: maximum scattering angle table
                                           [20][22] */
    float **tsmin,                   /* I: minimum scattering angle table
                                           [20][22] */
    float **nbfic,                   /* I: communitive number of azimuth angles
                                           [20][22] */
    float **nbfi,                    /* I: number of azimuth angles [20][22] */
    float tts[22],                   /* I: sun angle table */
    int32 indts[22],
    float **ttv,                     /* I: view angle table [20][22] */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float tauray[NSR_BANDS],         /* I: molecular optical thickness coeff */
    double ogtransa1[NSR_BANDS],     /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS],     /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS],     /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],      /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],      /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],      /* I: ozone transmission coeff */
    float rotoa,                     /* I: top of atmosphere reflectance */
    float *roslamb,                  /* O: lambertian surface reflectance */
    float *tgo,                      /* O: other gaseous transmittance */
    float *roatm,                    /* O: atmospheric reflectance */
    float *ttatmg,
    float *satm,                     /* O: spherical albedo */
    float *xrorayp,                  /* O: molecular reflectance */
    float *next                      /* O: ???? */
);

void local_chand
(
    float xphi,    /* I: azimuthal difference between sun and observation
                         (deg) */
    float xmuv,    /* I: cosine of observation zenith angle */
    float xmus,    /* I: cosine of solar zenith angle */
    float xtau,    /* I: molecular optical depth */
    float *xrray   /* O: molecular reflectance, 0.0 to 1.0 */
);

void comptg
(
    int iband,                   /* I: band index (0-based) */
    float xts,                   /* I: solar zenith angle */
    float xtv,                   /* I: observation zenith angle */
    float xmus,                  /* I: cosine of solar zenith angle */
    float xmuv,                  /* I: cosine of observation zenith angle */
    float uoz,                   /* I: total column ozone */
    float uwv,                   /* I: total column water vapor (precipital
                                       water vapor) */
    float atm_pres,              /* I: pressure at sea level */
    double ogtransa1[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS], /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS], /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],  /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],  /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],  /* I: ozone transmission coeff */
    float *tgoz,                 /* O: ozone transmission */
    float *tgwv,                 /* O: water vapor transmission */
    float *tgwvhalf,             /* O: water vapor transmission, half content */
    float *tgog                  /* O: other gases transmission */
);

void compsalb
(
    int ip1,            /* I: index variable for surface pressure */
    int ip2,            /* I: index variable for surface pressure */
    int iaot1,          /* I: index variable for AOT */
    int iaot2,          /* I: index variable for AOT */
    float raot550nm,    /* I: nearest value of AOT */
    int iband,          /* I: band index (0-based) */
    float pres,         /* I: surface pressure */
    float tpres[7],     /* I: surface pressure table */
    float aot550nm[22], /* I: AOT look-up table */
    float ***sphalbt,   /* I: spherical albedo table [NSR_BANDS][7][22] */
    float ***normext,   /* I: aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm)
                              [NSR_BANDS][7][22] */
    float *satm,        /* O: spherical albedo */
    float *next         /* O: ????? */
);

void comptrans
(
    int ip1,            /* I: index variable for surface pressure */
    int ip2,            /* I: index variable for surface pressure */
    int iaot1,          /* I: index variable for AOT */
    int iaot2,          /* I: index variable for AOT */
    float xts,          /* I: zenith angle */
    float raot550nm,    /* I: nearest value of AOT */
    int iband,          /* I: band index (0-based) */
    float pres,         /* I: surface pressure */
    float tpres[7],     /* I: surface pressure table */
    float aot550nm[22], /* I: AOT look-up table */
    float ****transt,   /* I: transmission table
                              [NSR_BANDS][7][22][22] */
    float xtsstep,      /* I: zenith angle step value */
    float xtsmin,       /* I: minimum zenith angle value */
    float tts[22],      /* I: sun angle table */
    float *xtts         /* O: downward transmittance */
);

void comproatm
(
    int ip1,            /* I: index variable for surface pressure */
    int ip2,            /* I: index variable for surface pressure */
    int iaot1,          /* I: index variable for AOT */
    int iaot2,          /* I: index variable for AOT */
    float xts,          /* I: solar zenith angle (deg) */
    float xtv,          /* I: observation zenith angle (deg) */
    float xmus,         /* I: cosine of solar zenith angle */
    float xmuv,         /* I: cosine of observation zenith angle */
    float cosxfi,       /* I: cosine of azimuthal difference */
    float raot550nm,    /* I: nearest value of AOT */
    int iband,          /* I: band index (0-based) */
    float pres,         /* I: surface pressure */
    float tpres[7],     /* I: surface pressure table */
    float aot550nm[22], /* I: AOT look-up table */
    float ****rolutt,   /* I: intrinsic reflectance table
                              [NSR_BANDS][7][22][8000] */
    float **tsmax,      /* I: maximum scattering angle table [20][22] */
    float **tsmin,      /* I: minimum scattering angle table [20][22] */
    float **nbfic,      /* I: communitive number of azimuth angles [20][22] */
    float **nbfi,       /* I: number of azimuth angles [20][22] */
    float tts[22],      /* I: sun angle table */
    int32 indts[22],
    float **ttv,        /* I: view angle table [20][22] */
    float xtsstep,      /* I: solar zenith step value */
    float xtsmin,       /* I: minimum solar zenith value */
    float xtvstep,      /* I: observation step value */
    float xtvmin,       /* I: minimum observation value */
    int its,            /* I: index for the sun angle table */
    int itv,            /* I: index for the view angle table */
    float *roatm        /* O: atmospheric reflectance */
);

int readluts
(
    float **tsmax,              /* O: maximum scattering angle table [20][22] */
    float **tsmin,              /* O: minimum scattering angle table [20][22] */
    float **ttv,                /* O: view angle table [20][22] */
    float tts[22],              /* O: sun angle table */
    float **nbfic,              /* O: communitive number of azimuth angles
                                      [20][22] */
    float **nbfi,               /* O: number of azimuth angles [20][22] */
    int32 indts[22],            /* O: */
    float ****rolutt,           /* O: intrinsic reflectance table
                                      [NSR_BANDS][7][22][8000] */
    float ****transt,           /* O: transmission table
                                      [NSR_BANDS][7][22][22] */
    float ***sphalbt,           /* O: spherical albedo table
                                      [NSR_BANDS][7][22] */
    float ***normext,           /* O: ?????
                                      [NSR_BANDS][7][22] */
    float xtsstep,              /* I: solar zenith step value */
    float xtsmin,               /* I: minimum solar zenith value */
    char anglehdf[STR_SIZE],    /* I: angle HDF filename */
    char intrefnm[STR_SIZE],    /* I: intrinsic reflectance filename */
    char transmnm[STR_SIZE],    /* I: transmission filename */
    char spheranm[STR_SIZE]     /* I: spherical albedo filename */
);

int subaeroret
(
    int iband1,                      /* I: band 1 index (0-based) */
    int iband3,                      /* I: band 3 index (0-based) */
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xmus,                      /* I: cosine of solar zenith angle */
    float xmuv,                      /* I: cosine of observation zenith angle */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float cosxfi,                    /* I: cosine of azimuthal difference */
    float pres,                      /* I: surface pressure */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float erelc[NSR_BANDS],          /* I: band ratio variable */
    float troatm[NSR_BANDS],         /* I: atmospheric reflectance table */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /* I: intrinsic reflectance table
                                           [NSR_BANDS][7][22][8000] */
    float ****transt,                /* I: transmission table
                                           [NSR_BANDS][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float ***sphalbt,                /* I: spherical albedo table
                                           [NSR_BANDS][7][22] */
    float ***normext,                /* I: ????
                                           [NSR_BANDS][7][22] */
    float **tsmax,                   /* I: maximum scattering angle table
                                           [20][22] */
    float **tsmin,                   /* I: minimum scattering angle table
                                           [20][22] */
    float **nbfic,                   /* I: communitive number of azimuth angles
                                           [20][22] */
    float **nbfi,                    /* I: number of azimuth anglesi [20][22] */
    float tts[22],                   /* I: sun angle table */
    int32 indts[22],
    float **ttv,                     /* I: view angle table [20][22] */
    float tauray[NSR_BANDS],         /* I: molecular optical thickness coeff */
    double ogtransa1[NSR_BANDS],     /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS],     /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS],     /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],      /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],      /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],      /* I: ozone transmission coeff */
    float *raot,                     /* O: AOT reflectance */
    float *residual,                 /* O: model residual */
    float *snext                     /* O: ????? */
);

int subaeroret_residual
(
    int iband1,                      /* I: band 1 index (0-based) */
    int iband3,                      /* I: band 3 index (0-based) */
    double ros1,                     /* I: surface reflectance for band 1 */
    double ros3,                     /* I: surface reflectance for band 3 */
    float roslamb,                   /* I: lambertian surface reflectance */
    double pratio,                   /* I: targeted ratio between the surface
                                           reflectance in two bands */
    float raot550nm,                 /* I: nearest input value of AOT */

    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xmus,                      /* I: cosine of solar zenith angle */
    float xmuv,                      /* I: cosine of observation zenith angle */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float cosxfi,                    /* I: cosine of azimuthal difference */
    float pres,                      /* I: surface pressure */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float erelc[NSR_BANDS],          /* I: band ratio variable */
    float troatm[NSR_BANDS],         /* I: atmospheric reflectance table */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /* I: intrinsic reflectance table
                                           [NSR_BANDS][7][22][8000] */
    float ****transt,                /* I: transmission table
                                           [NSR_BANDS][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float ***sphalbt,                /* I: spherical albedo table
                                           [NSR_BANDS][7][22] */
    float ***normext,                /* I: ????
                                           [NSR_BANDS][7][22] */
    float **tsmax,                   /* I: maximum scattering angle table
                                           [20][22] */
    float **tsmin,                   /* I: minimum scattering angle table
                                           [20][22] */
    float **nbfic,                   /* I: communitive number of azimuth angles
                                           [20][22] */
    float **nbfi,                    /* I: number of azimuth anglesi [20][22] */
    float tts[22],                   /* I: sun angle table */
    int32 indts[22],
    float **ttv,                     /* I: view angle table [20][22] */
    float tauray[NSR_BANDS],         /* I: molecular optical thickness coeff */
    double ogtransa1[NSR_BANDS],     /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS],     /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS],     /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],      /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],      /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],      /* I: ozone transmission coeff */
    float *residual,                 /* O: model residual */
    float *snext                     /* O: ????? */
);

int memory_allocation_main
(
    int nlines,          /* I: number of lines in the scene */
    int nsamps,          /* I: number of samples in the scene */
    uint16 **qaband,     /* O: QA band for the input image, nlines x nsamps */
    int16 ***sband       /* O: output surface reflectance and brightness temp
                               bands */
);

int memory_allocation_sr
(
    int nlines,          /* I: number of lines in the scene */
    int nsamps,          /* I: number of samples in the scene */
    int16 **aerob1,      /* O: atmospherically corrected band 1 data
                               (TOA refl), nlines x nsamps */
    int16 **aerob2,      /* O: atmospherically corrected band 2 data
                               (TOA refl), nlines x nsamps */
    int16 **aerob4,      /* O: atmospherically corrected band 4 data
                               (TOA refl), nlines x nsamps */
    int16 **aerob5,      /* O: atmospherically corrected band 5 data
                               (TOA refl), nlines x nsamps */
    int16 **aerob7,      /* O: atmospherically corrected band 7 data
                               (TOA refl), nlines x nsamps */
    uint8 **cloud,       /* O: bit-packed value that represent clouds,
                               nlines x nsamps */
    float **twvi,        /* O: interpolated water vapor value,
                               nlines x nsamps */
    float **tozi,        /* O: interpolated ozone value, nlines x nsamps */
    float **tp,          /* O: interpolated pressure value, nlines x nsamps */
    float **tresi,       /* O: residuals for each pixel, nlines x nsamps */
    float **taero,       /* O: aerosol values for each pixel, nlines x nsamps */
    uint8 **lw_mask,     /* O: land/water mask data, nlines x nsamps */
    int16 ***dem,        /* O: CMG DEM data array [DEM_NBLAT][DEM_NBLON] */
    int16 ***andwi,      /* O: avg NDWI [RATIO_NBLAT][RATIO_NBLON] */
    int16 ***sndwi,      /* O: standard NDWI [RATIO_NBLAT][RATIO_NBLON] */
    int16 ***ratiob1,    /* O: mean band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 ***ratiob2,    /* O: mean band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 ***ratiob7,    /* O: mean band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 ***intratiob1, /* O: band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 ***intratiob2, /* O: band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 ***intratiob7, /* O: band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 ***slpratiob1, /* O: slope band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 ***slpratiob2, /* O: slope band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 ***slpratiob7, /* O: slope band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    uint16 ***wv,        /* O: water vapor values [CMG_NBLAT][CMG_NBLON] */
    uint8 ***oz,         /* O: ozone values [CMG_NBLAT][CMG_NBLON] */
    float *****rolutt,   /* O: intrinsic reflectance table
                               [NSR_BANDS][7][22][8000] */
    float *****transt,   /* O: transmission table
                               [NSR_BANDS][7][22][22] */
    float ****sphalbt,   /* O: spherical albedo table [NSR_BANDS][7][22] */
    float ****normext,   /* O: aerosol extinction coefficient at the current
                               wavelength (normalized at 550nm)
                               [NSR_BANDS][7][22] */
    float ***tsmax,      /* O: maximum scattering angle table [20][22] */
    float ***tsmin,      /* O: minimum scattering angle table [20][22] */
    float ***nbfic,      /* O: communitive number of azimuth angles [20][22] */
    float ***nbfi,       /* O: number of azimuth angles [20][22] */
    float ***ttv         /* O: view angle table [20][22] */
);

int read_auxiliary_files
(
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm,        /* I: auxiliary filename for ozone and water vapor */
    int16 **dem,        /* O: CMG DEM data array [DEM_NBLAT][DEM_NBLON] */
    int16 **andwi,      /* O: avg NDWI [RATIO_NBLAT][RATIO_NBLON] */
    int16 **sndwi,      /* O: standard NDWI [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob1,    /* O: mean band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob2,    /* O: mean band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **ratiob7,    /* O: mean band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob1, /* O: band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob2, /* O: band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **intratiob7, /* O: band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob1, /* O: slope band1 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob2, /* O: slope band2 ratio [RATIO_NBLAT][RATIO_NBLON] */
    int16 **slpratiob7, /* O: slope band7 ratio [RATIO_NBLAT][RATIO_NBLON] */
    uint16 **wv,        /* O: water vapor values [CMG_NBLAT][CMG_NBLON] */
    uint8 **oz          /* O: ozone values [CMG_NBLAT][CMG_NBLON] */
);

#endif
