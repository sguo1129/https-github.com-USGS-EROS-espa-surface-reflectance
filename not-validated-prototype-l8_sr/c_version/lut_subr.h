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
int invaero
(
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /*** I: intrinsic reflectance table
                                           [16][7][22][8000] */
    float pres,                      /* I: surface pressure */
    float tpres[7],                  /* I: surface pressure table */
    float ****transt,                /*** I: transmission table
                                           [16][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float ***sphalbt,                /*** I: spherical albedo table
                                           [16][7][22] */
    float **tsmax,                   /* I: [20][22] */
    float **tsmin,                   /* I: [20][22] */
    float **nbfic,                   /* I: [20][22] */
    float **nbfi,                    /* I: [20][22] */
    float tts[22],
    int32 indts[22],
    float **ttv,                     /* I: [20][22] */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float tauray[16],                /* I: molecular optical thickness coeff */
    float ogtransa0[16],             /* I: other gases transmission coeff */
    float ogtransa1[16],             /* I: other gases transmission coeff */
    float ogtransb0[16],             /* I: other gases transmission coeff */
    float ogtransb1[16],             /* I: other gases transmission coeff */
    float ogtransc0[16],             /* I: other gases transmission coeff */
    float ogtransc1[16],             /* I: other gases transmission coeff */
    float wvtransa[16],              /* I: water vapor transmission coeff */
    float wvtransb[16],              /* I: water vapor transmission coeff */
    float wvtransc[16],              /* I: water vapor transmission coeff */
    float oztransa[16],              /* I: ozone transmission coeff */
    float trotoa[16],                /* I: top of atmos reflectance table */
    float erelc[16],
    int iband1,                      /* I: band 1 index (0-based) */
    int iband2,                      /* I: band 2 index (0-based) */
    float *raot550nm,                /* O: nearest value of AOT */
    float *roslamb1,                 /* O: lambertian surface reflectance of
                                           band 1 */
    float *residual                  /* O: model residual */
);

int invaeroocean
(
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /*** I: intrinsic reflectance table
                                           [16][7][22][8000] */
    float pres,                      /* I: surface pressure */
    float tpres[7],                  /* I: surface pressure table */
    float ****transt,                /*** I: transmission table
                                           [16][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float ***sphalbt,                /*** I: spherical albedo table
                                           [16][7][22] */
    float **tsmax,                   /* I: [20][22] */
    float **tsmin,                   /* I: [20][22] */
    float **nbfic,                   /* I: [20][22] */
    float **nbfi,                    /* I: [20][22] */
    float tts[22],
    int32 indts[22],
    float **ttv,                     /* I: [20][22] */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float tauray[16],                /* I: molecular optical thickness coeff */
    float ogtransa0[16],             /* I: other gases transmission coeff */
    float ogtransa1[16],             /* I: other gases transmission coeff */
    float ogtransb0[16],             /* I: other gases transmission coeff */
    float ogtransb1[16],             /* I: other gases transmission coeff */
    float ogtransc0[16],             /* I: other gases transmission coeff */
    float ogtransc1[16],             /* I: other gases transmission coeff */
    float wvtransa[16],              /* I: water vapor transmission coeff */
    float wvtransb[16],              /* I: water vapor transmission coeff */
    float wvtransc[16],              /* I: water vapor transmission coeff */
    float oztransa[16],              /* I: ozone transmission coeff */
    float trotoa[16],                /* I: top of atmos reflectance table */
    float erelc[16],
    int iband1,                      /* I: band 1 index (0-based) */
    int iband2,                      /* I: band 2 index (0-based) */
    float *aot2,
    float *roslamb1,                 /* O: lambertian surface reflectance of
                                           band 1 */
    float *residual,                 /* O: model residual */
    float *angexp                    /* O: Angstrom exponent */
);

int atmcorocea2
(
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float aot2,
    int iband,                       /* I: band index (0-based) */
    float pres,                      /* I: surface pressure */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /*** I: intrinsic reflectance table
                                           [16][7][22][8000] */
    float ****transt,                /*** I: transmission table
                                           [16][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float ***sphalbt,                /*** I: spherical albedo table
                                           [16][7][22] */
    float **tsmax,                   /* I: [20][22] */
    float **tsmin,                   /* I: [20][22] */
    float **nbfic,                   /* I: [20][22] */
    float **nbfi,                    /* I: [20][22] */
    float tts[22],
    int32 indts[22],
    float **ttv,                     /* I: [20][22] */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float tauray[16],                /* I: molecular optical thickness coeff */
    float ogtransa0[16],             /* I: other gases transmission coeff */
    float ogtransa1[16],             /* I: other gases transmission coeff */
    float ogtransb0[16],             /* I: other gases transmission coeff */
    float ogtransb1[16],             /* I: other gases transmission coeff */
    float ogtransc0[16],             /* I: other gases transmission coeff */
    float ogtransc1[16],             /* I: other gases transmission coeff */
    float wvtransa[16],              /* I: water vapor transmission coeff */
    float wvtransb[16],              /* I: water vapor transmission coeff */
    float wvtransc[16],              /* I: water vapor transmission coeff */
    float oztransa[16],              /* I: ozone transmission coeff */
    float rotoa,                     /* I: top of atmosphere reflectance */
    float *roslamb,                  /* O: lambertian surface reflectance */
    float angexp,
    float *tgo,                      /* O: other gaseous transmittance */
    float *roatm,                    /* O: atmospheric reflectance */
    float *ttatmg,
    float *satm,                     /* O: spherical albedo */
    float *xrorayp                   /* O: molecular reflectance */
);

int atmcorlamb2
(
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float raot550nm,                 /* I: nearest value of AOT */
    int iband,                       /* I: band index (0-based) */
    float pres,                      /* I: surface pressure */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /*** I: intrinsic reflectance table
                                           [16][7][22][8000] */
    float ****transt,                /*** I: transmission table
                                           [16][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float ***sphalbt,                /*** I: spherical albedo table
                                           [16][7][22] */
    float **tsmax,                   /* I: [20][22] */
    float **tsmin,                   /* I: [20][22] */
    float **nbfic,                   /* I: [20][22] */
    float **nbfi,                    /* I: [20][22] */
    float tts[22],
    int32 indts[22],
    float **ttv,                     /* I: [20][22] */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float tauray[16],                /* I: molecular optical thickness coeff */
    float ogtransa0[16],             /* I: other gases transmission coeff */
    float ogtransa1[16],             /* I: other gases transmission coeff */
    float ogtransb0[16],             /* I: other gases transmission coeff */
    float ogtransb1[16],             /* I: other gases transmission coeff */
    float ogtransc0[16],             /* I: other gases transmission coeff */
    float ogtransc1[16],             /* I: other gases transmission coeff */
    float wvtransa[16],              /* I: water vapor transmission coeff */
    float wvtransb[16],              /* I: water vapor transmission coeff */
    float wvtransc[16],              /* I: water vapor transmission coeff */
    float oztransa[16],              /* I: ozone transmission coeff */
    float rotoa,                     /* I: top of atmosphere reflectance */
    float *roslamb,                  /* O: lambertian surface reflectance */
    float *tgo,                      /* O: other gaseous transmittance */
    float *roatm,                    /* O: atmospheric reflectance */
    float *ttatmg,
    float *satm,                     /* O: spherical albedo */
    float *xrorayp                   /* O: molecular reflectance */
);

void raycorlamb2
(
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    int iband,                       /* I: band index (0-based) */
    float pres,                      /* I: surface pressure */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float tauray[16],                /* I: molecular optical thickness coeff */
    float ogtransa0[16],             /* I: other gases transmission coeff */
    float ogtransa1[16],             /* I: other gases transmission coeff */
    float ogtransb0[16],             /* I: other gases transmission coeff */
    float ogtransb1[16],             /* I: other gases transmission coeff */
    float ogtransc0[16],             /* I: other gases transmission coeff */
    float ogtransc1[16],             /* I: other gases transmission coeff */
    float wvtransa[16],              /* I: water vapor transmission coeff */
    float wvtransb[16],              /* I: water vapor transmission coeff */
    float wvtransc[16],              /* I: water vapor transmission coeff */
    float oztransa[16],              /* I: ozone transmission coeff */
    float rotoa,                     /* I: top of atmosphere reflectance */
    float *roslamb,                  /* O: lambertian surface reflectance */
    float *tgo,                      /* O: other gaseous transmittance */
    float *roatm,                    /* O: atmospheric reflectance */
    float *ttatmg,
    float *satm,                     /* O: spherical albedo */
    float *xrorayp                   /* O: molecular reflectance */
);

int atmcorlamb
(
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float raot550nm,
    int iband,                       /* I: band index (0-based) */
    float pres,                      /* I: surface pressure */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /*** I: intrinsic reflectance table
                                           [16][7][22][8000] */
    float ****transt,                /*** I: transmission table
                                           [16][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float ***sphalbt,                /*** I: spherical albedo table
                                           [16][7][22] */
    float **tsmax,                   /* I: [20][22] */
    float **tsmin,                   /* I: [20][22] */
    float **nbfic,                   /* I: [20][22] */
    float **nbfi,                    /* I: [20][22] */
    float tts[22],
    int32 indts[22],
    float **ttv,                     /* I: [20][22] */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float tauray[16],                /* I: molecular optical thickness coeff */
    float ogtransa0[16],             /* I: other gases transmission coeff */
    float ogtransa1[16],             /* I: other gases transmission coeff */
    float ogtransb0[16],             /* I: other gases transmission coeff */
    float ogtransb1[16],             /* I: other gases transmission coeff */
    float ogtransc0[16],             /* I: other gases transmission coeff */
    float ogtransc1[16],             /* I: other gases transmission coeff */
    float wvtransa[16],              /* I: water vapor transmission coeff */
    float wvtransb[16],              /* I: water vapor transmission coeff */
    float wvtransc[16],              /* I: water vapor transmission coeff */
    float oztransa[16],              /* I: ozone transmission coeff */
    float rotoa,                     /* I: top of atmosphere reflectance */
    float *roslamb                   /* O: lambertian surface reflectance */
);

void local_csalbr
(
    float xtau,       /* I: molecular optical depth */
    float *xalb       /* O: atmospheric (Rayleigh) spherical albedo */
);

float fintexp3
(
    float xtau       /* I: molecular optical depth */
);

float fintexp1
(
    float xtau       /* I: molecular optical depth */
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

void comptransray
(
    float xtaur,   /* I: rayleigh optical depth for surface pressure */
    float xmus,    /* I: cosine of solar zenith angle */
    float *ttray   /* O: */
);

void comptg
(
    int iband,                   /* I: band index (0-based) */
    float xts,                   /* I: solar zenith angle */
    float xtv,                   /* I: observation zenith angle */
    float uoz,                   /* I: total column ozone */
    float uwv,                   /* I: total column water vapor (precipital
                                       water vapor) */
    float pres,                  /* I: surface pressure */
    float ogtransa0[16],         /* I: other gases transmission coeff */
    float ogtransa1[16],         /* I: other gases transmission coeff */
    float ogtransb0[16],         /* I: other gases transmission coeff */
    float ogtransb1[16],         /* I: other gases transmission coeff */
    float ogtransc0[16],         /* I: other gases transmission coeff */
    float ogtransc1[16],         /* I: other gases transmission coeff */
    float wvtransa[16],          /* I: water vapor transmission coeff */
    float wvtransb[16],          /* I: water vapor transmission coeff */
    float wvtransc[16],          /* I: water vapor transmission coeff */
    float oztransa[16],          /* I: ozone transmission coeff */
    float *tgoz,                 /* O: ozone transmission */
    float *tgwv,                 /* O: water vapor transmission */
    float *tgwvhalf,             /* O: water vapor transmission, half content */
    float *tgog                  /* O: other gases transmission */
);

void compsalb
(
    float raot550nm,                 /* I: nearest value of AOT */
    int iband,                       /* I: band index (0-based) */
    float pres,                      /* I: surface pressure */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ***sphalbt,                /*** I: spherical albedo table
                                           [16][7][22] */
    float *satm                      /* O: spherical albedo */
);

int comptrans
(
    float xts,                       /* I: solar zenith */
    float raot550nm,                 /* I: nearest value of AOT */
    int iband,                       /* I: band index (0-based) */
    float pres,                      /* I: surface pressure */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****transt,                /*** I: transmission table
                                           [16][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float tts[22],
    float *xtts                      /* O: downward transmittance */
);

int comproatm
(
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float raot550nm,                 /* I: nearest value of AOT */
    int iband,                       /* I: band index (0-based) */
    float pres,                      /* I: surface pressure */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /*** I: intrinsic reflectance table
                                           [16][7][22][8000] */
    float **tsmax,                   /* I: [20][22] */
    float **tsmin,                   /* I: [20][22] */
    float **nbfic,                   /* I: [20][22] */
    float **nbfi,                    /* I: [20][22] */
    float tts[22],
    int32 indts[22],
    float **ttv,                     /* I: [20][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float *roatm                     /* O: atmospheric reflectance */
);

int readluts
(
    float tauray[16],                /* O: molecular optical thickness coeff */
    float oztransa[16],              /* O: ozone transmission coeff */
    float wvtransa[16],              /* O: water vapor transmission coeff */
    float wvtransb[16],              /* O: water vapor transmission coeff */
    float wvtransc[16],              /* O: water vapor transmission coeff */
    float ogtransa0[16],             /* O: other gases transmission coeff */
    float ogtransa1[16],             /* O: other gases transmission coeff */
    float ogtransb0[16],             /* O: other gases transmission coeff */
    float ogtransb1[16],             /* O: other gases transmission coeff */
    float ogtransc0[16],             /* O: other gases transmission coeff */
    float ogtransc1[16],             /* O: other gases transmission coeff */
    float **tsmax,                   /* O: [20][22] */
    float **tsmin,                   /* O: [20][22] */
    float **ttv,                     /* O: [20][22] */
    float tts[22],                   /* O: */
    float **nbfic,                   /* O: [20][22] */
    float **nbfi,                    /* O: [20][22] */
    int32 indts[22],                 /* O: */
    float ****rolutt,                /*** O: intrinsic reflectance table
                                           [16][7][22][8000] */
    float ****transt,                /*** O: transmission table
                                           [16][7][22][22] */
    float ***sphalbt,                /*** O: spherical albedo table
                                           [16][7][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    char sbandname[16][STR_SIZE],    /* I: "band" names for the molecular
                                           optical thickness file */
    char tauraynm[STR_SIZE],     /* I: molecular optical thickness filename */
    char gscoefnm[STR_SIZE],     /* I: gaseous transmission coef filename */
    char anglehdf[STR_SIZE],     /* I: angle HDF filename */
    char intrefnm[STR_SIZE],     /* I: intrinsic reflectance filename */
    char transmnm[STR_SIZE],     /* I: transmission filename */
    char spheranm[STR_SIZE]      /* I: spherical albedo filename */
);

int subaeroret
(
    int iband1,                      /* I: band 1 index (0-based) */
    int iband3,                      /* I: band 3 index (0-based) */
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float pres,                      /* I: surface pressure */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float erelc[16],
    float troatm[16],
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /*** I: intrinsic reflectance table
                                           [16][7][22][8000] */
    float ****transt,                /*** I: transmission table
                                           [16][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float ***sphalbt,                /*** I: spherical albedo table
                                           [16][7][22] */
    float **tsmax,                   /* I: [20][22] */
    float **tsmin,                   /* I: [20][22] */
    float **nbfic,                   /* I: [20][22] */
    float **nbfi,                    /* I: [20][22] */
    float tts[22],
    int32 indts[22],
    float **ttv,                     /* I: [20][22] */
    float tauray[16],                /* I: molecular optical thickness coeff */
    float ogtransa0[16],             /* I: other gases transmission coeff */
    float ogtransa1[16],             /* I: other gases transmission coeff */
    float ogtransb0[16],             /* I: other gases transmission coeff */
    float ogtransb1[16],             /* I: other gases transmission coeff */
    float ogtransc0[16],             /* I: other gases transmission coeff */
    float ogtransc1[16],             /* I: other gases transmission coeff */
    float wvtransa[16],              /* I: water vapor transmission coeff */
    float wvtransb[16],              /* I: water vapor transmission coeff */
    float wvtransc[16],              /* I: water vapor transmission coeff */
    float oztransa[16],              /* I: ozone transmission coeff */
    float *raot,
    float *residual,                 /* O: model residual */
    int *nit                         /* O: number of iterations */
);

#endif
