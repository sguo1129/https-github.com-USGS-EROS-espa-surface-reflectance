/*****************************************************************************
FILE: lut_subr.c
  
PURPOSE: Contains functions for reading the look-up tables and doing some
of the coefficient computations for the surface reflectance application.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/
#include "lut_subr.h"
#include "hdf.h"
#include "mfhdf.h"

/******************************************************************************
MODULE:  atmcorlamb2

PURPOSE:  Lambertian atmospheric correction 2.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred doing the atmospheric corrections.
SUCCESS        Successful completion

NOTES:
    1. Standard sea level pressure is 1013 millibars.
******************************************************************************/
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
    float ***normext,                /* I: aerosol extinction coefficient at
                                           the current wavelength (normalized
                                           at 550nm) [NSR_BANDS][7][22] */
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
    float *ttatmg,                   /* O: total atmospheric transmission */
    float *satm,                     /* O: spherical albedo */
    float *xrorayp,                  /* O: molecular reflectance */
    float *next                      /* O: */
)
{
    char FUNC_NAME[] = "atmcorlamb2";   /* function name */
    char errmsg[STR_SIZE];  /* error message */
    float xttv;         /* upward transmittance */
    float xtts;         /* downward transmittance */
    float ttatm;        /* total transmission of the atmosphere */
    float tgog;         /* other gases transmission */
    float tgoz;         /* ozone transmission */
    float tgwv;         /* water vapor transmission */
    float tgwvhalf;     /* water vapor transmission, half content */
    float xtaur;        /* rayleigh optical depth for surface pressure */
    float atm_pres;     /* atmospheric pressure at sea level */
    int ip;             /* surface pressure looping variable */
    int ip1, ip2;       /* index variables for the surface pressure */
    int iaot;           /* aerosol optical thickness (AOT) looping variable */
    int iaot1, iaot2;   /* index variables for the AOT and spherical albedo
                           arrays */
    int its;            /* index for the sun angle table */
    int itv;            /* index for the view angle table */

    /* Get the pressure and AOT related values for the current surface pressure
       and AOT.  These indices are passed into several functions. */
    /* Look for the appropriate pressure index in the surface pressure table.
       Stop at the second to last item in the table, so that we have the last
       two elements to use as ip1 and ip2, if needed. */
    ip1 = 0;
    for (ip = 0; ip < 6; ip++)  /* 7 elements in the array, stop one short */
    {
        if (pres < tpres[ip])
            ip1 = ip;
    }
    ip2 = ip1 + 1;
      
    /* Look for the appropriate AOT index in the AOT table.
       Stop at the second to last item in the table, so that we have the last
       two elements to use as iaot1 and iaot2, if needed. */
    iaot1 = 0;
    for (iaot = 0; iaot < 21; iaot++) /* 22 elements in table, stop one short */
    {
        if (raot550nm > aot550nm[iaot])
            iaot1 = iaot;
    }
    iaot2 = iaot1 + 1;

    /* Determine the index in the view angle table */
    if (xtv <= xtvmin)
        itv = 0;
    else
        itv = (int) ((xtv - xtvmin) / xtvstep + 1.0);

    /* Determine the index in the sun angle table */
    if (xts <= xtsmin) 
        its = 0;
    else
        its = (int) ((xts - xtsmin) / xtsstep);
    if (its > 19)
    {
        sprintf (errmsg, "Solar zenith (xts) is too large: %f", xts);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* This routine returns variables for calculating roslamb */
    comproatm (ip1, ip2, iaot1, iaot2, xts, xtv, xmus, xmuv, cosxfi,
        raot550nm, iband, pres, tpres, aot550nm, rolutt, tsmax, tsmin, nbfic,
        nbfi, tts, indts, ttv, xtsstep, xtsmin, xtvstep, xtvmin, its, itv,
        roatm);

    /* Compute the transmission for the solar zenith angle */
    comptrans (ip1, ip2, iaot1, iaot2, xts, raot550nm, iband, pres, tpres,
        aot550nm, transt, xtsstep, xtsmin, tts, &xtts);

    /* Compute the transmission for the observation zenith angle */
    comptrans (ip1, ip2, iaot1, iaot2, xtv, raot550nm, iband, pres, tpres,
        aot550nm, transt, xtvstep, xtvmin, tts, &xttv);

    /* Compute total transmission (product downward by  upward) */
    ttatm = xtts * xttv;
    
    /* Compute spherical albedo */
    compsalb (ip1, ip2, iaot1, iaot2, raot550nm, iband, pres, tpres, aot550nm,
        sphalbt, normext, satm, next);

    atm_pres = pres * ONE_DIV_1013;
    comptg (iband, xts, xtv, xmus, xmuv, uoz, uwv, atm_pres, ogtransa1,
        ogtransb0, ogtransb1, wvtransa, wvtransb, oztransa, &tgoz, &tgwv,
        &tgwvhalf, &tgog);

    /* Compute rayleigh component (intrinsic reflectance, at p=pres).
       Pressure in the atmosphere is pres / 1013. */
    xtaur = tauray[iband] * atm_pres;
    local_chand (xfi, xmuv, xmus, xtaur, xrorayp);

    /* Perform atmospheric correction */
    *roslamb = rotoa / (tgog * tgoz);
    *roslamb = (*roslamb) - ((*roatm) - (*xrorayp)) * tgwvhalf - (*xrorayp);
    *roslamb /= ttatm * tgwv;
    *roslamb = (*roslamb) / (1.0 + (*satm) * (*roslamb));
    *tgo = tgog * tgoz;
    *roatm = ((*roatm) - (*xrorayp)) * tgwvhalf + (*xrorayp);
    *ttatmg = ttatm * tgwv;

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  local_chand

PURPOSE:  Computes the atm/molecular reflectance from 0.0 to 1.0, based on the
sun and observation angles.

RETURN VALUE:
Type = None

NOTES:
 1. Here's how the xfd value was originally calculated. Given that these
    are static values, the xfd itself can really be static.
    xdep = 0.0279  // depolarization factor
    xfd = xdep / (2.0 - xdep)
        = 0.014147355
    xfd = (1.0 - xfd) / (1.0 + 2.0 * xfd)
        = .985852645 / 1.02829471
        = .958725777
******************************************************************************/
void local_chand
(
    float xphi,    /* I: azimuthal difference between sun and observation
                         (deg) */
    float xmuv,    /* I: cosine of observation zenith angle */
    float xmus,    /* I: cosine of solar zenith angle */
    float xtau,    /* I: molecular optical depth */
    float *xrray   /* O: molecular reflectance, 0.0 to 1.0 */
)
{
    int i;                             /* looping variable */
    float pl[10];
    float fs0, fs1, fs2;
    float phios;
    float xcosf2, xcosf3;
    float xph1, xph2, xph3;
    float xitm;
    float xp1, xp2, xp3;
    float cfonc1, cfonc2, cfonc3;
    float xlntau;                      /* log molecular optical depth */
    float xitot1, xitot2, xitot3;
    float xmus2, xmuv2;                /* square of xmus and xmuv */

    /* constant vars */
    const float xfd = 0.958725777;
    const float as0[10] = {
         0.33243832, -6.777104e-02, 0.16285370, 1.577425e-03,
        -0.30924818, -1.240906e-02, -0.10324388, 3.241678e-02, 0.11493334,
        -3.503695e-02};
    const float as1[2] = {0.19666292, -5.439061e-02};
    const float as2[2] = {0.14545937, -2.910845e-02};

    phios = (180.0 - xphi) * DEG2RAD;
    xcosf2 = cos (phios);
    xcosf3 = cos (2.0 * phios);

    /* xmus and xmuv squared is used frequently */
    xmus2 = xmus * xmus;
    xmuv2 = xmuv * xmuv;

    xph1 = 1.0 + (3.0 * xmus2 - 1.0) * (3.0 * xmuv2 - 1.0) * xfd * 0.125;
    xph2 = -xmus * xmuv * sqrt(1.0 - xmus2) * sqrt(1.0 - xmuv2);
    xph2 = xph2 * xfd * 0.5 * 1.5;
    xph3 = (1.0 - xmus2) * (1.0 - xmuv2);
    xph3 = xph3 * xfd * 0.5 * 0.375;

    xitm = (1.0 - exp(-xtau * (1.0 / xmus + 1.0 / xmuv))) *
        xmus / (4.0 * (xmus + xmuv));
    xp1 = xph1 * xitm;
    xp2 = xph2 * xitm;
    xp3 = xph3 * xitm;

    xitm = (1.0 - exp(-xtau / xmus)) * (1.0 - exp(-xtau / xmuv));
    cfonc1 = xph1 * xitm;
    cfonc2 = xph2 * xitm;
    cfonc3 = xph3 * xitm;

    xlntau = log (xtau);
    pl[0] = 1.0;
    pl[1] = xlntau;
    pl[2] = xmus + xmuv;
    pl[3] = xlntau * pl[2];
    pl[4] = xmus * xmuv;
    pl[5] = xlntau * pl[4];
    pl[6] = xmus2 + xmuv2;
    pl[7] = xlntau * pl[6];
    pl[8] = xmus2 * xmuv2;
    pl[9] = xlntau * pl[8];

    fs0 = 0.0;
    for (i = 0; i < 10; i++)
        fs0 += pl[i] * as0[i];
    fs1 = pl[0] * as1[0] + pl[1] * as1[1];
    fs2 = pl[0] * as2[0] + pl[1] * as2[1];
    xitot1 = xp1 + cfonc1 * fs0 * xmus;
    xitot2 = xp2 + cfonc2 * fs1 * xmus;
    xitot3 = xp3 + cfonc3 * fs2 * xmus;

    *xrray = xitot1;
    *xrray += xitot2 * xcosf2 * 2.0;
    *xrray += xitot3 * xcosf3 * 2.0;
    *xrray /= xmus;
}


/******************************************************************************
MODULE:  comptg

PURPOSE:  Computes the transmission of the water vapor, ozone, and other gases.

RETURN VALUE:
Type = N/A

NOTES:
1. Standard sea level pressure is 1013 millibars.
******************************************************************************/
void comptg
(
    int iband,                   /* I: band index (0-based) */
    float xts,                   /* I: solar zenith angle */
    float xtv,                   /* I: view zenith angle */
    float xmus,                  /* I: cosine of solar zenith angle */
    float xmuv,                  /* I: cosine of view zenith angle */
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
)
{
    float a, b;  /* water vapor transmission coefficient */
    float m;     /* ozone transmission coefficient */
    float x;     /* water vapor transmission coefficient */

    /* Compute ozone transmission */
    m = 1.0 / xmus + 1.0 / xmuv;
    *tgoz = exp(oztransa[iband] * m * uoz);

    /* Compute water vapor transmission */
    a = wvtransa[iband];
    b = wvtransb[iband];

    x = m * uwv;
    if (x > 1.0E-06)
        *tgwv = exp(-a * exp(log(x) * b));
    else
        *tgwv = 1.0;

    /* Compute water vapor transmission half the content */
    x *= 0.5;
    if (x > 1.0E-06)
        *tgwvhalf = exp(-a * exp(log(x) * b));
    else
        *tgwvhalf = 1.0;

    /* Compute other gases transmission */
    *tgog = -(ogtransa1[iband] * atm_pres) *
        pow(m, exp(-(ogtransb0[iband] + ogtransb1[iband] * atm_pres)));
    *tgog = exp(*tgog);
}


/******************************************************************************
MODULE:  compsalb

PURPOSE:  Computes spherical albedo

RETURN VALUE:
Type = N/A

NOTES:
******************************************************************************/
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
    float *next         /* O: */
)
{
    float xtiaot1, xtiaot2;         /* spherical albedo trans value */
    float satm1, satm2;             /* spherical albedo value */
    float next1, next2;
    float dpres;                    /* pressure ratio */
    float deltaaot;                 /* AOT ratio */

    /* Compute the delta AOT */
    deltaaot = raot550nm - aot550nm[iaot1];
    deltaaot /= aot550nm[iaot2] - aot550nm[iaot1];

    /* Compute the spherical albedo */
    xtiaot1 = sphalbt[iband][ip1][iaot1];
    xtiaot2 = sphalbt[iband][ip1][iaot2];
    satm1 = xtiaot1 + (xtiaot2 - xtiaot1) * deltaaot;

    xtiaot1 = sphalbt[iband][ip2][iaot1];
    xtiaot2 = sphalbt[iband][ip2][iaot2];
    satm2 = xtiaot1 + (xtiaot2 - xtiaot1) * deltaaot;

    dpres = (pres - tpres[ip1]) / (tpres[ip2] - tpres[ip1]);
    *satm = satm1 + (satm2 - satm1) * dpres;

    /* Compute the normalized?? spherical albedo */
    xtiaot1 = normext[iband][ip1][iaot1];
    xtiaot2 = normext[iband][ip1][iaot2];
    next1 = xtiaot1 + (xtiaot2 - xtiaot1) * deltaaot;

    xtiaot1 = normext[iband][ip2][iaot1];
    xtiaot2 = normext[iband][ip2][iaot2];
    next2 = xtiaot1 + (xtiaot2 - xtiaot1) * deltaaot;

    dpres = (pres - tpres[ip1]) / (tpres[ip2] - tpres[ip1]);
    *next = next1 + (next2 - next1) * dpres;
}


/******************************************************************************
MODULE:  comptrans

PURPOSE:  Compute transmission

RETURN VALUE:
Type = none

NOTES:
1. This is called by subaeroret for both the solar zenith angle and the
   observation zenith angle.  Thus, xts is not specific to the solar zenith
   angle in this case.
2. This function is heavily dependent upon the input solar zenith and
   observation zenith angles.  At the current time, these are static values
   in the overall application.  Knowing that, speedup is achievable in this
   routine if these values never change.  However, the long-term goal is to
   change the main function to compute the solar zenith angle on a per-pixel
   basis.  Therefore we will leave this routine as-is.
3. This function is also dependent upon surface pressure and AOT.
******************************************************************************/
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
)
{
    char FUNC_NAME[] = "comptrans";   /* function name */
    char errmsg[STR_SIZE];            /* error message */
    float xtiaot1, xtiaot2;         /* spherical albedo trans value */
    float xtts1, xtts2;
    float xmts, xtranst;
    float dpres;                    /* pressure ratio */
    float deltaaot;                 /* AOT ratio */
    int its;                        /* index for the sun angle table */

    /* Determine the index in the sun angle table */
    if (xts <= xtsmin) 
        its = 0;
    else
        its = (int) ((xts - xtsmin) / xtsstep);
    if (its > 19)
    {
        sprintf (errmsg, "Zenith angle (xts) is too large: %f", xts);
        error_handler (true, FUNC_NAME, errmsg);
        return;
    }

    xmts = (xts - tts[its]) * 0.25;
    xtranst = transt[iband][ip1][iaot1][its];
    xtiaot1 = xtranst + (transt[iband][ip1][iaot1][its+1] - xtranst) * xmts;

    xtranst = transt[iband][ip1][iaot2][its];
    xtiaot2 = xtranst + (transt[iband][ip1][iaot2][its+1] - xtranst) * xmts;

    deltaaot = raot550nm - aot550nm[iaot1];
    deltaaot /= aot550nm[iaot2] - aot550nm[iaot1];
    xtts1 = xtiaot1 + (xtiaot2 - xtiaot1) * deltaaot;

    xtranst = transt[iband][ip2][iaot1][its];
    xtiaot1 = xtranst + (transt[iband][ip2][iaot1][its+1] - xtranst) * xmts;

    xtranst = transt[iband][ip2][iaot2][its];
    xtiaot2 = xtranst + (transt[iband][ip2][iaot2][its+1] - xtranst) * xmts;
    xtts2 = xtiaot1 + (xtiaot2 - xtiaot1) * deltaaot;

    dpres = (pres - tpres[ip1]) / (tpres[ip2] - tpres[ip1]);
    *xtts = xtts1 + (xtts2 - xtts1) * dpres;
}


/******************************************************************************
MODULE:  comproatm

PURPOSE:  Computes the atmospheric reflectance

RETURN VALUE:
Type = none

NOTES:
1. This function is heavily dependent upon the solar zenith and observation
   zenith angles.  At the current time, these are static values in the
   overall application.  Knowing that, speedup is achievable in this routine
   if these values never change.  However, the long-term goal is to change
   the main function to compute the solar zenith angle on a per-pixel basis.
   Therefore we will leave this routine as-is.
2. This function is also dependent upon surface pressure and AOT.
******************************************************************************/
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
)
{
    int isca;
    int iindex;
    float nbfic1, nbfic2, nbfic3, nbfic4;
    float nbfi1, nbfi2, nbfi3, nbfi4;
    float ro, rop1, rop2;           /* reflectance at p1 and p2 */
    float xtsmax;
    float cscaa;
    float scaa;                     /* scattering angle */
    float sca1;
    float sca2;
    float dpres;                    /* pressure ratio */
    float deltaaot;                 /* AOT ratio */
    float roinf;
    float rosup;
    float ro1, ro2, ro3, ro4;
    float roiaot1, roiaot2;
    float t, u;
    float logaot550nm[22] =
        {-4.605170186, -2.995732274, -2.302585093,
         -1.897119985, -1.609437912, -1.203972804,
         -0.916290732, -0.510825624, -0.223143551,
          0.000000000, 0.182321557, 0.336472237,
          0.470003629, 0.587786665, 0.693157181,
          0.832909123, 0.955511445, 1.098612289,
          1.252762969, 1.386294361, 1.504077397,
          1.609437912};

    cscaa = -xmus * xmuv - cosxfi * sqrt(1.0 - xmus * xmus) *
        sqrt(1.0 - xmuv * xmuv);
    scaa = acos(cscaa) * RAD2DEG;    /* vs / DEG2RAD */

    nbfic1 = nbfic[itv][its];
    nbfi1 = nbfi[itv][its];
    nbfic2 = nbfic[itv][its+1];
    nbfi2 = nbfi[itv][its+1];
    nbfic3 = nbfic[itv+1][its];
    nbfi3 = nbfi[itv+1][its];
    nbfic4 = nbfic[itv+1][its+1];
    nbfi4 = nbfi[itv+1][its+1];

    /* Compute for ip1, iaot1 */
    /* Interpolate point 1 (its,itv) vs scattering angle */
    xtsmax = tsmax[itv][its];
    if ((its != 0) && (itv != 0))
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi1)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi1 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv][its];
        }

        iindex = indts[its] + nbfic1 - nbfi1 + isca - 1;
        roinf = rolutt[iband][ip1][iaot1][iindex];
        rosup = rolutt[iband][ip1][iaot1][iindex+1];
        ro1 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its] + nbfic1 - nbfi1;
        roinf = rolutt[iband][ip1][iaot1][iindex];
        rosup = roinf;
        ro1 = roinf;
    }

    /* Interpolate point 2 (its+1,itv) vs scattering angle */
    xtsmax = tsmax[itv][its+1];
    if (itv != 0)
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi2)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi2 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv][its+1];
        }

        iindex = indts[its+1] + nbfic2 - nbfi2 + isca - 1;
        roinf = rolutt[iband][ip1][iaot1][iindex];
        rosup = rolutt[iband][ip1][iaot1][iindex+1];
        ro2 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its+1] + nbfic2 - nbfi2;
        roinf = rolutt[iband][ip1][iaot1][iindex];
        rosup = roinf;
        ro2 = roinf;
    }

    /* Interpolate point 3 (its,itv+1) vs scattering angle */
    xtsmax = tsmax[itv+1][its];
    if (its != 0)
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi3)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi3 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv+1][its];
        }

        iindex = indts[its] + nbfic3 - nbfi3 + isca - 1;
        roinf = rolutt[iband][ip1][iaot1][iindex];
        rosup = rolutt[iband][ip1][iaot1][iindex+1];
        ro3 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its] + nbfic3 - nbfi3;
        roinf = rolutt[iband][ip1][iaot1][iindex];
        rosup = roinf;
        ro3 = roinf;
    }

    /* Interpolate point 4 (its+1,itv+1) vs scattering angle */
    xtsmax = tsmax[itv+1][its+1];
    isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
    if (isca <= 0)
        isca = 1;
    if (isca + 1 < nbfi4)
    {
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = xtsmax - isca * 4.0;
    }
    else
    {
        isca = nbfi4 - 1;
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = tsmin[itv+1][its+1];
    }

    iindex = indts[its+1] + nbfic4 - nbfi4 + isca - 1;
    roinf = rolutt[iband][ip1][iaot1][iindex];
    rosup = rolutt[iband][ip1][iaot1][iindex+1];
    ro4 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);

    /* Note: t and u are used elsewhere through this function */
    t = (tts[its+1] - xts) / (tts[its+1] - tts[its]);
    u = (ttv[itv+1][its] - xtv) / (ttv[itv+1][its] - ttv[itv][its]);
    roiaot1 = ro1 * t * u + ro2 * u * (1.0 - t) + ro3 * (1.0 - u) * t +
        ro4 * (1.0 - u) * (1.0 - t);

    /* Compute for ip1, iaot2 */
    /* Interpolate point 1 (its,itv) vs scattering angle */
    xtsmax = tsmax[itv][its];
    if ((its != 0) && (itv != 0))
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi1)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi1 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv][its];
        }

        iindex = indts[its] + nbfic1 - nbfi1 + isca - 1;
        roinf = rolutt[iband][ip1][iaot2][iindex];
        rosup = rolutt[iband][ip1][iaot2][iindex+1];
        ro1 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its] + nbfic1 - nbfi1;
        roinf = rolutt[iband][ip1][iaot2][iindex];
        rosup = roinf;
        ro1 = roinf;
    }

    /* Interpolate point 2 (its+1,itv) vs scattering angle */
    xtsmax = tsmax[itv][its+1];
    if (itv != 0)
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi2)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi2 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv][its+1];
        }

        iindex = indts[its+1] + nbfic2 - nbfi2 + isca - 1;
        roinf = rolutt[iband][ip1][iaot2][iindex];
        rosup = rolutt[iband][ip1][iaot2][iindex+1];
        ro2 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its+1] + nbfic2 - nbfi2;
        roinf = rolutt[iband][ip1][iaot2][iindex];
        rosup = roinf;
        ro2 = roinf;
    }

    /* Interpolate point 3 (its,itv+1) vs scattering angle */
    xtsmax = tsmax[itv+1][its];
    if (its != 0)
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi3)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi3 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv+1][its];
        }

        iindex = indts[its] + nbfic3 - nbfi3 + isca - 1;
        roinf = rolutt[iband][ip1][iaot2][iindex];
        rosup = rolutt[iband][ip1][iaot2][iindex+1];
        ro3 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its] + nbfic3 - nbfi3;
        roinf = rolutt[iband][ip1][iaot2][iindex];
        rosup = roinf;
        ro3 = roinf;
    }

    /* Interpolate point 4 (its+1,itv+1) vs scattering angle */
    xtsmax = tsmax[itv+1][its+1];
    isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
    if (isca <= 0)
        isca = 1;
    if (isca + 1 < nbfi4)
    {
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = xtsmax - isca * 4.0;
    }
    else
    {
        isca = nbfi4 - 1;
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = tsmin[itv+1][its+1];
    }

    iindex = indts[its+1] + nbfic4 - nbfi4 + isca - 1;
    roinf = rolutt[iband][ip1][iaot2][iindex];
    rosup = rolutt[iband][ip1][iaot2][iindex+1];
    ro4 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);

    roiaot2 = ro1 * t * u + ro2 * u * (1.0 - t) + ro3 * (1.0 - u) * t +
        ro4 * (1.0 - u) * (1.0 - t);

    /* Interpolation as log of tau */
    /* Note: delaaot is calculated here and used later in this function */
    deltaaot = logaot550nm[iaot2] - logaot550nm[iaot1];
    deltaaot = (log (raot550nm) - logaot550nm[iaot1]) / deltaaot;
    ro = roiaot1 + (roiaot2 - roiaot1) * deltaaot;
    rop1 = ro;

    /* Compute for ip2, iaot1 */
    /* Interpolate point 1 (its,itv) vs scattering angle */
    xtsmax = tsmax[itv][its];
    if ((its != 0) && (itv != 0))
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi1)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi1 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv][its];
        }

        iindex = indts[its] + nbfic1 - nbfi1 + isca - 1;
        roinf = rolutt[iband][ip2][iaot1][iindex];
        rosup = rolutt[iband][ip2][iaot1][iindex+1];
        ro1 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its] + nbfic1 - nbfi1;
        roinf = rolutt[iband][ip2][iaot1][iindex];
        rosup = roinf;
        ro1 = roinf;
    }

    /* Interpolate point 2 (its+1,itv) vs scattering angle */
    xtsmax = tsmax[itv][its+1];
    if (itv != 0)
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi2)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi2 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv][its+1];
        }

        iindex = indts[its+1] + nbfic2 - nbfi2 + isca - 1;
        roinf = rolutt[iband][ip2][iaot1][iindex];
        rosup = rolutt[iband][ip2][iaot1][iindex+1];
        ro2 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its+1] + nbfic2 - nbfi2;
        roinf = rolutt[iband][ip2][iaot1][iindex];
        rosup = roinf;
        ro2 = roinf;
    }

    /* Interpolate point 3 (its,itv+1) vs scattering angle */
    xtsmax = tsmax[itv+1][its];
    if (its != 0)
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi3)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi3 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv+1][its];
        }

        iindex = indts[its] + nbfic3 - nbfi3 + isca - 1;
        roinf = rolutt[iband][ip2][iaot1][iindex];
        rosup = rolutt[iband][ip2][iaot1][iindex+1];
        ro3 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its] + nbfic3 - nbfi3;
        roinf = rolutt[iband][ip2][iaot1][iindex];
        rosup = roinf;
        ro3 = roinf;
    }

    /* Interpolate point 4 (its+1,itv+1) vs scattering angle */
    xtsmax = tsmax[itv+1][its+1];
    isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
    if (isca <= 0)
        isca = 1;
    if (isca + 1 < nbfi4)
    {
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = xtsmax - isca * 4.0;
    }
    else
    {
        isca = nbfi4 - 1;
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = tsmin[itv+1][its+1];
    }

    iindex = indts[its+1] + nbfic4 - nbfi4 + isca - 1;
    roinf = rolutt[iband][ip2][iaot1][iindex];
    rosup = rolutt[iband][ip2][iaot1][iindex+1];
    ro4 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);

    roiaot1 = ro1 * t * u + ro2 * u * (1.0 - t) + ro3 * (1.0 - u) * t +
        ro4 * (1.0 - u) * (1.0 - t);

    /* Compute for ip2, iaot2 */
    /* Interpolate point 1 (its,itv) vs scattering angle */
    xtsmax = tsmax[itv][its];
    if ((its != 0) && (itv != 0))
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi1)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi1 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv][its];
        }

        iindex = indts[its] + nbfic1 - nbfi1 + isca - 1;
        roinf = rolutt[iband][ip2][iaot2][iindex];
        rosup = rolutt[iband][ip2][iaot2][iindex+1];
        ro1 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its] + nbfic1 - nbfi1;
        roinf = rolutt[iband][ip2][iaot2][iindex];
        rosup = roinf;
        ro1 = roinf;
    }

    /* Interpolate point 2 (its+1,itv) vs scattering angle */
    xtsmax = tsmax[itv][its+1];
    if (itv != 0)
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi2)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi2 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv][its+1];
        }

        iindex = indts[its+1] + nbfic2 - nbfi2 + isca - 1;
        roinf = rolutt[iband][ip2][iaot2][iindex];
        rosup = rolutt[iband][ip2][iaot2][iindex+1];
        ro2 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its+1] + nbfic2 - nbfi2;
        roinf = rolutt[iband][ip2][iaot2][iindex];
        rosup = roinf;
        ro2 = roinf;
    }

    /* Interpolate point 3 (its,itv+1) vs scattering angle */
    xtsmax = tsmax[itv+1][its];
    if (its != 0)
    {
        isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
        if (isca <= 0)
            isca = 1;
        if (isca + 1 < nbfi3)
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
            isca = nbfi3 - 1;
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = tsmin[itv+1][its];
        }

        iindex = indts[its] + nbfic3 - nbfi3 + isca - 1;
        roinf = rolutt[iband][ip2][iaot2][iindex];
        rosup = rolutt[iband][ip2][iaot2][iindex+1];
        ro3 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);
    }
    else
    {
        sca1 = xtsmax;
        sca2 = xtsmax;
        iindex = indts[its] + nbfic3 - nbfi3;
        roinf = rolutt[iband][ip2][iaot2][iindex];
        rosup = roinf;
        ro3 = roinf;
    }

    /* Interpolate point 4 (its+1,itv+1) vs scattering angle */
    xtsmax = tsmax[itv+1][its+1];
    isca = (int) ((xtsmax - scaa) * 0.25 + 1);   /* * 0.25 vs / 4.0 */
    if (isca <= 0)
        isca = 1;
    if (isca + 1 < nbfi4)
    {
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = xtsmax - isca * 4.0;
    }
    else
    {
        isca = nbfi4 - 1;
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = tsmin[itv+1][its+1];
    }

    iindex = indts[its+1] + nbfic4 - nbfi4 + isca - 1;
    roinf = rolutt[iband][ip2][iaot2][iindex];
    rosup = rolutt[iband][ip2][iaot2][iindex+1];
    ro4 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);

    roiaot2 = ro1 * t * u + ro2 * u * (1.0 - t) + ro3 * (1.0 - u) * t +
        ro4 * (1.0 - u) * (1.0 - t);

    /* Interpolation as log of tau */
    ro = roiaot1 + (roiaot2 - roiaot1) * deltaaot;
    rop2 = ro;

    dpres = (pres - tpres[ip1]) / (tpres[ip2] - tpres[ip1]);
    *roatm = rop1 + (rop2 - rop1) * dpres;
}


/******************************************************************************
MODULE:  readluts

PURPOSE:  Reads the look-up tables and input atmospheric files

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the look-up tables or atmospheric files
SUCCESS        Successful completion

NOTES:
******************************************************************************/
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
    float ***normext,           /* O: aerosol extinction coefficient at the
                                      current wavelength (normalized at 550nm)
                                      [NSR_BANDS][7][22] */
    float xtsstep,              /* I: solar zenith step value */
    float xtsmin,               /* I: minimum solar zenith value */
    char anglehdf[STR_SIZE],    /* I: angle HDF filename */
    char intrefnm[STR_SIZE],    /* I: intrinsic reflectance filename */
    char transmnm[STR_SIZE],    /* I: transmission filename */
    char spheranm[STR_SIZE]     /* I: spherical albedo filename */
)
{
    char FUNC_NAME[] = "readluts";   /* function name */
    char errmsg[STR_SIZE];  /* error message */
    char tmpstr[STR_SIZE];  /* temporary string variable, not use */
    int i, j;               /* looping variables */
    int iband;              /* band looping variable */
    int iaot;               /* aerosol optical thickness (AOT) index */
    int ipres;              /* looping variable for pressure */
    int itau;               /* looping variable for molecular optical thick */
    int ival;               /* looping variable for LUT */
    int status;             /* return status of the HDF function */
    int start[3];           /* starting point to read SDS data */
    int edges[3];           /* number of values to read in SDS data */
    char fname[STR_SIZE];   /* filename to be read */
    float *rolut = NULL;    /* intrinsic reflectance read from HDF file
                               [8000*22*7] */
    float ttsr[22];        /* GAIL - should this be 21 instead?? */
    float xx;               /* temporary float values, not used */
    int sd_id;              /* file ID for the HDF file */
    int sds_id;             /* ID for the current SDS */
    int sds_index;          /* index for the current SDS */
    FILE *fp = NULL;        /* file pointer for reading ascii files */

    /* Initialize some variables */
    for (i = 0; i < 20; i++)
    {
        for (j = 0; j < 22; j++)
        {
            nbfic[i][j] = 0.0;
            tts[j] = xtsmin + xtsstep * j;
        }
    }

    /* Open as HDF file for reading */
    sd_id = SDstart (anglehdf, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", anglehdf);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the 2D bands from the angle HDF file */
    start[0] = 0;   /* lines */
    start[1] = 0;   /* samples */
    edges[0] = 20;  /* number of lines */
    edges[1] = 22;  /* number of samples */

    /* Find the TSMAX SDS */
    sds_index = SDnametoindex (sd_id, "TSMAX");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find TSMAX in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access TSMAX for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    edges[0] = 1;   /* number of lines */
    for (i = 0; i < 20; i++)
    {
        start[0] = i;   /* lines */
        status = SDreaddata (sds_id, start, NULL, edges, tsmax[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the TSMAX SDS");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to TSMAX SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the TSMIN SDS */
    sds_index = SDnametoindex (sd_id, "TSMIN");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find TSMIN in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access TSMIN for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    edges[0] = 1;   /* number of lines */
    for (i = 0; i < 20; i++)
    {
        start[0] = i;   /* lines */
        status = SDreaddata (sds_id, start, NULL, edges, tsmin[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the TSMIN SDS");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to TSMIN SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the TTV SDS */
    sds_index = SDnametoindex (sd_id, "TTV");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find TTV in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access TTV for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    edges[0] = 1;   /* number of lines */
    for (i = 0; i < 20; i++)
    {
        start[0] = i;   /* lines */
        status = SDreaddata (sds_id, start, NULL, edges, ttv[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the TTV SDS");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to TTV SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the NBFI SDS */
    sds_index = SDnametoindex (sd_id, "NBFI");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find NBFI in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access NBFI for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    edges[0] = 1;   /* number of lines */
    for (i = 0; i < 20; i++)
    {
        start[0] = i;   /* lines */
        status = SDreaddata (sds_id, start, NULL, edges, nbfi[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the NBFI SDS");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to NBFI SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the NBFIC SDS */
    sds_index = SDnametoindex (sd_id, "NBFIC");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find NBFIC in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access NBFIC for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    edges[0] = 1;   /* number of lines */
    for (i = 0; i < 20; i++)
    {
        start[0] = i;   /* lines */
        status = SDreaddata (sds_id, start, NULL, edges, nbfic[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the NBFIC SDS");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to NBFIC SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the INDTS SDS */
    sds_index = SDnametoindex (sd_id, "INDTS");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find INDTS in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access INDTS for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    start[0] = 0;   /* lines */
    start[1] = 0;   /* samples */
    edges[0] = 20;  /* number of lines */
    edges[1] = 22;  /* number of samples */
    status = SDreaddata (sds_id, start, NULL, edges, indts);
    if (status == -1)
    {
        sprintf (errmsg, "Reading data from the INDTS SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to INDTS SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the TTS SDS */
    sds_index = SDnametoindex (sd_id, "TTS");
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find TTS in the HDF file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access TTS for reading");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    start[0] = 0;   /* lines */
    start[1] = 0;   /* samples */
    edges[0] = 20;  /* number of lines */
    edges[1] = 22;  /* number of samples */
    status = SDreaddata (sds_id, start, NULL, edges, tts);
    if (status == -1)
    {
        sprintf (errmsg, "Reading data from the TTS SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the HDF SDS */
    status = SDendaccess (sds_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to TTS SDS");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the HDF file */
    status = SDend (sd_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to HDF file: %s", anglehdf);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* rolutt[8000][22][7] */
    rolut = calloc (8000*22*7, sizeof (float));
    if (rolut == NULL)
    {
        sprintf (errmsg, "Error allocating memory for rolut");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Begin read look up table (intrinsic reflectance) */
    /* Open as HDF file for reading */
    sd_id = SDstart (intrefnm, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", intrefnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    start[0] = 0;  /* left-most dimension */
    start[1] = 0;
    start[2] = 0;  /* right-most dimension */
    edges[0] = 8000;
    edges[1] = 22;
    edges[2] = 7;
    for (iband = 0; iband < NSR_BANDS; iband++)
    {
        /* Get the sds name */
        sprintf (fname, "NRLUT_BAND_%d", iband+1);

        /* Find the SDS */
        sds_index = SDnametoindex (sd_id, fname);
        if (sds_index == -1)
        {
            sprintf (errmsg, "Unable to find %s in the %s HDF file", fname,
                intrefnm);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    
        /* Open the current band as an SDS */
        sds_id = SDselect (sd_id, sds_index);
        if (sds_id < 0)
        {
            sprintf (errmsg, "Unable to access %s for reading", fname);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    
        /* Read the whole band, as-is, then rearrange the order later */
        status = SDreaddata (sds_id, start, NULL, edges, rolut);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the %s SDS", fname);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    
        /* Close the HDF SDS */
        status = SDendaccess (sds_id);
        if (status == -1)
        {
            sprintf (errmsg, "Ending access to %s SDS", fname);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Flip the look-up table value to allow the index values to be in
           the fastest-changing dimension and the band values to be in the
           slowest-changing dimension.  This allows the access for the data
           in the application to be most efficient.  rolut is read from the
           HDF file (per band) as 8000 x 22 x 7. */
        for (ipres = 0; ipres < 7; ipres++)
            for (itau = 0; itau < 22; itau++)
                for (ival = 0; ival < 8000; ival++)
                    rolutt[iband][ipres][itau][ival] =
                        rolut[ival*22*7 + itau*7 + ipres];
    }  /* for iband */

    /* Free the temporary rolut array */
    free (rolut);

    /* Close the HDF file */
    status = SDend (sd_id);
    if (status == -1)
    {
        sprintf (errmsg, "Ending access to HDF file: %s", intrefnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Begin read look up table (transmission) */
    fp = fopen (transmnm, "r");
    if (fp == NULL)
    {
        sprintf (errmsg, "Opening transmission coefficient file: %s", transmnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* 8 bands of data */
    for (iband = 0; iband < NSR_BANDS; iband++)
    {
        /* This first read contains information about the band and source of
           the data; ignore for now */
        if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
        {
            sprintf (errmsg, "Skipping data source in transmission data file");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* 7 pressure levels (1050.0 mb, 1013.0 mb, 900.0 mb, 800.0 mb,
           700.0, 600.0 mb, 500.0 mb) */
        for (ipres = 0; ipres < 7; ipres++)
        {
            /* This next read contains information about the pressure level of
               the data; ignore for now */
            if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
            {
                sprintf (errmsg, "Skipping pressure level in transmission data "
                    "file");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* 21 lines per pressure level */
            for (i = 0; i < 21; i++)   /* GAIL - array is size 22 vs. 21 */
            {
                /* Grab the first value in the line.  Basically this is a
                   repeat of the previous pressure level and band, as all the
                   pressure levels have the same values for the ttsr. */
                if (fscanf (fp, "%f ", &ttsr[i]) != 1)
                {
                    sprintf (errmsg, "Reading first transmission value from "
                        "transmission coefficient file: %s", transmnm);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                if (fabs (tts[i] - ttsr[i]) > 1.0E-5)
                {
                    sprintf (errmsg, "Problem with transmission LUT: %s",
                        transmnm);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                /* Grab the remaining 22 values in the line.  Store the iaot
                   in the more efficient order for processing, not necessarily
                   for reading. */
                for (iaot = 0; iaot < 22; iaot++)
                {
                    if (fscanf (fp, "%f", &transt[iband][ipres][iaot][i]) != 1)
                    {
                        sprintf (errmsg, "Reading transmission values from "
                            "transmission coefficient file: %s", transmnm);
                        error_handler (true, FUNC_NAME, errmsg);
                        return (ERROR);
                    }
                }
            }  /* for i */

            /* Clear out the EOL for the last line */
            if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
            {
                sprintf (errmsg, "Skipping EOL in last line");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
        }  /* for ipres */
    }  /* for iband */

    /* Close transmission file */
    fclose (fp);

    /* Begin read look up table (spherical albedo) */
    fp = fopen (spheranm, "r");
    if (fp == NULL)
    {
        sprintf (errmsg, "Opening spherical albedo coefficient file: %s",
            spheranm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* 8 bands of data */
    for (iband = 0; iband < NSR_BANDS; iband++)
    {
        /* This first read contains information about the source of the data;
           ignore for now */
        if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
        {
            sprintf (errmsg, "Skipping data source in spherical albedo data "
                "file");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* 7 pressure levels (1050.0 mb, 1013.0 mb, 900.0 mb, 800.0 mb,
           700.0, 600.0 mb, 500.0 mb) */
        for (ipres = 0; ipres < 7; ipres++)
        {
            /* This next read contains information about the pressure level of
               the data; ignore for now */
            if (fgets (tmpstr, sizeof (tmpstr), fp) == NULL)
            {
                sprintf (errmsg, "Skipping pressure level in spherical albedo "
                    "data file");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* 22 lines of spherical albedo information */
            for (iaot = 0; iaot < 22; iaot++)
            {
                if (fscanf (fp, "%f %f %f\n", &xx, &sphalbt[iband][ipres][iaot],
                    &normext[iband][ipres][iaot]) != 3)
                {
                    sprintf (errmsg, "Reading spherical albedo values from "
                        "spherical albedo coefficient file: %s", spheranm);
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
            }
        }
    }

    /* Close spherical albedo file */
    fclose (fp);

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  memory_allocation_main

PURPOSE:  Allocates memory for all the various arrays within the L8 surface
reflectance application for the main application.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred allocating memory
SUCCESS        Successful completion

NOTES:
  1. Memory is allocated for each of the input variables, so it is up to the
     calling routine to free this memory.
  2. Each array passed into this function is passed in as the address to that
     1D, 2D, nD array.
******************************************************************************/
int memory_allocation_main
(
    int nlines,          /* I: number of lines in the scene */
    int nsamps,          /* I: number of samples in the scene */
    uint16 **qaband,     /* O: QA band for the input image, nlines x nsamps */
    int16 ***sband       /* O: output surface reflectance and brightness temp
                               bands */
)
{
    char FUNC_NAME[] = "memory_allocation_main"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int i;                   /* looping variables */

    *qaband = calloc (nlines*nsamps, sizeof (uint16));
    if (*qaband == NULL)
    {
        sprintf (errmsg, "Error allocating memory for qaband");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Given that the QA band is its own separate array of uint16s, we need
       one less band for the signed image data */
    *sband = calloc (NBAND_TTL_OUT-1, sizeof (int16*));
    if (*sband == NULL)
    {
        sprintf (errmsg, "Error allocating memory for sband");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    for (i = 0; i < NBAND_TTL_OUT-1; i++)
    {
        (*sband)[i] = calloc (nlines*nsamps, sizeof (int16));
        if ((*sband)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for sband");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  memory_allocation_sr

PURPOSE:  Allocates memory for all the various arrays needed specifically for
the L8 surface reflectance corrections.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred allocating memory
SUCCESS        Successful completion

NOTES:
  1. Memory is allocated for each of the input variables, so it is up to the
     calling routine to free this memory.
  2. Each array passed into this function is passed in as the address to that
     1D, 2D, nD array.
******************************************************************************/
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
)
{
    char FUNC_NAME[] = "memory_allocation_sr"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    int i, j, k;             /* looping variables */

    *aerob1 = calloc (nlines*nsamps, sizeof (int16));
    if (*aerob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob1");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *aerob2 = calloc (nlines*nsamps, sizeof (int16));
    if (*aerob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob2");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *aerob4 = calloc (nlines*nsamps, sizeof (int16));
    if (*aerob4 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob4");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *aerob5 = calloc (nlines*nsamps, sizeof (int16));
    if (*aerob5 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob5");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *aerob7 = calloc (nlines*nsamps, sizeof (int16));
    if (*aerob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for aerob7");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *twvi = calloc (nlines*nsamps, sizeof (float));
    if (*twvi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for twvi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *tozi = calloc (nlines*nsamps, sizeof (float));
    if (*tozi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tozi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *tp = calloc (nlines*nsamps, sizeof (float));
    if (*tp == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tp");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *tresi = calloc (nlines*nsamps, sizeof (float));
    if (*tresi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tresi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *taero = calloc (nlines*nsamps, sizeof (float));
    if (*taero == NULL)
    {
        sprintf (errmsg, "Error allocating memory for taero");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *cloud = calloc (nlines*nsamps, sizeof (uint8));
    if (*cloud == NULL)
    {
        sprintf (errmsg, "Error allocating memory for cloud");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *lw_mask = calloc (nlines*nsamps, sizeof (uint8));
    if (*lw_mask == NULL)
    {
        sprintf (errmsg, "Error allocating memory for lw_mask");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Allocate memory for all the climate modeling grid files */
    *dem = calloc (DEM_NBLAT, sizeof (int16*));
    if (*dem == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the DEM");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    for (i = 0; i < DEM_NBLAT; i++)
    {
        (*dem)[i] = calloc (DEM_NBLON, sizeof (int16));
        if ((*dem)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the DEM");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    *andwi = calloc (RATIO_NBLAT, sizeof (int16*));
    if (*andwi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the andwi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *sndwi = calloc (RATIO_NBLAT, sizeof (int16*));
    if (*sndwi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the sndwi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ratiob1 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (*ratiob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob1");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ratiob2 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (*ratiob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob2");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ratiob7 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (*ratiob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the ratiob7");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *intratiob1 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (*intratiob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the intratiob1");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *intratiob2 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (*intratiob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the intratiob2");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *intratiob7 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (*intratiob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the intratiob7");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *slpratiob1 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (*slpratiob1 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the slpratiob1");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *slpratiob2 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (*slpratiob2 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the slpratiob2");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *slpratiob7 = calloc (RATIO_NBLAT, sizeof (int16*));
    if (*slpratiob7 == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the slpratiob7");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    for (i = 0; i < RATIO_NBLAT; i++)
    {
        (*andwi)[i] = calloc (RATIO_NBLON, sizeof (int16));
        if ((*andwi)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the andwi");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*sndwi)[i] = calloc (RATIO_NBLON, sizeof (int16));
        if ((*sndwi)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the sndwi");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*ratiob1)[i] = calloc (RATIO_NBLON, sizeof (int16));
        if ((*ratiob1)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the ratiob1");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*ratiob2)[i] = calloc (RATIO_NBLON, sizeof (int16));
        if ((*ratiob2)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the ratiob2");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*ratiob7)[i] = calloc (RATIO_NBLON, sizeof (int16));
        if ((*ratiob7)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the ratiob7");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*intratiob1)[i] = calloc (RATIO_NBLON, sizeof (int16));
        if ((*intratiob1)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the intratiob1");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*intratiob2)[i] = calloc (RATIO_NBLON, sizeof (int16));
        if ((*intratiob2)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the intratiob2");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*intratiob7)[i] = calloc (RATIO_NBLON, sizeof (int16));
        if ((*intratiob7)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the intratiob7");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*slpratiob1)[i] = calloc (RATIO_NBLON, sizeof (int16));
        if ((*slpratiob1)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the slpratiob1");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*slpratiob2)[i] = calloc (RATIO_NBLON, sizeof (int16));
        if ((*slpratiob2)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the slpratiob2");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*slpratiob7)[i] = calloc (RATIO_NBLON, sizeof (int16));
        if ((*slpratiob7)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the slpratiob7");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    *wv = calloc (CMG_NBLAT, sizeof (int16*));
    if (*wv == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the wv");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *oz = calloc (CMG_NBLAT, sizeof (uint8*));
    if (*oz == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the oz");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    for (i = 0; i < CMG_NBLAT; i++)
    {
        (*wv)[i] = calloc (CMG_NBLON, sizeof (int16));
        if ((*wv)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the wv");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*oz)[i] = calloc (CMG_NBLON, sizeof (uint8));
        if ((*oz)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for the oz");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* rolutt[NSR_BANDS][7][22][8000] and transt[NSR_BANDS][7][22][22] and
       sphalbt[NSR_BANDS][7][22] and normext[NSR_BANDS][7][22] */
    *rolutt = calloc (NSR_BANDS, sizeof (float***));
    if (*rolutt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for rolutt");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *transt = calloc (NSR_BANDS, sizeof (float***));
    if (*transt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for transt");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *sphalbt = calloc (NSR_BANDS, sizeof (float**));
    if (*sphalbt == NULL)
    {
        sprintf (errmsg, "Error allocating memory for sphalbt");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *normext = calloc (NSR_BANDS, sizeof (float**));
    if (*normext == NULL)
    {
        sprintf (errmsg, "Error allocating memory for normext");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    for (i = 0; i < NSR_BANDS; i++)
    {
        (*rolutt)[i] = calloc (7, sizeof (float**));
        if ((*rolutt)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for rolutt");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*transt)[i] = calloc (7, sizeof (float**));
        if ((*transt)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for transt");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*sphalbt)[i] = calloc (7, sizeof (float*));
        if ((*sphalbt)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for sphalbt");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*normext)[i] = calloc (7, sizeof (float*));
        if ((*normext)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for normext");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        for (j = 0; j < 7; j++)
        {
            (*rolutt)[i][j] = calloc (22, sizeof (float*));
            if ((*rolutt)[i][j] == NULL)
            {
                sprintf (errmsg, "Error allocating memory for rolutt");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            (*transt)[i][j] = calloc (22, sizeof (float*));
            if ((*transt)[i][j] == NULL)
            {
                sprintf (errmsg, "Error allocating memory for transt");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            (*sphalbt)[i][j] = calloc (22, sizeof (float));
            if ((*sphalbt)[i][j] == NULL)
            {
                sprintf (errmsg, "Error allocating memory for sphalbt");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            (*normext)[i][j] = calloc (22, sizeof (float));
            if ((*normext)[i][j] == NULL)
            {
                sprintf (errmsg, "Error allocating memory for normext");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            for (k = 0; k < 22; k++)
            {
                (*rolutt)[i][j][k] = calloc (8000, sizeof (float));
                if ((*rolutt)[i][j][k] == NULL)
                {
                    sprintf (errmsg, "Error allocating memory for rolutt");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                (*transt)[i][j][k] = calloc (22, sizeof (float));
                if ((*transt)[i][j][k] == NULL)
                {
                    sprintf (errmsg, "Error allocating memory for transt");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
            }  /* for k */
        }  /* for j */
    }  /* for i */

    /* tsmax[20][22] and float tsmin[20][22] and float nbfic[20][22] and
       nbfi[20][22] and float ttv[20][22] */
    *tsmax = calloc (20, sizeof (float*));
    if (*tsmax == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tsmax");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *tsmin = calloc (20, sizeof (float*));
    if (*tsmin == NULL)
    {
        sprintf (errmsg, "Error allocating memory for tsmin");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *nbfic = calloc (20, sizeof (float*));
    if (*nbfic == NULL)
    {
        sprintf (errmsg, "Error allocating memory for nbfic");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *nbfi = calloc (20, sizeof (float*));
    if (*nbfi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for nbfi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    *ttv = calloc (20, sizeof (float*));
    if (*ttv == NULL)
    {
        sprintf (errmsg, "Error allocating memory for ttv");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    for (i = 0; i < 20; i++)
    {
        (*tsmax)[i] = calloc (22, sizeof (float));
        if ((*tsmax)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for tsmax");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*tsmin)[i] = calloc (22, sizeof (float));
        if ((*tsmin)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for tsmin");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*nbfic)[i] = calloc (22, sizeof (float));
        if ((*nbfic)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for nbfic");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*nbfi)[i] = calloc (22, sizeof (float));
        if ((*nbfi)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for nbfi");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        (*ttv)[i] = calloc (22, sizeof (float));
        if ((*ttv)[i] == NULL)
        {
            sprintf (errmsg, "Error allocating memory for ttv");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  read_auxiliary_files

PURPOSE:  Reads the auxiliary files required for this application.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading one of the auxiliary files
SUCCESS        Successful completion

NOTES:
  1. It is assumed that memory has already been allocated for the input data
     arrays.
******************************************************************************/
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
)
{
    char FUNC_NAME[] = "read_auxiliary_files"; /* function name */
    char errmsg[STR_SIZE];   /* error message */
    char sds_name[STR_SIZE]; /* name of the SDS being read */
    int i, j;            /* looping variables */
    int status;          /* return status of the HDF function */
    int start[5];        /* starting point to read SDS data; handles up to
                            4D dataset */
    int edges[5];        /* number of values to read in SDS data; handles up to
                            4D dataset */
    int sd_id;           /* file ID for the HDF file */
    int sds_id;          /* ID for the current SDS */
    int sds_index;       /* index for the current SDS */

    /* Read the DEM */
    sd_id = SDstart (cmgdemnm, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", cmgdemnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "averaged elevation");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the DEM file %s", sds_name,
            cmgdemnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < DEM_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = DEM_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, dem[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the DEM file */
    status = SDend (sd_id);
    if (status != 0)
    {
        sprintf (errmsg, "Closing DEM file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the RATIO file */
    sd_id = SDstart (rationm, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", rationm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 6) */
    strcpy (sds_name, "average ndvi");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, andwi[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 14) */
    strcpy (sds_name, "standard ndvi");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, sndwi[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 21) */
    strcpy (sds_name, "slope ratiob9");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, slpratiob1[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 22) */
    strcpy (sds_name, "inter ratiob9");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, intratiob1[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 15) */
    strcpy (sds_name, "slope ratiob3");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, slpratiob2[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 16) */
    strcpy (sds_name, "inter ratiob3");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, intratiob2[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 27) */
    strcpy (sds_name, "slope ratiob7");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, slpratiob7[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name (SDS 28) */
    strcpy (sds_name, "inter ratiob7");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the RATIO file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = RATIO_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, intratiob7[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the RATIO file */
    status = SDend (sd_id);
    if (status != 0)
    {
        sprintf (errmsg, "Closing RATIO file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Compute the band ratios based on the averaged NDWI */
    for (i = 0; i < RATIO_NBLAT; i++)
    {
        for (j = 0; j < RATIO_NBLON; j++)
        {
            ratiob1[i][j] = (int16) (andwi[i][j] * slpratiob1[i][j] * 0.001 +
                intratiob1[i][j]);
            ratiob2[i][j] = (int16) (andwi[i][j] * slpratiob2[i][j] * 0.001 +
                intratiob2[i][j]);
            ratiob7[i][j] = (int16) (andwi[i][j] * slpratiob7[i][j] * 0.001 +
                intratiob7[i][j]);
        }
    }

    /* Read ozone and water vapor from the user-specified auxiliary file */
    sd_id = SDstart (auxnm, DFACC_RDONLY);
    if (sd_id < 0)
    {
        sprintf (errmsg, "Unable to open %s for reading as SDS", auxnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "Coarse Resolution Ozone");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the AUX file %s", sds_name,
            auxnm);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < CMG_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = CMG_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, oz[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find the SDS name */
    strcpy (sds_name, "Coarse Resolution Water Vapor");
    sds_index = SDnametoindex (sd_id, sds_name);
    if (sds_index == -1)
    {
        sprintf (errmsg, "Unable to find %s in the AUX file", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the current band as an SDS */
    sds_id = SDselect (sd_id, sds_index);
    if (sds_id < 0)
    {
        sprintf (errmsg, "Unable to access %s for reading", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Read the data one line at a time */
    for (i = 0; i < CMG_NBLAT; i++)
    {
        start[0] = i;  /* line */
        start[1] = 0;  /* sample */
        edges[0] = 1;
        edges[1] = CMG_NBLON;
        status = SDreaddata (sds_id, start, NULL, edges, wv[i]);
        if (status == -1)
        {
            sprintf (errmsg, "Reading data from the SDS: %s", sds_name);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Close the SDS */
    status = SDendaccess (sds_id);
    if (status < 0)
    {
        sprintf (errmsg, "Ending access to %s", sds_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the AUX file */
    status = SDend (sd_id);
    if (status != 0)
    {
        sprintf (errmsg, "Closing AUX file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Successful completion */
    return (SUCCESS);
}

