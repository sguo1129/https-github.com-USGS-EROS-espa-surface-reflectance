/*****************************************************************************
FILE: lut_subr.c
  
PURPOSE: Contains functions for reading the look-up tables and doing some
of the coefficient computations for the surface reflectance application.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
6/25/2014    Gail Schmidt     Conversion of the original FORTRAN code delivered
                              by Eric Vermote, NASA GSFC
7/15/2014    Gail Schmidt     Cleaned up some of the divide by constants so
                              that we could speed up the processing
7/21/2014    Gail Schmidt     Flipped the use of xmuv, xmus in the call to
                              local_chand to match the parameter order.  This
                              is a cosmetic change, however, since the nature
                              of the equations means them being flip-flopped
                              really doesn't matter.
7/22/2014    Gail Schmidt     Cleaned up unused ogtransa0, ogtransc0,
                              ogtransc1, wvtransc arrays.  Made the rest of
                              these transmission arrays doubles and hard-coded
                              their static values in this code vs. reading
                              from a static ASCII file.
7/29/2014    Gail Schmidt     Defined a static NSR_BANDS variable for the
                              variables that refer to the surface reflectance
                              band-related bands (ogtrans, wvtrans, tauray,
                              erelc, etc.)  These previously were of size 16.
8/7/2014     Gail Schmidt     Removed the unused functions (invaero,
                              invaeroocean, atmcorocea2, raycorlamb2,
                              atmcorlamb, local_csalbr, fintexp3, fintexp1,
                              comptransray).
8/14/2014    Gail Schmidt     Updated for v1.3 delivered by Eric Vermote

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
ERROR          Error occurred ???
SUCCESS        Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/25/2014    Gail Schmidt     Conversion of the original FORTRAN code delivered
                              by Eric Vermote, NASA GSFC

NOTES:
    1. Standard sea level pressure is 1013 millibars.
******************************************************************************/
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
    float *ttatmg,
    float *satm,                     /* O: spherical albedo */
    float *xrorayp,                  /* O: molecular reflectance */
    float *next                      /* O: ???? */
)
{
    char FUNC_NAME[] = "atmcorlamb2";   /* function name */
    char errmsg[STR_SIZE];  /* error message */
    int retval;             /* function return value */
    float xttv;    /* upward transmittance */
    float xtts;    /* downward transmittance */
    float ttatm;   /* total transmission of the atmosphere */
    float tgog;    /* other gases transmission */
    float tgoz;    /* ozone transmission */
    float tgwv;    /* water vapor transmission */
    float tgwvhalf;  /* water vapor transmission, half content */
    float xtaur;   /* rayleigh optical depth for surface pressure */
    float xphi;    /* azimuthal difference between sun and observation (deg) */
    float xmus;    /* cosine of solar zenith angle */
    float xmuv;    /* cosine of observation zenith angle */

    /* This routine returns variables for calculating roslamb */
    retval = comproatm (xts, xtv, xfi, raot550nm, iband, pres, tpres,
        aot550nm, rolutt, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, xtsstep,
        xtsmin, xtvstep, xtvmin, roatm);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Computing atmospheric reflectance.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    retval = comptrans (xts, raot550nm, iband, pres, tpres, aot550nm, transt,
        xtsstep, xtsmin, tts, &xtts);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Computing transmission.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    retval = comptrans (xtv, raot550nm, iband, pres, tpres, aot550nm, transt,
        xtvstep, xtvmin, tts, &xttv);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Computing transmission.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Compute total transmission (product downward by  upward) */
    ttatm = xtts * xttv;
    
    /* Compute spherical albedo */
    compsalb (raot550nm, iband, pres, tpres, aot550nm, sphalbt, normext, satm,
        next);

    comptg (iband, xts, xtv, uoz, uwv, pres, ogtransa1, ogtransb0, ogtransb1,
        wvtransa, wvtransb, oztransa, &tgoz, &tgwv, &tgwvhalf, &tgog);

    /* Compute rayleigh component (intrinsic reflectance, at p0) */
    xphi = xfi;
    xmus = cos (xts * DEG2RAD);
    xmuv = cos (xtv * DEG2RAD);

    /* Compute rayleigh component (intrinsic reflectance, at p=pres).
       Pressure in the atmosphere is pres / 1013. */
    xtaur = tauray[iband] * pres * ONE_DIV_1013;  /* vs / 1013.0 */
    local_chand (xphi, xmuv, xmus, xtaur, xrorayp);

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

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/27/2014    Gail Schmidt     Conversion of the original FORTRAN code delivered
                              by Eric Vermote, NASA GSFC

NOTES:
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
    static float xdep = 0.0279;        /* depolarization factor */
    float pl[10];
    float fs0, fs1, fs2;
    float phios;
    float xcosf1, xcosf2, xcosf3;
    float xbeta2;
    float xfd;
    float xph1, xph2, xph3;
    float xitm;
    float xp1, xp2, xp3;
    float cfonc1, cfonc2, cfonc3;
    float xlntau;                      /* log molecular optical depth */
    float xitot1, xitot2, xitot3;
    float as0[10] = {0.33243832, -6.777104e-02, 0.16285370, 1.577425e-03,
        -0.30924818, -1.240906e-02, -0.10324388, 3.241678e-02, 0.11493334,
        -3.503695e-02};
    float as1[2] = {0.19666292, -5.439061e-02};
    float as2[2] = {0.14545937, -2.910845e-02};

    phios = 180.0 - xphi;
    xcosf1 = 1.0;
    xcosf2 = cos (phios * DEG2RAD);
    xcosf3 = cos (2.0 * phios * DEG2RAD);
    xbeta2 = 0.5;

    /***** GAIL -- these are actually static .... so compute ahead of time ***/
    xfd = xdep / (2.0 - xdep);
    xfd = (1.0 - xfd) / (1.0 + 2.0 * xfd);

    xph1 = 1.0 + (3.0 * xmus * xmus -1.0) * (3.0 * xmuv * xmuv - 1.0) *
        xfd * 0.125;   /* vs / 8.0 */
    xph2 = -xmus * xmuv * sqrt(1.0 - xmus * xmus) * sqrt(1.0 - xmuv * xmuv);
    xph2 = xph2 * xfd * xbeta2 * 1.5;
    xph3 = (1.0 - xmus * xmus) * (1.0 - xmuv * xmuv);
    xph3 = xph3 * xfd * xbeta2 * 0.375;

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
    pl[6] = xmus * xmus + xmuv * xmuv;
    pl[7] = xlntau * pl[6];
    pl[8] = xmus * xmus * xmuv * xmuv;
    pl[9] = xlntau * pl[8];

    fs0 = 0.0;
    for (i = 0; i < 10; i++)
        fs0 += pl[i] * as0[i];
    fs1 = pl[0] * as1[0] + pl[1] * as1[1];
    fs2 = pl[0] * as2[0] + pl[1] * as2[1];
    xitot1 = xp1 + cfonc1 * fs0 * xmus;
    xitot2 = xp2 + cfonc2 * fs1 * xmus;
    xitot3 = xp3 + cfonc3 * fs2 * xmus;

    *xrray = xitot1 * xcosf1;
    *xrray += xitot2 * xcosf2 * 2.0;
    *xrray += xitot3 * xcosf3 * 2.0;
    *xrray /= xmus;
}


/******************************************************************************
MODULE:  comptg

PURPOSE:  Computes the transmission of the water vapor, ozone, and other gases.

RETURN VALUE:
Type = N/A

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/27/2014    Gail Schmidt     Conversion of the original FORTRAN code delivered
                              by Eric Vermote, NASA GSFC

NOTES:
    1. Standard sea level pressure is 1013 millibars.
******************************************************************************/
void comptg
(
    int iband,                   /* I: band index (0-based) */
    float xts,                   /* I: solar zenith angle */
    float xtv,                   /* I: observation zenith angle */
    float uoz,                   /* I: total column ozone */
    float uwv,                   /* I: total column water vapor (precipital
                                       water vapor) */
    float pres,                  /* I: surface pressure */
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
    float a, b;
    float a1, b0, b1;
    float m;
    float x;
    float p;                /* pressure in atmosphere */

    /* Compute ozone transmission */
    m = 1.0 / cos(xts * DEG2RAD) + 1.0 / cos(xtv * DEG2RAD);
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
    x = m * uwv * 0.5;
    if (x > 1.0E-06)
        *tgwvhalf = exp(-a * exp(log(x) * b));
    else
        *tgwvhalf = 1.0;

    /* Compute other gases transmission */
    a1 = ogtransa1[iband];
    b0 = ogtransb0[iband];
    b1 = ogtransb1[iband];
    p = pres * ONE_DIV_1013;   /* vs / 1013.0 */
    *tgog = -(a1 * p) * pow(m, exp(-(b0 + b1 * p)));
    *tgog = exp(*tgog);
}


/******************************************************************************
MODULE:  compsalb

PURPOSE:  Computes spherical albedo

RETURN VALUE:
Type = N/A

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/27/2014    Gail Schmidt     Conversion of the original FORTRAN code delivered
                              by Eric Vermote, NASA GSFC
8/14/2014    Gail Schmidt     Updated for v1.3 delivered by Eric Vermote

NOTES:
******************************************************************************/
void compsalb
(
    float raot550nm,                 /* I: nearest value of AOT */
    int iband,                       /* I: band index (0-based) */
    float pres,                      /* I: surface pressure */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ***sphalbt,                /* I: spherical albedo table
                                           [NSR_BANDS][7][22] */
    float ***normext,                /* I: aerosol extinction coefficient at
                                           the current wavelength (normalized
                                           at 550nm) [NSR_BANDS][7][22] */
    float *satm,                     /* O: spherical albedo */
    float *next                      /* O: ????? */
)
{
    int ip1, ip2, ip;               /* index variables for the pressure,
                                       AOT, and spherical albedo arrays */
    int iaot;                       /* aerosol optical thickness (AOT) index */
    int iaot1, iaot2;               /* index variables for the AOT and
                                       spherical albedo arrays */
    float xtiaot1, xtiaot2;         /* spherical albedo trans value */
    float satm1, satm2;             /* spherical albedo value */
    float next1, next2;             /* ???? */
    float dpres;                    /* pressure ratio */
    float deltaaot;                 /* AOT ratio */

    ip1 = 0;
    for (ip = 0; ip < 7; ip++)
    {
        if (pres < tpres[ip])
            ip1 = ip;
    }
    if (ip1 == 6)
        ip1 = 5;
    ip2 = ip1 + 1;
      
    iaot1 = 0;
    for (iaot = 0; iaot < 22; iaot++)
    {
        if (raot550nm > aot550nm[iaot])
            iaot1 = iaot;
    }
    if (iaot1 == 21)
        iaot1 = 20;
    iaot2 = iaot1 + 1;

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
Type = int
Value          Description
-----          -----------
ERROR          Error occurred computing transmission
SUCCESS        Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/25/2014    Gail Schmidt     Conversion of the original FORTRAN code delivered
                              by Eric Vermote, NASA GSFC

NOTES:
******************************************************************************/
int comptrans
(
    float xts,                       /* I: solar zenith */
    float raot550nm,                 /* I: nearest value of AOT */
    int iband,                       /* I: band index (0-based) */
    float pres,                      /* I: surface pressure */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****transt,                /* I: transmission table
                                           [NSR_BANDS][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float tts[22],                   /* I: sun angle table */
    float *xtts                      /* O: downward transmittance */
)
{
    char FUNC_NAME[] = "comptrans";   /* function name */
    char errmsg[STR_SIZE];            /* error message */
    int its;                        /* index for the sun angle table */
    int ip1, ip2, ip;               /* index variables for the pressure,
                                       AOT, and spherical albedo arrays */
    int iaot;                       /* aerosol optical thickness (AOT) index */
    int iaot1, iaot2;               /* index variables for the AOT and
                                       spherical albedo arrays */
    float xtiaot1, xtiaot2;         /* spherical albedo trans value */
    float xtts1, xtts2;
    float xmts, xtranst;
    float dpres;                    /* pressure ratio */
    float deltaaot;                 /* AOT ratio */

    ip1 = 0;
    for (ip = 0; ip < 7; ip++)
    {
        if (pres < tpres[ip])
            ip1 = ip;
    }
    if (ip1 == 6)
        ip1 = 5;
    ip2 = ip1 + 1;
      
    iaot1 = 0;
    for (iaot = 0; iaot < 22; iaot++)
    {
        if (raot550nm > aot550nm[iaot])
            iaot1 = iaot;
    }
    if (iaot1 == 21)
        iaot1 = 20;
    iaot2 = iaot1 + 1;

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

    xmts = (xts - tts[its]) * 0.25;    /* vs / 4.0 */
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

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  comproatm

PURPOSE:  Computes the atmospheric reflectance

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred computing the atmospheric reflectance
SUCCESS        Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/25/2014    Gail Schmidt     Conversion of the original FORTRAN code delivered
                              by Eric Vermote, NASA GSFC
7/1/2014     Gail Schmidt     Passed in xtsstep, xtsmin, xtvstep, and xtvmin
                              since they are defined in the main routine and
                              then duplicated here.

NOTES:
  1. GAIL -- this is likely the major time hog in the overall interpolation
******************************************************************************/
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
    float ****rolutt,                /* I: intrinsic reflectance table
                                           [NSR_BANDS][7][22][8000] */
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
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float *roatm                     /* O: atmospheric reflectance */
)
{
    char FUNC_NAME[] = "comproatm";   /* function name */
    char errmsg[STR_SIZE];            /* error message */
    int ip1, ip2, ip;               /* index variables for the pressure,
                                       AOT, and spherical albedo arrays */
    int iaot1, iaot2;               /* index variables for the AOT and
                                       spherical albedo arrays */
    int its;                        /* index for the sun angle table */
    int itv;                        /* index for the view angle table */
    int iaot;                       /* aerosol optical thickness (AOT) index */
    int isca;
    int iindex;
    float nbfic1, nbfic2, nbfic3, nbfic4;
    float nbfi1, nbfi2, nbfi3, nbfi4;
    float ro, rop1, rop2;           /* reflectance at p1 and p2 */
    float xmus;    /* cosine of solar zenith angle */
    float xmuv;    /* cosine of observation zenith angle */
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


    ip1 = 0;
    for (ip = 0; ip < 7; ip++)
    {
        if (pres < tpres[ip])
            ip1 = ip;
    }
    if (ip1 == 6)
        ip1 = 5;
    ip2 = ip1 + 1;

    if (xtv <= xtvmin)
        itv = 0;
    else
        itv = (int) ((xtv - xtvmin) / xtvstep + 1.0);

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

    xmuv = cos(xtv * DEG2RAD);
    xmus = cos(xts * DEG2RAD);
    cscaa = -xmus * xmuv - cos(xfi * DEG2RAD) * sqrt(1.0 - xmus * xmus) *
        sqrt(1.0 - xmuv * xmuv);
    scaa = acos(cscaa) * RAD2DEG;    /* vs / DEG2RAD */

    /* Determine lower and uper limit in aot */
    iaot1 = 0;
    for (iaot = 0; iaot < 22; iaot++)
    {
        if (raot550nm > aot550nm[iaot])
            iaot1 = iaot;
    }
    if (iaot1 == 21)
        iaot1 = 20;
    iaot2 = iaot1 + 1;

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
        if ((isca + 1) < nbfi[itv][its])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv][its] - 1;
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
        if ((isca + 1) < nbfi[itv][its+1])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv][its+1] - 1;
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
        if ((isca + 1) < nbfi[itv+1][its])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv+1][its] - 1;
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
    if ((isca + 1) < nbfi[itv+1][its+1])
    {
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = xtsmax - isca * 4.0;
    }
    else
    {
        isca = nbfi[itv+1][its+1] - 1;
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = tsmin[itv+1][its+1];
    }

    iindex = indts[its+1] + nbfic4 - nbfi4 + isca - 1;
    roinf = rolutt[iband][ip1][iaot1][iindex];
    rosup = rolutt[iband][ip1][iaot1][iindex+1];
    ro4 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);

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
        if ((isca + 1) < nbfi[itv][its])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv][its] - 1;
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
        if ((isca + 1) < nbfi[itv][its+1])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv][its+1] - 1;
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
        if ((isca + 1) < nbfi[itv+1][its])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv+1][its] - 1;
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
    if ((isca + 1) < nbfi[itv+1][its+1])
    {
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = xtsmax - isca * 4.0;
    }
    else
    {
        isca = nbfi[itv+1][its+1] - 1;
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = tsmin[itv+1][its+1];
    }

    iindex = indts[its+1] + nbfic4 - nbfi4 + isca - 1;
    roinf = rolutt[iband][ip1][iaot2][iindex];
    rosup = rolutt[iband][ip1][iaot2][iindex+1];
    ro4 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);

    t = (tts[its+1] - xts) / (tts[its+1] - tts[its]);
    u = (ttv[itv+1][its] - xtv) / (ttv[itv+1][its] - ttv[itv][its]);
    roiaot2 = ro1 * t * u + ro2 * u * (1.0 - t) + ro3 * (1.0 - u) * t +
        ro4 * (1.0 - u) * (1.0 - t);

    /* Interpolation as log of tau */
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
        if ((isca + 1) < nbfi[itv][its])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv][its] - 1;
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
        if ((isca + 1) < nbfi[itv][its+1])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv][its+1] - 1;
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
        if ((isca + 1) < nbfi[itv+1][its])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv+1][its] - 1;
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
    if ((isca + 1) < nbfi[itv+1][its+1])
    {
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = xtsmax - isca * 4.0;
    }
    else
    {
        isca = nbfi[itv+1][its+1] - 1;
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = tsmin[itv+1][its+1];
    }

    iindex = indts[its+1] + nbfic4 - nbfi4 + isca - 1;
    roinf = rolutt[iband][ip2][iaot1][iindex];
    rosup = rolutt[iband][ip2][iaot1][iindex+1];
    ro4 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);

    t = (tts[its+1] - xts) / (tts[its+1] - tts[its]);
    u = (ttv[itv+1][its] - xtv) / (ttv[itv+1][its] - ttv[itv][its]);
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
        if ((isca + 1) < nbfi[itv][its])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv][its] - 1;
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
        if ((isca + 1) < nbfi[itv][its+1])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv][its+1] - 1;
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
        if ((isca + 1) < nbfi[itv+1][its])
        {
            sca1 = xtsmax - (isca - 1) * 4.0;
            sca2 = xtsmax - isca * 4.0;
        }
        else
        {
 	        isca = nbfi[itv+1][its] - 1;
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
    if ((isca + 1) < nbfi[itv+1][its+1])
    {
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = xtsmax - isca * 4.0;
    }
    else
    {
        isca = nbfi[itv+1][its+1] - 1;
        sca1 = xtsmax - (isca - 1) * 4.0;
        sca2 = tsmin[itv+1][its+1];
    }

    iindex = indts[its+1] + nbfic4 - nbfi4 + isca - 1;
    roinf = rolutt[iband][ip2][iaot2][iindex];
    rosup = rolutt[iband][ip2][iaot2][iindex+1];
    ro4 = roinf + (rosup - roinf) * (scaa - sca1) / (sca2 - sca1);

    t = (tts[its+1] - xts) / (tts[its+1] - tts[its]);
    u = (ttv[itv+1][its] - xtv) / (ttv[itv+1][its] - ttv[itv][its]);
    roiaot2 = ro1 * t * u + ro2 * u * (1.0 - t) + ro3 * (1.0 - u) * t +
        ro4 * (1.0 - u) * (1.0 - t);

    /* Interpolation as log of tau */
    deltaaot = logaot550nm[iaot2] - logaot550nm[iaot1];
    deltaaot = (log (raot550nm) - logaot550nm[iaot1]) / deltaaot;
    ro = roiaot1 + (roiaot2 - roiaot1) * deltaaot;
    rop2 = ro;

    dpres = (pres - tpres[ip1]) / (tpres[ip2] - tpres[ip1]);
    *roatm = rop1 + (rop2 - rop1) * dpres;
       
    /* Successful completion */
    return (SUCCESS);
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

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/27/2014    Gail Schmidt     Conversion of the original FORTRAN code delivered
                              by Eric Vermote, NASA GSFC
7/1/2014     Gail Schmidt     Passed in xtsstep and xtsmin since they are
                              defined in the main routine and then duplicated
                              here.
7/22/2014    Gail Schmidt     Removed the read of the tauray file, since it's
                              a small, static ASCII file.  Just made it a
                              hard-coded array in the main function.  Ditto for
                              the gaseous transmission coefficient file.
8/14/2014    Gail Schmidt     Updated for v1.3 delivered by Eric Vermote

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

    for (i = 0; i < 20; i++)
    {
        start[0] = i;   /* lines */
        edges[0] = 1;   /* number of lines */
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

    for (i = 0; i < 20; i++)
    {
        start[0] = i;   /* lines */
        edges[0] = 1;   /* number of lines */
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

    for (i = 0; i < 20; i++)
    {
        start[0] = i;   /* lines */
        edges[0] = 1;   /* number of lines */
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

    for (i = 0; i < 20; i++)
    {
        start[0] = i;   /* lines */
        edges[0] = 1;   /* number of lines */
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

    for (i = 0; i < 20; i++)
    {
        start[0] = i;   /* lines */
        edges[0] = 1;   /* number of lines */
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

