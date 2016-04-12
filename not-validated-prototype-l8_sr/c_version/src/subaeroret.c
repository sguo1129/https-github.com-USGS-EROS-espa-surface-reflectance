/*****************************************************************************
FILE: subaeroret.c
  
PURPOSE: Contains functions for handling the aerosol inversion as well as the
aerosol interpolations.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/
#include "lut_subr.h"
#include "l8_sr.h"

/******************************************************************************
MODULE:  subaeroret

PURPOSE:  This subroutine reads the lookup table (LUT) and performs the aerosol
retrievals, using band ratios.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the LUT or doing the correction
SUCCESS        Successful completion

NOTES:
******************************************************************************/
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
)
{
    char FUNC_NAME[] = "subaeroret";   /* function name */
    char errmsg[STR_SIZE];       /* error message */
    int nit;                /* number of iterations */
    int iter;               /* looping variable for iterations */
    int iaot;               /* aerosol optical thickness (AOT) index */
    int retval;             /* function return value */
    bool flagn;             /* flag to start AOT convergence */
    float raot550nm=0.0;    /* nearest input value of AOT */
    float roslamb;          /* lambertian surface reflectance */
    double ros1, ros3;      /* surface reflectance for bands */
    double raot1, raot2;    /* AOT ratios that bracket the predicted ratio */
    float next;             /* ???? */
    float tgo;              /* other gaseous transmittance */
    float roatm;            /* atmospheric intrinsic reflectance */
    float ttatmg;           /* total atmospheric transmission */
    float satm;             /* spherical albedo */
    float xrorayp;          /* reflectance of the atmosphere due to molecular
                               (Rayleigh) scattering */
    double aratio1, aratio2;
    double pratio;          /* targeted ratio between the surface reflectance
                               in two bands */
    double eratio;
    double eaot;            /* estimate of AOT */
    double th1, th3;
    double peratio;
    double pros1, pros3;    /* predicted surface reflectance */

    /* Correct band 3 and band 1 with increasing AOT (using pre till ratio is
       equal to erelc[2]). pratio is the targeted ratio between the surface
       reflectance in band 4 (ros4) and band 2 (ros2). */
    iaot = 0;
    pratio = erelc[iband3] / erelc[iband1];
    aratio1 = 1000.0;
    aratio2 = 2000.0;
    ros1 = 1.0;
    ros3 = 1.0;
    raot1 = 0.0001;
    raot2 = 0.0;
    flagn = false;
    th1 = 0.01;
    th3 = 0.01;
    pros1 = 0.0;
    pros3 = 0.0;

    /* The ratio decreases as the AOT increases.  The exit conditions in this
       loop are when two values of AOT can be found that bracket the predicted
       ratio (pratio). */
    nit = 0;
    while ((iaot < 22) && (aratio1 > pratio) && (ros1 > th1) && (ros3 > th3) &&
        ((aratio1 - 0.01) < aratio2) && (nit < 30))
    {
        ros1 = -1.0;
        ros3 = -1.0;

        /* If flagn is set... start converge to the AOT bounds by dichotomy to
           increase the accuracy of the retrieval */
        if (!flagn)
            raot550nm = aot550nm[iaot];
        else
            raot550nm = (raot1 + aot550nm[iaot]) * 0.5;

        /* Loop until convergence.  Add a mechanism to stop the loop from being
           infinite by stopping at 50 iterations. */
        iter = 0;
        while ((ros1 < th1 || ros3 < th3) && (iter < 50))
        {
            if (iter > 0)
            {
                if (iaot >= 1)
                {
                    raot550nm = (raot550nm + aot550nm[iaot-1]) * 0.5;
                }
                else
                {
                    /* Inversion failed.  Compute the model residual with
                       what we have and then return */
                    *raot = raot550nm;
                    retval = subaeroret_residual (iband1, iband3, ros1, ros3,
                        roslamb, pratio, raot550nm, xts, xtv, xmus, xmuv, xfi,
                        cosxfi, pres, uoz, uwv, erelc, troatm, tpres, aot550nm,
                        rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
                        sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts,
                        indts, ttv, tauray, ogtransa1, ogtransb0, ogtransb1,
                        wvtransa, wvtransb, oztransa, residual, snext);
                    if (retval != SUCCESS)
                    {
                        sprintf (errmsg, "Computing the subaeroret model "
                            "residual");
                        error_handler (true, FUNC_NAME, errmsg);
                        return (ERROR);
                    }

                    return (SUCCESS);
                }
            }

            /* Atmospheric correction for band 3 */
            retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm,
                iband3, pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin,
                xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi,
                tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                ogtransb1, wvtransa, wvtransb, oztransa, troatm[iband3],
                &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Performing lambertian atmospheric correction "
                    "type 2.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            ros3 = roslamb;

            /* Atmospheric correction for band 1 */
            retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm,
                iband1, pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin,
                xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi,
                tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                ogtransb1, wvtransa, wvtransb, oztransa, troatm[iband1],
                &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Performing lambertian atmospheric correction "
                    "type 2.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            ros1 = roslamb;

            /* Keep count of the iterations */
            iter++;
        } /* end while */

        if ((iter > 1) || flagn)
            flagn = true;
        else
            iaot++;

        if ((ros1 > th1) && (ros3 > th3))
        {
            aratio2 = aratio1;
            raot2 = raot1;
            raot1 = raot550nm;
            aratio1 = ros3 / ros1;
        }

        /* Increment the number of iterations */
        nit++;
    }  /* end while */

    /* Once the two values of AOT (raot2 and raot1) that gives ratios that
       bracket the predicted ratio are found, they are used to estimate the
       AOT (eaot) using linear interpolation. */
    if ((aratio1 > pratio) && (aratio2 > pratio))
    {
        /* Early break out if the ratios are not valid */
        if (raot1 < raot2)
            raot550nm = raot1;
        else
            raot550nm = raot2;
        *raot = raot550nm;

        retval = subaeroret_residual (iband1, iband3, ros1, ros3, roslamb,
            pratio, raot550nm, xts, xtv, xmus, xmuv, xfi, cosxfi, pres, uoz,
            uwv, erelc, troatm, tpres, aot550nm, rolutt, transt, xtsstep,
            xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic,
            nbfi, tts, indts, ttv, tauray, ogtransa1, ogtransb0, ogtransb1,
            wvtransa, wvtransb, oztransa, residual, snext);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Computing the subaeroret model residual");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        return (SUCCESS);
    }

    /* Compute the estimated AOT */
    eaot = (aratio1 - pratio) * (raot2 - raot1) / (aratio1 - aratio2) + raot1;

    /* The estimated AOT is refined by performing an additional iteration of
       atmospheric correction */
    /* Atmospheric correction for band 3 */
    raot550nm = eaot;
    retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm, iband3,
        pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
        sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
        uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
        oztransa, troatm[iband3], &roslamb, &tgo, &roatm, &ttatmg, &satm,
        &xrorayp, &next);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Performing lambertian atmospheric correction "
            "type 2.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    ros3 = roslamb;

    /* Atmospheric correction for band 1 */
    retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm, iband1,
        pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
        sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
        uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
        oztransa, troatm[iband1], &roslamb, &tgo, &roatm, &ttatmg, &satm,
        &xrorayp, &next);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Performing lambertian atmospheric correction "
            "type 2.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    ros1 = roslamb;

    /* Compute the estimated AOT, depending on which ratio is appropriate */
    eratio = ros3 / ros1;
    if (fabs (eratio - aratio1) > fabs (eratio - aratio2))
    {
        raot2 = eaot;
        aratio2 = eratio;
        eaot = (aratio1 - pratio) * (raot2 - raot1)/(aratio1 - aratio2) + raot1;
    }
    else
    {
        raot1 = eaot;
        aratio1 = eratio;
        eaot = (aratio1 - pratio) * (raot2 - raot1)/(aratio1 - aratio2) + raot1;
    }

    /* Atmospheric correction for band 3 */
    raot550nm = eaot;
    retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm, iband3,
        pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
        sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
        uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
        oztransa, troatm[iband3], &roslamb, &tgo, &roatm, &ttatmg, &satm,
        &xrorayp, &next);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Performing lambertian atmospheric correction "
            "type 2.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    ros3 = roslamb;

    /* Atmospheric correction for band 1 */
    retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm, iband1,
        pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
        sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
        uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
        oztransa, troatm[iband1], &roslamb, &tgo, &roatm, &ttatmg, &satm,
        &xrorayp, &next);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Performing lambertian atmospheric correction "
            "type 2.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    ros1 = roslamb;

    /* The last step of the algorithm is done by making very small increases or
       decreases of the estimated AOT.  The value of the estimated AOT which
       provides the closest ratio to the predicted ratio is finally selected. */
    eratio = ros3 / ros1;
    raot550nm = eaot;
    if (raot550nm >= 0.01)
    {
        peratio = 1000.0;
        if (eratio > pratio)
        {
            while ((eratio > pratio) && (peratio > eratio))
            {  /* Increase the raot550nm */
                /* Atmospheric correction for band 3 */
                pros1 = ros1;
                pros3 = ros3;
                raot550nm += 0.005;
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, iband3, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, troatm[iband3], &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, &next);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
                ros3 = roslamb;

                /* Atmospheric correction for band 1 */
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, iband1, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, troatm[iband1], &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, &next);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
                ros1 = roslamb;
                peratio = eratio;
                eratio = ros3 / ros1;
            }

            if (fabs (eratio - pratio) > fabs (peratio - pratio))
            {
                raot550nm -= 0.005;
                eratio = peratio;
                ros1 = pros1;
                ros3 = pros3;
            }
        }
        else
        {
            peratio = 0.0;
            while ((eratio < pratio) && (peratio < eratio))
            {  /* Decrease the raot550nm */
                /* Atmospheric correction for band 3 */
                pros1 = ros1;
                pros3 = ros3;
                raot550nm -= 0.005;
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, iband3, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, troatm[iband3], &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, &next);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
                ros3 = roslamb;

                /* Atmospheric correction for band 1 */
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, iband1, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, troatm[iband1], &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, &next);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
                ros1 = roslamb;
                peratio = eratio;
                eratio = ros3 / ros1;
            }

            if (fabs (eratio - pratio) > fabs (peratio - pratio))
            {
                raot550nm += 0.005;
                eratio = peratio;
                ros1 = pros1;
                ros3 = pros3;
            }
        }  /* end else */
    }  /* if raot550nm */
    *raot = raot550nm;

    /* Compute the model residual */
    retval = subaeroret_residual (iband1, iband3, ros1, ros3, roslamb,
        pratio, raot550nm, xts, xtv, xmus, xmuv, xfi, cosxfi, pres, uoz, uwv,
        erelc, troatm, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin,
        xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts,
        indts, ttv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
        wvtransb, oztransa, residual, snext);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Computing the subaeroret model residual");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  subaeroret_residual

PURPOSE:  Computes the model residual for the aerosol inversion.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred computing the residual
SUCCESS        Successful completion

NOTES:
******************************************************************************/
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
)
{
    char FUNC_NAME[] = "subaeroret_residual";   /* function name */
    char errmsg[STR_SIZE];  /* error message */
    int iband;              /* looping variable for bands */
    int retval;             /* function return value */
    float next;             /* ???? */
    float tgo;
    float roatm;            /* atmospherice intrinsic reflectance */
    float ttatmg;
    float satm;             /* spherical albedo */
    float xrorayp;          /* molecular reflectance */

    /* Compute the model residual
       Note - Eric indicated the residual on Band 7 was not to be used, so
       stop the residual calculations with Band 6. */
    *residual = fabs (ros3 - ros1 * pratio);
    for (iband = 0; iband <= DN_BAND6; iband++)
    {
        if (erelc[iband] > 0.0)
        {
            retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm,
                iband, pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin,
                xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi,
                tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                ogtransb1, wvtransa, wvtransb, oztransa, troatm[iband],
                &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Performing lambertian atmospheric correction "
                    "type 2.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            *residual += fabs (roslamb - ros1 * (erelc[iband] / erelc[iband1]));
            if (iband == iband3)
                *snext = next;
        }
    }

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  aerosol_interpolation

PURPOSE:  Performs interpolation for the zero-valued aerosol pixels.

RETURN VALUE: N/A

NOTES:
1. Does not use the water, cloud, or cirrus pixels as part of the interpolation.
2. Aerosol values are modified inline in the taero input buffer.
3. This tends to produce blocks in the interpolated aerosols, resulting in a
   blocky output product in areas of water, etc.
******************************************************************************/
void aerosol_interpolation
(
    int nlines,       /* I: number of lines in the aerosol buffer */
    int nsamps,       /* I: number of samples in the aerosol buffer */
    uint8 *cloud,     /* I: bit-packed value that represent clouds,
                            nlines x nsamps */
    float *tresi,     /* I: residuals for each pixel, nlines x nsamps;
                            tresi < 0.0 flags water pixels and pixels with
                            high residuals */
    float *taero      /* I/O: aerosol values for each pixel, nlines x nsamps */
)
{
    int i, j, k, l;    /* looping variable for pixels */
    int win_pix;       /* current pixel in the line,sample window */

    int nbaot;         /* number of AOT pixels (non-cloud/water) for aerosol
                          interpolation */
    int step;          /* step value for aerosol interpolation */
    bool hole;         /* is this a hole in the aerosol retrieval area? */
    double aaot;       /* average of AOT */
    double sresi;      /* sum of 1 / residuals */

    /* Aerosol interpolation. Does not use water, cloud, or cirrus pixels. */
    printf ("Performing aerosol interpolation ...\n");
    hole = true;
    step = 10;
    while (hole && (step < 1000))
    {
        hole = false;
        for (i = 0; i < nlines; i += step)
        {
            for (j = 0; j < nsamps; j += step)
            {
                nbaot = 0;
                aaot = 0.0;
                sresi = 0.0;

                /* Check the step x step window around the current pixel */
                for (k = i - step*0.5; k <= i + step*0.5; k++)
                {
                    /* Make sure the line is valid */
                    if (k < 0 || k >= nlines)
                        continue;

                    win_pix = k * nsamps + j - step*0.5;
                    for (l = j - step*0.5; l <= j + step*0.5; l++, win_pix++)
                    {
                        /* Make sure the sample is valid */
                        if (l < 0 || l >= nsamps)
                            continue;

                        /* Check for clear pixels with positive residuals */
                        if ((tresi[win_pix] > 0) && (cloud[win_pix] == 0))
                        {
                            nbaot++;
                            aaot += taero[win_pix] / tresi[win_pix];
                            sresi += 1.0 / tresi[win_pix];
                        }
                    }
                }

                /* If pixels were found */
                if (nbaot != 0)
                {
                    aaot /= sresi;

                    /* Check the step x step window around the current pixel */
                    for (k = i - step*0.5; k <= i + step*0.5; k++)
                    {
                        /* Make sure the line is valid */
                        if (k < 0 || k >= nlines)
                            continue;

                        win_pix = k * nsamps + j - step*0.5;
                        for (l = j - step*0.5; l <= j + step*0.5;
                             l++, win_pix++)
                        {
                            /* Make sure the sample is valid */
                            if (l < 0 || l >= nsamps)
                                continue;

                            if ((tresi[win_pix] < 0) &&
                                (!btest (cloud[win_pix], CIR_QA)) &&
                                (!btest (cloud[win_pix], CLD_QA)) &&
                                (!btest (cloud[win_pix], WAT_QA)))
                            {
                                taero[win_pix] = aaot;
                                tresi[win_pix] = 1.0;
                            }
                        }  /* for l */
                    }  /* for k */
                }
                else
                {  /* this is a hole */
                    hole = true;
                }
            }  /* end for j */
        }  /* end for i */

        /* Modify the step value */
        step *= 2;
    }  /* end while */
}


#ifdef SINGLE
/******************************************************************************
MODULE:  aerosol_interpolation2 (single-threaded)

PURPOSE:  Performs interpolation for the zero-valued aerosol pixels using a
moving window approach, based on the original interpolation code.

RETURN VALUE: N/A

NOTES:
1. Does not use the water, cloud, or cirrus pixels as part of the interpolation.
2. A local copy of the aerosols is made and used so that the interpolated
   aerosol pixels are not used to interpolate other aerosol pixels.  Only the
   original aerosol values are used for the interpolation.
******************************************************************************/
void aerosol_interpolation2
(
    int nlines,       /* I: number of lines in the aerosol buffer */
    int nsamps,       /* I: number of samples in the aerosol buffer */
    uint8 *cloud,     /* I: bit-packed value that represent clouds,
                            nlines x nsamps */
    float *tresi,     /* I: residuals for each pixel, nlines x nsamps;
                            tresi < 0.0 flags water pixels and pixels with
                            high residuals */
    float *taero      /* I/O: aerosol values for each pixel, nlines x nsamps */
)
{
    int i, j, k, l;    /* looping variable for pixels */
    int win_pix;       /* current pixel in the line,sample window */
    int curr_pix;      /* current pixel in the scene */
    int nbaot;         /* number of AOT pixels (non-cloud/water) for aerosol
                          interpolation */
    int step;          /* step value for aerosol interpolation */
    bool hole;         /* is this a hole in the aerosol retrieval area? */
    double aaot;       /* average of AOT */
    double sresi;      /* sum of 1 / residuals */

    /* Aerosol interpolation. Does not use water, cloud, or cirrus pixels. */
    printf ("Performing aerosol moving window interpolation ...\n");
    for (i = 0; i < nlines; i++)
    {
        curr_pix = i * nsamps;
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            /* If the current pixel is a hole in the aerosols and it's not
               cirrus, cloud, or water, then look at filling it */
            if ((tresi[curr_pix] < 0) &&
                (!btest (cloud[curr_pix], CIR_QA)) &&
                (!btest (cloud[curr_pix], CLD_QA)) &&
                (!btest (cloud[curr_pix], WAT_QA)))
            {
                /* Look at the step x step window and use it to fill the pixel
                   if there are enough non-fill aerosols */
                hole = true;
                step = 10;
                while (hole && (step < 1000))
                {
                    /* Initialize our variables for this pixel */
                    hole = false;
                    nbaot = 0;
                    aaot = 0.0;
                    sresi = 0.0;

                    /* Check the step x step window around the current pixel */
                    for (k = i - step*0.5; k <= i + step*0.5; k++)
                    {
                        /* Make sure the line is valid */
                        if (k < 0 || k >= nlines)
                            continue;

                        win_pix = k * nsamps + j - step*0.5;
                        for (l = j - step*0.5; l <= j + step*0.5; l++,
                             win_pix++)
                        {
                            /* Make sure the sample is valid */
                            if (l < 0 || l >= nsamps)
                                continue;

                            /* Check for clear pixels with positive residuals */
                            if ((tresi[win_pix] > 0) && (cloud[win_pix] == 0))
                            {
                                nbaot++;
                                aaot += taero[win_pix] / tresi[win_pix];
                                sresi += 1.0 / tresi[win_pix];
                            }
                        }
                    }

                    /* If at least 30% of the window had valid AOT values */
                    if (nbaot > step * step * 0.3)
                    {
                        taero[curr_pix] = aaot / sresi;
                        tresi[curr_pix] = 1.0;
                    }
                    else
                    {
                        /* this is a hole, increment the step value */
                        hole = true;
                        step *= 2;
                    }
                }  /* end while */
            }  /* end if this is fill and not cloud, water, cirrus */
        }  /* end for j */
    }  /* end for i */

    /* Successful completion */
    return (SUCCESS);
}
#endif


/******************************************************************************
MODULE:  aerosol_interpolation_mw

PURPOSE:  Performs interpolation for the zero-valued aerosol pixels using a
Gaussian moving window approach.  Interpolated values are not used for filling
the remaining holes.  For each aerosol-fill value in the scene, aerosol values
are interpolated using an NxN window that increases by a value of 2 if the
initial NxN window does not fulfill the requirements of the algorithm.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred allocating memory or copying taero
SUCCESS        Successful completion

NOTES:
1. Does not use the water, cloud, or cirrus pixels as part of the interpolation.
2. A local copy of the aerosols is made and used so that the interpolated
   aerosol pixels are not used to interpolate other aerosol pixels.  Only the
   original aerosol values are used for the interpolation.
******************************************************************************/
int aerosol_interpolation_mw
(
    int nlines,       /* I: number of lines in the aerosol buffer */
    int nsamps,       /* I: number of samples in the aerosol buffer */
    uint8 *cloud,     /* I: bit-packed value that represent clouds,
                            nlines x nsamps */
    float *tresi,     /* I: residuals for each pixel, nlines x nsamps;
                            tresi < 0.0 flags water pixels and pixels with
                            high residuals */
    float *taero      /* I/O: aerosol values for each pixel, nlines x nsamps */
)
{
    char FUNC_NAME[] = "aerosol_interpolation_mw";   /* function name */
    char errmsg[STR_SIZE];  /* error message */
    int i, j, k, l;    /* looping variable for pixels */
    int win_pix;       /* current pixel in the line,sample window */
    int curr_pix;      /* current pixel in the scene */
    int nbaot;         /* number of AOT pixels (non-cloud/water) for aerosol
                          interpolation */
    int step;          /* step value for aerosol interpolation */
    bool hole;         /* is this a hole in the aerosol retrieval area? */
    double aaot;       /* average of AOT */
    double sresi;      /* sum of 1 / residuals */
    float *ltaero = NULL;  /* local copy of aerosols for filling */
    float *ltresi = NULL;  /* local copy of residuals for filling */

    /* Allocate space for local copy of the aerosols and residuals */
    ltaero = calloc (nlines*nsamps, sizeof (float));
    if (ltaero == NULL)
    {
        sprintf (errmsg, "Error allocating memory for a local copy of taero");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    ltresi = calloc (nlines*nsamps, sizeof (float));
    if (ltresi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for a local copy of tresi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Make a local copy of the aerosols and residuals which will remain a
       clean copy of the unfilled data to be used for interpolation */
    memcpy (ltaero, taero, nlines * nsamps * sizeof (float));
    if (ltaero == NULL)
    {
        sprintf (errmsg, "Error copying memory for a local copy of taero");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    memcpy (ltresi, tresi, nlines * nsamps * sizeof (float));
    if (ltresi == NULL)
    {
        sprintf (errmsg, "Error copying memory for a local copy of tresi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Aerosol interpolation. Does not use water, cloud, or cirrus pixels. */
    printf ("Performing multi-threaded aerosol moving window interpolation "
        "...\n");
#ifdef _OPENMP
    #pragma omp parallel for private (i, j, k, l, hole, step, nbaot, aaot, sresi, curr_pix, win_pix)
#endif
    for (i = 0; i < nlines; i++)
    {
        curr_pix = i * nsamps;
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            /* If the current pixel is a hole in the aerosols and it's not
               cirrus, cloud, or water, then look at filling it */
            if ((ltresi[curr_pix] < 0) &&
                (!btest (cloud[curr_pix], CIR_QA)) &&
                (!btest (cloud[curr_pix], CLD_QA)) &&
                (!btest (cloud[curr_pix], WAT_QA)))
            {
                /* Look at the step x step window and use it to fill the pixel
                   if there are enough non-fill aerosols */
                hole = true;
                step = 10;
                while (hole && (step < 1000))
                {
                    /* Initialize our variables for this pixel */
                    hole = false;
                    nbaot = 0;
                    aaot = 0.0;
                    sresi = 0.0;

                    /* Check the step x step window around the current pixel */
                    for (k = i - step*0.5; k <= i + step*0.5; k++)
                    {
                        /* Make sure the line is valid */
                        if (k < 0 || k >= nlines)
                            continue;

                        win_pix = k * nsamps + j - step*0.5;
                        for (l = j - step*0.5; l <= j + step*0.5; l++,
                             win_pix++)
                        {
                            /* Make sure the sample is valid */
                            if (l < 0 || l >= nsamps)
                                continue;

                            /* Check for clear pixels with positive residuals */
                            if ((ltresi[win_pix] > 0) && (cloud[win_pix] == 0))
                            {
                                nbaot++;
                                aaot += ltaero[win_pix] / ltresi[win_pix];
                                sresi += 1.0 / ltresi[win_pix];
                            }
                        }
                    }

                    /* If at least 30% of the window had valid AOT values */
                    if (nbaot > step * step * 0.3)
                    {
                        taero[curr_pix] = aaot / sresi;
                        tresi[curr_pix] = 1.0;
                    }
                    else
                    {
                        /* this is a hole, increment the step value */
                        hole = true;
                        step *= 2;
                    }
                }  /* end while */
            }  /* end if this is fill and not cloud, water, cirrus */
        }  /* end for j */
    }  /* end for i */

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  aerosol_interpolation_mw2

PURPOSE:  Performs interpolation for the zero-valued aerosol pixels using a
moving window approach.  Start with a smaller moving window, then grow the
window in the next iteration.  Interpolated pixels are used in-between
iterations through the scene.  Thus the first iteration only uses the original
aerosol values.  The next iteration uses the interpolated values from the
previous iteration. Etc.  Iterations through the scene are made for a constant
NxN window, then the window is increased for the next iteration.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred allocating memory or copying taero
SUCCESS        Successful completion

NOTES:
1. Does not use the water, cloud, or cirrus pixels as part of the interpolation.
2. A local copy of the aerosols is made and used so that the interpolated
   aerosol pixels are not used to interpolate other aerosol pixels, within each
   sweep of the scene with the current moving window.
******************************************************************************/
int aerosol_interpolation_mw2
(
    int nlines,       /* I: number of lines in the aerosol buffer */
    int nsamps,       /* I: number of samples in the aerosol buffer */
    uint8 *cloud,     /* I: bit-packed value that represent clouds,
                            nlines x nsamps */
    float *tresi,     /* I: residuals for each pixel, nlines x nsamps;
                            tresi < 0.0 flags water pixels and pixels with
                            high residuals */
    float *taero      /* I/O: aerosol values for each pixel, nlines x nsamps */
)
{
    char FUNC_NAME[] = "aerosol_interpolation_mw2";   /* function name */
    char errmsg[STR_SIZE];  /* error message */
    int i, j, k, l;    /* looping variable for pixels */
    int win_pix;       /* current pixel in the line,sample window */
    int curr_pix;      /* current pixel in the scene */
    int nbaot;         /* number of AOT pixels (non-cloud/water) for aerosol
                          interpolation */
    int step;          /* step value for aerosol interpolation */
    int half_step;     /* half the step value for windowing */
    int min_valid_nbaot;  /* minimum number of AOT pixels for this window to
                             be a valid fill */
    bool hole;         /* is this a hole in the aerosol retrieval area? */
    double aaot;       /* average of AOT */
    double s1resi;     /* sum of 1 / residuals */
    double sresi;      /* sum of residuals */
    float valid_percent;   /* percentage of window pixels needed to be non-fill
                              in order to use this window for interpolation */
    float *ltaero = NULL;  /* local copy of aerosols for filling */
    float *ltresi = NULL;  /* local copy of residuals for filling */

    /* Allocate space for local copy of the aerosols and residuals */
    ltaero = calloc (nlines*nsamps, sizeof (float));
    if (ltaero == NULL)
    {
        sprintf (errmsg, "Error allocating memory for a local copy of taero");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    ltresi = calloc (nlines*nsamps, sizeof (float));
    if (ltresi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for a local copy of tresi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Make a local copy of the aerosols and residuals which will remain a
       clean copy of the unfilled data to be used for interpolation.  It will
       get updated with each sweep of the moving window through the scene to
       be used in the next sweep. */
    memcpy (ltaero, taero, nlines * nsamps * sizeof (float));
    if (ltaero == NULL)
    {
        sprintf (errmsg, "Error copying memory for a local copy of taero");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    memcpy (ltresi, tresi, nlines * nsamps * sizeof (float));
    if (ltresi == NULL)
    {
        sprintf (errmsg, "Error copying memory for a local copy of tresi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Aerosol interpolation. Does not use water, cloud, or cirrus pixels. */
    printf ("Performing multi-threaded aerosol moving window interpolation "
        "with growing window sizes ...\n");

    /* Look at the step x step window and use it to fill the pixel if there
       are enough non-fill aerosols.  Increase the size of the window the next
       time through the loop until all pixels are filled. */
    hole = true;
    step = 5;
    valid_percent = 0.5;  /* require 50% of non-fill pixels to be valid */
    while (hole && step < 1000)
    {
        printf ("   *Sweep with window size of %d x %d\n", step, step);
        /* Initialize our variables for this iteration */
        hole = false;
        half_step = step * 0.5;
        min_valid_nbaot = step * step * valid_percent;
#ifdef _OPENMP
    #pragma omp parallel for private (i, j, k, l, nbaot, aaot, s1resi, sresi, curr_pix, win_pix)
#endif
        for (i = 0; i < nlines; i++)
        {
            curr_pix = i * nsamps;
            for (j = 0; j < nsamps; j++, curr_pix++)
            {
                /* If the current pixel is a hole in the aerosols and it's not
                   cirrus, cloud, or water, then look at filling it */
                if ((ltresi[curr_pix] < 0) &&
                    (!btest (cloud[curr_pix], CIR_QA)) &&
                    (!btest (cloud[curr_pix], CLD_QA)) &&
                    (!btest (cloud[curr_pix], WAT_QA)))
                {
                    /* Look at the step x step window and use it to fill the
                       pixel if there are enough non-fill aerosols */
                    nbaot = 0;
                    aaot = 0.0;
                    s1resi = 0.0;
                    sresi = 0.0;

                    /* Check the step x step window around the current pixel */
                    for (k = i - half_step; k <= i + half_step; k++)
                    {
                        /* Make sure the line is valid */
                        if (k < 0 || k >= nlines)
                            continue;

                        win_pix = k * nsamps + j - half_step;
                        for (l = j - half_step; l <= j + half_step; l++,
                             win_pix++)
                        {
                            /* Make sure the sample is valid */
                            if (l < 0 || l >= nsamps)
                                continue;

                            /* Check for clear pixels with positive residuals.
                               This skips fill pixels as well, since they have
                               a tresi value of 0.0. */
                            if ((ltresi[win_pix] > 0) && (cloud[win_pix] == 0))
                            {
                                nbaot++;
                                aaot += ltaero[win_pix] / ltresi[win_pix];
                                s1resi += 1.0 / ltresi[win_pix];
                                sresi += ltresi[win_pix];
                            }
                        }
                    }

                    /* If the minimum number of valid AOT values for this
                       window is met then use the interpolation otherwise
                       wait to fill this pixel during the next sweep. */
                    if (nbaot > min_valid_nbaot)
                    {
                        taero[curr_pix] = aaot / s1resi;
                        tresi[curr_pix] = sresi / nbaot;
                    }
                    else
                    {
                        /* This is a hole in the current sweep.  Need to
                           protect this for updates across each line. */
                        /* TODO -- consider a local hole variable that is
                           private and then updated after all the samples
                           have been updated.  That should provide speedup. */
#ifdef _OPENMP
                        #pragma omp atomic
#endif
                        hole |= true;
                    }
                }  /* end if this is fill and not cloud, water, cirrus */
            }  /* end for j */
        }  /* end for i */

        /* Copy the results from this sweep into the final results to either
           be returned or to be used in the next sweep */
        memcpy (ltaero, taero, nlines * nsamps * sizeof (float));
        memcpy (ltresi, tresi, nlines * nsamps * sizeof (float));

        /* If there was a hole in this sweep, then increment the window size
           and sweep again.  The increased window size allows us to keep the
           original aerosol values in the interpolation and minimizes the values
           that were previously interpolated.  Copy the results from this sweep
           into the local copy of taero and tresi so that the interpolated
           values are used in the next sweep. */
        if (hole)
        {
            /* Increase the window size */
            step *= 2;

            if (step > 100)
            {
                /* Drop to 30% of non-fill pixels to be required */
                valid_percent = 0.3;
            }
        }
    }  /* end while hole */

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  aerosol_salt_and_pepper_fill

PURPOSE:  Performs interpolation at the single pixel level where the aerosol
holes are not more than a couple of pixels.  This algorithm uses a 3x3 window
around the current pixel to interpolate the aerosol value.  Clouds, cirrus,
and water pixels are not used to interpolation.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred allocating memory or copying taero
SUCCESS        Successful completion

NOTES:
1. A local copy of the aerosols is made and used so that the interpolated
   aerosol pixels are not used to interpolate other aerosol pixels.  Only the
   original aerosol values are used for the interpolation.
******************************************************************************/
int aerosol_salt_and_pepper_fill
(
    int nlines,       /* I: number of lines in the aerosol buffer */
    int nsamps,       /* I: number of samples in the aerosol buffer */
    uint8 *cloud,     /* I: bit-packed value that represent clouds,
                            nlines x nsamps */
    float *tresi,     /* I: residuals for each pixel, nlines x nsamps;
                            tresi < 0.0 flags water pixels and pixels with
                            high residuals */
    float *taero      /* I/O: aerosol values for each pixel, nlines x nsamps */
)
{
    char FUNC_NAME[] = "aerosol_salt_and_pepper_fill";   /* function name */
    char errmsg[STR_SIZE];  /* error message */
    int i, j, k, l;    /* looping variable for pixels */
    int curr_pix;      /* current pixel in the image to be filled */
    int win_pix;       /* current pixel in the line,sample window */
    int nbaot;         /* number of AOT pixels (non-cloud/water) for aerosol
                          interpolation */
    double aaot;       /* average of AOT */
    double sresi;      /* sum of 1 / residuals */
    float *ltaero = NULL;  /* local copy of aerosols for filling */
    float *ltresi = NULL;  /* local copy of residuals for filling */

    /* Allocate space for local copy of the aerosols and residuals */
    ltaero = calloc (nlines*nsamps, sizeof (float));
    if (ltaero == NULL)
    {
        sprintf (errmsg, "Error allocating memory for a local copy of taero");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    ltresi = calloc (nlines*nsamps, sizeof (float));
    if (ltresi == NULL)
    {
        sprintf (errmsg, "Error allocating memory for a local copy of tresi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Make a local copy of the aerosols and residuals which will remain a
       clean copy of the unfilled data to be used for interpolation */
    memcpy (ltaero, taero, nlines * nsamps * sizeof (float));
    if (ltaero == NULL)
    {
        sprintf (errmsg, "Error copying memory for a local copy of taero");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    memcpy (ltresi, tresi, nlines * nsamps * sizeof (float));
    if (ltresi == NULL)
    {
        sprintf (errmsg, "Error copying memory for a local copy of tresi");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Perform salt and pepper aerosol interpolation to fill the small holes
       in the aerosol data */
    printf ("Performing aerosol salt and pepper fill ...\n");
#ifdef _OPENMP
    #pragma omp parallel for private (i, j, k, l, nbaot, aaot, sresi, curr_pix, win_pix)
#endif
    for (i = 0; i < nlines; i++)
    {
        /* Determine the current pixel in the buffer to start this line */
        curr_pix = i * nsamps;
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            /* If the current pixel is a hole in the aerosols and it's not
               cirrus, cloud, or water, then look at filling it */
            if ((ltresi[curr_pix] < 0) &&
                (!btest (cloud[curr_pix], CIR_QA)) &&
                (!btest (cloud[curr_pix], CLD_QA)) &&
                (!btest (cloud[curr_pix], WAT_QA)))
            {
                /* Initialize our variables for this pixel */
                nbaot = 0;
                aaot = 0.0;
                sresi = 0.0;

                /* Check the 3x3 window around the current pixel */
                for (k = i - 1; k <= i + 1; k++)
                {
                    /* Make sure the line is valid */
                    if (k < 0 || k >= nlines)
                        continue;
    
                    win_pix = k * nsamps + j - 1;
                    for (l = j - 1; l <= j + 1; l++, win_pix++)
                    {
                        /* Make sure the sample is valid */
                        if (l < 0 || l >= nsamps)
                            continue;
    
                        /* Check for clear pixels with positive residuals */
                        if ((ltresi[win_pix] > 0) && (cloud[win_pix] == 0))
                        {
                            /* Increment the count for this window and add
                               the current residuals and aerosols for the
                               interpolation */
                            nbaot++;
                            aaot += ltaero[win_pix] / ltresi[win_pix];
                            sresi += 1.0 / ltresi[win_pix];
                        }
                    }  /* for l */
                }  /* for k */

                /* If more than half the 3x3 window is non-zero and clear
                   pixels, then fill the current pixel from the window. The
                   pixel will be filled in the original buffer. */
                if (nbaot >= 5)
                {
                    taero[curr_pix] = aaot / sresi;
                    tresi[curr_pix] = 1.0;
                }
            }  /* end if taero == 0.0 */
        }  /* end for j */
    }  /* end for i */

    /* Free memory */
    free (ltaero);
    free (ltresi);

    return (SUCCESS);
}


/******************************************************************************
MODULE:  aerosol_interpolation_ndd

PURPOSE:  Performs interpolation for the zero-valued aerosol pixels using
Newton's Divided Difference Interpolating Polynomials. This function applies
an Nth order interpolation polynomial to fill the aerosols in the line
direction. Once a gap is found, the first half of the gap is interpolated in
the left to right direction.  Then the function will interpolated the second
half of the gap in the right to left direction.

RETURN VALUE: N/A

NOTES:
1. Does not use the water, cloud, or cirrus pixels as part of the interpolation.
2. Aerosol values are modified inline in the taero input buffer.
3. The algorithm used is from http://pegasus.cc.ucf.edu/~klee/EGN3420/Notes/Interp_NDD.pdf
******************************************************************************/
#define Nth_ORDER 10   /* N+1 pixels are needed for an Nth order polynomial */
void aerosol_interpolation_ndd
(
    int nlines,       /* I: number of lines in the aerosol buffer */
    int nsamps,       /* I: number of samples in the aerosol buffer */
    uint8 *cloud,     /* I: bit-packed value that represent clouds,
                            nlines x nsamps */
    float *tresi,     /* I: residuals for each pixel, nlines x nsamps;
                            tresi < 0.0 flags water pixels and pixels with
                            high residuals */
    float *taero      /* I/O: aerosol values for each pixel, nlines x nsamps */
)
{
    int i, j;           /* looping variable for pixels */
    int hole;           /* looping variable for the current hole */
    long curr_pix;      /* current pixel in the image to be filled */
    long hole_pix;      /* current pixel in the hole */
    int line_hole_pix;  /* current pixel in the hole within the local line */
    long interp_pix;    /* current pixel in the interpolation observations */
    int line_interp_pix;  /* current pixel in the interpolation observations
                             within the local line */
    int hole_size;      /* size/width of the current aerosol hole */
    int nvals;          /* number values used in the interpolation */
    float x[Nth_ORDER+1]; /* array of sample locations to be used for
                             interpolation */
    float interp_aot;   /* taero produced from interpolation at sample x */
    float aero_vals[Nth_ORDER+1]; /* array of aerosol values to be used for the
                                     interpolation */

    /* Aerosol interpolation. Does not use water, cloud, or cirrus pixels. */
    printf ("Performing NDDIP aerosol interpolation ...\n");
    for (i = 0; i < nlines; i++)
    {
        /* Determine the current pixel in the buffer to start this line */
        curr_pix = i * nsamps;

        /* Interpolate holes from left to right, filling in the left half of
           the hole */
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            /* If the current pixel is a hole in the aerosols and it's not
               cirrus, cloud, or water, then look at filling it */
            if ((tresi[curr_pix] < 0) &&
                (!btest (cloud[curr_pix], CIR_QA)) &&
                (!btest (cloud[curr_pix], CLD_QA)) &&
                (!btest (cloud[curr_pix], WAT_QA)))
            {
                /* Determine the width of the hole in this line */
                line_hole_pix = j + 1;
                hole_pix = curr_pix + 1;
                while ((line_hole_pix < nsamps) && (tresi[hole_pix] < 0) &&
                    (!btest (cloud[hole_pix], CIR_QA)) &&
                    (!btest (cloud[hole_pix], CLD_QA)) &&
                    (!btest (cloud[hole_pix], WAT_QA)))
                {
                    line_hole_pix++;
                    hole_pix++;
                }
                hole_size = line_hole_pix - j;

                /* Loop through the left half of the hole and interpolate each
                   pixel */
                line_hole_pix = j;
                hole_pix = curr_pix;
                for (hole = 0; hole < (int) roundf (hole_size*0.5);
                     hole++, hole_pix++, line_hole_pix++)
                {
                    /* Stack up the values to the left of the hole for our
                       Nth-order polynomial. Nth-order polynomial requires N+1
                       values. */
                    nvals = 0;
                    line_interp_pix = line_hole_pix-1;
                    interp_pix = hole_pix-1;
                    while (nvals < Nth_ORDER+1 && line_interp_pix >= 0)
                    {
                        /* Don't use clouds, cirrus, or water */
                        /* TODO: add fill to this test as well as the later test*/
                        if ((!btest (cloud[interp_pix], CIR_QA)) &&
                            (!btest (cloud[interp_pix], CLD_QA)) &&
                            (!btest (cloud[interp_pix], WAT_QA)))
                        {
//                            x[nvals] = (j - line_interp_pix) *
//                                (1 - tresi[interp_pix]);
                            x[nvals] = (j - line_interp_pix);
                            aero_vals[nvals] = taero[interp_pix];
                            nvals++;
                        }

                        /* Look at the next pixel */
                        line_interp_pix--;
                        interp_pix--;
                    }

                    /* Interpolate current pixel, using pixels to the left */
                    nddip (x, aero_vals, nvals, 0, &interp_aot);
                    taero[hole_pix] = interp_aot;
                    tresi[hole_pix] = 1.0;
                }  /* end for hole */

                /* Move the current sample to end of the hole and continue.
                   Skip the next pixel after the hole, because it's already
                   been tested and is known to not be a hole. */
                j += hole_size + 1;
                curr_pix += hole_size + 1;
            } /* end if this is a hole */
        }  /* end for j */

        /* Interpolate holes from right to left, filling in the right half of
           the hole. There's no need to determine the size of the hole since
           all holes have already been filled from left to right. The remaining
           holes will be filled from right to left. */
        curr_pix--;  /* end of the current line */
        for (j = nsamps-1; j >= 0; j--, curr_pix--)
        {
            /* If the current pixel is a hole in the aerosols and it's not
               cirrus, cloud, or water, then look at filling it */
            if ((tresi[curr_pix] < 0) &&
                (!btest (cloud[curr_pix], CIR_QA)) &&
                (!btest (cloud[curr_pix], CLD_QA)) &&
                (!btest (cloud[curr_pix], WAT_QA)))
            {
                /* Stack up the values to the right of the hole for our
                   Nth-order polynomial. Nth-order polynomial requires N+1
                   values. */
                nvals = 0;
                line_interp_pix = j + 1;
                interp_pix = curr_pix + 1;
                while (nvals < Nth_ORDER+1 && line_interp_pix < nsamps)
                {
                    /* Don't use clouds, cirrus, or water */
                    if ((!btest (cloud[interp_pix], CIR_QA)) &&
                        (!btest (cloud[interp_pix], CLD_QA)) &&
                        (!btest (cloud[interp_pix], WAT_QA)))
                    {
/*                        x[nvals] = (line_interp_pix - j) *
                            (1 - tresi[interp_pix]); */
                        x[nvals] = (line_interp_pix - j);
                        aero_vals[nvals] = taero[interp_pix];
                        nvals++;
                    }

                    /* Look at the next pixel */
                    line_interp_pix++;
                    interp_pix++;
                }

                /* Interpolate the current pixel, using pixels to the right */
                /* Call the NDDIP function wiht x, aero_vals, and nvals */
                nddip (x, aero_vals, nvals, 0, &interp_aot);
                taero[curr_pix] = interp_aot;
                tresi[curr_pix] = 1.0;
            } /* end if this is a hole */
        }  /* end for j */
    }  /* end for i */
}


/******************************************************************************
MODULE:  nddip (Newton's Divided Difference Interpolating Polynomial)

PURPOSE:  Performs Newton's divided different interpolation using the x and y
values provided for the nth-order polynomial.

RETURN VALUE: N/A

NOTES:
1. The algorithm used is from http://pegasus.cc.ucf.edu/~klee/EGN3420/Notes/Interp_NDD.pdf

Given n+1 data points (not in any particular order) then
fn(x) = b0 + b1(x-x0) + b2(x-x0)(x-x1) + ... + bn(x-x0)(x-x1)...(x-x{n-1})
where
  b0 = f[x0] = f(x0)
  b1 = f[x1, x0] = (f(x1) - f(x0)) / (x1 - x0)
  b2 = f[x2, x1, x0] = (f[x2, x1] - f[x1, x0]) / (x2 - x0) =
       ((f(x2) - f(x1)) / (x2-x1)) - ((f(x1) - f(x0)) / (x1 - x0))} / (x2 - x0)
  bn = f[xn, x{n-1}, ... x0] =
       (f[xn, x{n-1},.. x2, x1] - f[x{n-1}, x{n-2}, ... x1, x0]) / (xn - x0)

This can be broken into a tree like the following.  For n observations, there
are n bx values that need to be computed in n-1 different tree levels.  The
results from each level feed the computation in the next level.  And, at each
level the first value is the bx value.

         Level0       Level1           Level2               Level3
         ------       ------           ------               ------
         (b0)
  x0 --> f(x0)   \     (b1)
                  \ f[x1, x0]   \       (b2)
                  /              \ f[x2, x1, x0]   \         (b3)
  x1 --> f(x1)   X               /                  \ f[x3, x2, x1, x0]
                  \ f[x2, x1]   X                   /
                  /              \ f[x3, x2, x1]   /
  x2 --> f(x2)   X               /
                  \ f[x3, x2]   /
                  /
  x3 --> f(x3)   /


******************************************************************************/
void nddip
(
    float *x_arr,   /* I: array of input x values (0 to n) */
    float *y_arr,   /* I: array of input f(x) values (0 to n) */
    int nobs,       /* I: number of observations in each array; the polynomial
                          will be an (nobs-1) order polynomial */
    float x,        /* I: x value for the f(x) to be interpolated */
    float *y        /* I: interpolated f(x) value */
)
{
    int i;             /* looping variables */
    int j;             /* current NDDIP tree level */
    int n;             /* local variable for the number of observations */
    float f1, f2;      /* f(n) values */
    float p[Nth_ORDER+1];  /* polynomial coeffients */
    float local_x[Nth_ORDER+1];  /* local copy of the x observations */
    float local_y[Nth_ORDER+1];  /* local copy of the y observations */
 
    /* Make local copy of the x and y observations because they are
       overwritten */
    memcpy (&local_x[0], x_arr, nobs * sizeof (float));
    memcpy (&local_y[0], y_arr, nobs * sizeof (float));

    /* Initializations */
    j = 1;  /* start with level 1 which had nobs-1 computations; b0 is known */
    f1 = 1.0;
    f2 = local_y[0];  /* b0 */
    n = nobs;

    /* Loop until all the coefficients and computations have been completed
       for each of the observations */
    do
    {
        /* Calculate the jth-level in the tree of the polynomial for f(xn),
           which has n-1 calculations.  y[0] after this loop contains
           bx for the jth level (i.e. bj) */
        for (i = 0; i < n-1; i++)
        {
            p[i] = ((local_y[i+1] - local_y[i]) / (local_x[i+j] - local_x[i]));
            local_y[i] = p[i]; /* save to be used in next factor */
        }

        /* Compute the (k-x0)(k-x1)...(k-x{n-1}) for the current level in
           the tree */
        f1 = 1.0;
        for (i = 0; i < j; i++)
            f1 *= (x - local_x[i]);

        /* Add the bx * (k-x0)(k-x1)... item to the sum for current level j */
        f2 += (local_y[0] * f1);

        /* Decrement the number of functions in the level and increment j to
           add the next element in the current level */
        n--;
        j++;
    } while (n != 1);

    /* Get the final f(x) value */
    *y = f2;
}

