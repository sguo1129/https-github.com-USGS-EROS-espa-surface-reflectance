/*****************************************************************************
FILE: subaeroret.c
  
PURPOSE: Contains functions for handling the atmosperic corrections.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
6/25/2014    Gail Schmidt     Conversion of the original FORTRAN code delivered
                              by Eric Vermote, NASA GSFC
8/14/2014    Gail Schmidt     Updated for v1.3 delivered by Eric Vermote

NOTES:
*****************************************************************************/
#include "lut_subr.h"

/******************************************************************************
MODULE:  subaeroret

PURPOSE:  Main driver for the atmospheric correction.  This subroutine reads
the lookup table (LUT) and performs the atmospheric corrections.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the LUT or doing the correction
SUCCESS        Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
6/25/2014    Gail Schmidt     Conversion of the original FORTRAN code delivered
                              by Eric Vermote, NASA GSFC
7/22/2014    Gail Schmidt     Number of iterations doesn't need to be returned
                              from this routine since it isn't used
7/22/2014    Gail Schmidt     Cleaned up unused ogtransa0, ogtransc0,
                              ogtransc1, wvtransc arrays.  Made the rest of
                              these transmission arrays doubles and hard-coded
                              their static values in this code vs. reading
                              from a static ASCII file.
7/29/2014    Gail Schmidt     Defined a static NSR_BANDS variable for the
                              variables that refer to the surface reflectance
                              band-related bands (ogtrans, wvtrans, tauray,
                              erelc, etc.)  These previously were of size 16.
                              Only compute the residual for the NSR_BANDS as
                              well vs. doing the residual computation for 16
                              bands.
8/14/2014    Gail Schmidt     Updated for v1.3 delivered by Eric Vermote

NOTES:
******************************************************************************/
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
    int iband;              /* looping variable for bands */
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
    float tgo;
    float roatm;            /* atmospherice reflectance */
    float ttatmg;
    float satm;             /* spherical albedo */
    float xrorayp;          /* molecular reflectance */
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
    nit = 0;
    pros1 = 0.0;
    pros3 = 0.0;

    /* The ratio decreases as the AOT increases.  The exit conditions in this
       loop are when two values of AOT can be found that bracket the predicted
       ratio (pratio). */
    while ((iaot < 22) && (aratio1 > pratio) && (ros1 > th1) && (ros3 > th3) &&
        ((aratio1 - 0.01) < aratio2) && (nit < 30))
    {
        nit++;
        ros1 = -1.0;
        ros3 = -1.0;

        /* If flagn is set... start converge to the AOT bounds by dichotomy to
           increase the accuracy of the retrieval */
        if (!flagn)
            raot550nm = aot550nm[iaot];
        else
            raot550nm = (raot1 + aot550nm[iaot]) * 0.5;

        /* Loop until convergence */
        iter = 0;
        while ((ros1 < th1) || (ros3 < th3)) 	
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
                        roslamb, pratio, raot550nm, xts, xtv, xfi, pres,
                        uoz, uwv, erelc, troatm, tpres, aot550nm, rolutt,
                        transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                        normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv,
                        tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                        wvtransb, oztransa, residual, snext);
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
            iband = iband3;
            retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres, tpres,
                aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
                sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv,
                uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                wvtransb, oztransa, troatm[iband], &roslamb, &tgo, &roatm,
                &ttatmg, &satm, &xrorayp, &next);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Performing lambertian atmospheric correction "
                    "type 2.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            ros3 = roslamb;

            /* Atmospheric correction for band 1 */
            iband = iband1;
            retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres, tpres,
                aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
                sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv,
                uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                wvtransb, oztransa, troatm[iband], &roslamb, &tgo, &roatm,
                &ttatmg, &satm, &xrorayp, &next);
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
            pratio, raot550nm, xts, xtv, xfi, pres, uoz, uwv, erelc, troatm,
            tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
            sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv,
            tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
            oztransa, residual, snext);
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

    /* The estimated AOT is refined by recomputing by performing an additional
       iteration of atmospheric correction */
    /* Atmospheric correction for band 3 */
    raot550nm = eaot;
    iband = iband3;
    retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres, tpres,
        aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
        normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
        ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb, oztransa,
        troatm[iband], &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Performing lambertian atmospheric correction "
            "type 2.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    ros3 = roslamb;

    /* Atmospheric correction for band 1 */
    iband = iband1;
    retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres, tpres,
        aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
        normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
        ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb, oztransa,
        troatm[iband], &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
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
    iband = iband3;
    retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres, tpres,
        aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
        normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
        ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb, oztransa,
        troatm[iband], &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Performing lambertian atmospheric correction "
            "type 2.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    ros3 = roslamb;

    /* Atmospheric correction for band 1 */
    iband = iband1;
    retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres, tpres,
        aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
        normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
        ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb, oztransa,
        troatm[iband], &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
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
                iband = iband3;
                retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres,
                    tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
                    xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts,
                    indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                    ogtransb1, wvtransa, wvtransb, oztransa, troatm[iband],
                    &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
                ros3 = roslamb;

                /* Atmospheric correction for band 1 */
                iband = iband1;
                retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres,
                    tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
                    xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts,
                    indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                    ogtransb1, wvtransa, wvtransb, oztransa, troatm[iband],
                    &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
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
                iband = iband3;
                retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres,
                    tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
                    xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, 
                    indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                    ogtransb1, wvtransa, wvtransb, oztransa, troatm[iband],
                    &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
                ros3 = roslamb;

                /* Atmospheric correction for band 1 */
                iband = iband1;
                retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres,
                    tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
                    xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts,
                    indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                    ogtransb1, wvtransa, wvtransb, oztransa, troatm[iband],
                    &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next);
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
        pratio, raot550nm, xts, xtv, xfi, pres, uoz, uwv, erelc, troatm,
        tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
        sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv,
        tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
        oztransa, residual, snext);
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

PURPOSE:  Computes the model residual for the subaeroret function.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the LUT or doing the correction
SUCCESS        Successful completion

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
8/14/2014    Gail Schmidt     Original development

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
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
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
    float roatm;            /* atmospherice reflectance */
    float ttatmg;
    float satm;             /* spherical albedo */
    float xrorayp;          /* molecular reflectance */

/**** GAIL -- should this be DN_BAND7 vs. DN_BAND6?? */
    /* Compute the model residual */
    *residual = fabs (ros3 - ros1 * pratio);
    for (iband = 0; iband <= DN_BAND6; iband++)
    {
        if (erelc[iband] > 0.0)
        {
            retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres, tpres,
                aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
                sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv,
                uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                wvtransb, oztransa, troatm[iband], &roslamb, &tgo, &roatm,
                &ttatmg, &satm, &xrorayp, &next);
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

