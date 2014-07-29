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
    float erelc[16],                 /* I: band ratio variable */
    float troatm[16],                /* I: atmospheric reflectance table */
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
    double ogtransa1[16],            /* I: other gases transmission coeff */
    double ogtransb0[16],            /* I: other gases transmission coeff */
    double ogtransb1[16],            /* I: other gases transmission coeff */
    double wvtransa[16],             /* I: water vapor transmission coeff */
    double wvtransb[16],             /* I: water vapor transmission coeff */
    double oztransa[16],             /* I: ozone transmission coeff */
    float *raot,
    float *residual                  /* O: model residual */
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
    float tgo;
    float roatm;            /* atmospherice reflectance */
    float ttatmg;
    float satm;             /* spherical albedo */
    float xrorayp;          /* molecular reflectance */
	double aratio1, aratio2;
    double pratio;           /* targeted ratio between the surface reflectance
                               in bands */
    double eratio;
    double eaot;             /* estimate of AOT */
    double th1, th3;
    double peratio;
    double pros1, pros3;     /* predicted surface reflectance */

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
            if (iter > 0 && iaot == 0)
                raot550nm = raot550nm * 0.5;
            else if (iter > 0 && iaot > 0)   /* GAIL don't we also need to check iaot > 0? Do we even care about iter > 0?? */
                raot550nm = (raot550nm + aot550nm[iaot-1]) * 0.5;

            /* Atmospheric correction for band 3 */
            iband = iband3;
            retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres, tpres,
                aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
                sphalbt, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
                uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                wvtransb, oztransa, troatm[iband], &roslamb, &tgo, &roatm,
                &ttatmg, &satm, &xrorayp);
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
                sphalbt, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
                uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                wvtransb, oztransa, troatm[iband], &roslamb, &tgo, &roatm,
                &ttatmg, &satm, &xrorayp);
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
    {  /* early break out if the ratios are not valid */
        *residual = 999999.0;
        *raot = raot550nm;
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
        tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
        ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb, oztransa,
        troatm[iband], &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp);
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
        tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
        ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb, oztransa,
        troatm[iband], &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp);
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
        tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
        ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb, oztransa,
        troatm[iband], &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp);
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
        tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
        ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb, oztransa,
        troatm[iband], &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp);
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
                    xtvmin, sphalbt, tsmax, tsmin, nbfic, nbfi, tts, indts,
                    ttv, uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1,
                    wvtransa, wvtransb, oztransa, troatm[iband], &roslamb,
                    &tgo, &roatm, &ttatmg, &satm, &xrorayp);
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
                    xtvmin, sphalbt, tsmax, tsmin, nbfic, nbfi, tts, indts,
                    ttv, uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1,
                    wvtransa, wvtransb, oztransa, troatm[iband], &roslamb,
                    &tgo, &roatm, &ttatmg, &satm, &xrorayp);
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
                    xtvmin, sphalbt, tsmax, tsmin, nbfic, nbfi, tts, indts,
                    ttv, uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1,
                    wvtransa, wvtransb, oztransa, troatm[iband], &roslamb,
                    &tgo, &roatm, &ttatmg, &satm, &xrorayp);
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
                    xtvmin, sphalbt, tsmax, tsmin, nbfic, nbfi, tts, indts,
                    ttv, uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1,
                    wvtransa, wvtransb, oztransa, troatm[iband], &roslamb,
                    &tgo, &roatm, &ttatmg, &satm, &xrorayp);
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
    *residual = fabs (ros3 - ros1 * pratio);
    for (iband = 0; iband < 6; iband++)
    {
        if (erelc[iband] > 0.0)
        {
            retval = atmcorlamb2 (xts, xtv, xfi, raot550nm, iband, pres, tpres,
                aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
                sphalbt, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
                uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
                wvtransb, oztransa, troatm[iband], &roslamb, &tgo, &roatm,
                &ttatmg, &satm, &xrorayp);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Performing lambertian atmospheric correction "
                    "type 2.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            *residual += fabs (roslamb - ros1 * (erelc[iband] / erelc[iband1]));
        }
    }

    /* Successful completion */
    return (SUCCESS);
}

