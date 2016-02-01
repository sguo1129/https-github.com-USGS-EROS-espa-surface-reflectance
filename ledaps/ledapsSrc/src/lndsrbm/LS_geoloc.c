#define _GNU_SOURCE
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <math.h> 
#include "cproj.h"
#include "cproj_prototypes.h"
#include "gctp_prototypes.h"

/* NOTE: Most of these variables and functions were pulled from the GCTP
 * package to handle the projection transformations.  Also, this code is not
 * set up to mix and match projections.  If one projection is initialized,
 * then another projection is initialized, it will overwrite the first
 * projection information.
 */
static double scale_factor   = 0.9996;   /* scale factor */
double pixel_size;      /* pixel size */
double sin_orien;       /* sine of the orientation angle */
double cos_orien;       /* cosine of the orientation angle */
double ul_corner[2];    /* UL projection x,y for UL (not center) of the pixel */

/* Assigns values to the semimajor axis, semiminor axis, and radius of sphere.
   Initializes the forward and/or inverse mapping function. */
int LSsphdz (
    char *projection,     /* I: projection name string */
    float coordinates[8], /* I: general coordinate info */
    double *parm,         /* I: Projection parameters */
    double *radius,       /* O: Radius of the sphere */
    double corner[2]      /* I: UL x,y for UL corner */
)
{
    double r_major;   /* major axis */
    double r_minor;   /* minor axis */
    double orient;    /* orientation */
    int ret = 0;      /* return code */
    long isph;        /* spheroid code number also known as datum */
    long zone;        /* zone code */
    
    /* Initialize global variables for the mapping */
    zone = (long) coordinates[4];
    isph = (long) coordinates[5];
    orient = coordinates[6];
    pixel_size = coordinates[7];

    sin_orien = sin (orient);
    cos_orien = cos (orient);

    ul_corner[0] = corner[0];
    ul_corner[1] = corner[1];
    
    /* Initialize the variables using GCTP */
    ret = sphdz (isph, parm, &r_major, &r_minor, radius);
    if (ret != 0)
        return (ret);
    
    if (zone == 0)
        zone = 31L;
      
    /* Do the forward or inverse transformation setup, depending on whether
       inverse is specified */
#ifdef INV
    if (!strcmp (projection, "GCTP_UTM"))
        ret = utminvint(r_major, r_minor, scale_factor, zone);
    else if (!strcmp (projection, "GCTP_PS"))
        ret = psinvint(r_major, r_minor, parm[4], parm[5], parm[6], parm[7]);
    else if (!strcmp (projection, "GCTP_ALBERS"))
        ret = alberinvint(r_major, r_minor, parm[2], parm[3], parm[4], parm[5],
            parm[6], parm[7]);
#else  
    if (!strcmp (projection, "GCTP_UTM"))
        ret = utmforint(r_major, r_minor, scale_factor, zone);
    else if (!strcmp (projection, "GCTP_PS"))
        ret = psforint(r_major, r_minor, parm[4], parm[5], parm[6], parm[7]);
    else if (!strcmp (projection, "GCTP_ALBERS"))
        ret = alberforint(r_major, r_minor, parm[2], parm[3], parm[4], parm[5],
            parm[6], parm[7]);
#endif  
    
    return (ret);
}


/* Universal Transverse Mercator inverse equations--mapping line,sample to
   x,y to lat,long 
  -----------------------------------------------------------------------*/
int LSutminv (
    double s,      /* I: sample */
    double l,      /* I: line */
    double *lon,   /* O: longitude (degrees) */
    double *lat    /* O: latitude (degrees) */
)
{
    int ret = 0;             /* return value */
    double x, y;             /* x,y projection coords */
    double dl, dp, dy, dx;   /* delta line, sample and x,y values */
    
    /* Calculate the x,y from the line,sample */
    dl = (l + 0.5) * pixel_size;
    dp = (s + 0.5) * pixel_size;
    dy = (dp * sin_orien) - (dl * cos_orien);
    dx = (dp * cos_orien) + (dl * sin_orien);
    y = ul_corner[1] + dy;
    x = ul_corner[0] + dx;
    
    /* Do the inverse mapping */
    ret = utminv (x, y, lon, lat);
    if (ret != 0)
    {
        printf ("Error in the inverse UTM mapping");
        return (ret);
    }
    
    /* Convert lat/long to degrees */
    *lat *= R2D;
    *lon *= R2D;
    
    return (0);
}


/* Universal Transverse Mercator inverse equations--mapping lat,long to x,y
   to line,sample
  -------------------------------------------------------------------------*/
int LSutmfor (
    double *s,    /* O: sample */
    double *l,    /* O: line */
    double lon,   /* I: longitude (degrees) */
    double lat    /* I: latitude (degrees) */
)
{
    int ret = 0;     /* return value */
    double x, y;     /* x,y projection coords */
    double dy, dx;   /* delta x, y */
    
    /* Convert lat/long from degrees to radians */
    lat *= D2R;
    lon *= D2R;

    /* Do the forward mapping */
    ret = utmfor (lon, lat, &x, &y);
    if (ret != 0)
    {
        printf ("Error in the forward UTM mapping");
        return (ret);
    }
    
    /* Convert the x,y back to line,sample */
    x -= ul_corner[0];
    y -= ul_corner[1];
    
    dx = (x * cos_orien) + (y * sin_orien);
    dy = (x * sin_orien) - (y * cos_orien);
    
    *s = dx / pixel_size - 0.5;
    *l = dy / pixel_size - 0.5;
    
    return(0);
}


/* Polar Stereographic inverse equations--mapping line,sample to x,y to
   lat/long
  ---------------------------------------------------------------------*/
int LSpsinv (
    double s,       /* I: sample */
    double l,       /* I: line */
    double *lon,    /* O: longitude (degrees) */
    double *lat     /* O: latitude (degrees) */
)
{
    int ret = 0;             /* return value */
    double x, y;             /* x,y projection coords */
    double dl, dp, dy, dx;   /* delta line, sample and x,y values */

    /* Calculate the x,y from the line,sample */
    dl = (l + 0.5) * pixel_size;
    dp = (s + 0.5) * pixel_size;
    dy = (dp * sin_orien) - (dl * cos_orien);
    dx = (dp * cos_orien) + (dl * sin_orien);
    y = ul_corner[1] + dy;
    x = ul_corner[0] + dx;

    /* Do the inverse mapping */
    ret = psinv (x, y, lon, lat);
    if (ret != 0)
    {
        printf ("Error in the inverse PS mapping");
        return (ret);
    }

    /* Convert lat/long to degrees */
    *lat *= R2D;
    *lon *= R2D;
    
    return(0);
}

/* Polar Stereographic forward equations--mapping lat,long to x,y to
   line,sample
  ------------------------------------------------------------------*/
int LSpsfor(
    double *s,    /* O: sample */
    double *l,    /* O: line */
    double lon,   /* I: longitude (degrees) */
    double lat    /* I: latitude (degrees) */
)
{
    int ret = 0;    /* return value */
    double x, y;    /* x,y projection coords */
    double dy, dx;  /* delta x, y */

    /* Convert lat/long from degrees to radians */
    lat *= D2R;
    lon *= D2R;

    /* Do the forward mapping */
    ret = psfor (lon, lat, &x, &y);
    if (ret != 0)
    {
        printf ("Error in the forward PS mapping");
        return (ret);
    }

    /* Convert the x,y back to line,sample */
    x -= ul_corner[0];
    y -= ul_corner[1];
    
    dx = (x * cos_orien) + (y * sin_orien);
    dy = (x * sin_orien) - (y * cos_orien);
    
    *s = dx / pixel_size - 0.5;
    *l = dy / pixel_size - 0.5;

    return(0);
}


/* Albers inverse equations--mapping line,sample to x,y to lat/long
  ---------------------------------------------------------------------*/
int LSalbersinv (
    double s,       /* I: sample */
    double l,       /* I: line */
    double *lon,    /* O: longitude (degrees) */
    double *lat     /* O: latitude (degrees) */
)
{
    int ret = 0;             /* return value */
    double x, y;             /* x,y projection coords */
    double dl, dp, dy, dx;   /* delta line, sample and x,y values */

    /* Calculate the x,y from the line,sample */
    dl = (l + 0.5) * pixel_size;
    dp = (s + 0.5) * pixel_size;
    dy = (dp * sin_orien) - (dl * cos_orien);
    dx = (dp * cos_orien) + (dl * sin_orien);
    y = ul_corner[1] + dy;
    x = ul_corner[0] + dx;

    /* Do the inverse mapping */
    ret = alberinv (x, y, lon, lat);
    if (ret != 0)
    {
        printf ("Error in the inverse Albers mapping");
        return (ret);
    }

    /* Convert lat/long to degrees */
    *lat *= R2D;
    *lon *= R2D;
    
    return(0);
}

/* Albers forward equations--mapping lat,long to x,y to line,sample
  ------------------------------------------------------------------*/
int LSalbersfor(
    double *s,    /* O: sample */
    double *l,    /* O: line */
    double lon,   /* I: longitude (degrees) */
    double lat    /* I: latitude (degrees) */
)
{
    int ret = 0;    /* return value */
    double x, y;    /* x,y projection coords */
    double dy, dx;  /* delta x, y */

    /* Convert lat/long from degrees to radians */
    lat *= D2R;
    lon *= D2R;

    /* Do the forward mapping */
    ret = alberfor (lon, lat, &x, &y);
    if (ret != 0)
    {
        printf ("Error in the forward Albers mapping");
        return (ret);
    }

    /* Convert the x,y back to line,sample */
    x -= ul_corner[0];
    y -= ul_corner[1];
    
    dx = (x * cos_orien) + (y * sin_orien);
    dy = (x * sin_orien) - (y * cos_orien);
    
    *s = dx / pixel_size - 0.5;
    *l = dy / pixel_size - 0.5;

    return(0);
}
