#define _GNU_SOURCE
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <math.h> 

/*
 * cc -O0 -o LS_geoloc.o -c LS_geoloc.c 
 */

/* NOTE: Most of these variables and functions were pulled from the GCTP
 * package to handle the projection transformations.  Also, this code is not
 * set up to mix and match projections.  If one projection is initialized,
 * then another projection is initialized, it will overwrite the first
 * projection information.
 *
 * Modified on 3/12/2013, Gail Schmidt, USGS EROS
 * Updated to add support for the Polar Stereographic projection, in addition
 *   to the existing support for UTM.
 */
#define PI      3.14159265358979323846
#define HALF_PI 1.57079632679489661923
#define TWO_PI  6.28318530717958647692
#define EPSLN   1.0e-10
#define R2D     57.2957795131
#define D2R     1.745329251994328e-2
#define S2R     4.848136811095359e-6
#define SQUARE(x)       x * x  
#define MAX_VAL 4
#define MAXLONG 2147483647.
#define DBLLONG 4.61168601e18

#define DATMCT 20
static double major[20] = {6378206.4, 6378249.145, 6377397.155, 6378157.5,
                    6378388.0, 6378135.0, 6377276.3452, 6378145.0,
                    6378137.0, 6377563.396, 6377304.063, 6377340.189,
                    6378137.0, 6378155.0, 6378160.0, 6378245.0,
                    6378270.0, 6378166.0, 6378150.0, 6370997.0};

static double minor[20] = {6356583.8, 6356514.86955, 6356078.96284, 6356772.2,
                    6356911.94613, 6356750.519915, 6356075.4133,
                    6356759.769356, 6356752.31414, 6356256.91,
                    6356103.039, 6356034.448, 6356752.314245,
                    6356773.3205, 6356774.719, 6356863.0188,
                    6356794.343479, 6356784.283666, 6356768.337303,
                    6370997.0};
static double scale_factor   = 0.9996;
static double false_easting  = 0.0;
static double false_northing = 0.0;

static double r_major;          /* major axis                           */
static double r_minor;          /* minor axis                           */
static double scale_factor;     /* scale factor                         */
static double lon_center;       /* Center longitude (projection center) */
static double lat_origin;       /* center latitude                      */
static double e0,e1,e2,e3;      /* eccentricity constants               */
static double e,e4,es,esp;      /* eccentricity constants               */
static double ml0;              /* small value m                        */
static long ind;		        /* sphere flag value			*/
static double center_lon;		/* center longitude		*/
static double center_lat;		/* center latitude		*/
static double fac;			/* sign variable		*/
static double mcs;			/* small m			*/
static double tcs;			/* small t			*/

double pixel_size, sin_orien, cos_orien, ul_corner[2];

/* Local prototypes */
int utminvint(double r_maj, double r_min, double scale_fact, long zone);
int utmforint(double r_maj, double r_min, double scale_fact, long zone);
int psinvint (double r_maj, double r_min, double c_lon, double c_lat,
    double false_east, double false_north);
int psforint (double r_maj, double r_min, double c_lon, double c_lat,
    double false_east, double false_north);
double e0fn(double);
double e1fn(double);
double e2fn(double);
double e3fn(double);
double mlfn(double, double, double, double, double);
double e4fn(double x);
double msfnz (double eccent, double sinphi, double cosphi);
double tsfnz (double eccent, double phi, double sinphi);
double asinz(double);
double phi2z (double eccent, double ts, long *flag);
double adjust_lon (double x);


int LSsphdz(char *projection, float coordinates[8], double *parm,
    double *radius, double corner[2])
/*long isph;		/ * spheroid code number also known as datum	*/
/*long zone;		/ * zone code                            	*/
/*double *parm;		/ * projection parameters			*/
/*double *radius;	/ * radius					*/

/* this is meant to be visible from other routines */

{
double r_major;
double r_minor;
int ret;
long isph;		/* spheroid code number also known as datum	*/
long zone;		/* zone code                            	*/

double t_major;		/* temporary major axis				*/
double t_minor;		/* temporary minor axis				*/
long jsph;		/* spheroid code number				*/
double orient;		/* orientation   				*/

zone = (long)coordinates[4];
isph = (long)coordinates[5];
orient     = coordinates[6];
pixel_size = coordinates[7];
sin_orien = sin(orient);
cos_orien = cos(orient);

ul_corner[0] = corner[0];
ul_corner[1] = corner[1];


if (isph < 0)
   {
   t_major = fabs(parm[0]);
   t_minor = fabs(parm[1]);
   
   if (t_major  > 0.0) 
     {
     if (t_minor > 1.0)
        {
        r_major = t_major;
  	r_minor = t_minor;
	*radius = t_major;
        } 
     else
     if (t_minor > 0.0)
        {
        r_major = t_major;
	*radius = t_major;
        r_minor = (sqrt(1.0 - t_minor)) * t_major; 
        }
     else
        {
        r_major = t_major;
	*radius = t_major;
        r_minor = t_major;
        }
     }
   else
   if (t_minor > 0.0)	/* t_major = 0 */
     {
     r_major = major[0];
     *radius = major[0];
     r_minor = minor[0];
     }
   else
     {
     r_major = major[DATMCT - 1];
     *radius = major[DATMCT - 1];
     r_minor = 6370997.0;
     }
  }
else		/* isph >= 0 */
  {
  jsph = abs(isph);
  if (jsph > 19)
     {
     printf("Invalid spheroid selection\n");
     isph = 1;
     jsph = 0;
     }
  r_major = major[jsph];
  r_minor = minor[jsph];
  *radius = major[DATMCT - 1];

  if (zone == 0) zone = 31L;

  
/* Do the forward or inverse transformation setup, depending on whether
   inverse is specified */
#ifdef INV
  if (!strcmp (projection, "GCTP_UTM"))
    ret = utminvint(r_major, r_minor, scale_factor, zone);
  else if (!strcmp (projection, "GCTP_PS"))
    ret = psinvint(r_major, r_minor, parm[4], parm[5], parm[6], parm[7]);
#else  
  if (!strcmp (projection, "GCTP_UTM"))
    ret = utmforint(r_major, r_minor, scale_factor, zone);
  else if (!strcmp (projection, "GCTP_PS"))
    ret = psforint(r_major, r_minor, parm[4], parm[5], parm[6], parm[7]);
#endif  
  }
return(ret);
}


/* Initialize the UTM projection for forward transformations
  ----------------------------------------------------------*/
int utmforint(double r_maj, double r_min, double scale_fact, long zone)
/*double r_maj;			/ * major axis				*/
/*double r_min;			/ * minor axis				*/
/*double scale_fact;		/ * scale factor			*/
/*long zone;			/ * zone number				*/

/* this is meant to be private to this file */

{
double temp;			/* temporary variable			*/

if ((abs(zone) < 1) || (abs(zone) > 60))
   {
   printf("Illegal zone number\n");
   return(11);
   }

r_major = r_maj;
r_minor = r_min;
scale_factor = scale_fact;
lat_origin = 0.0;
lon_center = ((6 * abs(zone)) - 183) * D2R;
false_easting = 500000.0;
false_northing = (zone < 0) ? 10000000.0 : 0.0;

temp = r_minor / r_major;
es = 1.0 - SQUARE(temp);
e = sqrt(es);
e0 = e0fn(es);
e1 = e1fn(es);
e2 = e2fn(es);
e3 = e3fn(es);
ml0 = r_major * mlfn(e0, e1, e2, e3, lat_origin);
esp = es / (1.0 - es);

if (es < .00001)
   ind = 1;
else 
   ind = 0;

return(0);
}


/* Initialize the UTM projection for inverse transformations
  ----------------------------------------------------------*/
int utminvint(double r_maj, double r_min, double scale_fact, long zone)
/*double r_maj;			/ * major axis				*/
/*double r_min;			/ * minor axis				*/
/*double scale_fact;		/ * scale factor			*/
/*long zone;			/ * zone number				*/

/* this is meant to be private to this file */

{
double temp;			/* temporary variables			*/

if ((abs(zone) < 1) || (abs(zone) > 60))
   {
   printf("Illegal zone number\n");
   return(11);
   }

r_major = r_maj;
r_minor = r_min;
scale_factor = scale_fact;
lat_origin = 0.0;
lon_center = ((6 * abs(zone)) - 183) * D2R;
false_easting = 500000.0;
false_northing = (zone < 0) ? 10000000.0 : 0.0;

temp = r_minor / r_major;
es = 1.0 - SQUARE(temp);
e = sqrt(es);
e0 = e0fn(es);
e1 = e1fn(es);
e2 = e2fn(es);
e3 = e3fn(es);
ml0 = r_major * mlfn(e0, e1, e2, e3, lat_origin);
esp = es / (1.0 - es);

if (es < .00001)
   ind = 1;
else 
   ind = 0;

return(0);
}

/* Initialize the Polar Stereographic projection for forward transformations
  --------------------------------------------------------------------------*/
int psforint (double r_maj, double r_min, double c_lon, double c_lat,
    double false_east, double false_north)
//    double r_maj,				/* major axis			*/
//    double r_min,				/* minor axis			*/
//    double c_lon,				/* center longitude		*/
//    double c_lat,				/* center latitude		*/
//    double false_east,			/* x offset in meters		*/
//    double false_north			/* y offset in meters		*/

/* this is meant to be private to this file */

{
double temp;				/* temporary variable		*/
double con1;				/* temporary angle		*/
double sinphi;				/* sin value			*/
double cosphi;				/* cos value			*/

r_major = r_maj;
r_minor = r_min;
false_northing = false_north;
false_easting = false_east;
temp = r_minor / r_major;
es = 1.0 - SQUARE(temp);
e = sqrt(es);
e4 = e4fn(e);
center_lon = c_lon;
center_lat = c_lat;

if (c_lat < 0)
   fac = -1.0;
else
   fac = 1.0;
ind = 0;
if (fabs(fabs(c_lat) - HALF_PI) > EPSLN)
{
   ind = 1;
   con1 = fac * center_lat; 
   sincos(con1,&sinphi,&cosphi);
   mcs = msfnz(e,sinphi,cosphi);
   tcs = tsfnz(e,con1,sinphi);
}

return(0);
}

/* Initialize the Polar Stereographic projection for inverse transformations
  --------------------------------------------------------------------------*/
int psinvint (double r_maj, double r_min, double c_lon, double c_lat,
    double false_east, double false_north)
//    double r_maj,				/* major axis			*/
//    double r_min,				/* minor axis			*/
//    double c_lon,				/* center longitude		*/
//    double c_lat,				/* center latitude		*/
//    double false_east,			/* x offset in meters		*/
//    double false_north			/* y offset in meters		*/

/* this is meant to be private to this file */

{
double temp;				/* temporary variable		*/
double con1;				/* temporary angle		*/
double sinphi;				/* sin value			*/
double cosphi;				/* cos value			*/
double es;                     /* eccentricity squared         */

r_major = r_maj;
r_minor = r_min;
false_easting = false_east;
false_northing = false_north;
temp = r_minor / r_major;
es = 1.0 - SQUARE(temp);
e = sqrt(es);
e4 = e4fn(e);
center_lon = c_lon;
center_lat = c_lat;

if (c_lat < 0)
   fac = -1.0;
else
   fac = 1.0;
ind = 0;
if (fabs(fabs(c_lat) - HALF_PI) > EPSLN)
{
   ind = 1;
   con1 = fac * center_lat; 
   sincos(con1,&sinphi,&cosphi);
   mcs = msfnz(e,sinphi,cosphi);
   tcs = tsfnz(e,con1,sinphi);
}

return(0);
}

int sign(double x)    { if (x < 0.0) return(-1); else return(1);}

double e0fn(double x) {return(1.0-0.25*x*(1.0+x/16.0*(3.0+1.25*x)));}

double e1fn(double x) {return(0.375*x*(1.0+0.25*x*(1.0+0.46875*x)));}

double e2fn(double x) {return(0.05859375*x*x*(1.0+0.75*x));}

double e3fn(double x) {return(x*x*x*(35.0/3072.0));}

double mlfn(double e0, double e1, double e2, double e3, double phi) 
    {return(e0*phi-e1*sin(2.0*phi)+e2*sin(4.0*phi)-e3*sin(6.0*phi));}

double e4fn(double x)
{
 double con;
 double com;
 con = 1.0 + x;
 com = 1.0 - x;
 return (sqrt((pow(con,con))*(pow(com,com))));
}
 
/* function to eliminate roundoff errors in asin */
double asinz (double con)
{
 if (fabs(con) > 1.0)
   {
   if (con > 1.0)
     con = 1.0;
   else
     con = -1.0;
   }
 return(asin(con));
}

/* Function to compute the constant small m which is the radius of
   a parallel of latitude, phi, divided by the semimajor axis */
double msfnz (double eccent, double sinphi, double cosphi)
{
    double con;

    con = eccent * sinphi;
    return((cosphi / (sqrt (1.0 - con * con))));
}

/* Function to compute the constant small t for use in the forward
   computations in the Lambert Conformal Conic and the Polar
   Stereographic projections */
double tsfnz (double eccent, double phi, double sinphi)
{
    double con;
    double com;
  
    con = eccent * sinphi;
    com = .5 * eccent; 
    con = pow(((1.0 - con) / (1.0 + con)),com);
    return (tan(.5 * (HALF_PI - phi))/con);
}

/* Function to compute the latitude angle, phi2, for the inverse of the
   Lambert Conformal Conic and Polar Stereographic projections */
double phi2z (double eccent, double ts, long *flag)
{
    double eccnth;
    double phi;
    double con;
    double dphi;
    double sinpi;
    long i;

    *flag = 0;
    eccnth = .5 * eccent;
    phi = HALF_PI - 2 * atan(ts);
    for (i = 0; i <= 15; i++)
    {
        sinpi = sin(phi);
        con = eccent * sinpi;
        dphi = HALF_PI - 2 * atan(ts *(pow(((1.0 - con)/(1.0 + con)),eccnth))) 
             - phi;
        phi += dphi; 
        if (fabs(dphi) <= .0000000001)
            return(phi);
    }
    printf ("Error: Convergence error in phi2z-conv\n");
    *flag = 002;
    return(002);
}

/* Function to adjust a longitude angle to range from -180 to 180 radians
   added if statments */
double adjust_lon (double x)
{
    long count = 0;
    for(;;)
    {
        if (fabs(x)<=PI)
            break;
        else if (((long) fabs(x / PI)) < 2)
            x = x-(sign(x) *TWO_PI);
        else if (((long) fabs(x / TWO_PI)) < MAXLONG)
            x = x-(((long)(x / TWO_PI))*TWO_PI);
        else if (((long) fabs(x / (MAXLONG * TWO_PI))) < MAXLONG)
            x = x-(((long)(x / (MAXLONG * TWO_PI))) * (TWO_PI * MAXLONG));
        else if (((long) fabs(x / (DBLLONG * TWO_PI))) < MAXLONG)
            x = x-(((long)(x / (DBLLONG * TWO_PI))) * (TWO_PI * DBLLONG));
        else
            x = x-(sign(x) *TWO_PI);
        count++;
        if (count > MAX_VAL)
           break;
    }

    return(x);
}


/* Universal Transverse Mercator inverse equations--mapping line,sample to
   x,y to lat,long 
  -----------------------------------------------------------------------*/
int LSutminv(double s, double l, double *lon, double *lat)
//double x;		/* (I) X projection coordinate 		*/
//double y;		/* (I) Y projection coordinate 		*/
//double *lon;		/* (O) Longitude 				*/
//double *lat;		/* (O) Latitude 				*/

/* this is meant to be visible from other routines */

{
double x, y;
double con,phi;		/* temporary angles				*/
double delta_phi;	/* difference between longitudes		*/
long i;			/* counter variable				*/
double sin_phi, cos_phi, tan_phi;	/* sin cos and tangent values	*/
double c, cs, t, ts, n, r, d, ds;	/* temporary variables		*/
double f, h, g, temp;			/* temporary variables		*/
long max_iter = 6;			/* maximun number of iterations	*/
double dl, dp, dy, dx;

dl = (l + 0.5) * pixel_size;
dp = (s + 0.5) * pixel_size;
dy = (dp * sin_orien) - (dl * cos_orien);
dx = (dp * cos_orien) + (dl * sin_orien);
y = ul_corner[1] + dy;
x = ul_corner[0] + dx;

/* -- start of original GCTP utminv -- */
/* fortran code for spherical form 
   --------------------------------*/
if (ind != 0)
   {
   f = exp(x/(r_major * scale_factor));
   g = .5 * (f - 1/f);
   temp = lat_origin + y/(r_major * scale_factor);
   h = cos(temp);
   con = sqrt((1.0 - h * h)/(1.0 + g * g));
   *lat = asinz(con);
   if (temp < 0)
     *lat = -*lat;
   if ((g == 0) && (h == 0))
     {
     *lon = lon_center;
     return(0);
     }
   else
     {
     *lon = (atan2(g,h) + lon_center);
     return(0);
     }
   }

/* Inverse equations
  -----------------*/
x = x - false_easting;
y = y - false_northing;

con = (ml0 + y / scale_factor) / r_major;
phi = con;
for (i=0;;i++)
   {
   delta_phi = ((con + e1 * sin(2.0*phi) - e2 * sin(4.0*phi) + e3 * sin(6.0*phi))
               / e0) - phi;
   phi += delta_phi;
   if (fabs(delta_phi) <= EPSLN) break;
   if (i >= max_iter) 
      { 
      printf("Latitude failed to converge\n"); 
      return(95);
      }
   }
if (fabs(phi) < HALF_PI)
   {
   cos_phi = cos(phi);
   sin_phi = sin(phi);
   tan_phi = tan(phi);
   c    = esp * SQUARE(cos_phi);
   cs   = SQUARE(c);
   t    = SQUARE(tan_phi);
   ts   = SQUARE(t);
   con  = 1.0 - es * SQUARE(sin_phi); 
   n    = r_major / sqrt(con);
   r    = n * (1.0 - es) / con;
   d    = x / (n * scale_factor);
   ds   = SQUARE(d);
   *lat = phi - (n * tan_phi * ds / r) * (0.5 - ds / 24.0 * (5.0 + 3.0 * t + 
          10.0 * c - 4.0 * cs - 9.0 * esp - ds / 30.0 * (61.0 + 90.0 * t +
          298.0 * c + 45.0 * ts - 252.0 * esp - 3.0 * cs)));
   *lon = lon_center + (d * (1.0 - ds / 6.0 * (1.0 + 2.0 * t +
          c - ds / 20.0 * (5.0 - 2.0 * c + 28.0 * t - 3.0 * cs + 8.0 * esp +
          24.0 * ts))) / cos_phi);
   }
else
   {
   *lat = HALF_PI * sign(y);
   *lon = lon_center;
   }
/* -- end of original GCTP utminv -- */

/* Convert lat/long to degrees */
*lat *= R2D;
*lon *= R2D;

return(0);
}


/* Universal Transverse Mercator inverse equations--mapping lat,long to x,y
   to line,sample
  -------------------------------------------------------------------------*/
int LSutmfor(double *s, double *l, double lon, double lat)
//double *x;		/* (O) X projection coordinate 		*/
//double *y;		/* (O) Y projection coordinate 		*/
//double lon;		/* (I) Longitude 				*/
//double lat;		/* (I) Latitude 				*/

/* this is meant to be visible from other routines */
{
double x;		/* X projection coordinate */
double y;		/* Y projection coordinate */

double delta_lon;	/* Delta longitude (Given longitude - center 	*/
double sin_phi, cos_phi;/* sin and cos value				*/
double al, als;		/* temporary values				*/
double b;		/* temporary values				*/
double c, t, tq;	/* temporary values				*/
double con, n, ml;	/* cone constant, small m			*/
double dy, dx;

lat *= D2R;
lon *= D2R;

/* -- start of original GCTP utminv -- */
/* Forward equations
  -----------------*/
delta_lon = lon - lon_center;
sincos(lat, &sin_phi, &cos_phi);

/* This part was in the fortran code and is for the spherical form 
----------------------------------------------------------------*/
if (ind != 0)
  {
  b = cos_phi * sin(delta_lon);
  if ((fabs(fabs(b) - 1.0)) < .0000000001)
     {
     printf("Point projects into infinity, tm-for\n");
     return(93);
     }
  else
     {
     *s = .5 * r_major * scale_factor * log((1.0 + b)/(1.0 - b));
     con = acos(cos_phi * cos(delta_lon)/sqrt(1.0 - b*b));
     if (lat < 0)
        con = - con;
     *l = r_major * scale_factor * (con - lat_origin); 
     return(0);
     }
  }

al  = cos_phi * delta_lon;
als = SQUARE(al);
c   = esp * SQUARE(cos_phi);
tq  = tan(lat);
t   = SQUARE(tq);
con = 1.0 - es * SQUARE(sin_phi);
n   = r_major / sqrt(con);
ml  = r_major * mlfn(e0, e1, e2, e3, lat);

x  = scale_factor * n * al * (1.0 + als / 6.0 * (1.0 - t + c + als / 20.0 *
      (5.0 - 18.0 * t + SQUARE(t) + 72.0 * c - 58.0 * esp))) + false_easting;

y  = scale_factor * (ml - ml0 + n * tq * (als * (0.5 + als / 24.0 *
      (5.0 - t + 9.0 * c + 4.0 * SQUARE(c) + als / 30.0 * (61.0 - 58.0 * t
      + SQUARE(t) + 600.0 * c - 330.0 * esp))))) + false_northing;
/* -- end of original GCTP utminv -- */

x -= ul_corner[0];
y -= ul_corner[1];

dx = ( (x) * cos_orien) + ( (y) * sin_orien);
dy = ( (x) * sin_orien) - ( (y) * cos_orien);

*s = dx/pixel_size - 0.5;
*l = dy/pixel_size - 0.5;

return(0);
}


/* Polar Stereographic inverse equations--mapping line,sample to x,y to
   lat/long
  ---------------------------------------------------------------------*/
int LSpsinv(double s, double l, double *lon, double *lat)
//    double x,			/* (O) X projection coordinate 	*/
//    double y,			/* (O) Y projection coordinate 	*/
//    double *lon,			/* (I) Longitude 		*/
//    double *lat			/* (I) Latitude 		*/

/* this is meant to be visible from other routines */

{
double x, y;
double rh;			/* height above ellipsiod	*/
double ts;			/* small value t		*/
double temp;			/* temporary variable		*/
long   flag;			/* error flag			*/
double dl, dp, dy, dx;

dl = (l + 0.5) * pixel_size;
dp = (s + 0.5) * pixel_size;
dy = (dp * sin_orien) - (dl * cos_orien);
dx = (dp * cos_orien) + (dl * sin_orien);
y = ul_corner[1] + dy;
x = ul_corner[0] + dx;

/* -- start of original GCTP psinv --*/
flag = 0;
x = (x - false_easting) * fac;
y = (y - false_northing) *fac;
rh = sqrt(x * x + y * y);
if (ind != 0)
  ts = rh * tcs/(r_major * mcs);
else
  ts = rh * e4 / (r_major * 2.0);
*lat = fac * phi2z(e,ts,&flag);
if (flag != 0)
   return(flag);
if (rh == 0)
   *lon = fac * center_lon;
else
   {
   temp = atan2(x, -y);
   *lon = adjust_lon(fac *temp + center_lon);
   }
/* -- end of original GCTP psinv --*/

/* Convert lat/long to degrees */
*lat *= R2D;
*lon *= R2D;

return(0);
}

/* Polar Stereographic forward equations--mapping lat,long to x,y to
   line,sample
  ------------------------------------------------------------------*/
int LSpsfor(double *s, double *l, double lon, double lat)
//double *x;		/* (O) X projection coordinate 		*/
//double *y;		/* (O) Y projection coordinate 		*/
//double lon;		/* (I) Longitude 				*/
//double lat;		/* (I) Latitude 				*/

/* this is meant to be visible from other routines */

{
double x;		/* X projection coordinate */
double y;		/* Y projection coordinate */
double con1;			/* adjusted longitude		*/
double con2;			/* adjusted latitude		*/
double rh;			/* height above ellipsoid	*/
double sinphi;			/* sin value			*/
double ts;			/* value of small t		*/
double dy, dx;

lat *= D2R;
lon *= D2R;

/* -- start of original GCTP psfor --*/
con1 = fac * adjust_lon(lon - center_lon);
con2 = fac * lat;
sinphi = sin(con2);
ts = tsfnz(e,con2,sinphi);
if (ind != 0)
   rh = r_major * mcs * ts / tcs;
else
   rh = 2.0 * r_major * ts / e4;
x = fac * rh * sin(con1) + false_easting;
y = -fac * rh * cos(con1) + false_northing;;
/* -- end of original GCTP psfor --*/

x -= ul_corner[0];
y -= ul_corner[1];

dx = ( (x) * cos_orien) + ( (y) * sin_orien);
dy = ( (x) * sin_orien) - ( (y) * cos_orien);

*s = dx/pixel_size - 0.5;
*l = dy/pixel_size - 0.5;

return(0);
}
