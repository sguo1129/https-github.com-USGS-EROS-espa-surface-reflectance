#ifndef _COMMON_H_
#define _COMMON_H_

#include "hdf.h"
#include "mfhdf.h"
typedef char byte;

/* Surface reflectance version */
#define SR_VERSION "0.10.0"

/* How many lines of data should be processed at one time */
#define PROC_NLINES 1000

/* For angle conversions -
   degrees to radians = PI/180
   radians to degrees = 180/PI */
#define DEG2RAD 0.017453293
#define RAD2DEG 57.29577951

/* For divisions - to minimize processing time */
#define ONE_DIV_1013 0.000987166
#define ONE_DIV_8500 0.000117647

/* Number of bands corrected to surface reflectance (bands 1-7).  The
   atmospheric correction variables store information for 8 bands, so we will
   go with that for the array size. */
#define NSR_BANDS 8

/* L8 Level-1 products have 8 reflectance bands (bands 1-7, and 9),
   2 thermal bands (band 10 and 11), 1 pan band (band 8), and 1 QA band
   (band 12) */
#define NBAND_REFL_MAX 8
#define NBAND_THM_MAX 2
#define NBAND_PAN_MAX 1
#define NBAND_QA_MAX 1
#define NBAND_TTL_MAX (NBAND_REFL_MAX + NBAND_THM_MAX + NBAND_PAN_MAX + NBAND_QA_MAX)

/* L8 surface reflectance products have 8 reflectance bands, 2 thermal bands, 
   0 pan bands, and 2 QA bands (pre-collection has 2 and collection has 1) */
#define NBAND_REFL_OUT 8
#define NBAND_THM_OUT 2
#define NBAND_PAN_OUT 0
#define NBAND_QA_OUT 2
#define NBAND_TTL_OUT (NBAND_REFL_OUT + NBAND_THM_OUT + NBAND_PAN_OUT + NBAND_QA_OUT)

/* CMG and DEM files are lat/long images where each pixel represents 0.05 deg x
   0.05 deg */
/* DEM information */
#define DEM_NBLAT 3600
#define DEM_NBLON 7200

/* Ratio file information */
#define RATIO_NBLAT 3600
#define RATIO_NBLON 7200

/* Ozone and water vapor information */
#define CMG_NBLAT 3600
#define CMG_NBLON 7200

/* Define the input products to be processed.  NOTE: DN_TTL should be the same
   as NBAND_TTL_MAX. */
typedef enum {DN_BAND1=0, DN_BAND2, DN_BAND3, DN_BAND4, DN_BAND5, DN_BAND6,
    DN_BAND7, DN_BAND8, DN_BAND9, DN_BAND10, DN_BAND11, DN_QA, DN_TTL}
    Mydn_band_t;

/* Define the output products to be processed. NOTE: SR_TTL should be the same
   as NBAND_TTL_OUT.  Pre-collection products use the cloud and ipflag bands,
   but collection products have a single QA band.  The cloud QA contains all
   the QA info and the ipflag band does not exist. */
typedef enum {SR_BAND1=0, SR_BAND2, SR_BAND3, SR_BAND4, SR_BAND5, SR_BAND6,
    SR_BAND7, SR_BAND9, SR_BAND10, SR_BAND11, SR_CLOUD, SR_IPFLAG, SR_TTL}
    Mysr_band_t;

/* Bit location of weight for cloudmask QA. Bit 4 is used as a temporary
   bit location as well as the first aerosol bit. */
typedef enum {
  CIR_QA=0,      /* cirrus cloud bit            = 1 */
  CLD_QA=1,      /* cloud bit                   = 2 */
  CLDA_QA=2,     /* adjacent cloud bit          = 4 */
  CLDS_QA=3,     /* cloud shadow bit            = 8 */
  CLDT_QA=4,     /* temporary cloud shadow bit, = 16
                    used and cleared before the aerosol is masked for the
                    actual QA (can't use the cloud shadow bit in this case) */
  AERO1_QA=4,    /* these two AERO bits mark the amount of aerosols and = 16 */
  AERO2_QA=5,    /* reflect the level of atmospheric correction made    = 32 */
  SNOW_QA=6,     /* snow bit = 64 (reserved for later) */
  WAT_QA=7,      /* water bit = 128 (deep water set via DEM and internal water
                    flag) */
/*#------------------ Used for Collection Products ------------------------#*/
  AERO_RETRIEVAL_NDVI_QA=8,   /* aerosol retrieval failed due to the NDVI test,
                                 bit = 256 */
  AERO_RETRIEVAL_RESID_QA=9,  /* aerosol retrieval failed due to the model
                                 residual test, bit = 512 */
  AERO_INTERP_QA=10           /* aerosol was interpolated bit = 1024 */
} Cloudqa_t;

/* Class values of ipflag (interpolation flag) QA */
typedef enum {
  IPFLAG_CLEAR=0,          /* IPFLAG is clear */
  IPFLAG_INTERP=1,         /* aerosol was interpolated */
  IPFLAG_NDVI_FAIL=2,      /* NDVI test failed */
  IPFLAG_RESIDUAL_FAIL=3,  /* residual test failed */
  IPFLAG_FILL=5,           /* fill value */
  IPFLAG_WATER=6           /* IPFLAG indicates water */
} Ipflag_t;

/* Satellite type definitions, mainly to allow future satellites to be
   supported if needed */
typedef enum {
  SAT_NULL = -1,
  SAT_LANDSAT_8 = 0, 
  SAT_MAX
} Sat_t;

/* Instrument type definition */
typedef enum {
  INST_NULL = -1,
  INST_OLI_TIRS = 0, 
  INST_OLI, 
  INST_MAX
} Inst_t;

/* World Reference System (WRS) type definition */
typedef enum {
  WRS_NULL = -1,
  WRS_1 = 0, 
  WRS_2,
  WRS_MAX
} Wrs_t;

typedef struct {
  int nlines;
  int nsamps;
  double pixsize[2];
} Img_coord_info_t;

#endif
