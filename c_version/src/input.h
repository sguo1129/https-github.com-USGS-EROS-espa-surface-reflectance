#ifndef INPUT_H
#define INPUT_H

#include <stdlib.h>
#include <stdio.h>
#include "input.h"
#include "date.h"
#include "common.h"
#include "espa_metadata.h"
#include "error_handler.h"
#include "raw_binary_io.h"

#define INPUT_FILL (0)
#define ANGLE_FILL (-999.0)
#define WRS_FILL (-1)
#define GAIN_BIAS_FILL (-999.0)

/* Structure for the input metadata */
typedef struct {
    Sat_t sat;               /* satellite */
    Inst_t inst;             /* instrument */
    Date_t acq_date;         /* acqsition date/time (scene center) */
    bool time_fill;          /* acqsition time fill; true = fill value (0h) */
    Date_t prod_date;        /* production date */
    float sun_zen;           /* solar zenith angle (degrees; scene center) */
    float sun_az;            /* solar azimuth angle (degrees; scene center) */
    Wrs_t wrs_sys;           /* WRS system */
    int ipath;               /* WRS path number */
    int irow;                /* WRS row number */
    uint16 fill;             /* fill value */
    int iband[NBAND_REFL_MAX];     /* reflectance band numbers */
    int iband_th[NBAND_THM_MAX];   /* thermal band numbers */
    int iband_pan[NBAND_PAN_MAX];  /* pan band numbers */
    int iband_qa[NBAND_QA_MAX];    /* QA band numbers */
    int iband_lw;                  /* land/water mask band number */
    bool gain_set;                 /* are the gains and biases set? */
    float gain[NBAND_REFL_MAX];    /* reflectance band gain */
    float gain_th[NBAND_THM_MAX];  /* thermal band gain */
    float gain_pan[NBAND_PAN_MAX]; /* pan band gain */
    float bias[NBAND_REFL_MAX];    /* reflectance band bias */
    float bias_th[NBAND_THM_MAX];  /* thermal band bias */
    float bias_pan[NBAND_PAN_MAX]; /* pan band bias */
} Input_meta_t;

/* Structure for the input data */
typedef struct {
    Input_meta_t meta;        /* input metadata */
    int nband;                /* number of reflectance bands */
    int nband_th;             /* number of thermal bands */
    int nband_pan;            /* number of pan bands */
    int nband_qa;             /* number of QA bands */
    Img_coord_info_t size;    /* input file size */
    Img_coord_info_t size_th; /* input thermal file size */
    Img_coord_info_t size_pan;/* input pan file size */
    Img_coord_info_t size_qa; /* input QA file size */
    Img_coord_info_t size_lw; /* input land/water mask file size */
    float scale_factor;       /* scale factor for reflectance bands */
    float scale_factor_th;    /* scale factor for thermal bands */
    float scale_factor_pan;   /* scale factor for pan bands */
    char *file_name[NBAND_REFL_MAX];  /* name of the input reflectance files */
    char *file_name_th[NBAND_THM_MAX];   /* name of the input thermal files */
    char *file_name_pan[NBAND_PAN_MAX];  /* name of the input pan files */
    char *file_name_qa[NBAND_QA_MAX];    /* name of the input QA files */
    char *file_name_lw;        /* name of the input land/water mask file */
    bool open[NBAND_REFL_MAX]; /* flag to indicate whether the specific input
                                  file is open for access; 'true' = open, 
                                  'false' = not open */
    bool open_th[NBAND_THM_MAX];  /* thermal band open flag */
    bool open_pan[NBAND_PAN_MAX]; /* pan band open flag */
    bool open_qa[NBAND_QA_MAX];   /* QA band open flag */
    bool open_lw;                 /* land/water band open flag */
    FILE *fp_bin[NBAND_REFL_MAX]; /* pointer for reflectance binary files */
    FILE *fp_bin_th[NBAND_THM_MAX];  /* pointer for thermal binary files */
    FILE *fp_bin_pan[NBAND_PAN_MAX]; /* pointer for pan binary files */
    FILE *fp_bin_qa[NBAND_QA_MAX];   /* pointer for QA binary files */
    FILE *fp_bin_lw;                 /* pointer for land/water binary file */
} Input_t;

/* Prototypes */
Input_t *open_input
(
    Espa_internal_meta_t *metadata      /* I: input metadata */
);

void close_input
(
    Input_t *this    /* I: pointer to input data structure */
);

void free_input
(
    Input_t *this    /* I: pointer to input data structure */
);

int get_input_refl_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current refl band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint16 *out_arr  /* O: output array to populate */
);

int get_input_th_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current thermal band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint16 *out_arr  /* O: output array to populate */
);

int get_input_pan_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current pan band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint16 *out_arr  /* O: output array to populate */
);

int get_input_qa_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iband,       /* I: current QA band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint16 *out_arr  /* O: output array to populate */
);

int get_input_lw_lines
(
    Input_t *this,   /* I: pointer to input data structure */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    uint8 *out_arr   /* O: output array to populate */
);

int get_xml_input
(
    Espa_internal_meta_t *metadata,  /* I: XML metadata */
    Input_t *this                    /* O: data structure for the input file */
);


#endif
