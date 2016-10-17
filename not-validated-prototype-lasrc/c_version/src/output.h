#ifndef OUTPUT_H
#define OUTPUT_H

#include "common.h"
#include "input.h"

/* Define some of the constants to use in the output data products */
#define FILL_VALUE -9999
#define CLOUD_FILL_VALUE 0
#define SCALE_FACTOR 0.0001
#define MULT_FACTOR 10000.0
#define SCALE_FACTOR_TH 0.1
#define MULT_FACTOR_TH 10.0
#define MIN_VALID -2000
#define MAX_VALID 16000
#define MIN_VALID_TH 1500
#define MAX_VALID_TH 3500

/* Structure for the 'output' data type */
typedef struct {
  bool open;            /* Flag to indicate whether output file is open;
                           'true' = open, 'false' = not open */
  bool process_collection; /* Is this scene a Collection product? */
  Inst_t inst;          /* instrument */
  int nband;            /* Number of output bands */
  int nlines;           /* Number of output lines */
  int nsamps;           /* Number of output samples */
  Espa_internal_meta_t metadata;  /* Metadata container to hold the band
                           metadata for the output bands; global metadata
                           won't be valid */
  FILE *fp_bin[NBAND_TTL_OUT];  /* File pointer for binary files; see common.h
                           for the bands and order of bands in the output */
} Output_t;

/* Prototypes */
Output_t *open_output
(
    Espa_internal_meta_t *in_meta,  /* I: input metadata structure */
    Input_t *input,                 /* I: input band data structure */
    bool toa,                       /* I: set this structure up for the TOA
                                          bands vs. the SR bands */
    bool process_collection         /* I: should this scene be processed as a
                                          collection product, which affects
                                          the output of QA bands */
);

int close_output
(
    Output_t *this,   /* I/O: Output data structure to close */
    bool toa          /* I: output structure is for TOA bands vs. SR bands */
);

int free_output
(
    Output_t *this    /* I/O: Output data structure to free */
);

int put_output_lines
(
    Output_t *this,    /* I: Output data structure; buf contains the line to
                             be written */
    void *buf,         /* I: buffer to be written */
    int iband,         /* I: current band to be written (0-based) */
    int iline,         /* I: current line to be written (0-based) */
    int nlines,        /* I: number of lines to be written */
    int nbytes         /* I: number of bytes per pixel in this band */
);

int get_output_lines
(
    Output_t *this,  /* I: pointer to output data structure */
    int iband,       /* I: current band to read (0-based) */
    int iline,       /* I: current line to read (0-based) */
    int nlines,      /* I: number of lines to read */
    int nbytes,      /* I: number of bytes per pixel in this band */
    void *buf        /* I: pointer to the buffer to be returned */
);

char *upper_case_str
(
    char *str    /* I: string to be converted to upper case */
);

#endif
