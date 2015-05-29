#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "lndcal.h"
#include "keyvalue.h"
#include "const.h"
#include "param.h"
#include "input.h"
#include "lut.h"
#include "output.h"
#include "cal.h"
#include "bool.h"
#include "error.h"
#include "util.h"

#include <time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

/* Type definitions */

#define NSDS NBAND_CAL_MAX

/* Functions */
/* !Revision:
 *
 * revision 1.0.0 9/12/2012  Gail Schmidt, USGS
 * - modified the application to write the thermal band QA bits to the
 *   output thermal product
 *
 * revision 1.0.1 2/20/2013  Gail Schmidt, USGS
 * - modified the application to calculate and write the bounding coords
 *
 * revision 1.1.2 3/22/2013  Gail Schmidt, USGS
 * - modified the application to write the UL and LR lat/long to the metadata
 *
 * revision 2.0.0 1/21/2014  Gail Schmidt, USGS
 * - modified application to use the ESPA internal raw binary file format
 * - removed any recalibration-related or DN map related code
 *
 * revision 2.0.1 8/5/2014  Gail Schmidt, USGS/EROS
 * - obtain the location of the ESPA schema file from an environment variable
 *   vs. the ESPA http site
 */

int main (int argc, const char **argv) {
  Param_t *param = NULL;
  Input_t *input = NULL;
  Lut_t *lut = NULL;
  Output_t *output = NULL;
  Output_t *output_th = NULL;
  int iline, isamp,oline, ib, jb, iz, val;
  unsigned char *line_in = NULL;
  unsigned char *line_in_thz = NULL;
  unsigned char *line_out_qa = NULL;
  int16 *line_out = NULL;
  int16 *line_out_th = NULL;
  int16 *line_out_thz = NULL;
  Cal_stats_t cal_stats;
  Cal_stats6_t cal_stats6;
  int nps,nls, nps6, nls6;
  int zoomx, zoomy;
  int i,odometer_flag=0;
  char msgbuf[1024];
  char envi_file[STR_SIZE]; /* name of the output ENVI header file */
  char *cptr=NULL;          /* pointer to the file extension */
  size_t input_psize;
  int qa_band = QA_BAND_NUM;
  int nband_refl = NBAND_REFL_MAX;
  int ifill, num_zero;
  int maxth=0;
  int mss_flag=0;
  Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
  Envi_header_t envi_hdr;   /* output ENVI header information */

  printf ("\nRunning lndcal ...\n");
  for (i=1; i<argc; i++)if ( !strcmp(argv[i],"-o") )odometer_flag=1;

  /* Read the parameters from the input parameter file */
  param = GetParam(argc, argv);
  if (param == (Param_t *)NULL) EXIT_ERROR("getting runtime parameters",
    "main");

  /* Validate the input metadata file */
  if (validate_xml_file (param->input_xml_file_name) != SUCCESS)
  {  /* Error messages already written */
      EXIT_ERROR("Failure validating XML file", "main");
  }

  /* Initialize the metadata structure */
  init_metadata_struct (&xml_metadata);

  /* Parse the metadata file into our internal metadata structure; also
     allocates space as needed for various pointers in the global and band
     metadata */
  if (parse_metadata (param->input_xml_file_name, &xml_metadata) != SUCCESS)
  {  /* Error messages already written */
    EXIT_ERROR("parsing XML file", "main");
  }

  /* Check to see if the gain and bias values were specified */
  if (!existRadGB (&xml_metadata))
    EXIT_ERROR("Gains and biases don't exist in XML file (TOA radiance gain "
      "and bias fields) for each band.  Make sure to utilize the latest LPGS "
      "MTL file for conversion to the ESPA internal raw binary format as the "
      "gains and biases should be in that file.", "main");
  
  /* Open input file */
  input = OpenInput (&xml_metadata);
  if (input == (Input_t *)NULL)
    EXIT_ERROR("setting up input from XML structure", "main");

  /* Get Lookup table */
  lut = GetLut(param, input->nband, input);
  if (lut == (Lut_t *)NULL) EXIT_ERROR("bad lut file", "main");

  nps6=  input->size_th.s;
  nls6=  input->size_th.l;
  nps =  input->size.s;
  nls =  input->size.l;
  zoomx= nint( (float)nps / (float)nps6 );
  zoomy= nint( (float)nls / (float)nls6 );

  for (ib = 0; ib < input->nband; ib++) 
    cal_stats.first[ib] = true;
  cal_stats6.first = true;
  if (input->meta.inst == INST_MSS)mss_flag=1; 

  /* Open the output files.  Raw binary band files will be be opened. */
  output = OpenOutput(&xml_metadata, input, param, lut, false /*not thermal*/,
    mss_flag);
  if (output == NULL) EXIT_ERROR("opening output file", "main");

  /* Allocate memory for the input buffer, enough for all reflectance bands */
  input_psize = sizeof(unsigned char);
  line_in = calloc (input->size.s * nband_refl, input_psize);
   if (line_in == NULL) 
     EXIT_ERROR("allocating input line buffer", "main");

  /* Create and open output thermal band, if one exists */
  if ( input->nband_th > 0 ) {
    output_th = OpenOutput (&xml_metadata, input, param, lut, true /*thermal*/,
      mss_flag);
    if (output_th == NULL)
      EXIT_ERROR("opening output therm file", "main");

    /* Allocate memory for the thermal input and output buffer, only holds
       one band */
    line_out_th = calloc(input->size_th.s, sizeof(int16));
    if (line_out_th == NULL) 
      EXIT_ERROR("allocating thermal output line buffer", "main");

    if (zoomx == 1) {
      line_out_thz = line_out_th;
      line_in_thz = line_in;
    }
    else {
      line_out_thz = calloc (input->size.s, sizeof(int16));
      if (line_out_thz == NULL) 
        EXIT_ERROR("allocating thermal zoom output line buffer", "main");
      line_in_thz = calloc (input->size.s, input_psize);
      if (line_in_thz == NULL) 
        EXIT_ERROR("allocating thermal zoom input line buffer", "main");
    }
  } else {
    printf("*** no output thermal file ***\n"); 
  }

  /* Allocate memory for output lines for both the image and QA data */
  line_out = calloc (input->size.s, sizeof (int16));
  if (line_out == NULL) 
    EXIT_ERROR("allocating output line buffer", "main");

  line_out_qa = calloc (input->size.s, sizeof(unsigned char));
  if (line_out_qa == NULL) 
    EXIT_ERROR("allocating qa output line buffer", "main");
  memset (line_out_qa, 0, input->size.s * sizeof(unsigned char));    

  /* Do for each THERMAL line */
  oline= 0;
  if (input->nband_th > 0) {
    ifill= (int)lut->in_fill;
    for (iline = 0; iline < input->size_th.l; iline++) {
      ib=0;
      if (!GetInputLineTh(input, iline, line_in))
        EXIT_ERROR("reading input data for a line", "main");

      if ( odometer_flag && ( iline==0 || iline ==(nls-1) || iline%100==0  ) ){ 
        if ( zoomy == 1 )
          printf("--- main loop BAND6 Line %d --- \r",iline); 
        else
          printf("--- main loop BAND6 Line in=%d out=%d --- \r",iline,oline); 
        fflush(stdout); 
      }

      memset(line_out_qa, 0, input->size.s*sizeof(unsigned char));    
      if (!Cal6(lut, input, line_in, line_out_th, line_out_qa, &cal_stats6,
        iline))
        EXIT_ERROR("doing calibration for a line", "main");

      if ( zoomx>1 ) {
        zoomIt(line_out_thz, line_out_th, nps/zoomx, zoomx );
        zoomIt8(line_in_thz, line_in, nps/zoomx, zoomx );
      }

      for ( iz=0; iz<zoomy; iz++ ) {
        for (isamp = 0; isamp < input->size.s; isamp++) {
          val= getValue(line_in_thz, isamp);
          if ( val> maxth) maxth=val;
          if ( val==ifill) line_out_qa[isamp] = lut->qa_fill; 
          else if ( val>=SATU_VAL6 ) line_out_qa[isamp] = ( 0x000001 << 6 ); 
        }

        if ( oline<nls ) {
          if (!PutOutputLine(output_th, ib, oline, line_out_thz)) {
            sprintf(msgbuf,"write thermal error ib=%d oline=%d iline=%d",ib,
              oline,iline);
            EXIT_ERROR(msgbuf, "main");
          }

          if (input->meta.inst != INST_MSS) 
            if (!PutOutputLine(output_th, ib+1, oline, line_out_qa)) {
	          sprintf(msgbuf,"write thermal QA error ib=%d oline=%d iline=%d",
                ib+1,oline,iline);
              EXIT_ERROR(msgbuf, "main");
            }
        }
        oline++;
      }
    } /* end loop for each thermal line */
  }
  if (odometer_flag) printf("\n");

  if (input->nband_th > 0)
    if (!CloseOutput(output_th))
      EXIT_ERROR("closing output thermal file", "main");

  /* Do for each REFLECTIVE line */
  ifill= (int)lut->in_fill;
  for (iline = 0; iline < input->size.l; iline++){
    /* Do for each band */

    if ( odometer_flag && ( iline==0 || iline ==(nls-1) || iline%100==0  ) )
     {printf("--- main reflective loop Line %d ---\r",iline); fflush(stdout);}

    memset(line_out_qa, 0, input->size.s*sizeof(unsigned char));
    
    for (ib = 0; ib < input->nband; ib++) {
      if (!GetInputLine(input, ib, iline, &line_in[ib*nps]))
        EXIT_ERROR("reading input data for a line", "main");
    }
    
    for (isamp = 0; isamp < input->size.s; isamp++){
      num_zero=0;
      for (ib = 0; ib < input->nband; ib++) {
        jb= (ib != 5 ) ? ib+1 : ib+2;
        val= getValue((unsigned char *)&line_in[ib*nps], isamp);
	    if ( val==ifill   )num_zero++;
        if ( val==SATU_VAL[ib] ) line_out_qa[isamp]|= ( 0x000001 <<jb ); 
      }
      /* Feng fixed bug by changing "|=" to "=" below (4/17/09) */
      if ( num_zero >  0 )line_out_qa[isamp] = lut->qa_fill; 
    }

    for (ib = 0; ib < input->nband; ib++) {
      if (!Cal(lut, ib, input, &line_in[ib*nps], line_out, line_out_qa,
        &cal_stats,iline))
        EXIT_ERROR("doing calibraton for a line", "main");

      if (!PutOutputLine(output, ib, iline, line_out))
        EXIT_ERROR("reading input data for a line", "main");
    } /* End loop for each band */
        
    if (input->meta.inst != INST_MSS) 
      if (!PutOutputLine(output, qa_band, iline, line_out_qa))
        EXIT_ERROR("writing qa data for a line", "main");
  } /* End loop for each line */

  if ( odometer_flag )printf("\n");

  for (ib = 0; ib < input->nband; ib++) {
    printf(
      " band %d rad min %8.5g max %8.4f  |  ref min  %8.5f max  %8.4f\n", 
      input->meta.iband[ib], cal_stats.rad_min[ib], cal_stats.rad_max[ib],
      cal_stats.ref_min[ib], cal_stats.ref_max[ib]);
  }

  if ( input->nband_th > 0 )
    printf(
      " band %d rad min %8.5g max %8.4f  |  tmp min  %8.5f max  %8.4f\n", 6,
      cal_stats6.rad_min,  cal_stats6.rad_max,
      cal_stats6.temp_min, cal_stats6.temp_max);

  /* Close input and output files */
  if (!CloseInput(input)) EXIT_ERROR("closing input file", "main");
  if (!CloseOutput(output)) EXIT_ERROR("closing input file", "main");

  /* Write the ENVI header for reflectance files */
  for (ib = 0; ib < output->nband; ib++) {
    /* Create the ENVI header file this band */
    if (create_envi_struct (&output->metadata.band[ib], &xml_metadata.global,
      &envi_hdr) != SUCCESS)
        EXIT_ERROR("Creating the ENVI header structure for this file.", "main");

    /* Write the ENVI header */
    strcpy (envi_file, output->metadata.band[ib].file_name);
    cptr = strchr (envi_file, '.');
    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
        EXIT_ERROR("Writing the ENVI header file.", "main");
  }

  /* Write the ENVI header for thermal files */
  for (ib = 0; ib < output_th->nband; ib++) {
    /* Create the ENVI header file this band */
    if (create_envi_struct (&output_th->metadata.band[ib], &xml_metadata.global,
      &envi_hdr) != SUCCESS)
        EXIT_ERROR("Creating the ENVI header structure for this file.", "main");

    /* Write the ENVI header */
    strcpy (envi_file, output_th->metadata.band[ib].file_name);
    cptr = strchr (envi_file, '.');
    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
        EXIT_ERROR("Writing the ENVI header file.", "main");
  }

  /* Append the reflective and thermal bands to the XML file */
  if (append_metadata (output->nband, output->metadata.band,
    param->input_xml_file_name) != SUCCESS)
    EXIT_ERROR("appending reflectance and QA bands", "main");
  if (input->nband_th > 0) {
    if (append_metadata (output_th->nband, output_th->metadata.band,
      param->input_xml_file_name) != SUCCESS)
      EXIT_ERROR("appending thermal and QA bands", "main");
  }

  /* Free the metadata structure */
  free_metadata (&xml_metadata);

  /* Free memory */
  if (!FreeParam(param)) 
    EXIT_ERROR("freeing parameter stucture", "main");

  if (!FreeInput(input)) 
    EXIT_ERROR("freeing input file stucture", "main");

  if (!FreeLut(lut)) 
    EXIT_ERROR("freeing lut file stucture", "main");

  if (!FreeOutput(output)) 
    EXIT_ERROR("freeing output file stucture", "main");

  free(line_out);
  line_out = NULL;
  free(line_in);
  line_in = NULL;
  free(line_out_qa);
  line_out_qa = NULL;
  free(line_out_th);
  line_out_th = NULL;
  if (zoomx != 1) {
    free(line_in_thz);
    free(line_out_thz);
  }
  line_in_thz = NULL;
  line_out_thz = NULL;

  /* All done */
  printf ("lndcal complete.\n");
  return (EXIT_SUCCESS);
}
