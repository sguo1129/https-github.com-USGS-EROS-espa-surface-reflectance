/*****************************************************************************
FILE: lndpm.c
  
PURPOSE: Contains functions for reading the input XML metadata file and
creating the parameter and metadata files needed for downstream processing
by LEDAPS applications.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/16/2014    Gail Schmidt     Modified to use ESPA internal file format.
                              The *.metadata.txt file is no longer needed by
                              lndcal nor is the MTL file and instead the XML
                              file is passed to both lndcal and lndsr.
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.
1/30/2014    Gail Schmidt     Reflective and thermal output files are no
                              longer needed as parameters in the lndcal*.txt
                              file, since outputs are handled in the XML file.
1/30/2014    Gail Schmidt     Removed any recalibration-related or DN map
                              related code
2/3/2014     Gail Schmidt     Reflective, thermal input and surface reflectance
                              output files are no longer needed as parameters
                              in the lndsr*.txt file, since inputs and outputs

NOTES:
  1. The XML metadata format written via this library follows the ESPA internal
     metadata format found in ESPA Raw Binary Format v1.0.doc.  The schema for
     the ESPA internal metadata format is available at
     http://espa.cr.usgs.gov/static/schema/espa_internal_metadata_v1_0.xsd.
*****************************************************************************/
#include <sys/stat.h>
#include "lndpm.h"

FILE *LOG_FP = NULL;     /* pointer to the log file */
int open_log (char name[]);
void LOG (char *fmt, ...);
int conv_date (int *mm, int *dd, int yyyy);
int find_file(char *path, char *name);

int main (int argc, char *argv[])
{
    char FUNC_NAME[] = "lndpm";    /* function name */
    char errmsg[STR_SIZE];         /* error message */
    char input_xml[STR_SIZE];      /* name of the input XML file */
    char tmpfile[STR_SIZE];        /* temporary storage of the XML file */
    char scene_name[STR_SIZE];     /* name of the scene */
    char lndcal_name[STR_SIZE];    /* name of the lndcal input file */
    char lndsr_name[STR_SIZE];     /* name of the lndsr input file */
    char dem[STR_SIZE];            /* name of DEM file */
    char ozone[STR_SIZE];          /* name of ozone file */
    char reanalysis[STR_SIZE];     /* name of NCEP file */
    char path_buf[DIR_BUF_SIZE];   /* path to the auxillary/cal file */
    char *anc_path = NULL;         /* path for LEDAPS ancillary data */
    char *token_ptr = NULL;        /* pointer used for obtaining scene name */
    char *file_ptr = NULL;         /* pointer used for obtaining file name */
    int year, month, day;          /* year, month, day of acquisition date */
    bool anc_missing = false;      /* is the ancillary data missing? */
    FILE *out = NULL;              /* pointer to the output parameter file */
    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */

    /* Open the log file */
    printf ("\nRunning lndpm ...\n");
    if (open_log ("lndpm : Landsat Metadata Parser") == ERROR)
    {
        sprintf (errmsg, "Create log report file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Check the command-line arguments and get the name of the XML file */
    if (argc != 2)
    {
        sprintf (errmsg, "Usage: lndpm <input_xml_file>");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    strcpy (input_xml, argv[1]);

    /* Validate the input metadata file */
    if (validate_xml_file (input_xml) != SUCCESS)
    {  /* Error messages already written */
        return (ERROR);
    }

    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (input_xml, &xml_metadata) != SUCCESS)
    {  /* Error messages already written */
        return (ERROR);
    }

    /* Get the scene name.  Strip off the path and file extension as it's
       assumed the XML filename is the scene name followed by the extension. */
    strcpy (tmpfile, input_xml);
    file_ptr = strrchr (tmpfile, '/');
    if (file_ptr != NULL)
        file_ptr++;
    else
        file_ptr = tmpfile;
    token_ptr = strtok (file_ptr, ".");
    sprintf (scene_name, "%s", token_ptr);
  
    /* Set up the names of the input files for downstream processing */
    sprintf (lndcal_name, "lndcal.%s.txt", scene_name);
    sprintf (lndsr_name, "lndsr.%s.txt", scene_name);

    /* Open the parameter file for lndcal for writing */
    out = fopen (lndcal_name, "w");
    if (out == NULL)
    {
        sprintf (errmsg, "Opening lndcal parameter file for writing: %s",
            lndcal_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Get the path for the LEDAPS auxillary products (NCEP, TOMS, DEM, etc.)
       from the ANC_PATH environment variable.  If it isn't defined, then
       assume the products are in the local directory. */
    anc_path = getenv ("ANC_PATH");
    if (anc_path == NULL)
    {
        anc_path = ".";
        sprintf (errmsg, "ANC_PATH environment variable isn't defined. It is "
            "assumed the LEDAPS ancillary products will be available from "
            "the local directory.");
        error_handler (false, FUNC_NAME, errmsg);
    }

    /* Write the parameter data to the lndcal parameter file */
    fprintf (out, "PARAMETER_FILE\n");
    fprintf (out, "XML_FILE = %s\n", input_xml);
    fprintf (out, "LEDAPSVersion = %s\n", LEDAPS_VERSION);
    fprintf (out, "END\n");
    fclose (out);

    /* Get year, month, and day from acquisition date */
    token_ptr = strtok (xml_metadata.global.acquisition_date, "-");
    sscanf (token_ptr, "%d", &year);
    token_ptr = strtok (NULL, "-");
    sscanf (token_ptr, "%d", &month);
    token_ptr = strtok (NULL, "-");
    sscanf (token_ptr, "%d", &day);
  
    /* Convert to day of year */
    if (conv_date (&month, &day, year) != SUCCESS)
    {
        sprintf (errmsg, "Not able to convert the month, day, year from the "
            "acquisition date to Julian DOY.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Find and prepare auxillary files */
    /* DEM file */
    strcpy (dem, "CMGDEM.hdf");
    strcpy (path_buf, anc_path);
    if (find_file (path_buf, dem))
    {
        strcpy (dem, path_buf);
        printf ("using DEM : %s\n", dem);
    }
    else
    {
        sprintf (errmsg, "Could not find DEM auxillary data: %s\n  Check "
            "ANC_PATH environment variable.", dem);
        error_handler (false, FUNC_NAME, errmsg);
        anc_missing = true;
    }

    /* TOMS ozone file */
    sprintf (ozone, "TOMS_%d%03d.hdf", year, day);
    strcpy (path_buf, anc_path);
    if (find_file (path_buf, ozone))
    {
        strcpy (ozone, path_buf);
        printf ("using TOMS : %s\n", ozone);
    }
    else
    {
        sprintf (errmsg, "Could not find TOMS auxillary data: %s\n  Check "
            "ANC_PATH environment variable.", ozone);
        error_handler (false, FUNC_NAME, errmsg);
        anc_missing = true;
    }
    
    /* NCEP file */
    sprintf (reanalysis, "REANALYSIS_%d%03d.hdf", year, day);
    strcpy (path_buf, anc_path);
    if (find_file (path_buf, reanalysis))
    {
        strcpy (reanalysis, path_buf);
        printf ("using REANALYSIS : %s\n", reanalysis);
    }
    else
    {
        sprintf (errmsg, "Could not find NCEP REANALYSIS auxillary data: %s\n"
            "  Check ANC_PATH environment variable.", reanalysis);
        error_handler (false, FUNC_NAME, errmsg);
        anc_missing = true;
    }

    /* Check to see if missing ancillary data */
    if (anc_missing)
    {
        sprintf (errmsg, "Verify the missing auxillary data products, then "
            "try reprocessing.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Open the parameter file for lndsr for writing */
    out = fopen (lndsr_name, "w");
    if (out == NULL)
    {
        sprintf (errmsg, "Opening lndsr parameter file for writing: %s",
            lndsr_name);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    fprintf (out, "PARAMETER_FILE\n");
    fprintf (out, "XML_FILE = %s\n", input_xml);
    fprintf (out, "DEM_FILE = %s\n", dem);
    if (fopen (ozone, "r") != NULL)
    {
        /* if ozone file doesn't exist then don't write it to the parameter file
           and instead use climatology estimation */
        fprintf (out, "OZON_FIL = %s\n", ozone);
    }
    fprintf (out, "PRWV_FIL = %s\n", reanalysis);
    fprintf (out, "LEDAPSVersion = %s\n", LEDAPS_VERSION);
    fprintf (out, "END\n");
    fclose (out);
    fclose (LOG_FP);

    /* Free the metadata structure */
    free_metadata (&xml_metadata);

    /* Successful completion */
    printf ("lndpm complete.\n");
    return (SUCCESS);
}


/******************************************************************************
MODULE:  conv_date

PURPOSE: Convert year-mm-dd to Julian day of year (DOY) or vice versa.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error with the conversion
SUCCESS         Successfully converted the date

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
             Feng Gao         Revised from program "jdoy.c" by MODIS LDOPE QA
                              team
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
  1. *mm and *dd are input/output varialbes
     Input julian day number:
        input *mm = 0;
        input *dd = julian day
        output in mm-dd-yyyy format
     Input in mm-dd-yyyy format
        output *mm = 0;
        output *dd = julian day number
******************************************************************************/
int conv_date
(
    int *mm,         /* I/O: month value */
    int *dd,         /* I/O: day of month or Julian day value */
    int yyyy         /* I: year value */
)
{
    char FUNC_NAME[] = "conv_date";    /* function name */
    char errmsg[STR_SIZE];             /* error message */
    int nm, im;        /* looping vars */
    int ndays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
                       /* number of days in each month */

    /* Check for leap year */
    if ((yyyy%400 == 0) || ((yyyy%4 == 0) && (yyyy%100 != 0))) ndays[1] = 29; 

    /* Do the conversion */
    if (*mm == 0)
    { /* Convert from Julian DOY to mm-dd-yyyy */
        if (*dd <= 0) {
            sprintf (errmsg, "Invalid input date: %d %d", *dd, yyyy);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        else
        {
            for (im = 0; ((im < 12) && (*dd > 0)); im++)
                *dd -= ndays[im];
            if ((im > 12) || ((im == 12) && (*dd > 0)))
            {
                sprintf (errmsg, "Invalid input date: %d %d", *dd, yyyy);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }
            else
            {
                *mm = im;
                *dd += ndays[*mm - 1];
            }
        }
    }
    else
    { /* Convert from mm-dd-yyyy to Julian DOY */
        if ((*mm <= 0) || (*dd <= 0))
        {
            sprintf (errmsg, "Invalid input date: %d-%d-%d", *mm, *dd, yyyy);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
        else
        {
            nm = *mm - 1;
            for (im = 0; im < nm; im++)
                *dd += ndays[im];
            *mm = 0;
        }
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  open_log

PURPOSE: Open the log file as a text file for appending and add a time stamp.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error opening log file
SUCCESS         Successfully opened log file

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
1. LOG_FP and logFile are global variables for all routines in this file to
   log to.
******************************************************************************/
int open_log
(
    char app_name[]     /* I: application name/info string */
)
{
    struct tm *currtime;    /* current time in GMT */
    time_t t;               /* current time variable */
    char str[STR_SIZE];     /* current time in HMS */
    char FUNC_NAME[] = "open_log";    /* function name */
    char errmsg[STR_SIZE];            /* error message */

    /* Open the log file for appending */
    LOG_FP = fopen (logFile, "a");
    if (LOG_FP == NULL)
    {
        sprintf (errmsg, "Cannot open log file for appending: %s", logFile);
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Log the time and application name or string */
    t = time(NULL);
    currtime = (struct tm *) gmtime (&t);
    strftime (str, 100, "%FT%H:%M:%SZ", currtime);
    LOG ("\n\n##################################\n");
    LOG ("Starting <%s> on %s\n", app_name, str);
    return (SUCCESS);
}


/******************************************************************************
MODULE:  LOG

PURPOSE: Log a message to the log file

RETURN VALUE: None

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
1. LOG_FP is a global variable for all routines in this file to log to.
******************************************************************************/
void LOG
(
    char *fmt, ...   /* format of log information */
)
{
    va_list ap;
    
    va_start (ap, fmt);
    vfprintf (LOG_FP, fmt, ap);
    va_end (ap);
}


/******************************************************************************
MODULE:  scan_dir

PURPOSE: Search the specified directory, recursively, for the specified
filename.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
non-zero        File is found, and path points to full path
zero            File is not found, and path is unchanged

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
******************************************************************************/
int scan_dir
(
    char *path,   /* I: directory/path to search; upon successfully finding
                        the file, path contains the location of the file */
    char *name    /* I: filename for which to search */
)
{
    char FUNC_NAME[] = "scan_dir";    /* function name */
    char errmsg[STR_SIZE];             /* error message */
    DIR *fd;                  /* pointer to directory */
    struct dirent *dirent_p;  /* pointer to directory */
    char *nbp;                /* pointer to the end of the path */
    int found = 0;            /* was the file found? */
    
    /* Add directory separator to end of directory name if not already there */
    nbp = path + strlen (path);
    if (*nbp != '/')
        *nbp++ = '/';
      
    if (nbp+MAXNAMLEN+2 >= path + DIR_BUF_SIZE)
    {
        sprintf (errmsg, "Path name too long -- cannot search: %s", path);
        error_handler (true, FUNC_NAME, errmsg);
        return (found);
    }

    /* Check the directory to see if it's accessible */
    if ((fd = opendir (path)) == NULL)
    {
        sprintf (errmsg, "Could not read directory: %s", path);
        error_handler (true, FUNC_NAME, errmsg);
        return (found);
    }

    /* Recursively search for the file while traversing through the directory
       structure */
    while ((dirent_p = readdir(fd)) != NULL)   /* search directory */
    {
        if (dirent_p->d_ino == 0)                /* slot not in use */
            continue;
        if (strcmp (dirent_p->d_name, ".") == 0   /* ignore current ... */
        || strcmp (dirent_p->d_name, "..") == 0)  /* and parent directory */
            continue;
        
        strcat (path, dirent_p->d_name);          /* check this path */
        if ((found = find_file (path, name)) != 0)
            break;                                /* found it */
        else
            *nbp = '\0';                          /* restore directory name */
    }

    /* Close the directory pointer and successful completion */
    closedir (fd);
    return (found);
}


/******************************************************************************
MODULE:  find_file

PURPOSE: Look in the current directory for the specified file.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
non-zero        File is found, and path points to full path
zero            File is not found, and path is unchanged

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/17/2014    Gail Schmidt     Updated to follow ESPA software guidelines an
                              use common ESPA funtions.

NOTES:
******************************************************************************/
int find_file
(
    char *path,   /* I: directory/path to search; upon successfully finding
                        the file, path contains the location of the file */
    char *name    /* I: filename for which to search */
)
{
    char FUNC_NAME[] = "find_file";    /* function name */
    char errmsg[STR_SIZE];             /* error message */
    struct stat stbuf;                 /* buffer for file/directory stat */
    char pbuf[DIR_BUF_SIZE] = {0};     /* path buffer */
    int found = 0;                     /* was the file found? */
    
    /* This is the path we are checking */
    strcpy (pbuf, path);

    /* Make sure the path exists */
    if (stat (pbuf, &stbuf) != 0)
    {
        sprintf (errmsg, "Can't stat directory: %s", pbuf);
        error_handler (true, FUNC_NAME, errmsg);
        return (found);
    }

    /* If this is a directory, then search it.  Otherwise it's a file so
       check it for our filename. */
    if ((stbuf.st_mode & S_IFMT) == S_IFDIR)
        found = scan_dir (pbuf, name);
    else
        found = (strcmp(pbuf + strlen(pbuf) - strlen(name), name) == 0);

    /* If the file was found remember the location of that file */
    if (found)
        strcpy (path, pbuf);

    return (found);
}

