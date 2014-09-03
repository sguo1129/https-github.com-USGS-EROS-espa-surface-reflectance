#! /usr/bin/env python
import sys
import os
import re
import commands
import datetime
from optparse import OptionParser

ERROR = 1
SUCCESS = 0

############################################################################
# Description: logIt logs the information to the logfile (if valid) or to
# stdout if the logfile is None.
#
# Inputs:
#   msg - message to be printed/logged
#   log_handler - log file handler; if None then print to stdout
#
# Returns: nothing
#
# Notes:
############################################################################
def logIt (msg, log_handler):
    if log_handler == None:
        print msg
    else:
        log_handler.write (msg + '\n')


#############################################################################
# Created on August 26, 2014 by Gail Schmidt, USGS/EROS
# Created Python script to run the Landsat 8 surface reflectance code based
# on the inputs specified by the user.  This script will determine the input
# auxiliary file needed for processing, based on the date of the Landsat 8
# input file.
#
# History:
#   Updated on {date} by {developer}, USGS/EROS
#   {description}
#
# Usage: do_l8_sr.py --help prints the help message
############################################################################
class SurfaceReflectance():

    def __init__(self):
        pass


    ########################################################################
    # Description: runSr will use the parameters passed for the input/output
    # files, logfile, and usebin.  If input/output files are None (i.e. not
    # specified) then the command-line parameters will be parsed for this
    # information.  The surface reflectance application is then executed to
    # generate the desired outputs on the specified input file.  If a log
    # file was specified, then the output from this application will be
    # logged to that file.
    #
    # Inputs:
    #   xml_infile - name of the input XML file
    #   process_sr - specifies whether the surface reflectance processing,
    #       should be completed.  True or False.  Default is True, otherwise
    #       the processing will halt after the TOA reflectance products are
    #       complete.
    #   logfile - name of the logfile for logging information; if None then
    #       the output will be written to stdout
    #   usebin - this specifies if the surface reflectance exe resides in the
    #       $BIN directory; if None then the surface reflectance exe is
    #       expected to be in the PATH
    #
    # Returns:
    #     ERROR - error running the surface reflectance application
    #     SUCCESS - successful processing
    #
    # Notes:
    #   1. The script obtains the path of the XML file and changes
    #      directory to that path for running the surface reflectance
    #      application.  If the XML file directory is not writable, then this
    #      script exits with an error.
    #   2. If the XML file is not specified and the information is
    #      going to be grabbed from the command line, then it's assumed all
    #      the parameters will be pulled from the command line.
    #######################################################################
    def runSr (self, xml_infile=None, process_sr=None, write_toa=False, \
        logfile=None, usebin=None):
        # if no parameters were passed then get the info from the
        # command line
        if xml_infile == None:
            # get the command line argument for the XML file
            parser = OptionParser()
            parser.add_option ("-i", "--xml", type="string",
                dest="xml",
                help="name of XML file", metavar="FILE")
            parser.add_option ("-s", "--process_sr", type="string",
                dest="process_sr",
                help="process the surface reflectance products; " + \
                     "True or False (default is True)  If False, then " + \
                     "processing will halt after the TOA reflectance " + \
                     "products are complete.")
            parser.add_option ("--write_toa", dest="write_toa", default=False,
                action="store_true",
                help="write the intermediate TOA reflectance products")
            parser.add_option ("--usebin", dest="usebin", default=False,
                action="store_true",
                help="use BIN environment variable as the location of " + \
                     "surface reflectance application")
            parser.add_option ("-l", "--logfile", type="string", dest="logfile",
                help="name of optional log file", metavar="FILE")
            (options, args) = parser.parse_args()
    
            # validate the command-line options
            usebin = options.usebin          # should $BIN directory be used
            logfile = options.logfile        # name of the log file

            # XML input file
            xml_infile = options.xml
            if xml_infile == None:
                parser.error ("missing input XML file command-line argument");
                return ERROR

            # surface reflectance options
            process_sr = options.process_sr
            write_toa = options.write_toa

        # open the log file if it exists; use line buffering for the output
        log_handler = None
        if logfile != None:
            log_handler = open (logfile, 'w', buffering=1)
        msg = 'Surface reflectance processing of Landsat file: %s' % xml_infile
        logIt (msg, log_handler)
        
        # should we expect the surface reflectance application to be in the
        # PATH or in the BIN directory?
        if usebin:
            # get the BIN dir environment variable
            bin_dir = os.environ.get('BIN')
            bin_dir = bin_dir + '/'
            msg = 'BIN environment variable: %s' % bin_dir
            logIt (msg, log_handler)
        else:
            # don't use a path to the surface reflectance application
            bin_dir = ""
            msg = 'L8 surface reflectance executable expected to be in the PATH'
            logIt (msg, log_handler)
        
        # make sure the XML file exists
        if not os.path.isfile(xml_infile):
            msg = "Error: XML file does not exist or is not accessible: " + \
                xml_infile
            logIt (msg, log_handler)
            return ERROR

        # use the base XML filename and not the full path.
        base_xmlfile = os.path.basename (xml_infile)
        msg = 'Processing XML file: %s' % base_xmlfile
        logIt (msg, log_handler)
        
        # get the path of the XML file and change directory to that location
        # for running this script.  save the current working directory for
        # return to upon error or when processing is complete.  Note: use
        # abspath to handle the case when the filepath is just the filename
        # and doesn't really include a file path (i.e. the current working
        # directory).
        mydir = os.getcwd()
        xmldir = os.path.dirname (os.path.abspath (xml_infile))
        if not os.access(xmldir, os.W_OK):
            msg = 'Path of XML file is not writable: %s.  Script needs ' + \
                'write access to the XML directory.' % xmldir
            logIt (msg, log_handler)
            return ERROR
        msg = 'Changing directories for surface reflectance processing: %s' % \
            xmldir
        logIt (msg, log_handler)
        os.chdir (xmldir)

        # pull the date from the XML filename to determine which auxiliary
        # file should be used for input.  Example: LC80410272013181LGN00.xml
        # uses L8ANC2013181.hdf_fused.
        aux_file = 'L8ANC' + base_xmlfile[9:16] + '.hdf_fused'
        print 'DEBUG: auxiliary file: %s' % aux_file

        # run surface reflectance algorithm, checking the return status.  exit
        # if any errors occur.
        process_sr_opt_str = "--process_sr=true "
        write_toa_opt_str = ""

        if process_sr == "False":
            process_sr_opt_str = "--process_sr=false "
        if write_toa:
            write_toa_opt_str = "--write_toa "

        cmdstr = "%sl8_sr --xml=%s --aux=%s %s%s--verbose" % \
            ("/media/sf_Software_SandBox/L8_SR/trunk/c_version/", \
             xml_infile, aux_file, process_sr_opt_str, \
             write_toa_opt_str)
##            (bin_dir, xml_infile, aux_file, process_sr_opt_str, \
##             write_toa_opt_str)
        print 'DEBUG: l8_sr command: %s' % cmdstr
        (status, output) = commands.getstatusoutput (cmdstr)
        logIt (output, log_handler)
        exit_code = status >> 8
        if exit_code != 0:
            msg = 'Error running l8_sr.  Processing will terminate.'
            logIt (msg, log_handler)
            os.chdir (mydir)
            return ERROR
        
        # successful completion.  return to the original directory.
        os.chdir (mydir)
        msg = 'Completion of surface reflectance.'
        logIt (msg, log_handler)
        if logfile != None:
            log_handler.close()
        return SUCCESS

######end of SurfaceReflectance class######

if __name__ == "__main__":
    sys.exit (SurfaceReflectance().runSr())
