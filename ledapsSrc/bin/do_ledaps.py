#! /usr/bin/python
##################################################################
# Created on November 13, 2012 by Gail Schmidt, USGS/EROS
# Created Python script so that lnd* application status values can
# be checked for successful completion and failures can be
# flagged.
#
# Usage: do_ledaps.py --help prints the help message
##################################################################
import sys
import os
import re
import commands
from optparse import OptionParser
ERROR = 1
SUCCESS = 0

def main ():
    # get the command line argument for the metadata file
    parser = OptionParser()
    parser.add_option ("-f", "--metfile", type="string", dest="metfile",
        help="name of Landsat MTL file", metavar="FILE")
    parser.add_option ("--usebin", dest="usebin", default=False,
        action="store_true",
        help="use BIN environment variable as the location of LEDAPS apps")
    (options, args) = parser.parse_args()
    meta_file = options.metfile        # Name of the metadata file
    if meta_file == None:
        parser.error ("missing meta_file command-line argument");
        return ERROR
    print 'LEDAPS processing of Landsat metadata file: %s' % meta_file
    
    # should we expect the lnd* applications to be in the PATH or in the BIN
    # directory?
    if options.usebin:
        # get the BIN dir environment variable
        bin_dir = os.environ.get('BIN')
        bin_dir = bin_dir + '/'
        print 'BIN environment variable: %s' % bin_dir
    else:
        # don't use a path to the lnd* applications
        bin_dir = ""
        print 'LEDAPS executables expected to be in the PATH'
    
    # parse the metadata filename, strip off the _MTL.txt or _MTL.met
    meta = re.sub('\.txt$', '', meta_file)
    meta = re.sub('\.met$', '', meta)
    meta = re.sub('_MTL', '', meta)
    print 'Processing meta basefile: %s' % meta
    
    # run LEDAPS modules, checking the return status of each module. exit if
    # any errors occur.
    cmdstr = "%slndpm %s" % (bin_dir, meta_file)
#    print 'lndpm command: %s' % cmdstr
    (status, output) = commands.getstatusoutput (cmdstr)
    print output
    exit_code = status >> 8
    if exit_code != 0:
        print 'Error running lndpm.  Processing will terminate.'
        return ERROR
    
    cmdstr = "%slndcal lndcal.%s.txt" % (bin_dir, meta)
#    print 'lndcal command: %s' % cmdstr
    (status, output) = commands.getstatusoutput (cmdstr)
    print output
    exit_code = status >> 8
    if exit_code != 0:
        print 'Error running lndcal.  Processing will terminate.'
        return ERROR
    
    cmdstr = "%slndsr lndsr.%s.txt" % (bin_dir, meta)
#    print 'lndsr command: %s' % cmdstr
    (status, output) = commands.getstatusoutput (cmdstr)
    print output
    exit_code = status >> 8
    if exit_code != 0:
        print 'Error running lndsr.  Processing will terminate.'
        return ERROR
    
    cmdstr = "%slndsrbm.ksh lndsr.%s.txt" % (bin_dir, meta)
#    print 'lndsrbm command: %s' % cmdstr
    (status, output) = commands.getstatusoutput (cmdstr)
    print output
    exit_code = status >> 8
    if exit_code != 0:
        print 'Error running lndsrbm.  Processing will terminate.'
        return ERROR
    
    # successful completion
    print 'Completion of LEDAPS.'
    return SUCCESS

if __name__ == "__main__":
    sys.exit (main())
