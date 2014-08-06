#! /usr/bin/env python

# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
# This is only coded to help me(Ron Dilley) run it for test SR output data.
# This code should never be the official execution script, in its current form.
# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
# TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO

import os
import sys
import glob
import subprocess
from cStringIO import StringIO
from argparse import ArgumentParser


# ============================================================================
def execute_cmd(cmd):
    '''
    Description:
      Execute a command line and return SUCCESS or ERROR

    Returns:
        output - The stdout and/or stderr from the executed command.
    '''

    output = ''
    proc = None
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, shell=True)
        output = proc.communicate()[0]

        if proc.returncode < 0:
            message = "Application terminated by signal [%s]" % cmd
            raise Exception(message)

        if proc.returncode != 0:
            message = "Application failed to execute [%s]" % cmd
            raise Exception(message)

        application_exitcode = proc.returncode >> 8
        if application_exitcode != 0:
            message = "Application [%s] returned error code [%d]" \
                      % (cmd, application_exitcode)
            raise Exception(message)

    finally:
        del proc

    return output
# END - execute_cmd


#=============================================================================
def build_argument_parser():
    '''
    Description:
      Build the command line argument parser.
    '''

    # Create a command line argument parser
    description = "Process an L8 scene to Surface Reflectance"
    parser = ArgumentParser(description=description)

    # Add parameters
    parser.add_argument('--input-dir', action='store',
                        dest='input_dir', required=True,
                        help="directory where a scene directory resides")

    parser.add_argument('--gen-sr',
                        action='store_true', dest='gen_sr', default=False,
                        help="generate the SR product")


    return parser
# END - build_argument_parser


#=============================================================================
if __name__ == '__main__':
    '''
    Description:
        Process an L8 scene to Surface Reflectance
    '''

    # Build the command line argument parser
    parser = build_argument_parser()

    # Parse the command line arguments
    args = parser.parse_args()

    scene_dirs = glob.glob(args.input_dir + '/LC*')

    for scene_dir in scene_dirs:
        sceneid = os.path.basename(scene_dir)
        gen_input_cmd = ['./processldcm.sh', args.input_dir, sceneid]
        gen_input_cmd = ' '.join(gen_input_cmd)
        output = execute_cmd(gen_input_cmd)
        if len(output) > 0:
            print "Scene ID [" + sceneid + "]"
            print output

        if args.gen_sr:
            gen_sr_cmd = ['./LDCMSR-v1.3', '<', sceneid + '.input', '>',
                          sceneid + '.log']

            output = ''
            try:
                gen_sr_cmd = ' '.join(gen_sr_cmd)
                output = execute_cmd(gen_sr_cmd)
                os.rename('correcteddata.hdf', 'sr-' + sceneid + '.hdf')
            except Exception, e:
                print str(e)
                pass

            if len(output) > 0:
                print "Scene ID [" + sceneid + "]"
                print output

    sys.exit(0)

