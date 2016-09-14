#! /usr/bin/env python
import sys
import os
import commands
import logging
from optparse import OptionParser
from espa import Metadata

ERROR = 1
SUCCESS = 0

class CommandFailed(Exception):
    """Exception to capture an external command failing."""
    pass

#############################################################################
# Check if a string is a float
#############################################################################
def is_float(string_to_check):
    """Check if the string is a float"""

    try:
        float(string_to_check)
        return True
    except ValueError:
        return False

#############################################################################
# Convert x/y to or from lat/lon
#
# The actual conversions are performed by the xy2geo and geo2xy tools, which
# are the tools supported by the command_name argument.
#
# This returns the row value (latitude or x) and column value (longitude or y).
#############################################################################
def convert_location(command_name, row, column, xml_filename, output_filename,
                     error_filename):
    """Convert x/y to latitude/longitude, or vice versa"""

    # Run tool to convert location
    cmdstr = ' '.join([command_name, xml_filename, str(column), str(row),
                       output_filename])
    commands.getstatusoutput(cmdstr)
    if os.path.exists(error_filename):
        message = 'Error running {0}. See file {1}'.format(command_name,
                                                           error_filename)
        raise CommandFailed(message)

    # Check that the output file exists.
    if not os.path.isfile(output_filename):
        message = ('Error running {0}.  Output file was not created. '
                   'See file {1}'.format(command_name, error_filename))
        raise CommandFailed(message)

    # Read the row and column from the file
    value_file = open(output_filename, "r")
    line = value_file.readline()
    value_file.close()
    line = line.strip()
    row_value = line.split()[8]
    column_value = line.split()[6]
    if not row_value:
        message = 'Row value in file {0} is missing. '.format(output_filename)
        raise IOError(message)
    if not column_value:
        message = 'Column value in file {0} is missing. '.format(
            output_filename)
        raise IOError(message)

    return float(row_value), float(column_value)

#############################################################################
# Retrieve air temperature values.
#
# Input: name of input ancillary file, xgrib and ygrib locations, temporary
# file for storing the air temperature values
#############################################################################
def get_air_temperatures(ancillary_filename, xgrib, ygrib,
                         air_temperature_filename):
    """Retrieve air temperature values from ancillary file"""

    # Read the air temperature values from the PRWV auxiliary file for
    # the center of the scene.  This uses the -t parameter which causes
    # the air temperature values to be sent to a data file
    cmdstr = ' '.join(['SDSreader3.0 -f ', ancillary_filename, ' -w \"',
                       str(ygrib), str(xgrib), ' 1 1 \" -t ',
                       air_temperature_filename])
    (status, output) = commands.getstatusoutput(cmdstr)
    exit_code = status >> 8
    if exit_code != 0:
        message = 'Error running SDSreader3.0.'
        raise CommandFailed(message)

    # Display the air temperature values
    air_temperature_file = open(air_temperature_filename, "r")
    air_temperature_lines = air_temperature_file.read()
    print air_temperature_lines.strip()
    air_temperature_file.close()

#############################################################################
# Retrieve scene center temperature.
#
# Input: Name of the surface reflectance text file to read
#
# Return: scene center temperature
#############################################################################
def get_center_temperature(scene_center_time, air_temperature_filename,
                           center_temperature_filename):
    """Retrieve the scene center temperature from a formatted text file"""

    # Initialize values
    temperature = -9999.0

    # Run tool to get the temperature at the scene center
    cmdstr = ' '.join(['comptemp', str(scene_center_time),
                       air_temperature_filename, center_temperature_filename])
    (status, output) = commands.getstatusoutput(cmdstr)
    exit_code = status >> 8
    if exit_code != 0:
        message = 'Error running comptemp.'
        raise CommandFailed(message)

    # Read the scene center temperature from the file
    center_temperature_file = open(center_temperature_filename, "r")
    temperature_line = center_temperature_file.readline()
    center_temperature_file.close()
    temperature_line = temperature_line.strip()
    if not temperature_line:
        message = 'Temperature value in file {0} is missing. '.format(
            center_temperature_filename)
        raise IOError(message)

    if is_float(temperature_line):
        temperature = float(temperature_line)
    else:
        message = 'Temperature value in file {0} should be a float. '.format(
            center_temperature_filename)
        raise FloatingPointError(message)

    return temperature

#############################################################################
# Retrieve and verify XML and ancillary filenames from the surface reflectance
# text file.
#
# Input: Name of the surface reflectance text file to read
#
# Return: XML filename, ancillary filename.
#############################################################################
def get_xml_and_ancillary_filenames(lndsr_filename):
    """Retrieve the XML and ancillary filenames from a formatted text file"""

    # Initialize values
    xml_filename = ''
    ancillary_filename = ''

    # Make sure the lndsr file exists
    if not os.path.exists(lndsr_filename):
        message = ('lndsr SR text file does not exist or is not '
                   'accessible: {0}'.format(lndsr_filename))
        raise IOError(message)

    # Make sure the lndsr file is really a file
    if not os.path.isfile(lndsr_filename):
        message = 'lndsr parameter is not a file: {0}'.format(lndsr_filename)
        raise IOError(message)

    # Find the XML filename and the reanalysis ancillary filename in the
    # lndsr file
    xml_filename = ''
    ancillary_filename = ''
    with open(lndsr_filename, 'r') as lndsr_file:
        lndsr_filename = lndsr_file.readlines()
        lndsr_file.close()

        # Process each line in the file
        for line in lndsr_filename:

            # Check for a match on the XML file.
            if line.startswith('XML_FILE'):
                columns = line.split()
                xml_filename = columns[2].rstrip()

            # Check for a match on the reanalysis file
            if line.startswith('PRWV_FIL'):
                columns = line.split()
                ancillary_filename = columns[2].rstrip()

    # Verify that the filenames were found
    if not xml_filename:
        message = 'Could not find the XML filename in: {0}'.format(
            lndsr_filename)
        raise IOError(message)
    if not ancillary_filename:
        message = 'Could not find the reanalysis filename in: {0}'.format(
            lndsr_filename)
        raise IOError(message)

    # Make sure the XML file exists
    if not os.path.exists(xml_filename):
        message = 'XML file does not exist or is not accessible: {0}'.format(
            xml_filename)
        raise IOError(message)

    # Make sure the ancillary file exists
    if not os.path.exists(ancillary_filename):
        message = 'Ancillary file does not exist or is not ' \
                  'accessible: {0}'.format(ancillary_filename)
        raise IOError(message)

    return xml_filename, ancillary_filename

#############################################################################
# Retrieve delta x and y values.
#
# Input: Scene center row and column, XML filename,
# intermediate filenames (scene_center_lat_long_filename,
# offset_lat_lon_filename, and xy_filename), and error filename (for xy2geo
# and geo2xy tools)
#
# Return: The delta x and delta y values are returned.  In addition,
# scene center latitude and longitude, and row and column, are returned.
#############################################################################
def get_deltas(center_row, center_column, xml_filename, center_lat_lon_filename,
               offset_lat_lon_filename, xy_filename, geoxy_error_filename):
    """Retrieve delta x and y, scene center latitude and longitude, row,
       and column"""

    # Obtain logger from logging using the module's name
    logger = logging.getLogger(__name__)

    # Initialize values
    delta_x = -9999.0
    delta_y = -9999.0
    scene_center_lat = -9999.0
    scene_center_lon = -9999.0
    cs_row = -9999.0
    cs_column = -9999.0

    # Run tool to convert scene center location to latitude and longitude
    try:
        scene_center_lat, scene_center_lon = convert_location(
            'xy2geo', center_row, center_column, xml_filename,
            center_lat_lon_filename, geoxy_error_filename)
    except (CommandFailed, IOError) as message:
        message = 'Failed to convert location for scene center.'
        logger.error(message)
        raise
    if scene_center_lat < -90.0 or scene_center_lat > 90.0:
        message = ('Scene center latitude {0} should be from -90.0 to '
                   '90.0'.format(scene_center_lat))
        raise ValueError(message)
    if scene_center_lon < -180.0 or scene_center_lon > 180.0:
        message = ('Scene center longitude {0} should be from -180.0 to '
                   '180'.format(scene_center_lon))
        raise ValueError(message)

    # Compute latitude and longitude of the point 100 pixels north from
    # the center
    try:
        offset_lat, offset_lon = convert_location('xy2geo', center_row - 100,
                                                  center_column, xml_filename,
                                                  offset_lat_lon_filename,
                                                  geoxy_error_filename)
    except (CommandFailed, IOError) as message:
        message = 'Failed to convert location for 100 pixel offset.'
        logger.error(message)
        raise
    if offset_lat < -90.0 or offset_lat > 90.0:
        message = ('Offset latitude {0} should be from -90.0 to 90.0. '
                   .format(offset_lat))
        raise ValueError(message)
    if offset_lon < -180.0 or offset_lon > 180.0:
        message = ('Offset longitude {0} should be from -180.0 to 180.0. '
                   .format(offset_lon))
        raise ValueError(message)

    # Now move to the longitude of the center of the image.
    try:
        cs_row, cs_column = convert_location('geo2xy', offset_lat,
                                             scene_center_lon, xml_filename,
                                             xy_filename, geoxy_error_filename)
    except (CommandFailed, IOError) as message:
        message = 'Failed to convert location for image center longitude.'
        logger.error(message)
        raise

    # Compute the deviation in pixels
    delta_x = float(cs_column) - center_column
    delta_y = center_row - float(cs_row)

    return delta_x, delta_y, scene_center_lat, scene_center_lon, cs_row, \
        cs_column

#############################################################################
# Retrieve specific metadata values from the XML format metadata file.
#
# Input: Name of the input metadata file to read
#
# Return: The metadata values are returned.
#############################################################################
def get_metadata(xml_filename):
    """Get various values from the specified metadata file"""

    # Initialize values
    lon1 = -9999.0
    lon2 = -9999.0
    lat1 = -9999.0
    lat2 = -9999.0
    number_lines = -9999
    number_samples = -9999
    scene_center_time = ''

    # Verify that the XML file exists
    if not os.path.exists(xml_filename):
        message = 'XML file does not exist or is not accessible: {0}'.format(
            xml_filename)
        raise IOError(message)

    # Extract values using ESPA metadata library
    espa_metadata = Metadata()
    espa_metadata.parse(xml_filename)
    global_metadata = espa_metadata.xml_object.global_metadata

    # Retrieve scene center time
    scene_center_time = str(global_metadata.scene_center_time)

    # Find sr_band1 and extract the number of lines and samples from it
    found_sr_band1 = 0
    for band in espa_metadata.xml_object.bands.band:
        if band.attrib['name'] == 'sr_band1':
            found_sr_band1 = 1
            number_lines = band.attrib['nlines']
            number_samples = band.attrib['nsamps']
            break

    if not found_sr_band1:
        message = 'Could not find XML data for surface reflectance band 1.'
        raise IOError(message)

    # Retrieve latitude and longitude values
    lon1 = global_metadata.bounding_coordinates.west
    lon2 = global_metadata.bounding_coordinates.east
    lat1 = global_metadata.bounding_coordinates.north
    lat2 = global_metadata.bounding_coordinates.south

    # Verify latitude and longitude values
    if lat1 < -90.0 or lat1 > 90.0:
        message = 'North latitude {0} should be from -90.0 to 90.0. '.format(
            lat1)
        raise ValueError(message)
    if lat2 < -90.0 or lat2 > 90.0:
        message = 'South latitude {0} should be from -90.0 to 90.0. '.format(
            lat2)
        raise ValueError(message)
    if lon1 < -180.0 or lon1 > 180.0:
        message = 'West longitude {0} should be from -180.0 to 180.0. '.format(
            lon1)
        raise ValueError(message)
    if lon2 < -180.0 or lon2 > 180.0:
        message = 'East longitude {0} should be from -180.0 to 180.0. '.format(
            lon2)
        raise ValueError(message)

    return scene_center_time, number_lines, number_samples, lon1, lon2, lat1, \
        lat2

#############################################################################
# Update the acquisition time from the metadata to scene center time in a
# format usable by the comptemp tool.
#
# Input: metadata acquisition time, longitude of scene center
#
# Return: scene center time
#############################################################################
def update_scene_center_time(scene_center_time, lonc):
    """Update the acquisition time to scene center time"""

    if scene_center_time == '00:00:00.000000Z':
        scene_center_time = 10.5 - lonc / 15
    else:
        scene_center_time_list = scene_center_time.split(':')
        hour = float(scene_center_time_list[0])
        minute = float(scene_center_time_list[1])
        scene_center_time = hour + minute / 60.0

    # Test the time for GMT
    scene_center_time_test = int(scene_center_time * 1000000)
    if scene_center_time_test < 0:
        scene_center_time += 24
        print 'We assume the scene center date is in GMT.'
    else:
        scene_center_time = int(100000 * scene_center_time) / 100000.0

    return scene_center_time

#############################################################################
# Created on May 3, 2016 by Ray Dittmeier, USGS/EROS
# Created Python script to replace lndsrbm.ksh.  Instead of parsing the output
# of a metadata dump tool, this looks up the metadata values directly.  Also,
# this reads file output from some other tools instead of standard output in
# case debug statements are sent to standard output.  Also, some verifications
# are added.
#
# History:
#
# Usage: lndsrbm.py --help prints the help message
############################################################################
class Lndsrbm():
    """Class to prepare for and run the lndsrbm application"""

    def __init__(self):
        pass

    ########################################################################
    # Description: lndsrbm will use the parameter passed for the lndsr_input.
    # If lndsr_input is None (i.e. not specified) then the command-line
    # parameters will be parsed for this information.  The lndsrbm tools and
    # application are then executed using values from the XML file specified
    # by the lndsr_input file.
    #
    # Inputs:
    #   lndsr_input - name of the surface reflectance lndsr text file
    #
    # Returns:
    #     ERROR - error running the lndsrbm application
    #     SUCCESS - successful processing
    #######################################################################
    def runLndsrbm(self, lndsr_input=None):
        """Prepare for and run the lndsrbm application"""

        # If no parameters were passed then get the info from the command
        # line
        if lndsr_input is None:

            # Get the command line argument for the lndsr file
            parser = OptionParser()
            parser.add_option('-f', '--lndsr',
                              type='string', dest='lndsr_input',
                              help='name of Landsat surface reflectance text '
                              'file',
                              metavar='FILE')
            (options, args) = parser.parse_args()

            # Validate the command-line option
            lndsr_input = options.lndsr_input  # name of the lndsr text file
            if lndsr_input is None:
                parser.error('missing lndsr SR text file command-line argument')
                return ERROR

        # Obtain logger from logging using the module's name
        logger = logging.getLogger(__name__)
        logger.info('Lndsrbm processing of surface reflectance text file: {0}'
                    .format(lndsr_input))

        # Retrieve the XML filename and the ancillary filename
        try:
            xml_filename, ancillary_filename \
                = get_xml_and_ancillary_filenames(lndsr_input)
        except IOError as message:
            message = ('Could not get XML and ancillary filenames. '
                       'Processing will terminate.')
            logger.error(message)
            raise
        logger.info('using ancillary data {0}'.format(ancillary_filename))

        # Retrieve the metadata values
        try:
            scene_center_time, number_lines, number_samples, lon1, lon2, lat1, \
                lat2 = get_metadata(xml_filename)
        except (IOError, ValueError) as message:
            message = 'Error running get_metadata.  Processing will terminate.'
            logger.error(message)
            raise
        logger.info('West bound: {0} East bound: {1} North bound {2} South '
                    'bound: {3}'.format(lon1, lon2, lat1, lat2))

        # Compute latitude and longitude of the center of the scene
        lonc = (lon1 + lon2) / 2
        latc = (lat1 + lat2) / 2
        logger.info('Center long: {0} Center lat: {1}'.format(lonc, latc))

        ygrib = int((90 - latc) * 73 / 180)
        xgrib = int((180 + lonc) * 144 / 360)
        logger.info('ygrib: {0} xgrib: {1}'.format(ygrib, xgrib))

        # Check and update the scene center time
        logger.info('Acquisition time: {0}'.format(scene_center_time))
        scene_center_time = update_scene_center_time(scene_center_time, lonc)
        logger.info('Scene time: {0}'.format(scene_center_time))

        # Get the air temperatures
        air_temperature_filename = 'tmp.airtemp'
        try:
            get_air_temperatures(ancillary_filename, xgrib, ygrib,
                                 air_temperature_filename)
        except CommandFailed as message:
            message = ('Error running get_air_temperatures.  Processing '
                       'will terminate.')
            logger.error(message)
            raise

        # Get the scene center temperature
        center_temperature_filename = 'scene_center_temperature.dat'
        try:
            temperature = get_center_temperature(
                scene_center_time, air_temperature_filename,
                center_temperature_filename)
        except (CommandFailed, IOError, FloatingPointError) as message:
            message = ('Error running get_center_temperature. '
                       'Processing will terminate.')
            logger.error(message)
            raise
        logger.info('tclear: {0}'.format(temperature))

        # Remove the geo_xy.ERROR file if it exists.  It is used to flag errors
        # with the xy2geo or geo2xy processing
        geoxy_error_filename = 'geo_xy.ERROR'
        try:
            os.remove(geoxy_error_filename)
            logger.info('Removing geoxy ERROR file: {0}'
                        .format(geoxy_error_filename))
        except OSError:
            pass

        # Compute scene orientation
        # Get row and column of image center
        center_row = float(number_lines) / 2
        center_column = float(number_samples) / 2
        logger.info('Center col/row: {0} {1}'.format(center_column, center_row))

        center_lat_lon_filename = 'scene_center_lat_lon.dat'
        offset_lat_lon_filename = 'offset_lat_lon.dat'
        xy_filename = 'xy.dat'
        try:
            delta_x, delta_y, scene_center_lat, scene_center_lon, cs_row, \
                cs_column = get_deltas(center_row, center_column, xml_filename,
                                       center_lat_lon_filename,
                                       offset_lat_lon_filename, xy_filename,
                                       geoxy_error_filename)
        except (CommandFailed, ValueError) as message:
            message = 'Could not get delta values. Processing will terminate.'
            logger.error(message)
            raise

        logger.info('Center lat/long: {0} {1}'
                    .format(scene_center_lat, scene_center_lon))
        logger.info('cscol/csrow: {0} {1}'.format(cs_column, cs_row))
        logger.info('delta x/y: {0} {1}'.format(delta_x, delta_y))

        # Update the cloud mask
        logger.info('Updating cloud mask')
        cmdstr = ' '.join(['lndsrbm --center_temp', str(temperature),
                           '--dx', str(delta_x), '--dy', str(delta_y), '--xml',
                           xml_filename])
        logger.info('{0}'.format(cmdstr))
        (status, output) = commands.getstatusoutput(cmdstr)
        logger.info(output)
        exit_code = status >> 8
        if exit_code != 0:
            message = 'Error running lndsrbm. Processing will terminate.'
            raise CommandFailed(message)

        # Clean up files
        cleanup_list = ['tmp.airtemp', center_temperature_filename,
                        center_lat_lon_filename, offset_lat_lon_filename,
                        xy_filename]
        for cleanup_filename in cleanup_list:
            try:
                os.remove(cleanup_filename)
                logger.info('Removing file: {0}'.format(cleanup_filename))
            except OSError:
                pass

        # Successful completion
        logger.info('Completion of lndsrbm.py.')

# ##### end of Lndsrbm class #####

if __name__ == "__main__":
    # Setup the default logger format and level. log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO)

    # Obtain logger from logging using the module's name
    logger = logging.getLogger(__name__)

    try:
        Lndsrbm().runLndsrbm()
    except Exception:
        logger.exception('lndsrbm.py failed')
        sys.exit(1)

    sys.exit(0)
