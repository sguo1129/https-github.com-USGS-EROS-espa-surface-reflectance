#! /usr/bin/env python

import unittest
import os

from lndsrbm import convert_location, get_air_temperatures, get_center_temperature, get_xml_and_ancillary_filenames, get_metadata, update_scene_center_time

class TestLndsrbm(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):

        # Clean up files left by the tests.
        air_temperature_filename = "tmp.airtemp"
        center_temperature_filename = "scene_center_temperature.dat"
        xy_filename = "xy.dat"
        center_lat_lon_filename = "scene_center_lat_lon.dat"
        offset_lat_lon_filename = "offset_lat_lon.dat"
        geoxy_error_filename = "geo_xy.ERROR"
        cleanup_list = [air_temperature_filename, center_temperature_filename,
                        xy_filename, center_lat_lon_filename,
                        offset_lat_lon_filename, geoxy_error_filename]
        for cleanup_filename in cleanup_list:
            try:
                os.remove(cleanup_filename)
            except OSError:
                pass


    def test_convert_location(self):
        """Tests converting from x/y to lat/lon and from lat/lon to x/y"""

        # Test converting from x/y to lat/lon.
        exe_dir = os.environ['BIN']
        center_row = 3490.5
        center_column = 3900.5
        xml_filename = "LE70450322002059EDC00.xml"
        center_lat_lon_filename = "scene_center_lat_lon.dat"
        geoxy_error_filename = "geo_xy.ERROR"
        if os.path.exists(geoxy_error_filename):
            os.remove(geoxy_error_filename)
        scene_center_lat, scene_center_lon = convert_location(
            exe_dir + "/xy2geo", center_row, center_column, xml_filename,
            center_lat_lon_filename, geoxy_error_filename)

        self.assertEqual(scene_center_lat, 40.374823)
        self.assertEqual(scene_center_lon, -122.857111)

        # Test converting from lat/lon to x/y.
        exe_dir = os.environ['BIN']
        offset_lat = 2.919583
        scene_center_lon = -70.151414
        xml_filename = "LE70450322002059EDC00.xml"
        xy_filename = "xy.dat"
        geoxy_error_filename = "geo_xy.ERROR"
        if os.path.exists(geoxy_error_filename):
            os.remove(geoxy_error_filename)
        cs_row, cs_column = convert_location(
            exe_dir + "/geo2xy", offset_lat, scene_center_lon, xml_filename,
            xy_filename, geoxy_error_filename)

        self.assertEqual(cs_row, 134892.664364)
        self.assertEqual(cs_column, 232997.908289)

        # Test converting from x/y to lat/lon with a failure.
        exe_dir = os.environ['BIN']
        center_row = 3490.5
        center_column = 3900.5
        xml_filename = "LE70050582002003EDC00.not_there"
        center_lat_lon_filename = "scene_center_lat_lon.dat"
        geoxy_error_filename = "geo_xy.ERROR"
        if os.path.exists(geoxy_error_filename):
            os.remove(geoxy_error_filename)
        status = 1
        try:
            scene_center_lat, scene_center_lon = convert_location(
                exe_dir + "/xy2geo", center_row, center_column, xml_filename,
                center_lat_lon_filename, geoxy_error_filename)
        except:
            # There should be an exception because the XML file isn't there.
            status = 0

        self.assertEqual(status, 0)

    def test_get_center_temperature(self):
        """Tests getting air temperature and scene center temperature"""

        # Test successful retrieval of air temperatures.
        ancillary_filename \
            = "/usr/local/ledaps/ANC/REANALYSIS/RE_2002/REANALYSIS_2002003.hdf"
        air_temperature_filename = "tmp.airtemp"
        xgrib = 43
        ygrib = 35
        get_air_temperatures(ancillary_filename, xgrib, ygrib,
                             air_temperature_filename)

        # Test successful retrieval of scene center temperature, based on
        # the air temperatures just retrieved.
        center_temperature_filename = "scene_center_temperature.dat"
        scene_center_time = 14.71666
        temperature = get_center_temperature(
            scene_center_time, air_temperature_filename,
            center_temperature_filename)

        self.assertEqual(temperature, 296.756409)


    def test_get_xml_and_ancillary_filenames(self):
        """Tests getting the XML filename and the ancillary filename"""

        # Test successful retrieval of XML and ancillary filenames.
        lndsr_filename = "lndsr.LE70450322002059EDC00.txt"
        xml_filename, ancillary_filename \
            = get_xml_and_ancillary_filenames(lndsr_filename)

        self.assertEqual(xml_filename, "LE70450322002059EDC00.xml")
        self.assertEqual(ancillary_filename,
            "/usr/local/ledaps/ANC/REANALYSIS/RE_2002/REANALYSIS_2002003.hdf")


    def test_update_scene_center_time(self):
        """Tests getting the scene center time in the required format"""

        # Test successful calculation and formatting of scene center time.
        scene_center_time = "14:43:14.7339290Z"
        lonc = -70.1525175
        scene_center_time = update_scene_center_time(scene_center_time, lonc)

        self.assertEqual(scene_center_time, 14.71666)


    def test_get_metadata(self):
        """Tests various elements that are retrieved with get_metadata()"""

        # Test successful retrieval of metadata fields.
        xml_filename = "LE70450322002059EDC00.xml"
        scene_center_time, number_lines, number_samples, lon1, lon2, lat1, \
            lat2 = get_metadata(xml_filename)

        self.assertEqual(scene_center_time, "18:40:25.5283427Z")
        self.assertEqual(number_lines, "7211")
        self.assertEqual(number_samples, "7931")
        self.assertEqual(lon1, -124.25339)
        self.assertEqual(lon2, -121.411018)
        self.assertEqual(lat1, 41.318504)
        self.assertEqual(lat2, 39.359287)

        # Test unsuccessful retrieval of metadata fields where the XML file
        # does not exist.
        xml_filename = "not_there.xml"
        status = 1
        try:
            scene_center_time, number_lines, number_samples, lon1, lon2, lat1, \
                lat2 = get_metadata(xml_filename)
        except IOError:
            # There should be an IOError exception since the file isn't there.
            status = 0

        self.assertEqual(status, 0)

if __name__ == '__main__':
    unittest.main(verbosity=2)

