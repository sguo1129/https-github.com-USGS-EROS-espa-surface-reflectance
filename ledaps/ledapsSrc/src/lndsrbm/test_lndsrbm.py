#! /usr/bin/env python

import unittest
import os

from lndsrbm import convert_location, get_air_temperatures, get_center_temperature, get_xml_and_ancillary_filenames, get_deltas, get_metadata, update_scene_center_time

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
        xml_filename = "LE70050582002003EDC00.xml"
        center_lat_lon_filename = "scene_center_lat_lon.dat"
        geoxy_error_filename = "geo_xy.ERROR"
        if os.path.exists(geoxy_error_filename):
            os.remove(geoxy_error_filename)
        scene_center_lat, scene_center_lon, status, status_string \
            = convert_location(
                exe_dir + "/xy2geo", center_row, center_column, xml_filename,
                center_lat_lon_filename, geoxy_error_filename)

        self.assertEqual(scene_center_lat, 2.892448)
        self.assertEqual(scene_center_lon, -70.151414)
        self.assertEqual(status, 0)
        self.assertEqual(status_string, "Success")

        # Test converting from lat/lon to x/y.
        exe_dir = os.environ['BIN']
        offset_lat = 2.919583
        scene_center_lon = -70.151414
        xml_filename = "LE70050582002003EDC00.xml"
        xy_filename = "xy.dat"
        geoxy_error_filename = "geo_xy.ERROR"
        if os.path.exists(geoxy_error_filename):
            os.remove(geoxy_error_filename)
        cs_row, cs_column, status, status_string \
            = convert_location(
                exe_dir + "/geo2xy", offset_lat, scene_center_lon, xml_filename,
                xy_filename, geoxy_error_filename)

        self.assertEqual(cs_row, 3390.501821)
        self.assertEqual(cs_column, 3900.600519)
        self.assertEqual(status, 0)
        self.assertEqual(status_string, "Success")

        # Test converting from x/y to lat/lon with a failure.
        exe_dir = os.environ['BIN']
        center_row = 3490.5
        center_column = 3900.5
        xml_filename = "LE70050582002003EDC00.ABC"
        center_lat_lon_filename = "scene_center_lat_lon.dat"
        geoxy_error_filename = "geo_xy.ERROR"
        if os.path.exists(geoxy_error_filename):
            os.remove(geoxy_error_filename)
        scene_center_lat, scene_center_lon, status, status_string \
            = convert_location(
                exe_dir + "/xy2geo", center_row, center_column, xml_filename,
                center_lat_lon_filename, geoxy_error_filename)
        status_string_list = status_string.split()

        self.assertEqual(scene_center_lat, -9999.0)
        self.assertEqual(scene_center_lon, -9999.0)
        self.assertEqual(status, 1)
        self.assertEqual(status_string_list[0], "Error")
        self.assertEqual(status_string_list[1], "running")


    def test_get_center_temperature(self):
        """Tests getting air temperature and scene center temperature"""

        # Test successful retrieval of air temperatures.
        ancillary_filename \
            = "/usr/local/ledaps/ANC/REANALYSIS/RE_2002/REANALYSIS_2002003.hdf"
        air_temperature_filename = "tmp.airtemp"
        xgrib = 43
        ygrib = 35
        status, status_string = get_air_temperatures(ancillary_filename,
                                                     xgrib, ygrib,
                                                     air_temperature_filename)
        self.assertEqual(status, 0)
        self.assertEqual(status_string, "Success")

        # Test successful retrieval of scene center temperature, based on
        # the air temperatures just retrieved.
        center_temperature_filename = "scene_center_temperature.dat"
        scene_center_time = 14.71666
        temperature, status, status_string = get_center_temperature(
            scene_center_time, air_temperature_filename,
            center_temperature_filename)

        self.assertEqual(temperature, 296.756409)
        self.assertEqual(status, 0)
        self.assertEqual(status_string, "Success")


    def test_get_xml_and_ancillary_filenames(self):
        """Tests getting the XML filename and the ancillary filename"""

        # Test successful retrieval of XML and ancillary filenames.
        lndsr_filename = "lndsr.LE70050582002003EDC00.txt"
        xml_filename, ancillary_filename, status, status_string \
            = get_xml_and_ancillary_filenames(lndsr_filename)

        self.assertEqual(xml_filename, "LE70050582002003EDC00.xml")
        self.assertEqual(ancillary_filename, 
            "/usr/local/ledaps/ANC/REANALYSIS/RE_2002/REANALYSIS_2002003.hdf")
        self.assertEqual(status, 0)
        self.assertEqual(status_string, "Success")


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
        xml_filename = "LE70050582002003EDC00.xml"
        scene_center_time, number_lines, number_samples, lon1, lon2, lat1, \
            lat2, status, status_string = get_metadata(xml_filename)

        self.assertEqual(scene_center_time, "14:43:14.7339290Z")
        self.assertEqual(number_lines, "6981")
        self.assertEqual(number_samples, "7801")
        self.assertEqual(lon1, -71.206259)
        self.assertEqual(lon2, -69.098776)
        self.assertEqual(lat1, 3.840657)
        self.assertEqual(lat2, 1.944488)
        self.assertEqual(status, 0)
        self.assertEqual(status_string, "Success")

        # Test unsuccessful retrieval of metadata fields where the XML file
        # does not exist.
        xml_filename = "not_there.xml"
        scene_center_time, number_lines, number_samples, lon1, lon2, lat1, \
            lat2, status, status_string = get_metadata(xml_filename)

        self.assertEqual(status, 1)
        self.assertEqual(status_string, "XML file does not exist or is not "
            "accessible: not_there.xml")

if __name__ == '__main__':
    unittest.main(verbosity=2)

