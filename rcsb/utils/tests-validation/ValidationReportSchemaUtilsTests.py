##
# File:    ValidationReportSchemaUtilsTests.py
# Author:  J. Westbrook
# Date:    3-Feb-2019
# Version: 0.001
#
# Update:
#
##
"""
Tests for wwPDB validation schema to dictionary translation tools.
"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.validation.ValidationReportSchemaUtils import ValidationReportSchemaUtils

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s')
logger = logging.getLogger()


class ValidationReportSchemaUtilsTests(unittest.TestCase):

    def setUp(self):
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), 'rcsb', 'mock-data')
        self.__xsdPath = os.path.join(HERE, 'test-data', 'wwpdb_validation_v002.xsd')
        self.__dictPath = os.path.join(HERE, 'test-output', 'vrpt_mmcif_ext.dic')
        #
        # This schema mapping file is used by the XML report data file reader.
        self.__schemaMapPath = os.path.join(HERE, 'test-output', 'vrpt_schemamap.json')

    def tearDown(self):
        pass

    def testProcessXsdSchema(self):
        vrsu = ValidationReportSchemaUtils()
        sObj = vrsu.readSchema(self.__xsdPath)
        logger.debug("Returns type %r" % type(sObj))
        logger.debug("Example length %d" % len(sObj))
        ok = vrsu.buildDictionary(self.__dictPath, sObj, schemaMapPath=self.__schemaMapPath)
        self.assertTrue(ok)
        #
        schemaMap = vrsu.fetchSchemaMap(self.__schemaMapPath)
        self.assertTrue('attributes' in schemaMap)
        self.assertTrue(len(schemaMap['attributes']) > 50)


def readValidationSchema():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ValidationReportSchemaUtilsTests("testProcessXsdSchema"))
    return suiteSelect


if __name__ == '__main__':

    if True:
        mySuite = readValidationSchema()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
