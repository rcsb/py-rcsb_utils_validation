##
# File:    ValidationReportReaderTests.py
# Author:  J. Westbrook
# Date:    3-Feb-2019
# Version: 0.001
#
# Update:
#
#
##
"""
Tests for various utilities for extracting data from wwPDB validation
report data files.

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.validation.ValidationReportReader import ValidationReportReader

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s')
logger = logging.getLogger()


class ValidationReportReaderTests(unittest.TestCase):

    def setUp(self):
        self.__mU = MarshalUtil()
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), 'rcsb', 'mock-data')
        self.__xsdPath = os.path.join(HERE, 'test-data', 'wwpdb_validation_v002.xsd')

        self.__exampleFileXray = os.path.join(self.__dirPath, 'MOCK_VALIDATION_REPORTS', 'cb', '1cbs',
                                              '1cbs_validation.xml.gz')
        self.__cifFileXray = os.path.join(HERE, 'test-output', '1cbs_validation.cif')
        self.__exampleFileNmr = os.path.join(self.__dirPath, 'MOCK_VALIDATION_REPORTS', 'dr', '6drg',
                                             '6drg_validation.xml.gz')
        self.__cifFileNmr = os.path.join(HERE, 'test-output', '6drg_validation.cif')

        self.__dictionaryMapPath = os.path.join(HERE, 'test-data', 'vrpt_dictmap.json')
        self.__dictionaryMap = self.__mU.doImport(self.__dictionaryMapPath, format="json")

    def tearDown(self):
        pass

    def testReadXrayValidationReport(self):
        dbu = ValidationReportReader(self.__dictionaryMap)
        cL = dbu.toCif(self.__exampleFileXray)
        ok = self.__mU.doExport(self.__cifFileXray, cL, format="mmcif")
        self.assertTrue(ok)

    def testReadNmrValidationReport(self):
        dbu = ValidationReportReader(self.__dictionaryMap)
        cL = dbu.toCif(self.__exampleFileNmr)
        ok = self.__mU.doExport(self.__cifFileNmr, cL, format="mmcif")
        self.assertTrue(ok)


def readValidationReport():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ValidationReportReaderTests("testReadXrayValidationReport"))
    suiteSelect.addTest(ValidationReportReaderTests("testReadNmrValidationReport"))
    return suiteSelect


if __name__ == '__main__':
    if True:
        mySuite = readValidationReport()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
