import unittest

from pyfakefs.fake_filesystem_unittest import TestCase

from src.ipcress.ipcress import (prettify_output)

class TestSlicer(TestCase):
    def test_prettify_output_defined(self):
        expected = 'test_cmd --pretty true'
        prettify_param = 'true'
        cmd = 'test_cmd'
        
        actual = prettify_output(prettify_param, cmd)

        self.assertEqual(actual, expected)       
 
    def test_prettify_output_undef(self):
        expected = 'test_cmd --pretty false'
        prettify_param = ''
        cmd = 'test_cmd'
        
        actual = prettify_output(prettify_param, cmd)

        self.assertEqual(actual, expected)       


if __name__ == '__main__':
    unittest.main()
