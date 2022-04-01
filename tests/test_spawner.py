import unittest

from runner.spawner import primer3_runner

class TestRunner(unittest.TestCase):
    @classmethod
    def test_runner(self):
        result = primer3_runner()
        self.assertTrue(result['PRIMER_LEFT_0_SEQUENCE'], 'GCATCAGTGAGTACAGCATGC')

if __name__ == '__main__':
    unittest.main()
