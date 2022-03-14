import unittest

from runner.spawner import primer3_runner

class TestRunner(unittest.TestCase):
    @classmethod
    def test_runner(self):
        self.assertTrue(primer3_runner(), 'test')

if __name__ == '__main__':
    unittest.main()
