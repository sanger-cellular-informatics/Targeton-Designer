from cli import primer_command

class TestCLIPrimerCommand(TestCase):
    def test_run_primer_command_success(self):
        test_slices_path = ''
        expected = 'LibAmpF'

        # act
        primer_command()

        # assert
        self.assertEqual(actual, expected)