import sys
import os

from unittest.mock import patch
from pathlib import Path

os.chdir(r'/home/ubuntu/lims2-webapp-filesystem/user/targeton-designer')
sys.path.insert(0, '')
sys.path.insert(0, 'src/')

from src.utils.arguments_parser import ParsedInputArguments
from src.cli import design_command

class DebugDesignCommand():
    def setUp(self):
        self.fasta_file_path = r"./tests/integration/fixtures/fasta_example.fa"
        self.bed_file_path = r"./tests/integration/fixtures/bed_example.bed"

    def debug_TDOutput(self):
        tmpdir = r"./tests/debug/"
        # Arrange
        # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
        with patch.object(sys, 'argv', ["./designer.sh", "design", "--bed", self.bed_file_path, "--fasta", self.fasta_file_path, "--dir", tmpdir]):
            parsed_input = ParsedInputArguments()
            args = parsed_input.get_args()

            # Act
            result = design_command(args)

            # Assert
            fields = result.fields()
            fields.remove('dir')
            for field in fields:
                field_value = getattr(result, field)
                path = Path(field_value)
                print(f"Checking file {field} -> {path.name}")
                

def main():
    test = DebugDesignCommand()
    test.setUp()
    test.debug_TDOutput()

if __name__ == '__main__':
    main()