import unittest
from io import StringIO
from os import path
import os
from pathlib import Path
import csv

import pandas as pd
import tempfile
from pyfakefs.fake_filesystem_unittest import TestCase
from unittest.mock import patch, Mock
from freezegun import freeze_time
from Bio import SeqIO

from src.utils.write_output_files import write_scoring_output, write_targeton_csv
from src.utils import write_output_files
from primer.slice_data import SliceData
from utils.write_output_files import export_retrieved_fasta


class TestWriteOutputFiles(TestCase):
    def setUp(self):
        self.setUpPyfakefs()
        self.fs.create_dir('test_dir')
        contents = (
            'region_1_1 ATCG GCTA 200 300\n'
            'region_1_2 ATCG GCTA 200 300\n'
            'region_2_1 ATCG GCTA 200 300\n'
        )
        self.fs.create_file('test_ipcress_input.txt', contents=contents)
        contents = (
            'chr1\t100\t250\tregion_1\t.\t+\n'
            'chr2\t100\t250\tregion_2\t.\t+\n'
        )
        self.fs.create_file('test.bed', contents=contents)
        self.slices = [
            SliceData('region_1', 100, 200, 'strand', 'chromosome', 'bases', 0, 0),
            SliceData('region_2', 300, 400, 'strand', 'chromosome', 'bases', 0, 0),
        ]

        self.slice_data = SliceData(
            name="TEST",
            start=100,
            end=200,
            strand="+",
            chromosome="7",
            bases="ATCGATCG",
            flanking_region=150,
            exclusion_region=0
        )

        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    @patch('builtins.print')
    def test_write_targeton_csv_success_with_timestamped_dir(self, mock_print):
        # arrange
        self.fs.create_dir('test_dir/td_310123')
        expected_dir = 'test_dir/td_310123'
        expected_file = 'test_dir/td_310123/targetons.csv'

        # act
        result = write_targeton_csv(
            'test_ipcress_input.txt', self.slices, 'test_dir/td_310123', dir_timestamped=True
        )

        # assert
        self.assertTrue(path.exists(expected_file))
        mock_print.assert_called_with('Targeton csv generated: test_dir/td_310123/targetons.csv')
        self.assertEqual(result.dir, expected_dir)
        self.assertEqual(result.csv, expected_file)

    @patch('builtins.print')
    @freeze_time('2023-01-31')
    def test_write_targeton_csv_success_with_non_timestamped_dir(self, mock_print):
        # arrange
        expected_dir = 'test_dir/td_20230131000000000000'
        expected_file = 'test_dir/td_20230131000000000000/targetons.csv'

        # act
        result = write_targeton_csv('test_ipcress_input.txt', self.slices, 'test_dir')

        # assert
        self.assertTrue(path.exists(expected_file))
        mock_print.assert_called_with(
            'Targeton csv generated: test_dir/td_20230131000000000000/targetons.csv'
        )
        self.assertEqual(result.dir, expected_dir)
        self.assertEqual(result.csv, expected_file)

    @patch('builtins.print')
    def test_write_targeton_csv_makes_csv_with_correct_contents(self, mock_print):
        # arrange
        expected = (
            'region_1_1,region_1\n'
            'region_1_2,region_1\n'
            'region_2_1,region_2\n'
        )

        # act
        write_targeton_csv('test_ipcress_input.txt', self.slices, 'test_dir', True)
        with open('test_dir/targetons.csv') as f:
            actual = f.read()

        # assert
        self.assertEqual(actual, expected)

    @patch('builtins.print')
    def test_write_scoring_output_success(self, mock_print):
        # arrange
        mock_scoring = Mock()
        expected_dir = ''
        expected_file = 'test_scoring.tsv'

        # act
        result = write_scoring_output(mock_scoring, 'test_scoring.tsv')

        # assert
        mock_scoring.save_mismatches.assert_called_with(expected_file)
        mock_print.assert_called_with('Scoring file saved: test_scoring.tsv')
        self.assertEqual(result.dir, expected_dir)
        self.assertEqual(result.tsv, expected_file)

    def test_export_to_csv(self):
        # arrange
        expected_file = 'test_scoring.tsv'
        fake_dir = 'fake_dir/test'
        self.fs.create_dir(fake_dir)
        mock_data = {'test': [1, 2, 3, 4, 5], 'test2': 'things'}
        mock_headers = ['test', 'test2']
        expected_delimiter = '\t'
        expected_file_path = Path(fake_dir) / expected_file

        # act
        result = write_output_files.export_to_csv(
            mock_data,
            fake_dir,
            expected_file,
            mock_headers,
            delimiter=expected_delimiter
        )

        sniffer = csv.Sniffer()
        with open(expected_file_path) as f:
            test_delimiter = sniffer.sniff(f.read()).delimiter
            f.seek(0)
            test_data = f.read()
        expected_read_data = f"test{expected_delimiter}test2\n[1, 2, 3, 4, 5]{expected_delimiter}things\n"
        # assert
        self.assertEqual(result, expected_file_path)
        self.assertEqual(test_delimiter, expected_delimiter)
        self.assertTrue(expected_file_path.exists())
        self.assertEqual(test_data, expected_read_data)

    def test_fasta_is_written_with_correct_header_and_sequence(self):

        fasta_path = export_retrieved_fasta(
            self.slice_data,
            self.temp_dir.name
        )

        self.assertTrue(os.path.exists(fasta_path))

        self.assertEqual(
            os.path.basename(fasta_path),
            "TEST_retrieved.fa"
        )

        expected_path = path.join(self.temp_dir.name, "TEST_retrieved.fa")
        expected_header = 'TEST:extended:GRCh38:7:100-200(+):150'

        result = export_retrieved_fasta(self.slice_data, self.temp_dir.name)

        self.assertEqual(result, expected_path)
        self.assertTrue(path.exists(expected_path))

        # verify fasta content
        records = list(SeqIO.parse(expected_path, 'fasta'))
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record.id, expected_header)
        self.assertEqual(str(record.seq), 'ATCGATCG')


    def test_empty_sequence_raises_value_error(self):

        self.slice_data.bases = ""

        with self.assertRaises(ValueError):
            export_retrieved_fasta(
                self.slice_data,
                self.temp_dir.name
            )

    def test_invalid_strand_raises_value_error(self):
        self.slice_data.strand = "?"

        with self.assertRaises(ValueError):
            export_retrieved_fasta(
                self.slice_data,
                self.temp_dir.name
            )



if __name__ == '__main__':
    unittest.main()
