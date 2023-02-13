import unittest
import json

from os import linesep
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch
from tempfile import TemporaryDirectory
from src.primer_designer import PrimerDesigner, Primer, PrimerPair, iterate_design, extract_primer_data, map_primer_data
from src.utils.write_output_files import DesignOutputData
from collections import defaultdict



### Test classes
## PrimerDesigner
class TestPrimerDesignerClass(TestCase):
    def setUp(self):
        self.scoring_output_tsv_path = r"scoring_output.tsv"
        self.p3_output_csv_path = r"p3_output.csv"
        self.scoring_data = (
            r"Targeton	Primer pair	A/B/Total	0	1	2	3	4	5	6	7	8	9	10	WGE format	Score"+linesep,
            r"exon1	exon1_2_LibAmp_0	A	1	0	0	0	0	0	0	0	0	0	0	{'0': 1, '1': 0, '2': 0, '3': 0, '4': 0, '5': 0, '6': 0, '7': 0, '8': 0, '9': 0, '10': 0}"+linesep,
            r"exon1	exon1_2_LibAmp_0	B	1	0	0	0	0	0	0	0	0	0	0	{'0': 1, '1': 0, '2': 0, '3': 0, '4': 0, '5': 0, '6': 0, '7': 0, '8': 0, '9': 0, '10': 0}"+linesep,
            r"exon1	exon1_2_LibAmp_0	Total	1	0	0	0	0	0	0	0	0	0	0	{'0': 1, '1': 0, '2': 0, '3': 0, '4': 0, '5': 0, '6': 0, '7': 0, '8': 0, '9': 0, '10': 0}	0.0"
        )
        self.primer_data = (
            r"primer,sequence,chr,primer_start,primer_end,tm,gc_percent,penalty,self_any_th,self_end_th,hairpin_th,end_stability"+linesep,
            r"exon1_2_LibAmpF_0,CTGTTCTGACAGTAGAAAGGCA,chr1,55,77,58.004800503683725,45.45454545454545,3.9951994963162747,1.588400990154355,0.0,46.57005916211301,4.75"+linesep,
            r"exon1_2_LibAmpR_0,AAGAATTTTCCCCAATGGTTGCT,chr1,242,265,59.347613464584356,39.130434782608695,3.652386535415644,0.0,0.0,34.76817642661916,3.91"
            )
        self.example_dict_primer_left = {
            'chromosome' : 'chr1',
            'chr_start' : '55',
            'chr_end' : '77',
            'seq' : 'CTGTTCTGACAGTAGAAAGGCA',
            'melting_temp' : '58.004800503683725'
        }
        self.example_dict_primer_right = {
            'chromosome' : 'chr1',
            'chr_start' : '242',
            'chr_end' : '265',
            'seq' : 'AAGAATTTTCCCCAATGGTTGCT',
            'melting_temp' : '59.347613464584356'
        }
        self.example_primer_pair_dict = {
            'pair'  : 'exon1_2_LibAmp_0',
            'score' : '0.0',
            'left' : self.example_dict_primer_left,
            'right' : self.example_dict_primer_right,
            'product_size' : 210
        }
        self.example_primer_pair = PrimerPair(self.example_primer_pair_dict)
        self.example_filled_primer_designer = PrimerDesigner()
        self.example_filled_primer_designer.from_dict(self.example_primer_pair_dict)
        self.example_flat_dict = [{
            'pair': 'exon1_2_LibAmp_0',
            'score': '0.0', 
            'product_size': 210, 
            'side': 'left', 
            'chromosome': 'chr1',
            'chr_start': '55',
            'chr_end': '77',
            'seq': 'CTGTTCTGACAGTAGAAAGGCA',
            'melting_temp': '58.004800503683725'
            },
            {
            'pair': 'exon1_2_LibAmp_0',
            'score': '0.0', 
            'product_size': 210, 
            'side': 'right',
            'chromosome': 'chr1',
            'chr_start': '242',
            'chr_end': '265',
            'seq': 'AAGAATTTTCCCCAATGGTTGCT',
            'melting_temp': '59.347613464584356'
            }
        ]
        self.example_primers = [
            {'primer': 'exon1_2_LibAmpF_0', 'sequence': 'CTGTTCTGACAGTAGAAAGGCA', 'chr': 'chr1',
                'primer_start': '55', 'primer_end': '77', 'tm': '58.004800503683725', 
                'gc_percent': '45.45454545454545', 'penalty': '3.9951994963162747',
                'self_any_th': '1.588400990154355', 'self_end_th': '0.0', 
                'hairpin_th': '46.57005916211301', 'end_stability': '4.75'},
            {'primer': 'exon1_2_LibAmpR_0', 'sequence': 'AAGAATTTTCCCCAATGGTTGCT', 'chr': 'chr1',
                'primer_start': '242', 'primer_end': '265', 'tm': '59.347613464584356',
                'gc_percent': '39.130434782608695', 'penalty': '3.652386535415644',
                'self_any_th': '0.0', 'self_end_th': '0.0', 'hairpin_th': '34.76817642661916',
                'end_stability': '3.91'}
        ]
        self.example_primer = {
            'primer': 'exon1_2_LibAmpF_0', 'sequence': 'CTGTTCTGACAGTAGAAAGGCA', 'chr': 'chr1',
            'primer_start': '55', 'primer_end': '77', 'tm': '58.004800503683725', 
            'gc_percent': '45.45454545454545', 'penalty': '3.9951994963162747',
            'self_any_th': '1.588400990154355', 'self_end_th': '0.0', 
            'hairpin_th': '46.57005916211301', 'end_stability': '4.75'
            }
        self.example_scoring = [
            {'Targeton': 'exon1', 'Primer pair': 'exon1_2_LibAmp_0', 'A/B/Total': 'A',
                '0': '1', '1': '0', '2': '0', '3': '0', '4': '0', '5': '0', '6': '0',
                '7': '0', '8': '0', '9': '0', '10': '0', 
                'WGE format': "{'0': 1, '1': 0, '2': 0, '3': 0, '4': 0, '5': 0, '6': 0, '7': 0, '8': 0, '9': 0, '10': 0}",
                'Score': None}, 
            {'Targeton': 'exon1', 'Primer pair': 'exon1_2_LibAmp_0', 'A/B/Total': 'B',
                '0': '1', '1': '0', '2': '0', '3': '0', '4': '0', '5': '0', '6': '0',
                '7': '0', '8': '0', '9': '0', '10': '0', 
                'WGE format': "{'0': 1, '1': 0, '2': 0, '3': 0, '4': 0, '5': 0, '6': 0, '7': 0, '8': 0, '9': 0, '10': 0}",
                'Score': None},
            {'Targeton': 'exon1', 'Primer pair': 'exon1_2_LibAmp_0', 'A/B/Total': 'Total',
                '0': '1', '1': '0', '2': '0', '3': '0', '4': '0', '5': '0', '6': '0', 
                '7': '0', '8': '0', '9': '0', '10': '0',
                'WGE format': "{'0': 1, '1': 0, '2': 0, '3': 0, '4': 0, '5': 0, '6': 0, '7': 0, '8': 0, '9': 0, '10': 0}",
                'Score': '0.0'}
        ]
        self.example_scoring_total_dict = {
            'Targeton': 'exon1', 'Primer pair': 'exon1_2_LibAmp_0', 'A/B/Total': 'Total',
            '0': '1', '1': '0', '2': '0', '3': '0', '4': '0', '5': '0', '6': '0', 
            '7': '0', '8': '0', '9': '0', '10': '0',
            'WGE format': "{'0': 1, '1': 0, '2': 0, '3': 0, '4': 0, '5': 0, '6': 0, '7': 0, '8': 0, '9': 0, '10': 0}",
            'Score': '0.0'}
        self.example_iter_pairs_dict = defaultdict(dict, {
            'exon1_2_LibAmp_0': {
                'F': {
                    'primer': 'exon1_2_LibAmpF_0', 'sequence': 'CTGTTCTGACAGTAGAAAGGCA', 'chr': 'chr1', 'primer_start': '55',
                    'primer_end': '77', 'tm': '58.004800503683725', 'gc_percent': '45.45454545454545', 'penalty': '3.9951994963162747',
                    'self_any_th': '1.588400990154355', 'self_end_th': '0.0', 'hairpin_th': '46.57005916211301', 'end_stability': '4.75'
                    }, 
                'score': '0.0',
                'R': {
                    'primer': 'exon1_2_LibAmpR_0', 'sequence': 'AAGAATTTTCCCCAATGGTTGCT', 'chr': 'chr1', 'primer_start': '242',
                    'primer_end': '265', 'tm': '59.347613464584356', 'gc_percent': '39.130434782608695', 'penalty': '3.652386535415644',
                    'self_any_th': '0.0', 'self_end_th': '0.0', 'hairpin_th': '34.76817642661916', 'end_stability': '3.91'
                    }
                }
            }
        )
        self.example_pair_key = 'exon1_2_LibAmp_0'
        
    def create_files(self,dir):
        scoring_path = dir / Path(self.scoring_output_tsv_path)
        p3_path = dir / Path(self.p3_output_csv_path)
        with open(scoring_path, 'w') as f:
            f.writelines(self.scoring_data)
        with open(p3_path, 'w') as f:
            f.writelines(self.primer_data)
        return p3_path, scoring_path

    def test_get_primer_pairs(self):
        # Arrange
        example_primer_designer = self.example_filled_primer_designer.copy()
        # Act 
        test_primer_pairs = example_primer_designer.get_primer_pairs()
        test_primer_pairs_dicts = [pair._asdict() for pair in test_primer_pairs]
        # Assert
        self.assertListEqual(test_primer_pairs_dicts, [self.example_primer_pair._asdict()])

    def test_append_pair(self):
        # Arrange
        # example_primer_designer = PrimerDesigner()
        # example_primer_designer.primer_pairs = [self.example_primer_pair]
        example_primer_designer = self.example_filled_primer_designer.copy()
        example_primer_list = [self.example_primer_pair, self.example_primer_pair]
        # Act
        example_primer_designer.append_pair(self.example_primer_pair)
        elements_equal = [pair._asdict()==example_pair._asdict() for pair,example_pair in zip(example_primer_designer.get_primer_pairs(),example_primer_list)]
        # Assert
        self.assertTrue(all(elements_equal))
        
    def test_get_fields(self):
        # Arrange
        example_fields = [
            'pair',
            'score',
            'left',
            'right',
            'product_size'
        ]
        example_fields.sort()
        # Act
        test_fields = self.example_filled_primer_designer.get_fields()
        test_fields.sort()
        # Assert
        self.assertListEqual(example_fields, test_fields)
    
    def test_dump_json(self):
        # Arrange
        with TemporaryDirectory() as dir:
            fn = Path('example_json_output')
            if not fn.suffix:
                fn = fn.with_suffix(r'.json')
            json_path = dir/fn
            # Act
            with open(json_path, 'w') as f:
                self.example_filled_primer_designer.dump_json(f, sort_keys=True, indent=4)
            
            with open(json_path, 'r') as f:
                loaded_json_primer_designer = PrimerDesigner()
                loaded_json_primer_designer.from_dict(json.load(f))
                                                            
            # Assert
            self.assertTrue(Path(json_path).is_file())
            self.assertGreater(Path(json_path).stat().st_size, 0)
            self.assertListEqual(loaded_json_primer_designer.to_list_dicts(), self.example_filled_primer_designer.to_list_dicts())
    
    def test_to_list_dict(self):
        # Arrange
        example_list_dicts = [self.example_primer_pair_dict]
        # Act
        test_list_dicts = self.example_filled_primer_designer.to_list_dicts()
        # Assert
        self.assertListEqual(test_list_dicts,example_list_dicts)
    
    def test_flatten(self):
        # Arrange/Act
        flat_dict = self.example_filled_primer_designer.flatten()
        # Assert
        self.assertListEqual(self.example_flat_dict,flat_dict)
        
    def test_copy(self):
        # Arrange
        example_primer_designer = self.example_filled_primer_designer
        # Act
        copied_primer_designer = example_primer_designer.copy()
        # Assert
        self.assertListEqual(copied_primer_designer.to_list_dicts(),example_primer_designer.to_list_dicts())
        
    def test_prepare_primer_designer(self):
        # Arrange
        with TemporaryDirectory() as tmpdir:
            p3_path, scoring_path = self.create_files(tmpdir)
            example_design_output_data = DesignOutputData(tmpdir)
            example_design_output_data.p3_csv = p3_path
            example_design_output_data.scoring_tsv = scoring_path
            # Act
            test_primer_designer = PrimerDesigner()
            test_primer_designer.prepare_primer_designer(example_design_output_data)
        # Assert
        self.assertListEqual(test_primer_designer.to_list_dicts(),self.example_filled_primer_designer.to_list_dicts())
        
    def test__init__(self):
        # Arrange
        with TemporaryDirectory() as tmpdir:
            p3_path, scoring_path = self.create_files(tmpdir)
            example_design_output_data = DesignOutputData(tmpdir)
            example_design_output_data.p3_csv = p3_path
            example_design_output_data.scoring_tsv = scoring_path
            # Act
            test_primer_designer = PrimerDesigner(data = example_design_output_data)
        # Assert
        self.assertListEqual(PrimerDesigner().get_primer_pairs(),[]) # Empty call
        self.assertListEqual(test_primer_designer.to_list_dicts(),self.example_filled_primer_designer.to_list_dicts()) # Data supplied call
        
    def test_build_pair_classes(self):
        # Arrange
        pairs = self.example_iter_pairs_dict
        test_primer_designer = PrimerDesigner()
        test_primer_designer.build_pair_classes(pairs)
        # Assert
        self.assertListEqual(test_primer_designer.to_list_dicts(),self.example_filled_primer_designer.to_list_dicts()) # Data supplied call
            
    def test_export_to_csv(self):
        # Arrange
        fn = 'test.csv'
        with TemporaryDirectory() as tmpdir:
            # Act
            path = Path(self.example_filled_primer_designer.export_to_csv(fn, tmpdir))
            # Assert
            self.assertTrue(path.is_file())
            self.assertGreater(path.stat().st_size, 0)
        
    def test_export_to_json(self):
        # Arrange
        fn = 'test.json'
        with TemporaryDirectory() as tmpdir:
            # Act
            path = Path(self.example_filled_primer_designer.export_to_json(fn, tmpdir))
            # Assert
            self.assertTrue(path.is_file())
            self.assertGreater(path.stat().st_size, 0)
    
    def test_write_output(self):
        # Arrange
        fn = 'test.json'
        with TemporaryDirectory() as tmpdir:
            # Act
            primer_designer_output = self.example_filled_primer_designer.write_output(existing_dir = tmpdir)
            path_csv = Path(primer_designer_output.json)
            path_json = Path(primer_designer_output.csv)
            # Assert
            self.assertTrue(path_csv.is_file())
            self.assertGreater(path_csv.stat().st_size, 0)
            self.assertTrue(path_json.is_file())
            self.assertGreater(path_json.stat().st_size, 0)
            
    def test_validate_input(self):
        # Assert
        # Missing all fields
        self.assertFalse(PrimerDesigner().validate_input(DesignOutputData('')))
        with TemporaryDirectory() as tmpdir:
            # Missing data fields
            example_design_output_data = DesignOutputData(tmpdir)
            self.assertFalse(PrimerDesigner().validate_input(example_design_output_data))
            # Missing scoring
            example_design_output_data.p3_csv = self.p3_output_csv_path
            self.assertFalse(PrimerDesigner().validate_input(example_design_output_data))
            # All present should pass
            example_design_output_data.scoring_tsv = self.scoring_output_tsv_path
            self.assertTrue(PrimerDesigner().validate_input(example_design_output_data))
            
### Functions
    def test_iterate_design(self):
        # Arrange
        example_pairs = self.example_iter_pairs_dict
        # Act
        test_pairs = iterate_design(self.example_primers, [self.example_scoring_total_dict])
        # Assert
        self.assertDictEqual(test_pairs, example_pairs)

    def test_extract_primer_data(self):
        # Arrange
        example_primer_data = self.example_dict_primer_left
        # Act
        # L == F
        test_primer_data = extract_primer_data(self.example_iter_pairs_dict[self.example_pair_key]['F'])
        # Assert
        self.assertDictEqual(test_primer_data, example_primer_data)
        
    def test_map_primer_data(self): 
        # Arrange
        example_pairs = self.example_iter_pairs_dict
        pairs = defaultdict(dict)
        # Act
        for primer in self.example_primers:
            test_pairs = map_primer_data(primer, [self.example_scoring_total_dict], pairs)
        # Assert
        self.assertDictEqual(test_pairs, example_pairs)

### Primer Pair class
class TestPrimerPairClass(TestCase):
    def setUp(self):
        self.example_dict_primer_left = {
            'chromosome' : 'chr1',
            'chr_start' : '55',
            'chr_end' : '77',
            'seq' : 'CTGTTCTGACAGTAGAAAGGCA',
            'melting_temp' : '58.004800503683725'
        }
        self.example_dict_primer_right = {
            'chromosome' : 'chr1',
            'chr_start' : '242',
            'chr_end' : '265',
            'seq' : 'AAGAATTTTCCCCAATGGTTGCT',
            'melting_temp' : '59.347613464584356'
        }
        self.example_primer_pair_dict = {
            'pair'  : 'exon1_2_LibAmp_0',
            'score' : '0.0',
            'left' : self.example_dict_primer_left,
            'right' : self.example_dict_primer_right,
            'product_size' : 210
        }
        self.example_primer_pair = PrimerPair(self.example_primer_pair_dict)
        self.example_fields = ['pair', 'score', 'left', 'right', 'product_size']
        self.example_product_size = 210
        
    def test__init__(self):
        # Arrange 
        example_primer_pair = self.example_primer_pair
        # Act 
        test_primer_pair = PrimerPair(self.example_primer_pair_dict)
        # Assert 
        self.assertDictEqual(test_primer_pair._asdict(), example_primer_pair._asdict())

    def test_get_paired_dict(self):
        # Arrange 
        example_paired_dict = self.example_primer_pair_dict
        # Act 
        test_paired_dict = self.example_primer_pair.get_paired_dict()
        # Assert 
        self.assertDictEqual(test_paired_dict, example_paired_dict)
    
    def test_get_fields(self):
        # Arrange 
        example_fields = self.example_fields
        # Act 
        test_fields = self.example_primer_pair.get_fields()
        # Assert 
        self.assertListEqual(test_fields, example_fields)
    
    def test__asdict(self):
        # Arrange 
        example_dict = self.example_primer_pair_dict
        # Act 
        test_dict = self.example_primer_pair._asdict()
        # Assert 
        self.assertDictEqual(test_dict, example_dict)
    
    def test_get_product_size(self):
        # Arrange 
        example_product_size = self.example_product_size
        # Act
        test_product_size = self.example_primer_pair.get_product_size()
        # Assert 
        self.assertEqual(test_product_size, example_product_size)
    
    def test_copy(self):
        # Arrange 
        example_copy = self.example_primer_pair
        # Act
        test_copy = self.example_primer_pair.copy()
        # Assert 
        self.assertDictEqual(test_copy._asdict(), example_copy._asdict())

### Primer class
class TestPrimerPairClass(TestCase):
    def setUp(self):
        self.example_dict_primer_left = {
            'chromosome' : 'chr1',
            'chr_start' : '55',
            'chr_end' : '77',
            'seq' : 'CTGTTCTGACAGTAGAAAGGCA',
            'melting_temp' : '58.004800503683725'
        }
        self.example_dict_primer_right = {
            'chromosome' : 'chr1',
            'chr_start' : '242',
            'chr_end' : '265',
            'seq' : 'AAGAATTTTCCCCAATGGTTGCT',
            'melting_temp' : '59.347613464584356'
        }
        self.example_primer = Primer(self.example_dict_primer_left)
        self.example_fields = ['chromosome', 'chr_start', 'chr_end', 'seq', 'melting_temp']
        self.example_item_chromosome = 'chr1'
        
    def test__init__(self):
        # Arrange 
        example_primer = self.example_primer
        # Act
        test_primer = Primer(self.example_dict_primer_left)
        # Assert
        self.assertDictEqual(test_primer._asdict(), example_primer._asdict())
    
    def test__getitem__(self):
        # Arrange 
        example_item = self.example_item_chromosome
        # Act
        test_item = self.example_primer.__getitem__('chromosome')
        # Assert
        self.assertEqual(test_item, example_item)
    
    def test_get_fields(self):
        # Arrange 
        example_fields = self.example_fields
        # Act
        test_fields = self.example_primer.get_fields()
        # Assert
        self.assertListEqual(test_fields, example_fields)
    
    def test__asdict(self):
        # Arrange 
        example_dict = self.example_dict_primer_left
        # Act
        test_dict = self.example_primer._asdict()
        # Assert
        self.assertDictEqual(test_dict, example_dict)
      
            
            
if __name__ == '__main__':
    unittest.main()
