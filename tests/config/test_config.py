from unittest import TestCase

from config.config import DesignerConfig


class TestDesignerConfigClass(TestCase):
    def setUp(self):
        self.config_path = 'tests/config_files/test_user_designer.config.json'
        self.default_config_path = 'tests/config_files/test_default_designer.config.json'

    def test_stringency_is_set(self):
        expected = [1, 0.5, 0.1]
        
        designer_config = DesignerConfig(args={'conf': self.config_path})
        
        self.assertEqual(designer_config.params['stringency_vector'], expected)

    def test_column_order_is_provided(self):
        expected = ["primer_type", "primer", "penalty", "stringency", "sequence", "primer_start", 
                    "primer_end", "tm", "gc_percent", "self_any_th", "self_end_th", "hairpin_th", 
                    "end_stability", "chromosome", "pre_targeton_start", "pre_targeton_end", 
                    "product_size", "targeton_id", "pair_uid"]
        
        designer_config = DesignerConfig(args={'conf': self.config_path})
        
        self.assertEqual(designer_config.params['csv_column_order'], expected)

    def test_read_config(self):
        expected = {'stringency_vector': [1, 0.5, 0.1],
                    'csv_column_order': ["primer_type", "primer", "penalty", "stringency", "sequence", "primer_start", 
                                         "primer_end", "tm", "gc_percent", "self_any_th", "self_end_th", "hairpin_th", 
                                         "end_stability", "chromosome", "pre_targeton_start", "pre_targeton_end", 
                                         "product_size", "targeton_id", "pair_uid"]}

        result = DesignerConfig.read_config(
            default_config_file=self.default_config_path, config_file=self.config_path
        )
        self.assertEqual(result, expected)

    def test_no_config_file_found(self):
        incorrect_path = 'tests/config/111.config.json'

        with self.assertRaises(FileNotFoundError):
            DesignerConfig.read_config(incorrect_path)

    def test_use_default_config(self):
        expected = {'stringency_vector': [1, 0.1],
                    'csv_column_order': ["primer_type", "primer", "penalty", "stringency", "sequence", "primer_start", 
                                         "primer_end", "tm", "gc_percent", "self_any_th", "self_end_th", "hairpin_th", 
                                         "end_stability", "chromosome", "pre_targeton_start", "pre_targeton_end", 
                                         "product_size", "targeton_id", "pair_uid"]}

        result = DesignerConfig.read_config(default_config_file=self.default_config_path)
        self.assertEqual(result, expected)

    def test_no_default_config_file_found(self):
        incorrect_default_path = 'tests/config/111.config.json'

        with self.assertRaises(FileNotFoundError):
            DesignerConfig.read_config(
                default_config_file=incorrect_default_path, config_file=self.config_path
            )
