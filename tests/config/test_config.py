from unittest import TestCase

from config.config import DesignerConfig, Primer3ParamsConfig


class TestDesignerConfigClass(TestCase):
    def setUp(self):
        self.config_path = 'tests/config/designer.config.json'
        self.default_config_path = 'tests/config/designer_default.config.json'

    def test_stringency_is_set(self):
        expected = [1, 0.5, 0.1]
        
        designer_config = DesignerConfig(self.config_path)
        
        self.assertEqual(designer_config.params['stringency_vector'], expected)

    def test_column_order_is_provided(self):
        expected = ["primer_type", "primer", "penalty", "stringency", "sequence", "primer_start", 
                    "primer_end", "tm", "gc_percent", "self_any_th", "self_end_th", "hairpin_th", 
                    "end_stability", "chromosome", "pre_targeton_start", "pre_targeton_end", 
                    "product_size", "targeton_id"]
        
        designer_config = DesignerConfig(self.config_path)
        
        self.assertEqual(designer_config.params['csv_column_order'], expected)

    def test_read_config(self):
        expected = {'stringency_vector': [1, 0.5, 0.1],
                    'csv_column_order': ["primer_type", "primer", "penalty", "stringency", "sequence", "primer_start", 
                                         "primer_end", "tm", "gc_percent", "self_any_th", "self_end_th", "hairpin_th", 
                                         "end_stability", "chromosome", "pre_targeton_start", "pre_targeton_end", 
                                         "product_size", "targeton_id"],
                    'filters': ["HAP1_variant", "duplicates"]
                    }
        
        designer_config = DesignerConfig(self.config_path)

        result = designer_config.read_config(
            default_config_file=self.default_config_path, config_file=self.config_path
        )
        self.assertEqual(result, expected)

    def test_no_config_file_found(self):
        incorrect_path = 'tests/config/111.config.json'

        with self.assertRaises(FileNotFoundError):
            config = DesignerConfig(incorrect_path)
            config.read_config(incorrect_path)

    def test_use_default_config(self):
        expected = {'stringency_vector': [1, 0.1],
                    'csv_column_order': ["primer_type", "primer", "penalty", "stringency", "sequence", "primer_start", 
                                         "primer_end", "tm", "gc_percent", "self_any_th", "self_end_th", "hairpin_th", 
                                         "end_stability", "chromosome", "pre_targeton_start", "pre_targeton_end", 
                                         "product_size", "targeton_id"],
                    'filters': ["HAP1_variant", "duplicates"]
                    }

        config = DesignerConfig()
        result = config.read_config(default_config_file=self.default_config_path)
        self.assertEqual(result, expected)

    def test_no_default_config_file_found(self):
        incorrect_default_path = 'tests/config/111.config.json'

        config = DesignerConfig(self.config_path)

        with self.assertRaises(FileNotFoundError):
            config.read_config(
                default_config_file=incorrect_default_path, config_file=self.config_path
            )


# class TestPrimer3ParamsClass(TestCase):
#     def setUp(self):
#         self.config_path = 'tests/config/primer3_params_simple.config.json'
#         self.default_config_path = 'tests/config/primer3_params_default.config.json'
#         self.p3_params_config = Primer3ParamsConfig(self.config_path)
#
#     def test_params_are_set(self):
#         expected = {'PRIMER_TASK': 'pick_cloning_primers'}
#
#         self.assertEqual(self.p3_params_config.params, expected)
#
#     def test_read_config(self):
#         expected = {'PRIMER_TASK': 'pick_cloning_primers'}
#
#         result = self.p3_params_config.read_config(
#             default_config_file=self.default_config_path, config_file=self.config_path
#         )
#
#         self.assertEqual(result, expected)
#
#     def test_no_config_file_found(self):
#         incorrect_path = 'tests/config/111.config.json'
#
#         with self.assertRaises(FileNotFoundError):
#             config = Primer3ParamsConfig(incorrect_path)
#             config.read_config(incorrect_path)
#
#     def test_use_default_config(self):
#         expected = {'PRIMER_TASK': 'generic'}
#
#         config = Primer3ParamsConfig()
#         result = config.read_config(default_config_file=self.default_config_path)
#
#         self.assertEqual(result, expected)
