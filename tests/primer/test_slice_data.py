#
# def test_read_input_fasta_valid_success(self):
#     # arrange
#     self.create_files()
#     fasta = '/fasta.fa'
#     expected = [SliceData('region1_1', '5', '10', '+', 'chr1',
#                           'GTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCT'
#                           'CACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCATGATTGGCTCCTTCTCTGTGG'
#                           'GCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCCATGC'
#                           'CCTCCT'
#                           )]
#
#     # act
#     actual = self.primer._read_input_fasta(fasta)
#
#     # assert
#     self.assertEqual(actual[0].name, expected[0].name)
#     self.assertEqual(actual[0].start, expected[0].start)
#     self.assertEqual(actual[0].end, expected[0].end)
#     self.assertEqual(actual[0].strand, expected[0].strand)
#     self.assertEqual(actual[0].bases, expected[0].bases)