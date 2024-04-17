# Example CLI commands with test data

### Primer Designer end-to-end (primer3)

Run **design** command with test input data:

`
./designer.sh design --bed ./examples/test_targeton_data.bed --fasta ../Homo_sapiens.GRCh38.dna.primary_assembly.fa --dir ./td_output
`

### Standalones 

Run only **slicer** with test input data: 

`
./designer.sh slicer --bed ./examples/test_targeton_data.bed --fasta ../Homo_sapiens.GRCh38.dna.primary_assembly.fa --dir ./td_output
`

Run only **primer** with test input data: 

`
./designer.sh primer --fasta ./examples/test_example_slice.fa --dir ./td_output
`

Generate **targeton csv** for primer scoring with test input data:

`
./designer.sh generate_targeton_csv --primers ./examples/test_ipcress_input.txt --bed ./examples/test_targeton_data.bed
`
