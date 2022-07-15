# Example CLI commands with test data

Run only **slicer** with test input data: 

`
./designer.sh slicer --bed ./examples/test_targeton_data.bed --fasta ../Homo_sapiens.GRCh38.dna.primary_assembly.fa --dir ./td_output
`

Run only **primer** with test input data: 

`
./designer.sh primer --fasta ./examples/test_slice_seqs.fa --dir ./td_output
`

Run only **ipcress** with test csv input: 

`
./designer.sh ipcress --dir ./td_output --fasta ../Homo_sapiens.GRCh38.dna.chromosome.1.fa --p3_csv ./examples/test_primer3_output.csv
`

Run only **ipcress** with test txt input: 

`
./designer.sh ipcress --dir ./td_output --fasta ../Homo_sapiens.GRCh38.dna.chromosome.1.fa --txt ./examples/test_ipcress_input.txt
`

Run Targeton Designer end-to-end (slicer-primer3-ipcress):

`
./designer.sh design --bed ./examples/test_targeton_data.bed --fasta ../Homo_sapiens.GRCh38.dna.primary_assembly.fa --dir ./td_output
`