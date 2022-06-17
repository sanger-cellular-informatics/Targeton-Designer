# Targeton Designer

Standalone targeton designer tool.

[[_TOC_]]

## Installation

Dependencies:

BedTools
```sh
sudo apt-get install bedtools
```

### Python Virtual Environment

Requirements:
 - Python3.8+
 - Python-venv

Install venv if you haven't already got it
```sh
sudo apt-get install python3-venv
```

Setting up Virtual Env:
```sh
python3 -m venv venv

source venv/bin/activate

pip3 install -r requirements.txt

deactivate
```

Run the tests:
```sh
venv/bin/activate

python3 -m unittest

deactivate
```

### Docker images

Upcoming feature in later releases

## Usage

### Python

Running Slicer tool:
```sh
python3 designer/slicer.py [-h] [-f5 FLANK_5] [-f3 FLANK_3] [-l LENGTH] [-o OFFSET] [--output_fasta OUTPUT_FASTA] [--output_slice_bed OUTPUT_SLICE_BED] bed fasta
```

Running Primer3:
```sh
python3 runner/primer3_runner.py [--seq INPUT_FASTA] [--dir OUTPUT_FOLDER] [--bed INPUT_BED] 
```
The input fasta and bed files are intended to be sourced from the slicer tool. Examples of how these files are constructed can be found below.

Example command:
```sh
python3 runner/primer3_runner.py --seq slicer_output.fa --dir ./p3_output/ --ref genomic_reference_file.fna --bed slicer_output.bed
```

### Docker

Upcoming feature in later releases

## File formats
### Genomic Reference file
A Fasta file of latest GRCh38 genome. This is used for gathering the slice sequences and retrieving primer information. 
Either supply a local genome reference file or download one from EnsEMBL and point to it with the relevant parameters:
http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/

### Slicer Input BED File
A BED file containing the regions you wish to slice across. 

The chromosome column data must match your reference fasta file IDs. If youre reference had >chr1 then you must call chromosome 1 'chr1' in this column and vice-versa.

Note: BED effectively are applied tsv files so use tabs to separate the values. Headers are optional in BED file and can be a cause of issues if the headers aren't perfect. While the file format leaves strand as optional, it is highly suggested to provide it for primer analysis. Score isn't used but necessary for the file format to be read correctly.

| chrom | chromStart | chromEnd | name | score | strand |
| ----- | ---------- | -------- | ---- | ----- | ------ |
| 1 | 42931046 | 42931206 | ENSE00003571441_HG6	| 0	| - |
| 1	| 42929593 | 42929780 | ENSE00000769557_HG8 | 0 | - |

Raw file
```
1	42931046	42931206	ENSE00003571441_HG6	0	-
1	42929593	42929780	ENSE00000769557_HG8	0	-
```

More information can be found here: https://en.wikipedia.org/wiki/BED_(file_format)

### Primer3 BED Input File (Slicer BED output)
This file is outputted by the slicer tool with the --output_bed parameter and is necessary for the primer3 runner. Similar to the input BED, this file details the locations of the slicers and their genomic locations

| chrom | chromStart | chromEnd | name | score | strand |
| ----- | ---------- | -------- | ---- | ----- | ------ |
| 1 | 42930996 | 42931206 | ENSE00003571441_HG6_1 | 0 | - | 
| 1 | 42931001 | 42931211 | ENSE00003571441_HG6_2 | 0 | - | 
| 1 | 42931006 | 42931216 | ENSE00003571441_HG6_3 | 0 | - |
| 1 | 42931011 | 42931221 | ENSE00003571441_HG6_4 | 0 | - |
| 1 | 42931016 | 42931226 | ENSE00003571441_HG6_5 | 0 | - |

Raw file
```
1	42930996	42931206	ENSE00003571441_HG6_1	0	-
1	42931001	42931211	ENSE00003571441_HG6_2	0	-
1	42931006	42931216	ENSE00003571441_HG6_3	0	-
1	42931011	42931221	ENSE00003571441_HG6_4	0	-
1	42931016	42931226	ENSE00003571441_HG6_5	0	-
```

### Primer3 Fasta Input File (Slicer Fasta output)
Slicer output of the slice sequences, their IDs, coordinates plus increment

```
>ENSE00003571441_HG6_1::1:42930996-42931206(-)
GTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCATGATTGGCTCC
TTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCCATGCCCTCCT
>ENSE00003571441_HG6_2::1:42931001-42931211(-)
TGCAGGTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCATGATTG
GCTCCTTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCCATGCC
>ENSE00003571441_HG6_3::1:42931006-42931216(-)
CCCCCTGCAGGTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCAT
GATTGGCTCCTTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCC
>ENSE00003571441_HG6_4::1:42931011-42931221(-)
CATCTCCCCCTGCAGGTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGG
GGCATGATTGGCTCCTTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAG
>ENSE00003571441_HG6_5::1:42931016-42931226(-)
GGCTGCATCTCCCCCTGCAGGTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTG
TTGGGGGCATGATTGGCTCCTTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCT
```

### Primer3 Output BED file
Genomic locations of the primers and their names. Names are incremented from 0 and given F and R depending on whether they're 5' or 3'

| chrom | chromStart | chromEnd | name | score | strand |
| ----- | ---------- | -------- | ---- | ----- | ------ |
| 1 | 42931021 | 42931039 | ENSE00003571441_HG6_6_LibAmpR_0 | 0 | - |
| 1 | 42931210 | 42931230 | ENSE00003571441_HG6_6_LibAmpF_0 | 0 | - |
| 1 | 42931021 | 42931039 | ENSE00003571441_HG6_6_LibAmpR_1 | 0 | - | 
| 1 | 42931211 | 42931230 | ENSE00003571441_HG6_6_LibAmpF_1 | 0 | - |

Raw file
```
1	42931021	42931039	ENSE00003571441_HG6_6_LibAmpR_0	0	-
1	42931210	42931230	ENSE00003571441_HG6_6_LibAmpF_0	0	-
1	42931021	42931039	ENSE00003571441_HG6_6_LibAmpR_1	0	-
1	42931211	42931230	ENSE00003571441_HG6_6_LibAmpF_1	0	-
```

### Primer3 Output CSV file
Contains all of the extra information from Primer3 for the individual primers

| primer | sequence | tm | gc_percent | penalty | self_any_th | self_end_th | hairpin_th | end_stability |
| ------ | -------- | -- | ---------- | ------- | ----------- | ----------- | ---------- | ------------- | 
| ENSE00003571441_HG6_6_LibAmpR_0 | ACCCAGGCTGCATCTCCC | 61.41508744063151 | 66.66666666666667 | 3.4150874406315097 | 9.564684449038168 | 0.0 | 0.0 | 4.3 |
| ENSE00003571441_HG6_6_LibAmpF_0 | AGTGCCAGGACCTCTCCTAC | 60.32483047348552 | 60.0 | 0.32483047348551963 | 0.0 | 0.0 | 46.300612411542886 | 3.18 |
| ENSE00003571441_HG6_6_LibAmpR_1 | ACCCAGGCTGCATCTCCC | 61.41508744063151 | 66.66666666666667 | 3.4150874406315097 | 9.564684449038168 | 0.0 | 0.0 | 4.3 | 
ENSE00003571441_HG6_6_LibAmpF_1 | AGTGCCAGGACCTCTCCTA | 58.90293358584404 | 57.89473684210526 | 2.097066414155961 | 0.0 | 0.0 | 46.300612411542886 | 2.94 | 

Raw File
```
primer,sequence,tm,gc_percent,penalty,self_any_th,self_end_th,hairpin_th,end_stability
ENSE00003571441_HG6_6_LibAmpR_0,ACCCAGGCTGCATCTCCC,61.41508744063151,66.66666666666667,3.4150874406315097,9.564684449038168,0.0,0.0,4.3
ENSE00003571441_HG6_6_LibAmpF_0,AGTGCCAGGACCTCTCCTAC,60.32483047348552,60.0,0.32483047348551963,0.0,0.0,46.300612411542886,3.18
ENSE00003571441_HG6_6_LibAmpR_1,ACCCAGGCTGCATCTCCC,61.41508744063151,66.66666666666667,3.4150874406315097,9.564684449038168,0.0,0.0,4.3
ENSE00003571441_HG6_6_LibAmpF_1,AGTGCCAGGACCTCTCCTA,58.90293358584404,57.89473684210526,2.097066414155961,0.0,0.0,46.300612411542886,2.94
```
