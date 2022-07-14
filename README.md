# Targeton Designer

Standalone targeton designer tool.

[[_TOC_]]

## Installation

Dependencies:

BedTools
```sh
sudo apt-get update
sudo apt-get install bedtools
sudo apt-get install build-essential
```

### Python3

Check Python3 version
```sh
python3 --version
```

Update if less than python3.8
```sh
sudo apt-get install python3.8-dev

sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 1
sudo update-alternatives --config python3
```
Select 3.8 then check it's updated successfully
```sh
python3 --version
```

### Clone the repo
Pull down the Targeton Designer repo and cd into it
```sh
git clone https://gitlab.internal.sanger.ac.uk/sci/targeton-designer.git
cd targeton-designer
```

### Python Virtual Environment

Requirements:
 - Python3.8+
 - Python-venv

Install venv if you haven't already got it. First install base then update to 3.8
```sh
sudo apt-get install python3-venv
sudo apt-get install python3.8-venv
```

Setting up Virtual Env:
```sh
python3 -m venv venv

source venv/bin/activate

pip3 install --upgrade pip
pip3 install -r requirements.txt

deactivate
```

Run the tests:
```sh
source venv/bin/activate

python3 -m unittest

deactivate
```

### Exonerate iPCRess

Nathan Weeks has placed Exonerate onto GitHub along with maintenance tweaks. The repo can be found here:
https://github.com/nathanweeks/exonerate

Install glib
```sh
sudo apt-get install libglib2.0-dev
```

Installing exonerate on your VM.
```sh
git clone https://github.com/nathanweeks/exonerate.git
cd exonerate
autoreconf -vfi
./configure [YOUR_CONFIGURE_OPTIONS]
make
make check
sudo make install
```

### Docker images

Upcoming feature in later releases

## Usage

### Command Line

Check Designer Version:
```sh
./designer.sh version
```

#### Designer Workflow

Running full Designer Workflow:
```sh
./designer.sh design [-h] [--bed INPUT_BED] [--fasta REF_FASTA]
```

Example Command
```sh
./designer.sh design --bed example.bed --fasta example_genomic_ref.fasta
```

#### Slicer Tool

Running Slicer tool:
```sh
./designer.sh slicer [-h] [-f5 FLANK_5] [-f3 FLANK_3] [-l LENGTH] [-o OFFSET] [-d DIR] [--bed INPUT_BED] [--fasta REF_FASTA]
```

Example command:
```sh
./designer.sh slicer --bed example.bed --fasta example_genomic_ref.fa -d example_dir
```

#### Primer3 Runner

Running Primer3:
```sh
./designer.sh primer [--fasta INPUT_FASTA] [--dir OUTPUT_FOLDER] 
```
The input fasta and BED files are intended to be sourced from the slicer tool. Examples of how these files are constructed can be found below.

Example command:
```sh
./designer.sh primer --fasta slices.fa --dir p3_output
```

#### Exonerate iPCRess

Running Exonerate iPCRess:
```sh
./designer ipcress [--dir DIR] [--fasta REF_FASTA] [--primers IPCRESS_INPUT] [--p3_csv PRIMER3_OUTPUT_CSV] 
```
Supply either a standard iPCRess input file or point P3 CSV to the output csv of the Primer3 runner.

Example command:
```sh
./designer ipcress --dir example_dir --fasta example_genomic_ref.fa --primers example_ipcress_input.txt
or
./designer ipcress --dir example_dir --fasta example_genomic_ref.fa --p3_csv example_p3_output.csv
```

### Docker

Upcoming feature in later releases

## File formats
### Genomic Reference file
A Fasta file of latest GRCh38 genome. This is used for gathering the slice sequences and retrieving primer information. 
Either supply a local genome reference file or download one from EnsEMBL and point to it with the relevant parameters:
http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/

We've included a bash script to download the reference genome:
```sh
./download_reference_genome.sh
```

Otherwise you can download it manually here:
```sh
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.alt.fa.gz
gunzip Homo_sapiens.GRCh38.dna.alt.fa.gz 
```

### Slicer Input BED File
A BED file containing the regions you wish to slice across. 

The chromosome column data must match your reference fasta file IDs. If youre reference had >chr1 then you must call chromosome 1 'chr1' in this column and vice-versa.

Note: BED effectively are applied tsv files so use tabs to separate the values. Headers are optional in BED file and can be a cause of issues if they aren't perfect. Strand is required for the slicer to ensure sequences are output in the correct orientation. Score isn't used but the field must be present for the file format to be read correctly.

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

### Slicer BED output
BED file output with row for each slice. This file will also be used for running VaLiAnt.

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
Contains the slice sequences, with their IDs including an increment, coordinates and strand in the header

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

### Primer3 Output CSV file (iPCRess Primer3 CSV Input file)
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

### iPCRess Standard Input File
Standard file format input for iPCRess. If you wish to use iPCRess as a standalone, use this file format.
Otherwise point iPCRess at the Primer3 output CSV.
Space separated text file containing name, left & right primer, minimum & maximum amplicon lengths.

| name | left | right | min_amplicon_len | max_amplicon_len |
| ---- | ---- | ----- | ---------------- | ---------------- |
| ENSE00003571441_HG6_6_LibAmp_0 | AGTGCCAGGACCTCTCCTAC | GGGAGATGCAGCCTGGGT | 200 | 300 |

Raw file
```
ENSE00003571441_HG6_6_LibAmp_0 AGTGCCAGGACCTCTCCTAC GGGAGATGCAGCCTGGGT 200 300
```
