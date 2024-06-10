# Primer Designer

Standalone primer designer tool.

## Table of contents

- [Installation](#installation)
  - [Clone Repository](#clone-repository)
  - [Python Virtual Environment](#python-virtual-environment)
  - [Git Hooks (for devs)](#git-hooks-for-devs)
  - [Docker Images](#docker-images)
  - [Python debugger](#python-debugger)
- [Usage](#usage)
    -  [Primer Designer Tool](#primer-designer-tool) 
        -  [Designer Workflow (Primer3)](#designer-workflow-primer3)
        -  [Primer3 Runner](#primer3-runner)
    -  [Primer Scoring Tool](#primer-scoring-tool)
    -  [Slicer Tool](#slicer-tool)
    -  [Other Tools](#other-tools)
        - [Targeton CSV generation](#targeton-csv-generation) 
        - [Primer data collation and output to csv & Json (for benchling)](#primer-data-collation-and-output-to-csv--json-for-benchling) 
        - [Post primers to Benchling](#post-primers-to-benchling) 
- [File formats](#file-formats)
   - [Genomic Reference file](#genomic-reference-file) 
   - [Slicer Input BED File](#slicer-input-bed-file) 
   - [Slicer BED output](#slicer-bed-output) 
   - [Primer3 and Designer Fasta Input File (Slicer Fasta output)](#primer3-and-designer-fasta-input-file-slicer-fasta-output) 
   - [Primer3 Output BED file](#primer3-output-bed-file) 
   - [Primer3 Output CSV file](#primer3-output-csv-file) 
   - [Primer Designer Output example](#primer-designer-output-example) 
   - [Scoring Tool Input iPCRess file example](#scoring-tool-input-ipcress-file-example) 


## Installation

### Clone Repository
Pull down the Primer Designer repo and cd into it.
Recursively pull any submodules.
```sh
git clone --recurse-submodule https://gitlab.internal.sanger.ac.uk/sci/targeton-designer.git
cd primer-designer
```

### Python Virtual Environment

Requirements:
 - Python3.8+
 - Python-venv

Run
```sh
make
make install
make setup-venv
```
```make``` sets up the git hooks that run unittests and pycodestyle on /src and /tests on ```git push```.
```make install``` installs dependencies below.
```make setup-venv``` creates a venv at ./venv and installs requirements.txt(s)

Dependencies:

Build-essential, BedTools and Python (3.8), Python-venv (3.8)
Change ```python``` command to point to Python (3.8), ubuntu expects python3 to be a specific version for compatibility.
```sh
sudo apt-get update \
&& sudo apt-get -y install build-essential bedtools python3.8-dev python3.8-venv \
&& sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.8 2  \
&& sudo update-alternatives --config python
```

Check Python3 (base) and Python (updated) version
```sh
python3 --version
python --version
```

Setting up Virtual Env:
```sh
python -m venv venv

source venv/bin/activate

pip install -U pip wheel setuptools 
pip install -r requirements.txt
pip install -r sge-primer-scoring/requirements.txt

deactivate
```

Run the tests:
```sh
source venv/bin/activate

python -m unittest discover --start-directory ./tests --top-level-directory .
cd sge-primer-scoring
python -m unittest discover --start-directory ./tests --top-level-directory .

deactivate
```

### Git Hooks (for devs)

Located in  .githooks/
Follows standard Git hook methodology: https://git-scm.com/docs/githooks
Currently just "pre-push" that is run on ```git push```
To run either invoke ```make``` or:
```sh
git config core.hooksPath .githooks
chmod +x .githooks/*
```

### Docker images

Upcoming feature in later releases

### Python debugger
To debug with a local debugger, insert at the top of the file:
```
import sys, os
os.chdir(r'/home/ubuntu/lims2-webapp-filesystem/user/primer-designer')
sys.path.insert(0, '')
sys.path.insert(0, 'src/')
```
This allows src and submodules inside src to be found.

To debug with vscode, make sure the cwd in the debugger settings are pointed at primer-designer.
Additionally, make sure the interpreter is pointed at the correct virtual environment (venv/bin/python).

## Usage

Make designer.sh executable
```sh
chmod +x ./designer.sh
```

Check Designer Version:
```sh
./designer.sh version
```

### Primer Designer Tool

#### Designer Workflow (Primer3)

Running full Designer Workflow:
```sh
./designer.sh design [-h] [--fasta SLICE_FASTA] [--primer3_params PRIMER_CONFIG_JSON]
```

Example Command
```sh
./designer.sh design --fasta tests/integration/fixtures/fasta_example.fa
```

#### Primer3 Runner

Primer runner uses 2 types of config:
##### Primer3 parameters config: 
Defines configuration for Primer3. 
To set up user defined configuration for primer3, you can pass ``--primer3_params`` argument with file path to the ``primer`` command.
You can use ``/config/primer3_example.config.json`` file as a template. 

If there is no config file for primer3, it will run with the default configuration.

##### Designer config
Is used to specify parameters specific for the Primer Designer tool. 
 Designer tool config contains a stringency vector to specify a set of stringencies to run Primer3, and a vector to specify the column order of the csv output file.  Available columns are: 'primer_type', 'primer', 'penalty', 'stringency', 'sequence', 'primer_start', 'primer_end', 'tm', 'gc_percent', 'self_any_th', 'self_end_th', 'hairpin_th', 'end_stability', 'chromosome', 'pre_targeton_start', 'pre_targeton_end' and 'product_size'.
To set up user defined configuration for primer3, you can pass ``--conf`` argument with file path to the ``primer`` command.


To set up user defined configuration for primer3, you can create file named ``primer3.config.json`` inside of the ``/config`` folder. 
You can use ``/config/primer3_example.config.json`` file as a template. 

If there is no config file for primer3, it will run with the default configuration.

Running Primer3:
```sh
./designer.sh primer [--fasta SLICE_FASTA] [--dir OUTPUT_FOLDER] [--primer3_params PRIMER_CONFIG_JSON] [--conf DESIGNER_CONFIG_JSON]
```
The input fasta and BED files are intended to be sourced from the slicer tool. Examples of how these files are constructed can be found below.

Example command:
```sh
./designer.sh primer --fasta slice.fa --dir p3_output
```

To use kmer lists for primer generation: 
1. Set config parameters (example config/primer3_kmer.config.json)
2. Provide 2 files with kmers: homo_sapiens_11.list and homo_sapiens_16.list

##### These kmer lists can be downloaded using:

```sh
chmod +x ./download_kmer_lists.sh
./download_kmer_lists.sh
```

##### Using Filters from Configuration file:

To use the `filters` add filter names with respect to their values in `designer.config.json` file under `filters` as follows:

```
{
    "stringency_vector": [...],
    "filters": {
        "duplicates": true,
        "HAP1_variant": true
    },
    "csv_column_order": ["primer_type", "primer", "penalty", ...]
}
```

Remember to use exact names as mentioned above. By default filter `duplicates` is applied when running primer designer command.


##### Primer Pairs Ranking:

You can use the ranking functionality to rank the primer pairs. Add the column names which you need to rank by as given example below.Add the following `ranking_criteria` in your `json` configuration file:

```
{
    "stringency_vector": [...],
    "ranking_criteria":["stringency", "product_size"],
    "csv_column_order": ["primer_type", "primer", "penalty", ...]
}
```

### Primer Scoring Tool

Running primer scoring:
```sh
./designer.sh scoring [--ipcress_file IPCRESS_FILE] [--scoring_mismatch SCORING_MISMATCH] [--output_tsv OUTPUT_TSV] [--targeton_csv TARGETON_CSV] 
```

Example command:
```sh
./designer.sh scoring --ipcress_file example_ipcress_file.txt --scoring_mismatch 4 --output_tsv example_output.tsv
```
Example command with targeton csv:
```sh
./designer.sh scoring --ipcress_file example_ipcress_file.txt --scoring_mismatch 4 --output_tsv example_targeton_output.tsv --targeton_csv example_targetons.csv
```
For more information and example files see the [Primer Scoring repo](https://gitlab.internal.sanger.ac.uk/sci/sge-primer-scoring).

### Slicer Tool

Running Slicer tool:
```sh
./designer.sh slicer [-h] [-f5 FLANK_5] [-f3 FLANK_3] [-l LENGTH] [-o OFFSET] [-d DIR] [--bed INPUT_BED] [--fasta REF_FASTA]
```

Example command:
```sh
./designer.sh slicer --bed example.bed --fasta example_genomic_ref.fa -d example_dir
```

### Other Tools

#### Targeton CSV generation

To generate the targeton CSV used in primer scoring:
```sh
./designer.sh generate_targeton_csv [--primers PRIMERS] [--bed BED] [--dir DIR]
```
Example command:
```sh
./designer.sh generate_targeton_csv --primers example_ipcress_input.txt --bed example.bed --dir example_dir
```
Please note primer pair names in the iPCRess input file must be prefixed by the corresponding region name in the BED file.

#### Primer data collation and output to csv & Json (for benchling)

To collate the primer and scoring data and output to CSV & JSON file:
```sh
./designer.sh collate_primer_data [--p3_csv Primer3_output.csv] [--score_tsv scoring_output.tsv] [--dir DIR]
```

Examples of the output can be found below.

Is also run as part of the design command.

#### Post primers to Benchling

To post the top 3 primer pairs for each targeton from the Primer Designer JSON output:
```sh
./designer.sh post_primers [--primer_json PRIMER_JSON]
```
Example command:
```sh
./designer.sh post_primers --primer_json primer_designer.json
```
A message will be printed if there are less than 3 primer pairs for a particular targeton. Please note some fields on Benchling will have to be updated manually for now.

## File formats

### Genomic Reference file
A Fasta file of latest GRCh38 genome. This is used for gathering the slice sequences and retrieving primer information. 
Either supply a local genome reference file or download one from EnsEMBL and point to it with the relevant parameters:
http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/

We've included a bash script to download the reference genome:
```sh
./download_reference_genome.sh
```

Otherwise, you can download it manually here:
```sh
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
```

### Slicer Input BED File
A BED file containing the regions you wish to slice across. 

The chromosome column data must match your reference fasta file IDs. If your reference had >chr1 then you must call chromosome 1 'chr1' in this column and vice-versa.

Note: BED effectively are applied tsv files so use tabs to separate the values. Headers are optional in BED file and can be a cause of issues if they aren't perfect. Strand is required for the slicer to ensure sequences are output in the correct orientation. Score isn't used but the field must be present for the file format to be read correctly.

| chrom | chromStart | chromEnd | name                 | score | strand |
|-------|------------|----------|----------------------|-------|--------|
| 1     | 42931046   | 42931206 | ENSE00003571441_HG6	 | 0	    | -      |
| 1	    | 42929593   | 42929780 | ENSE00000769557_HG8  | 0     | -      |

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

### Primer3 and Designer Fasta Input File (Slicer Fasta output)
Contains the slice sequence, with its ID, coordinates and strand in the header.
If multiple slices are provided in this file, only the first slice will be processed. The remaining slices
will be ignored, and the user will be notified.

```
>ENSE00003571441_HG6_1::1:42930996-42931206(-)
GTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCATGATTGGCTCC
TTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCCATGCCCTCCT
```

### Primer3 Output BED file
Genomic locations of the primers and their names. Names are incremented from 0 and given F and R depending on whether they're 5' or 3'

| chrom | chromStart | chromEnd | name                            | score | strand |
|-------|------------|----------|---------------------------------|-------|--------|
| 1     | 42931021   | 42931039 | ENSE00003571441_HG6_6_LibAmpR_0 | 0     | -      |
| 1     | 42931210   | 42931230 | ENSE00003571441_HG6_6_LibAmpF_0 | 0     | -      |
| 1     | 42931021   | 42931039 | ENSE00003571441_HG6_6_LibAmpR_1 | 0     | -      | 
| 1     | 42931211   | 42931230 | ENSE00003571441_HG6_6_LibAmpF_1 | 0     | -      |

Raw file
```
1	42931021	42931039	ENSE00003571441_HG6_6_LibAmpR_0	0	-
1	42931210	42931230	ENSE00003571441_HG6_6_LibAmpF_0	0	-
1	42931021	42931039	ENSE00003571441_HG6_6_LibAmpR_1	0	-
1	42931211	42931230	ENSE00003571441_HG6_6_LibAmpF_1	0	-
```

### Primer3 Output CSV file
It contains all the additional information from Primer3 for the individual primers. Column order can be specified through the Designer tool config. 

| primer_type | primer                          | penalty | stringency | sequence             | primer_start | primer_end | tm     | gc_percent | self_any_th | self_end_th | hairpin_th | end_stability | chromosome | pre_targeton_start | pre_targeton_end | product_size | targeton_id | pair_uid                             |
|-------------|---------------------------------|---------|------------|----------------------|--------------|------------|--------|------------|-------------|-------------|------------|---------------|------------|--------------------|------------------|--------------|-------------|--------------------------------------|
| LibAmp      | ENSE00003571441_HG6_6_LibAmpR_0 | 3.415   | 1          | ACCCAGGCTGCATCTCCC   | 42931021     | 42931039   | 61.415 | 66.667     | 9.565       | 0.0         | 0.0        | 4.3           | 1          | 42931021           | 42931231         | 170          | ENSE        | 5720fdd6&#8209;095a&#8209;11ef&#8209;8b15&#8209;fa163e1ef5d5 |
| LibAmp      | ENSE00003571441_HG6_6_LibAmpF_0 | 0.324   | 1          | AGTGCCAGGACCTCTCCTAC | 42931210     | 42931230   | 60.324 | 60.0       | 0.0         | 0.0         | 46.301     | 3.18          | 1          | 42931021           | 42931231         | 170          | ENSE        | 5720fdd6&#8209;095a&#8209;11ef&#8209;8b15&#8209;fa163e1ef5d5 |
| LibAmp      | ENSE00003571441_HG6_6_LibAmpR_1 | 3.415   | 1          | ACCCAGGCTGCATCTCCC   | 42931021     | 42931039   | 61.415 | 66.667     | 9.565       | 0.0         | 0.0        | 4.3           | 1          | 42931021           | 42931231         | 171          | ENSE        | 5720fdd7&#8209;095a&#8209;11ef&#8209;8b15&#8209;fa163e1ef5d5 |
| LibAmp      | ENSE00003571441_HG6_6_LibAmpF_1 | 2.097   | 1          | AGTGCCAGGACCTCTCCTA  | 42931211     | 42931230   | 58.902 | 57.895     | 0.0         | 0.0         | 46.301     | 2.94          | 1          | 42931021           | 42931231         | 171          | ENSE        | 5720fdd7&#8209;095a&#8209;11ef&#8209;8b15&#8209;fa163e1ef5d5 |

Raw File
```
primer_type,primer,penalty,stringency,sequence,primer_start,primer_end,tm,gc_percent,self_any_th,self_end_th,hairpin_th,end_stability,chromosome,pre_targeton_start,pre_targeton_end,product_size,targeton_id,pair_id
LibAmp,ENSE00003571441_HG6_6_LibAmpR_0,3.415,1,ACCCAGGCTGCATCTCCC,42931021,42931039,61.415,66.667,9.565,0.0,0.0,4.3,1,42931021,42931231,170,ENSE,5720fdd6-095a-11ef-8b15-fa163e1ef5d5
LibAmp,ENSE00003571441_HG6_6_LibAmpF_0,0.324,1,AGTGCCAGGACCTCTCCTAC,42931210,42931230,60.324,60.0,0.0,0.0,46.301,3.18,1,42931021,42931231,170,ENSE,5720fdd6-095a-11ef-8b15-fa163e1ef5d5
LibAmp,ENSE00003571441_HG6_6_LibAmpR_1,3.415,1,ACCCAGGCTGCATCTCCC,42931021,42931039,61.415,66.667,9.565,0.0,0.0,4.3,1,42931021,42931231,171,ENSE,5720fdd7-095a-11ef-8b15-fa163e1ef5d5
LibAmp,ENSE00003571441_HG6_6_LibAmpF_1,2.097,1,AGTGCCAGGACCTCTCCTA,42931211,42931230,58.902,57.895,0.0,0.0,46.301,2.94,1,42931021,42931231,171,ENSE,5720fdd7-095a-11ef-8b15-fa163e1ef5d5
```

### Primer Designer Output example
Creates a data structure of LIST_PAIRS->PRIMER_PAIR->PRIMER

Currently output as both a JSON and CSV.

As a nested structure JSON.
```
[
    {
        "left": {
            "chr_end": "77",
            "chr_start": "55",
            "chromosome": "chr1",
            "gc_content": "45.45454545454545",
            "melting_temp": "58.004800503683725",
            "seq": "CTGTTCTGACAGTAGAAAGGCA"
        },
        "pair": "exon1_2_LibAmp_0",
        "product_size": 210,
        "right": {
            "chr_end": "265",
            "chr_start": "242",
            "chromosome": "chr1",
            "gc_content": "39.130434782608695",
            "melting_temp": "59.347613464584356",
            "seq": "AAGAATTTTCCCCAATGGTTGCT"
        },
        "score": "0.0",
        "targeton": "exon1",
        "version": "01"
    },
    {
        "left": {
            "chr_end": "77",
            "chr_start": "55",
            "chromosome": "chr1",
            "gc_content": "45.45454545454545",
            "melting_temp": "58.004800503683725",
            "seq": "CTGTTCTGACAGTAGAAAGGCA"
        },
        "pair": "exon1_2_LibAmp_1",
        "product_size": 210,
        "right": {
            "chr_end": "265",
            "chr_start": "243",
            "chromosome": "chr1",
            "gc_content": "40.90909090909091",
            "melting_temp": "57.98020807087107",
            "seq": "AAGAATTTTCCCCAATGGTTGC"
        },
        "score": "0.0",
        "targeton": "exon1",
        "version": "01"
    },
    {
        "left": {
            "chr_end": "78",
            "chr_start": "55",
            "chromosome": "chr1",
            "gc_content": "43.47826086956522",
            "melting_temp": "58.426100173219595",
            "seq": "CTGTTCTGACAGTAGAAAGGCAT"
        },
        "pair": "exon1_2_LibAmp_2",
        "product_size": 210,
        "right": {
            "chr_end": "265",
            "chr_start": "242",
            "chromosome": "chr1",
            "gc_content": "39.130434782608695",
            "melting_temp": "59.347613464584356",
            "seq": "AAGAATTTTCCCCAATGGTTGCT"
        },
        "score": "0.0",
        "targeton": "exon1",
        "version": "01"
    },
    {
        "left": {
            "chr_end": "78",
            "chr_start": "55",
            "chromosome": "chr1",
            "gc_content": "43.47826086956522",
            "melting_temp": "58.426100173219595",
            "seq": "CTGTTCTGACAGTAGAAAGGCAT"
        },
        "pair": "exon1_2_LibAmp_3",
        "product_size": 210,
        "right": {
            "chr_end": "265",
            "chr_start": "243",
            "chromosome": "chr1",
            "gc_content": "40.90909090909091",
            "melting_temp": "57.98020807087107",
            "seq": "AAGAATTTTCCCCAATGGTTGC"
        },
        "score": "0.0",
        "targeton": "exon1",
        "version": "01"
    }
]
```

Flattened to a CSV table.
```
version,pair,score,targeton,product_size,side,chromosome,chr_start,chr_end,seq,melting_temp,gc_content
01,exon1_2_LibAmp_0,0.0,exon1,210,left,chr1,55,77,CTGTTCTGACAGTAGAAAGGCA,58.004800503683725,45.45454545454545
01,exon1_2_LibAmp_0,0.0,exon1,210,right,chr1,242,265,AAGAATTTTCCCCAATGGTTGCT,59.347613464584356,39.130434782608695
01,exon1_2_LibAmp_1,0.0,exon1,210,left,chr1,55,77,CTGTTCTGACAGTAGAAAGGCA,58.004800503683725,45.45454545454545
01,exon1_2_LibAmp_1,0.0,exon1,210,right,chr1,243,265,AAGAATTTTCCCCAATGGTTGC,57.98020807087107,40.90909090909091
01,exon1_2_LibAmp_2,0.0,exon1,210,left,chr1,55,78,CTGTTCTGACAGTAGAAAGGCAT,58.426100173219595,43.47826086956522
01,exon1_2_LibAmp_2,0.0,exon1,210,right,chr1,242,265,AAGAATTTTCCCCAATGGTTGCT,59.347613464584356,39.130434782608695
01,exon1_2_LibAmp_3,0.0,exon1,210,left,chr1,55,78,CTGTTCTGACAGTAGAAAGGCAT,58.426100173219595,43.47826086956522
01,exon1_2_LibAmp_3,0.0,exon1,210,right,chr1,243,265,AAGAATTTTCCCCAATGGTTGC,57.98020807087107,40.90909090909091
```

### Scoring Tool Input iPCRess file example
Two files, stnd and err.
stnd: Space separated text file. Sequence_id contains the chromosome and description can be 'forward', 'revcomp', 'single_A' or 'single_B'.
| sequence_id | experiment_id | product_length | primer_5 | pos_5 | mismatch_5 | primer_3 | pos_3 | mismatch_3 | description |
| ----------- | ------------- | -------------- | -------- | ----- | ---------- | -------- | ----- | ---------- | ----------- |
| 19:filter(unmasked) | ID0001 | 259 | A | 44907726 | 0 | B | 44907967 | 0 | forward |

Raw file
```
ipcress: chr1:filter(unmasked) exon1_2_LibAmp_0 210 A 55 0 B 242 0 forward
ipcress: chr1:filter(unmasked) exon1_2_LibAmp_1 210 A 55 0 B 243 0 forward
ipcress: chr1:filter(unmasked) exon1_2_LibAmp_2 210 A 55 0 B 242 0 forward
ipcress: chr1:filter(unmasked) exon1_2_LibAmp_3 210 A 55 0 B 243 0 forward
-- completed ipcress analysis
```

err: is a .txt error output from ipcress or if it fails early system.run().
```
** Message: 14:55:57.646: Loaded [1] experiments
** Message: 14:55:58.141: Loaded [1] experiments
** Message: 14:55:58.649: Loaded [1] experiments
** Message: 14:55:59.167: Loaded [1] experiments
```

#### Targeton CSV example
A row for each primer pair and the region it originated from.

| Primer pair      | Region |
|------------------|--------|
| exon1_2_LibAmp_0 | exon1  |

Raw file
```
exon1_2_LibAmp_0,exon1
exon1_2_LibAmp_1,exon1
exon1_2_LibAmp_2,exon1
exon1_2_LibAmp_3,exon1
```
