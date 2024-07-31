# Primer Designer

A standalone primer designer tool that takes a sequence of a target region of the genome (i.e., a targeton) 
and identifies the appropriate LibAmp primers at both ends of the targeton. 
The Primer3 Designer tool includes filtering and ranking of primers.

## Table of contents

1. [Installation](#1-installation)
   1. [Clone Repository](#11-clone-repository)
   2. [Install Dependencies](#12-install-dependencies)
   3. [Setup Environment](#13-setup-environment)
      1. [Setting up Python Virtual Environment](#131-setting-up-python-virtual-environment)
      2. [Downloading kmer lists for primer generation](#132-downloading-kmer-lists-for-primer-generation)
      3. [Make the Designer Script Executable and Check Version](#133-make-the-designer-script-executable-and-check-version)
   4. [Running unit tests](#14-running-unit-tests)
2. [Usage](#2-usage)
    1. [Quick start with primer command](#21-quick-start-with-primer-command)
    2.  [Primer Designer Tool](#22-primer-designer-tool) 
        1.  [Primer3 Runner](#221-primer3-runner)
            1. [Primer3 parameters config](#2211-primer3-parameters-config)
            2. [Designer config](#2212-designer-config)
            3. [Running Primer3](#2213-running-primer3)
            4. [Applying filters from the designer config file](#2214-applying-filters-from-the-designer-config-file)
            5. [Applying ranking from the designer config file](#2215-applying-ranking-from-the-designer-config-file)
            6. [Specifying column order through the designer config file](#2216-specifying-column-order-through-the-designer-config-file)
            7. [Using the designer config file to set command-line arguments](#2217-using-the-designer-config-file-to-set-command-line-arguments)
3. [File formats](#3-file-formats)
   1. [Primer3 and Designer Fasta Input File (Slicer Fasta output)](#31-primer3-and-designer-fasta-input-file-slicer-fasta-output) 
   2. [Primer3 Output BED file](#32-primer3-output-bed-file) 
   3. [Primer3 Output CSV file](#33-primer3-output-csv-file) 
   4. [Primer3 Output Optimal Primer Pairs CSV file](#34-primer3-output-optimal-primer-pairs-csv-file) 
   5. [Primer3 Output Discarded Primer Pairs CSV file](#35-primer3-output-discarded-primer-pairs-csv-file) 
   6. [Genomic Reference file](#36-genomic-reference-file) 
   7. [Slicer Input BED File](#37-slicer-input-bed-file) 
   8. [Slicer BED output](#38-slicer-bed-output)  
   9. [Scoring Tool Input iPCRess file example](#39-scoring-tool-input-ipcress-file-example) 
   10. [Primer Designer Output example](#310-primer-designer-output-example)
4. [For Developers](#4-for-developers)
   1. [Git Hooks](#41-git-hooks)
   2. [Python debugger](#42-python-debugger)
5. [Tools and commands no longer in use](#5-tools-and-commands-no-longer-in-use)
   1. [Designer Workflow (Primer3)](#51-designer-workflow-primer3)
   2. [Primer Scoring Tool](#52-primer-scoring-tool)
   3. [Slicer Tool](#53-slicer-tool)
   4. [Targeton CSV generation](#54-targeton-csv-generation)
   5. [Primer data collation and output to CSV and JSON (for Benchling)](#55-primer-data-collation-and-output-to-csv-and-json-for-benchling) 
   6. [Post primers to Benchling](#56-post-primers-to-benchling) 
6. [Upcoming releases](#6-upcoming-releases)


## 1. Installation

### 1.1 Clone Repository
Clone the Primer Designer repository and `cd` into it by using the following command.
Recursively pull any submodules.
```sh
git clone --recurse-submodule https://gitlab.internal.sanger.ac.uk/sci/targeton-designer.git
cd targeton-designer
```

Please note that the repository on GitLab is currently named `targeton-designer` due to the project's previous name, Targeton Designer. 
The project has since been renamed to Primer Designer, and the repository name will be updated accordingly in the future. 
We apologise for any confusion this may cause.

### 1.2 Install Dependencies

Build-essential, BedTools and Python (3.8), Python-venv (3.8).
Change ```python``` command to point to Python (3.8), Ubuntu expects python3 to be a specific version for compatibility.
```sh
sudo apt-get update \
&& sudo apt-get -y install build-essential bedtools python3.8-dev python3.8-venv \
&& sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.8 2  \
&& sudo update-alternatives --config python \
&& sudo apt-get install python3-pip \
&& sudo apt install make
```

Check Python3 (base) and Python (updated) version
```sh
python3 --version
python --version
```

### 1.3 Setup Environment

Requirements:
 - Python3.8+
 - Python-venv

Run

```sh
sudo make install
```

```make install``` installs dependencies and setup environment for running primer designer commands.


##### 1.3.1 Setting up Python Virtual Environment
Install Python virtual environment using the following command:
```
sudo pip3 install virtualenv 
```

Then, by using the following command, create a virtual environment and activate it:

```
python -m venv venv
source venv/bin/activate
```

To deactivate the virtual environment type the following command and hit enter:
```
deactivate
```

After creating the virtual environment, you need to install python packages from `requirements.txt` using the following commands:

```
pip install -U pip wheel setuptools 
pip install -r requirements.txt
pip install -r sge-primer-scoring/requirements.txt
```

##### 1.3.2 Downloading kmer lists for primer generation

`kmer` files are required in order to use the Primer3 masker function, which prevents the design of primers using repetitive template regions.

The `kmer` files are needed in the project root directory to run the `primer` command and to run the unit tests. 
(see [Primer3 Manual](https://primer3.org/manual.html#PRIMER_MASK_KMERLIST_PATH)).

1. Set config parameters (in your user-defined designer config file, see [below](
#2212-designer-config))
2. Provide 2 files with kmers: `homo_sapiens_11.list` and `homo_sapiens_16.list`

These kmer lists can be downloaded using:

```sh
chmod +x ./download_kmer_lists.sh
./download_kmer_lists.sh
```

##### 1.3.3 Make the Designer Script Executable and Check Version
```sh
chmod +x ./designer.sh
```

Check the Designer Version:
```sh
./designer.sh version
```

### 1.4 Running unit tests
```sh
python -m unittest discover --start-directory ./tests --top-level-directory .
```

## 2. Usage

### 2.1 Quick start with primer command
To run the `primer` command with minimal configuration, you can follow the steps below to get started. 
In the project directory there are example files you can use to run the `primer` command. 
For example, you will need the test `fasta` file from `./examples` as follows:

```
./designer.sh primer --fasta ./examples/test_example_slice.fa
```

By running the above command, you will see output files are generated in the `td_output` folder in the necessary file formats. 

### 2.2 Primer Designer Tool

#### 2.2.1 Primer3 Runner

Primer3 runner uses 2 types of config:
- Primer3 parameters config
- Designer config

##### 2.2.1.1 Primer3 parameters config

This defines the configuration settings for Primer3 (see the [Primer3 Manual](https://primer3.org/manual.html) for more details).

The default configuration can be found in `config/default_primer3.config.json` (which should NOT be moved, deleted or edited). 
This file contains the default configuration that will be applied if no user config file is provided.

To provide your own config parameters, please copy the default file into a new file, rename it, and edit it. 
You can pass your own user Primer3 config file to the `primer` command using the  `--primer3_params` argument indicating the file path.

##### 2.2.1.2 Designer config

This specifies parameters specific to the Primer Designer tool. You can specify:
1. A vector of different stringencies to be applied when running Primer3,
2. Any filters that should be applied to the list of primer pairs provided by Primer3 (see [below](
#2214-applying-filters-from-the-designer-config-file)),
3. The order in which primer pairs should be ranked in the output CSV files (see [below](
#2215-applying-ranking-from-the-designer-config-file)),
4. The column order for the output CSV files.

Note that, in the Designer config, `stringency_vector` corresponds to 'PRIMER_MASK_FAILURE_RATE' in the [Primer3 Manual](https://primer3.org/manual.html#PRIMER_MASK_FAILURE_RATE). 
This means that a value of 0.1 will apply more stringent settings for the masking algorithm than a value of 1.

The default configuration can be found in `config/default_designer.config.json` (which should NOT be moved, deleted or edited). 
This file contains the default configuration that will be applied if no user config file is provided.

 To provide your own config parameters, please copy the default file into a new file, rename it, and edit it. 
 You can pass your own user designer config file to the `primer` command using the  `--conf` argument indicating the file path. 
 For any fields missing in the user-defined config file, the default settings from `config/default_designer.config.json` will be applied.

##### 2.2.1.3 Running Primer3

```sh
./designer.sh primer [--fasta SLICE_FASTA] [--dir OUTPUT_FOLDER] [--primer3_params PRIMER_CONFIG_JSON] [--conf DESIGNER_CONFIG_JSON]
```

Example command:
```sh
./designer.sh primer --fasta slice.fa --dir p3_output
```

##### 2.2.1.4 Applying filters from the designer config file

You can set your own filtering parameters using your user designer config file (see [above](
#2212-designer-config)).
 * The `duplicates` filter will discard any duplicated primer pairs that have an equivalent pair with a lower primer mask failure rate (see [above](#2212-designer-config)).
 * The `HAP1_variant` filter will discard all primer pairs with at least one primer containing SNPs (variants) that differ between the HAP1 genome and the GRCh38 reference genome.
 These filters can be turned on (`true`) or off (`false`) as follows:

```
{
  ...
  "filters": {
    "duplicates": true,
    "HAP1_variant": true
  },
  ...
}
```

Remember to use the exact names mentioned above. If a filter name is missing, it will not be applied (e.g., if only `"filters": {"duplicates": true}`, `"HAP1_variant"` will not be applied).

If no user config file is passed, then the default `config/default_designer.config.json` will be applied. 
If a user config file is passed, but it does not contain a `filters` key, then the filters from `config/default_designer.config.json` will be applied. 
If the user config file contains a `filters` key and no filters are defined, i.e. `"filters": {},`, then the `duplicates` filter will be applied by default.

##### 2.2.1.5 Applying ranking from the designer config file

You can set your own ranking parameters using your user designer config file (see [above](
#2212-designer-config)). Ranking parameters can be turned on (`true`) or off (`false`) as follows:

```
{
  ...
  "ranking": {
    "stringency": true,
    "product_size": true
  },
  ...
}
```

The order specified in the config file will be retained for ranking. 
In this example, ranking will be applied first by the primer mask failure rate value (in ascending order) 
and then by product size (in descending order), i.e., primer pairs with the same primer mask failure rate value will be ranked according to their product size.

If a name is missing, it will not be used for ranking. If no user config file is passed, then the default `config/default_designer.config.json` will be applied. 
If a user config file is passed, but it does not contain a `ranking` key, then the ranking parameters from `config/default_designer.config.json` will be applied. 
If the user config file contains a `ranking` key and no ranking is defined, i.e. `"ranking": {},`, no ranking will be applied.

##### 2.2.1.6 Specifying column order through the designer config file

Column order can be specified through the user designer config file:

```
{
  ...
  "csv_column_order": ["primer_type",
                        "primer",
                        "penalty",
                        "stringency",
                        "sequence",
                        "primer_start",
                        "primer_end",
                        "tm",
                        "gc_percent",
                        "self_any_th",
                        "self_end_th",
                        "hairpin_th",
                        "end_stability",
                        "chromosome",
                        "pre_targeton_start",
                        "pre_targeton_end",
                        "product_size"]
}
```

All available columns are indicated in the example above. Note that any columns with names missing in the user designer config file will not be present in the output CSV files.

##### 2.2.1.7 Using the designer config file to set command-line arguments

Command-line arguments (`--dir`, `--fasta`, and `--primer3_params`) can, alternatively, be specified in the user designer config file. 
Note that, as no command line arguments are specified in `config/default_designer.config.json`, 
if no user designer config file is provided, `--fasta` (mandatory parameter), `--dir` and `--primer3_params` 
(both optional) will have to be passed as command-line arguments.

```
{
  "stringency_vector": [...],
  "filters": {...}
  "ranking": {...},
  "csv_column_order": [...],
  "dir": "/path/to/output/directory",
  "fasta": "/path/to/fasta/file.fa",
  "primer3_params": "/path/to/primer3/configuration/file"
}
```
Once you add the above configuration to `custom_config.json` file, you will be able to run the following commands:

```
./designer.sh design --conf custom_config.json
```
**or**

```
./designer.sh primer --conf custom_config.json
```

**Note:** Where these arguments are specified both in the command line and in the user designer config file, the parameters specified in the command line will take precedence.

## 3. File formats

File formats for data input and output files used by the Primer Designer Tool are listed and described with examples below.

### 3.1 Primer3 and Designer Fasta Input File (Slicer Fasta output)
Contains the slice sequence, with its ID, coordinates and strand in the header.
If multiple slices are provided in this file, only the first slice will be processed. 
The remaining slices will be ignored, and the user will be notified.

```
>STEQ::7:44490254-44490755(+)
CCGCGCTTCAAATTACTGAAGCCATTCTCACAAGCTCAACCCCAGGACACCAGGAAAAGGAGGAAACAGGCTGGGAGAGCTTGGAGGAGCGGGCGCCAGGAGTCAGGGCAGGCCGGGGGCGGCGGCTCCCCAGACCGCAGGCCCGCCCGCCTCACCTGCAGCACCAAGGCCTGCGCGGCCCCGGGCGGGAAGCCCATGCGGTCCGATGGGTGGCGCAGCAAGCGATAGAAGTCTGTCTTGCGGTAGAGGAAGCCAAAGAGAACGCGCAGGAAATCCTGGACGTTGCCCACGTGCTGCAGGATGCCCAAAAGGGCCTGGTCATACAGCTCGGCCGCCCCTGTCTCCATGTCGCCTCCCGCCCTAGGTACGCTTCACACACACAGCGCCGCCTCAGACCTGCCGACTGGCCACTTCCGGCGTCCGCAGCCAACGGCTCCGCCGGACGGCGCGGCTGCGGAACTTCCGGTCCGTGCTCGTTCGGCCCCCGCGGCCCCGGGCTGTT
```

### 3.2 Primer3 Output BED file
Genomic locations of the primers and their names. Names are incremented from 0 and given F and R depending on whether they reflect the positive (5'-3') or negative (3'-5') strand, respectively. 
See example showing the first 10 rows of the file.

Note that all sample output files from the primer designer command (i.e., this BED file and the two following CSV files) 
have been produced using the settings from the default config files and the FASTA file above.

| chrom | chromStart | chromEnd | name           | score | strand |
|-------|------------|----------|----------------|-------|--------|
| 7     | 44490309   | 44490328 | STEQ_LibAmpF_0 | 0     | +      |
| 7     | 44490404   | 44490422 | STEQ_LibAmpR_0 | 0     | -      |
| 7     | 44490309   | 44490328 | STEQ_LibAmpF_1 | 0     | +      |
| 7     | 44490403   | 44490421 | STEQ_LibAmpR_1 | 0     | -      |
| 7     | 44490293   | 44490312 | STEQ_LibAmpF_0 | 0     | +      |
| 7     | 44490404   | 44490422 | STEQ_LibAmpR_0 | 0     | -      |
| 7     | 44490293   | 44490312 | STEQ_LibAmpF_1 | 0     | +      |
| 7     | 44490403   | 44490421 | STEQ_LibAmpR_1 | 0     | -      |
| 7     | 44490294   | 44490313 | STEQ_LibAmpF_2 | 0     | +      |
| 7     | 44490404   | 44490422 | STEQ_LibAmpR_2 | 0     | -      |
| 7     | 44490294   | 44490313 | STEQ_LibAmpF_3 | 0     | +      |

Raw file (`p3_output.bed`)
```
7       44490309        44490328        STEQ_LibAmpF_0  0       +
7       44490404        44490422        STEQ_LibAmpR_0  0       -
7       44490309        44490328        STEQ_LibAmpF_1  0       +
7       44490403        44490421        STEQ_LibAmpR_1  0       -
7       44490293        44490312        STEQ_LibAmpF_0  0       +
7       44490404        44490422        STEQ_LibAmpR_0  0       -
7       44490293        44490312        STEQ_LibAmpF_1  0       +
7       44490403        44490421        STEQ_LibAmpR_1  0       -
7       44490294        44490313        STEQ_LibAmpF_2  0       +
7       44490404        44490422        STEQ_LibAmpR_2  0       -
7       44490294        44490313        STEQ_LibAmpF_3  0       +
```

### 3.3 Primer3 Output CSV file
It contains all the additional information from Primer3 for the individual primers. 
Column order can be specified through the Designer tool config. 
In this example, primer pairs have been ranked first by `stringency` and then by `product_size`, according to the default config file. 
See example showing the first 10 rows of the file (plus, in the table, the first 4 rows with a different stringency).

| primer_type | primer          | penalty | stringency | sequence             | primer_start | primer_end | tm     | gc_percent | self_any_th | self_end_th | hairpin_th | end_stability | chromosome | pre_targeton_start | pre_targeton_end | product_size | targeton_id | pair_uid                             |
|-------------|-----------------|---------|------------|----------------------|--------------|------------|--------|------------|-------------|-------------|------------|---------------|------------|--------------------|------------------|--------------|-------------|--------------------------------------|
| LibAmp      | STEQ_LibAmpF_18 | 0.894   | 0.1        | TTCTCACAAGCTCAACCCCA | 44490278     | 44490297   | 59.157 | 50.0       | 0.0         | 0.0         | 0.0        | 4.96          | 7          | 44490254           | 44490755         | 145          | STEQ        | 8fb6334a-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpR_18 | 2.056   | 0.1        | CCTTGGTGCTGCAGGTGAG  | 44490404     | 44490422   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.51          | 7          | 44490254           | 44490755         | 145          | STEQ        | 8fb6334a-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpF_6  | 0.455   | 0.1        | TCTCACAAGCTCAACCCCAG | 44490279     | 44490298   | 59.602 | 55.0       | 0.0         | 0.0         | 0.0        | 4.45          | 7          | 44490254           | 44490755         | 144          | STEQ        | 8fb6333e-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpR_6  | 2.056   | 0.1        | CCTTGGTGCTGCAGGTGAG  | 44490404     | 44490422   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.51          | 7          | 44490254           | 44490755         | 144          | STEQ        | 8fb6333e-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpF_10 | 0.455   | 0.1        | TCTCACAAGCTCAACCCCAG | 44490279     | 44490298   | 59.602 | 55.0       | 0.0         | 0.0         | 0.0        | 4.45          | 7          | 44490254           | 44490755         | 144          | STEQ        | 8fb63342-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpR_10 | 2.172   | 0.1        | CCTTGGTGCTGCAGGTGA   | 44490405     | 44490422   | 59.886 | 61.111     | 17.124      | 0.0         | 0.0        | 4.02          | 7          | 44490254           | 44490755         | 144          | STEQ        | 8fb63342-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpF_19 | 0.894   | 0.1        | TTCTCACAAGCTCAACCCCA | 44490278     | 44490297   | 59.157 | 50.0       | 0.0         | 0.0         | 0.0        | 4.96          | 7          | 44490254           | 44490755         | 144          | STEQ        | 8fb6334b-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpR_19 | 2.056   | 0.1        | CTTGGTGCTGCAGGTGAGG  | 44490403     | 44490421   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.86          | 7          | 44490254           | 44490755         | 144          | STEQ        | 8fb6334b-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpF_7  | 0.455   | 0.1        | TCTCACAAGCTCAACCCCAG | 44490279     | 44490298   | 59.602 | 55.0       | 0.0         | 0.0         | 0.0        | 4.45          | 7          | 44490254           | 44490755         | 143          | STEQ        | 8fb6333f-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpR_7  | 2.056   | 0.1        | CTTGGTGCTGCAGGTGAGG  | 44490403     | 44490421   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.86          | 7          | 44490254           | 44490755         | 143          | STEQ        | 8fb6333f-443c-11ef-91cc-fa163e1eb62c |
| ...         | ...             | ...     | ...        | ...                  | ...          | ...        | ...    | ...        | ...         | ...         | ...        | ...           | ...        | ...                | ...              | ...          | ...         | ...                                  |
| LibAmp      | STEQ_LibAmpF_0  | 0.22    | 0.5        | AAAGGAGGAAACAGGCTGGG | 44490309     | 44490328   | 59.887 | 55.0       | 20.341      | 0.0         | 0.0        | 4.45          | 7          | 44490254           | 44490755         | 114          | STEQ        | 8fb63324-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpR_0  | 2.056   | 0.5        | CCTTGGTGCTGCAGGTGAG  | 44490404     | 44490422   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.51          | 7          | 44490254           | 44490755         | 114          | STEQ        | 8fb63324-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpF_6  | 0.22    | 0.5        | AAAGGAGGAAACAGGCTGGG | 44490309     | 44490328   | 59.887 | 55.0       | 20.341      | 0.0         | 0.0        | 4.45          | 7          | 44490254           | 44490755         | 114          | STEQ        | 8fb6332a-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpR_6  | 2.172   | 0.5        | CCTTGGTGCTGCAGGTGA   | 44490405     | 44490422   | 59.886 | 61.111     | 17.124      | 0.0         | 0.0        | 4.02          | 7          | 44490254           | 44490755         | 114          | STEQ        | 8fb6332a-443c-11ef-91cc-fa163e1eb62c |
| ...         | ...             | ...     | ...        | ...                  | ...          | ...        | ...    | ...        | ...         | ...         | ...        | ...           | ...        | ...                | ...              | ...          | ...         | ...                                  |


Raw File (`p3_output.csv`)
```
primer_type,primer,penalty,stringency,sequence,primer_start,primer_end,tm,gc_percent,self_any_th,self_end_th,hairpin_th,end_stability,chromosome,pre_targeton_start,pre_targeton_end,product_size,targeton_id,pair_uid
LibAmp,STEQ_LibAmpF_18,0.894,0.1,TTCTCACAAGCTCAACCCCA,44490278,44490297,59.157,50.0,0.0,0.0,0.0,4.96,7,44490254,44490755,145,STEQ,8fb6334a-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpR_18,2.056,0.1,CCTTGGTGCTGCAGGTGAG,44490404,44490422,60.97,63.158,17.124,0.0,0.0,3.51,7,44490254,44490755,145,STEQ,8fb6334a-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpF_6,0.455,0.1,TCTCACAAGCTCAACCCCAG,44490279,44490298,59.602,55.0,0.0,0.0,0.0,4.45,7,44490254,44490755,144,STEQ,8fb6333e-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpR_6,2.056,0.1,CCTTGGTGCTGCAGGTGAG,44490404,44490422,60.97,63.158,17.124,0.0,0.0,3.51,7,44490254,44490755,144,STEQ,8fb6333e-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpF_10,0.455,0.1,TCTCACAAGCTCAACCCCAG,44490279,44490298,59.602,55.0,0.0,0.0,0.0,4.45,7,44490254,44490755,144,STEQ,8fb63342-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpR_10,2.172,0.1,CCTTGGTGCTGCAGGTGA,44490405,44490422,59.886,61.111,17.124,0.0,0.0,4.02,7,44490254,44490755,144,STEQ,8fb63342-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpF_19,0.894,0.1,TTCTCACAAGCTCAACCCCA,44490278,44490297,59.157,50.0,0.0,0.0,0.0,4.96,7,44490254,44490755,144,STEQ,8fb6334b-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpR_19,2.056,0.1,CTTGGTGCTGCAGGTGAGG,44490403,44490421,60.97,63.158,17.124,0.0,0.0,3.86,7,44490254,44490755,144,STEQ,8fb6334b-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpF_7,0.455,0.1,TCTCACAAGCTCAACCCCAG,44490279,44490298,59.602,55.0,0.0,0.0,0.0,4.45,7,44490254,44490755,143,STEQ,8fb6333f-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpR_7,2.056,0.1,CTTGGTGCTGCAGGTGAGG,44490403,44490421,60.97,63.158,17.124,0.0,0.0,3.86,7,44490254,44490755,143,STEQ,8fb6333f-443c-11ef-91cc-fa163e1eb62c
```

### 3.4 Primer3 Output Optimal Primer Pairs CSV file
It contains the top 3 optimal primer pairs from the previous CSV file (`p3_output.csv`).

| primer_type | primer          | penalty | stringency | sequence             | primer_start | primer_end | tm     | gc_percent | self_any_th | self_end_th | hairpin_th | end_stability | chromosome | pre_targeton_start | pre_targeton_end | product_size | targeton_id | pair_uid                             |
|-------------|-----------------|---------|------------|----------------------|--------------|------------|--------|------------|-------------|-------------|------------|---------------|------------|--------------------|------------------|--------------|-------------|--------------------------------------|
| LibAmp      | STEQ_LibAmpF_18 | 0.894   | 0.1        | TTCTCACAAGCTCAACCCCA | 44490278     | 44490297   | 59.157 | 50.0       | 0.0         | 0.0         | 0.0        | 4.96          | 7          | 44490254           | 44490755         | 145          | STEQ        | 8fb6334a-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpR_18 | 2.056   | 0.1        | CCTTGGTGCTGCAGGTGAG  | 44490404     | 44490422   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.51          | 7          | 44490254           | 44490755         | 145          | STEQ        | 8fb6334a-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpF_6  | 0.455   | 0.1        | TCTCACAAGCTCAACCCCAG | 44490279     | 44490298   | 59.602 | 55.0       | 0.0         | 0.0         | 0.0        | 4.45          | 7          | 44490254           | 44490755         | 144          | STEQ        | 8fb6333e-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpR_6  | 2.056   | 0.1        | CCTTGGTGCTGCAGGTGAG  | 44490404     | 44490422   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.51          | 7          | 44490254           | 44490755         | 144          | STEQ        | 8fb6333e-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpF_10 | 0.455   | 0.1        | TCTCACAAGCTCAACCCCAG | 44490279     | 44490298   | 59.602 | 55.0       | 0.0         | 0.0         | 0.0        | 4.45          | 7          | 44490254           | 44490755         | 144          | STEQ        | 8fb63342-443c-11ef-91cc-fa163e1eb62c |
| LibAmp      | STEQ_LibAmpR_10 | 2.172   | 0.1        | CCTTGGTGCTGCAGGTGA   | 44490405     | 44490422   | 59.886 | 61.111     | 17.124      | 0.0         | 0.0        | 4.02          | 7          | 44490254           | 44490755         | 144          | STEQ        | 8fb63342-443c-11ef-91cc-fa163e1eb62c |


Raw File (`optimal_primer_pairs.csv`)
```
primer_type,primer,penalty,stringency,sequence,primer_start,primer_end,tm,gc_percent,self_any_th,self_end_th,hairpin_th,end_stability,chromosome,pre_targeton_start,pre_targeton_end,product_size,targeton_id,pair_uid
LibAmp,STEQ_LibAmpF_18,0.894,0.1,TTCTCACAAGCTCAACCCCA,44490278,44490297,59.157,50.0,0.0,0.0,0.0,4.96,7,44490254,44490755,145,STEQ,8fb6334a-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpR_18,2.056,0.1,CCTTGGTGCTGCAGGTGAG,44490404,44490422,60.97,63.158,17.124,0.0,0.0,3.51,7,44490254,44490755,145,STEQ,8fb6334a-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpF_6,0.455,0.1,TCTCACAAGCTCAACCCCAG,44490279,44490298,59.602,55.0,0.0,0.0,0.0,4.45,7,44490254,44490755,144,STEQ,8fb6333e-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpR_6,2.056,0.1,CCTTGGTGCTGCAGGTGAG,44490404,44490422,60.97,63.158,17.124,0.0,0.0,3.51,7,44490254,44490755,144,STEQ,8fb6333e-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpF_10,0.455,0.1,TCTCACAAGCTCAACCCCAG,44490279,44490298,59.602,55.0,0.0,0.0,0.0,4.45,7,44490254,44490755,144,STEQ,8fb63342-443c-11ef-91cc-fa163e1eb62c
LibAmp,STEQ_LibAmpR_10,2.172,0.1,CCTTGGTGCTGCAGGTGA,44490405,44490422,59.886,61.111,17.124,0.0,0.0,4.02,7,44490254,44490755,144,STEQ,8fb63342-443c-11ef-91cc-fa163e1eb62c
```

### 3.5 Primer3 Output Discarded Primer Pairs CSV file
It contains primer pairs discarded during filtering, with discard reason in column discard_reason (last column). 
The column order is specified through the Designer tool config. See example showing the first 10 rows of the file.

| primer_type | primer         | penalty | stringency | sequence             | primer_start | primer_end | tm     | gc_percent | self_any_th | self_end_th | hairpin_th | end_stability | chromosome | pre_targeton_start | pre_targeton_end | product_size | targeton_id | pair_uid                             | discard reason                         |
|-------------|----------------|---------|------------|----------------------|--------------|------------|--------|------------|-------------|-------------|------------|---------------|------------|--------------------|------------------|--------------|-------------|--------------------------------------|----------------------------------------|
| LibAmp      | STEQ_LibAmpF_0 | 0.22    | 1.0        | AAAGGAGGAAACAGGCTGGG | 44490309     | 44490328   | 59.887 | 55.0       | 20.341      | 0.0         | 0.0        | 4.45          | 7          | 44490254           | 44490755         | 114          | STEQ        | 8fb63310-443c-11ef-91cc-fa163e1eb62c | has duplicate with a higher stringency |
| LibAmp      | STEQ_LibAmpR_0 | 2.056   | 1.0        | CCTTGGTGCTGCAGGTGAG  | 44490404     | 44490422   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.51          | 7          | 44490254           | 44490755         | 114          | STEQ        | 8fb63310-443c-11ef-91cc-fa163e1eb62c | has duplicate with a higher stringency |
| LibAmp      | STEQ_LibAmpF_1 | 0.22    | 1.0        | AAAGGAGGAAACAGGCTGGG | 44490309     | 44490328   | 59.887 | 55.0       | 20.341      | 0.0         | 0.0        | 4.45          | 7          | 44490254           | 44490755         | 113          | STEQ        | 8fb63311-443c-11ef-91cc-fa163e1eb62c | has duplicate with a higher stringency |
| LibAmp      | STEQ_LibAmpR_1 | 2.056   | 1.0        | CTTGGTGCTGCAGGTGAGG  | 44490403     | 44490421   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.86          | 7          | 44490254           | 44490755         | 113          | STEQ        | 8fb63311-443c-11ef-91cc-fa163e1eb62c | has duplicate with a higher stringency |
| LibAmp      | STEQ_LibAmpF_2 | 0.328   | 1.0        | CCCCAGGACACCAGGAAAAG | 44490293     | 44490312   | 60.251 | 60.0       | 0.0         | 0.0         | 31.239     | 2.27          | 7          | 44490254           | 44490755         | 130          | STEQ        | 8fb63312-443c-11ef-91cc-fa163e1eb62c | has duplicate with a higher stringency |
| LibAmp      | STEQ_LibAmpR_2 | 2.056   | 1.0        | CCTTGGTGCTGCAGGTGAG  | 44490404     | 44490422   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.51          | 7          | 44490254           | 44490755         | 130          | STEQ        | 8fb63312-443c-11ef-91cc-fa163e1eb62c | has duplicate with a higher stringency |
| LibAmp      | STEQ_LibAmpF_2 | 0.328   | 0.5        | CCCCAGGACACCAGGAAAAG | 44490293     | 44490312   | 60.251 | 60.0       | 0.0         | 0.0         | 31.239     | 2.27          | 7          | 44490254           | 44490755         | 130          | STEQ        | 8fb63326-443c-11ef-91cc-fa163e1eb62c | has duplicate with a higher stringency |
| LibAmp      | STEQ_LibAmpR_2 | 2.056   | 0.5        | CCTTGGTGCTGCAGGTGAG  | 44490404     | 44490422   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.51          | 7          | 44490254           | 44490755         | 130          | STEQ        | 8fb63326-443c-11ef-91cc-fa163e1eb62c | has duplicate with a higher stringency |
| LibAmp      | STEQ_LibAmpF_3 | 0.328   | 1.0        | CCCCAGGACACCAGGAAAAG | 44490293     | 44490312   | 60.251 | 60.0       | 0.0         | 0.0         | 31.239     | 2.27          | 7          | 44490254           | 44490755         | 129          | STEQ        | 8fb63313-443c-11ef-91cc-fa163e1eb62c | has duplicate with a higher stringency |
| LibAmp      | STEQ_LibAmpR_3 | 2.056   | 1.0        | CTTGGTGCTGCAGGTGAGG  | 44490403     | 44490421   | 60.97  | 63.158     | 17.124      | 0.0         | 0.0        | 3.86          | 7          | 44490254           | 44490755         | 129          | STEQ        | 8fb63313-443c-11ef-91cc-fa163e1eb62c | has duplicate with a higher stringency |


Raw File (`discarded_pairs.csv`)
```
primer_type,primer,penalty,stringency,sequence,primer_start,primer_end,tm,gc_percent,self_any_th,self_end_th,hairpin_th,end_stability,chromosome,pre_targeton_start,pre_targeton_end,product_size,targeton_id,pair_uid,discard_reason
t,pre_targeton_end,product_size,targeton_id,pair_uid,discard_reason
LibAmp,STEQ_LibAmpF_0,0.22,1.0,AAAGGAGGAAACAGGCTGGG,44490309,44490328,59.887,55.0,20.341,0.0,0.0,4.45,7,44490254,44490755,114,STEQ,8fb63310-443c-11ef-91cc-fa163e1eb62c,has duplicate with a higher stringency
LibAmp,STEQ_LibAmpR_0,2.056,1.0,CCTTGGTGCTGCAGGTGAG,44490404,44490422,60.97,63.158,17.124,0.0,0.0,3.51,7,44490254,44490755,114,STEQ,8fb63310-443c-11ef-91cc-fa163e1eb62c,has duplicate with a higher stringency
LibAmp,STEQ_LibAmpF_1,0.22,1.0,AAAGGAGGAAACAGGCTGGG,44490309,44490328,59.887,55.0,20.341,0.0,0.0,4.45,7,44490254,44490755,113,STEQ,8fb63311-443c-11ef-91cc-fa163e1eb62c,has duplicate with a higher stringency
LibAmp,STEQ_LibAmpR_1,2.056,1.0,CTTGGTGCTGCAGGTGAGG,44490403,44490421,60.97,63.158,17.124,0.0,0.0,3.86,7,44490254,44490755,113,STEQ,8fb63311-443c-11ef-91cc-fa163e1eb62c,has duplicate with a higher stringency
LibAmp,STEQ_LibAmpF_2,0.328,1.0,CCCCAGGACACCAGGAAAAG,44490293,44490312,60.251,60.0,0.0,0.0,31.239,2.27,7,44490254,44490755,130,STEQ,8fb63312-443c-11ef-91cc-fa163e1eb62c,has duplicate with a higher stringency
LibAmp,STEQ_LibAmpR_2,2.056,1.0,CCTTGGTGCTGCAGGTGAG,44490404,44490422,60.97,63.158,17.124,0.0,0.0,3.51,7,44490254,44490755,130,STEQ,8fb63312-443c-11ef-91cc-fa163e1eb62c,has duplicate with a higher stringency
LibAmp,STEQ_LibAmpF_2,0.328,0.5,CCCCAGGACACCAGGAAAAG,44490293,44490312,60.251,60.0,0.0,0.0,31.239,2.27,7,44490254,44490755,130,STEQ,8fb63326-443c-11ef-91cc-fa163e1eb62c,has duplicate with a higher stringency
LibAmp,STEQ_LibAmpR_2,2.056,0.5,CCTTGGTGCTGCAGGTGAG,44490404,44490422,60.97,63.158,17.124,0.0,0.0,3.51,7,44490254,44490755,130,STEQ,8fb63326-443c-11ef-91cc-fa163e1eb62c,has duplicate with a higher stringency
LibAmp,STEQ_LibAmpF_3,0.328,1.0,CCCCAGGACACCAGGAAAAG,44490293,44490312,60.251,60.0,0.0,0.0,31.239,2.27,7,44490254,44490755,129,STEQ,8fb63313-443c-11ef-91cc-fa163e1eb62c,has duplicate with a higher stringency
LibAmp,STEQ_LibAmpR_3,2.056,1.0,CTTGGTGCTGCAGGTGAGG,44490403,44490421,60.97,63.158,17.124,0.0,0.0,3.86,7,44490254,44490755,129,STEQ,8fb63313-443c-11ef-91cc-fa163e1eb62c,has duplicate with a higher stringency
```

### 3.6 Genomic Reference file
A fasta file of latest GRCh38 genome. This is used for gathering the slice sequences and retrieving primer information. 
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


\
\
\
The following file formats are generated by tools and commands that are no longer in use.
You can refer to the [Tools and commands no longer in use](#5-tools-and-commands-no-longer-in-use) section for more details.

### 3.7 Slicer Input BED File
A BED file containing the regions you wish to slice across. 

The chromosome column data must match your reference fasta file IDs. If your reference had >chr1 then you must call chromosome 1 'chr1' in this column and vice-versa.

Note: BED effectively are applied tsv files so use tabs to separate the values. 
Headers are optional in the BED file and can be a cause of issues if they aren't perfect. 
Strand is required for the slicer to ensure sequences are output in the correct orientation. 
Score isn't used but the field must be present for the file format to be read correctly.

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

### 3.8 Slicer BED output
BED file output with a row for each slice. This file will also be used for running VaLiAnt.

| chrom | chromStart | chromEnd | name                  | score | strand |
|-------|------------|----------|-----------------------|-------|--------|
| 1     | 42930996   | 42931206 | ENSE00003571441_HG6_1 | 0     | -      | 
| 1     | 42931001   | 42931211 | ENSE00003571441_HG6_2 | 0     | -      | 
| 1     | 42931006   | 42931216 | ENSE00003571441_HG6_3 | 0     | -      |
| 1     | 42931011   | 42931221 | ENSE00003571441_HG6_4 | 0     | -      |
| 1     | 42931016   | 42931226 | ENSE00003571441_HG6_5 | 0     | -      |

Raw file
```
1	42930996	42931206	ENSE00003571441_HG6_1	0	-
1	42931001	42931211	ENSE00003571441_HG6_2	0	-
1	42931006	42931216	ENSE00003571441_HG6_3	0	-
1	42931011	42931221	ENSE00003571441_HG6_4	0	-
1	42931016	42931226	ENSE00003571441_HG6_5	0	-
```

### 3.9 Scoring Tool Input iPCRess file example
Two files, stnd and err.
stnd: Space separated text file. Sequence_id contains the chromosome and description can be 'forward', 'revcomp', 'single_A' or 'single_B'.

| sequence_id         | experiment_id | product_length | primer_5 | pos_5    | mismatch_5 | primer_3 | pos_3    | mismatch_3 | description |
|---------------------|---------------|----------------|----------|----------|------------|----------|----------|------------|-------------|
| 19:filter(unmasked) | ID0001        | 259            | A        | 44907726 | 0          | B        | 44907967 | 0          | forward     |

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

### 3.10 Primer Designer Output example
The output generated by the Primer Designer using `design` command is identical to that of the `primer` command. 
This consistency arises because the `design` command workflow follows the same procedures as the `primer` command workflow.

## 4. For Developers

### 4.1 Git Hooks

Located in  .githooks/
Follows standard Git hook methodology: https://git-scm.com/docs/githooks
Currently just "pre-push" that is run on ```git push```
To run either invoke ```make``` or:
```sh
git config core.hooksPath .githooks
chmod +x .githooks/*
```

### 4.2 Python debugger
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

## 5. Tools and commands no longer in use

The following tools and commands are no longer in use because they are no longer part of the main Primer Designer Workflow. 
However, note that they can still be run as separate tools. We are currently not providing maintenance for these tools.

#### 5.1 Designer Workflow (Primer3)

Running full Designer Workflow:
```sh
./designer.sh design [-h] [--fasta SLICE_FASTA] [--primer3_params PRIMER_CONFIG_JSON]
```

Example Command
```sh
./designer.sh design --fasta examples/fasta_example.fa
```

### 5.2 Primer Scoring Tool

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

### 5.3 Slicer Tool

Running Slicer tool:
```sh
./designer.sh slicer [-h] [-f5 FLANK_5] [-f3 FLANK_3] [-l LENGTH] [-o OFFSET] [-d DIR] [--bed INPUT_BED] [--fasta REF_FASTA]
```

Example command:
```sh
./designer.sh slicer --bed example.bed --fasta example_genomic_ref.fa -d example_dir
```

### 5.4 Targeton CSV generation

To generate the targeton CSV used in primer scoring:
```sh
./designer.sh generate_targeton_csv [--primers PRIMERS] [--bed BED] [--dir DIR]
```
Example command:
```sh
./designer.sh generate_targeton_csv --primers example_ipcress_input.txt --bed example.bed --dir example_dir
```
Please note primer pair names in the iPCRess input file must be prefixed by the corresponding region name in the BED file.

### 5.5 Primer data collation and output to CSV and JSON (for Benchling)

To collate the primer and scoring data and output to CSV & JSON file:
```sh
./designer.sh collate_primer_data [--p3_csv Primer3_output.csv] [--score_tsv scoring_output.tsv] [--dir DIR]
```

Examples of the output can be found below.

This is also run as part of the design command.

### 5.6 Post primers to Benchling

To post the top 3 primer pairs for each targeton from the Primer Designer JSON output:
```sh
./designer.sh post_primers [--primer_json PRIMER_JSON]
```
Example command:
```sh
./designer.sh post_primers --primer_json primer_designer.json
```
A message will be printed if there are less than 3 primer pairs for a particular targeton. 
Please note that some fields on Benchling will have to be updated manually for now.

## 6. Upcoming releases

Docker image in upcoming feature in later release.
