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

Running Python3:
```sh
python3 runner/cmd.py --seqs targeton_file.csv --dir ./tmp_folder --ref genomic_reference_file.fna
```

### Docker

Upcoming feature in later releases

## File formats

### Slicer Input BED File

### Primer3 Fasta Input File

### Targeton Slices Input File 
**To be replaced by fa**
```sh
id,seq
ENSE00000893952,TCCACACAGGATGCCAGGCCAAGGTGGAGCAAGCGGTGGAGACAGAGCCGGAGCCCGAGCTGCGCCAGCAGACCGAGTGGCAGAGCGGCCAGCGCTGGGAACTGGCACTGGGTCGCTTTTGGGATTACCTGCGCTGGGTGCAGACACTGTCTGAGCAGGTGCAGGAGGAGCTGCTCAGCTCCCAGGTCACCCAGGAACTGAGGTGAGTGTCC
```

### Genomic Reference File

Either supply a local genome reference file or download one from EnsEMBL:
http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/
