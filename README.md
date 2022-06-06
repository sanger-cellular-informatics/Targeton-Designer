# Targeton Designer

Standalone targeton designer tool.

[[_TOC_]]

## Installation

Dependencies:

BedTools
```sh
apt-get install bedtools
```

### Python Virtual Environment

Requirements:
 - Python3.8+
 - Python-venv

Install venv if you haven't already got it
```sh
apt install python3.8-venv
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

Install docker-compose - https://docs.docker.com/compose/install/

## Usage

### Python

Running Slicer tool:
```sh
python3 designer/slicer.py bed bedfile fasta fasta_file -l 200
```

Running Python3:
```sh
python3 runner/cmd.py --seqs targeton_file.csv --dir ./tmp_folder --ref genomic_reference_file.fna
```

### Docker

Run Primer3:
docker run --rm Primer3 [-h] [-s SEQ] [-d DIR] [-r REF]
	Required args:
		'--seq' - Path to sequences file
		'--dir' - Output folder location
		'--ref' - Reference file path

Run Slicer:
docker run --rm Slicer [-h] [-f5 FLANK_5] [-f3 FLANK_3] [-l LENGTH] [-o OFFSET]
                 bed fasta

Get sequence slices for regions in BED file according to parameters specified

positional arguments:
  bed                   BED file containing regions of interest
  fasta                 FASTA file to retrieve sequences from

optional arguments:
  -h, --help            show this help message and exit
  -f5 FLANK_5, --flank_5 FLANK_5
                        how far to extend region at 5' end
  -f3 FLANK_3, --flank_3 FLANK_3
                        how far to extend region at 3' end
  -l LENGTH, --length LENGTH
                        length of each slice
  -o OFFSET, --offset OFFSET
                        offset between each slice

## File formats
