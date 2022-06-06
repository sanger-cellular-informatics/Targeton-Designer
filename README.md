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

Docker Installation:
https://docs.docker.com/get-started/

Docker Install
```sh
sudo apt-get update
sudo apt-get install \
    ca-certificates \
    curl \
    gnupg \
    lsb-releasesudo mkdir -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpgsudo apt install docker.iosudo apt install docker-compose
```

Download ref fasta
```sh
wget http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

Create output folder
```sh
mkdir p3_output
```

Local Install
```sh
git clone https://gitlab.internal.sanger.ac.uk/sci/targeton-designer.git
 
cd targeton_designer
 
pip install -r requirements.txt
```

Run primer3 image
```sh
sudo docker-compose build
```

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

**Run Slicer**
```sh
docker run --rm Slicer [-h] [-f5 FLANK_5] [-f3 FLANK_3] [-l LENGTH] [-o OFFSET] bed fasta
```

Slicer Arguments

Get sequence slices for regions in BED file according to parameters specified

positional arguments:
  bed                   BED file containing regions of interest
  fasta                 FASTA file to retrieve sequences from

Optional Slicer Arguments:
  - -h, --help            show this help message and exit
  - -f5 FLANK_5, --flank_5 FLANK_5
                        how far to extend region at 5' end
  - -f3 FLANK_3, --flank_3 FLANK_3
                        how far to extend region at 3' end
  - -l LENGTH, --length LENGTH
                        length of each slice
  - -o OFFSET, --offset OFFSET
                        offset between each slice

**Run Primer3**
```sh
docker run --rm Primer3 [-h] [-s SEQ] [-d DIR] [-r REF]
```

Required Primer3 args:
- --seq - Path to sequences file
- --dir - Output folder location
- --ref - Reference file path


## File formats

### Slicer Input BED File

### Fasta Input File

### Targeton Slices Input File 
**To be replaced by fa**
```sh
id,seq
ENSE00000893952,TCCACACAGGATGCCAGGCCAAGGTGGAGCAAGCGGTGGAGACAGAGCCGGAGCCCGAGCTGCGCCAGCAGACCGAGTGGCAGAGCGGCCAGCGCTGGGAACTGGCACTGGGTCGCTTTTGGGATTACCTGCGCTGGGTGCAGACACTGTCTGAGCAGGTGCAGGAGGAGCTGCTCAGCTCCCAGGTCACCCAGGAACTGAGGTGAGTGTCC
```

### Genomic Reference File

Either supply a local genome reference file or download one from EnsEMBL:
http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/
