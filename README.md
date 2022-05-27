# Targeton Designer

Standalone targeton designer tool

Install docker-compose - https://docs.docker.com/compose/install/

Running tests:
python3 -m unittest


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