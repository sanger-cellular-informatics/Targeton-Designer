#!/usr/bin/env python3
import designer.slicer as slicer
import runner.primer3_runner as primer3

import sys

def create_slicer_params(input_params):
    return  {
        'input_fasta': input_params['input_fasta'],
        'input_bed': input_params['input_bed'],
        'output_fasta': input_params['output_fasta'],
        'output_bed': input_params['output_bed'],
        'flank_5': input_params['flank_5'],
        'flank_3': input_params['flank_3'],
        'length': input_params['length'],
        'offset': input_params['offset'],
    }

def create_primer3_params(input_params):
    return  {
        'seq': input_params['output_fasta'],
        'bed': input_params['output_bed'],
        'dir': input_params['dir'],
    }

def main(params):
    try:
        slicer_params = create_slicer_params(params)
        slicer.main(slicer_params)

    except Exception:
        print('Slicer failed')
        return ''

    try:
        primer3_params = create_primer3_params(params)
        primer3.main(primer3_params)
    except Exception:
        print('Primer3 failed')
        return ''

    print('Designed successfully')
    return ''

if __name__ == '__main__':
    parsed_args = slicer.parse_args(sys.argv[1:])
    main(vars(parsed_args))
