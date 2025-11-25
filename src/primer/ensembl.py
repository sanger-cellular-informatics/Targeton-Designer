import time

import requests


def get_seq_from_ensembl_by_coords(chromosome: str, start: int, end: int, strand: str = '+'):
    ens_strand = '-1' if strand == '-' else '1'

    url = (
        f'https://rest.ensembl.org/sequence/region/human/'
        f'{chromosome}:{start}..{end}:{ens_strand}'
    )
    headers = {'Content-type': 'text/plain'}
    response = requests.get(url, headers=headers)
    time.sleep(0.1)
    if response.status_code == requests.codes.ok:
        return response.text
    else:
        raise requests.exceptions.RequestException
