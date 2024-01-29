from collections import defaultdict

import requests

from utils.file_system import parse_json


def parse_primer_json(primer_json: str) -> defaultdict:
    primer_data = parse_json(primer_json)
    targeton_primer_pairs = defaultdict(list)
    for primer_pair in primer_data:
        targeton_primer_pairs[primer_pair['targeton']].append(primer_pair)
    return targeton_primer_pairs


def filter_primer_pairs(targeton_primer_pairs: defaultdict) -> list:
    top_primer_data = []
    for targeton in targeton_primer_pairs.keys():
        primer_pairs = targeton_primer_pairs[targeton]
        if len(primer_pairs) < 3:
            print(f'Only {len(primer_pairs)} primer pair(s) for targeton: {targeton}')
        primer_pairs.sort(key=lambda pair: pair['score'])
        top_primer_pairs = primer_pairs[:3]
        top_primer_data.extend(top_primer_pairs)
    return top_primer_data


def post_primer_data(primer_data: list) -> None:
    response = requests.post(
        'https://sge-service.link:8081/libamp', json=primer_data, headers={'Content-Type': 'application/json'}
    )
    if response.status_code == 201:
        print(f'Successfully posted primers!')
    else:
        print(f'Issue with post request: {response.status_code} {response.reason}')


def post_primer_pairs(primer_json: str) -> None:
    targeton_primer_pairs = parse_primer_json(primer_json)
    top_primer_pairs = filter_primer_pairs(targeton_primer_pairs)
    post_primer_data(top_primer_pairs)
