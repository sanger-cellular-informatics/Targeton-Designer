import json

import requests


def contain_variant(chromosome: str, start: int, end: int) -> bool:
    url = f"https://z4ell7ogh5.execute-api.eu-west-2.amazonaws.com/prod?chromosome={chromosome}&start={start}&end={end}"
    headers = {'Content-type': 'application/json'}

    response = requests.get(url, headers=headers)

    if response.status_code == requests.codes.ok:
        variants_found = json.loads(response.text)["variants"]

        return len(variants_found) > 0
    else:
        raise requests.exceptions.RequestException
