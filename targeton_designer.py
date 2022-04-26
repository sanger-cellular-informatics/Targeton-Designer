import redis
from flask import Flask

app = Flask(__name__)
cache = redis.Redis(host='redis', port=6379)

from runner.primer3 import primer3_runner
from runner.exonerate import generate_ipcress_input 

@app.route('/')
def start():
    design_input = {
        'SEQUENCE_ID': 'ENSE00000893952',
        'SEQUENCE_TEMPLATE': 'TCCACACAGGATGCCAGGCCAAGGTGGAGCAAGCGGTGGAGACAGAGCCGGAGCCCGAGCTGCGCCAGCAGACCGAGTGGCAGAGCGGCCAGCGCTGGGAACTGGCACTGGGTCGCTTTTGGGATTACCTGCGCTGGGTGCAGACACTGTCTGAGCAGGTGCAGGAGGAGCTGCTCAGCTCCCAGGTCACCCAGGAACTGAGGTGAGTGTCC',
    }

    primers = primer3_runner(design_input)
    generate_ipcress_input(primers)
    
    return primers 
