import redis
from flask import Flask

app = Flask(__name__)
cache = redis.Redis(host='redis', port=6379)

from runner.primer3 import primer3_runner
from runner.exonerate import run_exonerate

@app.route('/')
def start():
    primers = primer3_runner()
    #run_exonerate()
    return primers   
