import redis
from flask import Flask

app = Flask(__name__)
cache = redis.Redis(host='redis', port=6379)

from runner.primer3 import primer3_runner

@app.route('/')
def start():
    primers = primer3_runner()
    return primers   
