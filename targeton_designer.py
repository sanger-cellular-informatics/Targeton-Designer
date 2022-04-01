import redis
from flask import Flask

app = Flask(__name__)
cache = redis.Redis(host='redis', port=6379)

from runner.spawner import primer3_runner

@app.route('/')
def hello():
    return primer3_runner()
