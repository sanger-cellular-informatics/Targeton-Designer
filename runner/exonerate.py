import docker

def run_exonerate():
    client = docker.from_env()
    client.containers.run('exonerate','--version')
