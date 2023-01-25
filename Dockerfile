# syntax=docker/dockerfile:1
FROM python:3.8-slim-bullseye as base

WORKDIR /

ENV PYTHONUNBUFFERED: 1

RUN uname -r

RUN apt-get update 
RUN apt-get install -y build-essential
RUN apt-get install -y libz-dev
RUN apt-get install -y bedtools

COPY requirements.txt requirements.txt
COPY sge-primer-scoring/requirements.txt scoring_requirements.txt

RUN pip3 install -U pip wheel setuptools 
RUN pip3 install -r requirements.txt
RUN pip3 install -r scoring_requirements.txt

COPY . .

FROM base as dev
ENTRYPOINT [ "python3", "./src/cli.py" ]

FROM base as test
CMD ["python3", "-m", "unittest"]


# TODO: Flesh out loading data into Docker Container ticket
