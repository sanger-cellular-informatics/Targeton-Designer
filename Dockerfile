# syntax=docker/dockerfile:1
FROM python:3.8-slim-buster as base

WORKDIR /

ENV PYTHONUNBUFFERED: 1

RUN apt-get update 
RUN apt-get install -y build-essential
RUN apt-get install -y libz-dev
RUN apt-get install -y bedtools

COPY requirements.txt requirements.txt

RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt

COPY . .

FROM base as dev
ENTRYPOINT [ "python3", "./src/cli.py" ]

FROM base as test
CMD ["python3", "-m", "unittest"]

