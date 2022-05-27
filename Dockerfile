# syntax=docker/dockerfile:1
FROM python:3.10-slim-buster AS slicer

WORKDIR /

ENV PYTHONUNBUFFERED: 1

RUN apt-get update 
RUN apt-get install -y build-essential
RUN apt-get install -y libz-dev 
RUN apt-get install -y bedtools

COPY designer/requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY designer .

ENTRYPOINT [ "python3", "./slicer.py" ]

FROM python:3.10-slim-buster AS primer3

WORKDIR /

ENV PYTHONUNBUFFERED: 1

RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y libz-dev 

COPY runner/requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY runner .
COPY cmd.py .

ENTRYPOINT [ "python3", "./cmd.py" ]

