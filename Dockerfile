FROM python:3.8.0 as base

WORKDIR /
COPY Makefile Makefile
RUN apt-get update
RUN make install

COPY requirements.txt requirements.txt
COPY sge-primer-scoring/requirements.txt sge-primer-scoring/requirements.txt
RUN make setup-venv

COPY src src
COPY tests tests
COPY sge-primer-scoring sge-primer-scoring
