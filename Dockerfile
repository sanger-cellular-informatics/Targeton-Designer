FROM python:3.8.0 as base

WORKDIR /

FROM base as build
COPY Makefile Makefile
COPY requirements.txt requirements.txt
COPY sge-primer-scoring sge-primer-scoring
RUN ls sge-primer-scoring
COPY src src
COPY tests tests
RUN apt-get update
RUN make install
RUN make setup-venv
