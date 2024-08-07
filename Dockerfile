FROM python:3.8.0 as base

WORKDIR /

COPY . /

RUN apt-get update

RUN make install
RUN make setup-venv
