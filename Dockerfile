FROM python:3.8.0

WORKDIR /

COPY . /

RUN apt-get update && \
    apt-get -y install build-essential && \
    apt-get -y install bedtools && \
    apt-get -y install curl && \
    apt-get -y install python3-setuptools

RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir -r sge-primer-scoring/requirements.txt