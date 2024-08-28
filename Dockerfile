FROM python:3.8.0

ARG USER_ID=1000

RUN apt-get update && \
    apt-get -y install build-essential && \
    apt-get -y install bedtools && \
    apt-get -y install curl && \
    apt-get -y install python3-setuptools

USER root

RUN groupadd -g $USER_ID primerdesigner && \
    useradd -m -u $USER_ID -g primerdesigner primerdesigner && \
    usermod -u $USER_ID primerdesigner

WORKDIR /targeton-designer

COPY --chown=primerdesigner:primerdesigner . /targeton-designer

RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir -r sge-primer-scoring/requirements.txt

USER primerdesigner