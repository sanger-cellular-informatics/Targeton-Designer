FROM python:3.8.0

ARG USER_ID=1001
ARG USER_NAME=primerdesigner

RUN apt-get update && \
    apt-get -y install build-essential && \
    apt-get -y install bedtools && \
    apt-get -y install curl && \
    apt-get -y install python3-setuptools

USER root

RUN groupadd -g $USER_ID $USER_NAME && \
    useradd -m -u $USER_ID -g $USER_NAME $USER_NAME && \
    usermod -u $USER_ID $USER_NAME

WORKDIR /targeton-designer

COPY --chown=$USER_NAME:$USER_NAME . /targeton-designer

RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir -r sge-primer-scoring/requirements.txt

USER $USER_NAME