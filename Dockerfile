FROM python:3.8.0 as base

RUN apt-get update && \
    apt-get -y install build-essential && \
    apt-get -y install bedtools && \
    apt-get -y install curl && \
    apt-get -y install python3-setuptools

RUN groupadd -g 1001 primerdesigner && \
    useradd -m -u 1001 -g primerdesigner primerdesigner

WORKDIR /targeton-designer

COPY . /targeton-designer

RUN mkdir /targeton-designer/td_output
RUN chown -hR primerdesigner:primerdesigner /targeton-designer/logs /targeton-designer/td_output && \
    chmod 770 /targeton-designer/td_output

RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir -r sge-primer-scoring/requirements.txt

USER primerdesigner