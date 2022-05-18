# syntax=docker/dockerfile:1
FROM python:3.10-slim-buster

WORKDIR ./

ENV PYTHONUNBUFFERED: 1
ENV FLASK_APP=./targeton_designer.py
ENV FLASK_RUN_HOST=0.0.0.0

RUN apt-get update && apt-get install -y build-essential
RUN apt-get install -y libz-dev bedtools
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

EXPOSE 5000

COPY . .

CMD ["flask", "run"]

