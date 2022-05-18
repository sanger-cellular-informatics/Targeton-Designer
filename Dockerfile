# syntax=docker/dockerfile:1
FROM python:3.10-alpine
WORKDIR ./
ENV PYTHONUNBUFFERED: 1
ENV FLASK_APP=./targeton_designer.py
ENV FLASK_RUN_HOST=0.0.0.0
RUN apk add --no-cache gcc musl-dev linux-headers make
RUN apk update && apk add zlib-dev 
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
RUN mv bedtools.static.binary bedtools
RUN chmod a+x bedtools
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt
EXPOSE 5000
COPY . .
CMD ["flask", "run"]

