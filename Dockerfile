# syntax=docker/dockerfile:1
FROM python:3.10-alpine
WORKDIR ./
ENV PYTHONUNBUFFERED: 1
ENV FLASK_APP=./targeton_designer.py
ENV FLASK_RUN_HOST=0.0.0.0
RUN apk add --no-cache gcc musl-dev linux-headers make
RUN apk update && apk add zlib-dev 
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
RUN tar -zxvf bedtools-2.29.1.tar.gz
RUN cd bedtools2 && make
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt
EXPOSE 5000
COPY . .
CMD ["flask", "run"]

