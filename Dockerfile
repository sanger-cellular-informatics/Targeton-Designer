# syntax=docker/dockerfile:1
FROM python:3.10-alpine
WORKDIR /runner
ENV FLASK_APP=targeton_designer.py
ENV FLASK_RUN_HOST=0.0.0.0
RUN apk add --no-cache gcc musl-dev linux-headers make
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt
EXPOSE 5000
COPY . .
CMD ["flask", "run"]

