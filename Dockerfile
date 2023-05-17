FROM python:3.8.0 as base

WORKDIR /

FROM base as build
COPY Makefile Makefile
COPY requirements.txt requirements.txt
COPY sge-primer-scoring sge-primer-scoring
COPY src src
COPY tests tests
RUN apt-get update
RUN make install
RUN make setup-venv

FROM base as unittest
ENV DOCKER_ENV=${DOCKER_ENV:-unittest}
RUN ls -l
CMD [ "sh", "-c", "make test" ]


# TODO: Flesh out loading data into Docker Container ticket
