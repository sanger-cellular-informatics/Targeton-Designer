FROM python:3.8.0 as base

WORKDIR /

COPY Makefile Makefile
COPY requirements.txt requirements.txt
COPY sge-primer-scoring/requirements.txt sge-primer-scoring/scoring_requirements.txt
COPY src src
COPY tests tests

RUN apt-get update

RUN make install
RUN make setup-venv

COPY . .

FROM base as unittest
ENV DOCKER_ENV=${DOCKER_ENV:-unittest}
CMD [ "sh", "-c", "make test" ]


# TODO: Flesh out loading data into Docker Container ticket
