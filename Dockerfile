FROM python:3.8.0 as base

WORKDIR /

COPY Makefile Makefile
COPY requirements.txt requirements.txt
COPY sge-primer-scoring/requirements.txt sge-primer-scoring/requirements.txt
COPY sge-primer-scoring/src sge-primer-scoring/src
COPY src src
COPY tests tests
COPY download_kmer_lists.sh download_kmer_lists.sh

RUN apt-get update

RUN make install
RUN make setup-venv

FROM base as unittest
ENV DOCKER_ENV=${DOCKER_ENV:-unittest}
CMD [ "sh", "-c", "make test" ]
