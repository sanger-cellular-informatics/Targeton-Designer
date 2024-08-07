#!/bin/sh

docker run -v $(pwd)/kmer/:/kmer --rm -it --entrypoint bash $1