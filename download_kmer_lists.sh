#!/bin/sh

mkdir kmer

curl -o kmer/homo_sapiens_16.list https://primer3.ut.ee/lists/homo_sapiens_16.list
curl -o kmer/homo_sapiens_11.list https://primer3.ut.ee/lists/homo_sapiens_11.list

ls kmer
