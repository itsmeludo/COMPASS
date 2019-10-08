#!/bin/bash

cat ../data/COMPASS.fasta.xz.part* > ../data/COMPASS.fasta.xz

xz --decompress ../data/COMPASS.fasta.xz
