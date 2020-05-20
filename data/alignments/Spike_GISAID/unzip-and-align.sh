#!/bin/bash

unzip spikeprot*.zip

sed 's, ,_,g' -i spikeprot*.fasta

mafft --thread 4 --op 10 spikeprot*.fasta > spike_GISAID_aligned.fasta