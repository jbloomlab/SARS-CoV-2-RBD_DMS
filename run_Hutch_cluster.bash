#!/bin/bash

snakemake \
    -j 999 \
    --cluster-config cluster.yaml \
    --cluster "sbatch -c {cluster.cpus} --mem={cluster.mem} -t {cluster.time} -J {cluster.name}" \
    --latency-wait 30 \
    --use-envmodules \
    --use-conda \
    --conda-prefix ./env
