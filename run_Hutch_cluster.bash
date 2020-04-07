#!/bin/bash -i
# bash shebang has -i for reason here: https://stackoverflow.com/a/55507956

set -e

echo "Using following version of conda:"
conda --version

conda deactivate

if [ ! -d "./env" ]
then
    echo "Building conda environment..."
    conda env create -f environment.yml -p ./env
fi

echo "Activating conda environment..."
conda activate ./env

echo "Running snakemake..."
snakemake \
    --use-conda \
    --conda-prefix ./env \
    -j 999 \
    --cluster-config cluster.yaml \
    --cluster "sbatch -p {cluster.partition} -c {cluster.cpus} --mem={cluster.mem} -t {cluster.time} -J {cluster.name}" \
    --latency-wait 30

echo "De-activating conda environment..."
conda deactivate

echo "Script complete."
