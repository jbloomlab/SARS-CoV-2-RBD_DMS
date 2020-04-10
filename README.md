# Deep mutational scanning of SARS-CoV-2 RBD
Analysis of deep mutational scanning of barcoded codon variants of SARS-CoV-2 RBD

Study and analysis by Tyler Starr, Allie Greaney, and [Jesse Bloom](https://research.fhcrc.org/bloom/en.html).

## Summary of workflow and results
For a summary of the workflow and results, [click here](results/summary/summary.md).
Reading this summary is the best way to understand the analysis.

## Running the analysis
The analysis consists of three components, all of which are contained in this repository:

 1. Instructions to build the computing environment.

 2. The computer code itself.

 3. The required input data.

First, set up the computing environment, which is done via `conda`.
Ensure you have `conda` installed; if not install it via Miniconda as described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation).
The environment is specified in [environment.yml](environment.yml).
If you have not previously built the conda environment, then build the environment specified in [environment.yml](environment.yml) to `./env` with:

    conda env create -f environment.yml -p ./env

Then activate it with:

    conda activate ./env

If you've previously built the environment into `./env`, just do the activation step.

Now you can run the entire analysis.
The analysis consists primarily of a series of Jupyter notebooks and R scripts in to the top-level directory along with some additional code in [Snakefile](Snakefile).
You can run the analysis by using [Snakemake](https://snakemake.readthedocs.io) to run [Snakefile](Snakefile), specifying the conda environment in `./env`, as in:

    snakemake --use-conda --conda-prefix ./env

However, you probably want to use the server to help with computationally intensive parts of the analysis.
To run using the cluster configuration for the Fred Hutch server, simply run the bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash), which executes [Snakefile](Snakefile) in a way that takes advantage of the Hutch server resources.
This bash script also automates the environment building steps above, so really all you have to do is run this script.
You likely want to submit [run_Hutch_cluster.bash](run_Hutch_cluster.bash) itself to the cluster (since it takes a while to run) with:

    sbatch -t 7-0 run_Hutch_cluster.bash

## Configuring the analysis
The configuration for the analysis is specifed in [config.yaml](config.yaml).
This file defines key variables for the analysis, and should be relatively self-explanatory.
You should modify the analysis by changing this configuration file; do **not** hard-code crucial experiment-specific variables within the Jupyter notebooks or `Snakefile`.

The input files pointed to by [config.yaml](config.yaml) are in the [./data/](data) subdirectory.
See the [./data/README.md](./data/README.md) file for details.

  - [data/amplicon.gb](data/amplicon.gb): the amplicon being sequenced by PacBio
  - [data/feature_parse_specs.yaml](data/feature_parse_specs.yaml): how to parse the amplicon when handling the PacBio data
  - [data/PacBio_runs.csv](data/PacBio_runs.csv): list of the PacBio runs used to call the variants
  - [data/barcode_runs.csv](data/barcode_runs.csv): list of the Illumina runs used to count the barcodes for different samples
  - [data/experimental_schematic.jpg](data/experimental_schematic.jpg): schematic of experiment

## Cluster configuration
There is a cluster configuration file [cluster.yaml](cluster.yaml) that configures [Snakefile](Snakefile) for the Fred Hutch cluster, as recommended by the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).
The [run_Hutch_cluster.bash](run_Hutch_cluster.bash) script uses this configuration to run [Snakefile](Snakefile).
If you are using a different cluster than the Fred Hutch one, you may need to modify the cluster configuration file.

## Notebooks that perform the analysis
The Jupyter notebooks that perform most of the analysis are in this top-level directory with the extension `*.ipynb`.
These notebooks read the key configuration values from [config.yaml](config.yaml).

There is also a [./scripts/](scripts) subdirectory with related scripts.

The notebooks need to be run in the order described in [the workflow and results summary](results/summary/summary.md).
This will occur automatically if you run them via [Snakefile](Snakefile) as described above.

## Results
Results are placed in the [./results/](results) subdirectory.
Most results are not tracked as they are very large.
Ones that are tracked:

  - [./results/summary](./results/summary): contains the Markdown output of the Jupyter notebooks as well as the top-level summary in [./results/summary/summary.md](./results/summary/summary.md).

## Updating the conda environment
The [environment.yml](environment.yml) file contains a fully pinned conda environment.
An environment without all of the versions pinned is in [environment_unpinned.yml](environment_unpinned.yml).
If you need to update the environment, the suggested way to do it is add the new requirement to [environment_unpinned.yml](environment_unpinned.yml), then build that environment and finally export the pinned version:

    conda env create -f environment_unpinned.yml --prefix ./env

Then activate the environment with:

    conda activate ./env

Finally, export the pinned version with:

    conda env export --prefix ./env > environment.yml
