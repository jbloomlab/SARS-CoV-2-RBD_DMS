# Deep mutational scanning of SARS-CoV-2 RBD
Analysis of deep mutational scanning of barcoded codon variants of SARS-CoV-2 RBD

Study and analysis by Tyler Starr, Allie Greaney, and [Jesse Bloom](https://research.fhcrc.org/bloom/en.html). Code modified from NIH45-46_DMS analysis by Tyler Starr and Jesse Bloom.

## Summary of workflow and results
For a summary of the workflow and results, [click here](results/summary/summary.md).
Reading this summary is the best way to understand the analysis.

## Running the analysis
The analysis consists of three components, all of which are contained in this repository:

 1. Instructions to build the computing environment to run the code.

 2. The computer code itself.

 3. The required input data for the computer code.

First, set up the computing environment, which is done via `conda`.
Ensure you have `conda` installed; if not install it via Miniconda as described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation).
Next, create the environment specified in [environment.yml](environment.yml).
To do this, begin by deactivating any existing conda environments:

    conda deactivate

The next step differs depending on whether you have built the conda environment previously or not (i.e., do you have already have a subdirectory called `./env`?).
If you have not previously built the conda environment, then build the environment specified in [environment.yml](environment.yml) to `./env` with:

    conda env create -f environment.yml -p ./env

Then activate it with:

    conda activate ./env

If you've previously built the environment into `./env`, just do the activation step.

Now you can run the entire analysis.
The analysis consists primarily of a series of Jupyter notebooks in to the top-level directory along with some additional code in [Snakefile](Snakefile).
You can run the analysis by using [Snakemake](https://snakemake.readthedocs.io) to run [Snakefile](Snakefile), specifying the conda environment in `./env`, as in:

    snakemake --use-conda --conda-prefix ./env

However, you probably want to use the server to help with computationally intensive parts of the analysis.
To run using the cluster configuration for the Fred Hutch server, simply run the bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash), which executes [Snakefile](Snakefile) in a way that takes advantage of the Hutch server resources.
This bash script also automates the environment building steps above, so really all you have to do is run this script.
You likely want to submit [run_Hutch_cluster.bash](run_Hutch_cluster.bash) itself to the cluster (since it takes a while to run) with:

    sbatch -t 7-0 run_Hutch_cluster.bash

Note that whenever you are done working on the project, you should deactivate the conda environment with:

    conda deactivate

## Configuring the analysis
The configuration for the analysis is specifed in [config.yaml](config.yaml), as recommended by the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).
This file defines key variables for the analysis, and should be relatively self-explanatory.
You should modify the analysis by changing this configuration file; do **not** hard-code crucial experiment-specific variables within the Jupyter notebooks or `Snakefile`.

Key input files pointed to by [config.yaml](config.yaml) are in the [./data/](data) subdirectory:

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
These notebooks read the key configuration values from [config.yaml](config.yaml); this avoids having to hardcode values throughout the notebooks.

There is also a [./scripts/](scripts) subdirectory with related scripts.

The notebooks need to be run in the order described in [the workflow and results summary](results/summary/summary.md).
This will occur automatically if you run them all via [Snakefile](Snakefile) as described above.
You have to run them by [Snakefile](Snakefile) to get the [nice Markdown summary](results/summary/summary.md).

## Results
Results are placed in the [./results/](results) subdirectory.
Most results are not tracked as they are very large.
Ones that are tracked:

  - [./results/summary](./results/summary): contains the Markdown output of the Jupyter notebooks as well as the top-level summary in [./results/summary/summary.md](./results/summary/summary.md).

## Creating / updating the conda environment
If you need to create a new conda environment, follow these steps to build an environment to `./env` and then export it to an `environment.yml` file:
 
1. Ensure you have `conda` installed. If not, you can install it via Miniconda as described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation).

2. Deactivate any existing `conda` environment with:

        conda deactivate

3. Create the environment in the desired location. To create an environment in the location `./env` for Python 3.7 and the `pip` Python installation tool, use:
 
        conda create --prefix ./env python=3.7 pip
        
4. Activate the `conda` environment with:

        conda activate ./env

5. Install any additional `conda` packages we need. Here we install `pbccs` (the PacBio `ccs` program), the `sra-toolkit` (which includes `fastq-dump`), and `minimap2` with:

        conda install pbccs sra-tools minimap2

6. Use `pip` to install Python dependencies from `PyPI`:

        pip install jupyter jupyter_contrib_nbextensions nbdime numpy pandas plotnine pyyaml snakemake

7. Use `pip` to install Python dependencies only available on GitHub as [described here](https://www.reddit.com/r/Python/comments/2crput/how_to_install_with_pip_directly_from_github/). If these are in private GitHub repos, `git` must be configured to access your credentials. This creates a `./src` directory since they are installed with `-e`:

        pip install -e git+https://github.com/jbloomlab/dms_variants.git#egg=dms_variants
        pip install -e git+https://github.com/jbloomlab/alignparse.git#egg=alignparse

8. Export the environment to an `environment.yml` file with:

        conda env export --prefix ./env > environment.yml

9. If you have installed any Python dependencies from GitHub as in step 7 above, they will **not** be in the created `environment.yml` file. You can see how they should be listed with `pip freeze` and then manually add the lines for them to `environment.yml`.
