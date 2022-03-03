# Deep mutational scanning of SARS-CoV-2 RBD
Analysis of deep mutational scanning of barcoded codon variants of SARS-CoV-2 RBD

Study and analysis by Tyler Starr, Allie Greaney, [Jesse Bloom](https://research.fhcrc.org/bloom/en.html), and co-authors.

For the full paper, see [Deep mutational scanning of SARS-CoV-2 receptor binding domain reveals constraints on folding and ACE2 binding](https://doi.org/10.1016/j.cell.2020.08.012).

**NOTE added 3 March, 2022. We have an improved deep mutational scanning dataset that we recommend using in place of this 2020 dataset. Please find the updated data [here](https://jbloomlab.github.io/SARS-CoV-2-RBD_DMS_variants/RBD-heatmaps/), the preprint describing the new datasets [here](https://www.biorxiv.org/content/10.1101/2022.02.24.481899v1), and the new GitHub repository [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_variants).**

## Summary of workflow and results
For a summary of the workflow and links to key results files, [click here](results/summary/summary.md).
Reading this summary is the best way to understand the analysis.

## Running the analysis
The analysis consists of three components, all of which are contained in this repository:

 1. Instructions to build the computing environment.

 2. The computer code itself.

 3. The required input data.

First, set up the computing environment, which is partially done via `conda`.
Ensure you have `conda` installed; if not install it via Miniconda as described [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/#regular-installation).
The environment is specified in [environment.yml](environment.yml).
If you have not previously built the conda environment, then build the environment specified in [environment.yml](environment.yml) to `./env` with:

    conda env create -f environment.yml -p ./env

Then activate it with:

    conda activate ./env

If you've previously built the environment into `./env`, just do the activation step.

Setting up the `conda` environment above installs everything to run all parts of the analysis **except** the `R` markdown notebooks.
For those, the pipeline currently uses the Fred Hutch computing cluster module `R/3.6.1-foss-2018b` as specified in `Snakefile`.
That module is not packaged with this repo, so if you aren't on the Fred Hutch cluster you'll have to create a similar `R` environment yourself (all the `R` packages are listed at the beginning of their output in the [summary results](results/summary/summary.md).

Now you can run the entire analysis.
The analysis consists primarily of a series of Jupyter notebooks and R markdown in to the top-level directory along with some additional code in [Snakefile](Snakefile).
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

Note that the raw sequencing data are on the SRA in BioProject [PRJNA639956](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA639956) as well as on the Hutch cluster.

## Cluster configuration
There is a cluster configuration file [cluster.yaml](cluster.yaml) that configures [Snakefile](Snakefile) for the Fred Hutch cluster, as recommended by the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).
The [run_Hutch_cluster.bash](run_Hutch_cluster.bash) script uses this configuration to run [Snakefile](Snakefile).
If you are using a different cluster than the Fred Hutch one, you may need to modify the cluster configuration file.

## Notebooks that perform the analysis
The Jupyter notebooks and R markdown dscripts that perform most of the analysis are in this top-level directory with the extension `*.ipynb` or `*.Rmd`.
These notebooks read the key configuration values from [config.yaml](config.yaml).

There is also a [./scripts/](scripts) subdirectory with related scripts.

The notebooks need to be run in the order described in [the workflow and results summary](results/summary/summary.md).
This will occur automatically if you run them via [Snakefile](Snakefile) as described above.

## Results
Results are placed in the [./results/](results) subdirectory.
Many of the files created in this subdirectory are not tracked in the `git` repo as they are very large.
However, key results files are tracked as well as a summary that shows the code and results.
Click [here](./results/summary/summary.md) to see that summary.

The large results files are tracked via [git-lfs](https://git-lfs.github.com/).
This requires `git-lfs` to be installed, which it is in the `conda` environment specified by [environment.yml](environment.yml).
The following commands were then run:

    git lfs install

You may need to run this if you are tracking these files and haven't installed `git-lfs` in your user account.
Then the large results files were added for tracking with:

    git lfs track results/variants/codon_variant_table.csv
    git lfs track results/counts/variant_counts.csv
    git lfs track results/expression_meanFs/expression_meanFs.csv
    git lfs track results/expression_meanFs/expression_meanFs_homologs.csv
    git lfs track results/binding_Kds/binding_Kds.csv
    git lfs track results/binding_Kds/binding_Kds_homologs.csv
    git lfs track results/single_mut_effects/single_mut_effects.csv
    git lfs track results/single_mut_effects/homolog_effects.csv
    git lfs track results/structure_function/6m0j_b-factor-mean-bind.pdb
    git lfs track results/structure_function/6m0j_b-factor-mean-expr.pdb

The one exception is the `interactive_heatmap.html` which is placed in `./docs/_includes/` so that it can be rendered using GitHub Pages.

## Updating the conda environment
The [environment.yml](environment.yml) file contains a fully pinned conda environment.
An environment without all of the versions pinned is in [environment_unpinned.yml](environment_unpinned.yml).
If you need to update the environment, the suggested way to do it is add the new requirement to [environment_unpinned.yml](environment_unpinned.yml), then build that environment and finally export the pinned version:

    conda env create -f environment_unpinned.yml --prefix ./env

Then activate the environment with:

    conda activate ./env

Finally, export the pinned version with:

    conda env export --prefix ./env > environment.yml
