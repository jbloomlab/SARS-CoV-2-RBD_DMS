# Input data
This directory contains input data for the analysis.

## Basic information about sequences and alignments

These files are used for the basic processing of the deep sequencing data to call variants by barcode and count barcodes:

   - [PacBio_amplicons.gb](PacBio_amplicons.gb): the amplicons being sequenced by PacBio.
     Note that there are a variety of possible amplicons because in addition to the SARS-CoV-2 RBD there are RBDs from a variety of other viral strains.

   - [feature_parse_specs.yaml](feature_parse_specs.yaml): how to parse the amplicon when handling the PacBio data.

   - [PacBio_runs.csv](PacBio_runs.csv): list of the PacBio runs used to call the variants.

   - [barcode_runs.csv](barcode_runs.csv): list of the Illumina runs used to count the barcodes for different samples. This file is built from [barcode_runs_orig-names.csv](barcode_runs_orig-names.csv) by the Jupyter notebook [build_barcode_runs.ipynb](build_barcode_runs.ipynb).

   - [RBD_sites.csv](RBD_sites.csv): gives site and residue information for SARS-CoV-2, including alignment of the RBD integer numbering with the Spike numbering for SARS-CoV-2 RBD, alignment to SARS-CoV, and structural annotations as detailed below.

## Alignments of different Spikes / RBDs
The [./alignments/](alignments) subdirectory contains alignments and phylogenetic trees of Spike and RBD constructs from human, bat, and other mammal sarbecovirus isolates. Sequences were aligned both as nucleotide gene sequecnes and as translated amino acid sequences. The script [./alignments/alignment-and-phylogeny.sh](alignments/alignment-and-phylogeny.sh) was used to align sequences with `mafft` and infer phylogenetic trees with `raxml`. The unaligned sequence inputs to this script are found in [./alignments/unaligned-sequences/](alignments/unaligned-sequences/), with each fasta header providing common names for sarbecovirus isolates and their genome accession number in NCBI or GISAID. Sequence selection includes:
   - The curated set of 30 sequences from [Letko, Marzi, and Munster (2020)](https://www.nature.com/articles/s41564-020-0688-y), which contained all known unique RBD sequences at the time of its publication. See [Extended Data Fig 1](https://www.nature.com/articles/s41564-020-0688-y/figures/6) for more information.
   - Bat sarbecovirus isolate WIV16, which was removed from the Letko dataset because it has an identical RBD amino acid sequence to WIV1.
   - RaTG13 and RmYN02, two newly described bat CoV sequences that are the most closely related strains currently known to SARS-CoV-2 (at least at the whole-genome level).
   - Two recent pangolin CoV sequences from [Lam _et al._ (2020)](https://www.nature.com/articles/s41586-020-2169-0), including the infectious virus isolated from a pangolin seized in Guanxi, and the consensus sequence from two closely related isolates from pangolins seized in Guangdong. The consensus sequence was constructed from an alignment of the GD pangolin CoV in GISAID accssion EPI\_ISL\_410544, and the reconstructed metagenomic viral genome reported in the [Lam _et al._ (2020)](https://www.nature.com/articles/s41586-020-2169-0) supplement.
   - BtKY72, a bat coronavirus isolated in Kenya which, along with the Bulgaria sample BM48-31 is used to root the phylogeny

## Input data from previous studies
The [./lit-measurements](lit-measurements) subdirectory contains raw data on sarbecovirus functional properties from previous experimental studies, which we will in comparison to our deep mutational scanning measurements. Files include:
   - [./lit-measurements/Letko\_2020\_data\_1e.csv](lit-measurements/Letko_2020_data_1e.csv), which reports quantitative measurements of the ability of VSV pseudoviruses containing SARS-CoV-1 Spike as a chimera with various sarbecovirus RBDs to mediate entry into BHK (baby hamster kidney) cells expressing human ACE2. This data is represented in [Figure 1e](https://www.nature.com/articles/s41564-020-0688-y/figures/1) of [Letko _et al.__ 2020](https://www.nature.com/articles/s41564-020-0688-y).
   - [./lit-measurements/Shang\_2020\_data\_3b.csv](lit-measurements/Shang_2020_data_3b.csv), which reports quantitative measurement of the ability of retrovirus pseudoviruses bearing Spike from SARS-CoV-2 or RaTG13 to enter HEK293T cells expressing human ACE2. This data is represented in Figure 3b of [Shang _et al._ 2020](https://www.nature.com/articles/s41586-020-2179-y).

## Isogenic titrations
The [./isogenic_titrations/](isogenic_titrations) subdirectory contains processed data from isogenic titration experiments use to fit single-variant titration curves.

## Structures
The [./structures/](structures) directory contains information relevant for structural analyses.
Specifically, [./structures/Ab-bound/](structures/Ab-bound) and [./structures/ACE2-bound](./structures/ACE2-bound) contain PDB files of SARS-CoV and SARS-CoV-2 spike and RBD proteins, either free, bound to ACE2, or antibodies. Each subdirectory contains the `.pdb` files, and a `.pse` session illustrating structural alignment of the component `.pdb` files.

The [./structures/annotate-structural-elements.R](structures/annotate-structural-elements.R) script analyzes the various PDB files to create strucutral annotations, including relative solvent accessibility, ACE2 contact residues, and antibody epitope sites. This script was run to create the file [./RBD_sites.csv](RBD_sites.csv)
