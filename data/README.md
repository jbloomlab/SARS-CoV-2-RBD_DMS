# Input data
This directory contains input data for the analysis.

## Basic information about sequences and alignments

These files are used for the basic processing of the deep sequencing data to call variants by barcode and count barcodes:

   - [data/PacBio_amplicons.gb](data/PacBio_amplicons.gb): the amplicons being sequenced by PacBio.
     Note that there are a variety of possible amplicons because in addition to the SARS-CoV-2 RBD there are RBDs from a variety of other viral strains.

   - [data/feature_parse_specs.yaml](data/feature_parse_specs.yaml): how to parse the amplicon when handling the PacBio data.

   - [data/PacBio_runs.csv](data/PacBio_runs.csv): list of the PacBio runs used to call the variants.

   - [data/barcode_runs.csv](data/barcode_runs.csv): list of the Illumina runs used to count the barcodes for different samples.

   - [RBD_sites.csv](RBD_sites.csv): gives site and residue information for SARS-CoV-2, including alignment of the RBD integer numbering with the Spike numbering for SARS-CoV-2 RBD, alignment to SARS-CoV, and structural annotations as detailed below.

## Alignments of different Spikes / RBDs
The [./alignments/](alignments) subdirectory contains alignments of Spike and RBD constructs from human, bat, and other mammal SARSr-CoV isolates.

## Isogenic titrations
The [./isogenic_titrations/](isogenic_titrations) subdirectory contains processed data from isogenic titration experiments use to fit single-variant titration curves.

## Structures
The [./structures/](structures) directory contains information relevant for structural analyses.
Specifically, [./structures/Ab-bound/](structures/Ab-bound) and [./structures/ACE2-bound](./structures/ACE2-bound) contain PDB files of SARS-CoV and SARS-CoV-2 spike and RBD proteins, either free, bound to ACE2, or antibodies. Each subdirectory contains the `.pdb` files, and a `.pse` session illustrating structural alignment of the component `.pdb` files.

The `./structures/annotate-structural-elements.R` script analyzes the various PDB files to create strucutral annotations, including relative solvent accessibility, ACE2 contact residues, and antibody epitope sites. This script was run to create the file `./RBD_sites.csv`


`PacBio_amplicon.gb` encodes the PacBio amplicon to map for library variants of the SARS-CoV-2 RBD

