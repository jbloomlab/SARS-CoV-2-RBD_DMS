# Summary of data directory components

## `./alignments` directory
Alignments of Spike and RBD constructs from human, bat, and other mammal SARSr-CoV isolates

## `./isogenic_titrations` directory
Processed data from isogenic titration experiments use to fit single-variant titration curves

## `./Pacbio_amplicon_controls` directory
Genbank files of the 11 internal control plasmids from bat, pangolin, and SARS-CoV RBDs that are spiked into the SARS-CoV_2 RBD variant library. These constructs are the expected PacBio sequencing constructs as submitted.

## `./structures` directory
the `./structures/Ab-bound` and `./structures/ACE2-bound` sub-directories contain PDB files of SARS-CoV and SARS-CoV-2 spike and RBD proteins, either free, bound to ACE2, or antibodies. Each subdirectory contains the `.pdb` files, and a `.pse` session illustrating structural alignment of the component `.pdb` files.

The `./structures/annotate-structural-elements.R` script analyzes the various PDB files to create strucutral annotations, including relative solvent accessibility, ACE2 contact residues, and antibody epitope sites. This script was run to create the file `./RBD_sites.csv`

## Other files

`barcode_runs.csv` and `PacBio_runs.csv` list Illumina and PacBio sequenecing runs, respectively, and required info for the analyses launced in the Snakemake pipeline

`feature_parse_specs.yaml` gives parsing specifications for `alignparse`

`PacBio_amplicon.gb` encodes the PacBio amplicon to map for library variants of the SARS-CoV-2 RBD

`RBD_sites.csv` gives site and residue information for SARS-CoV-2, including alignment of the RBD integer numbering with the Spike numbering for SARS-CoV-2 RBD, alignment to SARS-CoV, and structural annotations as described above
