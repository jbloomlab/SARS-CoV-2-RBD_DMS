# Summary

Analysis run by [Snakefile](../../Snakefile)
using [this config file](../../config.yaml).
See the [README in the top directory](../../README.md)
for details.

Here is the DAG of the computational workflow:
![dag.svg](dag.svg)

Here is the Markdown output of each Jupyter notebook in the
workflow:

1. [Process PacBio CCSs](process_ccs.md).

2. [Build variants from CCSs](build_variants.md).
   Creates a [codon variant table](../variants/codon_variant_table.csv)
   linking barcodes to the mutations in the variants.

3. [Count variants by barcode](count_variants.md).
   Creates a [variant counts file](../counts/variant_counts.csv)
   giving counts of each barcoded variant in each condition.

4. [QC analysis of sequencing counts](analyze_counts.md).

5. [Computation of ACE2-binding *K*<sub>D</sub>](compute_binding_Kd.md).
   Creates files giving the ACE2-binding of each barcoded variant
   [of SARS-CoV-2 RBD](../binding_Kds/binding_Kds.csv) and of
   [the homologs](../binding_Kds/binding_Kds_homologs.csv).

6. [Computation of expression mean fluorescence](compute_expression_meanF.md).
   Creates files giving the expression of each barcoded variant
   [of SARS-CoV-2 RBD](../expression_meanFs/expression_meanFs.csv) and of
   [the homologs](../expression_meanFs/expression_meanFs_homologs.csv).

7. [Global epistasis decomposition of binding effects](global_epistasis_binding.md).

8. [Global epistasis decomposition of expression effects](global_epistasis_expression.md).

9. [Calculation of final single mutant effects on binding and expression](single_mut_effects.md).
   Creates files giving the estimated expression and ACE2-binding of
   [single mutants to SARS-CoV-2 RBD](../single_mut_effects/single_mut_effects.csv)
   and [the homologs](../single_mut_effects/homolog_effects.csv).

10. [Structure-function analysis of mutational effects](structure_function.md).

11. [Characterization of amino acid variants in sequenced SARS-CoV-2 isolates](circulating_variants.md).

12. [Analysis of mutational tolerance within RBD antibody epitopes](antibody_epitopes.md)

13. [Analysis of amino acid diversity in natural sarbecovirus isolates](sarbecovirus_diversity.md)