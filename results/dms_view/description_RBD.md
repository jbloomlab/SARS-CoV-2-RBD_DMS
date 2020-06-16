## Deep mutational scanning of SARS-CoV-2 receptor binding domain

Tyler N. Starr, Allison J. Greaney, Sarah K. Hilton, Katharine H.D. Crawford,  ..., David Veesler, Jesse D. Bloom

For background, see our paper [here]().

## What data are shown here?

We are showing the effects of all mutations to the SARS-CoV-2 RBD as measured via deep mutational scanning. Menu bars can alternate the view between our measurements of ACE2 receptor binding affinity and RBD expression (a correlate of protein stability). Within each phenotype, different metrics can be visualized at the per-site and per-mutation level:

   - per-site `entropy` and `n_effective` describe the degree of constraint reflected in the amino acid preferences illustrated in logoplots as in Figure S4 of the manuscript. Lower `entropy` or `n_effective` indicates stronger mutational constraint.
   - per-site `mean-effect`, `min-effect`, and `max-effect` are statistics on the raw measurement values of mutational effects on binding and expression, as reflected in the heatmaps in Figure 3. Lower values of these metrics indicate stronger mutational constraint.
   - per-mutation `preference` is a relative preference value useful for visualization in logoplots, as the sum of preferences for all 20 amino acids at a site sum to 1. Larger values (taller letters) indicate more preferred amino acids.
   - per-mutation `delta_effect` is the raw measurement value of mutational effect on binding and expression, either delta-log<sub>10</sub>(_K_<sub>D,app</sub>) for binding, or delta-log(MFI) for expression. The height of the letter indicates the magnitude of mutational effect, with letters below the line indicating deleterious effects.

When sites are selected, their location is indicated on the ACE2-bound RBD crystal structure (PDB [6M0J](https://www.rcsb.org/structure/6M0J), from [Lan _et al._ 2020](https://www.nature.com/articles/s41586-020-2180-5))

For `dms-view` renderings of specific residue sets discussed in the paper, see the links on [this page](https://jbloomlab.github.io/SARS-CoV-2-RBD-DMS/structures)

## Where can I find more information?

   - Link to [paper]()
   - Link to an [interactive heatmap display](https://jbloomlab.github.io/SARS-CoV-2-RBD_DMS/) of the underlying data
   - Link to [GitHub repo](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS)
