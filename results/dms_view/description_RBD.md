# Deep mutational scanning of SARS-CoV-2 receptor binding domain reveals constraints on folding and ACE2 binding

Tyler N. Starr, Allison J. Greaney, Sarah K. Hilton, Katharine H.D. Crawford1,  ..., David Veesler, Jesse D. Bloom


## Abstract
The receptor binding domain (RBD) of the SARS-CoV-2 spike mediates viral attachment to ACE2 receptor, and so is a major determinant of host range and a dominant target of neutralizing antibodies. Here we experimentally measure how all amino acid mutations to the RBD affect expression of folded protein and its affinity for ACE2. Most mutations are deleterious for expression and ACE2 binding, and we identify highly constrained regions on the RBD’s surface that may be good targets for broadly neutralizing antibodies. But a substantial number of mutations are well tolerated or even enhance ACE2 binding, including at interface residues that vary across SARS-related coronaviruses. However, we find no evidence that these ACE2-affinity enhancing mutations have been selected in current SARS-CoV-2 pandemic isolates. We present an interactive visualization and open analysis pipeline to facilitate use of our comprehensive dataset for vaccine design and functional annotation of mutations observed during viral surveillance.

## What data is shown here?

We are showing the effects of all mutations to the SARS-CoV-2 RBD as measured in our [study](). Menu bars can alternate the view between our measurements of ACE2 receptor binding affinity and RBD expression (a correlate of protein stability). Within each phenotype, different metrics can be visualized at the per-site and per-mutation level:

   - per-site entropy and n\_effective describe the degree of constraint reflected in the amino acid preferences illustrated in logoplots as in Figure S4 of the manuscript.
   - per-site mean-effect, min-effect, and max-effect are statistics on the raw measurement values of mutational effects on binding and expression, as reflected in the heatmaps in Figure 3.
   - per-mutation preference is a relative preference value useful for visualization in logoplots, as the sum of preferences for all 20 amino acids at a site sum to 1.
   - per-mutation delta\_effect is the raw measurement value of mutational effect on binding and expression, either delta-log<sub>10</sub>(_K_<sub>D,app</sub>) for binding, or delta-log(MFI) for expression.

When sites are selected, their location is indicated on the ACE2-bound RBD crystal structure (PDB [6M0J](https://www.rcsb.org/structure/6M0J), from [Lan _et al._] 2020(https://www.nature.com/articles/s41586-020-2180-5))
   
## Where can I find more information?

   - Link to [paper]()
   - Link to [GitHub repo](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS)