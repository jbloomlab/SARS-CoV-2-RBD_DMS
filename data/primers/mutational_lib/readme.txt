21 March 2020 TNS
create mutagenesis primers for SARS-CoV-2 RBD

copied ../../2019/NIH45-46_primers/create_NNSprimers_py3.py to this directory

created SARS-CoV-2_RBD.txt file, which contains the mutagenic amplicon sequence for SARS-CoV-2 RBD as amplified from p2649, with the desired codons to be mutagenized in upper case

python create_NNSprimers_py3.py SARS-CoV-2_RBD.txt SARS-CoV-2_RBD 331 SARS-CoV-2_RBD_NNSprimers.txt

Going to order as IDT oPools (all oligos already pooled) rather than separate primers on plates

Duplicated output into forward and reverse primer pools

oPools only allow N and K degenerate bases. To keep NNS, I just duplicated each primer, and synthesize it as an NNG and an NNC separately (GNN and CNN for the reverse primer)

Collated into the IDT order form for upload, ordered two pools for forward and reverse primer sets at 1pmol each primer scale
