#8 April 2020
#TNS

#script to read in PDB files and annotate structural contacts with receptor and Abs

setwd("/fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARS-CoV-2-RBD_DMS")

library(bio3d)

#read in file that gives concordance between numbering of SARS-CoV-2 and SARS-CoV RBDs
sites_alignment <- read.csv(file="data/RBD_sites.csv",stringsAsFactors=F)

#read in PDB of 6m0j, gives SARS-CoV-2 RBD (chain E) bound to huACE2 (chain A)
RBD2_ACE2 <- read.pdb(file="data/structures/ACE2-bound/6m0j.pdb")
RBD2_ACE2_atoms <- RBD2_ACE2$atom
#read in same PDB but with the ACE2 atoms removed, to calculate buried surface area
RBD2 <- read.pdb(file="data/structures/ACE2-bound/6m0j_no-ACE2.pdb")
RBD2_atoms <- RBD2$atom



