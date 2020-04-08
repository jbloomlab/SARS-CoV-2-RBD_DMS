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

#calculate RSA of free and bound RBD, distance to nearest ACE2 contact residue, for each RBD position
#installed dssp executable to ./src with wget ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.0.4-linux-amd64 -O ./src/dssp and chmod a+x ./src/dssp
dssp_RBD2_ACE2 <- dssp(RBD2_ACE2, exefile="~/bloom_j/computational_notebooks/tstarr/2020/SARS-CoV-2-RBD_DMS/src/dssp")
dssp_RBD2 <- dssp(RBD2, exefile="~/bloom_j/computational_notebooks/tstarr/2020/SARS-CoV-2-RBD_DMS/src/dssp")

#seems like the dssp output does not automatically link chain, resi number, and the SASA and SSE stats? Annoying.
RBD2_ACE2_dssp_summary <- data.frame(chain=rep(NA,length(dssp_RBD2_ACE2$sse)),resi=rep(NA,length(dssp_RBD2_ACE2$sse)),SSE=rep(NA,length(dssp_RBD2_ACE2$sse)),SASA=rep(NA,length(dssp_RBD2_ACE2$sse)))

for(i in 1:length(dssp_RBD2_ACE2$sse)){
  RBD2_ACE2_dssp_summary$resi[i] <- as.numeric(strsplit(names(dssp_RBD2_ACE2$sse)[i],split="_")[[1]][1])
  RBD2_ACE2_dssp_summary$chain[i] <- strsplit(names(dssp_RBD2_ACE2$sse)[i],split="_")[[1]][2]
  RBD2_ACE2_dssp_summary$SSE[i] <- dssp_RBD2_ACE2$sse[i]
  RBD2_ACE2_dssp_summary$SASA[i]<- dssp_RBD2_ACE2$acc[i]
}

#repeat to get SASA in structure without ACE2
RBD2_dssp_summary <- data.frame(chain=rep(NA,length(dssp_RBD2$sse)),resi=rep(NA,length(dssp_RBD2$sse)),SASA=rep(NA,length(dssp_RBD2$sse)))
for(i in 1:length(dssp_RBD2$sse)){
  RBD2_dssp_summary$resi[i] <- as.numeric(strsplit(names(dssp_RBD2$sse)[i],split="_")[[1]][1])
  RBD2_dssp_summary$chain[i] <- strsplit(names(dssp_RBD2$sse)[i],split="_")[[1]][2]
  RBD2_dssp_summary$SASA[i]<- dssp_RBD2$acc[i]
}


for(i in 1:nrow(sites_alignment)){
  if(sites_alignment$site_SARS2[i] %in% RBD2_ACE2_dssp_summary[RBD2_ACE2_dssp_summary$chain=="E","resi"]){
    sites_alignment$sse[i] <- RBD2_ACE2_dssp_summary[RBD2_ACE2_dssp_summary$chain=="E" & RBD2_ACE2_dssp_summary$resi == sites_alignment[i,"site_SARS2"],"SSE"]
    sites_alignment$SASA_bound[i] <- RBD2_ACE2_dssp_summary[RBD2_ACE2_dssp_summary$chain=="E" & RBD2_ACE2_dssp_summary$resi == sites_alignment[i,"site_SARS2"],"SASA"]
    sites_alignment$SASA_unbound[i] <- RBD2_dssp_summary[RBD2_dssp_summary$chain=="E" & RBD2_dssp_summary$resi == sites_alignment[i,"site_SARS2"],"SASA"]
  }else{
    sites_alignment$sse[i] <- NA
    sites_alignment$SASA_bound[i] <- NA
    sites_alignment$SASA_unbound[i] <- NA
  }
}

#calculate RSA by normalizing SASA by maximum SASA of different amino acid types. Max from Tien et al.
max_SASA <- data.frame(AA=c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"),max_SASA=c(121,265,187,187,148,214,214,97,216,195,191,230,203,228,154,143,163,264,255,165))

for(i in 1:nrow(sites_alignment)){
  sites_alignment$RSA_bound[i] <- sites_alignment$SASA_bound[i]/max_SASA[max_SASA$AA==sites_alignment$amino_acid_SARS2[i],"max_SASA"]
  sites_alignment$RSA_unbound[i] <- sites_alignment$SASA_unbound[i]/max_SASA[max_SASA$AA==sites_alignment$amino_acid_SARS2[i],"max_SASA"]
}

#annotate as contact residue if within (5 Angstrom) from ACE2 atoms

#compare to annotated contact residues from the 6m0j paper -- how did they define contacts?