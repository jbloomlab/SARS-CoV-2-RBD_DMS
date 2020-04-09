#8 April 2020
#TNS

#script to read in PDB files and annotate structural contacts with receptor and Abs

setwd("/fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARS-CoV-2-RBD_DMS")

library(bio3d)

#read in file that gives concordance between numbering of SARS-CoV-2 and SARS-CoV RBDs
sites_alignment <- read.csv(file="data/RBD_sites.csv",stringsAsFactors=F)

#annotate residues as belonging to the RBM (SARS-CoV-2 residues N437-Y508)
sites_alignment$RBM <- F
sites_alignment[sites_alignment$site_SARS2 %in% seq(437,508),"RBM"] <- T

#annotate residues of paired disulfides: (C336/C361, C379/C432, C391/C525, C480/C488)
sites_alignment$disulfide <- NA
sites_alignment[sites_alignment$site_SARS2 %in% c(336, 361),"disulfide"] <- 1
sites_alignment[sites_alignment$site_SARS2 %in% c(379, 432),"disulfide"] <- 2
sites_alignment[sites_alignment$site_SARS2 %in% c(391, 525),"disulfide"] <- 3
sites_alignment[sites_alignment$site_SARS2 %in% c(480, 488),"disulfide"] <- 4

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

#annotate as contact residue if within (4 Angstrom) from ACE2 atoms
#use binding.site function from bio3d to determine interaction residues as <4 Angstrom distance
RBD2_ACE2_contacts <- binding.site(RBD2_ACE2, a.inds=atom.select(RBD2_ACE2, chain="E"), b.inds=atom.select(RBD2_ACE2, chain="A"),cutoff=4, hydrogens=F)
#this gives exact same set of contact residues as the 6m0j paper
sites_alignment$SARS2_ACE2_contact <- F
sites_alignment[sites_alignment$site_SARS2 %in% RBD2_ACE2_contacts$resno,"SARS2_ACE2_contact"] <- T

#repeat to annotate SARS1 RBD ACE2 contacts with 2ajf structure, chain A (ACE2) and chain E (SARS-CoV RBD)
RBD1_ACE2 <- read.pdb(file="data/structures/ACE2-bound/2ajf.pdb")
RBD1_ACE2_contacts <- binding.site(RBD1_ACE2, a.inds=atom.select(RBD1_ACE2, chain="E"), b.inds=atom.select(RBD1_ACE2, chain="A"),cutoff=4, hydrogens=F)
sites_alignment$SARS1_ACE2_contact <- F
sites_alignment[sites_alignment$site_SARS1 %in% RBD1_ACE2_contacts$resno & !is.na(sites_alignment$site_SARS1),"SARS1_ACE2_contact"] <- T

#'key 5 sites' from Fang Li, etc. biochemical and human/civet adaptation experiments: SARS2 residues L455, F486, Q493, S494, N501, (Y505 also highlighted in the Anderson review)
sites_alignment$SARS1_key_adaptation <- F
sites_alignment[sites_alignment$site_SARS2 %in% c(455, 486, 493, 494, 501),"SARS1_key_adaptation"] <- T

#next, annotate residues as being in epitopes for known SARS1 and SARS2 antibodies
#CR3022 Fab bound to SARS2 RBD, pdb 6w41 (chain C = RBD, chains H and L = Fab)
RBD2_CR3022 <- read.pdb(file="data/structures/Ab-bound/CR3022_6w41.pdb")
RBD2_CR3022_contacts <- binding.site(RBD2_CR3022, a.inds=atom.select(RBD2_CR3022, chain="C"), b.inds=atom.select(RBD2_CR3022, chain=c("H","L")),cutoff=4, hydrogens=F)
sites_alignment$epitope_CR3022 <- F
sites_alignment[sites_alignment$site_SARS2 %in% RBD2_CR3022_contacts$resno,"epitope_CR3022"] <- T

#VHH-72 camelid nanobody bound to SARS1 RBD, pdb 6waq (chain B = RBD, chain A = nanobody)
RBD1_VHH72 <- read.pdb(file="data/structures/Ab-bound/VHH-72_6waq.pdb")
RBD1_VHH72_contacts <- binding.site(RBD1_VHH72, a.inds=atom.select(RBD1_VHH72, chain="B"), b.inds=atom.select(RBD1_VHH72, chain=c("A")),cutoff=4, hydrogens=F)
sites_alignment$epitope_VHH72 <- F
sites_alignment[sites_alignment$site_SARS1 %in% RBD1_VHH72_contacts$resno & !is.na(sites_alignment$site_SARS1),"epitope_VHH72"] <- T

#S230 Fab bound to SARS1 prefusion trimer, pdb 6nb6 (chain C = spike monomer, interacts with chains H_L Fab, also interaction chains B/M+I) and 6nb7 (chain A=spike monomer, chain H and L=Fab; other trimers are B/D+E and C/G+I chain interactions)
#this is an EM model with lower resolution -- only backbone atoms appear to be modeled. Therefore, do 8A cutoff for more standard Calpha cutoff. Keep only interactoins that appear in majority of modeled protomer interactions?
RBD1_S230_1 <- read.pdb(file="data/structures/Ab-bound/S230_6nb6.pdb")
RBD1_S230_1_contacts <- binding.site(RBD1_S230_1, a.inds=atom.select(RBD1_S230_1, chain="C"), b.inds=atom.select(RBD1_S230_1, chain=c("H","L")),cutoff=8, hydrogens=F)
RBD1_S230_1_contacts2 <- binding.site(RBD1_S230_1, a.inds=atom.select(RBD1_S230_1, chain="B"), b.inds=atom.select(RBD1_S230_1, chain=c("M","I")),cutoff=8, hydrogens=F)

RBD1_S230_2 <- read.pdb(file="data/structures/Ab-bound/S230_6nb7.pdb")
RBD1_S230_2_contacts <- binding.site(RBD1_S230_2, a.inds=atom.select(RBD1_S230_2, chain="A"), b.inds=atom.select(RBD1_S230_2, chain=c("H","L")),cutoff=8, hydrogens=F)
RBD1_S230_2_contacts2 <- binding.site(RBD1_S230_2, a.inds=atom.select(RBD1_S230_2, chain="B"), b.inds=atom.select(RBD1_S230_2, chain=c("D","E")),cutoff=8, hydrogens=F)
RBD1_S230_2_contacts3 <- binding.site(RBD1_S230_2, a.inds=atom.select(RBD1_S230_2, chain="C"), b.inds=atom.select(RBD1_S230_2, chain=c("G","I")),cutoff=8, hydrogens=F)

#keep only residues that are contacts in >3 of the 5 modeled protomers
contacts <- unique(RBD1_S230_1_contacts$resno,RBD1_S230_1_contacts2$resno,RBD1_S230_2_contacts$resno,RBD1_S230_2_contacts2$resno,RBD1_S230_2_contacts3$resno)
for(i in 1:length(contacts)){
  if(sum(contacts[i] == c(RBD1_S230_1_contacts$resno,RBD1_S230_1_contacts2$resno,RBD1_S230_2_contacts$resno,RBD1_S230_2_contacts2$resno,RBD1_S230_2_contacts3$resno))<4){
    contacts[i] <- NA
  }
}

sites_alignment$epitope_S230 <- F
sites_alignment[sites_alignment$site_SARS1 %in% contacts  & !is.na(sites_alignment$site_SARS1),"epitope_S230"] <- T

#m396 Fab bound to SARS1 RBD, pdb 2dd8 (chain S = RBD, chain H+L = Fab)
RBD1_m396 <- read.pdb(file="data/structures/Ab-bound/m396_2dd8.pdb")
RBD1_m396_contacts <- binding.site(RBD1_m396, a.inds=atom.select(RBD1_m396, chain="S"), b.inds=atom.select(RBD1_m396, chain=c("H","L")),cutoff=4, hydrogens=F)
sites_alignment$epitope_m396 <- F
sites_alignment[sites_alignment$site_SARS1 %in% RBD1_m396_contacts$resno & !is.na(sites_alignment$site_SARS1),"epitope_m396"] <- T

#F26G19 Fab bound to SARS1 RBD, pdb 3bgf (chain S = RBD, chain H+L = Fab)
RBD1_F26G19 <- read.pdb(file="data/structures/Ab-bound/F26G19_3bgf.pdb")
RBD1_F26G19_contacts <- binding.site(RBD1_F26G19, a.inds=atom.select(RBD1_F26G19, chain="S"), b.inds=atom.select(RBD1_F26G19, chain=c("H","L")),cutoff=4, hydrogens=F)
sites_alignment$epitope_F26G19 <- F
sites_alignment[sites_alignment$site_SARS1 %in% RBD1_F26G19_contacts$resno & !is.na(sites_alignment$site_SARS1),"epitope_F26G19"] <- T

#80R scFv bound to SARS1 RBD, pdb 2ghw (chain A = RBD, chain B = scFv)
RBD1_80R <- read.pdb(file="data/structures/Ab-bound/80R_2ghw.pdb")
RBD1_80R_contacts <- binding.site(RBD1_80R, a.inds=atom.select(RBD1_80R, chain="A"), b.inds=atom.select(RBD1_80R, chain=c("B")),cutoff=4, hydrogens=F)
sites_alignment$epitope_80R <- F
sites_alignment[sites_alignment$site_SARS1 %in% RBD1_80R_contacts$resno & !is.na(sites_alignment$site_SARS1),"epitope_80R"] <- T

#save RBD_sites.csv with structural annotations for downstream data analysis
write.csv(sites_alignment,file="data/RBD_sites.csv",row.names=F)



