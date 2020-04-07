#21 November 2019
#TNS

#script to analyze outputs of global epistasis analyses for binding and expression phenotypes

setwd("~/bloom_j/computational_notebooks/tstarr/2019/NIH45-46_DMS")

#import R libraries
library(yaml)
library(ggplot2)
library(data.table)
#library(Rpdb)
library(bio3d)
#library(seqinr)

#read the configuration file
config <- read_yaml("config.yaml")

#make output files/directories
#make output directory
if(!file.exists(config$analyze_mut_effects_dir)){
  dir.create(file.path(config$analyze_mut_effects_dir))
}

#read in PDB of 3u7y, NIH45-45 fab bound to gp120, make dataframe of atoms and coordinates, pasting together "resid" and "insert" for residue numbering
pdb <- read.pdb(file="data/3u7y.pdb")
pdb_atoms <- pdb$atom
pdb_atoms[!is.na(pdb_atoms$insert),"resno"] <- paste(pdb_atoms[!is.na(pdb_atoms$insert),"resno"],pdb_atoms[!is.na(pdb_atoms$insert),"insert"],sep="")
#read in same PDB but without the gp120 molecule, to calculate buried surface area
pdb_fab <- read.pdb(file="data/3u7y_fab-only.pdb")
pdb_fab_atoms <- pdb_fab$atom
pdb_fab_atoms[!is.na(pdb_fab_atoms$insert),"resno"] <- paste(pdb_fab_atoms[!is.na(pdb_fab_atoms$insert),"resno"],pdb_fab_atoms[!is.na(pdb_fab_atoms$insert),"insert"],sep="")

#read in file that gives concordance between numbering on the scFv construct and the Fab numbering (heavy/light chain, as well as loop ABC numbering scheme, both per the pdb file I use as well as Kabat nomenclature)
sites_alignment <- read.csv(file="data/NIH45-46_sites.csv",stringsAsFactors=F)

#read in tables of per-barcode measured and predicted phenotypes
bc_bind <- data.table(read.csv(file="results/global_epistasis_binding/globalepistasis_binding_predictions_wide.csv",stringsAsFactors = F))
bc_expr <- data.table(read.csv(file="results/global_epistasis_expression/globalepistasis_expression_predictions_wide.csv",stringsAsFactors = F))
#can reproduce plots in the global epistasis notebook, e.g.
#libA, Cauchy
#png("results/analyze_mut_effects/libA2-bind_ddG-v-latent-phenotype_with_red-line.png",bg="transparent",width=330,height=360)
plot(bc_bind[library=="libA2",latent_phenotype_Cauchy_A],bc_bind[library=="libA2",ddG],pch=16,col="#00000005",xlab="global epistasis latent phenotype",ylab="Tite-seq ddG (kcal/mol)")
points(bc_bind[library=="libA2",latent_phenotype_Cauchy_A][order(bc_bind[library=="libA2",latent_phenotype_Cauchy_A])],bc_bind[library=="libA2",predicted_phenotype_Cauchy_A][order(bc_bind[library=="libA2",latent_phenotype_Cauchy_A])],type="l",col="red",lwd=3)
#dev.off()
#libB, Cauchy
#png("results/analyze_mut_effects/libB2-bind_ddG-v-latent-phenotype_with_red-line.png",bg="transparent",width=330,height=360)
plot(bc_bind[library=="libB2",latent_phenotype_Cauchy_B],bc_bind[library=="libB2",ddG],pch=16,col="#00000005",xlab="global epistasis latent phenotype",ylab="Tite-seq ddG (kcal/mol)")
points(bc_bind[library=="libB2",latent_phenotype_Cauchy_B][order(bc_bind[library=="libB2",latent_phenotype_Cauchy_B])],bc_bind[library=="libB2",predicted_phenotype_Cauchy_B][order(bc_bind[library=="libB2",latent_phenotype_Cauchy_B])],type="l",col="red",lwd=3)
#dev.off()
#libA, Gaussian
plot(bc_bind[library=="libA2",latent_phenotype_Gaussian_A],bc_bind[library=="libA2",ddG],pch=16,col="#00000005",xlab="global epistasis latent phenotype",ylab="Tite-seq ddG (kcal/mol)")
points(bc_bind[library=="libA2",latent_phenotype_Gaussian_A][order(bc_bind[library=="libA2",latent_phenotype_Gaussian_A])],bc_bind[library=="libA2",predicted_phenotype_Gaussian_A][order(bc_bind[library=="libA2",latent_phenotype_Gaussian_A])],type="l",col="red",lwd=3)
#libB, Gaussian
plot(bc_bind[library=="libB2",latent_phenotype_Gaussian_B],bc_bind[library=="libB2",ddG],pch=16,col="#00000005",xlab="global epistasis latent phenotype",ylab="Tite-seq ddG (kcal/mol)")
points(bc_bind[library=="libB2",latent_phenotype_Gaussian_B][order(bc_bind[library=="libB2",latent_phenotype_Gaussian_B])],bc_bind[library=="libB2",predicted_phenotype_Gaussian_B][order(bc_bind[library=="libB2",latent_phenotype_Gaussian_B])],type="l",col="red",lwd=3)
#expression phenotype
#png("results/analyze_mut_effects/libA2-expr_ML_meanF-v-latent-phenotype_with_red-line.png",bg="transparent",width=330,height=360)
plot(bc_expr[library=="libA2",latent_phenotype_Cauchy_A],bc_expr[library=="libA2",ML_meanF],pch=16,col="#00000005",xlab="global epistasis latent phenotype",ylab="Sort-seq ML meanF (a.u)")
points(bc_expr[library=="libA2",latent_phenotype_Cauchy_A][order(bc_expr[library=="libA2",latent_phenotype_Cauchy_A])],bc_expr[library=="libA2",predicted_phenotype_Cauchy_A][order(bc_expr[library=="libA2",latent_phenotype_Cauchy_A])],type="l",col="red",lwd=3)
#dev.off()
#libB, Cauchy
#png("results/analyze_mut_effects/libB2-expr_ML_meanF-v-latent-phenotype_with_red-line.png",bg="transparent",width=330,height=360)
plot(bc_expr[library=="libB2",latent_phenotype_Cauchy_B],bc_expr[library=="libB2",ML_meanF],pch=16,col="#00000005",xlab="global epistasis latent phenotype",ylab="Sort-seq ML meanF (a.u.)")
points(bc_expr[library=="libB2",latent_phenotype_Cauchy_B][order(bc_expr[library=="libB2",latent_phenotype_Cauchy_B])],bc_expr[library=="libB2",predicted_phenotype_Cauchy_B][order(bc_expr[library=="libB2",latent_phenotype_Cauchy_B])],type="l",col="red",lwd=3)
#dev.off()


#binding measurements: read in coefficients from each model, look at correlations between replicates. For Cauchy, compare to mutant effects from model without global epistasis correction
#observed scale, Cauchy likelihood
betas_bind_observed_Cauchy_A <- read.csv('results/global_epistasis_binding/Cauchy-predicted-effects_binding_A.csv',stringsAsFactors=F)
betas_bind_observed_Cauchy_B <- read.csv('results/global_epistasis_binding/Cauchy-predicted-effects_binding_B.csv',stringsAsFactors=F)
betas_bind_observed_Cauchy <- merge(betas_bind_observed_Cauchy_A,betas_bind_observed_Cauchy_B,by=c("site","mutation","wildtype","mutant"),all=T,sort=T,suffixes=c("_A","_B"));rm(betas_bind_observed_Cauchy_A);rm(betas_bind_observed_Cauchy_B)
#pdf(file="results/analyze_mut_effects/correlation-betas_bind-observed.pdf",bg="transparent",width=4,height=4.5,useDingbats=F)
plot(betas_bind_observed_Cauchy$effect_A[!betas_bind_observed_Cauchy$site %in% seq(126,140)],betas_bind_observed_Cauchy$effect_B[!betas_bind_observed_Cauchy$site %in% seq(126,140)],pch=16,col="#00000067",xlab="rep A beta (observed scale)",ylab="rep B beta (observed scale)",main="R-squared 0.9422");summary(lm(betas_bind_observed_Cauchy$effect_B[!betas_bind_observed_Cauchy$site %in% seq(126,140)]~betas_bind_observed_Cauchy$effect_A[!betas_bind_observed_Cauchy$site %in% seq(126,140)]));abline(lm(betas_bind_observed_Cauchy$effect_B[!betas_bind_observed_Cauchy$site %in% seq(126,140)]~betas_bind_observed_Cauchy$effect_A[!betas_bind_observed_Cauchy$site %in% seq(126,140)]),col="red")
#dev.off()
#R-squared 0.9422
betas_bind_observed_Cauchy[betas_bind_observed_Cauchy$effect_B<1 & betas_bind_observed_Cauchy$effect_A>6,]
betas_bind_observed_Cauchy[betas_bind_observed_Cauchy$effect_B>5 & betas_bind_observed_Cauchy$effect_A<2,]

#latent scale, Cauchy likelihood
betas_bind_latent_Cauchy_A <- read.csv('results/global_epistasis_binding/Cauchy-latent-effects_binding_A.csv',stringsAsFactors=F)
betas_bind_latent_Cauchy_B <- read.csv('results/global_epistasis_binding/Cauchy-latent-effects_binding_B.csv',stringsAsFactors=F)
betas_bind_latent_Cauchy <- merge(betas_bind_latent_Cauchy_A,betas_bind_latent_Cauchy_B,by=c("site","mutation","wildtype","mutant"),all=T,sort=T,suffixes=c("_A","_B"));rm(betas_bind_latent_Cauchy_A);rm(betas_bind_latent_Cauchy_B)
#pdf(file="results/analyze_mut_effects/correlation-betas_bind-latent.pdf",bg="transparent",width=4,height=4.5)
plot(betas_bind_latent_Cauchy$effect_A[!betas_bind_latent_Cauchy$site %in% seq(126,140)],betas_bind_latent_Cauchy$effect_B[!betas_bind_latent_Cauchy$site %in% seq(126,140)],pch=16,col="#00000067",xlab="rep A beta (latent scale)",ylab="rep B beta (latent scale)",main="R-squared 0.947");summary(lm(betas_bind_latent_Cauchy$effect_B[!betas_bind_latent_Cauchy$site %in% seq(126,140)]~betas_bind_latent_Cauchy$effect_A[!betas_bind_latent_Cauchy$site %in% seq(126,140)]));abline(lm(betas_bind_latent_Cauchy$effect_B[!betas_bind_latent_Cauchy$site %in% seq(126,140)]~betas_bind_latent_Cauchy$effect_A[!betas_bind_latent_Cauchy$site %in% seq(126,140)]),col="red")
#dev.off()
#R-squared 0.9472
#can see the G54W (site 55 in scFv) is indeed beneficial, but there are some stronger ones?
points(betas_bind_latent_Cauchy[betas_bind_latent_Cauchy$mutation=="G55W","effect_A"],betas_bind_latent_Cauchy[betas_bind_latent_Cauchy$mutation=="G55W","effect_B"],pch=20,col="green")

#without global epistasis correction
betas_bind_observed_nonepistatic_Cauchy_A <- read.csv('results/global_epistasis_binding/nonepistatic-Cauchy-predicted-effects_binding_A.csv',stringsAsFactors=F)
betas_bind_observed_nonepistatic_Cauchy_B <- read.csv('results/global_epistasis_binding/nonepistatic-Cauchy-predicted-effects_binding_B.csv',stringsAsFactors=F)
betas_bind_observed_nonepistatic_Cauchy <- merge(betas_bind_observed_nonepistatic_Cauchy_A,betas_bind_observed_nonepistatic_Cauchy_B,by=c("site","mutation","wildtype","mutant"),all=T,sort=T,suffixes=c("_A","_B"));rm(betas_bind_observed_nonepistatic_Cauchy_A);rm(betas_bind_observed_nonepistatic_Cauchy_B)
#pdf(file="results/analyze_mut_effects/correlation-betas_bind-nonepistatic.pdf",bg="transparent",width=4,height=4.5)
plot(betas_bind_observed_nonepistatic_Cauchy$effect_A[!betas_bind_observed_nonepistatic_Cauchy$site %in% seq(126,140)],betas_bind_observed_nonepistatic_Cauchy$effect_B[!betas_bind_observed_nonepistatic_Cauchy$site %in% seq(126,140)],pch=16,col="#00000067",xlab="rep A beta (no epistasis correction)",ylab="rep B beta (no epistasis correction)",main="R-squared 0.927");summary(lm(betas_bind_observed_nonepistatic_Cauchy$effect_B[!betas_bind_observed_nonepistatic_Cauchy$site %in% seq(126,140)]~betas_bind_observed_nonepistatic_Cauchy$effect_A[!betas_bind_observed_nonepistatic_Cauchy$site %in% seq(126,140)]));abline(lm(betas_bind_observed_nonepistatic_Cauchy$effect_B[!betas_bind_observed_nonepistatic_Cauchy$site %in% seq(126,140)]~betas_bind_observed_nonepistatic_Cauchy$effect_A[!betas_bind_observed_nonepistatic_Cauchy$site %in% seq(126,140)]),col="red")
#dev.off()
#R-squared 0.9272
plot(rowMeans(betas_bind_observed_nonepistatic_Cauchy[,c("effect_A","effect_B")],na.rm=T),rowMeans(betas_bind_observed_Cauchy[,c("effect_A","effect_B")]),pch=16,col="#00000067")
plot(rowMeans(betas_bind_observed_nonepistatic_Cauchy[,c("effect_A","effect_B")],na.rm=T),rowMeans(betas_bind_latent_Cauchy[,c("effect_A","effect_B")]),pch=16,col="#00000067")


#observed scale, Gaussian likelihood
betas_bind_observed_Gaussian_A <- read.csv('results/global_epistasis_binding/Gaussian-predicted-effects_binding_A.csv',stringsAsFactors=F)
betas_bind_observed_Gaussian_B <- read.csv('results/global_epistasis_binding/Gaussian-predicted-effects_binding_B.csv',stringsAsFactors=F)
betas_bind_observed_Gaussian <- merge(betas_bind_observed_Gaussian_A,betas_bind_observed_Gaussian_B,by=c("site","mutation","wildtype","mutant"),all=T,sort=T,suffixes=c("_A","_B"));rm(betas_bind_observed_Gaussian_A);rm(betas_bind_observed_Gaussian_B)
plot(betas_bind_observed_Gaussian$effect_A,betas_bind_observed_Gaussian$effect_B,pch=16,col="#00000067");summary(lm(betas_bind_observed_Gaussian$effect_B~betas_bind_observed_Gaussian$effect_A));abline(lm(betas_bind_observed_Gaussian$effect_B~betas_bind_observed_Gaussian$effect_A),col="red")
#R-squared 0.9168

#latent scale, Gaussian likelihood
betas_bind_latent_Gaussian_A <- read.csv('results/global_epistasis_binding/Gaussian-latent-effects_binding_A.csv',stringsAsFactors=F)
betas_bind_latent_Gaussian_B <- read.csv('results/global_epistasis_binding/Gaussian-latent-effects_binding_B.csv',stringsAsFactors=F)
betas_bind_latent_Gaussian <- merge(betas_bind_latent_Gaussian_A,betas_bind_latent_Gaussian_B,by=c("site","mutation","wildtype","mutant"),all=T,sort=T,suffixes=c("_A","_B"));rm(betas_bind_latent_Gaussian_A);rm(betas_bind_latent_Gaussian_B)
plot(betas_bind_latent_Gaussian$effect_A,betas_bind_latent_Gaussian$effect_B,pch=16,col="#00000067");summary(lm(betas_bind_latent_Gaussian$effect_B~betas_bind_latent_Gaussian$effect_A));abline(lm(betas_bind_latent_Gaussian$effect_B~betas_bind_latent_Gaussian$effect_A),col="red")
#R-squared 0.8981

#Cauchy likelihood seems to give better results -- is R-squared an appropriate metric for deciding to go ahead with Cauchy?


#expression measurements: read in coefficients from each model, look at correlations between replicates
#observed scale, Cauchy likelihood
betas_expr_observed_Cauchy_A <- read.csv('results/global_epistasis_expression/Cauchy-predicted-effects_expression_A.csv',stringsAsFactors=F)
betas_expr_observed_Cauchy_B <- read.csv('results/global_epistasis_expression/Cauchy-predicted-effects_expression_B.csv',stringsAsFactors=F)
betas_expr_observed_Cauchy <- merge(betas_expr_observed_Cauchy_A,betas_expr_observed_Cauchy_B,by=c("site","mutation","wildtype","mutant"),all=T,sort=T,suffixes=c("_A","_B"));rm(betas_expr_observed_Cauchy_A);rm(betas_expr_observed_Cauchy_B)
plot(betas_expr_observed_Cauchy$effect_A[!betas_expr_observed_Cauchy$site %in% seq(126,140)],betas_expr_observed_Cauchy$effect_B[!betas_expr_observed_Cauchy$site %in% seq(126,140)],pch=16,col="#00000067");summary(lm(betas_expr_observed_Cauchy$effect_B[!betas_expr_observed_Cauchy$site %in% seq(126,140)]~betas_expr_observed_Cauchy$effect_A[!betas_expr_observed_Cauchy$site %in% seq(126,140)]));abline(lm(betas_expr_observed_Cauchy$effect_B[!betas_expr_observed_Cauchy$site %in% seq(126,140)]~betas_expr_observed_Cauchy$effect_A[!betas_expr_observed_Cauchy$site %in% seq(126,140)]),col="red")
#excluding STOP variants?
#pdf(file="results/analyze_mut_effects/correlation-betas_expr-observed.pdf",bg="transparent",width=4,height=4.5)
plot(betas_expr_observed_Cauchy$effect_A[!betas_expr_observed_Cauchy$site %in% seq(126,140) & betas_expr_observed_Cauchy$mutant!="*"],betas_expr_observed_Cauchy$effect_B[!betas_expr_observed_Cauchy$site %in% seq(126,140) & betas_expr_observed_Cauchy$mutant!="*"],pch=16,col="#00000067",xlab="rep A expression beta (observed scale)",ylab="reb B expression beta (observed scale)",main="R-squared 0.908");summary(lm(betas_expr_observed_Cauchy$effect_B[!betas_expr_observed_Cauchy$site %in% seq(126,140) & betas_expr_observed_Cauchy$mutant!="*"]~betas_expr_observed_Cauchy$effect_A[!betas_expr_observed_Cauchy$site %in% seq(126,140) & betas_expr_observed_Cauchy$mutant!="*"]));abline(lm(betas_expr_observed_Cauchy$effect_B[!betas_expr_observed_Cauchy$site %in% seq(126,140) & betas_expr_observed_Cauchy$mutant!="*"]~betas_expr_observed_Cauchy$effect_A[!betas_expr_observed_Cauchy$site %in% seq(126,140) & betas_expr_observed_Cauchy$mutant!="*"]),col="red")
#dev.off()
#R-squared 0.9142, 0.9078 w/o STOP

#latent scale, Cauchy likelihood
betas_expr_latent_Cauchy_A <- read.csv('results/global_epistasis_expression/Cauchy-latent-effects_expression_A.csv',stringsAsFactors=F)
betas_expr_latent_Cauchy_B <- read.csv('results/global_epistasis_expression/Cauchy-latent-effects_expression_B.csv',stringsAsFactors=F)
betas_expr_latent_Cauchy <- merge(betas_expr_latent_Cauchy_A,betas_expr_latent_Cauchy_B,by=c("site","mutation","wildtype","mutant"),all=T,sort=T,suffixes=c("_A","_B"));rm(betas_expr_latent_Cauchy_A);rm(betas_expr_latent_Cauchy_B)
plot(betas_expr_latent_Cauchy$effect_A[!betas_expr_latent_Cauchy$site %in% seq(126,140)],betas_expr_latent_Cauchy$effect_B[!betas_expr_latent_Cauchy$site %in% seq(126,140)],pch=16,col="#00000067",xlab="mutant deltaE, replicate A",ylab="mutant deltaE, replicate B");summary(lm(betas_expr_latent_Cauchy$effect_B[!betas_expr_latent_Cauchy$site %in% seq(126,140)]~betas_expr_latent_Cauchy$effect_A[!betas_expr_latent_Cauchy$site %in% seq(126,140)]));abline(lm(betas_expr_latent_Cauchy$effect_B[!betas_expr_latent_Cauchy$site %in% seq(126,140)]~betas_expr_latent_Cauchy$effect_A[!betas_expr_latent_Cauchy$site %in% seq(126,140)]),col="red")
#excluding STOP variants?
#pdf(file="results/analyze_mut_effects/correlation-betas_expr-latent.pdf",bg="transparent",width=4,height=4.5)
plot(betas_expr_latent_Cauchy$effect_A[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*"],betas_expr_latent_Cauchy$effect_B[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*"],pch=16,col="#00000067",xlab="rep A expression beta (latent scale)",ylab="reb B expression beta (latent scale)",main="R-squared 0.857");summary(lm(betas_expr_latent_Cauchy$effect_B[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*"]~betas_expr_latent_Cauchy$effect_A[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*"]));abline(lm(betas_expr_latent_Cauchy$effect_B[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*"]~betas_expr_latent_Cauchy$effect_A[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*"]),col="red")
#dev.off()
#zoom in on unextrapolated region
plot(betas_expr_latent_Cauchy$effect_A[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*" & betas_expr_latent_Cauchy$effect_A> -2 & betas_expr_latent_Cauchy$effect_B > -2],betas_expr_latent_Cauchy$effect_B[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*" & betas_expr_latent_Cauchy$effect_A> -2 & betas_expr_latent_Cauchy$effect_B > -2],pch=16,col="#00000067",xlab="mutant deltaE, replicate A",ylab="mutant deltaE, replicate B");summary(lm(betas_expr_latent_Cauchy$effect_B[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*" & betas_expr_latent_Cauchy$effect_A> -2 & betas_expr_latent_Cauchy$effect_B > -2]~betas_expr_latent_Cauchy$effect_A[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*" & betas_expr_latent_Cauchy$effect_A> -2 & betas_expr_latent_Cauchy$effect_B > -2]));abline(lm(betas_expr_latent_Cauchy$effect_B[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*" & betas_expr_latent_Cauchy$effect_A> -2 & betas_expr_latent_Cauchy$effect_B > -2]~betas_expr_latent_Cauchy$effect_A[!betas_expr_latent_Cauchy$site %in% seq(126,140) & betas_expr_latent_Cauchy$mutant!="*" & betas_expr_latent_Cauchy$effect_A> -2 & betas_expr_latent_Cauchy$effect_B > -2]),col="red")
#R-squared 0.9414, 0.8571 w/o STOP, 0.8859 w/in "bulk" region of distribution (>-2 deltaE)

#without global epistasis correction
betas_expr_observed_nonepistatic_Cauchy_A <- read.csv('results/global_epistasis_expression/nonepistatic-Cauchy-predicted-effects_expression_A.csv',stringsAsFactors=F)
betas_expr_observed_nonepistatic_Cauchy_B <- read.csv('results/global_epistasis_expression/nonepistatic-Cauchy-predicted-effects_expression_B.csv',stringsAsFactors=F)
betas_expr_observed_nonepistatic_Cauchy <- merge(betas_expr_observed_nonepistatic_Cauchy_A,betas_expr_observed_nonepistatic_Cauchy_B,by=c("site","mutation","wildtype","mutant"),all=T,sort=T,suffixes=c("_A","_B"));rm(betas_expr_observed_nonepistatic_Cauchy_A);rm(betas_expr_observed_nonepistatic_Cauchy_B)
plot(betas_expr_observed_nonepistatic_Cauchy$effect_A[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140)],betas_expr_observed_nonepistatic_Cauchy$effect_B[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140)],pch=16,col="#00000067");summary(lm(betas_expr_observed_nonepistatic_Cauchy$effect_B[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140)]~betas_expr_observed_nonepistatic_Cauchy$effect_A[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140)]));abline(lm(betas_expr_observed_nonepistatic_Cauchy$effect_B[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140)]~betas_expr_observed_nonepistatic_Cauchy$effect_A[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140)]),col="red")
#excluding STOP variants?
#pdf(file="results/analyze_mut_effects/correlation-betas_expr-nonepistatic.pdf",bg="transparent",width=4,height=4.5)
plot(betas_expr_observed_nonepistatic_Cauchy$effect_A[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140) & betas_expr_observed_nonepistatic_Cauchy$mutant!="*"],betas_expr_observed_nonepistatic_Cauchy$effect_B[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140) & betas_expr_observed_nonepistatic_Cauchy$mutant!="*"],pch=16,col="#00000067",xlab="rep A expression beta (no epistasis correction)",ylab="reb B expression beta (no epistasis correction)",main="R-squared 0.831");summary(lm(betas_expr_observed_nonepistatic_Cauchy$effect_B[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140) & betas_expr_observed_nonepistatic_Cauchy$mutant!="*"]~betas_expr_observed_nonepistatic_Cauchy$effect_A[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140) & betas_expr_observed_nonepistatic_Cauchy$mutant!="*"]));abline(lm(betas_expr_observed_nonepistatic_Cauchy$effect_B[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140) & betas_expr_observed_nonepistatic_Cauchy$mutant!="*"]~betas_expr_observed_nonepistatic_Cauchy$effect_A[!betas_expr_observed_nonepistatic_Cauchy$site %in% seq(126,140) & betas_expr_observed_nonepistatic_Cauchy$mutant!="*"]),col="red")
#dev.off()
#R-squared 0.8796, 0.8313 w/o STOP
plot(rowMeans(betas_expr_observed_nonepistatic_Cauchy[,c("effect_A","effect_B")],na.rm=T),rowMeans(betas_expr_observed_Cauchy[,c("effect_A","effect_B")]),pch=16,col="#00000067")
plot(rowMeans(betas_expr_observed_nonepistatic_Cauchy[,c("effect_A","effect_B")],na.rm=T),rowMeans(betas_expr_latent_Cauchy[,c("effect_A","effect_B")]),pch=16,col="#00000067")


#observed scale, Gaussian likelihood
betas_expr_observed_Gaussian_A <- read.csv('results/global_epistasis_expression/Gaussian-predicted-effects_expression_A.csv',stringsAsFactors=F)
betas_expr_observed_Gaussian_B <- read.csv('results/global_epistasis_expression/Gaussian-predicted-effects_expression_B.csv',stringsAsFactors=F)
betas_expr_observed_Gaussian <- merge(betas_expr_observed_Gaussian_A,betas_expr_observed_Gaussian_B,by=c("site","mutation","wildtype","mutant"),all=T,sort=T,suffixes=c("_A","_B"));rm(betas_expr_observed_Gaussian_A);rm(betas_expr_observed_Gaussian_B)
plot(betas_expr_observed_Gaussian$effect_A[!betas_expr_observed_Gaussian$site %in% seq(126,140)],betas_expr_observed_Gaussian$effect_B[!betas_expr_observed_Gaussian$site %in% seq(126,140)],pch=16,col="#00000067");summary(lm(betas_expr_observed_Gaussian$effect_B[!betas_expr_observed_Gaussian$site %in% seq(126,140)]~betas_expr_observed_Gaussian$effect_A[!betas_expr_observed_Gaussian$site %in% seq(126,140)]));abline(lm(betas_expr_observed_Gaussian$effect_B[!betas_expr_observed_Gaussian$site %in% seq(126,140)]~betas_expr_observed_Gaussian$effect_A[!betas_expr_observed_Gaussian$site %in% seq(126,140)]),col="red")
#excluding STOP
plot(betas_expr_observed_Gaussian$effect_A[!betas_expr_observed_Gaussian$site %in% seq(126,140) & betas_expr_observed_Gaussian$mutant!="*"],betas_expr_observed_Gaussian$effect_B[!betas_expr_observed_Gaussian$site %in% seq(126,140) & betas_expr_observed_Gaussian$mutant!="*"],pch=16,col="#00000067");summary(lm(betas_expr_observed_Gaussian$effect_B[!betas_expr_observed_Gaussian$site %in% seq(126,140) & betas_expr_observed_Gaussian$mutant!="*"]~betas_expr_observed_Gaussian$effect_A[!betas_expr_observed_Gaussian$site %in% seq(126,140) & betas_expr_observed_Gaussian$mutant!="*"]));abline(lm(betas_expr_observed_Gaussian$effect_B[!betas_expr_observed_Gaussian$site %in% seq(126,140) & betas_expr_observed_Gaussian$mutant!="*"]~betas_expr_observed_Gaussian$effect_A[!betas_expr_observed_Gaussian$site %in% seq(126,140) & betas_expr_observed_Gaussian$mutant!="*"]),col="red")
#R-squared  0.907, 0.8596

#latent scale, Gaussian likelihood
betas_expr_latent_Gaussian_A <- read.csv('results/global_epistasis_expression/Gaussian-latent-effects_expression_A.csv',stringsAsFactors=F)
betas_expr_latent_Gaussian_B <- read.csv('results/global_epistasis_expression/Gaussian-latent-effects_expression_B.csv',stringsAsFactors=F)
betas_expr_latent_Gaussian <- merge(betas_expr_latent_Gaussian_A,betas_expr_latent_Gaussian_B,by=c("site","mutation","wildtype","mutant"),all=T,sort=T,suffixes=c("_A","_B"));rm(betas_expr_latent_Gaussian_A);rm(betas_expr_latent_Gaussian_B)
plot(betas_expr_latent_Gaussian$effect_A[!betas_expr_latent_Gaussian$site %in% seq(126,140)],betas_expr_latent_Gaussian$effect_B[!betas_expr_latent_Gaussian$site %in% seq(126,140)],pch=16,col="#00000067");summary(lm(betas_expr_latent_Gaussian$effect_B[!betas_expr_latent_Gaussian$site %in% seq(126,140)]~betas_expr_latent_Gaussian$effect_A[!betas_expr_latent_Gaussian$site %in% seq(126,140)]));abline(lm(betas_expr_latent_Gaussian$effect_B[!betas_expr_latent_Gaussian$site %in% seq(126,140)]~betas_expr_latent_Gaussian$effect_A[!betas_expr_latent_Gaussian$site %in% seq(126,140)]),col="red")
#excluding STOP variants?
plot(betas_expr_latent_Gaussian$effect_A[!betas_expr_latent_Gaussian$site %in% seq(126,140) & betas_expr_latent_Gaussian$mutant!="*"],betas_expr_latent_Gaussian$effect_B[!betas_expr_latent_Gaussian$site %in% seq(126,140) & betas_expr_latent_Gaussian$mutant!="*"],pch=16,col="#00000067",xlab="mutant deltaE, replicate A",ylab="mutant deltaE, replicate B");summary(lm(betas_expr_latent_Gaussian$effect_B[!betas_expr_latent_Gaussian$site %in% seq(126,140) & betas_expr_latent_Gaussian$mutant!="*"]~betas_expr_latent_Gaussian$effect_A[!betas_expr_latent_Gaussian$site %in% seq(126,140) & betas_expr_latent_Gaussian$mutant!="*"]));abline(lm(betas_expr_latent_Gaussian$effect_B[!betas_expr_latent_Gaussian$site %in% seq(126,140) & betas_expr_latent_Gaussian$mutant!="*"]~betas_expr_latent_Gaussian$effect_A[!betas_expr_latent_Gaussian$site %in% seq(126,140) & betas_expr_latent_Gaussian$mutant!="*"]),col="red")
#R-squared 0.845, 0.7146 w/o STOP (something nonlinear going on)

#Cauchy seems better, once again

#keep Cauchy betas data.frames, combine latent and observed measurements (though I want to use latent scale); remove others
rm(betas_bind_latent_Gaussian);rm(betas_bind_observed_Gaussian);rm(betas_expr_latent_Gaussian);rm(betas_expr_observed_Gaussian);rm(betas_bind_observed_nonepistatic_Cauchy);rm(betas_expr_observed_nonepistatic_Cauchy)
names(betas_bind_latent_Cauchy)[names(betas_bind_latent_Cauchy)=="effect_A"] <- "bind_latent_A"
names(betas_bind_latent_Cauchy)[names(betas_bind_latent_Cauchy)=="effect_B"] <- "bind_latent_B"
names(betas_bind_latent_Cauchy)[names(betas_bind_latent_Cauchy)=="effect_A"] <- "bind_latent_A"
names(betas_bind_latent_Cauchy)[names(betas_bind_latent_Cauchy)=="effect_B"] <- "bind_latent_B"
names(betas_bind_observed_Cauchy)[names(betas_bind_observed_Cauchy)=="effect_A"] <- "bind_observed_A"
names(betas_bind_observed_Cauchy)[names(betas_bind_observed_Cauchy)=="effect_B"] <- "bind_observed_B"
names(betas_bind_observed_Cauchy)[names(betas_bind_observed_Cauchy)=="effect_A"] <- "bind_observed_A"
names(betas_bind_observed_Cauchy)[names(betas_bind_observed_Cauchy)=="effect_B"] <- "bind_observed_B"
names(betas_expr_latent_Cauchy)[names(betas_expr_latent_Cauchy)=="effect_A"] <- "expr_latent_A"
names(betas_expr_latent_Cauchy)[names(betas_expr_latent_Cauchy)=="effect_B"] <- "expr_latent_B"
names(betas_expr_latent_Cauchy)[names(betas_expr_latent_Cauchy)=="effect_A"] <- "expr_latent_A"
names(betas_expr_latent_Cauchy)[names(betas_expr_latent_Cauchy)=="effect_B"] <- "expr_latent_B"
names(betas_expr_observed_Cauchy)[names(betas_expr_observed_Cauchy)=="effect_A"] <- "expr_observed_A"
names(betas_expr_observed_Cauchy)[names(betas_expr_observed_Cauchy)=="effect_B"] <- "expr_observed_B"
names(betas_expr_observed_Cauchy)[names(betas_expr_observed_Cauchy)=="effect_A"] <- "expr_observed_A"
names(betas_expr_observed_Cauchy)[names(betas_expr_observed_Cauchy)=="effect_B"] <- "expr_observed_B"

betas_bind <- merge(betas_bind_latent_Cauchy,betas_bind_observed_Cauchy,by=c("mutation","wildtype","site","mutant"),all=T,sort=F)
betas_expr <- merge(betas_expr_latent_Cauchy,betas_expr_observed_Cauchy,by=c("mutation","wildtype","site","mutant"),all=T,sort=F)

betas <- merge(betas_bind,betas_expr, by=c("site","mutation","wildtype","mutant"),sort=T,all=T)
betas <- data.table(betas)

rm(betas_bind_latent_Cauchy);rm(betas_bind_observed_Cauchy);rm(betas_expr_latent_Cauchy);rm(betas_expr_observed_Cauchy)

#how would replicates correlate if I just took the mutational effects from single mutant barcodes? also calculate number of times a beta was sampled across single-mut barcodes or all barcodes
#for each barcode that is determined in a replicate, output separated mutations
bc_bind[,aa_subs_list := list(strsplit(aa_substitutions,split=" ")),by=.(library,barcode)]
bc_expr[,aa_subs_list := list(strsplit(aa_substitutions,split=" ")),by=.(library,barcode)]
betas[,n_bc_bind_A := sum(unlist(lapply(bc_bind[library=="libA2" & !is.na(ddG),aa_subs_list], function(x) mutation %in% x))),by=mutation]
betas[,n_bc_bind_B := sum(unlist(lapply(bc_bind[library=="libB2" & !is.na(ddG),aa_subs_list], function(x) mutation %in% x))),by=mutation]
betas[,n_bc_expr_A := sum(unlist(lapply(bc_expr[library=="libA2" & !is.na(ML_meanF),aa_subs_list], function(x) mutation %in% x))),by=mutation]
betas[,n_bc_expr_B := sum(unlist(lapply(bc_expr[library=="libB2" & !is.na(ML_meanF),aa_subs_list], function(x) mutation %in% x))),by=mutation]
#last calculation takes a while, save betas table for interim work until finalizing notebook
#save(betas,file="results/analyze_mut_effects/betas.temp.Rda")
load("results/analyze_mut_effects/betas.temp.Rda")
for(i in 1:nrow(betas)){
  print(i)
  ddGs_A <- bc_bind[aa_substitutions==betas[i,"mutation"] & library=="libA2",ddG]
  ddGs_B <- bc_bind[aa_substitutions==betas[i,"mutation"] & library=="libB2",ddG]
  meanFs_A <- bc_expr[aa_substitutions==betas[i,"mutation"] & library=="libA2",ML_meanF]
  meanFs_B <- bc_expr[aa_substitutions==betas[i,"mutation"] & library=="libB2",ML_meanF]
  betas$n_bc_1mut_bind_A[i] <- sum(!is.na(ddGs_A))
  betas$n_bc_1mut_bind_B[i] <- sum(!is.na(ddGs_B))
  betas$n_bc_1mut_expr_A[i] <- sum(!is.na(meanFs_A))
  betas$n_bc_1mut_expr_B[i] <- sum(!is.na(meanFs_B))
  betas$bc_1mut_bind_A[i] <- mean(ddGs_A,na.rm=T)
  betas$bc_1mut_bind_B[i] <- mean(ddGs_B,na.rm=T)
  betas$bc_1mut_expr_A[i] <- mean(meanFs_A,na.rm=T)
  betas$bc_1mut_expr_B[i] <- mean(meanFs_B,na.rm=T)
}

sum(!is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_bind_A"]))/nrow(betas[betas$mutant!="*" & !betas$site %in% seq(126,140) & wildtype!=mutant ,])
sum(!is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_bind_B"]))/nrow(betas[betas$mutant!="*" & !betas$site %in% seq(126,140) & wildtype!=mutant,])
sum(!is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_bind_A"]) | !is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_bind_B"]))/nrow(betas[betas$mutant!="*" & !betas$site %in% seq(126,140) & wildtype!=mutant,])
sum(!is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_bind_A"]) & !is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_bind_B"]))/nrow(betas[betas$mutant!="*" & !betas$site %in% seq(126,140) & wildtype!=mutant,])
#97.0% single mutants directly measured in titration repA, 94.2% in repB, 99.05% in at least one rep, 92.2% in both
sum(!is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_expr_A"]))/nrow(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),])
sum(!is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_expr_B"]))/nrow(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),])
sum(!is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_expr_A"]) | !is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_expr_B"]))/nrow(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),])
sum(!is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_expr_A"]) & !is.na(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),"bc_1mut_expr_B"]))/nrow(betas[betas$mutant!="*" & wildtype!=mutant & !betas$site %in% seq(126,140),])
#96.5% single mutants directly measured in expression repA, 97.9% in repB, 99.3% in at least one rep, 95.1% in both

#pdf(file="results/analyze_mut_effects/correlations_direct-single-mutant-bc-mean-ddGs.pdf",bg="transparent",width=4,height=4.5)
plot(betas[mutant!="*" & !site %in% seq(126,140),bc_1mut_bind_A],betas[mutant!="*" & !site %in% seq(126,140),bc_1mut_bind_B],pch=16,col="#00000067",xlab="single mutant bcs ddG, rep A",ylab="single mutant bcs ddG, rep B",main="R-squared 0.9490 (92.2% of mutants)");model <- lm(betas[mutant!="*" & !site %in% seq(126,140),bc_1mut_bind_B] ~ betas[mutant!="*" & !site %in% seq(126,140),bc_1mut_bind_A]);summary(model);abline(model,col="red")
#R-squared 0.9490 for direct measurements of single mutants
#actually marginally better than the global epistasis betas correlation? (though, I should check correlation when excluding same NA betas as are missing here:)
plot(betas[mutant!="*" & !site %in% seq(126,140) & !is.na(bc_1mut_bind_A) & !is.na(bc_1mut_bind_B),bind_latent_A],betas[mutant!="*" & !site %in% seq(126,140) & !is.na(bc_1mut_bind_A) & !is.na(bc_1mut_bind_B),bind_latent_B],pch=16,col="#00000067",xlab="single mutant bcs betas bind (latent scale), rep A",ylab="single mutant bcs betas bind (latent scale), rep B",main="R-squared 0.9511");model <- lm(betas[mutant!="*" & !site %in% seq(126,140) & !is.na(bc_1mut_bind_A) & !is.na(bc_1mut_bind_B),bind_latent_B] ~ betas[mutant!="*" & !site %in% seq(126,140) & !is.na(bc_1mut_bind_A) & !is.na(bc_1mut_bind_B),bind_latent_A]);summary(model);abline(model,col="red")
plot(betas[mutant!="*" & !site %in% seq(126,140) & !is.na(bc_1mut_bind_A) & !is.na(bc_1mut_bind_B),bind_observed_A],betas[mutant!="*" & !site %in% seq(126,140) & !is.na(bc_1mut_bind_A) & !is.na(bc_1mut_bind_B),bind_observed_B],pch=16,col="#00000067",xlab="single mutant bcs betas bind (observed scale), rep A",ylab="single mutant bcs betas bind (observed scale), rep B",main="R-squared 0.9501");model <- lm(betas[mutant!="*" & !site %in% seq(126,140) & !is.na(bc_1mut_bind_A) & !is.na(bc_1mut_bind_B),bind_observed_B] ~ betas[mutant!="*" & !site %in% seq(126,140) & !is.na(bc_1mut_bind_A) & !is.na(bc_1mut_bind_B),bind_observed_A]);summary(model);abline(model,col="red")

#calculate means from duplicate global epistasis coefficients
betas$bind_latent <- rowMeans(betas[,c("bind_latent_A","bind_latent_B")],na.rm=T)
betas$bind_observed <- rowMeans(betas[,c("bind_observed_A","bind_observed_B")],na.rm=T)
betas$expr_latent <- rowMeans(betas[,c("expr_latent_A","expr_latent_B")],na.rm=T)
betas$expr_observed <- rowMeans(betas[,c("expr_observed_A","expr_observed_B")],na.rm=T)

#plot final betas versus the average of single-mutant bc values
plot(betas[mutant!="*" & !site %in% seq(126,140),(bc_1mut_bind_A+bc_1mut_bind_B)/2],betas[mutant!="*" & !site %in% seq(126,140),bind_observed],pch=16,col="#00000067",xlab="ddGmut from single mutant barcodes",ylab="ddGmut from global epistasis model",main="R-squared 0.9831")
abline(a=0,b=1,lty=2,col="red")
model <- lm(betas[mutant!="*" & !site %in% seq(126,140),bind_observed] ~ betas[mutant!="*" & !site %in% seq(126,140),(bc_1mut_bind_A+bc_1mut_bind_B)/2]);summary(model)
#dev.off()

#coverage among determined mutants (excluding linker sites):
table(betas_expr[!(betas_expr$site %in% seq(126,140)),"site"])[table(betas_expr[!(betas_expr$site %in% seq(126,140)),"site"])<21]
#no missing mutations for expression, between the two replicates
table(betas_bind[!(betas_bind$site %in% seq(126,140)),"site"])[table(betas_bind[!(betas_bind$site %in% seq(126,140)),"site"])<20]
#for binding, missing one mutation each at sites 22 (K), 99 (Q), 144 (E)
#number of backgrounds on which a mutation is sampled, for binding and expression, pooling variants with measurements in libraries A and B
#pdf(file="results/analyze_mut_effects/mutation-coverage-ecdf.pdf",width=4.5,height=4,useDingbats=F,bg="transparent")
plot(ecdf(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_bind_A+n_bc_bind_B]),verticals=T,pch=NA,col.01line=NA,main=paste("pooled replicates, binding scores, median =",median(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_bind_A+n_bc_bind_B])),xlab="number of barcodes",ylab="fraction mutations found < X times")
abline(v=median(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_bind_A+n_bc_bind_B]),lty=2);abline(h=0.5,lty=2)
plot(ecdf(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_expr_A+n_bc_expr_B]),verticals=T,pch=NA,col.01line=NA,main=paste("pooled replicates, expression scores, median =",median(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_expr_A+n_bc_expr_B])),xlab="number of barcodes",ylab="fraction mutations found < X times")
abline(v=median(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_expr_A+n_bc_expr_B]),lty=2);abline(h=0.5,lty=2)
#number of single-mutant backgrounds on which a mutation is sampled
plot(ecdf(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_1mut_bind_A+n_bc_1mut_bind_B]),verticals=T,pch=NA,col.01line=NA,main=paste("pooled replicates 1mut-barcodes, binding scores, median =",median(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_1mut_bind_A+n_bc_1mut_bind_B])),xlab="number of single-mutant barcodes",ylab="fraction mutations found < X times",xlim=c(0,100))
abline(v=median(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_1mut_bind_A+n_bc_1mut_bind_B]),lty=2);abline(h=0.5,lty=2)
plot(ecdf(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_1mut_expr_A+n_bc_1mut_expr_B]),verticals=T,pch=NA,col.01line=NA,main=paste("pooled replicates 1mut-barcodes, expression scores, median =",median(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_1mut_expr_A+n_bc_1mut_expr_B])),xlab="number of single-mutant barcodes",ylab="fraction mutations found < X times")
abline(v=median(betas[!site %in% seq(126,140) & mutant!="*" & wildtype!=mutant,n_bc_1mut_expr_A+n_bc_1mut_expr_B]),lty=2);abline(h=0.5,lty=2)
#dev.off()

#as an example, plot of correlation of mutational ddGs with censoring for number of barcode backgrounds (like a coverage stat); as expected, correlation gets better when censoring for mutations found across more barcodes
min <- 10;plot(betas[!site %in% seq(126,140) & n_bc_bind_A > min & n_bc_bind_B > min,bind_observed_A],betas[!site %in% seq(126,140) & n_bc_bind_A > min & n_bc_bind_B > min,bind_observed_B],pch=16);model<-lm(betas[!site %in% seq(126,140) & n_bc_bind_A > min & n_bc_bind_B > min,bind_observed_B]~betas[!site %in% seq(126,140) & n_bc_bind_A > min & n_bc_bind_B > min,bind_observed_A]);abline(model);summary(model);nrow(betas[!site %in% seq(126,140) & n_bc_bind_A > min & n_bc_bind_B > min & !is.na(bind_observed_A) & !is.na(bind_observed_B),])


#correlation between mutational effects on stability, binding
plot(betas[mutant!="*",expr_latent],betas[mutant!="*",bind_latent],pch=16,col="#00000033")
plot(betas[mutant!="*",expr_observed],betas[mutant!="*",bind_observed],pch=16,col="#00000067")
plot(betas[mutant!="*",expr_latent],betas[mutant!="*",bind_observed],pch=16,col="#00000067")


nsites <- 244
sites <- data.frame(site=1:nsites)
for(i in 1:nrow(sites)){
  sites$AA[i] <- sites_alignment[sites_alignment$site_scFv==sites$site[i], "amino_acid_scFv"]
  sites$kabat_site[i] <- sites_alignment[sites_alignment$site_scFv==sites$site[i], "site_Kabat"]
  sites$kabat_chain[i] <- sites_alignment[sites_alignment$site_scFv==sites$site[i], "chain_Kabat"]
  sites$kabat_annotation[i] <- sites_alignment[sites_alignment$site_scFv==sites$site[i], "Kabat_annotation"]
  sites$pdb_site[i] <- sites_alignment[sites_alignment$site_scFv==sites$site[i], "site_3U7Y"]
  sites$pdb_chain[i] <- sites_alignment[sites_alignment$site_scFv==sites$site[i], "chain_3U7Y"]
  sites$mean_bind[i] <- mean(betas[site==sites$site[i],bind_latent],na.rm=TRUE)
  sites$mean_bind_A[i] <- mean(betas[site==sites$site[i],bind_latent_A],na.rm=TRUE)
  sites$mean_bind_B[i] <- mean(betas[site==sites$site[i],bind_latent_B],na.rm=TRUE)
  sites$mean_bind_observed[i] <- mean(betas[site==sites$site[i],bind_observed],na.rm=TRUE)
  sites$mean_bind_observed_A[i] <- mean(betas[site==sites$site[i],bind_observed_A],na.rm=TRUE)
  sites$mean_bind_observed_B[i] <- mean(betas[site==sites$site[i],bind_observed_B],na.rm=TRUE)
  sites$max_bind[i] <- max(betas[site==sites$site[i] & wildtype != mutant,bind_observed],na.rm=TRUE)
  sites$min_bind[i] <- min(betas[site==sites$site[i] & wildtype != mutant,bind_observed],na.rm=TRUE)
  sites$mean_expr[i] <- mean(betas[site==sites$site[i] & mutant!="*",expr_latent],na.rm=TRUE)
  sites$mean_expr_A[i] <- mean(betas[site==sites$site[i] & mutant!="*",expr_latent_A],na.rm=TRUE)
  sites$mean_expr_B[i] <- mean(betas[site==sites$site[i] & mutant!="*",expr_latent_B],na.rm=TRUE)
  sites$mean_expr_observed[i] <- mean(betas[site==sites$site[i] & mutant!="*",expr_observed],na.rm=TRUE)
  sites$mean_expr_observed_A[i] <- mean(betas[site==sites$site[i] & mutant!="*",expr_observed_A],na.rm=TRUE)
  sites$mean_expr_observed_B[i] <- mean(betas[site==sites$site[i] & mutant!="*",expr_observed_B],na.rm=TRUE)
  sites$max_expr[i] <- max(betas[site==sites$site[i] & wildtype != mutant & mutant!="*",expr_latent],na.rm=TRUE)
  sites$min_expr[i] <- min(betas[site==sites$site[i] & wildtype != mutant & mutant!="*",expr_latent],na.rm=TRUE)
}

#add pdb, kabat site numbers to betas dataframe
for(i in 1:nrow(betas)){
  betas$kabat_site[i]<- sites[sites$site==betas$site[i],"kabat_site"]
  betas$kabat_chain[i]<- sites[sites$site==betas$site[i],"kabat_chain"]
  betas$kabat_annotation[i]<- sites[sites$site==betas$site[i],"kabat_annotation"]
  betas$pdb_site[i] <- sites[sites$site==betas$site[i],"pdb_site"]
  betas$pdb_chain[i] <- sites[sites$site==betas$site[i],"pdb_chain"]
}


#to complement correlation in mutational effects above, what is the correlation in mean mutational effect per site?
#pdf(file="results/analyze_mut_effects/correlation-mean-beta-per-site_bind-latent.pdf",bg="transparent",width=4,height=4.5,useDingbats=F)
plot(sites$mean_bind_A[!sites$site %in% seq(126,140)],sites$mean_bind_B[!sites$site %in% seq(126,140)],pch=16,col="#00000067",xlab="rep A bind, mean beta per site (latent scale)",ylab="rep B bind, mean beta per site (latent scale)",main="R-squared 0.9904");model <- lm(sites$mean_bind_B[!sites$site %in% seq(126,140)]~sites$mean_bind_A[!sites$site %in% seq(126,140)]);summary(model);abline(model,col="red")
#dev.off()
#R-squared 0.9901 (latent effects on binding)
#pdf(file="results/analyze_mut_effects/correlation-mean-beta-per-site_bind-observed.pdf",bg="transparent",width=4,height=4.5,useDingbats=F)
plot(sites$mean_bind_observed_A[!sites$site %in% seq(126,140)],sites$mean_bind_observed_B[!sites$site %in% seq(126,140)],pch=16,col="#00000067",xlab="rep A bind, mean beta per site (observed scale)",ylab="rep B bind, mean beta per site (observed scale)",main="R-squared 0.9865");model <- lm(sites$mean_bind_observed_B[!sites$site %in% seq(126,140)]~sites$mean_bind_observed_A[!sites$site %in% seq(126,140)]);summary(model);abline(model,col="red")
#dev.off()
#R-squared 0.9865 (observed scale effects on binding)
#pdf(file="results/analyze_mut_effects/correlation-mean-beta-per-site_expr-latent.pdf",bg="transparent",width=4,height=4.5)
plot(sites$mean_expr_A[!sites$site %in% seq(126,140)],sites$mean_expr_B[!sites$site %in% seq(126,140)],pch=16,col="#00000067",xlab="rep A expr, mean beta per site (latent scale)",ylab="rep B expr, mean beta per site (latent scale)",main="R-squared 0.9736");model <- lm(sites$mean_expr_B[!sites$site %in% seq(126,140)]~sites$mean_expr_A[!sites$site %in% seq(126,140)]);summary(model);abline(model,col="red")
#dev.off()
#R-squared 0.9736 (latent effects on expression)
#pdf(file="results/analyze_mut_effects/correlation-mean-beta-per-site_expr-observed.pdf",bg="transparent",width=4,height=4.5)
plot(sites$mean_expr_observed_A[!sites$site %in% seq(126,140)],sites$mean_expr_observed_B[!sites$site %in% seq(126,140)],pch=16,col="#00000067",xlab="rep A expr, mean beta per site (observed scale)",ylab="rep B expr, mean beta per site (observed scale)",main="R-squared 0.9830");model <- lm(sites$mean_expr_observed_B[!sites$site %in% seq(126,140)]~sites$mean_expr_observed_A[!sites$site %in% seq(126,140)]);summary(model);abline(model,col="red")
#dev.off()
#R-squared 0.9830 (observfed scale effects on expression)

plot(sites$mean_bind_observed,sites$mean_expr,pch=16,col="#00000050")


#histogram of ddGs
hist(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140),bind_observed],col="gray50",breaks=30,xlim=c(-1,8),xlab="ddGmut (kcal/mol)",main="")
#want to have an overlay of the expected noise in the sampling distribution of mutational effects -- do a pooled SEM from all duplicate measurements (near zero, since larger error for larger coefficients)
#pooled standard deviation SEM = sqrt[sum({s_i - mean(s_i)}^2)/N]
pooled_SEM_bind <- sqrt(sum(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & !is.na(bind_observed_A) & !is.na(bind_observed_B) & bind_observed < 1,bind_observed_A-bind_observed]^2,betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & !is.na(bind_observed_A) & !is.na(bind_observed_B) & bind_observed < 1,bind_observed_B-bind_observed]^2)/(nrow(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & !is.na(bind_observed_A) & !is.na(bind_observed_B) & bind_observed < 1,])*2))
sqrt(sum(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & !is.na(bind_observed_A) & !is.na(bind_observed_B) & bind_observed < 1,bind_observed_A-bind_observed_B]^2)/(2*nrow(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & !is.na(bind_observed_A) & !is.na(bind_observed_B) & bind_observed < 1,])))/sqrt(2)
#can use this pooled_SEM_bind to simulate a neutral sampling distribution on any histogram I'm wanting to show

#want to look at 'beneficial' mutations -- ascribe a p-value of being beneficial by comparing to the error distribution, then use FDR to select candidates
betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140),p_ben := pnorm(bind_observed,mean=0,sd=pooled_SEM_bind)]
betas$FDR_ben <- p.adjust(betas$p_ben, method="fdr")

#pdf(file="results/analyze_mut_effects/hist_betas-bind-observed.pdf",width=5,height=4,useDingbats=F)
hist(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140),bind_observed],col="gray50",breaks=30,xlim=c(-1,8),xlab="ddGmut (kcal/mol)",main="")
#curve(100*dnorm(x,mean=0,sd=pooled_SEM_bind),add=T,col="darkblue",lwd=2)
#instead of curve, simulate draws from the distribution and plot as histogram
set.seed(1981);hist(rnorm(1000,mean=0,sd=pooled_SEM_bind),add=T,col="#0000ff50",breaks=4)
#abline(v=max(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & FDR_ben<0.1,bind_observed]),lty=2)
#dev.off()

#expression
#pdf(file="results/analyze_mut_effects/hist_betas-expr-latent.pdf",width=5,height=4,useDingbats=F,bg="transparent")
hist(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140),expr_latent],col="gray50",breaks=20,xlab="deltaE (a.u.)",main="",xlim=c(-14,1))
hist(betas[wildtype != mutant & mutant == "*" & !site %in% seq(126,140),expr_latent],col="red",breaks=40,xlab="deltaE (a.ul)",main="",add=T)
#dev.off()

#look at the most beneficial mutations for binding -- get sites in pdb numbering to visualize on structure. also look to see if it's possibly an expression thing
hist(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & bind_observed<1,bind_observed],col="gray50",xlab="ddGmut (kcal/mol)",main="")
#curve(100*dnorm(x,mean=0,sd=pooled_SEM_bind),add=T,col="darkblue",lwd=2)
set.seed(1981);hist(rnorm(500,mean=0,sd=pooled_SEM_bind),add=T,col="#0000ff50",breaks=4)
abline(v=betas[betas$mutation=="G55W","bind_observed"],col="green",lty=2) #previously described beneficial mutation
abline(v=max(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & FDR_ben<0.1,bind_observed]),lty=2) #10% FDR line

beneficial_muts <- betas[bind_observed <= max(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & FDR_ben<0.1,bind_observed]) & !is.na(betas$bind_observed) & !(betas$site %in% seq(126,140)),]
beneficial_muts[,c("mutation","kabat_site","kabat_chain","pdb_site","bind_observed","expr_latent","bc_1mut_bind_A","bc_1mut_bind_B","FDR_ben")]

#pdf(file="results/analyze_mut_effects/scatter_zoom-on-beneficial_expr-latent_v_ddG-observed_colored-validation-beneficial-muts.pdf",width=4,height=4.5,bg="transparent",useDingbats=F)
plot(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140),bind_observed],betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140),expr_latent],pch=16,col="#00000066",xlim=c(1,-1),ylim=c(-1,1),xlab="ddG_mut (kcal/mol)",ylab="delta-expression")
points(beneficial_muts$bind_observed,beneficial_muts$expr_latent,pch=16,col="blue")
abline(h=0,lty=2);abline(v=0,lty=2)
points(beneficial_muts[mutation %in% c("E28P","L51I","G55N","G55W","D77G","D77N","S168P","S170P"),bind_observed],beneficial_muts[mutation %in% c("E28P","L51I","G55N","G55W","D77G","D77N","S168P","S170P"),expr_latent],pch=16,col="green")#points for mutants I want to validate
#dev.off()

#are these beneficial mutations unlikely to occur with the SHM machinery?
#read in mutabilities from SHM simulations
aa_mutabilities <- read.csv(file="results/SHM-propensities/mutation-probabilites.csv",stringsAsFactors = F)
for(i in 1:nrow(betas)){
  if(!is.na(betas$kabat_chain[i])){
    betas$mature_mutability[i] <- aa_mutabilities[aa_mutabilities$mutation==betas[i,mutation],"mature_mutability"]    
  }else{
    betas$mature_mutability[i] <- NA  
  }
}

#interesting, is there a trend here???
plot(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140), mature_mutability], betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140), bind_observed],pch=16,col="#00000015",log="x")

hist(log10(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & FDR_ben>0.1,mature_mutability]+1e-7),col="gray50")
hist(log10(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & FDR_ben<0.1,mature_mutability]+1e-7),col="green",add=T)

hist(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & FDR_ben>0.1,mature_mutability],col="gray50")
hist(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & FDR_ben<0.1,mature_mutability],col="green",add=F)

betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & FDR_ben<0.1,.(mutation,mature_mutability)]

#the ones I'm validating
betas[mutation %in% c("E28P","L51I","G55N","G55W","D77G","D77N","S168P","S170P"),.(mutation,mature_mutability)]

wilcox.test(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & FDR_ben>0.1,mature_mutability],betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & FDR_ben<0.1,mature_mutability])
wilcox.test(betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & mutation %in% c("E28P","L51I","G55N","G55W","D77G","D77N","S168P","S170P"),mature_mutability],betas[wildtype != mutant & mutant != "*" & !site %in% seq(126,140) & !(mutation %in% c("E28P","L51I","G55N","G55W","D77G","D77N","S168P","S170P")),mature_mutability])
#not significantly less likely to occur with SHM machinery than the typical mutation, though I might note that they are mostly <0.01 probability in a 3% SHM bout so they're definitely not mutationally privileged

#what about comparing on the sites scale?
site_mutabilities <- read.csv(file="results/SHM-propensities/sites_alignment.csv")
sites$naive_mutability <- NA
sites$mature_mutability <- NA
for(i in 1:nrow(sites)){
  sites$naive_mutability[i] <- site_mutabilities[site_mutabilities$site_scFv==sites[i,"site"] & !is.na(site_mutabilities$site_scFv),"mutability_naive"]
  sites$mature_mutability[i] <- site_mutabilities[site_mutabilities$site_scFv==sites[i,"site"] & !is.na(site_mutabilities$site_scFv),"mutability_mature"]
}

plot(sites[!(sites$site %in% seq(126,140)),"mature_mutability"],sites[!(sites$site %in% seq(126,140)),"mean_bind_observed"],pch=16)
plot(sites[!(sites$site %in% seq(126,140)),"naive_mutability"],sites[!(sites$site %in% seq(126,140)),"mean_bind_observed"],pch=16)

wilcox.test(sites[sites$site %in% beneficial_muts$site & !is.na(sites$site) & !(sites$site %in% seq(126,140)),"mature_mutability"],sites[!(sites$site %in% beneficial_muts$site) & !is.na(sites$site) & !(sites$site %in% seq(126,140)),"mature_mutability"])
#no difference in average mutabilities

#biochemical correlates of mutational effects
#distance to interface (change in RSA when liganded)
#RSA and relative amount of RSA buried at gp120 interface
#use bio3d package, installed dssp executable to ./src with wget ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2.0.4-linux-amd64 -O ./src/dssp and chmod a+x ./src/dssp
dssp <- dssp(pdb, exefile="~/bloom_j/computational_notebooks/tstarr/2019/NIH45-46_DMS/src/dssp")
dssp_fab <- dssp(pdb_fab, exefile="~/bloom_j/computational_notebooks/tstarr/2019/NIH45-46_DMS/src/dssp")

#annotate sites df with secondary structural element (sse), SASA_bound, and SASA_unbound
resi <- pdb_atoms[pdb_atoms$elety=="CA" & pdb_atoms$resid != "FLC",]
for(i in 1:nrow(resi)){
  resi$SSE[i] <- dssp$sse[i]
  resi$SASA_bound[i] <- dssp$acc[i]
}
resi_fab <- pdb_fab_atoms[pdb_fab_atoms$elety=="CA" & pdb_fab_atoms$resid != "FLC",]
for(i in 1:nrow(resi_fab)){
  resi_fab$SASA_unbound[i] <- dssp_fab$acc[i]
}
#calculate RSA by normalizing SASA by maximum SASA of different amino acid types. Max from Tien et al.
max_SASA <- data.frame(AA=c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"),max_SASA=c(121,265,187,187,148,214,214,97,216,195,191,230,203,228,154,143,163,264,255,165))

for(i in 1:nrow(sites)){
  if(!is.na(sites$pdb_site[i])){
    sites$SSE[i] <- resi[resi$resno==sites$pdb_site[i] & resi$chain==sites$pdb_chain[i],"SSE"]
    sites$SASA_bound[i] <- resi[resi$resno==sites$pdb_site[i] & resi$chain==sites$pdb_chain[i],"SASA_bound"]
    sites$SASA_unbound[i] <- resi_fab[resi_fab$resno==sites$pdb_site[i] & resi_fab$chain==sites$pdb_chain[i],"SASA_unbound"]
    sites$RSA_bound[i] <- sites$SASA_bound[i]/max_SASA[max_SASA$AA==sites$AA[i],"max_SASA"]
    sites$RSA_unbound[i] <- sites$SASA_unbound[i]/max_SASA[max_SASA$AA==sites$AA[i],"max_SASA"]
    sites$buried_SASA[i] <- sites$SASA_unbound[i] - sites$SASA_bound[i]
    sites$buried_RSA[i] <- sites$RSA_unbound[i] - sites$RSA_bound[i]
  }else{
    sites$SSE[i] <- NA
    sites$SASA_bound[i] <- NA
    sites$SASA_unbound[i] <- NA
    sites$RSA_bound[i] <- NA
    sites$RSA_unbound[i] <- NA
    sites$buried_SASA[i] <- NA
    sites$buried_RSA[i] <- NA
  }
}

plot(sites$SASA_unbound,sites$SASA_bound,pch=16,col="#00000050")
plot(sites$RSA_unbound,sites$RSA_bound,pch=16,col="#00000050")

#plot bound and unbound as points on same plot?
#pdf(file="results/analyze_mut_effects/scatter_mean-ddG_v_RSA_bound-and-unbound_latent-scale.pdf",bg="transparent",width=4,height=4.5,useDingbats=F)
plot(sites$RSA_unbound,sites$mean_bind,pch=16,col="gray85",ylim=rev(range(sites$mean_bind,na.rm=T)), ylab="mean ddG of mutations at site (kcal/mol)",xlab="relative solvent accessibility")
points(sites$RSA_bound,sites$mean_bind,pch=16)
#dev.off()
#observed scale
#pdf(file="results/analyze_mut_effects/scatter_mean-ddG_v_RSA_bound-and-unbound_observed-scale.pdf",bg="transparent",width=4,height=4.5,useDingbats=F)
plot(sites$RSA_unbound,sites$mean_bind_observed,pch=16,col="gray85",ylim=rev(range(sites$mean_bind_observed,na.rm=T)), ylab="mean ddG of mutations at site (kcal/mol)",xlab="relative solvent accessibility")
points(sites$RSA_bound,sites$mean_bind_observed,pch=16)
#dev.off()
plot(sites$SASA_bound,sites$mean_expr,pch=16,col="#00000050")
plot(sites$SASA_unbound,sites$mean_expr,pch=16,col="#00000050")

#plot mutational effects to b factors on structures, visualize
#output pdbs with replaced b factors per http://thegrantlab.org/bio3d/tutorials/structure-analysis)
pdb_mean_bind <- pdb
pdb_mean_bind$atom$b <- NA
for(i in 1:nrow(pdb_mean_bind$atom)){
  if(!is.na(pdb_mean_bind$atom$insert[i])){res <- paste(pdb_mean_bind$atom$resno[i],pdb_mean_bind$atom$insert[i],sep="")}else{res <- as.character(pdb_mean_bind$atom$resno[i])}
  chain <- pdb_mean_bind$atom$chain[i]
  mean_bind <- sites[sites$pdb_site==res & sites$pdb_chain==chain & !is.na(sites$pdb_site),"mean_bind_observed"]
  if(length(mean_bind)>0){pdb_mean_bind$atom$b[i] <- mean_bind}else{pdb_mean_bind$atom$b[i] <- 0}
}
#save pdb
write.pdb(pdb=pdb_mean_bind,file=paste(config$analyze_mut_effects_dir,"/3u7y_b-factor-mean-ddG.pdb",sep=""), b = pdb_mean_bind$atom$b)

#want to correlate mutational constraint on the antibody position with the mutation constraint at contacting sites on the HIV gp120 protein
#use entropy from mutational preferences as determined by Haddox et al. PLOS Pathogens, who performed DMS on the lab-adapted clade B HIV Env from LAI (can also do with the transmitted founder stains from clade A in Haddox et al. eLife paper)
LAI_prefs <- read.csv(file="data/gp120_alignments/Haddox_LAI-preferences.csv",stringsAsFactors = F)

#take all pairs of contacting residues (heavy atoms <5A distnace), see how constraint on the HIV residue correlates with constraint on the Ab residue
calc.dist <- function(x1,y1,z1,x2,y2,z2){
  return(sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2))
}

nsites <- 244
pdb_dists <- expand.grid(site1=1:nsites,site2=1:nsites)
pdb_dists <- pdb_dists[order(pdb_dists$site1, pdb_dists$site2),]
pdb_dists <- pdb_dists[pdb_dists$site1 < pdb_dists$site2,]
for(i in 1:nrow(pdb_dists)){
  print(i)
  pdb_dists$site1_PDB[i] <- sites_alignment[sites_alignment$site_scFv==pdb_dists$site1[i],"site_3U7Y"]
  pdb_dists$site1_chain[i] <- sites_alignment[sites_alignment$site_scFv==pdb_dists$site1[i],"chain_3U7Y"]
  pdb_dists$site2_PDB[i] <- sites_alignment[sites_alignment$site_scFv==pdb_dists$site2[i],"site_3U7Y"]
  pdb_dists$site2_chain[i] <- sites_alignment[sites_alignment$site_scFv==pdb_dists$site2[i],"chain_3U7Y"]
}

for(i in 1:nrow(pdb_dists)){
  print(i)
  if(!is.na(pdb_dists[i,"site1_PDB"]) & !is.na(pdb_dists[i,"site2_PDB"])){
    atoms1 <- pdb_atoms[pdb_atoms$chain==pdb_dists[i,"site1_chain"] & pdb_atoms$resno==pdb_dists[i,"site1_PDB"],]
    atoms2 <- pdb_atoms[pdb_atoms$chain==pdb_dists[i,"site2_chain"] & pdb_atoms$resno==pdb_dists[i,"site2_PDB"],]
    if(nrow(atoms1)>0 & nrow(atoms2)>0){
      pdb_dists$CA_dist[i] <- calc.dist(atoms1[atoms1$elety=="CA","x"],atoms1[atoms1$elety=="CA","y"],atoms1[atoms1$elety=="CA","z"],
                                        atoms2[atoms2$elety=="CA","x"],atoms2[atoms2$elety=="CA","y"],atoms2[atoms2$elety=="CA","z"])
      all_dists <- c()
      for(i1 in 1:nrow(atoms1)){for(i2 in 1:nrow(atoms2)){
        all_dists <- c(all_dists,calc.dist(atoms1[i1,"x"],atoms1[i1,"y"],atoms1[i1,"z"],atoms2[i2,"x"],atoms2[i2,"y"],atoms2[i2,"z"]))}}
      pdb_dists$min_dist[i] <- min(all_dists)
    }else{pdb_dists$CA_dist[i] <- NA;pdb_dists$min_dist[i] <- NA}
  }else{pdb_dists$CA_dist[i] <- NA;pdb_dists$min_dist[i] <- NA}
}

contacts <- data.frame(scFv_site=NA,pdb_site=NA,pdb_chain=NA,scFv_AA=NA,gp120_site=NA,gp120_AA=NA)
gp120_atoms <- pdb_atoms[pdb_atoms$chain=="G" & pdb_atoms$type != "HETATM",]
cutoff <- 5 #minimum distance between heavy atoms on residues
for(i in 1:nrow(sites)){
  print(i)
  if(!is.na(sites$pdb_site[i])){
    scFv_atoms <- pdb_atoms[pdb_atoms$chain==sites$pdb_chain[i] & pdb_atoms$resno==sites$pdb_site[i],]
    for(j in 1:nrow(scFv_atoms)){
      all_dists <- sapply(1:nrow(gp120_atoms), function(x) calc.dist(scFv_atoms[j,"x"],scFv_atoms[j,"y"],scFv_atoms[j,"z"],gp120_atoms[x,"x"],gp120_atoms[x,"y"],gp120_atoms[x,"z"]))
      if(length(all_dists[all_dists<cutoff])>0){
        gp120_contacts <- gp120_atoms[which(all_dists<cutoff),"resno"]
        for(k in 1:length(gp120_contacts)){
          contacts <- rbind(contacts, unlist(c(sites[i,c("site","pdb_site","pdb_chain","AA")],gp120_contacts[k],LAI_prefs[LAI_prefs$site==gp120_contacts[k],"wildtype"])))
        }
      }
    }
  }
}
contacts <- contacts[-1,]
#annotate with functional constraint on gp120, Fab residue from DMS datasets
for(i in 1:nrow(contacts)){
  contacts$scFv_constraint[i] <- sites[sites$site==contacts[i,"scFv_site"],"mean_bind_observed"]
  contacts$gp120_constraint[i] <- LAI_prefs[LAI_prefs$site==contacts[i,"gp120_site"],"entropy"]
}

contacts_unique <- contacts[!duplicated(contacts[]),]

plot(contacts_unique$gp120_constraint,contacts_unique$scFv_constraint,pch=16,col="#00000067")
summary(lm(contacts_unique$scFv_constraint~contacts_unique$gp120_constraint)) #no significant trend

#notice interaction between W110 (average ddG>6) and areas near the N276 glycan (absent in this gp120) -- Haddox et al. noted that N glycans are not preserved in this DMS dataset despite being conserved among many Envs. How does this plot compare if I use entropy from frequencies in natural alignment instead of from DMS preferences?
Env_align <- read.fasta(file="data/gp120_alignments/Haddox_group-M-align_AA_PDB-numbering.fasta")
Env_entropy <- entropy(Env_align)$H#.10 #uses "10 letter" alphabet (collapses biochemically similar AAs) to compute entropy
for(i in 1:nrow(LAI_prefs)){
  print(i)
  if(!is.na(LAI_prefs$alignment_index[i])){
    LAI_prefs$alignment_entropy[i] <- Env_entropy[LAI_prefs[i,"alignment_index"]]
  }else{
    LAI_prefs$alignment_entropy[i] <- NA
  }
}
plot(LAI_prefs$alignment_entropy,LAI_prefs$entropy,pch=16)

for(i in 1:nrow(contacts)){
  contacts$gp120_conservation[i] <- LAI_prefs[LAI_prefs$site==contacts[i,"gp120_site"],"alignment_entropy"]
}

contacts_unique <- contacts[!duplicated(contacts[]),]

plot(contacts_unique$gp120_conservation,contacts_unique$scFv_constraint,pch=16,col="#00000067")
summary(lm(contacts_unique$scFv_constraint~contacts_unique$gp120_conservation)) #no significant trend, though I see now later down when I boil down to single interaction per residue -- that it comes down to the loop D interaction, where constrained scFv sites contact variable HIV sites, which track to escape mutations

#make PDB where gp120 b factor is alignment conservation
pdb_Ab_gp120_conservation <- pdb_Ab_gp120_constraints
for(i in 1:nrow(pdb_Ab_gp120_conservation$atom)){
  if(pdb_Ab_gp120_conservation$atom$chain[i]=="G"){
    if(!is.na(pdb_Ab_gp120_conservation$atom$insert[i])){res <- paste(pdb_Ab_gp120_conservation$atom$resno[i],pdb_Ab_gp120_conservation$atom$insert[i],sep="")}else{res <- as.character(pdb_Ab_gp120_conservation$atom$resno[i])}
    site_entropy <- LAI_prefs[LAI_prefs$site==res,"alignment_entropy"]
    if(length(site_entropy)>0){if(!is.na(site_entropy)){pdb_Ab_gp120_conservation$atom$b[i] <- site_entropy}else{pdb_Ab_gp120_conservation$atom$b[i] <- 0}}
  }
}
#save pdb
write.pdb(pdb=pdb_Ab_gp120_conservation,file=paste(config$analyze_mut_effects_dir,"/3u7y_b-factor-mean-ddG-and-LAI-alignment-entropy.pdb",sep=""), b = pdb_Ab_gp120_conservation$atom$b)

####instead of looking for overall trend, take just the residues that are burying RSA at the interface -- does immune focusing explain which of these "matter" for binding? (which are mutationally constrained)
buried <- sites[sites$buried_RSA>0.0 & !is.na(sites$buried_RSA),] #sites that bury any SASA at interface
#for each of these residues, what's its closest gp120 residue?
for(i in 1:nrow(buried)){
  print(i)
  scFv_atoms <- pdb_atoms[pdb_atoms$chain==buried$pdb_chain[i] & pdb_atoms$resno==buried$pdb_site[i],]
  min_dist <- 100
  min_res <- NA
  for(j in 1:nrow(scFv_atoms)){
    all_dists <- sapply(1:nrow(gp120_atoms), function(x) calc.dist(scFv_atoms[j,"x"],scFv_atoms[j,"y"],scFv_atoms[j,"z"],gp120_atoms[x,"x"],gp120_atoms[x,"y"],gp120_atoms[x,"z"]))
    if(min(all_dists) < min_dist){
      min_dist <- min(all_dists)
      min_res <- gp120_atoms[which(all_dists==min(all_dists)),"resno"]
    }
  }
  buried$closest_gp120_dist[i] <- min_dist
  buried$closest_gp120_resi[i] <- min_res
}
for(i in 1:nrow(buried)){
  buried$closest_gp120_constraint[i] <- LAI_prefs[LAI_prefs$site==buried$closest_gp120_resi[i],"entropy"]
  buried$closest_gp120_conservation[i] <- LAI_prefs[LAI_prefs$site==buried$closest_gp120_resi[i],"alignment_entropy"]
}
plot(buried$closest_gp120_dist,buried$mean_bind_observed,pch=16,ylab="mean ddG",xlab="distance to nearest gp120 site (Angstrom)",ylim=rev(range(buried$mean_bind_observed,na.rm=T)))
#pdf(file="results/analyze_mut_effects/mean-ddG_v_alignment-entropy-of-nearest-gp120-contact_latent-scale.pdf",bg="transparent",width=4,height=4.5,useDingbats=F)
plot(buried$closest_gp120_conservation,buried$mean_bind,pch=16,ylab="mean ddG of mutations to contact residue",xlab="variability (alignment entropy) of nearest gp120 site",ylim=rev(range(buried$mean_bind,na.rm=T)))
#dev.off()
#pdf(file="results/analyze_mut_effects/mean-ddG_v_alignment-entropy-of-nearest-gp120-contact_observed-scale.pdf",bg="transparent",width=4,height=4.5,useDingbats=F)
plot(buried$closest_gp120_conservation,buried$mean_bind_observed,pch=16,ylab="mean ddG of mutations to contact residue",xlab="variability (alignment entropy) of nearest gp120 site",ylim=rev(range(buried$mean_bind_observed,na.rm=T)))
#dev.off()
buried[buried$closest_gp120_conservation>0.5 & buried$mean_bind_observed > 2,]
#three outlier here: N100E_hc with 281, W100F_hc with site 279, and Y91_lc with site 278
plot(buried$closest_gp120_conservation,buried$min_bind,pch=16,ylab="min ddG of mutations to contact residue",xlab="variability (alignment entropy) of nearest gp120 site",ylim=rev(range(buried$min_bind,na.rm=T)))
plot(buried$closest_gp120_constraint,buried$mean_bind,pch=16)

#look at Adam Immunity paper, Lynch J Virol, Bar et al. NEJM 2016 (VRC01 infusion and rebound) for escape mutations -- how do they map to this info?
#Adam sees VRC01 MAP highlights regions in the D loop, including especially site 279 (strongest constrained contact with W100F)
#Lynch validates four escape mutations in the donor 45 viral evolution, three of the four are in the D loop and are these exact three sites that map to non-epitope-focused contacts (278, 279, 281)!
#Bar et al. NEJM, the largest changes in several patients are within this loop D -- notably site 278 (T->S change, maintaining glycan), 279 (D279N), 281 (V281A), 282 (K282R) ;NO changes within the CD4bs, which is where the epitope focused contacts are

#read in alignments of VRC01-class antibodies from within the 45-46 clonal family and from independent lineages.
#germline reconstruction of the NIH45-46 naive antibody, from Bonsignori et al. 2018. Look at sites and states of SHM
intra_naive_align <- read.fasta(file="data/Ab_alignments/germline_align_NIH45-46-numbering.fasta")
#in "sites" df, annotate with naive amino acid at each site, whehter it mutated (SHM), and whether it's mutated in minVRC01 (Jardine et al.)
for(i in 1:nrow(sites)){
  sites$naive_AA[i] <- intra_naive_align$ali["VRC01_UCA_Bonsignori_MK032222/1-241",i]
  sites$minVRC01_AA[i] <- intra_naive_align$ali["min_VRC01",i]
  if(sites$naive_AA[i]==sites$AA[i] & !is.na(sites$naive_AA[i])){
    sites$SHM[i] <- F
  }else if(sites$naive_AA[i]!=sites$AA[i] & !is.na(sites$naive_AA[i])){
    sites$SHM[i] <- T
  }else if(is.na(sites$naive_AA[i])){
    sites$SHM[i] <- NA
  }
  if(sites$naive_AA[i] != "-"){
    sites$naive_ddG[i] <- betas[site==sites$site[i] & mutant==sites$naive_AA[i],bind_latent]
    sites$naive_ddG_observed[i] <- betas[site==sites$site[i] & mutant==sites$naive_AA[i],bind_observed]
  }else{
    sites$naive_ddG[i] <- NA
    sites$naive_ddG_observed[i] <- NA
  }
  if(i %in% 126:140){
    sites$naive_AA[i] <- NA
    sites$SHM[i] <- NA
    sites$naive_ddG[i] <- NA
    sites$naive_ddG_observed[i] <- NA
  }
}
#pdf(file="results/analyze_mut_effects/violin-plots_mean-ddG-of-mutations_SHM-v-non-SHM-sites_latent-scale.pdf",bg="transparent",width=4,height=4.5,useDingbats=F)
ggplot(sites[!is.na(sites$SHM),], aes(x=SHM,y=mean_bind))+
  geom_violin(trim=F)+coord_flip()+
  geom_dotplot(binaxis='y',stackdir='center',dotsize=0.4)
  #geom_jitter(shape=16,position=position_jitter(0.2))
  #stat_summary(fun.y=median,geom="point",shape=23,size=2)
wilcox.test(sites[!is.na(sites$SHM) & sites$SHM==T,"mean_bind"],sites[!is.na(sites$SHM) & sites$SHM==F,"mean_bind"]) #SHMs accrued at more tolerant positions
#dev.off()

#pdf(file="results/analyze_mut_effects/violin-plots_mean-ddG-of-mutations_SHM-v-non-SHM-sites_observed-scale.pdf",bg="transparent",width=4,height=4.5,useDingbats=F)
ggplot(sites[!is.na(sites$SHM),], aes(x=SHM,y=mean_bind_observed))+
  geom_violin(trim=F)+coord_flip()+
  geom_dotplot(binaxis='y',stackdir='center',dotsize=0.4)
#geom_jitter(shape=16,position=position_jitter(0.2))
#stat_summary(fun.y=median,geom="point",shape=23,size=2)
wilcox.test(sites[!is.na(sites$SHM) & sites$SHM==T,"mean_bind_observed"],sites[!is.na(sites$SHM) & sites$SHM==F,"mean_bind_observed"]) #SHMs accrued at more tolerant positions
#dev.off()

#pdf(file="results/analyze_mut_effects/violin-plots_mean-naive-mutabilities_SHM-v-non-SHM-sites.pdf",bg="transparent",width=4,height=4.5,useDingbats=F)
ggplot(sites[!is.na(sites$SHM),], aes(x=SHM,y=naive_mutability))+
  geom_violin(trim=F)+coord_flip()+
  geom_dotplot(binaxis='y',stackdir='center',dotsize=0.4)
#geom_jitter(shape=16,position=position_jitter(0.2))
#stat_summary(fun.y=median,geom="point",shape=23,size=2)
wilcox.test(sites[!is.na(sites$SHM) & sites$SHM==T,"naive_mutability"],sites[!is.na(sites$SHM) & sites$SHM==F,"naive_mutability"]) #SHMs are not more inherently mutable sites?
#dev.off()

#pdf(file="results/analyze_mut_effects/violin-plots_mean-mature-mutabilities_SHM-v-non-SHM-sites.pdf",bg="transparent",width=4,height=4.5,useDingbats=F)
ggplot(sites[!is.na(sites$SHM),], aes(x=SHM,y=mature_mutability))+
  geom_violin(trim=F)+coord_flip()+
  geom_dotplot(binaxis='y',stackdir='center',dotsize=0.4)
#geom_jitter(shape=16,position=position_jitter(0.2))
#stat_summary(fun.y=median,geom="point",shape=23,size=2)
wilcox.test(sites[!is.na(sites$SHM) & sites$SHM==T,"mature_mutability"],sites[!is.na(sites$SHM) & sites$SHM==F,"mature_mutability"]) #SHMs are not more inherently mutable sites?
#dev.off()


#same, but actual mutations at these positions
betas[,site_SHM := sites[sites$site==site,"SHM"],by=mutation]
ggplot(betas[!is.na(site_SHM),], aes(x=site_SHM,y=bind_observed))+
  geom_violin(trim=F)+coord_flip()+
  geom_dotplot(binaxis='y',stackdir='center',dotsize=0.01)
wilcox.test(betas[!is.na(site_SHM) & site_SHM==T,bind_observed],betas[!is.na(site_SHM) & site_SHM==F,bind_observed]) #SHMs accrued at more tolerant positions


#pdf(file="results/analyze_mut_effects/histogram_ddG-of-reversion-to-naive-state_latent-scale.pdf",width=4.5,height=4.5,useDingbats=F)
hist(sites[!is.na(sites$SHM) & sites$SHM==T,"naive_ddG"],breaks=20,freq=T,col="gray40",ylim=c(0,30),ylab="number somatic hypermutations",xlab="ddG of reversion to naive state",main=paste("sum = ",sum(sites[!is.na(sites$SHM) & sites$SHM==T,"naive_ddG"],na.rm=T),sep=""))
sum(sites[!is.na(sites$SHM) & sites$SHM==T,"naive_ddG"],na.rm=T)
#dev.off()
#set.seed(1042);hist(betas[betas$site %in% 126:140,"bind_observed"][sample(1:89,65)],add=T,col="#0000ff50",breaks=5) #perhaps an example null distribution to compare to?
#many SHMs have small or neutral effects; few have large effects on affinity

#pdf(file="results/analyze_mut_effects/histogram_ddG-of-reversion-to-naive-state_observed-scale.pdf",width=4.5,height=4.5,useDingbats=F)
hist(sites[!is.na(sites$SHM) & sites$SHM==T,"naive_ddG_observed"],breaks=20,freq=T,col="gray40",ylim=c(0,30),ylab="number somatic hypermutations",xlab="ddG of reversion to naive state",main=paste("sum = ",sum(sites[!is.na(sites$SHM) & sites$SHM==T,"naive_ddG_observed"],na.rm=T),sep=""))
set.seed(285802);hist(rnorm(40,mean=0,sd=pooled_SEM_bind),add=T,col="#0000ff50",breaks=3)
#hist(sites[sites$minVRC01_AA == sites$AA & sites$minVRC01_AA != sites$naive_AA & !is.na(sites$minVRC01_AA),"naive_ddG_observed"],add=T,col="green",breaks=20)
#dev.off()

#compare to differences with minVRC01? See in the histogram above tha tthe three most deleterious reversions are all still retained in minVRC01, like we'd expect

#next, look at diversity across VRC01-class Abs from within and the d45 lineage as well as from other donors -- is there more incompatibilities inter-lineage?
intra_VRC01class <- read.fasta(file="data/Ab_alignments/intra_align_NIH45-46-numbering.fasta") #alignment of sequences from the NIH45-46 lineage. Includes probe-isolated and functionally validated NGS sequences
inter_VRC01class <- read.fasta(file="data/Ab_alignments/inter_VRC01-class_align_NIH45-46-numbering.fasta") #alignment of sequences of VRC01-class antibodies from other donors

#alignment entropy at each site in each of these alignments -- compare to mean ddG of site in NIH45-46
intra_entropy <- entropy(intra_VRC01class)$H
inter_entropy <- entropy(inter_VRC01class)$H
for(i in 1:nrow(sites)){
  sites$intra_entropy[i] <- intra_entropy[i]
  sites$inter_entropy[i] <- inter_entropy[i]
  if(sites$site[i] %in% c(126:140,242:244)){
    sites$intra_entropy[i] <- NA
    sites$inter_entropy[i] <- NA
  }
}
plot(sites$inter_entropy,sites$intra_entropy,pch=16)
plot(sites$intra_entropy,sites$mean_bind_observed,pch=16,ylab="mean ddG (kcal/mol)",xlab="conservation (alignment entropy), within lineage")
plot(sites$inter_entropy,sites$mean_bind_observed,pch=16,ylab="mean ddG (kcal/mol)",xlab="conservation (alignment entropy), between lineages")

intra_subs <- vector()
for(j in 1:ncol(intra_VRC01class$ali)){
  for(i in 2:nrow(intra_VRC01class$ali)){
    if(intra_VRC01class$ali[i,j] != "-"){
      if(intra_VRC01class$ali[i,j] != intra_VRC01class$ali[1,j]){
        intra_subs <- c(intra_subs, unlist(paste(intra_VRC01class$ali[1,j],j,intra_VRC01class$ali[i,j],sep="")))
      }
    }
  }
}
intra_subs_unique <- unique(intra_subs)

inter_subs <- vector()
for(j in 1:ncol(inter_VRC01class$ali)){
  for(i in 2:nrow(inter_VRC01class$ali)){
    if(inter_VRC01class$ali[i,j] != "-"){
      if(inter_VRC01class$ali[i,j] != inter_VRC01class$ali[1,j]){
        inter_subs <- c(inter_subs, unlist(paste(inter_VRC01class$ali[1,j],j,inter_VRC01class$ali[i,j],sep="")))
      }
    }
  }
}
inter_subs_unique <- unique(inter_subs)

intra_subs <- data.frame(mutation=intra_subs,stringsAsFactors = F)
inter_subs <- data.frame(mutation=inter_subs,stringsAsFactors = F)

for(i in 1:nrow(intra_subs)){
  intra_subs$wildtype[i] <- strsplit(intra_subs$mutation[i],split="")[[1]][1]
  intra_subs$mutant[i] <- strsplit(intra_subs$mutation[i],split="")[[1]][length(strsplit(intra_subs$mutation[i],split="")[[1]])]
  intra_subs$site[i] <- as.numeric(paste(strsplit(intra_subs$mutation[i],split="")[[1]][-c(1,length(strsplit(intra_subs$mutation[i],split="")[[1]]))],collapse=""))
  intra_subs$pdb_site[i] <- sites[sites$site==intra_subs$site[i],"pdb_site"]
  intra_subs$pdb_chain[i] <- sites[sites$site==intra_subs$site[i],"pdb_chain"]
  if(length(betas[mutation==intra_subs$mutation[i],bind_observed])>0){
    intra_subs$ddG[i] <- betas[mutation==intra_subs$mutation[i],bind_latent]
    intra_subs$ddG_observed[i] <- betas[mutation==intra_subs$mutation[i],bind_observed]
  }else{
    intra_subs$ddG[i] <- NA
    intra_subs$ddG_observed[i] <- NA
  }
}

for(i in 1:nrow(inter_subs)){
  inter_subs$wildtype[i] <- strsplit(inter_subs$mutation[i],split="")[[1]][1]
  inter_subs$mutant[i] <- strsplit(inter_subs$mutation[i],split="")[[1]][length(strsplit(inter_subs$mutation[i],split="")[[1]])]
  inter_subs$site[i] <- as.numeric(paste(strsplit(inter_subs$mutation[i],split="")[[1]][-c(1,length(strsplit(inter_subs$mutation[i],split="")[[1]]))],collapse=""))
  inter_subs$pdb_site[i] <- sites[sites$site==inter_subs$site[i],"pdb_site"]
  inter_subs$pdb_chain[i] <- sites[sites$site==inter_subs$site[i],"pdb_chain"]
  if(length(betas[mutation==inter_subs$mutation[i],bind_observed])>0){
    inter_subs$ddG[i] <- betas[mutation==inter_subs$mutation[i],bind_latent]
    inter_subs$ddG_observed[i] <- betas[mutation==inter_subs$mutation[i],bind_observed]
  }else{
    inter_subs$ddG[i] <- NA
    inter_subs$ddG_observed[i] <- NA
  }
}

#combine dataframes with indicator as to whether it's intra or inter alignment
intra_subs$alignment <- "intra"
inter_subs$alignment <- "inter"
subs <- rbind(intra_subs,inter_subs)

boxplot(subs$ddG_observed ~ subs$alignment + subs$pdb_chain,notch=T)

ggplot(subs, aes(x=alignment,y=ddG_observed))+
  geom_violin(trim=F)+coord_flip()+
  geom_dotplot(binaxis='y',stackdir='center',dotsize=0.01)
wilcox.test(subs[!is.na(subs$alignment) & subs$alignment=="intra","ddG_observed"],subs[!is.na(subs$alignment) & subs$alignment=="inter","ddG_observed"])

#pdf(file="results/analyze_mut_effects/hist_intra-lineage-differences_ddG_latent-scale.pdf",width=5,height=4.5,useDingbats=F)
hist(intra_subs$ddG,ylab="number amino acid differences",xlab="ddG (kcal/mol)",main="differences within the NIH45-46 lineage",breaks=30,col="gray40")
#dev.off()
unique(intra_subs[intra_subs$ddG>1.4,]) #looks like a lot of contact residues? maybe site 33H an intrinsic thing

#pdf(file="results/analyze_mut_effects/hist_intra-lineage-differences_ddG_observed-scale.pdf",width=5,height=4.5,useDingbats=F)
hist(intra_subs$ddG_observed,ylab="number amino acid differences",xlab="ddG (kcal/mol)",main="differences within the NIH45-46 lineage",breaks=30,col="gray40")
set.seed(19714);hist(rnorm(2500,mean=0,sd=pooled_SEM_bind),add=T,col="#0000ff50",breaks=3)
#dev.off()
unique(intra_subs[intra_subs$ddG_observed>1,]) #looks like a lot of contact residues? maybe site 33H an intrinsic thing
t.test(intra_subs$ddG_observed);t.test(intra_subs$ddG_observed[intra_subs$ddG_observed < 1.4]) #shifted toward deleterious effects, even when excluding outliers

#pdf(file="results/analyze_mut_effects/hist_inter-lineage-differences_ddG_latent-scale.pdf",width=5,height=4.5,useDingbats=F)
hist(inter_subs$ddG,ylab="number amino acid differences",xlab="ddG (kcal/mol)",main="differences among heterologous VRC01-class Abs",breaks=30,col="gray40")
#dev.off()
unique(inter_subs[inter_subs$ddG>1.4,]) #looks like a lot of contact residues?

#pdf(file="results/analyze_mut_effects/hist_inter-lineage-differences_ddG_observed-scale.pdf",width=5,height=4.5,useDingbats=F)
hist(inter_subs$ddG_observed,ylab="number amino acid differences",xlab="ddG (kcal/mol)",main="differences among heterologous VRC01-class Abs",breaks=30,col="gray40")
set.seed(47194);hist(rnorm(800,mean=0,sd=pooled_SEM_bind),add=T,col="#0000ff50",breaks=2)
#dev.off()
unique(inter_subs[inter_subs$ddG_observed>1,]) #looks like a lot of contact residues?
t.test(inter_subs$ddG_observed);t.test(inter_subs$ddG_observed[inter_subs$ddG_observed < 1.4]) #shifted toward deleterious effects, even when excluding outliers



#make dataframe of Calpha and minimum atomic distances between each pair of residues in 3u7y NIH45-46 crystal structure
nsites <- 244
pdb_dists <- expand.grid(site1=1:nsites,site2=1:nsites)
pdb_dists <- pdb_dists[order(pdb_dists$site1, pdb_dists$site2),]
pdb_dists <- pdb_dists[pdb_dists$site1 < pdb_dists$site2,]
for(i in 1:nrow(pdb_dists)){
  pdb_dists$site1_PDB[i] <- sites_alignment[sites_alignment$site_scFv==pdb_dists$site1[i],"site_3U7Y"]
  pdb_dists$site1_chain[i] <- sites_alignment[sites_alignment$site_scFv==pdb_dists$site1[i],"chain_3U7Y"]
  pdb_dists$site2_PDB[i] <- sites_alignment[sites_alignment$site_scFv==pdb_dists$site2[i],"site_3U7Y"]
  pdb_dists$site2_chain[i] <- sites_alignment[sites_alignment$site_scFv==pdb_dists$site2[i],"chain_3U7Y"]
}

calc.dist <- function(x1,y1,z1,x2,y2,z2){
  return(sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2))
}

for(i in 1:nrow(pdb_dists)){
  print(i)
  if(!is.na(pdb_dists[i,"site1_PDB"]) & !is.na(pdb_dists[i,"site2_PDB"])){
    atoms1 <- pdb_atoms[pdb_atoms$chain==pdb_dists[i,"site1_chain"] & pdb_atoms$resno==pdb_dists[i,"site1_PDB"],]
    atoms2 <- pdb_atoms[pdb_atoms$chain==pdb_dists[i,"site2_chain"] & pdb_atoms$resno==pdb_dists[i,"site2_PDB"],]
    if(nrow(atoms1)>0 & nrow(atoms2)>0){
      pdb_dists$CA_dist[i] <- calc.dist(atoms1[atoms1$elety=="CA","x"],atoms1[atoms1$elety=="CA","y"],atoms1[atoms1$elety=="CA","z"],
                                        atoms2[atoms2$elety=="CA","x"],atoms2[atoms2$elety=="CA","y"],atoms2[atoms2$elety=="CA","z"])
      all_dists <- c()
      for(i1 in 1:nrow(atoms1)){for(i2 in 1:nrow(atoms2)){
        all_dists <- c(all_dists,calc.dist(atoms1[i1,"x"],atoms1[i1,"y"],atoms1[i1,"z"],atoms2[i2,"x"],atoms2[i2,"y"],atoms2[i2,"z"]))}}
      pdb_dists$min_dist[i] <- min(all_dists)
    }else{pdb_dists$CA_dist[i] <- NA;pdb_dists$min_dist[i] <- NA}
  }else{pdb_dists$CA_dist[i] <- NA;pdb_dists$min_dist[i] <- NA}
}
pdb_dists <- data.table(pdb_dists)

#can use these distances to check out possible epistatic candidates?
#line to put in a (pdb)site/chain and get residues less than a certain distance
site <- 33; chain <- "H"; distance <- 6
head(pdb_dists)
pdb_dists[((site2_chain==chain & site2_PDB==site) | (site1_chain==chain & site1_PDB==site)) & min_dist < distance,]

#e.g. list of differences between NIH45-46 and 3BNC60
diffs_NIH_3BNC60 <- vector()
for(j in 1:ncol(inter_VRC01class$ali)){
  if(inter_VRC01class$ali[9,j] != "-"){
    if(inter_VRC01class$ali[9,j] != inter_VRC01class$ali[1,j]){
      diffs_NIH_3BNC60 <- c(diffs_NIH_3BNC60, unlist(paste(inter_VRC01class$ali[1,j],j,inter_VRC01class$ali[9,j],sep="")))
    }
  }
}
diffs_NIH_3BNC60 <- data.frame(mutation=diffs_NIH_3BNC60,stringsAsFactors = F)
for(i in 1:nrow(diffs_NIH_3BNC60)){
  diffs_NIH_3BNC60$wildtype[i] <- strsplit(diffs_NIH_3BNC60$mutation[i],split="")[[1]][1]
  diffs_NIH_3BNC60$mutant[i] <- strsplit(diffs_NIH_3BNC60$mutation[i],split="")[[1]][length(strsplit(diffs_NIH_3BNC60$mutation[i],split="")[[1]])]
  diffs_NIH_3BNC60$site[i] <- as.numeric(paste(strsplit(diffs_NIH_3BNC60$mutation[i],split="")[[1]][-c(1,length(strsplit(diffs_NIH_3BNC60$mutation[i],split="")[[1]]))],collapse=""))
  diffs_NIH_3BNC60$pdb_site[i] <- sites[sites$site==diffs_NIH_3BNC60$site[i],"pdb_site"]
  diffs_NIH_3BNC60$pdb_chain[i] <- sites[sites$site==diffs_NIH_3BNC60$site[i],"pdb_chain"]
  if(length(betas[betas$mutation==diffs_NIH_3BNC60$mutation[i],"bind_observed"])>0){
    diffs_NIH_3BNC60$ddG[i] <- betas[betas$mutation==diffs_NIH_3BNC60$mutation[i],"bind_latent"]
    diffs_NIH_3BNC60$ddG_observed[i] <- betas[betas$mutation==diffs_NIH_3BNC60$mutation[i],"bind_observed"]
  }else{
    diffs_NIH_3BNC60$ddG[i] <- NA
    diffs_NIH_3BNC60$ddG_observed[i] <- NA
  }
}
diffs_NIH_3BNC60[diffs_NIH_3BNC60$pdb_chain=="H","pdb_site"]

#e.g. list of differences between NIH45-46 and VRC01
diffs_NIH_VRC01 <- vector()
for(j in 1:ncol(intra_VRC01class$ali)){
  if(intra_VRC01class$ali[2,j] != "-"){
    if(intra_VRC01class$ali[2,j] != intra_VRC01class$ali[1,j]){
      diffs_NIH_VRC01 <- c(diffs_NIH_VRC01, unlist(paste(intra_VRC01class$ali[1,j],j,intra_VRC01class$ali[2,j],sep="")))
    }
  }
}
diffs_NIH_VRC01 <- data.frame(mutation=diffs_NIH_VRC01,stringsAsFactors = F)
for(i in 1:nrow(diffs_NIH_VRC01)){
  diffs_NIH_VRC01$wildtype[i] <- strsplit(diffs_NIH_VRC01$mutation[i],split="")[[1]][1]
  diffs_NIH_VRC01$mutant[i] <- strsplit(diffs_NIH_VRC01$mutation[i],split="")[[1]][length(strsplit(diffs_NIH_VRC01$mutation[i],split="")[[1]])]
  diffs_NIH_VRC01$site[i] <- as.numeric(paste(strsplit(diffs_NIH_VRC01$mutation[i],split="")[[1]][-c(1,length(strsplit(diffs_NIH_VRC01$mutation[i],split="")[[1]]))],collapse=""))
  diffs_NIH_VRC01$pdb_site[i] <- sites[sites$site==diffs_NIH_VRC01$site[i],"pdb_site"]
  diffs_NIH_VRC01$pdb_chain[i] <- sites[sites$site==diffs_NIH_VRC01$site[i],"pdb_chain"]
  if(length(betas[betas$mutation==diffs_NIH_VRC01$mutation[i],"bind_observed"])>0){
    diffs_NIH_VRC01$ddG[i] <- betas[betas$mutation==diffs_NIH_VRC01$mutation[i],"bind_latent"]
    diffs_NIH_VRC01$ddG_observed[i] <- betas[betas$mutation==diffs_NIH_VRC01$mutation[i],"bind_observed"]
  }else{
    diffs_NIH_VRC01$ddG[i] <- NA
    diffs_NIH_VRC01$ddG_observed[i] <- NA
  }
}
diffs_NIH_VRC01[diffs_NIH_VRC01$pdb_chain=="H","pdb_site"]

#looking at effects by site, mutation
#order mutant as a factor for grouping by biochemical group
betas$mutant <- factor(betas$mutant, levels=c("*","C","P","G","V","M","L","I","A","F","W","Y","T","S","N","Q","E","D","H","K","R"))
#make factors for kabat numberings to keep ordering on display
betas$kabat_annotation <- factor(betas$kabat_annotation,levels=unique(betas$kabat_annotation))
betas$kabat_index <- factor(paste(betas$kabat_chain,betas$kabat_site,sep=""),levels=unique(paste(betas$kabat_chain,betas$kabat_site,sep="")))
#add character vector indicating widltype to use as plotting symbols for wt
betas[,wildtype_indicator := as.character(NA)]
betas[mutant==wildtype,wildtype_indicator := "x"]

#violin plots by annotations: FWR/CDR, HC/LC, interface/scaffold, buried/surface exposed, diagnostic of VRC01-class, sampled in lineages, etc.
boxplot(betas$bind_observed ~ betas$kabat_annotation,notch=T,pch=16)
boxplot(betas$bind_observed ~ betas$kabat_chain,notch=T)
wilcox.test(betas[kabat_chain=="H",bind_observed],betas[kabat_chain=="L",bind_observed])

##relationship between mutational effects on binding and expression
plot(betas$bind_observed,betas$expr_observed,pch=16,col="#00000020")
plot(betas$bind_latent,betas$expr_observed,pch=16,col="#00000020")
plot(betas$bind_observed,betas$expr_latent,pch=16,col="#00000020",ylim=c(-6,1))
plot(betas$bind_latent,betas$expr_latent,pch=16,col="#00000020",ylim=c(-6,1),xlab="binding betas, latent scale",ylab="expression betas, latent scale")
#notice that the coupling between mutant effects on binding and exprssion goes to binding of about +2, which is where the global epistasis curve for binding inflects -- need to think about how to explain this but intuitively it seems to be saying that mutational effects on binding are not decoupled from expression with Titeseq?

#correlation between mutational effect on expression and the response variable from the titration fit?
plot(bc_expr$latent_phenotype_Cauchy_A,bc_bind$response,pch=16,col="#00000005",xlab="latent phenotype, expression",ylab="Titration curve fit plateau height",xlim=c(-3,1))
plot(bc_expr$predicted_phenotype_Cauchy_A,bc_bind$response,pch=16,col="#00000005",xlab="predicted phenotype, expression",ylab="Titration curve fit plateau height")
plot(bc_expr$latent_phenotype_Cauchy_B,bc_bind$response,pch=16,col="#00000005",xlab="latent phenotype, expression",ylab="Titration curve fit plateau height",xlim=c(-3,1))
plot(bc_expr$predicted_phenotype_Cauchy_B,bc_bind$response,pch=16,col="#00000005",xlab="predicted phenotype, expression",ylab="Titration curve fit plateau height")

#biochemical determinants of per-site mutational effects
#in addition to heatmaps (below), I want to do something similar to the analysis of Abriata et al. PLoS ONE 2015, where they test a suite of biochemical characteristics to see which provides the best linear fit with the site-specific mutational effects
AA_props <- read.csv(file="data/2015_Abriata_AA-characteristics.csv")
#two options are, to make a multiple linear regression with a subset of these characteristics that are the least autocorrelated, or as they do in Abriata, test each correlation and keep the best correlation (but also require RMSD < 1 or whatever to avoid outliers driving the correlation)
#see what individual examples look like
#W50h, N58h, R71h, W100Fh (100B in VRC01?), W67l, CDRL3 (positions 91 and 96 side chains

dt <- betas[kabat_index=="H58" & mutant!="*",]
col <- 12;x <- AA_props[AA_props$Residue==as.character(dt[,mutant]),col];y <- dt[,bind_observed];model <- lm(y~x);summary(model);plot(x,y,ylab="ddG",xlab=names(AA_props)[col],ylim=rev(range(y,na.rm=T)),pch=as.character(dt$mutant));abline(model,lty=2,col="red")
col <- "Volume";x <- AA_props[AA_props$Residue==as.character(dt[,mutant]),col];y <- dt[,bind_observed];model <- lm(y~x);summary(model);plot(x,y,ylab="ddG",xlab=names(AA_props)[col],ylim=rev(range(y,na.rm=T)),pch=as.character(dt$mutant));abline(model,lty=2,col="red")
col <- "Hydrophobicity";x <- AA_props[AA_props$Residue==as.character(dt[,mutant]),col];y <- dt[,bind_observed];model <- lm(y~x);summary(model);plot(x,y,ylab="ddG",xlab=names(AA_props)[col],ylim=rev(range(y,na.rm=T)),pch=as.character(dt$mutant));abline(model,lty=2,col="red")
col <- "Isoelectric.point";x <- AA_props[AA_props$Residue==as.character(dt[,mutant]),col];y <- dt[,bind_observed];model <- lm(y~x);summary(model);plot(x,y,ylab="ddG",xlab=names(AA_props)[col],ylim=rev(range(y,na.rm=T)),pch=as.character(dt$mutant));abline(model,lty=2,col="red")

#good fits have RSE around 1-1.25 (see site 230), R-squared ~0.5
#decent fits have RSE 1.5-2 (if there are highly del mutations), R-squared 0.46

#for now, keep anything with RSE < 2 and R-squared >0.3, and only keep multiple properties if they're within 0.05 R-squared to the "best" correlate

for(i in 1:nrow(sites)){ #for each site
  print(i)
  if(!is.na(sites$kabat_chain[i])){
    properties <- vector() #list of properties that match fit quality criteria
    RSEs <- vector()
    Rsquareds <- vector()
    slopes <- vector()
    for(j in 2:ncol(AA_props)){ #for each AA property
      model <- lm(betas[mutant != "*" & site==sites[i,"site"],bind_observed] ~ AA_props[AA_props$Residue==as.character(betas[mutant != "*" & site==sites[i,"site"],mutant]),j])
      RSE <- summary(model)$sigma
      Rsquared <- summary(model)$r.squared
      slope <- summary(model)$coefficient[2,1]
      if(RSE < 2 & Rsquared > 0.3){
        properties <- c(properties,colnames(AA_props)[j])
        RSEs <- c(RSEs,RSE)
        Rsquareds <- c(Rsquareds, Rsquared)
        slopes <- c(slopes, slope)
      }
    }
    properties <- properties[which(Rsquareds > (max(Rsquareds)-0.05))] #only keep properties within 5% variance explained compared to best property
    RSEs <- RSEs[which(Rsquareds > (max(Rsquareds)-0.05))]
    slopes <- slopes[which(Rsquareds > (max(Rsquareds)-0.05))]
    Rsquareds <- Rsquareds[which(Rsquareds > (max(Rsquareds)-0.05))]
    sites$n_correlates[i] <- length(properties)
    if(sites$n_correlates[i]>0){
      sites$correlates[i] <- list(properties)
      sites$correlate_RSEs[i] <- list(RSEs)
      sites$correlate_Rsquareds[i] <- list(Rsquareds)
      sites$correlate_slopes[i] <- list(slopes)
      sites$best_correlate[i] <- properties[which(Rsquareds==max(Rsquareds))]
      sites$best_Rsquared[i] <- Rsquareds[which(Rsquareds==max(Rsquareds))]
      sites$best_slope[i] <- slopes[which(Rsquareds==max(Rsquareds))]
    }else{
      sites$correlates[i] <- NA
      sites$correlate_RSEs[i] <- NA
      sites$correlate_Rsquareds[i] <- NA
      sites$correlate_slopes[i] <- NA
      sites$best_correlate[i] <- NA
      sites$best_Rsquared[i] <- NA
      sites$best_slope[i] <- NA
    }
  }else{
    if(is.na(sites$kabat_chain[i])){
      sites$n_correlates[i] <- NA
      sites$correlates[i] <- NA
      sites$correlate_RSEs[i] <- NA
      sites$correlate_Rsquareds[i] <- NA
      sites$correlate_slopes[i] <- NA
    }
  }
}

dt <- betas[site==110 & mutant!="*",];sites[sites$site==110,]
col <- "Volume";x <- AA_props[AA_props$Residue==as.character(dt[,mutant]),col];y <- dt[,bind_observed];model <- lm(y~x);summary(model);plot(x,y,ylab="ddG",xlab=names(AA_props)[col],ylim=rev(range(y,na.rm=T)),pch=as.character(dt$mutant));abline(model,lty=2,col="red")

#make a data frame that captures the biochemical parameters (and the direction of their effect) for each site, similar to the betas frame, for plotting heatmaps; make it wide, then convert to long dataframe for plotiting
site_correlates <- sites[!is.na(sites$kabat_chain),c("site","AA","kabat_site","kabat_chain","kabat_annotation","pdb_site","pdb_chain")]
for(i in 2:ncol(AA_props)){ #add columns of 0 for each property
  site_correlates[,names(AA_props)[i]] <- rep(0,nrow(site_correlates))
}
#if helix propensity isn't in the list of factors that are ever occupied, don't bother plotting
if(!("Helix_propensity" %in% unique(unlist(sites$correlates)))){
  site_correlates$Helix_propensity <- NULL
}
for(i in 1:nrow(site_correlates)){
  properties <- sites[sites$site==site_correlates[i,"site"],"correlates"][[1]]
  slopes <- sites[sites$site==site_correlates[i,"site"],"correlate_slopes"][[1]]
  if(sum(!is.na(properties))>0){
    for(j in properties){
      if(slopes[which(properties==j)]>0){
        site_correlates[i,j] <- 1
      }else if(slopes[which(properties==j)]<0){
        site_correlates[i,j] <- -1
      }
    }
  }
}

site_correlates <- melt(site_correlates,id.vars=c("site","AA","kabat_site","kabat_chain","kabat_annotation","pdb_site","pdb_chain"),variable.name="property",value.name="indicator")
site_correlates$kabat_index <- factor(paste(site_correlates$kabat_chain,site_correlates$kabat_site,sep=""),levels=unique(paste(site_correlates$kabat_chain,site_correlates$kabat_site,sep="")))
site_correlates <- data.table(site_correlates)

#output heat maps of biochemical correlates. I don't know how to properly arrange these on top of each other easily, easier to just manually align post-hoc?
#pdf(file="results/analyze_mut_effects/heatmap_HC_biochem-correlates.pdf",width=7.5,height=4,bg="transparent")
ggplot(site_correlates[kabat_chain=="H",],aes(kabat_index,property))+geom_tile(aes(fill=factor(indicator)))+
  scale_fill_manual(values=c("#7378B9","#FFFFFF","#F48365"),breaks=c("-1","0","1"))+
  labs(x="site (Kabat numbering)",y="property")+theme_classic(base_size=4.5)+
  coord_equal()+
  theme(axis.text.x = element_text(angle=45,hjust=1))
#dev.off()

#pdf(file="results/analyze_mut_effects/heatmap_LC_biochem-correlates.pdf",width=7.5,height=4,bg="transparent")
ggplot(site_correlates[kabat_chain=="L",],aes(kabat_index,property))+geom_tile(aes(fill=factor(indicator)))+
  scale_fill_manual(values=c("#7378B9","#FFFFFF","#F48365"),breaks=c("-1","0","1"))+
  labs(x="site (Kabat numbering)",y="property")+theme_classic(base_size=4.5)+
  coord_equal()+
  theme(axis.text.x = element_text(angle=45,hjust=1))
#dev.off()

# #multiple linear regression using the parameters of Volume, log.solubility, Hydrophobicity, Isoelectric point, and helix propensity which have <0.3 self-correlation amongst themselves
# #center and normalize the three parameters from AA_props
# AA_props_normalize <- data.frame(Residue=AA_props[,"Residue"])
# AA_props_normalize$Volume <- (AA_props$Volume-mean(AA_props$Volume))/sd(AA_props$Volume)
# AA_props_normalize$Log.solubility <- (AA_props$Log.Solubility..m.at.20ºC..-mean(AA_props$Log.Solubility..m.at.20ºC..))/sd(AA_props$Log.Solubility..m.at.20ºC..)
# AA_props_normalize$Hydrophobicity <- (AA_props$Hydrophobicity-mean(AA_props$Hydrophobicity))/sd(AA_props$Hydrophobicity)
# AA_props_normalize$Isoelectric.point <- (AA_props$Isoelectric.point-mean(AA_props$Isoelectric.point))/sd(AA_props$Isoelectric.point)
# AA_props_normalize$Helix.propensity <- (AA_props$P.helix.-mean(AA_props$P.helix.))/sd(AA_props$P.helix.)
# 
# dt <- betas[site==50 & mutant!="*",]
# ddG <- dt[,bind_observed]; Volume <- AA_props_normalize[AA_props_normalize$Residue==as.character(dt[,mutant]),"Volume"]; Solubility <- AA_props_normalize[AA_props_normalize$Residue==as.character(dt[,mutant]),"Log.solubility"]; Hydrophobicity <- AA_props_normalize[AA_props_normalize$Residue==as.character(dt[,mutant]),"Hydrophobicity"]; Isoelectric_point <- AA_props_normalize[AA_props_normalize$Residue==as.character(dt[,mutant]),"Isoelectric.point"]; Helix <- AA_props_normalize[AA_props_normalize$Residue==as.character(dt[,mutant]),"Helix.propensity"];  model <- lm(ddG ~ Volume+Solubility+Hydrophobicity+Isoelectric_point+Helix);summary(model)
# summary(model)
# plot(model)
# 
# col <- "Volume";x <- AA_props[AA_props$Residue==as.character(dt[,mutant]),col];y <- dt[,bind_observed];model <- lm(y~x);summary(model);plot(x,y,ylab="ddG",xlab=names(AA_props)[col],ylim=rev(range(y,na.rm=T)),pch=as.character(dt$mutant));abline(model,lty=2,col="red")
# col <- "Hydrophobicity";x <- AA_props[AA_props$Residue==as.character(dt[,mutant]),col];y <- dt[,bind_observed];model <- lm(y~x);summary(model);plot(x,y,ylab="ddG",xlab=names(AA_props)[col],ylim=rev(range(y,na.rm=T)),pch=as.character(dt$mutant));abline(model,lty=2,col="red")
# col <- "Isoelectric.point";x <- AA_props[AA_props$Residue==as.character(dt[,mutant]),col];y <- dt[,bind_observed];model <- lm(y~x);summary(model);plot(x,y,ylab="ddG",xlab=names(AA_props)[col],ylim=rev(range(y,na.rm=T)),pch=as.character(dt$mutant));abline(model,lty=2,col="red")
# 
# ggplot(betas[mutant!="*" & site %in% c(72),],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
#   scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
#   labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
#   coord_equal()+
#   geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
#   theme(axis.text.x = element_text(angle=45,hjust=1))


#heatmaps of mutational effects, also annotate with the "best" biochemical correlate?
#add best correlate column
for(i in 1:nrow(betas)){
  if(betas[i,mutant]=="A"){
    betas$best_correlate[i] <- as.character(sites[sites$site==betas$site[i],"best_correlate"])
  }else{
    betas$best_correlate[i] <- as.character(NA)
  }
}


#all sites, binding
#pdf(file="results/analyze_mut_effects/heatmap_all-positions.pdf",width=12,height=4,bg="transparent")
ggplot(betas[mutant!="*" & !(site %in% seq(126,140)),],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=4)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=1,color="gray50")
  theme(axis.text.x = element_text(angle=45,hjust=1))
#dev.off()
#HC, binding
#pdf(file="results/analyze_mut_effects/heatmap_HC.pdf",width=7.5,height=4,bg="transparent")
ggplot(betas[kabat_chain=="H" & mutant!="*" & !(site %in% seq(126,140)),],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=4.5)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=1,color="gray50")+
  #geom_text(aes(x=kabat_index,y=length(unique(mutant))+0.7,label=best_correlate,angle=45))+
  theme(axis.text.x = element_text(angle=45,hjust=1))
#dev.off()

#LC, binding
#pdf(file="results/analyze_mut_effects/heatmap_LC.pdf",width=7.5,height=4,bg="transparent")
ggplot(betas[kabat_chain=="L" & mutant!="*" & !(site %in% seq(126,140)),],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=5)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=1,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))
#dev.off()
#all sites, expression
ggplot(betas[!(site %in% seq(126,140)),],aes(kabat_index,mutant))+geom_tile(aes(fill=expr_latent))+
  scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-15,2),values=c(0,7/17,11/17,15/17,16/17,1),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=4)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=1,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))
#anything below 7 equally red
#HC, expression
ggplot(betas[kabat_chain=="H" & !(site %in% seq(126,140)),],aes(kabat_index,mutant))+geom_tile(aes(fill=expr_latent))+
  scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-15,2),values=c(0,7/17,11/17,15/17,16/17,1),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=4)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=1,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))
#anything below 7 equally red
#LC, expression
ggplot(betas[kabat_chain=="L" & !(site %in% seq(126,140)),],aes(kabat_index,mutant))+geom_tile(aes(fill=expr_latent))+
  scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-15,2),values=c(0,7/17,11/17,15/17,16/17,1),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=4)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=1,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))
#anything below 7 equally red


#do heat maps look more valuable faceted on just "interesting" residues? e.g. contact residues, those known to be important in the literature, etc.
#note, on heatmaps, can I add an X for the "wildtype" state?
#CDRs
ggplot(betas[mutant!="*" & kabat_annotation %in% c("CDRH1","CDRH2","CDRH3","CDRL1","CDRL2","CDRL3"),],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#CDRs, one at a time
#CDRH1 residues
ggplot(betas[mutant!="*" & kabat_annotation=="CDRH1",],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#CDRH2 residues
ggplot(betas[mutant!="*" & kabat_annotation=="CDRH2",],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#CDRH3 residues
ggplot(betas[mutant!="*" & kabat_annotation=="CDRH3",],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#CDRL1 residues
ggplot(betas[mutant!="*" & kabat_annotation=="CDRL1",],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#CDRL2 residues
ggplot(betas[mutant!="*" & kabat_annotation=="CDRL2",],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#CDRL3 residues
ggplot(betas[mutant!="*" & kabat_annotation=="CDRL3",],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#sites of beneficial mutations
ggplot(betas[mutant!="*" & site %in% c(28,51,55,77,168,170),],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#"diagnostic" VRC01-class residues...
#West et al. PNAS focus on structural roles of residues W50h, N58h, R71h, W100Fh (100B in VRC01?), W67l, CDRL3 (positions 91 and 96 side chains)
#note that W67 thought to perhaps contact N276 glycan? Absent in my construct
#Diskin 2011, Scharf 2016 mention G54 makes backbone contact (and has many beneficial mutations available)
#Wiehe et al. 2018 say that reversions of improbable mutations affect neutralization: E16A, E28T (HC), I21L, W67S, N72T, Y71F (LC)
#F112 (kabat 100H), seen in Havenar-Davenaugh among naive VRC01-class Abs to be enriched in F
ggplot(betas[mutant!="*" & (kabat_index %in% c("H50","H54","H58","H71","H100F","L67") | kabat_annotation %in% c("CDRL3")),],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(betas[mutant!="*" & kabat_index %in% c("H8", "H9", "H10","H11"),],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#minimal SHMs (minVRC01 from Jardine et al.)
ggplot(betas[mutant!="*" & site %in% c(32,33,34,35,52,54,57,58,62,74,75,76,77,165,168,169,189,204,205,207),],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))



#sites that bury surface area at interface
ggplot(betas[mutant!="*" & site %in% sites[sites$buried_RSA>0,"site"],],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#sites that are buried in free structure
ggplot(betas[mutant!="*" & site %in% sites[sites$RSA_unbound==0,"site"],],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#"sensitive" positions
ggplot(betas[mutant!="*" & site %in% sites[sites$mean_bind_observed>2,"site"],],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#sites that changed along lineage
ggplot(betas[mutant!="*" & site %in% sites[sites$SHM==T & !is.na(sites$SHM),"site"],],aes(kabat_index,mutant))+geom_tile(aes(fill=bind_observed))+
  scale_fill_gradientn(colours=c("#383C6C","#7378B9","#FFFFFF","#F48365","#A94E35"),limits=c(-1,8),values=c(0,0.5/9,1/9,5/9,9/9),na.value="gray")+
  labs(x="site (Kabat numbering)",y="mutant")+theme_classic(base_size=8)+
  coord_equal()+
  geom_text(aes(label=wildtype_indicator),size=3,color="gray50")+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#structural distribution of mutational effects
#pdf(file="results/analyze_mut_effects/line-plot_mean-ddG_and_mean-expr_per_site_latent-scales.pdf",width=9,height=3.5,bg="transparent")
par(mar=c(5,5,2,5))
plot(sites$site, sites$mean_bind,type="l",xlab="site",ylab="mean ddG",ylim=c(5.5,-1),col="magenta")
par(new=T)
plot(sites$site, sites$mean_expr,type="l",col="green",axes=F,xlab=NA,ylab=NA,ylim=c(-6,1.2))
axis(side=4)
mtext(side=4, line=3, "mean change in expression (ML meanF)")
points(sites$site, rep(0.75,nsites),pch=as.character(sites$AA),cex=0.25)
#dev.off()

#pdf(file="results/analyze_mut_effects/line-plot_mean-ddG_and_mean-expr_per_site_observed-scales.pdf",width=9,height=3.5,bg="transparent")
par(mar=c(5,5,2,5))
plot(sites$site, sites$mean_bind_observed,type="l",xlab="site",ylab="mean ddG",ylim=c(6.5,-1),col="magenta",lwd=2)
par(new=T)
plot(sites$site, sites$mean_expr_observed,type="l",col="green",axes=F,xlab=NA,ylab=NA,ylim=c(-4,0.5))
axis(side=4)
mtext(side=4, line=3, "mean change in expression (ML meanF)")
points(sites$site, rep(0.35,nsites),pch=as.character(sites$AA),cex=0.25)
#dev.off()

#pdf(file="results/analyze_mut_effects/line-plot_mean-ddG-observed_and_mean-expr-latent_per_site_mixed-scales.pdf",width=9,height=3.5,bg="transparent")
par(mar=c(5,5,2,5))
plot(sites$site, sites$mean_bind_observed,type="l",xlab="site",ylab="mean ddG",ylim=c(6.5,-1),col="magenta",lwd=2)
par(new=T)
plot(sites$site, sites$mean_expr,type="l",col="green",axes=F,xlab=NA,ylab=NA,ylim=c(-4,0.5))
axis(side=4)
mtext(side=4, line=3, "mean change in expression (latent scale)")
points(sites$site, rep(0.35,nsites),pch=as.character(sites$AA),cex=0.25)
#dev.off()


#pairwise epistatic deviations -- for double mutant barcodes, compare actual measurement to that predicted from the component beta coefficients. (And plot versus 3D distance of residue pair)
#get bc_expr data for just double mutants
bc_expr_dbl <- bc_expr[n_aa_substitutions==2 & variant_class == ">1 nonsynonymous",]
bc_expr_dbl[,mutation1 := aa_subs_list[[1]][1],by=.(library,barcode)]
bc_expr_dbl[,mutation2 := aa_subs_list[[1]][2],by=.(library,barcode)]
bc_expr_dbl[,site1 := as.numeric(paste(unlist(strsplit(mutation1,split=""))[-c(1,length(unlist(strsplit(mutation1,split=""))))],collapse="")),by=.(library,barcode)]
bc_expr_dbl[,site2 := as.numeric(paste(unlist(strsplit(mutation2,split=""))[-c(1,length(unlist(strsplit(mutation2,split=""))))],collapse="")),by=.(library,barcode)]

#calculate actual dE from ML_meanF versus average wildtype value
bc_expr_dbl[library=="libA2", deltaE := ML_meanF - 9.490705]
bc_expr_dbl[library=="libB2", deltaE := ML_meanF - 9.427806]

#calculate expected dE from the expr_latent effects of the two component single mutations
bc_expr_dbl[,deltaE_mut1 := betas[mutation==mutation1,expr_latent],by=.(library,barcode)]
bc_expr_dbl[,deltaE_mut2 := betas[mutation==mutation2,expr_latent],by=.(library,barcode)]
bc_expr_dbl[,deltaE_predicted := deltaE_mut1 + deltaE_mut2,by=.(library,barcode)]

#how does this correlate with the latent phenotype estimated for the double mutant, and the "observed" phenotype?
plot(bc_expr_dbl$deltaE_predicted,bc_expr_dbl$latent_phenotype_Cauchy_A,pch=16,col="#00000033")
plot(bc_expr_dbl$deltaE_predicted,bc_expr_dbl$predicted_phenotype_Cauchy_A,pch=16,col="#00000033")
#correlates with latent, as it should (not 1:1 because these are from one replicate versus the other, betas are averaged)

#distribution of epistasis
plot(bc_expr_dbl$deltaE_predicted,bc_expr_dbl$deltaE,pch=16,col="#00000015")
#can see lots of weirdness with the censoring of the observable deltaE range compared to the predicted (latent scale)

plot(bc_expr_dbl$latent_phenotype_Cauchy_A, bc_expr_dbl$deltaE,pch=16,col="#00000015")
#actually, confused what to do here. Do with binding on observed, scale, seems more clear waht to do
#add residue pair distances to df

#plot deviance versus residue pair distance


#get bc_expr data for just double mutants
bc_bind_dbl <- bc_bind[n_aa_substitutions==2 & variant_class == ">1 nonsynonymous",]
bc_bind_dbl[,mutation1 := aa_subs_list[[1]][1],by=.(library,barcode)]
bc_bind_dbl[,mutation2 := aa_subs_list[[1]][2],by=.(library,barcode)]
bc_bind_dbl[,site1 := as.numeric(paste(unlist(strsplit(mutation1,split=""))[-c(1,length(unlist(strsplit(mutation1,split=""))))],collapse="")),by=.(library,barcode)]
bc_bind_dbl[,site2 := as.numeric(paste(unlist(strsplit(mutation2,split=""))[-c(1,length(unlist(strsplit(mutation2,split=""))))],collapse="")),by=.(library,barcode)]

#calculate expected ddG from the bind_observed effects of the two component single mutations
bc_bind_dbl[,ddG_mut1 := betas[mutation==mutation1,bind_observed],by=.(library,barcode)]
bc_bind_dbl[,ddG_mut2 := betas[mutation==mutation2,bind_observed],by=.(library,barcode)]
bc_bind_dbl[,ddG_predicted := ddG_mut1 + ddG_mut2,by=.(library,barcode)]

#how does this correlate with the latent phenotype estimated for the double mutant, and the "observed" phenotype?
plot(bc_bind_dbl$ddG_predicted,bc_bind_dbl$latent_phenotype_Cauchy_A,pch=16,col="#00000033")
plot(bc_bind_dbl$ddG_predicted,bc_bind_dbl$predicted_phenotype_Cauchy_A,pch=16,col="#00000033")

#distribution of epistasis
plot(bc_bind_dbl$ddG_predicted,bc_bind_dbl$ddG,pch=16,col="#00000015")
abline(h=7,lty=2);abline(v=7,lty=2)
abline(a=0,b=1,lty=2,col="red")

#calculate dddG (difference between predicted and actual ddG -- so, negative means double mutant performed better than expected from component singles)
bc_bind_dbl[,dddG := ddG - ddG_predicted,by=.(library,barcode)]
#probably will want to censor some of these downstream to make sure either predicted or observed ddGs aren't >6 or 8 kcal/mol (upper detectable limit)
hist(bc_bind_dbl$dddG)

#add residue pair distances to df
bc_bind_dbl$CA_dist <- sapply(1:nrow(bc_bind_dbl), function(x) return(pdb_dists[site1==bc_bind_dbl[x,site1] & site2==bc_bind_dbl[x,site2],CA_dist]))
bc_bind_dbl$min_dist <- sapply(1:nrow(bc_bind_dbl), function(x) return(pdb_dists[site1==bc_bind_dbl[x,site1] & site2==bc_bind_dbl[x,site2],min_dist]))

#plot deviance versus residue pair distance
plot(bc_bind_dbl$min_dist, bc_bind_dbl$dddG, pch=16,col="#00000020")
plot(bc_bind_dbl[ddG_predicted < 7 & ddG < 7, min_dist], bc_bind_dbl[ddG_predicted < 7 & ddG < 7, dddG], pch=16,col="#00000010")
plot(bc_bind_dbl[ddG_predicted < 7 | ddG < 7, min_dist], bc_bind_dbl[ddG_predicted < 7 | ddG < 7, dddG], pch=16,col="#00000010")


bc_bind_dbl[site1==35 & site2==109,]

#save Rsession (grabnode limits session to being open for 7 days)
save.image(file="results/analyze_mut_effects/2020_03_25_session-reload.Rdata")
#load(file="results/analyze_mut_effects/2020_03_25_session-reload.Rdata")

#save betas dataframe as csv
write.csv(betas, file="results/analyze_mut_effects/betas.csv")
#write sites dataframe as csv (leaving out some columns)
write.csv(sites[,c("site","AA","kabat_site","kabat_chain","kabat_annotation","pdb_site","pdb_chain","mean_bind","mean_bind_observed","mean_expr","mean_expr_observed","SSE","RSA_bound","RSA_unbound","naive_AA","minVRC01_AA")], file="results/analyze_mut_effects/sites.csv")

####################################
#older analysis used latent/predicted phenotypes to devise a 'deviance'?

#repeat for binding data
bc_bind[,n_mut := length(strsplit(aa_substitutions,split=" ")[[1]]),by=barcode]
bc_bind_dbl <- bc_bind[n_mut==2 & variant_class == ">1 nonsynonymous",]
for(i in 1:nrow(bc_bind_dbl)){
  subs <- bc_bind_dbl[i,aa_substitutions]
  subs <- unlist(strsplit(subs,split=" "))
  mut1 <- unlist(strsplit(subs[1],split=""))
  bc_bind_dbl[i,"site1"] <- as.numeric(paste(mut1[c(-1,-length(mut1))],collapse=""))
  mut2 <- unlist(strsplit(subs[2],split=""))
  bc_bind_dbl[i,"site2"] <- as.numeric(paste(mut2[c(-1,-length(mut2))],collapse=""))
  bc_bind_dbl[i,"wt1"] <- as.character(mut1[1])
  bc_bind_dbl[i,"wt2"] <- as.character(mut2[1])
  bc_bind_dbl[i,"mut1"] <- as.character(mut1[length(mut1)])
  bc_bind_dbl[i,"mut2"] <- as.character(mut2[length(mut2)])
}

plot(c(bc_bind_dbl[library=="libA2",predicted_phenotype_Cauchy_A],bc_bind_dbl[library=="libB2",predicted_phenotype_Cauchy_B]),c(bc_bind_dbl[library=="libA2",ddG],bc_bind_dbl[library=="libB2",ddG]),pch=16,col="#00000010")

hist(c(bc_bind_dbl[library=="libA2",ddG]-bc_bind_dbl[library=="libA2",predicted_phenotype_Cauchy_A],bc_bind_dbl[library=="libB2",ddG]-bc_bind_dbl[library=="libB2",predicted_phenotype_Cauchy_B]),col="#ff000025",breaks=20)
hist(c(bc_bind[variant_class %in% c("synonymous","wildtype") & library=="libA2",ddG]-bc_bind[variant_class %in% c("synonymous","wildtype") & library=="libA2",predicted_phenotype_Cauchy_A],bc_bind[variant_class %in% c("synonymous","wildtype") & library=="libB2",ddG]-bc_bind[variant_class %in% c("synonymous","wildtype") & library=="libB2",predicted_phenotype_Cauchy_B]),add=T,col="#0000ff25",breaks=20)  
hist(c(bc_bind[variant_class %in% c("1 nonsynonymous") & library=="libA2",ddG]-bc_bind[variant_class %in% c("1 nonsynonymous") & library=="libA2",predicted_phenotype_Cauchy_A],bc_bind[variant_class %in% c("1 nonsynonymous") & library=="libB2",ddG]-bc_bind[variant_class %in% c("1 nonsynonymous") & library=="libB2",predicted_phenotype_Cauchy_B]),add=T,col="#0000ff25",breaks=20)  
#doubles seem to have excess deviance relative to singles/WT

#add residue pair distances to df
bc_bind_dbl$CA_dist <- sapply(1:nrow(bc_bind_dbl), function(x) return(pdb_dists[site1==bc_bind_dbl[x,site1] & site2==bc_bind_dbl[x,site2],CA_dist]))
bc_bind_dbl$min_dist <- sapply(1:nrow(bc_bind_dbl), function(x) return(pdb_dists[site1==bc_bind_dbl[x,site1] & site2==bc_bind_dbl[x,site2],min_dist]))
bc_bind_dbl[library=="libA2", deviance := ddG - predicted_phenotype_Cauchy_A]
bc_bind_dbl[library=="libB2", deviance := ddG - predicted_phenotype_Cauchy_B]
bc_bind_dbl[, abs_deviance := abs(deviance)]

plot(bc_bind_dbl$CA_dist, bc_bind_dbl$deviance,pch=16,col="#00000020",xlab="CA distance",ylab="dddG (kcal/mol)")
plot(bc_bind_dbl$min_dist, bc_bind_dbl$deviance,pch=16,col="#00000010",xlab="residue pair distance",ylab="dddG (measured - predicted binding, kcal/mol)",main="binding phenotype")
points(bc_bind_dbl[site1==35 & site2==109,min_dist],bc_bind_dbl[site1==35 & site2==109,deviance],col="green")


#just the "negative epistasis" (actual measured ddG is better than predicted)
plot(bc_bind_dbl[deviance>0, min_dist],bc_bind_dbl[deviance>0, deviance],pch=16,col="#00000020")
bc_bind_dbl[deviance > 6,]

#just the "positive epistasis" (actual measured ddG is better than predicted)
plot(bc_bind_dbl[deviance<0, min_dist],bc_bind_dbl[deviance<0, deviance],pch=16,col="#00000020")
bc_bind_dbl[deviance < -6,]
bc_bind_dbl[deviance < -2 & min_dist < 10,]

#look at the epistatic interactions I want to validate -- do I have dbl mutant barcodes between sites 35/109, and 37/232?
bc_bind_dbl[site1==35 & site2==109,]
bc_expr_dbl[site1==35 & site2==109,]
bc_bind_dbl[site1==37 & site2==232 & !is.na(ddG),]
bc_expr_dbl[site1==37 & site2==232 & !is.na(ML_meanF),]

betas[mutation=="F232V",.(bind_observed,expr_observed)]

#test: for positive and all deviations -- variance versus distance -- expect higher scatter at low distance indicative of epistasis



#2d heatmap (position 1 by position 2, plot mean abs deviance at site-pairs on upper left, distance on lower right -- can I see "regions" of epistasis?)


#translating mutational effects onto a "preferences" scale for ExpCM analyses



