Structure function analysis of mutational effects
================
Tyler Starr
5/11/2020

This notebook analyzes our single mutant effects on binding and
expression in light of structural features of the RBD.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","bio3d","gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

#read in config file
config <- read_yaml("config.yaml")

#read in file giving concordance between RBD numbering and SARS-CoV-2 Spike numbering
RBD_sites <- data.table(read.csv(file="data/RBD_sites.csv",stringsAsFactors=F))

#make output directories
if(!file.exists(config$structure_function_dir)){
  dir.create(file.path(config$structure_function_dir))
}
#make output directory
if(!file.exists(config$dms_view_dir)){
  dir.create(file.path(config$dms_view_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 14.04.6 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /app/easybuild/software/OpenBLAS/0.2.18-GCC-5.4.0-2.26-LAPACK-3.6.1/lib/libopenblas_prescottp-r0.2.18.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3     bio3d_2.3-4       forcats_0.4.0    
    ##  [4] stringr_1.4.0     dplyr_0.8.3       purrr_0.3.2      
    ##  [7] readr_1.3.1       tidyr_0.8.3       tibble_2.1.3     
    ## [10] ggplot2_3.2.0     tidyverse_1.2.1   data.table_1.12.2
    ## [13] yaml_2.2.0        knitr_1.23       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.1       cellranger_1.1.0 pillar_1.4.2     compiler_3.6.1  
    ##  [5] tools_3.6.1      digest_0.6.20    lubridate_1.7.4  jsonlite_1.6    
    ##  [9] evaluate_0.14    nlme_3.1-140     gtable_0.3.0     lattice_0.20-38 
    ## [13] pkgconfig_2.0.2  rlang_0.4.0      cli_1.1.0        rstudioapi_0.10 
    ## [17] parallel_3.6.1   haven_2.1.1      xfun_0.7         withr_2.1.2     
    ## [21] xml2_1.2.0       httr_1.4.0       hms_0.4.2        generics_0.0.2  
    ## [25] grid_3.6.1       tidyselect_0.2.5 glue_1.3.1       R6_2.4.0        
    ## [29] readxl_1.3.1     rmarkdown_1.13   modelr_0.1.4     magrittr_1.5    
    ## [33] backports_1.1.4  scales_1.0.0     htmltools_0.3.6  rvest_0.3.4     
    ## [37] assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3    lazyeval_0.2.2  
    ## [41] munsell_0.5.0    broom_0.5.2      crayon_1.3.4

## Setup

Read in tables of variant effects on binding and expression, for single
mutations to the SARS-CoV-2 RBD and for a panel of homolog RBDs from the
sarbecovirus clade.

``` r
homologs <- data.table(read.csv(file=config$homolog_effects_file,stringsAsFactors = F))
mutants <- data.table(read.csv(file=config$single_mut_effects_file,stringsAsFactors = F))

#rename mutants site indices to prevent shared names with RBD_sites, simplifying some downstream calculations that cross-index these tables
setnames(mutants, "site_RBD", "RBD_site");setnames(mutants, "site_SARS2", "SARS2_site")

#add color column to homologs, by clade
homologs$clade_color <- as.character(NA); homologs[clade=="Clade 1",clade_color := "#EF4136"]; homologs[clade=="Clade 2",clade_color := "#009444"]; homologs[clade=="Clade 3",clade_color := "#EE2A7B"]; homologs[clade=="SARS-CoV-2",clade_color := "#2E3192"]
```

## General structural constraints on RBD affinity and stability

Let’s investigate how the RBD structure influences mutational effects on
expression and binding. First, we compute the *mean* effect of mutations
on binding and expression for each RBD site, as well as the best (max)
and worst (min) mutational effects on these two measurements (excluding
nonsense and synonymous mutants).

``` r
RBD_sites[,mean_bind := mean(mutants[SARS2_site==site_SARS2 & wildtype != mutant & mutant != "*",bind_avg],na.rm=T),by=site_SARS2]
RBD_sites[,max_bind := max(mutants[SARS2_site==site_SARS2 & wildtype != mutant & mutant != "*",bind_avg],na.rm=T),by=site_SARS2]
RBD_sites[,min_bind := min(mutants[SARS2_site==site_SARS2 & wildtype != mutant & mutant != "*",bind_avg],na.rm=T),by=site_SARS2]

RBD_sites[,mean_expr := mean(mutants[SARS2_site==site_SARS2 & wildtype != mutant & mutant != "*",expr_avg],na.rm=T),by=site_SARS2]
RBD_sites[,max_expr := max(mutants[SARS2_site==site_SARS2 & wildtype != mutant & mutant != "*",expr_avg],na.rm=T),by=site_SARS2]
RBD_sites[,min_expr := min(mutants[SARS2_site==site_SARS2 & wildtype != mutant & mutant != "*",expr_avg],na.rm=T),by=site_SARS2]
```

First, let’s see how mutational effects on binding and expression
correlate at the level of individual mutations and at the site-level
mean effects of mutation. We can see below that for many mutations and
sites, mutational effects on expression and binding are related,
indicating that stability is a generic constraint on mutational effects
on ACE2-binding. However, there are a handful of positions that deviate
from this pattern, being tolerant to mutation with respect to expression
despite being quite sensitive to mutation with respect to binding.
Below, we output the sites of these mutations, which we can visualize on
the ACE2-bound RBD structure using `dms-view`, linked [here for
sites](https://dms-view.github.io/?pdb-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2F6m0j.pdb&markdown-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2FBloomLab_rbd.md&data-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2Fresults%2FBloomLab2020_rbd.csv&condition=natural+frequencies&site_metric=site_entropy&mutation_metric=mut_frequency&selected_sites=447%2C449%2C456%2C473%2C476%2C487%2C489%2C496%2C500%2C502%2C505)
and [here for
mutations](https://dms-view.github.io/?pdb-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2F6m0j.pdb&markdown-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2FBloomLab_rbd.md&data-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2Fresults%2FBloomLab2020_rbd.csv&condition=natural+frequencies&site_metric=site_entropy&mutation_metric=mut_frequency&selected_sites=443%2C455%2C456%2C475%2C487%2C489%2C496%2C498%2C500%2C501%2C502%2C505)
that fall in this binding-specific defective class. We can see that
these sites exhibiting binding-specific mutational sensitivity are at
the ACE2-contact interface, or in the case of one mutation (S443N),
perhaps second shell posititions that are still ACE2-proximal. This is
consistent with these positions having binding constraints independent
of stability because of their direct interaction with ACE2.

``` r
par(mfrow=c(1,2))
plot(RBD_sites$mean_expr,RBD_sites$mean_bind,pch=19,col="#00000050",xlab="mean mutational effect on expression",ylab="mean mutational effect on binding",main="binding versus expression effects,\naverage per site")

plot(mutants$expr_avg,mutants$bind_avg,pch=19,col="#00000050",xlab="mutational effect on expression",ylab="mutational effect on binding",main="binding versus expression effects,\nper mutant")
```

<img src="structure_function_files/figure-gfm/correlation_mut_expression_binding-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structure_function_dir,"/correlation_expression_v_binding.pdf",sep="")))

#output sites and mutations with seemingly binding-specific detrimental effects
RBD_sites[mean_expr > -1 & mean_bind < -1,site_SARS2]
```

    ##  [1] 443 447 449 456 473 475 476 487 489 496 500 502 505

``` r
mutants[expr_avg > -0.5 & bind_avg < -2,mutation]
```

    ##  [1] "S443N" "L455D" "L455E" "F456A" "F456E" "F456G" "F456N" "F456Q"
    ##  [9] "F456R" "F456S" "A475D" "N487C" "N487E" "N487K" "N487L" "N487M"
    ## [17] "N487Q" "N487R" "Y489A" "Y489C" "Y489E" "Y489I" "Y489K" "Y489L"
    ## [25] "Y489N" "Y489P" "Y489Q" "Y489R" "Y489S" "Y489T" "Y489V" "G496D"
    ## [33] "G496E" "G496F" "Q498K" "T500I" "N501D" "N501K" "N501R" "G502A"
    ## [41] "G502C" "G502D" "G502E" "G502F" "G502H" "G502I" "G502K" "G502L"
    ## [49] "G502M" "G502N" "G502P" "G502Q" "G502R" "G502S" "G502T" "G502V"
    ## [57] "G502W" "G502Y" "Y505A" "Y505C" "Y505D" "Y505E" "Y505G" "Y505I"
    ## [65] "Y505K" "Y505L" "Y505M" "Y505Q" "Y505R" "Y505S" "Y505T" "Y505V"

The relative solvent accessibility (RSA) of an amino acid residue is
known to be a dominant factor influencing its tolerance to mutation.
Let’s see how RSA of a position is related to its mutational
sensitivity for binding. We use RSA in two different structural contexts
– the free RBD structure, and the RBD structure when complexed with
ACE2. We can see that mutational sensitivity of a position with respect
to binding is better described by RSA in the ACE2-bound RBD complex.
This explains our observation above – some positions are mutationally
sensitive, because they are buried in the isolated RBD, and so mutations
destabilize the core fold and thereby hamper binding, while others are
mutationally sensitive not because they are buried in the core fold and
are sensitive to destabilizing mutations, but because they are buried at
the ACE2 interface. These two factors combine to explain our overall
observed patterns of mutational tolerance.

<img src="structure_function_files/figure-gfm/mean_mut_effect_versus_RSA-1.png" style="display: block; margin: auto;" />

To further visualize site-wise mutational sensitivity on the 3D
structure, let’s output `.pdb` files for the ACE2-bound RBD structure in
which we replace the B factor column with metrics of interest for each
site: 1) the mean effect of mutation on binding, 2) the mean effect of
mutation on expression, 3) the max effect of any mutation on binding,
and 4) the max effect of any mutation on expression. In PyMol, we can
then visualize spheres at each Calpha, colored by spectrum from low
(yellow) to high (blue) for each metric by manually executing the
following commands in a PyMol session in which one of the output `pdb`
files is opened:

    hide all; show cartoon
    color warmpink, chain A; color gray80, chain E
    set sphere_scale, 0.6
    create RBD_CA, chain E and name ca
    hide cartoon, RBD_CA; show spheres, RBD_CA
    spectrum b, yellow green blue, RBD_CA

``` r
pdb <- read.pdb(file="data/structures/ACE2-bound/6m0j.pdb")
```

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
#color by mean effect on binding
pdb_mean_bind <- pdb
pdb_mean_bind$atom$b <- NA
for(i in 1:nrow(pdb_mean_bind$atom)){
  res <- pdb_mean_bind$atom$resno[i]
  chain <- pdb_mean_bind$atom$chain[i]
  mean_bind <- RBD_sites[site_SARS2==res & chain_6M0J==chain & !is.na(chain_6M0J),mean_bind]
  if(length(mean_bind)>0){pdb_mean_bind$atom$b[i] <- mean_bind}else{pdb_mean_bind$atom$b[i] <- 0}
}
#save pdb
write.pdb(pdb=pdb_mean_bind,file=paste(config$structure_function_dir,"/6m0j_b-factor-mean-bind.pdb",sep=""), b = pdb_mean_bind$atom$b)

#color by max effect on binding
pdb_max_bind <- pdb
pdb_max_bind$atom$b <- NA
for(i in 1:nrow(pdb_max_bind$atom)){
  res <- pdb_max_bind$atom$resno[i]
  chain <- pdb_max_bind$atom$chain[i]
  max_bind <- RBD_sites[site_SARS2==res & chain_6M0J==chain & !is.na(chain_6M0J),max_bind]
  if(length(max_bind)>0){pdb_max_bind$atom$b[i] <- max_bind}else{pdb_max_bind$atom$b[i] <- 0}
}
#save pdb
write.pdb(pdb=pdb_max_bind,file=paste(config$structure_function_dir,"/6m0j_b-factor-max-bind.pdb",sep=""), b = pdb_max_bind$atom$b)

#color by mean effect on expression
pdb_mean_expr <- pdb
pdb_mean_expr$atom$b <- NA
for(i in 1:nrow(pdb_mean_expr$atom)){
  res <- pdb_mean_expr$atom$resno[i]
  chain <- pdb_mean_expr$atom$chain[i]
  mean_expr <- RBD_sites[site_SARS2==res & chain_6M0J==chain & !is.na(chain_6M0J),mean_expr]
  if(length(mean_expr)>0){pdb_mean_expr$atom$b[i] <- mean_expr}else{pdb_mean_expr$atom$b[i] <- 0}
}
#save pdb
write.pdb(pdb=pdb_mean_expr,file=paste(config$structure_function_dir,"/6m0j_b-factor-mean-expr.pdb",sep=""), b = pdb_mean_expr$atom$b)

#color by max effect on expression
pdb_max_expr <- pdb
pdb_max_expr$atom$b <- NA
for(i in 1:nrow(pdb_max_expr$atom)){
  res <- pdb_max_expr$atom$resno[i]
  chain <- pdb_max_expr$atom$chain[i]
  max_expr <- RBD_sites[site_SARS2==res & chain_6M0J==chain & !is.na(chain_6M0J),max_expr]
  if(length(max_expr)>0){pdb_max_expr$atom$b[i] <- max_expr}else{pdb_max_expr$atom$b[i] <- 0}
}
#save pdb
write.pdb(pdb=pdb_max_expr,file=paste(config$structure_function_dir,"/6m0j_b-factor-max-expr.pdb",sep=""), b = pdb_max_expr$atom$b)
```

Does mutational tolerance with respect to binding and/or expression
differ systematically between positions in the core-RBD versus the RBM
loops? (Hypothesis is: expression is more constrained in the core RBD,
while binding is more constrained in the RBM) – analysis upcoming

## Distribution of functional effects of mutatioin

Let’s look at the distribution of single-mutant effects on our two
phenotypes, and compare the fraction of mutations that are within the
window defined by known functional RBD homologs for these two
phenotypes.

For the binding plot on the left, the intermediate blue point on the
x-scale is RaTG13, which *can* promote huACE2-mediated cell entry in an
in vitro cellular infection assay (though less efficiently than
SARS-CoV-2) according to [Shang et
al. 2020](https://www.nature.com/articles/s41586-020-2179-y/figures/3),
though whether this is sufficient to enable efficient viral replication
in more complex models is uncertain. For the cluster of homologs near 0,
the farthest-left point is LYRa11, which according to [Letko et
al. 2020](https://www.nature.com/articles/s41564-020-0688-y/figures/1)
can also promote huACE2-mediated cellular entry, though less efficiently
than SARS-CoV-1 and other bat CoV isolates such as WIV1/16 (identical
RBDs). Therefore, these two points define a window of affinites that can
at least support in vitro cellular infection – but in reality, the
window of possible “neutrality” with regards to actual human infectivity
is perhaps better set by the remaining four points with delta
log<sub>10</sub>(*K*<sub>A,app</sub>) values \~ 0 – these four points
are the RBDs from SARS-CoV-1, WIV1/16, SARS-CoV-2, and GD-Pangolin RBDs,
in that rank-order. Taken together, we identify 1732 single mutants
(45.55%) whose affinity effects are within a neutral window that
potentially enables huACE2-mediated infectivity (SARS-CoV-1 cutoff), and
3188 single mutants (83.85%) whose affinity is potentially sufficient to
enable huACE2-mediated in vitro cellular infectivity (RaTG13 cutoff).
Taken together, this suggests a quite large sequence space of RBD
diversity that is consistent with huACE2 binding and entry. (Though, of
course, our isolated RBD-only affinity measurements may have more
complex constraints in the context of full-Spike trimer.)

For expression, our range set by homologs may be improperly aligned to
the SARS-CoV-2 range – not least, because all other RBD homologs were
observed to have very slightly higher expression than SARS-CoV-2, so
direct interpretaiton of DFE versus homologs is a bit more complicated
(could be artefactual because the SARS-CoV-2 variants went through the
mutagenesis protocol, though we eliminated the clear outliers that
likely bear unseen mutations outside the RBD sequence. SARS-CoV-2 is
still posited to have lower expression than the other homologs, but this
elimination of low-expression wildtype/synonymous variants may have
further distanced it from the true “library-average” wildtype expression
and thereby make more of the single mutants look deleterious? Tough
problem to solve…)

<img src="structure_function_files/figure-gfm/DFE_bind-1.png" style="display: block; margin: auto;" />

## Exploratory heatmaps

Next, let’s make heatmaps of per-amino acid mutational effects on
binding and expression. We first make these heatmaps for all mutations
at all sites, colored by delta log<sub>10</sub>(*K*<sub>A,app</sub>). (I
would eventually like to flesh out these heatmaps by providing
additional indicator variables as ‘heatmap’ style rows on the top –
things could include a color scale for RSA, conservation, an asterisk or
something to indicate contact residues, indicators for NLGS or
disulfides, an indicator for strucutral contacts in full length Spike
trimer, etc. I also want to think about whether there are better color
scales to use, including the “divergence” in scale between the blue
(goes from white to bluest in a 0.5-unit scale for binding, 1-unit scale
for expression), versus red (0 to reddest in a 5-unit scale). Happy for
additional things to think about to refine these plots):

<img src="structure_function_files/figure-gfm/heatmap_binding_all-1.png" style="display: block; margin: auto;" />
And next, the same heat map, colored by delta mean fluorescence
(expression). I altered the scale such that anything less than -4 gets
the darkest red color – the lowest missense mutant score is -4.07, only
nonsense mutants push the scale down to -5, so this helps the relative
red scale be better calibrated between the binding and expression
measurements w.r.t. single missense mutations. (once again, happy for
thoughts.):

<img src="structure_function_files/figure-gfm/heatmap_expression_all-1.png" style="display: block; margin: auto;" />

To make more sense of these heatmaps, we now zoom in and group residues
in special structural and functional classes. First, let’s visualize
heatmaps of paired cysteines that form disulfides. The following
heatmaps reorder sites by cysteine pair. (I would love to add a line at
the top grouping the paired cysteines, or create a small gap between
each pair of columns, but I am gg-inept and so this will do for now.)

We can see that cysteines within a disulfide pair have similiar
sensitivities to mutation and even similar biochemical preferences, but
there are differing patterns of sensitivity for binding and expression.
Disulfide pair 4 is the most sensitive to mutation with respect to
binding, whereas it only is moderately important for expression – this
is the disulfide pair within the RBM loop that stabilizes regions of the
ACE2 interface, consistent with its exacerbated importance for binding
versus expression. The remaining three disulfides are in the core RBD
domain, and show similar sensitivities for expression and binding – pair
3 is tolerant to mutation (and may even support replacement with polar
amino acids), pair 1 is moderately sensitive, and pair 2 is most
sensitive – though hydrophobic mutations have noticeably decreased
deleterious effect for binding, which is interesting to see.

This trend is consistent with a series of Cys-\>Ala mutations made in
SARS-CoV-2 RBD in [Wong et
al. JBC 2004](https://www.jbc.org/content/279/5/3197.short) (I report
the numbers converted to SARS-CoV-2 numbering): mutations at SARS-CoV-1
RBD positions 379 and 432 severely impaired expression, mutations at
361, 480, and 488 had only mild effects on expression but strongly
decreased ACE2-binding, and mutations at 336 and 391 had little effect
on expression or binding (this is in 331-524 RBD construct, so lacking
cysteine 525).

<img src="structure_function_files/figure-gfm/heatmap_bind_expr_disulfide-1.png" style="display: block; margin: auto;" />

Next, let’s look at the N-linked glycosylation sites. In particular, we
look at mutational effects at NLGS motifs consisting of N331 and T333,
and N343 and T345, as well as N370 and A372, which in SARS-CoV-1 is an
NST NLGS motif. A prior paper, [Chen et
al. 2014](https://www.tandfonline.com/doi/full/10.4161/hv.27464),
showed that in SARS-CoV-1 RBD produced in Pichia yeast, yields were
lowered when progressively knocking out each of these glycans,
suggesting they are important for expression. (They didn’t do individual
knockouts, it seems.) Finally, on the righthand side, we look at the
effects of mutations i+2 from surface asparagines, to see whether
potential glycan knockins (mutations to S or T) have any interesting
effects.

As expected the two NLGS in the SARS-CoV-2 RBD are important for
stability. We can see that mutations to the focal asparagine or the N+2
threonine are uniiversally deleterious with regards to expression, with
the exception of T-\>S mutations which retain the NLGS. Overall, the
N343 glycan appears more important to RBD stability than the N331
glycan. Re-introducing the SARS-CoV-1 NLGS at site 370 has just a mildly
deleterious effect on expression. None of these three glycan mutants
have major impacts on binding, consistent with their distance from the
ACE2-interface.

At other non-glycosylated asparagines, I don’t see many strong effects
when introducing putative NLGS motifs with an i+2 mutation to S/T.
Glycosylation of N354, N360, N448, N450, N481, and N487 may have mildly
beneficial effects on expression, and glycosylation may be \~neutral or
tolerated with only minor detrimental effects with regards to expression
at sites 388, 394, 437, 439, 460, and 501. These sites are [visualized
on the RBD
structure](https://dms-view.github.io/?pdb-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2F6m0j.pdb&markdown-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2FBloomLab_rbd.md&data-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2Fresults%2FBloomLab2020_rbd.csv&condition=natural+frequencies&site_metric=site_entropy&mutation_metric=mut_frequency&selected_sites=354%2C360%2C388%2C394%2C437%2C439%2C448%2C450%2C460%2C481%2C487%2C501),
illustrating that quite a few surfaces on the RBD could potentially be
masked with glycan introductions, including different portions of the
RBM if wanting to target different epitopes with vaccine immunogens or
probes for isolating mAbs for particular epitope surfaces.

In contrast to novel glycans that could be tolerated in an engineered
RBD construct with respect to stability/expression, we can see that
introduction of an NLGS to N501 with mutation of site 503 to T or S may
have a detrimental effect on biinding, consistent with this interface
residue making key ACE2-interactions (site 501 is one of the “key
contacts” from the SARS-CoV-1 literature). Knockin of a possible N439
NLGS (via T/S mutations at site 441) also has a mild deleterious effect
on affinity, specific to the T/S mutations at this i+2 position; site
439 is sort of “second shell” from the ACE2 interface, but close enough
to imagine a glycan could impact affinity. There are other positions
where i+2 S/T mutations have large deleterious effects on binding
affinity (putatively glycosylating N422 (sort of buried, so might not
actually be glycosylated) and N487 (interface\!)), but other amino acid
mutations at these positions also have strong deleterious effects, so it
is hard to know whether the effect of the T/S mutants stems from the
addition of a glycan, or loss of the wildtype amino acid at this i+2
residue independent of the glycan effect.

<img src="structure_function_files/figure-gfm/heatmap_bind_expr_NLGS-1.png" style="display: block; margin: auto;" />

Here are the heat maps, zoomed in on residues in the “Receptor Binding
Motif”, the section of the RBM that extends out from the core alpha+beta
fold and contains the ACE2-contact residues. This motif is contiguous in
space, so there’s really no point to this heat map beyond the full-RBD
map shown above, once I figure out how to add e.g. a line at the top
across these positions to incidicate RBM. I have some interpretation of
expression/binding tradeoff from this heatmap, but this point is better
made in the subsequent figure focusing only on contact positions, so I
will elaborate it there. Haven’t looked into this RBM map in too much
more detail beyond that, so there might be something else interesting,
or this figure could definitely get the axe as it probably doesn’t add
much.

<img src="structure_function_files/figure-gfm/heatmap_RBM-1.png" style="display: block; margin: auto;" />

Next, let’s look at expression and binding effects for annotated contact
residues. Below, we are looking at the 19 residues that form ACE2
contacts in the SARS-CoV-2 (6m0j) or SARS-CoV-1 (2ajf) ACE2-bound RBD
crystal structures, where we annotated a residue as a contact if it
contains non-hydrogen atoms within 4 Angstroms of ACE2. 14 residues were
annotated as contacts in both structures, 3 contacts (417, 446, 475) are
unique to the SARS-CoV-2 structure, and 2 (439, 503) are unique to
SARS-CoV-1. We indicate the wildtype SARS-CoV-2 and -1 amino acids with
an “x” and “o”, respectively.

Putting this reduced display of binding and expression next to each
other is really interesting, as it highlights some potential
binding-expression tradeoffs at positions 417, 449, 455, 486, 502, and
505, which are visualized in `dms-view`
[here](https://dms-view.github.io/?pdb-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2F6m0j.pdb&markdown-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2FBloomLab_rbd.md&data-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2Fresults%2FBloomLab2020_rbd.csv&condition=natural+frequencies&site_metric=site_entropy&mutation_metric=mut_frequency&selected_sites=417%2C449%2C455%2C486%2C502%2C505).

For four of these residues (449, 455, 486, 505), binding prefers the
wildtype hydrophobic state, though this exposed hydrophobic amino acid
is evidently detrimental to expression, as polar amino acid mutations
improve expression. Site 417 is the opposite, consistent with its more
buried position a bit further from ACE2 – for expression, keeping this
amino acid hydrophobic would be preferred, but mutations to K or R can
enhance affinity, presumably because their long side chains can snorkel
out and make contact with ACE2 interface. Finally, site 502 is a
glycine, which accomodates the ACE2 loop bearing the critical K31
“hotspot” residue, which would clash sterically if mutating site 502
to any amino acid with a side chain, even though most polar side chains
at this position would improve exression. These positions are generally
solvent-accessible even in the full-Spike trimer RBD “down”
conformation, but if we want to follow up on this we should really look
more into whether full Spike trimer would ameliorate any of these
potential detrimental expression effects of surface-exposed hydrophobic
residues.

Within this heatmap is specific information about the “key” contact
sites that are often discussed in the literature. These sites primarily
came from observations about amino acid changes relating civet- and
human-adapted SARS-CoV-1 sequences, and other experimental evolution and
structural studies in the SARS-CoV-1 side of the tree. These sites are
455, 486, 493, 494, 501. May be worth extra glances when looking at
these exploratory data, though other positions may be more important
here in the SARS-CoV-2 side of the tree.

<img src="structure_function_files/figure-gfm/heatmap_ACE2_contacts-1.png" style="display: block; margin: auto;" />

Next, are heatmaps for all positions that differ in amino acid identity
between SARS-CoV-1 and -2, with the same “x” and “o” indicators for the
SARS-CoV-2 and -1 wildtype state. The positions of all of the variable
amino acids between these virus RBDs is shown in `dms-view`
[here](https://dms-view.github.io/?pdb-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2F6m0j.pdb&markdown-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2FBloomLab_rbd.md&data-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2Fresults%2FBloomLab2020_rbd.csv&condition=natural+frequencies&site_metric=site_entropy&mutation_metric=mut_frequency&selected_sites=346%2C348%2C354%2C357%2C372%2C373%2C384%2C393%2C402%2C403%2C406%2C417%2C430%2C434%2C438%2C439%2C441%2C443%2C444%2C445%2C446%2C452%2C455%2C456%2C458%2C459%2C460%2C462%2C470%2C471%2C472%2C473%2C474%2C475%2C476%2C477%2C478%2C481%2C482%2C484%2C485%2C486%2C490%2C493%2C494%2C498%2C499%2C501%2C503%2C519%2C529).

Overall, there are not any *major* incompatibilities that emerge – that
is, no amino acids found in SARS-CoV-1 have huge negative effects if
introduced into the SARS-CoV-2 background – the only quite minor
affinity defects that are evident seem to be horizontal swap mutations
K417V, L455V, II472P, A475P, G476D, G482P, F486L, and S494D, along with
slight expression defects caused by S438T and N439R. These
negative-effect SARS-CoV-1 mutations are highlighted in `dms-view`
[here](https://dms-view.github.io/?pdb-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2F6m0j.pdb&markdown-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2FBloomLab_rbd.md&data-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2Fresults%2FBloomLab2020_rbd.csv&condition=natural+frequencies&site_metric=site_entropy&mutation_metric=mut_frequency&selected_sites=417%2C455%2C472%2C475%2C476%2C482%2C486%2C494),
where we can see that they cluster in two parts of the RBM – the lateral
loop containing the RBM disulfide, and the medial portion of the RBM
loop.

On the other hand, many mutations to SARS-CoV-1 amino acid identities
seem to have *positive* effects on binding and expression. Most notably
for expression, are mutations N354E, R357K (mild), K417V (discussed
above, also in context of binding/expression trade-off), and L452K, and
for binding, beneficial effects are also seen for N460K, Q498V (quite
substantial positive effect, direct interface residue), N501T (also
quite substantial, interface, site of key adaptation in SARS-CoV-1),
V503I, and H519N. These positions are highlighted on the structure
[here](https://dms-view.github.io/?pdb-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2F6m0j.pdb&markdown-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2FBloomLab_rbd.md&data-url=https%3A%2F%2Fraw.githubusercontent.com%2Fdms-view%2FSARS-CoV-2%2Fmaster%2Fdata%2FSpike%2FBloomLab2020%2Fresults%2FBloomLab2020_rbd.csv&condition=natural+frequencies&site_metric=site_entropy&mutation_metric=mut_frequency&selected_sites=354%2C357%2C417%2C452%2C460%2C498%2C501%2C503%2C519),
where we can see that many are distant from the ACE2 interface
(primarily the expression observations). Interestingly, the beneficial
*binding* SARS-CoV-1 residues are on the *other* lateral end of the RBM
compared to the detrimental swaps highlighted above. Therefore, it seems
like, of the interface differences between SARS-CoV-1 and -2, one “edge”
of the RBM saddle is better optimized for ACE2-binding in SARS-CoV-1,
while the middle and the other edge are more optimized in SARS-CoV-2. I
am actually surprised to see that there are seemingly so many
*beneficial* swaps to SARS-CoV-1 identities (should more rigorously
compare in number and magnitude with the *detrimental* swaps), in light
of the fact that SARS-CoV-2 has repeatedly been demonstrated to have
overall tighter ACE2-binding affinity than SARS-CoV-1 – and this
suggests that the SARS-CoV-2 RBD could bind ACE2 even more tightly if
there was pressure to do so (although we don’t actually know what the
impact of tighter affinity is with regards to transmissibility,
pathogenicity, tropism etc. which is a very important caveat.)

<img src="structure_function_files/figure-gfm/heatmap_SARS_CoV_2_1_diff-1.png" style="display: block; margin: auto;" />

Next, let’s look at the differences between RaTG13 and SARS-CoV-2 – in
contrast to SARS-CoV-1, which has just slightly weaker ACE2-binding
affinity, RaTG13 has affinity a couple orders of magnitude lower than
SARS-CoV-2 – so, comparing their states might identify some of the key
adaptations that occurred from their common ancestor on the SARS-CoV-2
lineage. However, the branch to RaTG13 is relatively long, so many of
the differences are more likely to be RaTG13-specific substitutions,
rather than a reflection of the ancestral state. This is supported when
also mapping the GD-Pangolin state at these positions, which more
frequently matches the SARS-CoV-2 identity. A more rigorous analysis
will need to do proper ancestral sequence reconstruction to polarize the
substitutions that occurred on the focal branch leading to SARS-CoV-2,
given the complicated phylogenetic relationships in the RBD in this
clade. For the “by-eye” parsimony reconstructions, in the heatmaps
below, “o” continues to mark the WT SARS-CoV-2 amino acid, “^”
represents RaTG13, and “\#” represents GD-Pangolin.

The main thing to note, is that site 501 seems to explain a large degree
of difference between RaTG13 and SARS-CoV-2 affinity – from the
GD-Pangolin sequence, it appears that the difference may be an N501D
substitution on the RaTG13 branch (so not part of the adaptation from
ancestor to SARS-CoV-2 – GX-Pangolin is a T here, so perhaps not helpful
in the by-eye polarization). Site 501 was a key site of adaptation
between civet- and human-adapted strains in SARS-CoV-1, so it’s
interesting to see it cropping up again here, regardless of the
polarization of the change (and different amino acid states than the S/T
SARS-CoV-1 stuff). Mutations to the RaTG13 amino acid are also mildly
deleterious w.r.t. binding for R403T (GD-Pangolin has R, GX-Pangolin K),
Y449F (GD-Pangolin has Y, GX-Pangolin Y), F486L (GD-Pangolin has F,
GX-Pangolin L…), and Y505H (GD-Pangolin has Y, GX-Pangolin Y). So, from
parsimony, it does seem most of these ‘deleterious’ amino acids are
*derived* substitutions on the RaTG13 lineage, instead of reflecting
‘adaptive’ changes on the SARS-CoV-2 ancestral lineage.

<img src="structure_function_files/figure-gfm/heatmap_SARS_CoV_2_RaTG13_diff-1.png" style="display: block; margin: auto;" />

Based on heatmap gazing along with structures and prior literature, I
would propose the following validation mutants. As you can see, several
positions are prioritized in these panels – sites 455, and 501 are
positions of interest from prior literature on SARS-CoV-1 adaptation;
site 502 is the second most constrained position w.r.t. to binding in
our dataset, with the most sensitive position (G431) being more
constrained by stability/expression affects than binding itself per se.
Site 498 exhibits lots of differences amongst the relevant strains I’ve
been looking at (SARS-CoV-2 versus -1, RaTG13, GD-Pangolin) and has
large variation in functional effects of mutation, and so is clearly a
position of interest for our data. All of these positions are in the RBM
and direct or near-direct ACE2 contact positions.

For yeast display validation of beneficial mutation effects and
expression/binding tradeoff, I would propose to validate the following
six mutations:

| mutation | RBD\_site | bind\_lib1 | bind\_lib2 | bind\_avg | expr\_lib1 | expr\_lib2 | expr\_avg | SARS1\_indicator | RaTG13\_indicator | GDPangolin\_indicator |
| :------- | --------: | ---------: | ---------: | --------: | ---------: | ---------: | --------: | :--------------- | :---------------- | :-------------------- |
| Q498H    |       168 |       0.30 |       0.31 |      0.30 |       0.15 |       0.16 |      0.16 |                  |                   | \#                    |
| Q498Y    |       168 |       0.01 |       0.30 |      0.16 |       0.19 |     \-0.07 |      0.06 | o                | ^                 |                       |
| N501D    |       171 |     \-2.40 |     \-2.44 |    \-2.42 |       0.09 |       0.07 |      0.08 |                  | ^                 |                       |
| N501F    |       171 |       0.22 |       0.36 |      0.29 |     \-0.03 |     \-0.16 |    \-0.10 |                  |                   |                       |
| N501T    |       171 |     \-0.14 |       0.33 |      0.10 |     \-0.35 |     \-0.16 |    \-0.25 | o                |                   |                       |
| G502D    |       172 |     \-4.29 |     \-4.03 |    \-4.16 |       0.62 |       0.66 |      0.64 |                  |                   |                       |

For pseudovirus growth assays, I would propose to validate the following
six mutations:

| mutation | RBD\_site | bind\_lib1 | bind\_lib2 | bind\_avg | expr\_lib1 | expr\_lib2 | expr\_avg | SARS1\_indicator | RaTG13\_indicator | GDPangolin\_indicator |
| :------- | --------: | ---------: | ---------: | --------: | ---------: | ---------: | --------: | :--------------- | :---------------- | :-------------------- |
| L455Y    |       125 |     \-1.48 |     \-1.52 |    \-1.50 |     \-0.16 |     \-0.03 |    \-0.10 | o                |                   |                       |
| A475P    |       145 |     \-1.52 |     \-1.72 |    \-1.62 |     \-1.44 |     \-1.34 |    \-1.39 | o                |                   |                       |
| Q498Y    |       168 |       0.01 |       0.30 |      0.16 |       0.19 |     \-0.07 |      0.06 | o                | ^                 |                       |
| N501D    |       171 |     \-2.40 |     \-2.44 |    \-2.42 |       0.09 |       0.07 |      0.08 |                  | ^                 |                       |
| N501Q    |       171 |     \-0.09 |     \-0.04 |    \-0.06 |       0.02 |       0.26 |      0.14 |                  |                   |                       |
| G502D    |       172 |     \-4.29 |     \-4.03 |    \-4.16 |       0.62 |       0.66 |      0.64 |                  |                   |                       |

## Output for dms-view visualization of mutational effects

Let’s output data in a format for import into `dms-view`. We want to
remove nonsense mutants, rename columns per the `dms-view` input column
specifications, and move our bind and expression effect measurements
into a long data frame format. We will start with our delta binding and
expression measurements relative to wildtype, and add a rescaled metric
that resembles our more traditional logo plot visualizations.

We retitle the columns already in our data table to the header names
specified in the `dms-view` input file guidelines. We then convert to
“long” data frame format, and incorporate the desired mutation- and
site-level mutation measurements. We treat binding and expression as
separate “conditions”, for which we will show site-wise metrics of mean,
max, and min effect of mutations at the site on the condition phenotype,
and the mutations will be illustrated as their delta-phenotype compared
to wildtype, where phenotype is log<sub>10</sub>(*K*<sub>A,app</sub>)
for binding, or mean fluorescence for expression. Finally, we output a
different mutation-level metric, which takes the exponentiated value of
each delta metric, and normalizes it by the sum of all exponentiated
values at that position. This should scale mutational effects from 0-1,
and displays logo plots in a more straightforward manner where large
letters indicate preferred amino acids, and small letters indicate
deleterious mutational effects.

``` r
dmsview_output_bind <- mutants[mutant != "*",.(site=RBD_site,label_site=SARS2_site,wildtype,mutation=mutant,protein_chain="E",protein_site=SARS2_site,condition="ACE2-binding",mut_delta_effect=bind_avg)]
dmsview_output_expr <- mutants[mutant != "*",.(site=RBD_site,label_site=SARS2_site,wildtype,mutation=mutant,protein_chain="E",protein_site=SARS2_site,condition="expression",mut_delta_effect=expr_avg)]

for(i in 1:nrow(dmsview_output_bind)){
  dmsview_output_bind$mut_relative_effect[i] <- exp(dmsview_output_bind$mut_delta_effect[i])/sum(exp(dmsview_output_bind[label_site==dmsview_output_bind[i,label_site],mut_delta_effect]))
}
for(i in 1:nrow(dmsview_output_expr)){
  dmsview_output_expr$mut_relative_effect[i] <- exp(dmsview_output_expr$mut_delta_effect[i])/sum(exp(dmsview_output_expr[label_site==dmsview_output_expr[i,label_site],mut_delta_effect]))
}

dmsview_output <- rbind(dmsview_output_bind,dmsview_output_expr)

dmsview_output[condition=="ACE2-binding", c("site_mean_effect", "site_max_effect", "site_min_effect") := RBD_sites[site_SARS2==label_site,c("mean_bind","max_bind","min_bind")],by=c("label_site","mutation","condition")]
dmsview_output[condition=="expression", c("site_mean_effect", "site_max_effect", "site_min_effect") := RBD_sites[site_SARS2==label_site,c("mean_expr","max_expr","min_expr")],by=c("label_site","mutation","condition")]

write.csv(dmsview_output,file=paste(config$dms_view_dir,"/dms-view_table.csv",sep=""))
```