Circulating SARS-CoV-2 RBD variants
================
Tyler Starr
5/12/2020

This notebook analyzes RBD variants that have been sampled in isolates
within the current SARS-CoV-2 pandemic.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra","bio3d","seqinr")
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

#make output directory
if(!file.exists(config$circulating_variants_dir)){
  dir.create(file.path(config$circulating_variants_dir))
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
    ##  [1] seqinr_3.4-5      bio3d_2.3-4       gridExtra_2.3    
    ##  [4] forcats_0.4.0     stringr_1.4.0     dplyr_0.8.3      
    ##  [7] purrr_0.3.2       readr_1.3.1       tidyr_0.8.3      
    ## [10] tibble_2.1.3      ggplot2_3.2.0     tidyverse_1.2.1  
    ## [13] data.table_1.12.2 yaml_2.2.0        knitr_1.23       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.5 xfun_0.7         haven_2.1.1      lattice_0.20-38 
    ##  [5] colorspace_1.4-1 generics_0.0.2   htmltools_0.3.6  rlang_0.4.0     
    ##  [9] pillar_1.4.2     glue_1.3.1       withr_2.1.2      modelr_0.1.4    
    ## [13] readxl_1.3.1     munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0
    ## [17] rvest_0.3.4      evaluate_0.14    parallel_3.6.1   broom_0.5.2     
    ## [21] Rcpp_1.0.1       scales_1.0.0     backports_1.1.4  jsonlite_1.6    
    ## [25] hms_0.4.2        digest_0.6.20    stringi_1.4.3    grid_3.6.1      
    ## [29] ade4_1.7-13      cli_1.1.0        tools_3.6.1      magrittr_1.5    
    ## [33] lazyeval_0.2.2   crayon_1.3.4     pkgconfig_2.0.2  MASS_7.3-51.4   
    ## [37] xml2_1.2.0       lubridate_1.7.4  assertthat_0.2.1 rmarkdown_1.13  
    ## [41] httr_1.4.0       rstudioapi_0.10  R6_2.4.0         nlme_3.1-140    
    ## [45] compiler_3.6.1

## Setup

Read in tables of variant effects on binding and expression for single
mutations to the SARS-CoV-2 RBD.

``` r
mutants <- data.table(read.csv(file=config$single_mut_effects_file,stringsAsFactors = F))

#rename mutants site indices to prevent shared names with RBD_sites, simplifying some downstream calculations that cross-index these tables
setnames(mutants, "site_RBD", "RBD_site");setnames(mutants, "site_SARS2", "SARS2_site")
```

The 5/16/2020 GISAID report on Receptor binding surveillance for current
sequences lists the following RBD variants (number of times seen). I
believe this report focuses on just RBM variants, not all variants
observed in the rest of the RBD. I also think this is just mutants added
in the most recent data update? Not sure:

  - N439K (94x Scotland)
  - K444R (1x Spain)
  - V445I (1x England)
  - G446S (1x England)
  - G446V (2x Australia)
  - L455F (1x England)
  - F456L (1x USA)
  - A475V (2x USA)
  - G476S (9x USA, 1x Belgium)
  - T478I (44x England)
  - V483A (29x USA, 1x England)
  - V483F (4x Spain)
  - V483I (2x England)
  - F490L (1x Australia, 1x USA)
  - F490S (1x England)
  - S494P (3x USA, 1x each England, Spain, India, Sweden)
  - Y495N (1x Luxembourg)
  - V503F (1x USA)

First, let’s output the effects of these mutations in our summary
mutation table.

``` r
variants <- data.frame(mutation=c("N439K","K444R","V445I","G446S","G446V","L455F","F456L","A475V","G476S","T478I","V483A","V483F","V483I","F490L","F490S","S494P","Y495N","V503F"),
nobs=c(94,1,1,1,2,1,1,2,10,44,30,4,2,2,1,7,1,1),
ncountry=c(1,1,1,1,1,1,1,1,2,1,2,1,1,2,1,5,1,1))

for(i in 1:nrow(variants)){
variants$bind_lib1[i] <- mutants[mutation==variants$mutation[i], bind_lib1]
variants$bind_lib2[i] <- mutants[mutation==variants$mutation[i], bind_lib2]
variants$bind_avg[i] <- mutants[mutation==variants$mutation[i], bind_avg]
variants$expr_lib1[i] <- mutants[mutation==variants$mutation[i], expr_lib1]
variants$expr_lib2[i] <- mutants[mutation==variants$mutation[i], expr_lib2]
variants$expr_avg[i] <- mutants[mutation==variants$mutation[i], expr_avg]
}

kable(variants)
```

| mutation | nobs | ncountry | bind\_lib1 | bind\_lib2 | bind\_avg | expr\_lib1 | expr\_lib2 | expr\_avg |
| :------- | ---: | -------: | ---------: | ---------: | --------: | ---------: | ---------: | --------: |
| N439K    |   94 |        1 |       0.11 |     \-0.02 |      0.04 |     \-0.33 |     \-0.36 |    \-0.35 |
| K444R    |    1 |        1 |     \-0.08 |     \-0.04 |    \-0.06 |     \-0.01 |     \-0.17 |    \-0.09 |
| V445I    |    1 |        1 |     \-0.20 |       0.18 |    \-0.01 |       0.08 |     \-0.16 |    \-0.04 |
| G446S    |    1 |        1 |     \-0.19 |     \-0.22 |    \-0.20 |     \-0.25 |     \-0.55 |    \-0.40 |
| G446V    |    2 |        1 |     \-0.26 |     \-0.28 |    \-0.27 |     \-0.46 |     \-0.50 |    \-0.48 |
| L455F    |    1 |        1 |     \-0.16 |     \-0.22 |    \-0.19 |     \-0.05 |     \-0.05 |    \-0.05 |
| F456L    |    1 |        1 |     \-0.05 |     \-0.18 |    \-0.11 |     \-0.40 |     \-0.58 |    \-0.49 |
| A475V    |    2 |        1 |     \-0.11 |     \-0.17 |    \-0.14 |     \-0.10 |     \-0.32 |    \-0.21 |
| G476S    |   10 |        2 |     \-0.02 |     \-0.03 |    \-0.02 |     \-0.05 |     \-0.07 |    \-0.06 |
| T478I    |   44 |        1 |     \-0.05 |     \-0.02 |    \-0.04 |     \-0.14 |     \-0.18 |    \-0.16 |
| V483A    |   30 |        2 |       0.00 |     \-0.05 |    \-0.03 |       0.01 |       0.17 |      0.09 |
| V483F    |    4 |        1 |       0.04 |     \-0.06 |    \-0.01 |     \-0.29 |     \-0.33 |    \-0.31 |
| V483I    |    2 |        1 |     \-0.01 |     \-0.06 |    \-0.03 |     \-0.23 |     \-0.06 |    \-0.14 |
| F490L    |    2 |        2 |     \-0.04 |     \-0.17 |    \-0.11 |     \-0.38 |     \-0.32 |    \-0.35 |
| F490S    |    1 |        1 |       0.02 |     \-0.02 |      0.00 |     \-0.07 |     \-0.12 |    \-0.10 |
| S494P    |    7 |        5 |     \-0.05 |       0.06 |      0.00 |     \-0.02 |     \-0.02 |    \-0.02 |
| Y495N    |    1 |        1 |     \-1.36 |     \-1.50 |    \-1.43 |     \-2.16 |     \-2.03 |    \-2.09 |
| V503F    |    1 |        1 |       0.06 |     \-0.07 |    \-0.01 |     \-0.45 |     \-0.32 |    \-0.38 |

Let’s visualize the positions with circulating variants in our
*favorite* exploratory heatmaps\!

<img src="circulating_variants_files/figure-gfm/heatmap_circulating_variants-1.png" style="display: block; margin: auto;" />

Let’s expand our sequence set to all Spike sequences available on
GISAID. On the EpiCoV page, under downloads, one of the pre-made options
is an alignment of all Spike sequences isolated thus far, which is
updated each day. We load in this alignment using the `read.fasta`
function of the `bio3d` package, and trim the alignment to RBD residues.
(This function trims fasta header names at the first space, which causes
premature cuttiing off e.g. of “Hong Kong” sequences – in the future, we
might want to include a script in the data/alignments/Spike\_GISAID
folder which could fill spaces with underscores, remove any sequences
whose header matches “pangolin”, “bat”, etc.). We remove sequecnes from
bat and pangolin isolates, and then iterate through the alignment and
save any observed mutations. Finally, we add counts of observations for
each mutation in our mutations data table

``` r
alignment <- bio3d::read.fasta(file="data/alignments/Spike_GISAID/spikeprot0516.fasta", rm.dup=F)

alignment_RBD <- alignment; alignment_RBD$ali <- alignment$ali[,RBD_sites$site_SARS2]

#quick but dumb way to check that the first sequence entry matches our reference RBD sequence
# for(i in 1:length(alignment_RBD$ali[1,])){
#   print(alignment_RBD$ali[1,i] == RBD_sites[i,amino_acid_SARS2])
# }

#remove bat, pangolin sequences
remove <- grep("bat",alignment_RBD$id);  alignment_RBD$ali <- alignment_RBD$ali[-remove,]; alignment_RBD$id <- alignment_RBD$id[-remove]
remove <- grep("pangolin",alignment_RBD$id);  alignment_RBD$ali <- alignment_RBD$ali[-remove,]; alignment_RBD$id <- alignment_RBD$id[-remove]

#output all mutation differences from the WT/reference RBD sequence
#I do this by iterating over rows and columns of the alignment matrix which is STUPID but effective
variants_vec <- c()
isolates_vec <- c()
for(j in 1:ncol(alignment_RBD$ali)){
  #print(i)
  for(i in 1:nrow(alignment_RBD$ali)){
    if(alignment_RBD$ali[i,j] != alignment_RBD$ali[1,j] & !(alignment_RBD$ali[i,j] %in% c("X","-"))){
      variants_vec <- c(variants_vec, paste(alignment_RBD$ali[1,j],j,alignment_RBD$ali[i,j],sep=""))
      isolates_vec <- c(isolates_vec, alignment_RBD$id[i])
    }
  }
}

mutants[,nobs:=sum(mutation_RBD == variants_vec),by=mutation]
```

We see 19985 amino acid polymorphisims among sequences uploaded in
GISAID, which represent 3209 of our 3819 measured missense mutants.
(This, of course, is not making any attempt to deal with sequencing
artefacts, TC adaptations, etc.). In the table below, we can see that
most of these mutations are observed only one or a few times, so these
may either be sequencing errors or unfit genotypes. There is sort of a
second “mode” of observed mutation counts around \~35 – I wonder if the
statistics of sequencing errors could somehow explain these two modes?

We plot each mutations delta-log<sub>10</sub>(*K*<sub>A,app</sub>)
versus the number of times it is observed in the circulating Spike
alignment, with nobs=0 points in red to distinguish from \>0. We can see
that many of the mutations that are observed \<25 times are highly
deleterious, though some of the mutations with \~40 observations are
also highly deleterious. The plot on the right zooms in on the range
closer to neutral, to better see any distribution \>0. There is no
significant difference in median
delta-log<sub>10</sub>(*K*<sub>A,app</sub>) for mutations that are
observed more than 25 times versus those observed 1 to 25 times (P-value
0.7278863, Wilcoxon rank-sum test). Similarly, there is no significant
difference in median delta-expression for mutations observed more than
25 times versus those observed 1 to 25 times (P-value 0.8070109,
Wilcoxon rank-sum test). Similarly, there is not a significantly
different median effect on binding or expression for mutations observed
\>0 versus 0 times (excluding nonsense mutants – P-value 0.9606674 and
0.6399633 for binding and expression, respectively. P-values are also
non-significant for \>25 nobs versus 0 nobs comparisons).

``` r
table(mutants$nobs)
```

    ## 
    ##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    ## 1012  584  527  505  364  259  187  116  130  107   78   62   34   23   16 
    ##   15   16   17   18   19   20   21   24   33   34   35   36   37   38   39 
    ##    5    5    4    2    1    1    1    1    7   20   39   32   31   12   11 
    ##   40   41   42   43   44   45   46   47   49   95 
    ##   11    9   10    5    4    1    2    1    1    1

<img src="circulating_variants_files/figure-gfm/scatter_circulating_variants_nobs-1.png" style="display: block; margin: auto;" />

This preliminary anlaysis, though perhaps with its caveats given its
simplified parsing of mutations in the alignment, shows the overall
distribution of functional effects of observed mutations does not differ
from the overall distribution of functional effects of mutation,
suggesting some combination of a) purifying selection on RBD being
relatively weak (could be consistent with big population expansion into
a susceptible population with little inter-strain competition) and b)
these variants are deleterious tip sequences that are stochastically
sampled but are in the process of being purged (i.e. normal observation
that viral tips often have deleterious mutations even though purifying
selection along the trunk lineage is retained – hence why all of the
observed homolog RBDs have a tight distribution of expression scores
even though so many mutations reduce affinity, including those sampled
in SARS-CoV-2 isolates). Anyway, if these are points we wanted to
pursue, we should think about how we should properly be approaching
these conclusions.

In addition to touching on the statement that purifying selection might
not be super strong on the RBD, we also can see that there isn’t really
an enrichment of affinity- (or stability) enhancing mutations among the
\>20,000 sequenced isolates – if we show that plenty of single-nt
mutations in the SARS-CoV-2 RBD *could* enhance ACE2-affinity, I think
we could make a somewhat interesting conclusiona round the idea that
there is no strong selective benefit for enhanced ACE2-affinity, at
least so far at this stage in the pandemic. (And that sufficient
adaptation to huACE2 occurred prior to the origins of the current
pandemic).

Let’s look not only at the overall DFE, but the distribution of
funcitonal effects of single-nt changes to the SARS-CoV-2 RBD genetic
sequence. The RBD\_sites table already gives the codon sequence for each
site, so we can quickly output

``` r
#define a function that takes a character of three nucleotides (a codon), and outputs all amino acids that can be accessed via single-nt mutation of that codon
get.codon.muts <- function(codon){
  nt <- c("a","c","g","t")
  codon_split <- strsplit(codon,split="")[[1]]
  codon_muts <- vector()
  for(i in nt[nt!=codon_split[1]]){
    codon_muts <- c(codon_muts,translate(c(i,codon_split[2:3])))
  }
  for(i in nt[nt!=codon_split[2]]){
    codon_muts <- c(codon_muts,translate(c(codon_split[1],i,codon_split[3])))
  }
  for(i in nt[nt!=codon_split[3]]){
    codon_muts <- c(codon_muts,translate(c(codon_split[1:2],i)))
  }
  return(codon_muts)
}

#make mutants table that is restricted to single nt mutants of the SARS-CoV-2 sequence
mutants_1mut <- copy(mutants)
mutants_1mut[,codon:=RBD_sites[site_SARS2==SARS2_site,codon_SARS2],by=mutation]
mutants_1mut[,singlemut := mutant %in% get.codon.muts(codon),by=mutation]
mutants_1mut[singlemut==F,bind_avg:=NA]
mutants_1mut[singlemut==F,expr_avg:=NA]
```

Below is a heatmap for all mutations accessible in single nucleotide
changes from the SARS-CoV-2 WT reference sequence, with others grayed
out. We can see that several affinity-enhancing mutations (including at
least two that we are currently planning to validiate), are accessible
with single-nt mutations. Therefore, affinity-enhancing mutations, if
selectively beneficial, would be readily accessible via mutation.
However, there are some positions where beneficial mutations are
possible, but none are available from single-nt mutations (e.g. site
455, perhaps others if I dig in more…)

<img src="circulating_variants_files/figure-gfm/heatmap_binding_single_muts-1.png" style="display: block; margin: auto;" />

As has been seen in other studies, attributed to the conservative nature
of the genetic code with regards to codon mutational neighbors tending
toward conservative biochemical properties, we can see that amino acid
changes introduced with single nt changes have a smaller median
deleterious effect than all amino acid changes agnostic of codon
structure.

``` r
wilcox.test(mutants_1mut$bind_avg,mutants$bind_avg)
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  mutants_1mut$bind_avg and mutants$bind_avg
    ## W = 3152525, p-value = 9.005e-14
    ## alternative hypothesis: true location shift is not equal to 0

``` r
median(mutants_1mut[mutant!=wildtype,bind_avg],na.rm=T);median(mutants[mutant!=wildtype,bind_avg],na.rm=T)
```

    ## [1] -0.22

    ## [1] -0.34
