Circulating SARS-CoV-2 RBD variants
================
Tyler Starr
5/12/2020

This notebook analyzes RBD variants that have been sampled in isolates
within the current SARS-CoV-2 pandemic.

## Setup

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

Read in tables of variant effects on binding and expression for single
mutations to the SARS-CoV-2 RBD.

``` r
mutants <- data.table(read.csv(file=config$single_mut_effects_file,stringsAsFactors = F))

#rename mutants site indices to prevent shared names with RBD_sites, simplifying some downstream calculations that cross-index these tables
setnames(mutants, "site_RBD", "RBD_site");setnames(mutants, "site_SARS2", "SARS2_site")
```

## Analyzing amino acid diversity in GISAID Spike sequences

We constructed an alignment of all Spike sequences available on GISAID.
On the EpiCoV page, under downloads, one of the pre-made options is a
fasta of all Spike sequences isolated thus far, which is updated each
day. I have downloaded this file, unzipped, replaced spaces in fasta
headers with underscores, and aligned sequences. We load in this
alignment using the `read.fasta` function of the `bio3d` package, and
trim the alignment to RBD residues. We remove sequecnes from non-human
isolates (e.g. bat, pangolin, “environment”, mink, cat, TIGER), and then
iterate through the alignment and save any observed mutations. We then
filter mutations based on rationale below, and add counts of filtered
observations for each mutation as an ‘nobs’ colum in our overall mutants
data table.

We filter out any mutations that were *only* observed on sequences with
large numbers of missing characters – from my initial pass, I saw some
singleton amino acid variants which would require \>1 nt change, and
these were only found in a single sequence with many X amino acid
characters (the first half of the sequence was well determined, but the
second half was all X’s, with the annotated “differences” being within
short stretches between Xs with determiined amino acids), which made me
realize I needed to be careful not only of sequences rich in gap “-”
characters, but also ambiguous “X” characters. However, I didn’t want to
remove all sequences with undetermined characters off the bat, because
another pattern I saw is that for isolates bearing the N439K mutation,
\>10 are well determined across the entire RBD, but \~80 have many X
characters (in part of the RBD that is *not* near the N439K sequence
call). So, my preference would be to believe a mutation observed in an
X-rich sequence *if the variant in the sequence is observed in at least
one variant that does not contain an X*, but not believe mutations that
are *only* observed in X-rich sequences. (I noticed this issue with
N439K, but this is not the only mutation like this which is observed on
0X sequences at least once but other times on sequences with X
characters.) That is the filtering I therefore do below. This is
basically the limits of my genomic dataset sleuthing. Is there anything
else we should be doing to assess validity of observed mutants,
particularly e.g. could N439K simply be a biased sequencing error that
emerges in these Scotland sequencing samples? Would love ideas or help
in better filtering amino acid variants to retain.

``` r
alignment <- bio3d::read.fasta(file="data/alignments/Spike_GISAID/spike_GISAID_aligned.fasta", rm.dup=T)

#remove non-human samples
keep <- grep("Human",alignment$id);  alignment$ali <- alignment$ali[keep,]; alignment$id <- alignment$id[keep]

#remove columns that are gaps in first reference sequence
alignment$ali <- alignment$ali[,alignment$ali[1,]!="-"]

alignment_RBD <- alignment; alignment_RBD$ali <- alignment$ali[,RBD_sites$site_SARS2]

#check that the first sequence entry matches our reference RBD sequence
stopifnot(sum(!(alignment_RBD$ali[1,] == RBD_sites[,amino_acid_SARS2]))==0)

#remove sequences that are >5% gaps, as the amino acid calls may be generally unreliable
remove <- c()
for(i in 1:nrow(alignment_RBD$ali)){
  if(sum(alignment_RBD$ali[i,]=="-") > 0.05*ncol(alignment_RBD$ali)){remove <- c(remove,i)}
}

alignment_RBD$ali <- alignment_RBD$ali[-remove,];alignment_RBD$id <- alignment_RBD$id[-remove]

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

#remove any mutations that are *only* observed in X-rich sequences of dubious quality (keep counts in X-rich sequences if they were observed in at least one higher quality isolate)
#make a data frame that gives each observed mutation, the isolate it was observed in, and the number of X characters in that sequence
variants <- data.frame(isolate=isolates_vec,mutation=variants_vec)
for(i in 1:nrow(variants)){
  variants$number_X[i] <- sum(alignment_RBD$ali[which(alignment_RBD$id == variants[i,"isolate"]),]=="X")
}
#filter the sequence set for mutations observed in at least one X=0 background
variants_filtered <- data.frame(mutation=unique(variants[variants$number_X==0,"mutation"])) #only keep variants observed in at least one sequence with 0 X
for(i in 1:nrow(variants_filtered)){
  variants_filtered$n_obs[i] <- sum(variants$mutation == variants_filtered$mutation[i]) #but keep counts for any sequence with observed filtered muts
}

mutants[,nobs:=0]
for(i in 1:nrow(mutants)){
  if(mutants$mutation_RBD[i] %in% variants_filtered$mutation){
    mutants$nobs[i] <- variants_filtered[variants_filtered$mutation==mutants$mutation_RBD[i],"n_obs"]
  }
}
```

We see 368 amino acid polymorphisims within the 27827 sequences uploaded
in GISAID, which represents 89 of our 3819 measured missense mutants. In
the table below, we can see that many of these mutations are observed
only one or a few times, so there may still be unaccounted for
sequencinig artifacts, which we tried to account for at least minimally
with some filtering above.

``` r
kable(table(mutants[mutant!=wildtype & mutant!="*",nobs]),col.names=c("mutation count","frequency"))
```

| mutation count | frequency |
| :------------- | --------: |
| 0              |      3730 |
| 1              |        53 |
| 2              |        15 |
| 3              |         7 |
| 4              |         2 |
| 5              |         1 |
| 6              |         1 |
| 7              |         1 |
| 8              |         1 |
| 9              |         2 |
| 10             |         1 |
| 12             |         1 |
| 22             |         1 |
| 30             |         1 |
| 44             |         1 |
| 94             |         1 |

We plot each mutations experimental phenotype versus the number of times
it is observed in the circulating Spike alignment, for binding (top) and
expression (bottom), with the righthand plots simply zooming in on the
region surrounding zero for better visualization. We can see that some
of the mutations that are observed just one or a couple of times are
highly deleterious, and anything sampled more than a handful of times
exhibits \~neutral or perhaps small positive binding and/or expression
effects.

<img src="circulating_variants_files/figure-gfm/scatter_circulating_variants_nobs-1.png" style="display: block; margin: auto;" />

Here are tables giving mutations that were seen \>20 times, and those
seen any number of times with measured binding effects \>0.05. We are
currently slotted to validate the effect of N439K in both yeast display
and pseudovirus/mammalian experimental assays, and V367F, T478I, and
V483A in yeast display.

| mutation | expr\_lib1 | expr\_lib2 | expr\_avg | bind\_lib1 | bind\_lib2 | bind\_avg | nobs |
| :------- | ---------: | ---------: | --------: | ---------: | ---------: | --------: | ---: |
| V367F    |         NA |       0.74 |      0.74 |       0.02 |       0.13 |      0.07 |   22 |
| N439K    |     \-0.33 |     \-0.36 |    \-0.35 |       0.11 |     \-0.02 |      0.04 |   94 |
| T478I    |     \-0.14 |     \-0.18 |    \-0.16 |     \-0.05 |     \-0.02 |    \-0.04 |   44 |
| V483A    |       0.01 |       0.17 |      0.09 |       0.00 |     \-0.05 |    \-0.03 |   30 |

| mutation | expr\_lib1 | expr\_lib2 | expr\_avg | bind\_lib1 | bind\_lib2 | bind\_avg | nobs |
| :------- | ---------: | ---------: | --------: | ---------: | ---------: | --------: | ---: |
| G339D    |       0.21 |       0.40 |      0.30 |       0.04 |       0.07 |      0.06 |    1 |
| V367F    |         NA |       0.74 |      0.74 |       0.02 |       0.13 |      0.07 |   22 |
| K378R    |       0.07 |       0.24 |      0.16 |       0.03 |       0.13 |      0.08 |    1 |
| E406Q    |       0.06 |       0.02 |      0.04 |       0.08 |       0.06 |      0.07 |    1 |
| Q414A    |       0.35 |       0.24 |      0.30 |       0.07 |       0.24 |      0.16 |    1 |
| S477N    |       0.02 |       0.10 |      0.06 |       0.02 |       0.09 |      0.06 |    1 |
| Y508H    |       0.13 |       0.14 |      0.14 |       0.10 |       0.05 |      0.07 |    1 |

Let’s visualize the positions with interesting circulating variants in
our *favorite* exploratory heatmaps\! Below, we output maps first for
positions with circulating variants observed \>20 times (left two maps),
and second for those positions with circulating variants of at least
\>0.05 effect on binding (right two maps). These maps show the
SARS-CoV-2 wildtype state with an “x” indicator, the SARS-CoV-1 state
with an “o”, and “^” marks any amino acid variants observed at least one
time in the GISAID sequences.

<img src="circulating_variants_files/figure-gfm/heatmap_circulating_variants-1.png" style="display: block; margin: auto;" />

## Strength of selection among circulating variants

To characterize the effect of selection, we can compare the distribution
of functional effects of mutations at different n\_obs cutoffs, compared
to the “raw” distribution of functional effects – in this case, we
should look at the DFE of only those amino acid mutations that can be
introduced with single nucleotide mutations given the SARS-CoV-2
reference nt sequence.

First, we add a column to our mutants data frame indicating whether a
mutation is accessible by single nucleotide mutations.

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

mutants[,SARS2_codon:=RBD_sites[site_SARS2==SARS2_site,codon_SARS2],by=mutation]
mutants[,singlemut := mutant %in% get.codon.muts(SARS2_codon),by=mutation]
```

Are any of our observed GISAID mutations \>1nt changes? Below, we see
one mutation with a single observation that would require 2+ nucleotide
changes from the wildtype SARS-CoV-2 codon (or other synonymous codons
encoding the SARS-CoV-2 amino acid). (My first-pass analysis here showed
four such multi-nt mutations, three of which had obviously suspicious
mutation calls which caused me to go back to update my filtering
protocol. Perhaps this remaining mutation points to something additional
we could add.)

Is there anything suspicious about the isolates bearing this mutation?
This mutation is observed in the sequence
Spike|hCoV-19/USA/AZ-TG271878/2020|2020-03-20|EPI\_ISL\_426542|Original|hCoV-19^^Arizona|Human|AZ\_SPHL|TGen\_North|Lemmer|USA.
This sequence is observed to have 2 amino acid variants (D75V, Q84A),
both of which are unique to this sequence (and the former of which has
moderately deleterious effects on binding and expression). There do not
appear to be any ambiguous nts in this RBD sequence, so I’m not sure if
there is reason to filter it from our dataset besides ad hoc (we could
filter out singleton sequences containing \>1 amino acid mutation if we
think that is a generally unreliable class?).

``` r
kable(mutants[singlemut==F & nobs>0,.(mutation_RBD,mutation,expr_lib1,expr_lib2,expr_avg,bind_lib1,bind_lib2,bind_avg,nobs,SARS2_codon)])
```

| mutation\_RBD | mutation | expr\_lib1 | expr\_lib2 | expr\_avg | bind\_lib1 | bind\_lib2 | bind\_avg | nobs | SARS2\_codon |
| :------------ | :------- | ---------: | ---------: | --------: | ---------: | ---------: | --------: | ---: | :----------- |
| Q84A          | Q414A    |       0.35 |       0.24 |       0.3 |       0.07 |       0.24 |      0.16 |    1 | caa          |

Below is a heatmap of binding effects for all mutations accessible in
single nucleotide changes from the SARS-CoV-2 WT reference sequence,
with others grayed out. We can see that several affinity-enhancing
mutations (including at least two that we are currently planning to
validiate), are accessible with single-nt mutations. Therefore,
affinity-enhancing mutations, if selectively beneficial, would be
readily accessible via mutation. However, there are probably some
positions where beneficial mutations are possible, but none are
available from single-nt mutations which we could dig into if
interesting.

<img src="circulating_variants_files/figure-gfm/heatmap_binding_single_muts-1.png" style="display: block; margin: auto;" />

As has been seen in other studies, the genetic code is structured in a
conservative way such that neighboring codons exhibit similiar
biochemical properties, which causes single-nt amino acid changes to
have less deleterious effects than multiple-nt-mutant codon changes. We
see this trend in our data as well, further illustrating why we should
compare circulating variants to the single-nt-mutant amino acid effects
as the “raw” distribution of functional effects. (Median mutational
effect of single-nt mutants = -0.22; median mutational effect of all
amino acid muts = -0.34; P-value 2.3^{-4}, Wilcoxon rank-sum test.)

To illustrate how selection acts on circulating variants, let’s compare
the functional effects of observed mutations to those observed zero
times. The violin plots below show the distribution of functional
effects on binding (left) and expression (right), comparing single-nt
amino acid mutations with 0 observed counts in GISAID versus
increasingly stringent GISAID count cutoffs. We can see a bias among
circulating mutants for both binding and expression effects that are
visually by eye higher than expected by random mutation alone. This
suggests that purifying selection is removing deleterious RBD mutations
that affect traits correlated with our measured binding and expression
phenotypes.

How do we want to quantify this? The simplest statistical approach would
be permutations – draw random sets of mutation from the single-nt amino
acid mutation pool, and compare the observed e.g. median mutational
effect of mutations actually observed \>X times with the permuted sets.
This will surely show significant effects of purifying selection (median
higher in actual observed mutant set than randomly selected sets), but
is there something more interesting to be done than simply confirm the
significant shift in median? And, can we say anything with meat on it
about whether there is indeed a lack of selection for affinity-enhancing
mutations? (That is, is there a way to quantitatively compare the
probability of observing a stronger affinity-enhancing mutation than we
see in our observed set?) These are the things I want to think about
here\!

<img src="circulating_variants_files/figure-gfm/circulating_variant_DFEs-1.png" style="display: block; margin: auto;" />
