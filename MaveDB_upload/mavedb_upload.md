Output MaveDB uploads
================
Tyler Starr
6/18/2020

  - [MaveDB organization](#mavedb-organization)
  - [Per-barcode score set
    generation](#per-barcode-score-set-generation)
  - [Per-mutation score set
    generation](#per-mutation-score-set-generation)

This notebook outputs DMS measurements in the format required for upload
to [MaveDB](https://www.mavedb.org)

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","seqinr")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

#read in config file
config <- read_yaml("../config.yaml")
```

## MaveDB organization

I defined an “Experiment Set” on MaveDB, which contains two
“Experiments”: “Effects on binding of mutation in the SARS-CoV-2
RBD” and “Effects on expression of mutations in the SARS-CoV-2 RBD”.
Within each of these experiments, we will upload two “Score Sets”,
consisting of our measurements
(delta-log<sub>10</sub>(*K*<sub>D,app</sub>) for binding and
delta-log(MFI) for expression), at the level of individual barcodes, and
our final single-mutation-effect dataset.

## Per-barcode score set generation

Let’s read in the per-barcode binding and expression measurements.

``` r
bc_bind <- data.table(read.csv(file=paste("../",config$Titeseq_Kds_file,sep="")))
head(bc_bind)
```

    ##    library     target          barcode variant_call_support avgcount
    ## 1:    lib1 SARS-CoV-2 AAAAAAAAAATTGTAA                    1     0.00
    ## 2:    lib1 SARS-CoV-2 AAAAAAAAACTTAAAT                    2     4.30
    ## 3:    lib1 SARS-CoV-2 AAAAAAAAACTTCAAT                    5    74.50
    ## 4:    lib1 SARS-CoV-2 AAAAAAAACAAGCAGA                    6   146.32
    ## 5:    lib1 SARS-CoV-2 AAAAAAAACAATATAA                    1     2.48
    ## 6:    lib1 SARS-CoV-2 AAAAAAAACAGGTTGC                    4    47.48
    ##    log10Ka delta_log10Ka log10SE response baseline nMSR    variant_class
    ## 1:      NA            NA      NA       NA       NA   NA >1 nonsynonymous
    ## 2:      NA            NA      NA       NA       NA   NA >1 nonsynonymous
    ## 3:    8.72         -2.05    0.12     1.64     1.13 0.01 >1 nonsynonymous
    ## 4:   10.35         -0.42    0.05     2.87     1.01 0.00  1 nonsynonymous
    ## 5:      NA            NA      NA       NA       NA   NA >1 nonsynonymous
    ## 6:    6.00         -4.76    2.00     1.50     1.15 0.01 >1 nonsynonymous
    ##              aa_substitutions n_aa_substitutions
    ## 1:                 Y91L K199Y                  2
    ## 2: N13S L60P K94N S147T C150Y                  5
    ## 3:     A22C R127G E141D L188V                  4
    ## 4:                       N13F                  1
    ## 5:        C6K T15W K94Y V103W                  4
    ## 6:           V71K P149L N157T                  3

``` r
bc_expr <- data.table(read.csv(file=paste("../",config$expression_sortseq_file,sep="")))
head(bc_expr)
```

    ##    library     target          barcode variant_call_support total_count
    ## 1:    lib1 SARS-CoV-2 AAAAAAAAAATTGTAA                    1        0.00
    ## 2:    lib1 SARS-CoV-2 AAAAAAAAACTTAAAT                    2       64.71
    ## 3:    lib1 SARS-CoV-2 AAAAAAAAACTTCAAT                    5      117.96
    ## 4:    lib1 SARS-CoV-2 AAAAAAAACAAGCAGA                    6      244.34
    ## 5:    lib1 SARS-CoV-2 AAAAAAAACAATATAA                    1       95.35
    ## 6:    lib1 SARS-CoV-2 AAAAAAAACAGGTTGC                    4      212.43
    ##    ML_meanF delta_ML_meanF var_ML_meanF    variant_class
    ## 1:       NA             NA           NA >1 nonsynonymous
    ## 2:     7.45          -3.01         0.04 >1 nonsynonymous
    ## 3:     7.92          -2.54         0.03 >1 nonsynonymous
    ## 4:     8.93          -1.53         0.01  1 nonsynonymous
    ## 5:     6.21          -4.25         0.03 >1 nonsynonymous
    ## 6:     7.73          -2.73         0.02 >1 nonsynonymous
    ##              aa_substitutions n_aa_substitutions
    ## 1:                 Y91L K199Y                  2
    ## 2: N13S L60P K94N S147T C150Y                  5
    ## 3:     A22C R127G E141D L188V                  4
    ## 4:                       N13F                  1
    ## 5:        C6K T15W K94Y V103W                  4
    ## 6:           V71K P149L N157T                  3

We need to rename the genotype column to `hgvs_pro` and conform to [HGVS
guidelines](http://varnomen.hgvs.org/recommendations/general/). For
protein mutations, we:

  - Count from the first mutated position as position 1 (I believe we
    can add the offset for spike numbering later on)
  - Convert to three-letter amino acid codes
  - For single mutants, add the `p.` designator to specify amino acid
    changes, e.g. `p.Asn1Gln`
  - For multiple mutants, format as `p.[Asn1Gln;Ile2Leu]`
  - Stop mutants stay as single letter `*` character

We then rename our functional score for each table to be `score`

We include as optional columns for the binding data:

  - `avg_count` (currently named `avgcount`), `barcode`, and `library`

We include as optional columns for the expression data:

  - `count` (currently named `total_count`), `barcode`, and `library`

Below, we make these modifications and output the resulting `.csv` files

``` r
#modify aaa() function to not modify * to Stp, as HGVS still prefers * for stop even in three letter code
aaa_stop <- function(x){if(x=="*"){return(x)}else{return(aaa(x))}}

#function to return HGVS nomenclature mutants from our one letter code notation
hgvs_naming <- function(x){
  if(x==""){
    return("")
  }else{
    muts <- strsplit(as.character(x),split=" ")[[1]]
    if(length(muts)==1){
      split <- strsplit(muts,split="")[[1]]
      split[1] <- aaa_stop(split[1])
      split[length(split)] <- aaa_stop(split[length(split)])
      paste <- paste("p.",paste(split,collapse=""),sep="")
      return(paste)
    }else if(length(muts) > 1){
      for(i in 1:length(muts)){
        split <- strsplit(muts[i],split="")[[1]]
        split[1] <- aaa_stop(split[1])
        split[length(split)] <- aaa_stop(split[length(split)])
        paste <- paste(split,collapse="")
        muts[i] <- paste
      }
      paste("p.[",paste(muts, collapse=";"),"]",sep="")
    }
  }
}

bc_bind[, hgvs_pro := hgvs_naming(aa_substitutions),by="aa_substitutions"]

setnames(bc_bind, "delta_log10Ka", "score"); setnames(bc_bind, "avgcount", "avg_count")

bc_bind <- bc_bind[,.(hgvs_pro, score, avg_count, library)]

bc_bind[library=="lib1",library:="1"]
bc_bind[library=="lib2",library:="2"]

#empty unmutated sequence not allowed. Try "p.="
bc_bind[hgvs_pro=="",hgvs_pro:="p.="]

head(bc_bind)
```

    ##                                              hgvs_pro score avg_count
    ## 1:                             p.[Tyr91Leu;Lys199Tyr]    NA      0.00
    ## 2: p.[Asn13Ser;Leu60Pro;Lys94Asn;Ser147Thr;Cys150Tyr]    NA      4.30
    ## 3:         p.[Ala22Cys;Arg127Gly;Glu141Asp;Leu188Val] -2.05     74.50
    ## 4:                                         p.Asn13Phe -0.42    146.32
    ## 5:            p.[Cys6Lys;Thr15Trp;Lys94Tyr;Val103Trp]    NA      2.48
    ## 6:                   p.[Val71Lys;Pro149Leu;Asn157Thr] -4.76     47.48
    ##    library
    ## 1:       1
    ## 2:       1
    ## 3:       1
    ## 4:       1
    ## 5:       1
    ## 6:       1

``` r
write.csv(bc_bind,file="score-set_binding_bc.csv",row.names=F)
```

``` r
bc_expr[, hgvs_pro := hgvs_naming(aa_substitutions),by="aa_substitutions"]

setnames(bc_expr, "delta_ML_meanF", "score"); setnames(bc_expr, "total_count", "count")

bc_expr <- bc_expr[,.(hgvs_pro, score, count, library)]

bc_expr[library=="lib1",library:="1"]
bc_expr[library=="lib2",library:="2"]

#empty unmutated sequence not allowed. Try "p.="
bc_expr[hgvs_pro=="",hgvs_pro:="p.="]

head(bc_expr)
```

    ##                                              hgvs_pro score  count library
    ## 1:                             p.[Tyr91Leu;Lys199Tyr]    NA   0.00       1
    ## 2: p.[Asn13Ser;Leu60Pro;Lys94Asn;Ser147Thr;Cys150Tyr] -3.01  64.71       1
    ## 3:         p.[Ala22Cys;Arg127Gly;Glu141Asp;Leu188Val] -2.54 117.96       1
    ## 4:                                         p.Asn13Phe -1.53 244.34       1
    ## 5:            p.[Cys6Lys;Thr15Trp;Lys94Tyr;Val103Trp] -4.25  95.35       1
    ## 6:                   p.[Val71Lys;Pro149Leu;Asn157Thr] -2.73 212.43       1

``` r
write.csv(bc_expr,file="score-set_expression_bc.csv",row.names=F)
```

## Per-mutation score set generation

We’ll perform similar actions to output score sets consisting of our
final estimates of single-mutant effects.

``` r
muts <- data.table(read.csv(file=paste("../",config$single_mut_effects_file,sep=""),stringsAsFactors=F))
muts <- muts[mutant != wildtype,]
head(muts)
```

    ##    site_RBD site_SARS2 wildtype mutant mutation mutation_RBD bind_lib1
    ## 1:        1        331        N      A    N331A          N1A     -0.05
    ## 2:        1        331        N      C    N331C          N1C     -0.08
    ## 3:        1        331        N      D    N331D          N1D      0.00
    ## 4:        1        331        N      E    N331E          N1E      0.02
    ## 5:        1        331        N      F    N331F          N1F     -0.03
    ## 6:        1        331        N      G    N331G          N1G     -0.06
    ##    bind_lib2 bind_avg expr_lib1 expr_lib2 expr_avg
    ## 1:     -0.02    -0.03     -0.14     -0.08    -0.11
    ## 2:     -0.10    -0.09     -1.56     -0.97    -1.26
    ## 3:      0.07     0.03     -0.75     -0.12    -0.44
    ## 4:     -0.02     0.00     -0.39     -0.24    -0.31
    ## 5:     -0.16    -0.10     -0.83     -0.57    -0.70
    ## 6:     -0.02    -0.04     -0.21     -0.29    -0.25

We need to rename the `mutation_RBD` column to `hgvs_pro` and conform to
[HGVS guidelines](http://varnomen.hgvs.org/recommendations/general/).
For protein mutations, we:

  - Count from the first mutated position as position 1 (I believe we
    can add the offset for spike numbering later on)
  - Convert to three-letter amino acid codes
  - Add the `p.` designator to specify amino acid changes,
    e.g. `p.Asn1Gln`
  - Stop mutants stay as single letter `*` character

We then split the table into binding and expression tables, and rename
our functional score for each table to be `score`

We include as optional columns for the both datasets:

  - `library` (so, collapse to long form instead of wide), and we’ll
    also include the `average` with each row as the average of the two
    measurements

Below, we make these modifications and output the resulting `.csv` files

``` r
muts[,hgvs_pro := paste("p.",aaa_stop(wildtype),site_RBD,aaa_stop(mutant),sep=""),by=mutation]
mut_expr <- muts

mut_bind_lib1 <- muts[,.(hgvs_pro,bind_lib1,bind_avg)]
mut_bind_lib1[,library:="1"]
mut_bind_lib2 <- muts[,.(hgvs_pro,bind_lib2,bind_avg)]
mut_bind_lib2[,library:="2"]

mut_expr_lib1 <- muts[,.(hgvs_pro,expr_lib1,expr_avg)]
mut_expr_lib1[,library:="1"]
mut_expr_lib2 <- muts[,.(hgvs_pro,expr_lib2,expr_avg)]
mut_expr_lib2[,library:="2"]

setnames(mut_bind_lib1, c("bind_lib1", "bind_avg"), c("score", "average"))
setnames(mut_bind_lib2, c("bind_lib2", "bind_avg"), c("score", "average"))
setnames(mut_expr_lib1, c("expr_lib1", "expr_avg"), c("score", "average"))
setnames(mut_expr_lib2, c("expr_lib2", "expr_avg"), c("score", "average"))

mut_bind <- rbind(mut_bind_lib1, mut_bind_lib2)

head(mut_bind)
```

    ##     hgvs_pro score average library
    ## 1: p.Asn1Ala -0.05   -0.03       1
    ## 2: p.Asn1Cys -0.08   -0.09       1
    ## 3: p.Asn1Asp  0.00    0.03       1
    ## 4: p.Asn1Glu  0.02    0.00       1
    ## 5: p.Asn1Phe -0.03   -0.10       1
    ## 6: p.Asn1Gly -0.06   -0.04       1

``` r
write.csv(mut_bind,file="score-set_binding_mut.csv",row.names=F)
```

``` r
mut_expr <- rbind(mut_expr_lib1, mut_expr_lib2)

head(mut_expr)
```

    ##     hgvs_pro score average library
    ## 1: p.Asn1Ala -0.14   -0.11       1
    ## 2: p.Asn1Cys -1.56   -1.26       1
    ## 3: p.Asn1Asp -0.75   -0.44       1
    ## 4: p.Asn1Glu -0.39   -0.31       1
    ## 5: p.Asn1Phe -0.83   -0.70       1
    ## 6: p.Asn1Gly -0.21   -0.25       1

``` r
write.csv(mut_expr,file="score-set_expression_mut.csv",row.names=F)
```
