Titer Spike pseudovirus by flow cytometry
================

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))
#list of packages to install/load
packages = c("ggplot2", "data.table", "tidyverse", "dplyr", "broom", "gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))
#make results directory
if(!file.exists("results")){
 dir.create(file.path("results"))
}
```

## Experiment: Titer Spike point mutant pseudovirus by flow cytometry

**2020-06-01** Mutants tested: Neutral in bulk yeast display
experiments: \* wildtype \* N439K

Deleterious in bulk yeast display experiments: \* N501D \* L455Y \*
C432D \* G502D

Beneficial in bulk yeast display experiments: \* N501F \* Q498Y

### Define colorblind-friendly palette

``` r
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")
# The palette with black
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")
```

### Read in data from flow cytometry

``` r
dt <- read.csv(file="pseudovirus_titers.csv", stringsAsFactors=F)
head(dt)
```

    ##   cells single_cells single_cells2 single_cells2_FITC_geomean
    ## 1  87.4         94.4          85.2                       3950
    ## 2  84.5         93.9          80.9                       4269
    ## 3  89.2         96.1          89.4                       3560
    ## 4  85.5         95.7          82.1                       5634
    ## 5  83.4         95.8          80.2                       6432
    ## 6  84.1         96.3          80.8                       5089
    ##   FITC_perc_polygon_gate FITC_perc_linegate WELL.ID plasmid genotype
    ## 1                   87.3              11.60     C01    2736 wildtype
    ## 2                   91.4              11.40     C02    2736 wildtype
    ## 3                   86.7               6.08     C03    2736 wildtype
    ## 4                   95.3              22.20     C04    2759    Q498Y
    ## 5                   97.0              27.40     C05    2759    Q498Y
    ## 6                   94.6              17.70     C06    2759    Q498Y
    ##   expectation n_cells uL_virus virus_per_mL_stringent virus_per_mL_relaxed
    ## 1     neutral   12700       50                  29500               222000
    ## 2     neutral   12700       50                  29000               232000
    ## 3     neutral   12700       50                  15400               220000
    ## 4  beneficial   12700       50                  56400               242000
    ## 5  beneficial   12700       50                  69600               246000
    ## 6  beneficial   12700       50                  45000               240000

### Plot titers

``` r
dt <- dt %>% mutate(genotype = factor(genotype, 
                                levels=c("none", "VSV_G", "wildtype", 
                                         "N439K", "Q498Y", "N501F", 
                                         "N501D", "G502D", "L455Y", "C432D")),
              uL_virus = factor(uL_virus),
              yeast_display = factor(expectation,
                                   levels=c("neutral", "beneficial", "deleterious", "none"))
              )

ggplot(dt %>%
         filter(uL_virus == "50" | genotype=="none")
       ,
       aes(genotype, virus_per_mL_stringent/exp(mean(log(dt[dt$genotype=="wildtype" & dt$uL_virus=="50","virus_per_mL_stringent"]))))
       ) +
  geom_point(size=3.5,
             shape=16,
             position=position_jitter(width = 0.1, height = 0, seed = 0), 
             aes(fill=yeast_display)) +
  geom_hline(yintercept=25/exp(mean(log(dt[dt$genotype=="wildtype" & dt$uL_virus=="50","virus_per_mL_stringent"]))), color="grey", linetype="dashed") + #Limit of detection
  scale_y_log10(lim=c(8e-4, 1e1)) +
  xlab("Spike gene") +
  ylab("relative titer") +
  theme_bw()  +
  scale_fill_manual(name="effect on ACE2 binding \nin yeast display",
                         breaks=c("neutral", "beneficial", "deleterious", "none"),
                         labels=c("neutral", "increased affinity", "decreased affinity", "no Spike"),
                      values=c(cbPalette[1:3], "black")) +
  ggtitle("Spike-pseudotyped virus") +
  theme(plot.title = element_text(hjust = 0.5))
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](pseudovirus_titer_files/figure-gfm/unnamed-chunk-2-1.svg)<!-- -->

``` r
ggsave(
  "./results/stringent_titers.pdf",
  scale = 1,
  width = 8,
  height = 4,
  useDingbats=F
)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
