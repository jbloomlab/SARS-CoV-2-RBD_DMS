# Qualitative QC analysis of counts in each sample
This Python Jupyter notebook performs a qualitative analysis of the counts for the variants in each sample.
It then computes functional scores for each variant.
These are heuristic functional scores based on overall enrichment, not the more detailed bin-based scores generated elsewhere.

## Set up analysis

### Import Python modules.
This notebook primarily makes use of the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package, and uses [plotnine](https://github.com/has2k1/plotnine) for ggplot2-like plotting syntax:


```python
import collections
import os
import warnings

import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.plotnine_themes

from IPython.display import display, HTML, Image

import numpy

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://github.com/has2k1/plotnine) theme to the gray-grid one defined in [dms_variants](https://jbloomlab.github.io/dms_variants):


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using dms_variants version {dms_variants.__version__}")
print(f"Using alignparse version {alignparse.__version__}")
```

    Using dms_variants version 0.6.0
    Using alignparse version 0.1.3


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
os.makedirs(config['func_scores_dir'], exist_ok=True)
```

Read information about the samples:


```python
sample_info = pd.read_csv(config['barcode_runs'])
```

Get lists of the *SortSeq* and *TiteSeq* samples:


```python
sortseq_samples = sample_info.query('sample_type == "SortSeq"')['sample'].unique().tolist()
titeseq_samples = sample_info.query('sample_type == "TiteSeq"')['sample'].unique().tolist()
```

## Initialize codon-variant table
Initialize [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) from wildtype gene sequence and the variant counts CSV file.
We will then use the plotting functions of this variant table to analyze the counts per sample:


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(
                seqsfile=config['amplicons'],
                feature_parse_specs=config['feature_parse_specs'])
geneseq = targets.get_target(config['primary_target']).get_feature('gene').seq
print(f"Read gene of {len(geneseq)} nts for {config['primary_target']} from {config['amplicons']}")
      
print('Initializing CodonVariantTable from gene sequence and ' +
      config['variant_counts_file'])
      
variants = dms_variants.codonvarianttable.CodonVariantTable.from_variant_count_df(
                geneseq=geneseq,
                variant_count_df_file=config['variant_counts_file'],
                primary_target=config['primary_target'])
      
print('Done initializing CodonVariantTable.')
```

    Read gene of 603 nts for SARS-CoV-2 from data/PacBio_amplicons.gb
    Initializing CodonVariantTable from gene sequence and results/counts/variant_counts.csv
    Done initializing CodonVariantTable.


## Analyze counts for samples
Analyze the variant counts per cell for each sample:


```python
counts = (
    variants.variant_count_df
    .groupby(['library', 'sample'])
    .aggregate({'count': 'sum'})
    .merge(sample_info, on=['library', 'sample'], validate='one_to_one')
    .merge(variants.n_variants_df(samples=None)
                   .groupby('library')
                   .aggregate(n_variants=pd.NamedAgg('count', 'sum')),
           on='library', validate='many_to_one')
    .assign(counts_per_cell=lambda x: x['count'] / x['number_cells'])
    .drop(columns=['R1', 'date'])
    )

counts.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>library</th>
      <th>sample</th>
      <th>count</th>
      <th>sample_type</th>
      <th>sort_bin</th>
      <th>concentration</th>
      <th>number_cells</th>
      <th>n_variants</th>
      <th>counts_per_cell</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>lib1</td>
      <td>SortSeq_bin1</td>
      <td>22911930</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>6600000</td>
      <td>99648</td>
      <td>3.471505</td>
    </tr>
    <tr>
      <th>1</th>
      <td>lib1</td>
      <td>SortSeq_bin2</td>
      <td>11262091</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>3060000</td>
      <td>99648</td>
      <td>3.680422</td>
    </tr>
    <tr>
      <th>2</th>
      <td>lib1</td>
      <td>SortSeq_bin3</td>
      <td>9598488</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>2511000</td>
      <td>99648</td>
      <td>3.822576</td>
    </tr>
    <tr>
      <th>3</th>
      <td>lib1</td>
      <td>SortSeq_bin4</td>
      <td>11220871</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>2992000</td>
      <td>99648</td>
      <td>3.750291</td>
    </tr>
    <tr>
      <th>4</th>
      <td>lib1</td>
      <td>TiteSeq_01_bin1</td>
      <td>2700607</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>1098088</td>
      <td>99648</td>
      <td>2.459372</td>
    </tr>
  </tbody>
</table>
</div>



For the *SortSeq* arm plot the average counts per variant, the average cells per variant, and the average counts per cell for each library / sample:


```python
p = (
    ggplot(counts.query('sample in @sortseq_samples')
                 .assign(counts_per_variant=lambda x: x['count'] / x['n_variants'],
                         cells_per_variant=lambda x: x['number_cells'] / x['n_variants'])
                 .melt(id_vars=['library', 'sample'],
                       value_vars=['counts_per_variant', 'cells_per_variant', 'counts_per_cell'],
                       var_name='statistic'),
           aes('sample', 'value')) +
    geom_bar(stat='identity') +
    facet_grid('statistic ~ library', scales='free_y') +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.4 * len(sortseq_samples) * counts['library'].nunique(), 6),
          panel_grid_major_x=element_blank()) +
    xlab('')
    )

_ = p.draw()
```


![png](analyze_counts_files/analyze_counts_23_0.png)


For the *TiteSeq* arm, plot cells per variant aggregated across all bins for each concentration:


```python
p = (
    ggplot(counts.query('sample in @titeseq_samples')
                 .groupby(['library', 'concentration', 'n_variants'])
                 .aggregate({'number_cells': 'sum'})
                 .reset_index()
                 .assign(cells_per_variant=lambda x: x['number_cells'] / x['n_variants'],
                         concentration=lambda x: pd.Categorical(x['concentration'].astype(int))),
           aes('concentration', 'cells_per_variant')) +
    geom_bar(stat='identity') +
    facet_grid('~ library', scales='nrow') +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.3 * counts['concentration'].nunique() * counts['library'].nunique(), 2),
          panel_grid_major_x=element_blank())
    )

_ = p.draw()
```


![png](analyze_counts_files/analyze_counts_25_0.png)


Also for the *TiteSeq* arm, plot counts per cell for each sample (bin and concentration):


```python
lower_limit = 0.1  # plot values less than this as this

p = (
    ggplot(counts.query('sample in @titeseq_samples')
                 .assign(sort_bin=lambda x: pd.Categorical(x['sort_bin']),
                         concentration=lambda x: pd.Categorical(x['concentration'].astype(int)),
                         counts_per_cell=lambda x: numpy.clip(x['counts_per_cell'],
                                                              lower_limit, None)
                         ),
           aes('concentration', 'counts_per_cell', color='sort_bin')) +
    geom_point(size=2) +
    facet_wrap('~ library') +
    scale_color_manual(values=CBPALETTE) +
    scale_y_log10() +
    geom_hline(yintercept=lower_limit, color=CBPALETTE[-1], linetype='dotted') +
    theme(axis_text_x=element_text(angle=90),
          figure_size=(0.3 * counts['concentration'].nunique() * counts['library'].nunique(), 2),
          panel_grid_major_x=element_blank()
          )
    )

_ = p.draw()
```


![png](analyze_counts_files/analyze_counts_27_0.png)


## Mutations per variant
Average number of mutations per gene among all variants of the primary target for *SortSeq* samples:


```python
p = variants.plotNumCodonMutsByType(variant_type='all',
                                    orientation='v',
                                    libraries=variants.libraries,
                                    samples=sortseq_samples,
                                    heightscale=0.8)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
```


![png](analyze_counts_files/analyze_counts_29_0.png)


Mutation frequency across gene among all variants for primary target for *SortSeq* samples only:


```python
p = variants.plotMutFreqs('all',
                          'codon',
                          orientation='v',
                          libraries=variants.libraries,
                          samples=sortseq_samples)
_ = p.draw()
```


![png](analyze_counts_files/analyze_counts_31_0.png)


Now plot the average mutation frequency for the *TiteSeq* samples.
To make the plot not too big, chunk the plots out for smaller subsets:


```python
i = 0
chunksize = 8

while i < len(titeseq_samples):
    p = variants.plotNumCodonMutsByType(
                                    variant_type='all',
                                    orientation='v',
                                    libraries=variants.libraries,
                                    samples=titeseq_samples[i: i + chunksize],
                                    heightscale=0.9,
                                    widthscale=1.1)
    p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
    _ = p.draw()
    i += chunksize
```


![png](analyze_counts_files/analyze_counts_33_0.png)



![png](analyze_counts_files/analyze_counts_33_1.png)



![png](analyze_counts_files/analyze_counts_33_2.png)



![png](analyze_counts_files/analyze_counts_33_3.png)



![png](analyze_counts_files/analyze_counts_33_4.png)



![png](analyze_counts_files/analyze_counts_33_5.png)



![png](analyze_counts_files/analyze_counts_33_6.png)



![png](analyze_counts_files/analyze_counts_33_7.png)

