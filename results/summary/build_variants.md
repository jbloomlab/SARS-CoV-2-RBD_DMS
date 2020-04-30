# Build barcode variant table
This Python Jupyter notebook builds consensus sequences for barcoded variants from the mutations called in the processed PacBio CCSs.
It then uses these consensus sequences to build a codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://github.com/has2k1/plotnine) for ggplot2-like plotting.

The analysis uses the Bloom lab's [alignparse](https://jbloomlab.github.io/alignparse) and [dms_variants](https://jbloomlab.github.io/dms_variants) packages:


```python
import os
import re
import warnings

import alignparse
import alignparse.consensus
import alignparse.targets

import dms_variants
import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import numpy

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://github.com/has2k1/plotnine) theme to the gray grid one defined in [dms_variants](https://jbloomlab.github.io/dms_variants):


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.1.3
    Using dms_variants version 0.6.0


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

Make output directories if needed:


```python
os.makedirs(config['variants_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Build variant consensus sequences
We now use the processed CCSs to build a consensus sequence (really, list of mutations) for each barcoded variant in each library.

### Get processed CCSs
Read the CSV file with the processed CCSs into a data frame:


```python
processed_ccs = pd.read_csv(config['processed_ccs_file'], na_filter=None)

nlibs = processed_ccs['library'].nunique()  # number of unique libraries

ntargets = processed_ccs['target'].nunique()  # number of unique targets

print(f"Read {len(processed_ccs)} CCSs from {nlibs} libraries and {ntargets} targets.")
```

    Read 999925 CCSs from 2 libraries and 12 targets.


Here are the first few lines:


```python
display(HTML(processed_ccs.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>query_name</th>
      <th>gene_mutations</th>
      <th>gene_accuracy</th>
      <th>barcode</th>
      <th>barcode_accuracy</th>
      <th>target</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib1</td>
      <td>m54228_200414_185350/4391490/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>ATTTCCATCAAATTAA</td>
      <td>1.0</td>
      <td>BM48-31</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>m54228_200414_185350/5243143/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>GCGCTGTACCCTCTGG</td>
      <td>1.0</td>
      <td>BM48-31</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>m54228_200414_185350/5243413/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>TTCAATTAAGTCCGCC</td>
      <td>1.0</td>
      <td>BM48-31</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>m54228_200414_185350/5505876/ccs</td>
      <td></td>
      <td>1.0</td>
      <td>GTTTCACTGTGATATA</td>
      <td>1.0</td>
      <td>BM48-31</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>m54228_200414_185350/5636186/ccs</td>
      <td>T355C</td>
      <td>1.0</td>
      <td>CAGCTAATACTTCTCC</td>
      <td>1.0</td>
      <td>BM48-31</td>
    </tr>
  </tbody>
</table>


Overall statistics on number of total CCSs and number of unique barcodes:


```python
display(HTML(
    processed_ccs
    .groupby(['target', 'library'])
    .aggregate(total_CCSs=('barcode', 'size'),
               unique_barcodes=('barcode', 'nunique'))
    .assign(avg_CCSs_per_barcode=lambda x: x['total_CCSs'] / x['unique_barcodes'])
    .round(2)
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>total_CCSs</th>
      <th>unique_barcodes</th>
      <th>avg_CCSs_per_barcode</th>
    </tr>
    <tr>
      <th>target</th>
      <th>library</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">BM48-31</th>
      <th>lib1</th>
      <td>288</td>
      <td>110</td>
      <td>2.62</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>459</td>
      <td>121</td>
      <td>3.79</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">GD-Pangolin</th>
      <th>lib1</th>
      <td>338</td>
      <td>122</td>
      <td>2.77</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>452</td>
      <td>129</td>
      <td>3.50</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU3-1</th>
      <th>lib1</th>
      <td>269</td>
      <td>108</td>
      <td>2.49</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>432</td>
      <td>123</td>
      <td>3.51</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">LYRa11</th>
      <th>lib1</th>
      <td>252</td>
      <td>86</td>
      <td>2.93</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>329</td>
      <td>99</td>
      <td>3.32</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">RaTG13</th>
      <th>lib1</th>
      <td>328</td>
      <td>114</td>
      <td>2.88</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>489</td>
      <td>123</td>
      <td>3.98</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">Rf1</th>
      <th>lib1</th>
      <td>298</td>
      <td>110</td>
      <td>2.71</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>404</td>
      <td>123</td>
      <td>3.28</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">Rp3</th>
      <th>lib1</th>
      <td>250</td>
      <td>98</td>
      <td>2.55</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>388</td>
      <td>112</td>
      <td>3.46</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SARS-CoV</th>
      <th>lib1</th>
      <td>371</td>
      <td>126</td>
      <td>2.94</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>490</td>
      <td>131</td>
      <td>3.74</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SARS-CoV-2</th>
      <th>lib1</th>
      <td>425538</td>
      <td>104240</td>
      <td>4.08</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>566406</td>
      <td>102987</td>
      <td>5.50</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WIV16</th>
      <th>lib1</th>
      <td>305</td>
      <td>105</td>
      <td>2.90</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>423</td>
      <td>122</td>
      <td>3.47</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">ZC45</th>
      <th>lib1</th>
      <td>323</td>
      <td>126</td>
      <td>2.56</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>468</td>
      <td>143</td>
      <td>3.27</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">ZXC21</th>
      <th>lib1</th>
      <td>265</td>
      <td>107</td>
      <td>2.48</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>360</td>
      <td>114</td>
      <td>3.16</td>
    </tr>
  </tbody>
</table>


### Filter processed CCSs
We have the PacBio `ccs` program's estimated "accuracy" for both the barcode and the gene sequence for each processed CCS.
We will filter the CCSs to only keep ones of sufficiently high accuracy.

First, we want to plot the accuracies.
It is actually visually easier to look at the error rate, which is one minus the accuracy.
Because we want to plot on a log scale (which can't show error rates of zero), we define a *error_rate_floor*, and set all error rates less than this to that value:


```python
error_rate_floor = 1e-7  # error rates < this set to this
if error_rate_floor >= config['max_error_rate']:
    raise ValueError('error_rate_floor must be < max_error_rate')

processed_ccs = (
    processed_ccs
    .assign(barcode_error=lambda x: numpy.clip(1 - x['barcode_accuracy'],
                                               error_rate_floor, None),
            gene_error=lambda x: numpy.clip(1 - x['gene_accuracy'],
                                            error_rate_floor, None)
            )
    )
```

Now plot the error rates, drawing a dashed vertical line at the threshold separating the CCSs we retain for consensus building versus those that we discard:


```python
_ = (
 ggplot(processed_ccs
        .melt(value_vars=['barcode_error', 'gene_error'],
              var_name='feature_type', value_name='error rate'),
        aes('error rate')) +
 geom_histogram(bins=25) +
 geom_vline(xintercept=config['max_error_rate'],
            linetype='dashed',
            color=CBPALETTE[1]) +
 facet_wrap('~ feature_type') +
 theme(figure_size=(4.5, 2)) +
 ylab('number of CCSs') +
 scale_x_log10()
 ).draw()
```


![png](build_variants_files/build_variants_22_0.png)


Flag the CCSs to retain, and indicate how many we are retaining and purging due to the accuracy filter:


```python
processed_ccs = (
    processed_ccs
    .assign(retained=lambda x: ((x['gene_error'] < config['max_error_rate']) &
                                (x['barcode_error'] < config['max_error_rate'])))
    )
```

Here are number of retained CCSs:


```python
_ = (
 ggplot(processed_ccs.assign(xlabel=lambda x: x['target'] + ', ' + x['library'])
                     .groupby(['xlabel', 'retained'])
                     .size()
                     .rename('count')
                     .reset_index(),
        aes('xlabel', 'count', color='retained', label='count')) +
 geom_point(size=3) +
 geom_text(va='bottom', size=7, ha='center',format_string='{:.3g}', nudge_y=0.2) +
 theme(figure_size=(0.5 * nlibs * ntargets, 3),
       panel_grid_major_x=element_blank(),
       axis_text_x=element_text(angle=90),
       ) +
 scale_y_log10(name='number of CCSs') +
 xlab('') +
 scale_color_manual(values=CBPALETTE[1:])
 ).draw()
```


![png](build_variants_files/build_variants_26_0.png)


### Sequences per barcode
How many times is each barcode sequenced?
This is useful to know for thinking about building the barcode consensus.

First, plot the distribution of the number of times each **barcode** is observed among the retained CCSs:


```python
max_count = 8 # in plot, group all barcodes with >= this many counts

p = (
 ggplot(
    processed_ccs
     .query('retained')
     .groupby(['library', 'barcode'])
     .size()
     .rename('nseqs')
     .reset_index()
     .assign(nseqs=lambda x: numpy.clip(x['nseqs'], None, max_count)),
    aes('nseqs')) +
 geom_bar() +
 facet_wrap('~ library', nrow=1) +
 theme(figure_size=(1.75 * nlibs, 2),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of barcodes') +
 xlab('CCSs for barcode')
 )

_ = p.draw()
```


![png](build_variants_files/build_variants_28_0.png)


Now we plot the distribution of the number of **sequences** with barcodes that are observed a given number of time (again among retained CCSs).
This plot differs from the one above because if a barcode is observed multiple times it is given one count in the barcode-count plot above (which tallies barcodes) but multiple counts in the sequence-count plot below:


```python
_ = (
 ggplot(
    processed_ccs
     .query('retained')
     .assign(barcode_counts=1)
     .assign(barcode_counts=lambda x: numpy.clip(
        x.groupby(['library', 'barcode']).transform('count'), None, max_count)),
    aes('barcode_counts')) +
 geom_bar() +
 facet_wrap('~ library', nrow=1) +
 theme(figure_size=(1.75 * nlibs, 2),
       panel_grid_major_x=element_blank(),) +
 ylab('number of CCSs') +
 xlab('CCSs per barcode')
 ).draw()
```


![png](build_variants_files/build_variants_30_0.png)


### Empirical accuracy of CCSs
We want to directly estimate the accuracy of the gene-barcode link rather than relying on the PacBio `ccs` accuracy, which doesn't include inaccuracies due to things like strand exchange or the same barcode on different sequences.

One way to do this is to examine instances when we have multiple sequences for the same barcode. 
We can calculate the empirical accuracy of the sequences by looking at all instances of multiple sequences of the same barcode and determining how often they are identical.
This calculation is performed by `alignparse.consensus.empirical_accuracy` using the equations described in the docs for that function.

We will do this four for sets of sequences:

 1. All of the CCSs retained above.
 2. CCSs retained by applying a PacBio `ccs` accuracy filter 10-fold more stringent than the one above.
    The rationale is that if this improves the concordance (real accuracy) of the CCSs substantially then maybe we should make the accuracy filter more stringent.
 3. Like (1) but excluding all CCSs with indels.
    the rationale is that we only really care about substitutions, and will exclude sequences with indels anyway.
 4. Like (2) but excluding all CCSs with indels.
 
First, we annotate the sequences with the number of indels and whether they have an indel to enable categorization into the aforementioned sets:


```python
processed_ccs = alignparse.consensus.add_mut_info_cols(processed_ccs,
                                                       mutation_col='gene_mutations',
                                                       n_indel_col='n_indels')

processed_ccs = processed_ccs.assign(has_indel=lambda x: x['n_indels'] > 0)
```

Plot how many sequences have indels:


```python
_ = (
 ggplot(processed_ccs,
        aes('retained', fill='has_indel')) +
 geom_bar(position='dodge') +
 geom_text(aes(label='..count..'), stat='count', va='bottom', size=7,
           position=position_dodge(width=0.9), format_string='{:.2g}') +
 theme(figure_size=(2.5 * nlibs, 3),
       panel_grid_major_x=element_blank(),
       ) +
 ylab('number of CCSs') +
 scale_fill_manual(values=CBPALETTE[1:]) +
 facet_wrap('~ library', nrow=1)
 ).draw()
```


![png](build_variants_files/build_variants_34_0.png)


Now get the empirical accuracy for each of the CCS groups mentioned above:


```python
high_acc = config['max_error_rate'] / 10
empirical_acc = []

for desc, query_str in [
        ('retained', 'retained'),
        ('retained, no indel', 'retained and not has_indel'),
        ('10X accuracy',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc})"),
        ('10X accuracy, no indel',
         f"(gene_error < {high_acc}) and (barcode_error < {high_acc}) and not has_indel")
        ]:
    # get just CCSs in that category
    df = processed_ccs.query(query_str)
    
    # compute empirical accuracy
    empirical_acc.append(
        alignparse.consensus.empirical_accuracy(df,
                                                mutation_col='gene_mutations')
        .assign(description=desc)
        .merge(df
               .groupby('library')
               .size()
               .rename('number_CCSs')
               .reset_index()
               )
        )

# make description categorical to preserve order, and annotate as "actual"
# the category ("retained, no indel") that we will use for building variants.
empirical_acc = (
    pd.concat(empirical_acc, ignore_index=True, sort=False)
    .assign(description=lambda x: pd.Categorical(x['description'],
                                                 x['description'].unique(),
                                                 ordered=True),
            actual=lambda x: numpy.where(x['description'] == 'retained, no indel',
                                         True, False),
            )
    )
```

Display table of the empirical accuracies:


```python
display(HTML(empirical_acc.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>accuracy</th>
      <th>description</th>
      <th>number_CCSs</th>
      <th>actual</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib1</td>
      <td>0.984716</td>
      <td>retained</td>
      <td>408296</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>0.984339</td>
      <td>retained</td>
      <td>531571</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>0.998122</td>
      <td>retained, no indel</td>
      <td>390947</td>
      <td>True</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>0.998340</td>
      <td>retained, no indel</td>
      <td>509301</td>
      <td>True</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>0.989409</td>
      <td>10X accuracy</td>
      <td>395183</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>0.990143</td>
      <td>10X accuracy</td>
      <td>507336</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>0.998138</td>
      <td>10X accuracy, no indel</td>
      <td>380115</td>
      <td>False</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>0.998347</td>
      <td>10X accuracy, no indel</td>
      <td>488936</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


Plot the empirical accuracies, using a different color to show the category that we will actually use:


```python
p = (
    ggplot(empirical_acc,
           aes('description', 'accuracy', color='actual', label='accuracy')
           ) +
    geom_point(size=3) +
    geom_text(va='bottom', size=9, format_string='{:.3g}', nudge_y=0.003) +
    facet_wrap('~ library') +
    theme(figure_size=(1.75 * nlibs, 2.25),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank(),
          ) +
    xlab('') +
    scale_y_continuous(name='empirical accuracy', limits=(0.95, 1.005)) +
    scale_color_manual(values=CBPALETTE, guide=False)
    )

plotfile = os.path.join(config['figs_dir'], 'empirical_CCS_accuracy.pdf')
print(f"Saving plot to {plotfile}")
_ = p.draw()
```

    Saving plot to results/figures/empirical_CCS_accuracy.pdf



![png](build_variants_files/build_variants_40_1.png)


The above analysis shows that if we exclude sequences with indels (which we plan to do among our consensus sequences), then the accuracy of each CCS is around 99%. 
We do **not** get notably higher empirical accuracy by imposing a more stringent filter from the PacBio `ccs` program, indicating that the major sources of error are due to processes that are not modeled in this program's accuracy filter (perhaps strand exchange or barcode sharing).

Note that this empirical accuracy is for a **single** CCS.
When we build the consensus sequences for each barcode below, we will take the consensus of CCSs within a barcode.
So for barcodes with multiple CCSs, the actual accuracy of the consensus sequences will be higher than the empirical accuracy above due to capturing information from multiple CCSs.

### Consensus sequences for barcodes
We call the consensus sequence for each barcode using the simple method implemented in [alignparse.consensus.simple_mutconsensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html?highlight=simple_mutconsensus#alignparse.consensus.simple_mutconsensus).
The documentation for that function explains the method in detail, but basically it works like this:
 1. When there is just one CCS per barcode, the consensus is just that sequence.
 2. When there are multiple CCSs per barcode, they are used to build a consensus--however, the entire barcode is discarded if there are many differences between CCSs with the barcode, or high-frequency non-consensus mutations. The reason that barcodes are discarded in such cases as many differences between CCSs or high-frequency non-consensus mutations suggest errors such as barcode collisions or strand exchange.
 
First, call the consensus for each barcode including **all** retained sequences, even those with indels:


```python
consensus, dropped = alignparse.consensus.simple_mutconsensus(
                        processed_ccs.query('retained'),
                        group_cols=('library', 'barcode', 'target'),
                        mutation_col='gene_mutations',
                        )
```

Here are the first few lines of the data frame of consensus sequences for each barcode.
In addition to giving the library, barcode, target, and mutations, it also has a column indicating how many CCSs support the variant call:


```python
display(HTML(consensus.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>gene_mutations</th>
      <th>variant_call_support</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib1</td>
      <td>AAAAAAAAAATTGTAA</td>
      <td>SARS-CoV-2</td>
      <td>T271C A272T T273C A595T A597C</td>
      <td>1</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAAAAAACTTAAAT</td>
      <td>SARS-CoV-2</td>
      <td>A38G T39C T179C A180G G282C T439A G449A</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAAAAAACTTCAAT</td>
      <td>SARS-CoV-2</td>
      <td>G64T C65G A66C C379G T381C A423C C562G G564C</td>
      <td>5</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAAAAACAAGCAGA</td>
      <td>SARS-CoV-2</td>
      <td>A37T A38T T39C</td>
      <td>6</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAAAAACAATATAA</td>
      <td>SARS-CoV-2</td>
      <td>T16A G17A T18G A43T C44G T45G A280T G282C G307T T308G</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Since we retain variants with substitutions but ignore those with indels, add information about substitution mutations and number of indels:


```python
consensus = alignparse.consensus.add_mut_info_cols(
                    consensus,
                    mutation_col='gene_mutations',
                    sub_str_col='substitutions',
                    n_indel_col='number_of_indels',
                    overwrite_cols=True)

display(HTML(consensus.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>gene_mutations</th>
      <th>variant_call_support</th>
      <th>substitutions</th>
      <th>number_of_indels</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib1</td>
      <td>AAAAAAAAAATTGTAA</td>
      <td>SARS-CoV-2</td>
      <td>T271C A272T T273C A595T A597C</td>
      <td>1</td>
      <td>T271C A272T T273C A595T A597C</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAAAAAACTTAAAT</td>
      <td>SARS-CoV-2</td>
      <td>A38G T39C T179C A180G G282C T439A G449A</td>
      <td>2</td>
      <td>A38G T39C T179C A180G G282C T439A G449A</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAAAAAACTTCAAT</td>
      <td>SARS-CoV-2</td>
      <td>G64T C65G A66C C379G T381C A423C C562G G564C</td>
      <td>5</td>
      <td>G64T C65G A66C C379G T381C A423C C562G G564C</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAAAAACAAGCAGA</td>
      <td>SARS-CoV-2</td>
      <td>A37T A38T T39C</td>
      <td>6</td>
      <td>A37T A38T T39C</td>
      <td>0</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAAAAACAATATAA</td>
      <td>SARS-CoV-2</td>
      <td>T16A G17A T18G A43T C44G T45G A280T G282C G307T T308G</td>
      <td>1</td>
      <td>T16A G17A T18G A43T C44G T45G A280T G282C G307T T308G</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


Plot distribution of number of CCSs supporting each variant call (consensus), indicating whether or not there is an indel:


```python
max_variant_call_support = 6  # group variants with >= this much support

_ = (
 ggplot(consensus
        .assign(variant_call_support=lambda x: numpy.clip(x['variant_call_support'],
                                                          None,
                                                          max_variant_call_support),
                indel_state=lambda x: numpy.where(x['number_of_indels'] > 0,
                                                  'has indel', 'no indel')
                ),
        aes('variant_call_support')) +
 geom_bar() +
 ylab('number of variants') +
 facet_grid('indel_state ~ library') +
 theme(figure_size=(1.75 * nlibs, 3.5),
       panel_grid_major_x=element_blank(),
       ) 
 ).draw()
```


![png](build_variants_files/build_variants_49_0.png)


We see that most variant consensus sequences do **not** have indels, especially if we limit to the more "accurate" ones that have multiple CCSs supporting them.

We will ignore all consensus sequences with indels in the variant-barcode lookup table. 
We do this for two reasons:
 1. When there is just one CCS supporting a consensus, it is less likely to be accurate as indels are the main mode of PacBio error.
 2. For the purposes of our studies, we are interested in point mutations rather than indels anyway.
 
Here are number of valid consensus sequence (no indels) for each library and target:


```python
consensus = consensus.query('number_of_indels == 0')

lib_target_counts = (
    consensus
    .groupby(['library', 'target'])
    .size()
    .rename('consensus sequences')
    .reset_index()
    )

display(HTML(lib_target_counts.to_html(index=False)))

p = (ggplot(lib_target_counts.assign(xlabel=lambda x: x['target'] + ', ' + x['library']),
            aes('xlabel', 'consensus sequences')) +
     geom_point(size=3) +
     theme(figure_size=(0.5 * nlibs * ntargets, 1.75),
           axis_text_x=element_text(angle=90)) +
     xlab('') +
     scale_y_log10()
     )

_ = p.draw()
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>target</th>
      <th>consensus sequences</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib1</td>
      <td>BM48-31</td>
      <td>108</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>GD-Pangolin</td>
      <td>119</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>HKU3-1</td>
      <td>106</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>LYRa11</td>
      <td>84</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>RaTG13</td>
      <td>112</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>Rf1</td>
      <td>106</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>Rp3</td>
      <td>95</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>SARS-CoV</td>
      <td>124</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>98569</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>WIV16</td>
      <td>103</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>ZC45</td>
      <td>123</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>ZXC21</td>
      <td>107</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>BM48-31</td>
      <td>118</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>GD-Pangolin</td>
      <td>127</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>HKU3-1</td>
      <td>122</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>LYRa11</td>
      <td>94</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>RaTG13</td>
      <td>120</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>Rf1</td>
      <td>121</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>Rp3</td>
      <td>110</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>SARS-CoV</td>
      <td>125</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>SARS-CoV-2</td>
      <td>96515</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>WIV16</td>
      <td>116</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>ZC45</td>
      <td>141</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>ZXC21</td>
      <td>112</td>
    </tr>
  </tbody>
</table>



![png](build_variants_files/build_variants_51_1.png)


For the non-primary targets, we want to drop all barcodes with mutations:


```python
primary_target = config['primary_target']
print(f"Dropping variants with mutations for all targets except {primary_target}")

consensus = (
    consensus
    .assign(has_substitutions=lambda x: x['substitutions'].str.len().astype(bool))
    )

has_subs_by_target = (
        consensus
        .groupby(['target', 'library', 'has_substitutions'])
        .aggregate(n_barcodes=pd.NamedAgg('barcode', 'count'))
        .reset_index()
        )

display(HTML(has_subs_by_target
             .pivot_table(index=['target', 'library'],
                          columns='has_substitutions',
                          values='n_barcodes',
                          fill_value=0)
             .to_html()))

p = (ggplot(has_subs_by_target.assign(xlabel=lambda x: x['target'] + ', ' + x['library']),
            aes('xlabel', 'n_barcodes', color='has_substitutions')) +
     geom_point(size=3, alpha=0.7) +
     theme(figure_size=(0.5 * nlibs * ntargets, 1.75),
           axis_text_x=element_text(angle=90)) +
     xlab('') +
     scale_y_log10() +
     scale_color_manual(values=CBPALETTE)
     )

_ = p.draw()
```

    Dropping variants with mutations for all targets except SARS-CoV-2



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>has_substitutions</th>
      <th>False</th>
      <th>True</th>
    </tr>
    <tr>
      <th>target</th>
      <th>library</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">BM48-31</th>
      <th>lib1</th>
      <td>107</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>117</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">GD-Pangolin</th>
      <th>lib1</th>
      <td>119</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>127</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">HKU3-1</th>
      <th>lib1</th>
      <td>106</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>122</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">LYRa11</th>
      <th>lib1</th>
      <td>35</td>
      <td>49</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>41</td>
      <td>53</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">RaTG13</th>
      <th>lib1</th>
      <td>112</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>118</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">Rf1</th>
      <th>lib1</th>
      <td>55</td>
      <td>51</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>58</td>
      <td>63</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">Rp3</th>
      <th>lib1</th>
      <td>95</td>
      <td>0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>110</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SARS-CoV</th>
      <th>lib1</th>
      <td>122</td>
      <td>2</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>123</td>
      <td>2</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">SARS-CoV-2</th>
      <th>lib1</th>
      <td>4481</td>
      <td>94088</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>4643</td>
      <td>91872</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WIV16</th>
      <th>lib1</th>
      <td>102</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>115</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">ZC45</th>
      <th>lib1</th>
      <td>122</td>
      <td>1</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>140</td>
      <td>1</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">ZXC21</th>
      <th>lib1</th>
      <td>104</td>
      <td>3</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>111</td>
      <td>1</td>
    </tr>
  </tbody>
</table>



![png](build_variants_files/build_variants_53_2.png)


Print most common non-primary targets with mutations (useful for debugging):


```python
display(HTML(
    consensus
    .query('target != @primary_target')
    .query('has_substitutions')
    .groupby(['target', 'library', 'substitutions'])
    .aggregate(n_barcodes=pd.NamedAgg('barcode', 'count'))
    .reset_index()
    .sort_values('n_barcodes', ascending=False)
    .head(n=20)
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>substitutions</th>
      <th>n_barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Rf1</td>
      <td>lib2</td>
      <td>T506A</td>
      <td>60</td>
    </tr>
    <tr>
      <td>LYRa11</td>
      <td>lib2</td>
      <td>G64A</td>
      <td>52</td>
    </tr>
    <tr>
      <td>LYRa11</td>
      <td>lib1</td>
      <td>G64A</td>
      <td>48</td>
    </tr>
    <tr>
      <td>Rf1</td>
      <td>lib1</td>
      <td>T506A</td>
      <td>48</td>
    </tr>
    <tr>
      <td>ZXC21</td>
      <td>lib1</td>
      <td>T491C</td>
      <td>2</td>
    </tr>
    <tr>
      <td>RaTG13</td>
      <td>lib2</td>
      <td>G343T</td>
      <td>2</td>
    </tr>
    <tr>
      <td>BM48-31</td>
      <td>lib1</td>
      <td>T355C</td>
      <td>1</td>
    </tr>
    <tr>
      <td>SARS-CoV</td>
      <td>lib1</td>
      <td>T303C T384C</td>
      <td>1</td>
    </tr>
    <tr>
      <td>ZXC21</td>
      <td>lib1</td>
      <td>C299T</td>
      <td>1</td>
    </tr>
    <tr>
      <td>ZC45</td>
      <td>lib2</td>
      <td>C225T</td>
      <td>1</td>
    </tr>
    <tr>
      <td>ZC45</td>
      <td>lib1</td>
      <td>C225T</td>
      <td>1</td>
    </tr>
    <tr>
      <td>WIV16</td>
      <td>lib2</td>
      <td>C98T</td>
      <td>1</td>
    </tr>
    <tr>
      <td>WIV16</td>
      <td>lib1</td>
      <td>C98T</td>
      <td>1</td>
    </tr>
    <tr>
      <td>SARS-CoV</td>
      <td>lib2</td>
      <td>T303C T384C</td>
      <td>1</td>
    </tr>
    <tr>
      <td>SARS-CoV</td>
      <td>lib2</td>
      <td>C8T</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Rf1</td>
      <td>lib2</td>
      <td>C457T</td>
      <td>1</td>
    </tr>
    <tr>
      <td>SARS-CoV</td>
      <td>lib1</td>
      <td>C8T</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BM48-31</td>
      <td>lib2</td>
      <td>T355C</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Rf1</td>
      <td>lib2</td>
      <td>C300T</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Rf1</td>
      <td>lib2</td>
      <td>C165T</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Remove the non-primary targets with any substitutions:


```python
print(f"Culling the {len(consensus)} barcodes to remove mutated non-primary targets")

consensus = consensus.query('(target == @primary_target) or (has_substitutions == False)')

print(f"Retained {len(consensus)} barcodes after culling")
```

    Culling the 197577 barcodes to remove mutated non-primary targets
    Retained 197345 barcodes after culling


Are there any barcodes in the same library that are shared across targets?
If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    consensus
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(consensus)} barcodes:")
consensus = (
    consensus
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(consensus)} barcodes.")
```

    Here are duplicated barcodes:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib2</td>
      <td>AACTGAGTTAACGTCC</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>ATAACCAAAACATCCA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>lib2</td>
      <td>CTCAAGGTAACGATCC</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


    
    Removing the 3 duplicated barcodes.Started with 197345 barcodes:
    After removing duplicates, there are 197339 barcodes.


Below we write the retained consensus sequences to a CSV file that links the nucleotide mutations to the barcodes.
(The next section analyzes this variant table in detail, and provides have more precise information on the number of variants and relevant statistics):


```python
print(f"Writing nucleotide variants to {config['nt_variant_table_file']}")
      
(consensus
 [['target', 'library', 'barcode', 'substitutions', 'variant_call_support']]
 .to_csv(config['nt_variant_table_file'], index=False)
 )
      
print('Here are the first few lines of this file:')
display(HTML(
    pd.read_csv(config['nt_variant_table_file'], na_filter=None)
    .head()
    .to_html(index=False)
    ))
```

    Writing nucleotide variants to results/variants/nucleotide_variant_table.csv
    Here are the first few lines of this file:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>substitutions</th>
      <th>variant_call_support</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>SARS-CoV-2</td>
      <td>lib1</td>
      <td>AAAAAAAAAATTGTAA</td>
      <td>T271C A272T T273C A595T A597C</td>
      <td>1</td>
    </tr>
    <tr>
      <td>SARS-CoV-2</td>
      <td>lib1</td>
      <td>AAAAAAAAACTTAAAT</td>
      <td>A38G T39C T179C A180G G282C T439A G449A</td>
      <td>2</td>
    </tr>
    <tr>
      <td>SARS-CoV-2</td>
      <td>lib1</td>
      <td>AAAAAAAAACTTCAAT</td>
      <td>G64T C65G A66C C379G T381C A423C C562G G564C</td>
      <td>5</td>
    </tr>
    <tr>
      <td>SARS-CoV-2</td>
      <td>lib1</td>
      <td>AAAAAAAACAAGCAGA</td>
      <td>A37T A38T T39C</td>
      <td>6</td>
    </tr>
    <tr>
      <td>SARS-CoV-2</td>
      <td>lib1</td>
      <td>AAAAAAAACAATATAA</td>
      <td>T16A G17A T18G A43T C44G T45G A280T G282C G307T T308G</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


What happened to the barcodes that we "dropped" because we could not construct a reliable consensus?
The `dropped` data frame from [alignparse.consensus.simple_mutconsensus](https://jbloomlab.github.io/alignparse/alignparse.consensus.html?highlight=simple_mutconsensus#alignparse.consensus.simple_mutconsensus) has this information:


```python
display(HTML(dropped.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>target</th>
      <th>drop_reason</th>
      <th>nseqs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib1</td>
      <td>AAAAAAAAATGATTTC</td>
      <td>SARS-CoV-2</td>
      <td>subs diff too large</td>
      <td>6</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAAACCTTACAGAA</td>
      <td>SARS-CoV-2</td>
      <td>subs diff too large</td>
      <td>7</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAAATCCGACCCTG</td>
      <td>SARS-CoV-2</td>
      <td>minor subs too frequent</td>
      <td>11</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAAATTCCAGAGGC</td>
      <td>SARS-CoV-2</td>
      <td>minor subs too frequent</td>
      <td>6</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>AAAAACACGCCTGTAA</td>
      <td>SARS-CoV-2</td>
      <td>minor subs too frequent</td>
      <td>5</td>
    </tr>
  </tbody>
</table>


Summarize the information in this data frame on dropped barcodes with the plot below.
This plot shows several things.
First, we see that the total number of barcodes dropped is modest (just a few thousand per library) relative to the total number of barcodes per library (seen above to be on the order of hundreds of thousands).
Second, the main reason that barcodes are dropped is that there are CCSs within the same barcode with suspiciously large numbers of mutations relative to the consensus---which we use as a filter to discard the entire barcode as it could indicate strand exchange or some other issue.
In any case, the modest number of dropped barcodes indicates that there probably isn't much of a need to worry: 


```python
max_nseqs = 8  # plot together all barcodes with >= this many sequences

_ = (
 ggplot(
    dropped.assign(nseqs=lambda x: numpy.clip(x['nseqs'], None, max_nseqs)),
    aes('nseqs')) + 
 geom_bar() + 
 scale_x_continuous(limits=(1, None)) +
 xlab('number of sequences for barcode') +
 ylab('number of barcodes') +
 facet_grid('library ~ drop_reason') +
 theme(figure_size=(10, 1.5 * nlibs),
       panel_grid_major_x=element_blank(),
       )
 ).draw()
```


![png](build_variants_files/build_variants_65_0.png)


## Create barcode-variant table
We now create a [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) that stores and processes all the information about the variant consensus sequences.
Below we initialize such a table, and then analyze information about its composition.

### Initialize codon variant table
In order to initialize the codon variant table, we need two pieces of information:
  1. The wildtype gene sequence.
  2. The list of nucleotide mutations for each variant as determined in the consensus calling above.

Read "wildtype" gene sequence to which we made the alignments (in order to do this, initialize an `alignparse.Targets` and get the gene sequence from it):


```python
targets = alignparse.targets.Targets(seqsfile=config['amplicons'],
                                     feature_parse_specs=config['feature_parse_specs'])
geneseq = targets.get_target(config['primary_target']).get_feature('gene').seq

print(f"Read gene of {len(geneseq)} nts for {config['primary_target']} from {config['amplicons']}")
```

    Read gene of 603 nts for SARS-CoV-2 from data/PacBio_amplicons.gb


Now initialize the codon variant table using this wildtype sequence and our list of nucleotide mutations for each variant:


```python
variants = dms_variants.codonvarianttable.CodonVariantTable(
                barcode_variant_file=config['nt_variant_table_file'],
                geneseq=geneseq,
                primary_target=config['primary_target'],
                )
```

### Basic stats on variants
We now will analyze the variants.
In this call and in the plots below, we set `samples=None` as we aren't looking at variant counts in specific samples, but are simply looking at properties of the variants in the table.

Here are the number of variants for each target:


```python
display(HTML(
    variants
    .n_variants_df(samples=None)
    .pivot_table(index=['target'],
                 columns='library',
                 values='count')
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>lib1</th>
      <th>lib2</th>
      <th>all libraries</th>
    </tr>
    <tr>
      <th>target</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>BM48-31</th>
      <td>107</td>
      <td>117</td>
      <td>224</td>
    </tr>
    <tr>
      <th>GD-Pangolin</th>
      <td>119</td>
      <td>127</td>
      <td>246</td>
    </tr>
    <tr>
      <th>HKU3-1</th>
      <td>106</td>
      <td>121</td>
      <td>227</td>
    </tr>
    <tr>
      <th>LYRa11</th>
      <td>35</td>
      <td>40</td>
      <td>75</td>
    </tr>
    <tr>
      <th>RaTG13</th>
      <td>112</td>
      <td>118</td>
      <td>230</td>
    </tr>
    <tr>
      <th>Rf1</th>
      <td>55</td>
      <td>58</td>
      <td>113</td>
    </tr>
    <tr>
      <th>Rp3</th>
      <td>95</td>
      <td>110</td>
      <td>205</td>
    </tr>
    <tr>
      <th>SARS-CoV</th>
      <td>122</td>
      <td>123</td>
      <td>245</td>
    </tr>
    <tr>
      <th>SARS-CoV-2</th>
      <td>98569</td>
      <td>96512</td>
      <td>195081</td>
    </tr>
    <tr>
      <th>WIV16</th>
      <td>102</td>
      <td>115</td>
      <td>217</td>
    </tr>
    <tr>
      <th>ZC45</th>
      <td>122</td>
      <td>140</td>
      <td>262</td>
    </tr>
    <tr>
      <th>ZXC21</th>
      <td>104</td>
      <td>110</td>
      <td>214</td>
    </tr>
  </tbody>
</table>


Plot the number of variants supported by each number of CCSs:


```python
max_support = 10  # group variants with >= this much support

p = variants.plotVariantSupportHistogram(max_support=max_support,
                                         widthscale=1.1,
                                         heightscale=0.9)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
```


![png](build_variants_files/build_variants_73_0.png)


### Mutations per variant
Plot the number of barcoded variants with each number of amino-acid and codon mutations.
This is for the primary target only, and doesn't include the spiked-in secondary targets:


```python
max_muts = 7  # group all variants with >= this many mutations

for mut_type in ['aa', 'codon']:
    p = variants.plotNumMutsHistogram(mut_type, samples=None, max_muts=max_muts,
                                      widthscale=1.1,
                                      heightscale=0.9)
    p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
    _ = p.draw()
    plotfile = os.path.join(config['figs_dir'], f"n_{mut_type}_muts_per_variant.pdf")
    print(f"Saving plot to {plotfile}")
    p.save(plotfile)
```

    Saving plot to results/figures/n_aa_muts_per_variant.pdf
    Saving plot to results/figures/n_codon_muts_per_variant.pdf



![png](build_variants_files/build_variants_75_1.png)



![png](build_variants_files/build_variants_75_2.png)


Plot the frequencies of different codon mutation types among **all** variants (any number of mutations), again only for primary target:


```python
p = variants.plotNumCodonMutsByType(variant_type='all', samples=None,
                                    ylabel='mutations per variant',
                                    heightscale=0.8)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
plotfile = os.path.join(config['figs_dir'], f"avg_muts_per_variant.pdf")
print(f"Saving plot to {plotfile}")
p.save(plotfile)
```

    Saving plot to results/figures/avg_muts_per_variant.pdf



![png](build_variants_files/build_variants_77_1.png)


Variants supported by multiple PacBio CCSs should have fewer spurious mutations since sequencing errors are very unlikely to occur on two CCSs.
Below we plot the number of codon mutations per variant among variants with at least two CCSs supporting their call.
The difference in mutation rates here and in the plot above (that does not apply the `min_support=2` filter) gives some estimate of the frequency of mutations in our variants our spurious.
In fact, we see the numbers are very similar, indicating that few of the mutations are spurious:


```python
p = variants.plotNumCodonMutsByType(variant_type='all', samples=None,
                                    ylabel='mutations per variant', 
                                    min_support=2, heightscale=0.8)
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
```


![png](build_variants_files/build_variants_79_0.png)


### Completeness of mutation sampling
We examine how completely amino-acid mutations are sampled by the variants for the primary target, looking at single-mutant variants only and all variants.
The plot below shows that virtually every mutation is found in a variant in each library, even if we just look among the single mutants.
Things look especially good if we aggregate across libraries:


```python
for variant_type in ['all', 'single']:
    p = variants.plotCumulMutCoverage(variant_type, mut_type='aa', samples=None)
    _ = p.draw()
    plotfile = os.path.join(config['figs_dir'],
                            f"variant_cumul_{variant_type}_mut_coverage.pdf")
    print(f"Saving plot to {plotfile}")
    p.save(plotfile)
```

    Saving plot to results/figures/variant_cumul_all_mut_coverage.pdf
    Saving plot to results/figures/variant_cumul_single_mut_coverage.pdf



![png](build_variants_files/build_variants_81_1.png)



![png](build_variants_files/build_variants_81_2.png)


To get more quantitative information like that plotted above, we determine how many mutations are found 0, 1, or >1 times both among single and all mutants for the primary target:


```python
count_dfs = []
for variant_type in ['all', 'single']:
    count_dfs.append(
        variants.mutCounts(variant_type, mut_type='aa', samples=None)
        .assign(variant_type=variant_type)
        )
    
display(HTML(
    pd.concat(count_dfs)
    .assign(count=lambda x: (numpy.clip(x['count'], None, 2)
                             .map({0: '0', 1: '1', 2:'>1'}))
            )
    .groupby(['variant_type', 'library', 'count'])
    .aggregate(number_of_mutations=pd.NamedAgg(column='mutation', aggfunc='count'))
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th></th>
      <th>number_of_mutations</th>
    </tr>
    <tr>
      <th>variant_type</th>
      <th>library</th>
      <th>count</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="9" valign="top">all</th>
      <th rowspan="3" valign="top">lib1</th>
      <th>0</th>
      <td>24</td>
    </tr>
    <tr>
      <th>1</th>
      <td>15</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3981</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">lib2</th>
      <th>0</th>
      <td>24</td>
    </tr>
    <tr>
      <th>1</th>
      <td>17</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3979</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">all libraries</th>
      <th>0</th>
      <td>16</td>
    </tr>
    <tr>
      <th>1</th>
      <td>10</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3994</td>
    </tr>
    <tr>
      <th rowspan="9" valign="top">single</th>
      <th rowspan="3" valign="top">lib1</th>
      <th>0</th>
      <td>507</td>
    </tr>
    <tr>
      <th>1</th>
      <td>609</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>2904</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">lib2</th>
      <th>0</th>
      <td>491</td>
    </tr>
    <tr>
      <th>1</th>
      <td>572</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>2957</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">all libraries</th>
      <th>0</th>
      <td>189</td>
    </tr>
    <tr>
      <th>1</th>
      <td>254</td>
    </tr>
    <tr>
      <th>&gt;1</th>
      <td>3577</td>
    </tr>
  </tbody>
</table>


### Mutation frequencies along gene
We plot the frequencies of mutations along the gene among the variants for the primary target.
Ideally, this would be uniform.
We make the plot for both all variants and single-mutant / wildtype variants:


```python
for variant_type in ['all', 'single']:
    p = variants.plotMutFreqs(variant_type, mut_type='codon', samples=None)
    p.draw()
```


![png](build_variants_files/build_variants_85_0.png)



![png](build_variants_files/build_variants_85_1.png)


We can also use heat maps to examine the extent to which specific amino-acid or codon mutations are over-represented.
These heat maps are large, so we make them just for all variants and the merge of all libraries:


```python
for mut_type in ['aa', 'codon']:
    p = variants.plotMutHeatmap('all', mut_type, samples=None, #libraries='all_only',
                                widthscale=2)
    p.draw()
```


![png](build_variants_files/build_variants_87_0.png)



![png](build_variants_files/build_variants_87_1.png)


### Write codon-variant table
We write the codon variant table to a CSV file.
This table looks like this:


```python
display(HTML(
    variants.barcode_variant_df
    .head()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>BM48-31</td>
      <td>lib1</td>
      <td>AAAACTAGATGAATAA</td>
      <td>3</td>
      <td>BM48-31</td>
      <td>BM48-31</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>BM48-31</td>
      <td>lib1</td>
      <td>AAACCCTTTCAAGTAG</td>
      <td>1</td>
      <td>BM48-31</td>
      <td>BM48-31</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>BM48-31</td>
      <td>lib1</td>
      <td>AAAGAAGTAATAATTG</td>
      <td>2</td>
      <td>BM48-31</td>
      <td>BM48-31</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>BM48-31</td>
      <td>lib1</td>
      <td>AAAGGGGCGGACGGAC</td>
      <td>1</td>
      <td>BM48-31</td>
      <td>BM48-31</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>BM48-31</td>
      <td>lib1</td>
      <td>AACATTCTACCGACTG</td>
      <td>3</td>
      <td>BM48-31</td>
      <td>BM48-31</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


Note how this table differs from the nucleotide variant table we generated above and used to initialize the [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) in that it gives **codon** substitutions and associated amino-acid substitutions.

Write it to CSV file:


```python
print(f"Writing codon-variant table to {config['codon_variant_table_file']}")

variants.barcode_variant_df.to_csv(config['codon_variant_table_file'], index=False)
```

    Writing codon-variant table to results/variants/codon_variant_table.csv

