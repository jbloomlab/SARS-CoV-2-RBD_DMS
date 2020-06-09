# Logo plots showing effects of mutations
Creates logo plot visualizations and input file for [dms-view](https://dms-view.github.io/) visualizations.

## Imports and read data
Import Python modules:


```python
import math
import os

from dms_variants.constants import CBPALETTE

import matplotlib
import matplotlib.pyplot as plt

import numpy

import pandas as pd

from plotnine import *

import yaml

import logomaker_utils  # custom module that uses logomaker to create plots
```

    Matplotlib is building the font cache using fc-list. This may take a moment.


Get the Uni Sans (all caps) font set up.
Needed if we want to make plots where same letter is coded twice (for above and below line) and we want them to look like capital letters both thimes:


```python
font_files = matplotlib.font_manager.findSystemFonts('data/uni-sans')
font_list = matplotlib.font_manager.createFontList(font_files)
matplotlib.font_manager.fontManager.ttflist.extend(font_list)
```

    /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARS-CoV-2-RBD_DMS/env/lib/python3.7/site-packages/ipykernel_launcher.py:2: MatplotlibDeprecationWarning: 
    The createFontList function was deprecated in Matplotlib 3.2 and will be removed two minor releases later. Use FontManager.addfont instead.
      


Set `plotnine` theme:


```python
theme_set(theme_classic())
```

Read in the configuration file, and then read the input data files from that:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)

print(f"Reading single-mutant effects from {config['single_mut_effects_file']}")
mut_effects = pd.read_csv(config['single_mut_effects_file'])
mut_effects.head()
```

    Reading single-mutant effects from results/single_mut_effects/single_mut_effects.csv





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
      <th>site_RBD</th>
      <th>site_SARS2</th>
      <th>wildtype</th>
      <th>mutant</th>
      <th>mutation</th>
      <th>mutation_RBD</th>
      <th>bind_lib1</th>
      <th>bind_lib2</th>
      <th>bind_avg</th>
      <th>expr_lib1</th>
      <th>expr_lib2</th>
      <th>expr_avg</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>N331A</td>
      <td>N1A</td>
      <td>-0.05</td>
      <td>-0.02</td>
      <td>-0.03</td>
      <td>-0.14</td>
      <td>-0.08</td>
      <td>-0.11</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>C</td>
      <td>N331C</td>
      <td>N1C</td>
      <td>-0.08</td>
      <td>-0.10</td>
      <td>-0.09</td>
      <td>-1.56</td>
      <td>-0.97</td>
      <td>-1.26</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>N331D</td>
      <td>N1D</td>
      <td>0.00</td>
      <td>0.07</td>
      <td>0.03</td>
      <td>-0.75</td>
      <td>-0.12</td>
      <td>-0.44</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>N331E</td>
      <td>N1E</td>
      <td>0.02</td>
      <td>-0.02</td>
      <td>0.00</td>
      <td>-0.39</td>
      <td>-0.24</td>
      <td>-0.31</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>N331F</td>
      <td>N1F</td>
      <td>-0.03</td>
      <td>-0.16</td>
      <td>-0.10</td>
      <td>-0.83</td>
      <td>-0.57</td>
      <td>-0.70</td>
    </tr>
  </tbody>
</table>
</div>



Read annnotations of various properties of sites:


```python
site_annotations = pd.read_csv('data/RBD_sites.csv')

# make sure consistent with site and wildtype data in mut_effects data frame
pd.testing.assert_frame_equal(
        mut_effects[['site_RBD', 'site_SARS2', 'wildtype']].drop_duplicates().reset_index(drop=True),
        site_annotations[['site_RBD', 'site_SARS2', 'amino_acid_SARS2']].rename(columns={'amino_acid_SARS2': 'wildtype'}),
        )

# first few lines
site_annotations.head()
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
      <th>site_RBD</th>
      <th>amino_acid_SARS2</th>
      <th>site_SARS2</th>
      <th>amino_acid_SARS1</th>
      <th>site_SARS1</th>
      <th>chain_6M0J</th>
      <th>codon_SARS2</th>
      <th>amino_acid_RaTG13</th>
      <th>amino_acid_GD_Pangolin</th>
      <th>entropy</th>
      <th>...</th>
      <th>SARS1_key_adaptation</th>
      <th>epitope_CR3022</th>
      <th>epitope_VHH72</th>
      <th>epitope_S230</th>
      <th>epitope_m396</th>
      <th>epitope_F26G19</th>
      <th>epitope_80R</th>
      <th>epitope_B38</th>
      <th>epitope_S309</th>
      <th>buried_downRBD</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>N</td>
      <td>331</td>
      <td>N</td>
      <td>318.0</td>
      <td>NaN</td>
      <td>aat</td>
      <td>N</td>
      <td>N</td>
      <td>0.000000</td>
      <td>...</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>I</td>
      <td>332</td>
      <td>I</td>
      <td>319.0</td>
      <td>NaN</td>
      <td>att</td>
      <td>I</td>
      <td>I</td>
      <td>0.000000</td>
      <td>...</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>T</td>
      <td>333</td>
      <td>T</td>
      <td>320.0</td>
      <td>E</td>
      <td>aca</td>
      <td>T</td>
      <td>T</td>
      <td>0.000000</td>
      <td>...</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>N</td>
      <td>334</td>
      <td>N</td>
      <td>321.0</td>
      <td>E</td>
      <td>aac</td>
      <td>N</td>
      <td>N</td>
      <td>0.179256</td>
      <td>...</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>L</td>
      <td>335</td>
      <td>L</td>
      <td>322.0</td>
      <td>E</td>
      <td>ttg</td>
      <td>L</td>
      <td>L</td>
      <td>1.468456</td>
      <td>...</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 31 columns</p>
</div>



Make output directory:


```python
os.makedirs(config['logoplots_dir'], exist_ok=True)
```

## Compute letter heights
Here we explain how we use the experimental measurements to compute letter heights for plotting in the logo plots.

### Letter heights as "probabilities" of observing amino acids
For each site $r$, we first compute the probability $p_{r,a}$ of observing amino acid $a$ under some notion that the probability of observing an amino acid relates to its favorability on expression or ACE2-binding as measured in the experiments.
We compute $p_{r,a}$ over all non-stop amino acids observed at a site, so that $1 = \sum_a p_{r,a}$
Typically there are 20 amino acids observed at a site, although in some cases a few may not be observed.
If an amino acid $a$ is not observed at a site $r$, then we set $p_{r,a} = 0$ under the assumption that non-observed amino acids are generally deleterious (particularly for binding, where the non-observed ones are generally missing because they were eliminated in the pre-sort for expression).

One way to plot the data is to then just make the height of each letter at each site equal to $p_{r,a}$.
In this case, the stacks at all sites sum to one.

### Letter heights as "information content"
After computing the probabilities $p_{r,a}$, we can determine the "information content" at a site, which is how much these probabilities deviate from all being $\frac{1}{20}$ as would occur if there was no constraint at a height (note that there is no need to correct for background amino-acid content as is done for standard information-content sequence logos, as that was done when analyzing the experimental data to get the mutational effects).

The information content $I_r$ at site $r$ is just the total possible entropy if all amino acids were equally probable minus the actual entropy, so $I_r = \ln 20 + \sum_a p_{r,a} \ln p_{r,a}$.

The height of each letter in the information content scheme is then simply the information content times the probability of that amino acid, so $I_r \times p_{r,a}$.
Note that these letter heights are analogous to those in classic sequence logos used to show transcription-factor binding motifs.

### Computing probabilities from experimental measurements
For each site $r$ and amino-acid $a$, we have an experimental measurement of ACE2-binding or expression.
We will denote these measurements as $x_{r,a}$, using $x_{r,a}^{\rm{bind}}$ or $x_{r,a}^{
\rm{express}}$ when we need to distinguish measurements on ACE2-binding from measurements on RBD expression.

These measurements are all expressed as differences in the log measured value of the mutant relative to the wildtype.
Specifically, for ACE2-binding the measurements are $x_{r,a}^{\rm{bind}} = \log_{10} K_{A,\rm{app}}^{r,a} - \log_{10} K_{A,\rm{app}}^{\rm{wt}}$ where $K_{A,\rm{app}}^{r,a}$ is the apparent assocation constant for binding to ACE2 of the variant with amino-acid $a$ at site $r$, and $K_{A,\rm{app}}^{\rm{wt}}$ is the same assocation constant for wildtype.
For expression, the measurements are $x_{r,a}^{\rm{expr}} = \ln E_{r,a} - \ln E_{\rm{wt}}$ where $E_{r,a}$ is the flow cytometry measured expression for the variant with amino-acid $a$ at site $r$ and $E_{\rm{wt}}$ is the same measurement for wildtype.
So by definition, $x_{r,a} = 0$ when $a$ is equal to the wildtype amino acid at site $r$.

It therefore follows that if we want to compute the "probability" of observing each amino acid at a site under the idea that the probability is proportional to the measured value for that amino acid, we have $p_{r,a} = \frac{\exp\left(\alpha x_{r,a} \right)}{\sum_{a'}\exp\left(\alpha x_{r,a'} \right)}$ where $\alpha > 0$ is a constant.
Larger values of $\alpha$ make the $p_{r,a}$ values more strongly peaked on amino acids with more favorable measurements.

What is the "correct" choice for $\alpha$?
If we interpret the binding probabilities in a thermodynamic perspective, then the value of $\alpha$ should depend on two things: the choice of base when taking the log of the experimental measurement, and the temperature for which we want to compute the probabilities relative to the temperature at which the experimental measurements were made.
In particular, in a thermodynamic perspective, we should have $\alpha = \left(\ln b\right) \times \frac{T_{\rm{experiment}}}{T}$ where $b$ is the base used when taking the log of the experimental measurements, $T$ is the temperature at which we want to compute the probabilities, and $T_{\rm{experiment}}$ is the temperature at which the experiments were performed (where temperatures are in absolute units, e.g. Kelvin).
Assuming that we want to compute probabilities at the same temperature at which we performed the experiments (room temperature), then this gives values of $\alpha_{\rm{bind}} = \ln 10$ (basically correcting for the use of log base-10 when processing experimental data).

There is no natural thermodynamic interpretation for expression, so we arbitrarily choose $\alpha_{\rm{expr}}$ so that the range of the exponents in computing $p_{r,a}^{\rm{expr}}$ is the same as that for computing $p_{r,a}^{\rm{bind}}$, which gives $\alpha_{\rm{expr}} = \alpha_{\rm{bind}} \frac{\max x_{r,a}^{\rm{expr}} - \min x_{r,a}^{\rm{expr}}}{\max x_{r,a}^{\rm{bind}} - \min x_{r,a}^{\rm{bind}}}$.

Note, however, since in the end we are creating visualizations, it would also be reasonable to manually tune the values of $\alpha$ to adjust how much more favorable amino acids have larger letters (on a heat-map representation, this adjustment of $\alpha$ just corresponds to changing the log base or the color scale).
So the the actual calculations below, we add a manual scaling factor $>1$ to make things look a bit more peaked than in the thermodynamic interpretation:

### Actually compute letter heights
We now actually compute the letter heights use the definitions above.
In the final data frame, `prob_bind` / `prob_expr` are the letter heights defined by $p_{r,a}$ above, and `mut_info_bind` / `mut_info_expr` are the letter heights defined by $I_r \times p_{r,a}$.


```python
# get mutation effects without stop codons
nostop_mut_effects = mut_effects.query('(mutant != "*") and (wildtype != "*")')

# exponent scaling factors alpha for binding and expression
manual_scale_alpha = 1.4  # ad hoc additional scaling factor for exponents
alpha_bind = math.log(10) * manual_scale_alpha
alpha_expr = alpha_bind * ((nostop_mut_effects['expr_avg'].max() - nostop_mut_effects['expr_avg'].min()) /
                           (nostop_mut_effects['bind_avg'].max() - nostop_mut_effects['bind_avg'].min()))
print(f"alpha (exponent scale factor):\nbinding = {alpha_bind:.3f}\nexpression = {alpha_expr:.3f}")

def compute_site_info(p):
    p = p[p > 0]
    assert len(p), 'no non-zero values'
    return math.log(20) + sum(p * numpy.log(p))

letter_heights = (
    nostop_mut_effects
    .assign(nonorm_prob_bind=lambda x: numpy.exp(alpha_bind * x['bind_avg']).fillna(0),
            nonorm_prob_expr=lambda x: numpy.exp(alpha_expr * x['expr_avg']).fillna(0),
            sitesum_prob_bind=lambda x: x.groupby('site_RBD')['nonorm_prob_bind'].transform('sum'),
            sitesum_prob_expr=lambda x: x.groupby('site_RBD')['nonorm_prob_expr'].transform('sum'),
            prob_bind=lambda x: x['nonorm_prob_bind'] / x['sitesum_prob_bind'],
            prob_expr=lambda x: x['nonorm_prob_expr'] / x['sitesum_prob_expr'],
            site_info_bind=lambda x: x.groupby('site_RBD')['prob_bind'].transform(compute_site_info),
            site_info_expr=lambda x: x.groupby('site_RBD')['prob_expr'].transform(compute_site_info),
            mut_info_bind=lambda x: x['prob_bind'] * x['site_info_bind'],
            mut_info_expr=lambda x: x['prob_expr'] * x['site_info_expr'],
            )
    [['site_SARS2', 'site_RBD', 'wildtype', 'mutant',
      'prob_bind', 'prob_expr', 'site_info_bind', 'site_info_expr',
      'mut_info_bind', 'mut_info_expr', 'bind_avg', 'expr_avg']]
    )

assert not letter_heights[['mut_info_bind', 'mut_info_expr', 'prob_bind', 'prob_expr']].isnull().any(axis=None)

letter_heights.head().round(3)
```

    alpha (exponent scale factor):
    binding = 3.224
    expression = 3.173





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
      <th>site_SARS2</th>
      <th>site_RBD</th>
      <th>wildtype</th>
      <th>mutant</th>
      <th>prob_bind</th>
      <th>prob_expr</th>
      <th>site_info_bind</th>
      <th>site_info_expr</th>
      <th>mut_info_bind</th>
      <th>mut_info_expr</th>
      <th>bind_avg</th>
      <th>expr_avg</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>A</td>
      <td>0.050</td>
      <td>0.124</td>
      <td>0.011</td>
      <td>0.358</td>
      <td>0.001</td>
      <td>0.044</td>
      <td>-0.03</td>
      <td>-0.11</td>
    </tr>
    <tr>
      <th>1</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>C</td>
      <td>0.041</td>
      <td>0.003</td>
      <td>0.011</td>
      <td>0.358</td>
      <td>0.000</td>
      <td>0.001</td>
      <td>-0.09</td>
      <td>-1.26</td>
    </tr>
    <tr>
      <th>2</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>D</td>
      <td>0.061</td>
      <td>0.044</td>
      <td>0.011</td>
      <td>0.358</td>
      <td>0.001</td>
      <td>0.016</td>
      <td>0.03</td>
      <td>-0.44</td>
    </tr>
    <tr>
      <th>3</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>E</td>
      <td>0.055</td>
      <td>0.066</td>
      <td>0.011</td>
      <td>0.358</td>
      <td>0.001</td>
      <td>0.024</td>
      <td>0.00</td>
      <td>-0.31</td>
    </tr>
    <tr>
      <th>4</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>F</td>
      <td>0.040</td>
      <td>0.019</td>
      <td>0.011</td>
      <td>0.358</td>
      <td>0.000</td>
      <td>0.007</td>
      <td>-0.10</td>
      <td>-0.70</td>
    </tr>
  </tbody>
</table>
</div>



As a quick sanity check, plot the site information content $I_r$ for binding and expression at each site:


```python
p = (ggplot(letter_heights.melt(id_vars='site_SARS2',
                                value_vars=['site_info_bind', 'site_info_expr'],
                                value_name='site information content',
                                var_name='property')
            ) +
     aes('site_SARS2', 'site information content', color='property') +
     geom_step() +
     theme(figure_size=(8, 2.25)) +
     scale_color_manual(values=CBPALETTE[1: ])
     )

_ = p.draw()
```


![png](logoplots_of_muteffects_files/logoplots_of_muteffects_20_0.png)


## Draw logo plots

### Define coloring and highlighting scheme
First combine letter heights with key site annotations.
We want to plot as contact sites ones that contact ACE2 in **either** SARS-CoV-2 or SARS-CoV-1, so make a column that indicates this:


```python
logo_df = (
    letter_heights
    .merge(site_annotations.assign(ACE2_contact=lambda x: x['SARS2_ACE2_contact'] | x['SARS1_ACE2_contact'])
                           [['site_SARS2', 'amino_acid_SARS1', 'RSA_bound', 'RSA_unbound', 'ACE2_contact']],
           on='site_SARS2',
           validate='many_to_one')
    )

logo_df.head()
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
      <th>site_SARS2</th>
      <th>site_RBD</th>
      <th>wildtype</th>
      <th>mutant</th>
      <th>prob_bind</th>
      <th>prob_expr</th>
      <th>site_info_bind</th>
      <th>site_info_expr</th>
      <th>mut_info_bind</th>
      <th>mut_info_expr</th>
      <th>bind_avg</th>
      <th>expr_avg</th>
      <th>amino_acid_SARS1</th>
      <th>RSA_bound</th>
      <th>RSA_unbound</th>
      <th>ACE2_contact</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>A</td>
      <td>0.050009</td>
      <td>0.124117</td>
      <td>0.010774</td>
      <td>0.358381</td>
      <td>0.000539</td>
      <td>0.044481</td>
      <td>-0.03</td>
      <td>-0.11</td>
      <td>N</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
    </tr>
    <tr>
      <th>1</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>C</td>
      <td>0.041214</td>
      <td>0.003228</td>
      <td>0.010774</td>
      <td>0.358381</td>
      <td>0.000444</td>
      <td>0.001157</td>
      <td>-0.09</td>
      <td>-1.26</td>
      <td>N</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>D</td>
      <td>0.060681</td>
      <td>0.043553</td>
      <td>0.010774</td>
      <td>0.358381</td>
      <td>0.000654</td>
      <td>0.015609</td>
      <td>0.03</td>
      <td>-0.44</td>
      <td>N</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>E</td>
      <td>0.055087</td>
      <td>0.065794</td>
      <td>0.010774</td>
      <td>0.358381</td>
      <td>0.000594</td>
      <td>0.023579</td>
      <td>0.00</td>
      <td>-0.31</td>
      <td>N</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>F</td>
      <td>0.039907</td>
      <td>0.019085</td>
      <td>0.010774</td>
      <td>0.358381</td>
      <td>0.000430</td>
      <td>0.006840</td>
      <td>-0.10</td>
      <td>-0.70</td>
      <td>N</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
</div>



Now assign colors. 

We color amino acid letters:
 - dark blue if wildtype in SARS-CoV-2
 - green if wildtype in SARS-CoV-1 (but not SARS-CoV-2)
 - gray otherwise.
 
We highlight sites as follows (always using an alpha transparency of 0.25):
 - yellow if an ACE2 contact site (in SARS-CoV-1 or SARS-CoV-2)


```python
def assign_color(row):
    if row['wildtype'] == row['mutant']:
        return CBPALETTE[5]
    elif row['mutant'] == row['amino_acid_SARS1']:
        return CBPALETTE[3]
    else:
        return CBPALETTE[0]
    
def assign_highlight_color(row):
    if row['ACE2_contact'] is True:
        return CBPALETTE[4]
    else:
        return None
    
logo_df = (
    logo_df
    .assign(color=lambda x: x.apply(assign_color, axis=1),
            highlight_color=lambda x: x.apply(assign_highlight_color, axis=1),
            highlight_alpha=0.25,
            )
    )
```

### Draw logos for binding and expression separately
First we draw logo plots that show the amino acids preferred for binding and expression separately in either the probability representation.
For now, we don't show the information representation: if you want that, use `info_bind` / `info_expr` rather than `prob_bind` / `prob_expr` below:


```python
for height_col, ylabel in [
        ('prob_bind', 'binding to ACE2'),
        ('prob_expr', 'RBD expression'),
        ]:

    filename = os.path.join(config['logoplots_dir'], f"{height_col}_logo.pdf")
    print(f"\nPlotting {height_col} and saving to {filename}...")
    
    fig = logomaker_utils.line_wrapped_logo(
        logo_df,
        site_col='site_SARS2',
        letter_col='mutant',
        height_col=height_col,
        color_col='color',
        highlight_color_col='highlight_color',
        highlight_alpha_col='highlight_alpha',
        missing_letter='zero_height',
        sites_per_line=51,  # number of sites per line
        scaleheight=1.2,  # stretch plot vertically by this much
        style_xticks_kwargs={'spacing': 2,  # number every other residue  
                             },
        logo_kwargs={'font_name': 'Uni Sans',
                     'width': 0.9,  # width of letters
                     'vpad': 0.05,  # vertical padding between letters
                     'stack_order': 'big_on_top',  # big letters on top
                     },
        xlabel='site',
        ylabel=ylabel,
        )
    
    display(fig)
    fig.savefig(filename)
    plt.close(fig)
```

    
    Plotting prob_bind and saving to results/logoplots/prob_bind_logo.pdf...



![png](logoplots_of_muteffects_files/logoplots_of_muteffects_27_1.png)


    
    Plotting prob_expr and saving to results/logoplots/prob_expr_logo.pdf...



![png](logoplots_of_muteffects_files/logoplots_of_muteffects_27_3.png)


### Draw combined logo with binding and expression
We draw logos that show **both** binding and expression data on the same plot, with binding as capital letters above the line and binding as lower-case letters below the center line.

First, make a data frame suitable for this by concatenating the logo data frame, mutating the letters for binding to be lower case and re-naming columns to no longer specify binding versus expression.
We also make a column that colors letters differently if they are for expression or binding:


```python
combo_logo_df = pd.concat([
        # data for binding
        logo_df.drop(columns=[c for c in logo_df.columns if '_expr' in c])
               .rename(columns={c: c.replace('_bind', '') for c in logo_df.columns})
               .assign(property_type='binding'),
        # data for expression, make lowercase letters and negative heights
        logo_df.drop(columns=[c for c in logo_df.columns if '_bind' in c])
               .rename(columns={c: c.replace('_expr', '') for c in logo_df.columns})
               .assign(mutant=lambda x: x['mutant'].str.lower(),
                       prob=lambda x: -x['prob'],
                       mut_info=lambda x: -x['mut_info'],
                       property_type='expression'),
        ],
        ignore_index=True, sort=False)

combo_logo_df
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
      <th>site_SARS2</th>
      <th>site_RBD</th>
      <th>wildtype</th>
      <th>mutant</th>
      <th>prob</th>
      <th>site_info</th>
      <th>mut_info</th>
      <th>bind_avg</th>
      <th>expr_avg</th>
      <th>amino_acid_SARS1</th>
      <th>RSA_bound</th>
      <th>RSA_unbound</th>
      <th>ACE2_contact</th>
      <th>color</th>
      <th>highlight_color</th>
      <th>highlight_alpha</th>
      <th>property_type</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>A</td>
      <td>0.050009</td>
      <td>0.010774</td>
      <td>0.000539</td>
      <td>-0.03</td>
      <td>-0.11</td>
      <td>N</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>#999999</td>
      <td>None</td>
      <td>0.25</td>
      <td>binding</td>
    </tr>
    <tr>
      <th>1</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>C</td>
      <td>0.041214</td>
      <td>0.010774</td>
      <td>0.000444</td>
      <td>-0.09</td>
      <td>-1.26</td>
      <td>N</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>#999999</td>
      <td>None</td>
      <td>0.25</td>
      <td>binding</td>
    </tr>
    <tr>
      <th>2</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>D</td>
      <td>0.060681</td>
      <td>0.010774</td>
      <td>0.000654</td>
      <td>0.03</td>
      <td>-0.44</td>
      <td>N</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>#999999</td>
      <td>None</td>
      <td>0.25</td>
      <td>binding</td>
    </tr>
    <tr>
      <th>3</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>E</td>
      <td>0.055087</td>
      <td>0.010774</td>
      <td>0.000594</td>
      <td>0.00</td>
      <td>-0.31</td>
      <td>N</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>#999999</td>
      <td>None</td>
      <td>0.25</td>
      <td>binding</td>
    </tr>
    <tr>
      <th>4</th>
      <td>331</td>
      <td>1</td>
      <td>N</td>
      <td>F</td>
      <td>0.039907</td>
      <td>0.010774</td>
      <td>0.000430</td>
      <td>-0.10</td>
      <td>-0.70</td>
      <td>N</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>#999999</td>
      <td>None</td>
      <td>0.25</td>
      <td>binding</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>8035</th>
      <td>531</td>
      <td>201</td>
      <td>T</td>
      <td>s</td>
      <td>-0.054449</td>
      <td>0.008464</td>
      <td>-0.000461</td>
      <td>0.01</td>
      <td>0.02</td>
      <td>T</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>#999999</td>
      <td>None</td>
      <td>0.25</td>
      <td>expression</td>
    </tr>
    <tr>
      <th>8036</th>
      <td>531</td>
      <td>201</td>
      <td>T</td>
      <td>t</td>
      <td>-0.051101</td>
      <td>0.008464</td>
      <td>-0.000433</td>
      <td>0.00</td>
      <td>0.00</td>
      <td>T</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>#0072B2</td>
      <td>None</td>
      <td>0.25</td>
      <td>expression</td>
    </tr>
    <tr>
      <th>8037</th>
      <td>531</td>
      <td>201</td>
      <td>T</td>
      <td>v</td>
      <td>-0.042241</td>
      <td>0.008464</td>
      <td>-0.000358</td>
      <td>0.01</td>
      <td>-0.06</td>
      <td>T</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>#999999</td>
      <td>None</td>
      <td>0.25</td>
      <td>expression</td>
    </tr>
    <tr>
      <th>8038</th>
      <td>531</td>
      <td>201</td>
      <td>T</td>
      <td>w</td>
      <td>-0.039643</td>
      <td>0.008464</td>
      <td>-0.000336</td>
      <td>-0.02</td>
      <td>-0.08</td>
      <td>T</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>#999999</td>
      <td>None</td>
      <td>0.25</td>
      <td>expression</td>
    </tr>
    <tr>
      <th>8039</th>
      <td>531</td>
      <td>201</td>
      <td>T</td>
      <td>y</td>
      <td>-0.043603</td>
      <td>0.008464</td>
      <td>-0.000369</td>
      <td>-0.01</td>
      <td>-0.05</td>
      <td>T</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>False</td>
      <td>#999999</td>
      <td>None</td>
      <td>0.25</td>
      <td>expression</td>
    </tr>
  </tbody>
</table>
<p>8040 rows × 17 columns</p>
</div>



Now draw these logo plots with binding (upper-case letters) at top and expression (lower-case letters) at bottom.
We plot this only for the probabilities; if you want the information content then change `height_col='prob'` to `height_col='info'`:


```python
filename = os.path.join(config['logoplots_dir'], f"prob_logo.pdf")
print(f"\nPlotting and saving to {filename}...")
   
fig = logomaker_utils.line_wrapped_logo(
        combo_logo_df,
        site_col='site_SARS2',
        letter_col='mutant',
        height_col='prob',
        color_col='color',
        highlight_color_col='highlight_color',
        highlight_alpha_col='highlight_alpha',
        missing_letter='zero_height',
        sites_per_line=51,
        scaleheight=1.5,
        fade_letters_by_height=(0.5, 1),  # smaller letters more transparent
        style_xticks_kwargs={'spacing': 2,  # number every other residue  
                             },
        logo_kwargs={'font_name': 'Uni Sans',
                     'width': 0.9,  # width of letters
                     'vpad': 0.05,  # vertical padding between letters
                     'flip_below': False,  # flip letters below line
                     'stack_order': 'small_on_top',  # small letters closer to baseline
                     },
        xlabel='site',
        ylabel='RBD expression $\Longleftrightarrow$ binding to ACE2',
        all_letters=combo_logo_df['mutant'].unique()
        )
    
display(fig)
fig.savefig(filename)
plt.close(fig)
```

    
    Plotting and saving to results/logoplots/prob_logo.pdf...



![png](logoplots_of_muteffects_files/logoplots_of_muteffects_31_1.png)


## Zoomed logos
Make some zoomed logos:

### Look at cysteine pairs
Plot cysteine pairs for expression only.
We want a space between them, so make a dummy site label called `isite` that skips between each pair.


```python
cys_pairs = [(336, 361),
             (379, 432),
             (391, 525),
             (480, 488)]

site_to_isite = {}
for i, (cys1, cys2) in enumerate(cys_pairs):
    site_to_isite[cys1] = 3 * i
    site_to_isite[cys2] = 3 * i + 1
    
cys_df = (
    logo_df
    .query(f"site_SARS2 in {list(site_to_isite)}")
    .assign(isite=lambda x: x['site_SARS2'].map(site_to_isite))
    )

fig = logomaker_utils.line_wrapped_logo(
        cys_df,
        site_col='isite',
        letter_col='mutant',
        height_col='prob_bind',
        color_col='color',
        sitelabel_col='site_SARS2',
        missing_letter='zero_height',
        sites_per_line=10,
        scaleheight=1.3,
        scalewidth=1.6,
        logo_kwargs={'font_name': 'Uni Sans',
                     'width': 0.9,  # width of letters
                     'vpad': 0.05,  # vertical padding between letters
                     'flip_below': True,  # flip letters below line
                     'stack_order': 'big_on_top',  # small letters closer to baseline
                     },
        xlabel='site',
        )
    
display(fig)
plt.close(fig)
```


![png](logoplots_of_muteffects_files/logoplots_of_muteffects_34_0.png)


### Look at contact residues
Plot the ACE2 contact residues in the combined expression / binding logo.
Note how we subset on ACE2 contact sites and then make a dummy variable `isite`.
We don't put spaces between non-consecutive sites here.


```python
contact_sites = sorted(logo_df.query('ACE2_contact == True')['site_SARS2'].unique())
site_to_isite = {site: isite for isite, site in enumerate(contact_sites)}

contact_df = (
    combo_logo_df
    .query(f"site_SARS2 in {list(site_to_isite)}")
    .assign(isite=lambda x: x['site_SARS2'].map(site_to_isite))
    )

fig = logomaker_utils.line_wrapped_logo(
        contact_df,
        site_col='isite',
        letter_col='mutant',
        height_col='prob',
        color_col='color',
        sitelabel_col='site_SARS2',
        missing_letter='zero_height',
        sites_per_line=51,
        scaleheight=1.9,
        scalewidth=1.2,
        fade_letters_by_height=(0.5, 1),  # smaller letters more transparent
        logo_kwargs={'font_name': 'Uni Sans',
                     'width': 0.9,  # width of letters
                     'vpad': 0.05,  # vertical padding between letters
                     'flip_below': True,  # flip letters below line
                     'stack_order': 'small_on_top',  # small letters closer to baseline
                     },
        xlabel='site',
        ylabel='expression $\leftrightarrow$ binding',
        all_letters=combo_logo_df['mutant'].unique()
        )
    
display(fig)
plt.close(fig)
```


![png](logoplots_of_muteffects_files/logoplots_of_muteffects_36_0.png)


## Create input file for `dms-view`


```python
dms_view_df = (
    letter_heights
    .rename(columns={'site_RBD': 'site',
                     'site_SARS2': 'label_site',
                     'mutant': 'mutation',
                     'prob_bind': 'ACE2-binding',
                     'prob_expr': 'expression',
                     })
    .assign(protein_chain='E',
            protein_site=lambda x: x['label_site'],
            )
    .melt(id_vars=['site', 'label_site', 'wildtype', 'mutation',
                   'protein_chain', 'protein_site'],
          value_vars=['ACE2-binding', 'expression'],
          var_name='condition',
          value_name='mut_preference',
          )
    .merge(letter_heights
           .rename(columns={'bind_avg': 'ACE2-binding',
                            'expr_avg': 'expression',
                            'site_SARS2': 'label_site',
                            'mutant': 'mutation',
                            })
           .melt(id_vars=['label_site', 'mutation'],
                 value_vars=['ACE2-binding', 'expression'],
                 var_name='condition',
                 value_name='mut_delta_effect',
                 ),
           validate='one_to_one'
           )
    .assign(
        site_entropy=lambda x: x.groupby(['site', 'condition'])
                                ['mut_preference']
                                .transform(lambda p: -sum(p[p > 0] * 
                                                          numpy.log(p[p > 0]))),
        site_n_effective=lambda x: numpy.exp(x['site_entropy']),
        site_mean_effect=lambda x: x.groupby(['site', 'condition'])
                                    ['mut_delta_effect']
                                    .transform('mean'),
        site_min_effect=lambda x: x.groupby(['site', 'condition'])
                                    ['mut_delta_effect']
                                    .transform('min'),
        site_max_effect=lambda x: x.groupby(['site', 'condition'])
                                   ['mut_delta_effect']
                                   .transform('max'),
        )
    )

print(f"Writing to {config['dms_view_file']}")
os.makedirs(config['dms_view_dir'], exist_ok=True)
dms_view_df.to_csv(config['dms_view_file'], index=False)

dms_view_df
```

    Writing to results/dms_view/dms-view_table.csv





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
      <th>site</th>
      <th>label_site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>protein_chain</th>
      <th>protein_site</th>
      <th>condition</th>
      <th>mut_preference</th>
      <th>mut_delta_effect</th>
      <th>site_entropy</th>
      <th>site_n_effective</th>
      <th>site_mean_effect</th>
      <th>site_min_effect</th>
      <th>site_max_effect</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>E</td>
      <td>331</td>
      <td>ACE2-binding</td>
      <td>0.050009</td>
      <td>-0.03</td>
      <td>2.984958</td>
      <td>19.785669</td>
      <td>-0.0335</td>
      <td>-0.16</td>
      <td>0.06</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>C</td>
      <td>E</td>
      <td>331</td>
      <td>ACE2-binding</td>
      <td>0.041214</td>
      <td>-0.09</td>
      <td>2.984958</td>
      <td>19.785669</td>
      <td>-0.0335</td>
      <td>-0.16</td>
      <td>0.06</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>E</td>
      <td>331</td>
      <td>ACE2-binding</td>
      <td>0.060681</td>
      <td>0.03</td>
      <td>2.984958</td>
      <td>19.785669</td>
      <td>-0.0335</td>
      <td>-0.16</td>
      <td>0.06</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>E</td>
      <td>331</td>
      <td>ACE2-binding</td>
      <td>0.055087</td>
      <td>0.00</td>
      <td>2.984958</td>
      <td>19.785669</td>
      <td>-0.0335</td>
      <td>-0.16</td>
      <td>0.06</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>F</td>
      <td>E</td>
      <td>331</td>
      <td>ACE2-binding</td>
      <td>0.039907</td>
      <td>-0.10</td>
      <td>2.984958</td>
      <td>19.785669</td>
      <td>-0.0335</td>
      <td>-0.16</td>
      <td>0.06</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>8035</th>
      <td>201</td>
      <td>531</td>
      <td>T</td>
      <td>S</td>
      <td>E</td>
      <td>531</td>
      <td>expression</td>
      <td>0.054449</td>
      <td>0.02</td>
      <td>2.987268</td>
      <td>19.831436</td>
      <td>-0.0095</td>
      <td>-0.08</td>
      <td>0.08</td>
    </tr>
    <tr>
      <th>8036</th>
      <td>201</td>
      <td>531</td>
      <td>T</td>
      <td>T</td>
      <td>E</td>
      <td>531</td>
      <td>expression</td>
      <td>0.051101</td>
      <td>0.00</td>
      <td>2.987268</td>
      <td>19.831436</td>
      <td>-0.0095</td>
      <td>-0.08</td>
      <td>0.08</td>
    </tr>
    <tr>
      <th>8037</th>
      <td>201</td>
      <td>531</td>
      <td>T</td>
      <td>V</td>
      <td>E</td>
      <td>531</td>
      <td>expression</td>
      <td>0.042241</td>
      <td>-0.06</td>
      <td>2.987268</td>
      <td>19.831436</td>
      <td>-0.0095</td>
      <td>-0.08</td>
      <td>0.08</td>
    </tr>
    <tr>
      <th>8038</th>
      <td>201</td>
      <td>531</td>
      <td>T</td>
      <td>W</td>
      <td>E</td>
      <td>531</td>
      <td>expression</td>
      <td>0.039643</td>
      <td>-0.08</td>
      <td>2.987268</td>
      <td>19.831436</td>
      <td>-0.0095</td>
      <td>-0.08</td>
      <td>0.08</td>
    </tr>
    <tr>
      <th>8039</th>
      <td>201</td>
      <td>531</td>
      <td>T</td>
      <td>Y</td>
      <td>E</td>
      <td>531</td>
      <td>expression</td>
      <td>0.043603</td>
      <td>-0.05</td>
      <td>2.987268</td>
      <td>19.831436</td>
      <td>-0.0095</td>
      <td>-0.08</td>
      <td>0.08</td>
    </tr>
  </tbody>
</table>
<p>8040 rows × 14 columns</p>
</div>


