# Fit global epistasis models to DMS data

The [`dms_variants.globalepistasis`](https://jbloomlab.github.io/dms_variants/dms_variants.globalepistasis) module is based off the "global epistasis" concepts of [Otwinoski et
al](https://doi.org/10.1073/pnas.1804015115) and [Sailer and
Harms](https://doi.org/10.1534/genetics.116.195214) -- that there exists an underlying latent phenotype that mutations affect additively, and
then an observed (measured) phenotype that is a non-linear
function of the latent phenotype.


## Setup for analysis

Import Python modules / packages:


```python
import collections
import os
import itertools
import random
import tempfile
import time
import warnings

import pandas as pd

from plotnine import *

import dms_variants.binarymap
import dms_variants.codonvarianttable
import dms_variants.globalepistasis
import dms_variants.plotnine_themes
import dms_variants.simulate
from dms_variants.constants import CBPALETTE, CODONS_NOSTOP

import yaml
```

Set the gray grid [plotnine theme](https://plotnine.readthedocs.io/en/stable/generated/plotnine.themes.theme.html?highlight=themes) defined in [dms_variants.plotnine_themes](https://jbloomlab.github.io/dms_variants/dms_variants.plotnine_themes.html):


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using dms_variants version 0.6.0


Set pandas display options to show large chunks of Data Frames in this
example:


```python
pd.set_option('display.max_columns', 30)
pd.set_option('display.width', 500)
```

Hide warnings that clutter output:


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
os.makedirs(config['global_epistasis_binding_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Read in Tite-seq affinity scores by barcode
Read in log10Ka measurements. Rename log10Ka column to be func_score, and calculate variance from the ddG_log10Ka column. Remove rows with NaN for func_score. For empty aa_substitutions (wildtype), it is being replaced with NA. Make back to an empty string.


```python
df = pd.read_csv(config['Titeseq_Kds_file'])
df.rename(columns={'log10Ka':'func_score'},inplace=True)
df['func_score_var'] = df['log10SE']**2
func_scores = df[pd.notnull(df['func_score'])]
func_scores.fillna('',inplace=True)
func_scores.head()
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
      <th>target</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>avgcount</th>
      <th>func_score</th>
      <th>delta_log10Ka</th>
      <th>log10SE</th>
      <th>Kd</th>
      <th>Kd_SE</th>
      <th>response</th>
      <th>baseline</th>
      <th>nMSR</th>
      <th>variant_class</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>func_score_var</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAAACTTCAAT</td>
      <td>5</td>
      <td>74.50</td>
      <td>8.72</td>
      <td>-2.05</td>
      <td>0.12</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.64</td>
      <td>1.13</td>
      <td>0.01</td>
      <td>&gt;1 nonsynonymous</td>
      <td>A22C R127G E141D L188V</td>
      <td>4</td>
      <td>0.0144</td>
    </tr>
    <tr>
      <th>3</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAACAAGCAGA</td>
      <td>6</td>
      <td>146.32</td>
      <td>10.35</td>
      <td>-0.42</td>
      <td>0.05</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.87</td>
      <td>1.01</td>
      <td>0.00</td>
      <td>1 nonsynonymous</td>
      <td>N13F</td>
      <td>1</td>
      <td>0.0025</td>
    </tr>
    <tr>
      <th>5</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAACAGGTTGC</td>
      <td>4</td>
      <td>47.48</td>
      <td>6.00</td>
      <td>-4.76</td>
      <td>2.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.50</td>
      <td>1.15</td>
      <td>0.01</td>
      <td>&gt;1 nonsynonymous</td>
      <td>V71K P149L N157T</td>
      <td>3</td>
      <td>4.0000</td>
    </tr>
    <tr>
      <th>6</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAACATTAAAT</td>
      <td>6</td>
      <td>47.44</td>
      <td>10.15</td>
      <td>-0.61</td>
      <td>0.09</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.64</td>
      <td>1.09</td>
      <td>0.00</td>
      <td>&gt;1 nonsynonymous</td>
      <td>A18V T148S H189Y</td>
      <td>3</td>
      <td>0.0081</td>
    </tr>
    <tr>
      <th>8</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAACCTTACAA</td>
      <td>2</td>
      <td>18.78</td>
      <td>9.61</td>
      <td>-1.15</td>
      <td>0.19</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.00</td>
      <td>1.11</td>
      <td>0.01</td>
      <td>&gt;1 nonsynonymous</td>
      <td>T63D A89N</td>
      <td>2</td>
      <td>0.0361</td>
    </tr>
  </tbody>
</table>
</div>



## Fit global epistasis models
We now fit global epistasis models to the functional scores.
For background on these models, see [Otwinoski et al (2018)](https://www.pnas.org/content/115/32/E7550).
The models we fits are implemented in [dms_variants.globalepistasis](https://jbloomlab.github.io/dms_variants/dms_variants.globalepistasis.html), which has extensive documentation.

The primary model of interest is [MonotonicSplineEpistasis](https://jbloomlab.github.io/dms_variants/dms_variants.globalepistasis.html#dms_variants.globalepistasis.MonotonicSplineEpistasis), which assumes that the observed phenotype is a monotonic non-linear function of an underlying additive latent phenotype.
As a control, we also fit a [NoEpistasis](https://jbloomlab.github.io/dms_variants/dms_variants.globalepistasis.html#dms_variants.globalepistasis.NoEpistasis) model, which assumes that mutations simply contribute additively to the observed phenotype:

For the fitting, we first convert the set of functional scores into a binary representation using a [BinaryMap]((https://jbloomlab.github.io/dms_variants/dms_variants.binarymap.html#dms_variants.binarymap.BinaryMap)).
Then we create the model, fit it, and store it.


```python
# NBVAL_IGNORE_OUTPUT

models = {}  # store models, keyed by `(epistasistype, likelihoodtype, lib)`

for (lib), scores in func_scores.groupby(['library']):
   
    bmap = dms_variants.binarymap.BinaryMap(scores)
    
    for epistasistype, likelihoodtype, Model in [
            ('global epistasis', 'Gaussian', dms_variants.globalepistasis.MonotonicSplineEpistasisGaussianLikelihood),
            ('no epistasis', 'Gaussian', dms_variants.globalepistasis.NoEpistasisGaussianLikelihood),
            ('global epistasis', 'Cauchy', dms_variants.globalepistasis.MonotonicSplineEpistasisCauchyLikelihood),
            ('no epistasis', 'Cauchy', dms_variants.globalepistasis.NoEpistasisCauchyLikelihood),
            ]:
        print(f"Fitting {epistasistype} with {likelihoodtype} likelihood model to {lib}...", end=' ')
    
        start = time.time()
        model = Model(bmap)
        model.fit()  # do NOT change ftol in normal use, this is just for test
        print(f"fitting took {time.time() - start:.1f} sec.")
        models[(epistasistype, likelihoodtype, lib)] = model
```

    Fitting global epistasis with Gaussian likelihood model to lib1... fitting took 84.7 sec.
    Fitting no epistasis with Gaussian likelihood model to lib1... fitting took 3.2 sec.
    Fitting global epistasis with Cauchy likelihood model to lib1... fitting took 213.4 sec.
    Fitting no epistasis with Cauchy likelihood model to lib1... fitting took 16.7 sec.
    Fitting global epistasis with Gaussian likelihood model to lib2... fitting took 81.9 sec.
    Fitting no epistasis with Gaussian likelihood model to lib2... fitting took 3.3 sec.
    Fitting global epistasis with Cauchy likelihood model to lib2... fitting took 238.5 sec.
    Fitting no epistasis with Cauchy likelihood model to lib2... fitting took 21.3 sec.


Now we want to see which model fits the data better.
To do this, we get the log likelihood of each model along with the number of model parameters and use it to calculate the [AIC](https://en.wikipedia.org/wiki/Akaike_information_criterion).
Models with lower AIC are better, and below we see that the [MonotonicSplineEpistasis](https://jbloomlab.github.io/dms_variants/dms_variants.globalepistasis.html#dms_variants.globalepistasis.MonotonicSplineEpistasis) global epistasis model always fits the data much better. (Note, I previous tried alterning number of meshpoints, and the default 4 was the lowest AIC):


```python
# NBVAL_IGNORE_OUTPUT

logliks_df = (
    pd.DataFrame.from_records(
            [(epistasistype, likelihoodtype, lib, model.nparams, model.loglik) for
             (epistasistype, likelihoodtype, lib), model in models.items()],
            columns=['model', 'likelihood type', 'library',
                     'n_parameters', 'log_likelihood']
            )
    .assign(AIC=lambda x: 2 * x['n_parameters'] - 2 * x['log_likelihood'])
    .set_index(['library'])
    )

logliks_df.round(1)
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
      <th>model</th>
      <th>likelihood type</th>
      <th>n_parameters</th>
      <th>log_likelihood</th>
      <th>AIC</th>
    </tr>
    <tr>
      <th>library</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>lib1</th>
      <td>global epistasis</td>
      <td>Gaussian</td>
      <td>3801</td>
      <td>-53065.8</td>
      <td>113733.5</td>
    </tr>
    <tr>
      <th>lib1</th>
      <td>no epistasis</td>
      <td>Gaussian</td>
      <td>3795</td>
      <td>-80441.6</td>
      <td>168473.3</td>
    </tr>
    <tr>
      <th>lib1</th>
      <td>global epistasis</td>
      <td>Cauchy</td>
      <td>3801</td>
      <td>-39036.0</td>
      <td>85674.0</td>
    </tr>
    <tr>
      <th>lib1</th>
      <td>no epistasis</td>
      <td>Cauchy</td>
      <td>3795</td>
      <td>-60021.4</td>
      <td>127632.9</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>global epistasis</td>
      <td>Gaussian</td>
      <td>3798</td>
      <td>-49768.0</td>
      <td>107132.0</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>no epistasis</td>
      <td>Gaussian</td>
      <td>3792</td>
      <td>-73422.4</td>
      <td>154428.8</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>global epistasis</td>
      <td>Cauchy</td>
      <td>3798</td>
      <td>-32464.9</td>
      <td>72525.8</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>no epistasis</td>
      <td>Cauchy</td>
      <td>3792</td>
      <td>-54067.1</td>
      <td>115718.3</td>
    </tr>
  </tbody>
</table>
</div>




```python
# check to confirm the global epistasis model is better in all cases
assert (logliks_df
        .reset_index()
        .pivot_table(index=['library'],
                     values='AIC',
                     columns='model'
                     )
        .assign(global_better=lambda x: x['global epistasis'] < x['no epistasis'])
        ['global_better']
        .all()
        )
```

Because the [MonotonicSplineEpistasis](https://jbloomlab.github.io/dms_variants/dms_variants.globalepistasis.html#dms_variants.globalepistasis.MonotonicSplineEpistasis) global epistasis models fit so much better than the additive [NoEpistasis](https://jbloomlab.github.io/dms_variants/dms_variants.globalepistasis.html#dms_variants.globalepistasis.NoEpistasis) models, below we are just going to analyze the results from the global epistasis model.

First, we will examine how the model looks on all the actual variants used to fit the model.
We use [AbstractEpistasis.phenotypes_df](https://jbloomlab.github.io/dms_variants/dms_variants.globalepistasis.html#dms_variants.globalepistasis.AbstractEpistasis.phenotypes_df) to get a data frame of all the variants used to fit each global epistasis model along with their functional scores and the latent and observed phenotypes predicted by each model.


```python
# NBVAL_IGNORE_OUTPUT

variants_df = pd.concat(
        [model.phenotypes_df
         .assign(library=lib,
                 likelihoodtype=likelihoodtype,
                 )
         for (epistasistype, likelihoodtype, lib), model in models.items()
         if (epistasistype == 'global epistasis')],
        ignore_index=True, sort=False)

#predictionsfile = os.path.join(config['global_epistasis_binding_dir'], 'globalepistasis_binding_predictions.csv')
#variants_df.to_csv(predictionsfile, index=False)
#print(f"Writing predictions to {predictionsfile}")

variants_df.head().round(2)
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
      <th>aa_substitutions</th>
      <th>func_score</th>
      <th>func_score_var</th>
      <th>latent_phenotype</th>
      <th>observed_phenotype</th>
      <th>library</th>
      <th>likelihoodtype</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>A22C R127G E141D L188V</td>
      <td>8.72</td>
      <td>0.01</td>
      <td>-1.54</td>
      <td>9.06</td>
      <td>lib1</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>1</th>
      <td>N13F</td>
      <td>10.35</td>
      <td>0.00</td>
      <td>-0.62</td>
      <td>10.46</td>
      <td>lib1</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>2</th>
      <td>V71K P149L N157T</td>
      <td>6.00</td>
      <td>4.00</td>
      <td>-3.99</td>
      <td>6.64</td>
      <td>lib1</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>3</th>
      <td>A18V T148S H189Y</td>
      <td>10.15</td>
      <td>0.01</td>
      <td>-0.89</td>
      <td>10.14</td>
      <td>lib1</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>4</th>
      <td>T63D A89N</td>
      <td>9.61</td>
      <td>0.04</td>
      <td>-1.33</td>
      <td>9.43</td>
      <td>lib1</td>
      <td>Gaussian</td>
    </tr>
  </tbody>
</table>
</div>



Below we plot the relationships among the latent phenotype from the model, the observed phenotype from the model, and the measured functional score for all variants used to fit the model:


```python
for x, y in itertools.combinations(['latent_phenotype',
                                    'observed_phenotype',
                                    'func_score'],
                                   2):
    p = (
        ggplot(variants_df, aes(x, y)) +
        geom_point(alpha=0.05, size=0.5) +
        facet_grid('library ~ likelihoodtype', scales='free_y') +
        theme(figure_size=(2 * variants_df['likelihoodtype'].nunique(),
                           2 * variants_df['library'].nunique()),
              )
        )
#         plotfile = os.path.join(config['figs_dir'], f'{y}-v-{x}_by_{epistasistype}.pdf')
#         print(f"Saving to {plotfile}")
#         p.save(plotfile)
    _ = p.draw()
```


![png](global_epistasis_binding_files/global_epistasis_binding_24_0.png)



![png](global_epistasis_binding_files/global_epistasis_binding_24_1.png)



![png](global_epistasis_binding_files/global_epistasis_binding_24_2.png)


The two Cauchy fits looks good -- clearly picking up on the complete censoring of *K*<sub>A,app</sub> values >10<sup>6</sup> to the 10<sup>6</sup> boundary condition. The Gaussian likelihood model fails to find this boundary in the lib2 data, and in the lib1 data it seems to fit a boundary closer to 7 than to 6. I am curious to see how these curves change when I remove the SE/variance estimates, which get wonky in these values at this boundary (and are of uncertain value even for points within the dynamic range).

The Cauchy likelihood fits also have fewer points that are observed at this boundary condition in the experimental functional scores but given latent phenotypes near neutral (aka, likely "false positives"). Overall, these models that use the Cauchy liklelihood look better, which is consistent with what I've seen from previous Tite-seq results. We will confirm this later on when we evaluate the correlation in coefficients between replicates from these different models.

To get a better view into the "shape" of global epistasis, let's re-make the plots above but only showing single mutant barcodes.


```python
mask = (variants_df['aa_substitutions'].str.len() < 6) & (variants_df['aa_substitutions'].str.len() > 0)
single_variants_df = variants_df.loc[mask]
single_variants_df.head()
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
      <th>aa_substitutions</th>
      <th>func_score</th>
      <th>func_score_var</th>
      <th>latent_phenotype</th>
      <th>observed_phenotype</th>
      <th>library</th>
      <th>likelihoodtype</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>N13F</td>
      <td>10.35</td>
      <td>0.0025</td>
      <td>-0.620878</td>
      <td>10.456437</td>
      <td>lib1</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>10</th>
      <td>P7S</td>
      <td>10.67</td>
      <td>0.0081</td>
      <td>-0.082381</td>
      <td>10.731078</td>
      <td>lib1</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>21</th>
      <td>P149Q</td>
      <td>10.72</td>
      <td>0.0025</td>
      <td>-0.234197</td>
      <td>10.714161</td>
      <td>lib1</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>24</th>
      <td>D90Y</td>
      <td>10.57</td>
      <td>0.0049</td>
      <td>-0.622221</td>
      <td>10.455074</td>
      <td>lib1</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>26</th>
      <td>L5C</td>
      <td>10.69</td>
      <td>0.0036</td>
      <td>-0.257591</td>
      <td>10.706989</td>
      <td>lib1</td>
      <td>Gaussian</td>
    </tr>
  </tbody>
</table>
</div>




```python
for x, y in itertools.combinations(['latent_phenotype',
                                    'observed_phenotype',
                                    'func_score'],
                                   2):
    p = (
        ggplot(single_variants_df, aes(x, y)) +
        geom_point(alpha=0.05, size=0.5) +
        facet_grid('library ~ likelihoodtype', scales='free_y') +
        theme(figure_size=(2 * variants_df['likelihoodtype'].nunique(),
                           2 * variants_df['library'].nunique()),
              )
        )
#         plotfile = os.path.join(config['figs_dir'], f'{y}-v-{x}_by_{epistasistype}.pdf')
#         print(f"Saving to {plotfile}")
#         p.save(plotfile)
    _ = p.draw()
```


![png](global_epistasis_binding_files/global_epistasis_binding_27_0.png)



![png](global_epistasis_binding_files/global_epistasis_binding_27_1.png)



![png](global_epistasis_binding_files/global_epistasis_binding_27_2.png)



```python
#add predicted and latent phenotypes to table by barcode with additional columns used for interpretation, save output file
dt = pd.read_csv(config['Titeseq_Kds_file'])
dt[['aa_substitutions']] = dt[['aa_substitutions']].fillna(value='')

dt = models.get(('global epistasis', 'Gaussian', 'lib1')).add_phenotypes_to_df(df=dt, substitutions_col='aa_substitutions',latent_phenotype_col='latent_phenotype_Gaussian_1',observed_phenotype_col='predicted_phenotype_Gaussian_1',unknown_as_nan=True)
dt = models.get(('global epistasis', 'Cauchy', 'lib1')).add_phenotypes_to_df(df=dt, substitutions_col='aa_substitutions',latent_phenotype_col='latent_phenotype_Cauchy_1',observed_phenotype_col='predicted_phenotype_Cauchy_1',unknown_as_nan=True)
dt = models.get(('global epistasis', 'Gaussian', 'lib2')).add_phenotypes_to_df(df=dt, substitutions_col='aa_substitutions',latent_phenotype_col='latent_phenotype_Gaussian_2',observed_phenotype_col='predicted_phenotype_Gaussian_2',unknown_as_nan=True)
dt = models.get(('global epistasis', 'Cauchy', 'lib2')).add_phenotypes_to_df(df=dt, substitutions_col='aa_substitutions',latent_phenotype_col='latent_phenotype_Cauchy_2',observed_phenotype_col='predicted_phenotype_Cauchy_2',unknown_as_nan=True)

dt.to_csv(config['global_epistasis_binding_file'], index=False)
print(f"Writing predictions to {config['global_epistasis_binding_file']}")

dt.head().round(2)
```

    Writing predictions to results/global_epistasis_binding/global_epistasis_binding_predictions.csv





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
      <th>target</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>avgcount</th>
      <th>log10Ka</th>
      <th>delta_log10Ka</th>
      <th>log10SE</th>
      <th>Kd</th>
      <th>Kd_SE</th>
      <th>response</th>
      <th>baseline</th>
      <th>nMSR</th>
      <th>variant_class</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>latent_phenotype_Gaussian_1</th>
      <th>predicted_phenotype_Gaussian_1</th>
      <th>latent_phenotype_Cauchy_1</th>
      <th>predicted_phenotype_Cauchy_1</th>
      <th>latent_phenotype_Gaussian_2</th>
      <th>predicted_phenotype_Gaussian_2</th>
      <th>latent_phenotype_Cauchy_2</th>
      <th>predicted_phenotype_Cauchy_2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAAAATTGTAA</td>
      <td>1</td>
      <td>0.00</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&gt;1 nonsynonymous</td>
      <td>Y91L K199Y</td>
      <td>2</td>
      <td>-0.34</td>
      <td>10.67</td>
      <td>-0.33</td>
      <td>10.66</td>
      <td>-0.41</td>
      <td>10.66</td>
      <td>-0.33</td>
      <td>10.69</td>
    </tr>
    <tr>
      <th>1</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAAACTTAAAT</td>
      <td>2</td>
      <td>4.30</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&gt;1 nonsynonymous</td>
      <td>N13S L60P K94N S147T C150Y</td>
      <td>5</td>
      <td>-4.56</td>
      <td>6.61</td>
      <td>-4.76</td>
      <td>6.14</td>
      <td>-4.56</td>
      <td>6.54</td>
      <td>-4.82</td>
      <td>6.06</td>
    </tr>
    <tr>
      <th>2</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAAACTTCAAT</td>
      <td>5</td>
      <td>74.50</td>
      <td>8.72</td>
      <td>-2.05</td>
      <td>0.12</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.64</td>
      <td>1.13</td>
      <td>0.01</td>
      <td>&gt;1 nonsynonymous</td>
      <td>A22C R127G E141D L188V</td>
      <td>4</td>
      <td>-1.54</td>
      <td>9.06</td>
      <td>-1.38</td>
      <td>9.43</td>
      <td>-1.57</td>
      <td>9.07</td>
      <td>-1.23</td>
      <td>9.67</td>
    </tr>
    <tr>
      <th>3</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAACAAGCAGA</td>
      <td>6</td>
      <td>146.32</td>
      <td>10.35</td>
      <td>-0.42</td>
      <td>0.05</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.87</td>
      <td>1.01</td>
      <td>0.00</td>
      <td>1 nonsynonymous</td>
      <td>N13F</td>
      <td>1</td>
      <td>-0.62</td>
      <td>10.46</td>
      <td>-0.67</td>
      <td>10.38</td>
      <td>-0.63</td>
      <td>10.46</td>
      <td>-0.70</td>
      <td>10.36</td>
    </tr>
    <tr>
      <th>4</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAACAATATAA</td>
      <td>1</td>
      <td>2.48</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>&gt;1 nonsynonymous</td>
      <td>C6K T15W K94Y V103W</td>
      <td>4</td>
      <td>-7.48</td>
      <td>6.61</td>
      <td>-6.45</td>
      <td>6.12</td>
      <td>-5.29</td>
      <td>6.53</td>
      <td>-7.16</td>
      <td>6.03</td>
    </tr>
  </tbody>
</table>
</div>



## Repeat fits for pooled library measurements

Repeat the fits for all barcodes pooled together. There is slight variation in the average *K*<sub>D,app</sub> ascribed to wildtype genotypes in each library. To account for this, express the functional score in each library as the delta_log10Ka relative to the average of wildtype measurements in that library, thereby standardizing the small difference in mean WT between the two replicates and avoiding consequential artefacts in a joint fit.


```python
df_joint = pd.read_csv(config['Titeseq_Kds_file'])
df_joint.rename(columns={'delta_log10Ka':'func_score'},inplace=True)
df_joint['func_score_var'] = df_joint['log10SE']**2
func_scores_joint = df_joint[pd.notnull(df_joint['func_score'])]
func_scores_joint.fillna('',inplace=True)
func_scores_joint.head()
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
      <th>target</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>avgcount</th>
      <th>log10Ka</th>
      <th>func_score</th>
      <th>log10SE</th>
      <th>Kd</th>
      <th>Kd_SE</th>
      <th>response</th>
      <th>baseline</th>
      <th>nMSR</th>
      <th>variant_class</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>func_score_var</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAAACTTCAAT</td>
      <td>5</td>
      <td>74.50</td>
      <td>8.72</td>
      <td>-2.05</td>
      <td>0.12</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.64</td>
      <td>1.13</td>
      <td>0.01</td>
      <td>&gt;1 nonsynonymous</td>
      <td>A22C R127G E141D L188V</td>
      <td>4</td>
      <td>0.0144</td>
    </tr>
    <tr>
      <th>3</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAACAAGCAGA</td>
      <td>6</td>
      <td>146.32</td>
      <td>10.35</td>
      <td>-0.42</td>
      <td>0.05</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.87</td>
      <td>1.01</td>
      <td>0.00</td>
      <td>1 nonsynonymous</td>
      <td>N13F</td>
      <td>1</td>
      <td>0.0025</td>
    </tr>
    <tr>
      <th>5</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAACAGGTTGC</td>
      <td>4</td>
      <td>47.48</td>
      <td>6.00</td>
      <td>-4.76</td>
      <td>2.00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1.50</td>
      <td>1.15</td>
      <td>0.01</td>
      <td>&gt;1 nonsynonymous</td>
      <td>V71K P149L N157T</td>
      <td>3</td>
      <td>4.0000</td>
    </tr>
    <tr>
      <th>6</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAACATTAAAT</td>
      <td>6</td>
      <td>47.44</td>
      <td>10.15</td>
      <td>-0.61</td>
      <td>0.09</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.64</td>
      <td>1.09</td>
      <td>0.00</td>
      <td>&gt;1 nonsynonymous</td>
      <td>A18V T148S H189Y</td>
      <td>3</td>
      <td>0.0081</td>
    </tr>
    <tr>
      <th>8</th>
      <td>lib1</td>
      <td>SARS-CoV-2</td>
      <td>AAAAAAAACCTTACAA</td>
      <td>2</td>
      <td>18.78</td>
      <td>9.61</td>
      <td>-1.15</td>
      <td>0.19</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>2.00</td>
      <td>1.11</td>
      <td>0.01</td>
      <td>&gt;1 nonsynonymous</td>
      <td>T63D A89N</td>
      <td>2</td>
      <td>0.0361</td>
    </tr>
  </tbody>
</table>
</div>




```python
# NBVAL_IGNORE_OUTPUT

models_joint = {}  # store models, keyed by `(epistasistype, likelihoodtype)`

for (target), scores in func_scores_joint.groupby(['target']):
   
    bmap = dms_variants.binarymap.BinaryMap(scores)
    
    for epistasistype, likelihoodtype, Model in [
            ('global epistasis', 'Gaussian', dms_variants.globalepistasis.MonotonicSplineEpistasisGaussianLikelihood),
            ('no epistasis', 'Gaussian', dms_variants.globalepistasis.NoEpistasisGaussianLikelihood),
            ('global epistasis', 'Cauchy', dms_variants.globalepistasis.MonotonicSplineEpistasisCauchyLikelihood),
            ('no epistasis', 'Cauchy', dms_variants.globalepistasis.NoEpistasisCauchyLikelihood),
            ]:
        print(f"Fitting {epistasistype} with {likelihoodtype} likelihood model...", end=' ')
    
        start = time.time()
        model = Model(bmap)
        model.fit()  # do NOT change ftol in normal use, this is just for test
        print(f"fitting took {time.time() - start:.1f} sec.")
        models_joint[(epistasistype, likelihoodtype)] = model

```

    Fitting global epistasis with Gaussian likelihood model... fitting took 156.0 sec.
    Fitting no epistasis with Gaussian likelihood model... fitting took 6.1 sec.
    Fitting global epistasis with Cauchy likelihood model... fitting took 128.1 sec.
    Fitting no epistasis with Cauchy likelihood model... fitting took 30.1 sec.



```python
# NBVAL_IGNORE_OUTPUT

variants_df_joint = pd.concat(
        [model.phenotypes_df
         .assign(likelihoodtype=likelihoodtype)
         for (epistasistype, likelihoodtype), model in models_joint.items()
         if (epistasistype == 'global epistasis')],
        ignore_index=True, sort=False)

variants_df_joint.head().round(2)
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
      <th>aa_substitutions</th>
      <th>func_score</th>
      <th>func_score_var</th>
      <th>latent_phenotype</th>
      <th>observed_phenotype</th>
      <th>likelihoodtype</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>A22C R127G E141D L188V</td>
      <td>-2.05</td>
      <td>0.01</td>
      <td>-1.60</td>
      <td>-1.75</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>1</th>
      <td>N13F</td>
      <td>-0.42</td>
      <td>0.00</td>
      <td>-0.64</td>
      <td>-0.33</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>2</th>
      <td>V71K P149L N157T</td>
      <td>-4.76</td>
      <td>4.00</td>
      <td>-6.19</td>
      <td>-4.19</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>3</th>
      <td>A18V T148S H189Y</td>
      <td>-0.61</td>
      <td>0.01</td>
      <td>-0.92</td>
      <td>-0.67</td>
      <td>Gaussian</td>
    </tr>
    <tr>
      <th>4</th>
      <td>T63D A89N</td>
      <td>-1.15</td>
      <td>0.04</td>
      <td>-1.31</td>
      <td>-1.27</td>
      <td>Gaussian</td>
    </tr>
  </tbody>
</table>
</div>




```python
for x, y in itertools.combinations(['latent_phenotype',
                                    'observed_phenotype',
                                    'func_score'],
                                   2):
    p = (
        ggplot(variants_df_joint, aes(x, y)) +
        geom_point(alpha=0.05, size=0.5))
#        facet_grid('likelihoodtype', scales='free_y') +
#        theme(figure_size=(2 * variants_df_novar['likelihoodtype'].nunique()),
#              )
#        )
#         plotfile = os.path.join(config['figs_dir'], f'{y}-v-{x}_by_{epistasistype}.pdf')
#         print(f"Saving to {plotfile}")
#         p.save(plotfile)
    _ = p.draw()
```


![png](global_epistasis_binding_files/global_epistasis_binding_33_0.png)



![png](global_epistasis_binding_files/global_epistasis_binding_33_1.png)



![png](global_epistasis_binding_files/global_epistasis_binding_33_2.png)


We will look more at each of these models (Cauchy and Gaussian likelihoods, observed and latent scale measurements) in the next notebook, in which we evaluate coefficients, process our final phenotype map, and compare our measurements to some validation datasets.

## Output epistasis model parameters



```python
#lib1 models
models.get(('global epistasis', 'Gaussian', 'lib1')).single_mut_effects(phenotype='latent',standardize_range=False).to_csv('results/global_epistasis_binding/Gaussian-latent-effects_binding_1.csv',index=False)
models.get(('global epistasis', 'Gaussian', 'lib1')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/Gaussian-predicted-effects_binding_1.csv',index=False)
models.get(('global epistasis', 'Cauchy', 'lib1')).single_mut_effects(phenotype='latent',standardize_range=False).to_csv('results/global_epistasis_binding/Cauchy-latent-effects_binding_1.csv',index=False)
models.get(('global epistasis', 'Cauchy', 'lib1')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/Cauchy-predicted-effects_binding_1.csv',index=False)
models.get(('no epistasis', 'Gaussian', 'lib1')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/nonepistatic-Gaussian-predicted-effects_binding_1.csv',index=False)
models.get(('no epistasis', 'Cauchy', 'lib1')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/nonepistatic-Cauchy-predicted-effects_binding_1.csv',index=False)
#lib2 models
models.get(('global epistasis', 'Gaussian', 'lib2')).single_mut_effects(phenotype='latent',standardize_range=False).to_csv('results/global_epistasis_binding/Gaussian-latent-effects_binding_2.csv',index=False)
models.get(('global epistasis', 'Gaussian', 'lib2')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/Gaussian-predicted-effects_binding_2.csv',index=False)
models.get(('global epistasis', 'Cauchy', 'lib2')).single_mut_effects(phenotype='latent',standardize_range=False).to_csv('results/global_epistasis_binding/Cauchy-latent-effects_binding_2.csv',index=False)
models.get(('global epistasis', 'Cauchy', 'lib2')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/Cauchy-predicted-effects_binding_2.csv',index=False)
models.get(('no epistasis', 'Gaussian', 'lib2')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/nonepistatic-Gaussian-predicted-effects_binding_2.csv',index=False)
models.get(('no epistasis', 'Cauchy', 'lib2')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/nonepistatic-Cauchy-predicted-effects_binding_2.csv',index=False)
#joint models
models_joint.get(('global epistasis', 'Gaussian')).single_mut_effects(phenotype='latent',standardize_range=False).to_csv('results/global_epistasis_binding/Gaussian-latent-effects_binding_joint.csv',index=False)
models_joint.get(('global epistasis', 'Gaussian')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/Gaussian-predicted-effects_binding_joint.csv',index=False)
models_joint.get(('global epistasis', 'Cauchy')).single_mut_effects(phenotype='latent',standardize_range=False).to_csv('results/global_epistasis_binding/Cauchy-latent-effects_binding_joint.csv',index=False)
models_joint.get(('global epistasis', 'Cauchy')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/Cauchy-predicted-effects_binding_joint.csv',index=False)
models_joint.get(('no epistasis', 'Gaussian')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/nonepistatic-Gaussian-predicted-effects_binding_joint.csv',index=False)
models_joint.get(('no epistasis', 'Cauchy')).single_mut_effects(phenotype='observed',standardize_range=False).to_csv('results/global_epistasis_binding/nonepistatic-Cauchy-predicted-effects_binding_joint.csv',index=False)



```


```python
#! jupyter nbconvert --to markdown global_epistasis_binding.ipynb --output-dir ./results/summary/ --output global_epistasis_binding.md
```


```python

```
