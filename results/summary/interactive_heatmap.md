# Interactive heat maps of effects of mutations
Use [Altair](https://altair-viz.github.io/) to create interactive heat maps of mutational effects with relevant annotations.

## Imports and read data



```python
import math
import os
import altair as alt
import pandas as pd
import numpy as np
import yaml
alt.data_transformers.enable('default', max_rows=None)
```




    DataTransformerRegistry.enable('default')



Read in the configuration file, and then read the input data files from that:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
    
os.makedirs(config['figs_dir'], exist_ok=True)

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
site_annotations['SARS-CoV-1/2 differences'] = site_annotations[['amino_acid_SARS2', 'amino_acid_SARS1']].apply(lambda x: False if x[0] == x[1] else True, axis=1)
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
      <th>epitope_CR3022</th>
      <th>epitope_VHH72</th>
      <th>epitope_S230</th>
      <th>epitope_m396</th>
      <th>epitope_F26G19</th>
      <th>epitope_80R</th>
      <th>epitope_B38</th>
      <th>epitope_S309</th>
      <th>buried_downRBD</th>
      <th>SARS-CoV-1/2 differences</th>
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
      <td>True</td>
      <td>False</td>
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
      <td>True</td>
      <td>False</td>
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
      <td>True</td>
      <td>False</td>
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
      <td>True</td>
      <td>False</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 32 columns</p>
</div>



## Data frame for heat maps
This data frame needs to have the mutational effects along with the additional site annotations we'd like to show in interactive format:


```python
# mut effects to keep and their names
mut_effects_to_keep = {'site_SARS2': 'site',
                       'mutant': 'mutant',
                       'wildtype': 'wildtype',
                       'mutation': 'mutation',
                       'bind_avg': 'ACE2_binding',
                       'expr_avg': 'expression',
                      }

# site annotations to keep and their names
site_annotations_to_keep = {'site_SARS2': 'site',
                            'RSA_bound': 'RSA_bound',
                            'SARS2_ACE2_contact': 'ACE2 contact SARS-CoV-2',
                            'SARS1_ACE2_contact': 'ACE2 contact SARS-CoV-1',
                            'amino_acid_SARS1': 'SARS-CoV-1 aa',
                            'amino_acid_RaTG13': 'RaTG13 aa',
                            'amino_acid_GD_Pangolin': 'GD-Pangolin aa',
                            'SARS-CoV-1/2 differences': 'SARS-CoV-1/2 differences',
                            'epitope_B38': 'epitope B38',
                            'epitope_S309': 'epitope S309',
                            'epitope_CR3022': 'epitope CR3022',
                            'epitope_VHH72': 'epitope VHH-72'
                           }

# conditions is the data shown in the heatmaps
# targets are the values in the dropdown menus
conditions = ['expression', 'ACE2 binding']
targets = ['all sites', 'ACE2 contact SARS-CoV-2', 'ACE2 contact SARS-CoV-1', 'SARS-CoV-1/2 differences',
           'epitope B38', 'epitope S309', 'epitope CR3022', 'epitope VHH-72']

# merge the mutational effects and the site annotations
df = (mut_effects[list(mut_effects_to_keep)]
      .rename(columns=mut_effects_to_keep)
      .query('(mutant != "*") and (wildtype != "*")')
      .merge(site_annotations[list(site_annotations_to_keep)]
             .rename(columns=site_annotations_to_keep),
             on='site')
     )
# format RSA
df['RSA_bound'] = df['RSA_bound'].round(2)

# format 'missing data' codes
df['wildtype_code'] = (df[['wildtype', 'mutant']].apply(lambda x: 'x' if x[0] == x[1] else '', axis=1))
df['RSA_code'] = (df['RSA_bound'].apply(lambda x: '-' if np.isnan(x) else ''))

# we want one dropdown site subset to be all sites
df['all sites'] = True

# pull out data used for the zoom bar
zoom_bar_info = df[['site', 'ACE2 contact SARS-CoV-2']]

# melt data frame to allow degeneracy in the dropdown mneu 
df = pd.melt(df,
            id_vars=[x for x in df.columns.values if x not in targets],
            var_name='subset')
df = df.drop(df[df['value'] == False].index).drop(columns='value')
df = pd.merge(df, zoom_bar_info, on='site')

# clean up 
df = df.drop(columns=['wildtype']).drop_duplicates()
df = df.fillna(np.nan)
df.head()
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
      <th>site</th>
      <th>mutant</th>
      <th>mutation</th>
      <th>ACE2_binding</th>
      <th>expression</th>
      <th>RSA_bound</th>
      <th>SARS-CoV-1 aa</th>
      <th>RaTG13 aa</th>
      <th>GD-Pangolin aa</th>
      <th>wildtype_code</th>
      <th>RSA_code</th>
      <th>subset</th>
      <th>ACE2 contact SARS-CoV-2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>417</td>
      <td>A</td>
      <td>K417A</td>
      <td>-0.35</td>
      <td>0.03</td>
      <td>0.19</td>
      <td>V</td>
      <td>K</td>
      <td>R</td>
      <td></td>
      <td></td>
      <td>ACE2 contact SARS-CoV-2</td>
      <td>True</td>
    </tr>
    <tr>
      <th>20</th>
      <td>417</td>
      <td>C</td>
      <td>K417C</td>
      <td>-0.42</td>
      <td>-0.05</td>
      <td>0.19</td>
      <td>V</td>
      <td>K</td>
      <td>R</td>
      <td></td>
      <td></td>
      <td>ACE2 contact SARS-CoV-2</td>
      <td>True</td>
    </tr>
    <tr>
      <th>40</th>
      <td>417</td>
      <td>D</td>
      <td>K417D</td>
      <td>-1.04</td>
      <td>-0.33</td>
      <td>0.19</td>
      <td>V</td>
      <td>K</td>
      <td>R</td>
      <td></td>
      <td></td>
      <td>ACE2 contact SARS-CoV-2</td>
      <td>True</td>
    </tr>
    <tr>
      <th>60</th>
      <td>417</td>
      <td>E</td>
      <td>K417E</td>
      <td>-0.75</td>
      <td>-0.25</td>
      <td>0.19</td>
      <td>V</td>
      <td>K</td>
      <td>R</td>
      <td></td>
      <td></td>
      <td>ACE2 contact SARS-CoV-2</td>
      <td>True</td>
    </tr>
    <tr>
      <th>80</th>
      <td>417</td>
      <td>F</td>
      <td>K417F</td>
      <td>-0.13</td>
      <td>-0.11</td>
      <td>0.19</td>
      <td>V</td>
      <td>K</td>
      <td>R</td>
      <td></td>
      <td></td>
      <td>ACE2 contact SARS-CoV-2</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
</div>



## Create plot functions


```python
def RSA_bar(data, metric, missing_metric):
    """Create one layer heatmap of RSA data.
    Parameters
    ----------
    data :pandas.DataFrame
        Data frame with site, RSA value, and missing data code.
    metric : str
        Column in `data` with RSA values.
    missing_metric : str
        Column in `data` with codes for missing values.
    Returns
    -------
    altair.Chart
    """
    RSA_base = (alt.Chart(data)
        .encode(x=alt.X('site:O',
                        axis=alt.Axis(title=None,
                                      ticks=False,
                                      labels=False))))
    RSA = (RSA_base
               .mark_rect()
               .encode(color=alt.Color(metric,
                                       scale=alt.Scale(scheme='greys',
                                                       domain=[1.0, 0.0]),
                                       legend=None),
                      ).interactive()
              )

    RSA_missing = (RSA_base
                .mark_text(color='black')
                .encode(text=alt.Text(missing_metric)
                       )
               )

    RSA = ((RSA + RSA_missing).properties(title={'text': ['relative solvent accessibility'],
                                                  'subtitle': ['(white=1, black=0)']})
                .add_selection(subset_select)  # add dropdown menu
                .transform_filter(subset_select)  # add dropdown filtering
                .transform_filter(zoom_brush)
                )
    return RSA
```


```python
def zoom_bar(data, metric, title):
    """Create one layer heatmap for zoom bar.
    Parameters
    ----------
    data :pandas.DataFrame
        Data frame with site and metric value.
    metric : str
        Column in `data` with values to color by.
    title : str
        Title of the plot.
    Returns
    -------
    altair.Chart
    """
    zoom = (alt.Chart(data)
            .mark_rect()
            .encode(x='site:O',
                    color=alt.Color(metric, 
                                     scale=alt.Scale(domain=[True, False],  
                                                    range=['lightgrey',
                                                           'darkgrey']),
                                    legend=alt.Legend(orient='left',
                                                     labelFontSize=15,
                                                     titleFontSize=16,
                                                     title=title)))
            .add_selection(zoom_brush)
            .add_selection(subset_select)
            .transform_filter(subset_select)
            .properties(width=900,
                        title='zoom bar'))
    return zoom
```


```python
def DMS_heatmaps(data, metric):
    """Create main heatmap for one condition.
    The heatmap is the results of three layers.
    *heatmap* is the main DMS data
    *wildtype* marks wildtype data with an 'x'
    *nulls* creates grey cells for missing data.
    If you exclude nulls, missing data is white, 
    which is appropriate for some color schemes
    but not all.
    Parameters
    ----------
    data :pandas.DataFrame
        Main dataframe
    metric : str
        Column in `data` with values to color by.
    Returns
    -------
    altair.Chart
    """
    aa_order = ['R', 'K', 'H', 'D', 'E', 'Q', 'N', 'S', 'T', 'Y',
            'W', 'F', 'A', 'I', 'L', 'M', 'V', 'G', 'P', 'C', '*']
    tooltips = [c for c in data.columns if c not in
                {'site', 'wildtype', 'mutant', 'wildtype_code', 'subset', 'RSA_code'}]

    # everything is site v mutant
    base = (alt.Chart(data)
            .encode(x=alt.X('site:O',
                             axis=alt.Axis(titleFontSize=15)),
                    y=alt.Y('mutant:O',
                            sort=aa_order,
                            axis=alt.Axis(labelFontSize=12,
                                          titleFontSize=15))
                   )
           )
    heatmap = (base
               .mark_rect()
               .encode(color=alt.Color(metric,
                                       type='quantitative', 
                                       scale=alt.Scale(scheme='redblue',
                                                       domain=[data[metric].min(),
                                                               data[metric].max()],
                                                       domainMid=0),
                                       legend=alt.Legend(orient='left',
                                                         title='grey is n.d.',
                                                        gradientLength=100)),
                       stroke=alt.value('black'),
                       strokeWidth=alt.condition(cell_selector,
                                                 alt.value(2),
                                                 alt.value(0)),
                       tooltip=tooltips
                      )
              )
    
    wildtype = (base
                .mark_text(color='black')
                .encode(text=alt.Text('wildtype_code:N')
                       )
               )

    nulls = (base
             .mark_rect()
             .transform_filter(f"!isValid(datum.{metric})")
             .mark_rect(opacity=0.5)
             .encode(alt.Color(f'{metric}:N',
                               scale=alt.Scale(scheme='greys'),
                               legend=None)
                    )
            )
    
    return ((heatmap + nulls + wildtype)
            .interactive()
            .add_selection(subset_select)  # add dropdown menu
            .add_selection(cell_selector)  # mouse over highlighting
            .transform_filter(subset_select)  # add dropdown filtering
            .transform_filter(zoom_brush)  # add zoom bar filtering
            .properties(height=250, title=' '.join(metric.split('_'))))
```

## Create plots


```python
# SELECTIONS
# select the cell in the other heatmap
cell_selector = alt.selection_single(on='mouseover',
                                     empty='none')

# A dropdown filter
subset_dropdown = alt.binding_select(options=targets)
subset_select = alt.selection_single(fields=['subset'],
                                     bind=subset_dropdown,
                                     name="site",
                                     init={'subset': 'all sites'})
# zoom brush
zoom_brush = alt.selection_interval(encodings=['x'], mark=alt.BrushConfig(stroke='black',
                                                                          strokeWidth=2))

# create each plot
RSA = RSA_bar(df, 'RSA_bound', 'RSA_code')
expression = DMS_heatmaps(df, 'expression')
binding = DMS_heatmaps(df, 'ACE2_binding')
zoom = zoom_bar(df, 'ACE2 contact SARS-CoV-2', 'ACE2 contact')

# paste them together
chart = alt.vconcat(binding, expression, spacing=0)
chart = alt.vconcat(chart, RSA, spacing=0)

chart = (alt.vconcat(zoom, chart, spacing=0)
         .configure_title(anchor='start',
                          fontSize=20))
chart
```





<div id="altair-viz-b7a7352724fd45fea6c28a2f7eff3bd0"></div>
<script type="text/javascript">
  (function(spec, embedOpt){
    let outputDiv = document.currentScript.previousElementSibling;
    if (outputDiv.id !== "altair-viz-b7a7352724fd45fea6c28a2f7eff3bd0") {
      outputDiv = document.getElementById("altair-viz-b7a7352724fd45fea6c28a2f7eff3bd0");
    }
    const paths = {
      "vega": "https://cdn.jsdelivr.net/npm//vega@5?noext",
      "vega-lib": "https://cdn.jsdelivr.net/npm//vega-lib?noext",
      "vega-lite": "https://cdn.jsdelivr.net/npm//vega-lite@4.8.1?noext",
      "vega-embed": "https://cdn.jsdelivr.net/npm//vega-embed@6?noext",
    };

    function loadScript(lib) {
      return new Promise(function(resolve, reject) {
        var s = document.createElement('script');
        s.src = paths[lib];
        s.async = true;
        s.onload = () => resolve(paths[lib]);
        s.onerror = () => reject(`Error loading script: ${paths[lib]}`);
        document.getElementsByTagName("head")[0].appendChild(s);
      });
    }

    function showError(err) {
      outputDiv.innerHTML = `<div class="error" style="color:red;">${err}</div>`;
      throw err;
    }

    function displayChart(vegaEmbed) {
      vegaEmbed(outputDiv, spec, embedOpt)
        .catch(err => showError(`Javascript Error: ${err.message}<br>This usually means there's a typo in your chart specification. See the javascript console for the full traceback.`));
    }

    if(typeof define === "function" && define.amd) {
      requirejs.config({paths});
      require(["vega-embed"], displayChart, err => showError(`Error loading script: ${err.message}`));
    } else if (typeof vegaEmbed === "function") {
      displayChart(vegaEmbed);
    } else {
      loadScript("vega")
        .then(() => loadScript("vega-lite"))
        .then(() => loadScript("vega-embed"))
        .catch(showError)
        .then(() => displayChart(vegaEmbed));
    }
</script>




```python
print(f"Saving chart to {config['interactive_heatmap']}")
os.makedirs(os.path.dirname(config['interactive_heatmap']), exist_ok=True)

chart.save(config['interactive_heatmap'])
```

    Saving chart to docs/_includes/interactive_heatmap.html



```python

```