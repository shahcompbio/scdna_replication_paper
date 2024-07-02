---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: spectrum_scdna
    language: python
    name: spectrum_scdna
---

```python

import anndata as ad
import pandas as pd
import scgenome


```

```python

pdxs = [
    'SA1047',
    'SA1049',
    'SA1050',
    'SA1051',
    'SA1052',
    'SA1053',
    'SA1091',
    'SA1093',
    'SA1096',
    'SA1162',
    'SA1181',
    'SA1182',
    'SA1184',
    'SA501',
    'SA530',
    'SA604',
]

cell_lines = [
    'OV2295',
    'SA039',
    'SA1054',
    'SA1055',
    'SA1056',
    'SA1188',
    'SA1292',
    'SA906a',
    'SA906b',
]

cell_frac_rep = []

for dataset in pdxs:
    filename = f'/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz'
    cell_frac_rep.append(pd.read_csv(filename, usecols=['cell_id', 'clone_id', 'cell_frac_rep']).drop_duplicates(['cell_id', 'clone_id']).assign(dataset=dataset))

for dataset in cell_lines:
    filename = f'/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv'
    cell_frac_rep.append(pd.read_csv(filename, sep='\t', usecols=['cell_id', 'clone_id', 'cell_frac_rep']).drop_duplicates(['cell_id', 'clone_id']).assign(dataset=dataset))

filename = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/laks_flow/GM18507/cn_s_pyro_inferred_composite_filtered.tsv'
cell_frac_rep.append(pd.read_csv(filename, sep='\t', usecols=['cell_id', 'clone_id', 'cell_frac_rep']).drop_duplicates(['cell_id', 'clone_id']).assign(dataset='GM18507'))

filename = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/laks_flow/T47D/cn_s_pyro_inferred_composite_filtered.tsv'
cell_frac_rep.append(pd.read_csv(filename, sep='\t', usecols=['cell_id', 'clone_id', 'cell_frac_rep']).drop_duplicates(['cell_id', 'clone_id']).assign(dataset='T47D'))

cell_frac_rep = pd.concat(cell_frac_rep)

```

```python

# For both files, the fraction of replicated bins for that cell will be noted as cell_frac_rep with values ->1 representing
# fully replicated late S, and ->0 representing early S

cell_frac_rep

```

```python

features = pd.read_csv('/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/rt_model/features.csv.gz', low_memory=False)

features['clone'] = features['dataset'] + '_clone' + features['clone_id']
features = features.set_index('clone')
features['signature_NA'] = 1 - features.filter(regex='signature_.*', axis=1).sum(axis=1)
features['wgd'] = (features['ploidy'] > 2) * 1

features_ploidy = features['ploidy'].astype('int').astype('str').to_frame()

features_wgd = features['wgd'].to_frame()

features_signature = features.filter(regex='signature_.*', axis=1).melt(var_name='signature', value_name='indicator', ignore_index=False)
features_signature = features_signature[features_signature['indicator'] == 1]
features_signature['signature'] = features_signature['signature'].str.replace('signature_', '')
features_signature = features_signature.drop('indicator', axis=1)

features_type = features.filter(regex='type_.*', axis=1).melt(var_name='type', value_name='indicator', ignore_index=False)
features_type = features_type[features_type['indicator'] == 1]
features_type['type'] = features_type['type'].str.replace('type_', '')
features_type = features_type.drop('indicator', axis=1)

features_df = pd.concat([features_type, features_signature, features_ploidy, features_wgd], axis=1).fillna('N/A')
features_df.index.name = 'clone'

cell_frac_rep['clone'] = cell_frac_rep['dataset'] + '_clone' + cell_frac_rep['clone_id']

```

```python

cell_frac_rep = cell_frac_rep.merge(features_df.reset_index())
cell_frac_rep

```

```python

import seaborn as sns

plot_data = cell_frac_rep

remove = ['SA1162_cloneB', 'SA1162_cloneC', 'SA1162_cloneD', 'SA1162_cloneG', 'SA1052_cloneG']

plot_data = plot_data[~plot_data['clone'].isin(remove)]
# plot_data = plot_data[~plot_data['clone'].str.startswith('SA1091')]
plot_data = plot_data.merge(features.reset_index()[['clone', 'num_cells_s']])
plot_data['clone'] = plot_data['clone'] + ' n_s=' + plot_data['num_cells_s'].astype(str)
plot_data = plot_data[plot_data['num_cells_s'] >= 20]

sns.displot(data=plot_data.query('signature == "FBI"'), x='cell_frac_rep', col='clone', kind='hist', col_wrap=4, common_norm=False, stat='proportion', bins=20)

```

```python

sns.displot(data=plot_data, x='cell_frac_rep', col='signature', kind='kde', common_norm=False)#, common_norm=False, stat='proportion', bins=20)

```

```python

import matplotlib.pyplot as plt

plt.figure(figsize=(2,2))
sns.barplot(data=plot_data, x='signature', y='cell_frac_rep', color='0.5', width=0.4)
sns.despine()

plt.figure(figsize=(5,2))
sns.barplot(data=plot_data, x='type', y='cell_frac_rep', color='0.5', width=0.4)
sns.despine()

plt.figure(figsize=(3,2))
sns.barplot(data=plot_data, x='signature', hue='type', y='cell_frac_rep')
plt.legend(ncols=2, loc='upper left', bbox_to_anchor=(1., 1.))
sns.despine()

```

```python

import matplotlib.pyplot as plt

plt.figure(figsize=(1,2))
sns.barplot(data=plot_data, x='wgd', y='cell_frac_rep', color='0.5', width=0.4)
sns.despine()

```

```python

```
