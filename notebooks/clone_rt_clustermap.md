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
    display_name: Spectrum new scgenome
    language: python
    name: spectrum_newscgenome
---

```python

import pandas as pd
import matplotlib.pyplot as plt

clone_rt = pd.read_csv('/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/rt_model/clone_rt.csv.gz', low_memory=False)
features = pd.read_csv('/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/rt_model/features.csv.gz', low_memory=False)

features['clone'] = features['dataset'] + '_clone' + features['clone_id']
features = features.set_index('clone')
features['signature_NA'] = 1 - features.filter(regex='signature_.*', axis=1).sum(axis=1)
features['wgd'] = (features['ploidy'] > 2) * 1

clone_rt = clone_rt.set_index(['chr', 'start', 'end']).T
clone_rt.index = [a.replace('pseudobulk_', '').replace('_model_rep_state', '') for a in clone_rt.index]

remove = ['SA1162_cloneB', 'SA1162_cloneC', 'SA1162_cloneD', 'SA1162_cloneG', 'SA1052_cloneG']

clone_rt = clone_rt[~clone_rt.index.isin(remove)]

features = features.reindex(clone_rt.index)
features.index.name = 'clone'

```

```python
features
```

```python

columns = [
    'chr',
    'start',
    'end',
    'mcf7_rt',
    'bg02es_rt',
    'bj_rt',
    'gm06990_rt',
    'gm12801_rt',
    'gm12812_rt',
    'gm12813_rt',
    'gm12878_rt',
    'helas3_rt',
    'hepg2_rt',
    'huvec_rt',
    'imr90_rt',
    'k562_rt',
    'sknsh_rt',
    'nhek_rt',
]

encode_filename = '/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/simulation/D1.0/s_phase_cells_pyro_composite_filtered.tsv'
encode = pd.read_csv(encode_filename, sep='\t', low_memory=True, usecols=columns).drop_duplicates()
encode['start'] += 1
encode = encode.set_index(['chr', 'start', 'end']).T
encode.index

```

```python

features['num_cells_s'].hist(bins=100)

features = features[features['num_cells_s'] >= 20]

clone_rt = clone_rt.reindex(features.index)

```

```python

clone_rt.head()

```

```python

features_ploidy = features['ploidy'].astype('int').astype('str').to_frame()
features_ploidy

```

```python

features_wgd = features['wgd'].to_frame()
features_wgd.head()

```

```python

assert (features.filter(regex='signature_.*', axis=1).sum(axis=1) == 1).all()

features_signature = features.filter(regex='signature_.*', axis=1).melt(var_name='signature', value_name='indicator', ignore_index=False)
features_signature = features_signature[features_signature['indicator'] == 1]
features_signature['signature'] = features_signature['signature'].str.replace('signature_', '')
features_signature = features_signature.drop('indicator', axis=1)
features_signature = features_signature.reindex(clone_rt.index, fill_value='NA')
features_signature

```

```python

features_type = features.filter(regex='type_.*', axis=1).melt(var_name='type', value_name='indicator', ignore_index=False)
features_type = features_type[features_type['indicator'] == 1]
features_type['type'] = features_type['type'].str.replace('type_', '')
features_type = features_type.drop('indicator', axis=1)
features_type

```

```python

features_dataset = features['dataset'].to_frame()
features_dataset.head()

```

```python

import seaborn as sns
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc

# Get unique attributes and sort them
def get_colors(df, palettes):
    colors = df[palettes.keys()].copy()
    attribute_to_color = dict()
    for col in palettes.keys():
        unique_attrs = sorted(colors[col].unique())
        cmap = sns.color_palette(palette=palettes[col], n_colors=len(unique_attrs))
        attribute_to_color[col] = dict(zip(unique_attrs, cmap))
        colors[col] = colors[col].map(attribute_to_color[col])
    return colors, attribute_to_color

features_df = pd.concat([features_type, features_signature, features_dataset, features_ploidy, features_wgd], axis=1).fillna('N/A')
features_df.index.name = 'clone'

palettes = {
    'type': 'tab10',
    'signature': 'bright',
    'wgd': 'mako',
    'dataset': 'tab20',
}

col_colors, attribute_to_color = get_colors(features_df, palettes)

clone_rt_corr = clone_rt.T.corr()
clone_rt_dism = 1 - clone_rt_corr
linkage = hc.linkage(sp.distance.squareform(clone_rt_dism), method='average')
g = sns.clustermap(clone_rt_corr, row_linkage=linkage, col_linkage=linkage, cmap=sns.cm.rocket, col_colors=col_colors, figsize=(8, 8))
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_yticks([])

for feature in attribute_to_color.keys():
    plt.figure(figsize=(4, 1))
    ax = plt.gca()
    for attribute, color in attribute_to_color[feature].items():
        ax.bar(0, 0, color=color, label=attribute, linewidth=0)
    ax.legend(loc='center', ncols=7, bbox_to_anchor=(0.5, 0.5), title=feature)
    ax.axis('off')

```

```python

clone_rt_dism.mean(axis=1).hist(bins=100)

outliers = clone_rt_dism.mean(axis=1).rename('mean_dist')
outliers = outliers[outliers > 0.4]
outliers.to_frame()

```

```python

import umap
import seaborn as sns
import matplotlib.pyplot as plt

embedding = umap.UMAP(random_state=10).fit_transform(clone_rt.values)
plot_data = pd.DataFrame(index=clone_rt.index)
plot_data['UMAP1'] = embedding[:, 0]
plot_data['UMAP2'] = embedding[:, 1]

plot_data = plot_data.merge(features[['dataset', 'num_cells_s']], left_index=True, right_index=True)
plot_data = plot_data.merge(features_signature, left_index=True, right_index=True)
plot_data = plot_data.merge(features_type, left_index=True, right_index=True)
plot_data = plot_data.merge(features_ploidy, left_index=True, right_index=True)
plot_data = plot_data.merge(features_wgd, left_index=True, right_index=True)

plt.figure()
sns.scatterplot(x='UMAP1', y='UMAP2', data=plot_data, hue='type', s=plot_data['num_cells_s'] * 0.1)
plt.legend(ncols=4, loc='upper left', bbox_to_anchor=(1., 1.))

plt.figure()
sns.scatterplot(x='UMAP1', y='UMAP2', data=plot_data, hue='signature', s=plot_data['num_cells_s'] * 0.1)
plt.legend(ncols=4, loc='upper left', bbox_to_anchor=(1., 1.))

plt.figure()
sns.scatterplot(x='UMAP1', y='UMAP2', data=plot_data, hue='wgd', s=plot_data['num_cells_s'] * 0.1)
plt.legend(ncols=4, loc='upper left', bbox_to_anchor=(1., 1.))

plt.figure()
sns.scatterplot(x='UMAP1', y='UMAP2', data=plot_data, hue='dataset', s=plot_data['num_cells_s'] * 0.1)
plt.legend(ncols=4, loc='upper left', bbox_to_anchor=(1., 1.))

```

```python

# features = features[(features['type_HGSOC'] == 1) | (features['type_TNBC'] == 1)]
# clone_rt = clone_rt.loc[features.index]

```

```python

import numpy as np
from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

feature_cols = [
    'wgd',
    'type_GM18507',
    'type_HGSOC',
    'type_OV2295',
    'type_T47D',
    'type_TNBC',
    'type_hTERT',
    'signature_FBI',
    'signature_HRD',
    'signature_TD',
    'signature_NA',
]

scaler = StandardScaler()

selected_features = features[feature_cols]

observations = scaler.fit_transform(clone_rt.values[:, :])
covariates = scaler.fit_transform(selected_features.values)

n_comp = 3

# Apply CCA
cca = CCA(n_components=n_comp)
observations_c, covariates_c = cca.fit_transform(observations, covariates)

comp_corr = [np.corrcoef(observations_c[:, i], covariates_c[:, i])[1][0] for i in range(n_comp)]
plt.bar(range(n_comp), comp_corr, color='lightgrey', width = 0.8, edgecolor='k')

```

```python

import rcca

n_comp = 6

# Initialize a CCA object as regularized or sparse
cca = rcca.CCA(kernelcca=True, reg=1, numCC=n_comp)

# Train on data
# X and Y are the two datasets (e.g., features and covariates)
cca.train([observations, covariates])

```

```python

comp_corr = [np.corrcoef(cca.comps[0][:, i], cca.comps[1][:, i])[1][0] for i in range(n_comp)]
plt.bar(range(n_comp), comp_corr, color='lightgrey', width = 0.8, edgecolor='k')
_ = plt.xticks(range(n_comp))

```

```python

weights = pd.DataFrame(cca.ws[1], index=selected_features.columns)

weights[3].plot.bar()

```

```python

selected_features.groupby(['wgd', 'type_HGSOC', 'signature_HRD']).size()

```

```python

Y = clone_rt.copy()

chromosomes = [str(a) for a in range(1, 23)] + ['X']

plot_data = Y.T.groupby(level=0).mean().melt(ignore_index=False, var_name='clone', value_name='mean_rt').reset_index()
plot_data = plot_data.merge(features_df.reset_index())

plt.figure(figsize=(20, 4))
sns.barplot(x='chr', hue='type', y='mean_rt', order=chromosomes, data=plot_data)

plt.figure(figsize=(4, 4))
sns.barplot(x='type', y='mean_rt', data=plot_data)

```

<!-- #region jp-MarkdownHeadingCollapsed=true -->

# Clone rt histograms

<!-- #endregion -->

```python

plot_data = (
    clone_rt.T.melt(ignore_index=False, var_name='clone', value_name='rt')
        .reset_index().merge(features_df.reset_index()).merge(features[['num_cells_s']].reset_index()))
plot_data['clone'] = plot_data['clone'] + ' n_s=' + plot_data['num_cells_s'].astype(str)

sns.displot(data=plot_data, x='rt', col='clone', kind='hist', col_wrap=4, common_norm=False, stat='proportion', bins=20)

```

```python

sns.displot(
    data=plot_data[plot_data['clone'].str.startswith('SA530')],
    x='rt', col='clone', kind='hist', col_wrap=3, common_norm=False, stat='proportion', bins=20)

```


# Analysis of global shifts


```python

Y = clone_rt.copy()

plot_data = Y.T.groupby(level=0).mean().melt(ignore_index=False, var_name='clone', value_name='mean_rt').reset_index()
plot_data = Y.T.melt(ignore_index=False, var_name='clone', value_name='rt').reset_index()
plot_data = plot_data[plot_data['chr'] != 'X']
plot_data = plot_data.merge(features_df.reset_index())

plt.figure(figsize=(4, 2))
sns.barplot(x='type', y='rt', data=plot_data, color='0.5', width=0.4)
sns.despine()

plt.figure(figsize=(4, 2))
sns.barplot(x='wgd', y='rt', data=plot_data, color='0.5', width=0.4)
sns.despine()

plt.figure(figsize=(4, 2))
sns.barplot(x='signature', y='rt', data=plot_data, color='0.5', width=0.4)
sns.despine()

```

```python

merged_index = encode.columns.intersection(clone_rt.columns)

combined = pd.concat([
    encode.loc[:, merged_index],
    clone_rt.loc[:, merged_index],
], axis=0)

combined

```

```python

Y = combined.copy()
Y.values[:] = scaler.fit_transform(Y.T).T

chromosomes = [str(a) for a in range(1, 23)] + ['X']
chromosomes = ['1', '2', 'X']
width = 20
width = 5

plot_data = Y.T.groupby(level=0).mean().melt(ignore_index=False, var_name='clone', value_name='mean_rt').reset_index()
plot_data = plot_data.merge(features_df.reset_index(), how='left')
plot_data.loc[plot_data['clone'].str.endswith('_rt'), 'type'] = plot_data.loc[plot_data['clone'].str.endswith('_rt'), 'clone']

plt.figure(figsize=(width, 4))
sns.barplot(x='chr', hue='type', y='mean_rt', order=chromosomes, data=plot_data, palette='tab10')
plt.legend(ncols=2, loc='upper left', bbox_to_anchor=(1., 1.), title='type')
sns.despine()

plot_data2 = plot_data[~plot_data['clone'].str.endswith('_rt')]

plt.figure(figsize=(width, 4))
sns.barplot(x='chr', hue='type', y='mean_rt', order=chromosomes, data=plot_data2, palette='tab10')
plt.legend(ncols=2, loc='upper left', bbox_to_anchor=(1., 1.), title='type')
sns.despine()

plot_data2 = plot_data[plot_data['clone'].str.endswith('_rt')]
plot_data2 = plot_data[
    # (plot_data['clone'] == 'mcf7_rt') |
    # (plot_data['clone'] == 'helas3_rt') |
    # (plot_data['clone'] == 'k562_rt') |
    (plot_data['clone'] == 'gm12878_rt') |
    (plot_data['clone'] == 'bj_rt') |
    (~plot_data['clone'].str.endswith('_rt'))]

plot_data2['type'].unique()

plt.figure(figsize=(width, 4))
sns.barplot(x='chr', hue='type', y='mean_rt', order=chromosomes, data=plot_data2, palette='tab10')
plt.legend(ncols=2, loc='upper left', bbox_to_anchor=(1., 1.), title='type')
sns.despine()

plot_data2 = plot_data[(plot_data['clone'] == 'mcf7_rt') | (~plot_data['clone'].str.endswith('_rt'))]

plt.figure(figsize=(width, 4))
sns.barplot(x='chr', hue='type', y='mean_rt', order=chromosomes, data=plot_data2, palette='tab10')
plt.legend(ncols=2, loc='upper left', bbox_to_anchor=(1., 1.), title='type')
sns.despine()

```

```python

Y = combined.copy()
Y.values[:] = scaler.fit_transform(Y.T).T

chromosomes = [str(a) for a in range(1, 23)] + ['X']
chromosomes = ['1', '2', 'X']
width = 20
width = 5

types_order = ['gm12878_rt', 'bj_rt', 'hTERT', 'OV2295', 'HGSOC', 'TNBC', 'T47D', 'GM18507']

plot_data = Y.T.groupby(level=0).mean().melt(ignore_index=False, var_name='clone', value_name='mean_rt').reset_index()
plot_data = plot_data.merge(features_df.reset_index(), how='left')
plot_data.loc[plot_data['clone'].str.endswith('_rt'), 'type'] = plot_data.loc[plot_data['clone'].str.endswith('_rt'), 'clone']

plot_data2 = plot_data[
    (plot_data['clone'] == 'gm12878_rt') |
    (plot_data['clone'] == 'bj_rt') |
    (~plot_data['clone'].str.endswith('_rt'))]

plt.figure(figsize=(width, 4))
sns.barplot(x='chr', hue='type', y='mean_rt', order=chromosomes, hue_order=types_order, data=plot_data2, palette='tab10')
plt.legend(ncols=2, loc='upper left', bbox_to_anchor=(0.05, 0.5), title='type', fontsize=8)
sns.despine()

```

```python

plot_data2 = plot_data[(plot_data['clone'] == 'mcf7_rt') | (~plot_data['clone'].str.endswith('_rt'))]

sns.heatmap(plot_data2.groupby(['chr', 'type'])['mean_rt'].mean().unstack().T)

```

```python

Y = clone_rt.copy()
Y.values[:] = scaler.fit_transform(Y.T).T

chromosomes = [str(a) for a in range(1, 23)] + ['X']

plot_data = Y.T.groupby(level=0).mean().melt(ignore_index=False, var_name='clone', value_name='mean_rt').reset_index()
plot_data = plot_data.merge(features_df.reset_index())

plt.figure(figsize=(20, 4))
sns.barplot(x='chr', hue='type', y='mean_rt', order=chromosomes, data=plot_data, palette='tab10')
plt.legend(ncols=2, loc='upper left', bbox_to_anchor=(1., 1.), title='signature')
sns.despine()

plt.figure(figsize=(20, 4))
sns.barplot(x='chr', hue='wgd', y='mean_rt', order=chromosomes, data=plot_data)
plt.legend(ncols=2, loc='upper left', bbox_to_anchor=(1., 1.))
sns.despine()

plt.figure(figsize=(20, 4))
sns.barplot(x='chr', hue='signature', y='mean_rt', order=chromosomes, data=plot_data)
plt.legend(ncols=2, loc='upper left', bbox_to_anchor=(1., 1.))
sns.despine()

```

```python

Y = clone_rt.copy()
Y = Y.set_index(features_type.loc[clone_rt.index]['type'], append=True)

chromosomes = [str(a) for a in range(1, 23)] + ['X']

plot_data = Y.T.var().rename('var_rt').reset_index(level=1).rename_axis('clone').sort_values('var_rt')
plot_data['num_cells_s'] = features['num_cells_s']

sns.scatterplot(x='var_rt', y='num_cells_s', data=plot_data)

plot_data.tail(30)

```

```python

import pandas as pd
import statsmodels.api as sm
import tqdm

X = selected_features
X = sm.add_constant(X)

Y = clone_rt.copy()
# Y.values[:] = scaler.fit_transform(Y.T).T

betas = pd.DataFrame(index=Y.columns, columns=X.columns)
betas_min = pd.DataFrame(index=Y.columns, columns=X.columns)
betas_max = pd.DataFrame(index=Y.columns, columns=X.columns)
pvalues = pd.DataFrame(index=Y.columns, columns=X.columns)

for idx in tqdm.tqdm(range(Y.shape[1])):
    y = Y.values[:, idx]
    model = sm.OLS(y, X)
    results = model.fit()
    # results = model.fit_regularized()
    betas.iloc[idx] = results.params
    betas_min.iloc[idx] = results.conf_int()[0]
    betas_max.iloc[idx] = results.conf_int()[1]
    pvalues.iloc[idx] = results.pvalues

```

```python

chromosomes = [
    '1',
    '2',
    '3',
    '4',
    '5',
    '6',
    '7',
    '8',
    '9',
    '10',
    '11',
    '12',
    '13',
    '14',
    '15',
    '16',
    '17',
    '18',
    '19',
    '20',
    '21',
    '22',
    'X',
]

plt.figure(figsize=(3, 3))
((betas.loc[chromosomes] ** 2).sum() ** 0.5).plot.bar()
sns.despine()

plt.figure(figsize=(3, 3))
(((betas - betas.mean(axis=0)).loc[chromosomes] ** 2).sum() ** 0.5).plot.bar()
sns.despine()

plt.figure(figsize=(3, 3))
(pvalues.loc[chromosomes] < 0.05).sum().plot.bar()
sns.despine()


```

```python


chromosomes = [str(a) for a in range(1, 23)] + ['X']

plot_data = (betas_max - betas_min).groupby(level=0).mean().melt(ignore_index=False, var_name='covariate', value_name='mean_beta').reset_index()

plt.figure(figsize=(20, 4))
sns.barplot(x='chr', hue='covariate', y='mean_beta', order=chromosomes, data=plot_data, palette='tab10')
plt.legend(ncols=2, loc='upper left', bbox_to_anchor=(1., 1.))

```

```python

chromosomes = [str(a) for a in range(1, 23)] + ['X']
# chromosomes = ['1', '2', 'X']

plot_data = betas.groupby(level=0).mean().melt(ignore_index=False, var_name='covariate', value_name='mean_beta').reset_index()
plot_data

plt.figure(figsize=(20, 4))
sns.barplot(x='chr', hue='covariate', y='mean_beta', order=chromosomes, data=plot_data, palette='tab10')
plt.legend(ncols=2, loc='upper left', bbox_to_anchor=(1., 1.))

```

```python

import scgenome
import anndata as ad

var = betas.reset_index()[['chr', 'start', 'end']]
var['chr'] = var['chr'].astype('category')
data = betas.reset_index(drop=True).T

adata = ad.AnnData(
    betas.reset_index(drop=True).T,
    var=var,
    layers={
        'betas_min': betas_min.reset_index(drop=True).T,
        'betas_max': betas_max.reset_index(drop=True).T,
        'pvalue': pvalues.reset_index(drop=True).T,
    })
adata.obs

```

```python

from scgenome.cnplot import *

def plot_cell_cn_profile(ax, cn_data, value_field_name, cn_field_name=None, max_cn=13, chromosome=None, s=5,
                         squashy=False, rawy=False, cmap=None, color=None):
    """ Plot copy number profile on a genome axis

    Args:
        ax: matplotlib axis
        cn_data: copy number table
        value_field_name: column in cn_data to use for the y axis value
    
    Kwargs:
        cn_field_name: state column to color scatter points
        max_cn: max copy number for y axis
        chromosome: single chromosome plot
        s: size of scatter points
        squashy: compress y axis
        rawy: raw data on y axis

    The cn_data table should have the following columns (in addition to value_field_name and
    optionally cn_field_name):
        - chr
        - start
        - end
    """
    chromosome_info = refgenome.info.chromosome_info[['chr', 'chromosome_start', 'chromosome_end']].copy()
    chromosome_info['chr'] = pd.Categorical(chromosome_info['chr'], categories=cn_data['chr'].cat.categories)

    plot_data = cn_data.merge(chromosome_info)
    plot_data = plot_data[plot_data['chr'].isin(refgenome.info.chromosomes)]
    plot_data['start'] = plot_data['start'] + plot_data['chromosome_start']
    plot_data['end'] = plot_data['end'] + plot_data['chromosome_start']

    squash_coeff = 0.15
    squash_f = lambda a: np.tanh(squash_coeff * a)
    if squashy:
        plot_data[value_field_name] = squash_f(plot_data[value_field_name])

    plot_data = plot_data.sort_values('start')
    
    if cn_field_name is not None:
        if cmap is not None:
            ax.plot(
                plot_data['start'], plot_data[value_field_name],
                c=plot_data[cn_field_name],
                cmap=cmap, lw=1, color=color,
            )
        else:
            ax.plot(
                plot_data['start'], plot_data[value_field_name],
                c=plot_data[cn_field_name],
                cmap=get_cn_cmap(plot_data[cn_field_name].astype(int).values), lw=1, color=color,
            )
    else:
        ax.plot(
            plot_data['start'], plot_data[value_field_name], lw=1, color=color,
        )

    if chromosome is not None:
        chromosome_length = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_length']
        chromosome_start = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_start']
        chromosome_end = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_end']
        xticks = np.arange(0, chromosome_length, 2e7)
        xticklabels = ['{0:d}M'.format(int(x / 1e6)) for x in xticks]
        xminorticks = np.arange(0, chromosome_length, 1e6)
        ax.set_xlabel(f'chromosome {chromosome}')
        ax.set_xticks(xticks + chromosome_start)
        ax.set_xticklabels(xticklabels)
        ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(xminorticks + chromosome_start))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlim((chromosome_start, chromosome_end))

    else:
        ax.set_xlim((-0.5, refgenome.info.chromosome_info['chromosome_end'].max()))
        ax.set_xlabel('chromosome')
        ax.set_xticks([0] + list(refgenome.info.chromosome_info['chromosome_end'].values))
        ax.set_xticklabels([])
        ax.xaxis.tick_bottom()
        ax.yaxis.tick_left()
        ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(refgenome.info.chromosome_info['chromosome_mid']))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(refgenome.info.plot_chromosomes))

    if squashy and not rawy:
        yticks = np.array([0, 2, 4, 7, 20])
        yticks_squashed = squash_f(yticks)
        ytick_labels = [str(a) for a in yticks]
        ax.set_yticks(yticks_squashed)
        ax.set_yticklabels(ytick_labels)
        ax.set_ylim((-0.01, 1.01))
        ax.spines['left'].set_bounds(0, 1)
    elif not rawy:
        ax.set_ylim((-0.05 * max_cn, max_cn))
        ax.set_yticks(range(0, int(max_cn) + 1))
        ax.spines['left'].set_bounds(0, max_cn)

    return chromosome_info


def plot_cn_profile(
        adata: AnnData,
        obs_id: str,
        value_layer_name=None,
        state_layer_name=None,
        ax=None,
        max_cn=13,
        chromosome=None,
        s=5,
        squashy=False,
        rawy=False,
        color=None,
    ):
    """Plot scatter points of copy number across the genome or a chromosome.

    Parameters
    ----------
    adata : AnnData
        copy number data
    obs_id : str
        observation to plot
    value_layer_name : str, optional
        layer with values for y axis, None for X, by default None
    state_layer_name : str, optional
        layer with states for colors, None for no color by state, by default None
    ax : [type], optional
        existing axis to plot into, by default None
    max_cn : int, optional
        max copy number for y axis, by default 13
    chromosome : [type], optional
        single chromosome plot, by default None
    s : int, optional
        size of scatter points, by default 5
    squashy : bool, optional
        compress y axis, by default False
    rawy : bool, optional
        raw data on y axis, by default False

    Examples
    -------

    .. plot::
        :context: close-figs

        import scgenome
        adata = scgenome.datasets.OV2295_HMMCopy_reduced()
        scgenome.pl.plot_cn_profile(adata, 'SA922-A90554B-R27-C43', value_layer_name='copy', state_layer_name='state')

    TODO: missing return
    """

    cn_data = adata.var.copy()

    if value_layer_name is not None:
        cn_data['value'] = np.array(adata[[obs_id], :].layers[value_layer_name][0])
    else:
        cn_data['value'] = np.array(adata[[obs_id], :].X[0])

    cn_field_name = None
    if state_layer_name is not None:
        cn_data['state'] = np.array(adata[[obs_id], :].layers[state_layer_name][0])
        cn_field_name = 'state'

    if ax is None:
        ax = plt.gca()

    cn_data = cn_data.dropna(subset=['value'])

    plot_cell_cn_profile(
        ax, cn_data, 'value', cn_field_name=cn_field_name, max_cn=max_cn,
        chromosome=chromosome, s=s, squashy=squashy, rawy=rawy, color=color)

    return ax


```

```python

smoothed = np.zeros(adata.shape)

for idx in range(smoothed.shape[0]):
    smoothed[idx, :] = pd.Series(adata.X[idx, :]).rolling(20).mean()
    
adata.layers['smoothed'] = smoothed

```

```python

chromosome = 'X'

plt.figure(figsize=(10, 2))
# scgenome.pl.plot_cn_profile(adata, 'const', rawy=True, s=1, chromosome=chromosome)
# plot_cn_profile(adata, 'const', value_layer_name=None, rawy=True, s=1, chromosome=chromosome, color='k')
plot_cn_profile(adata, 'const', value_layer_name='smoothed', rawy=True, s=1, chromosome=chromosome, color='k')
# scgenome.pl.plot_cn_profile(adata, 'const', value_layer_name='betas_min', rawy=True, s=1, chromosome=chromosome)
# scgenome.pl.plot_cn_profile(adata, 'const', value_layer_name='betas_max', rawy=True, s=1, chromosome=chromosome)
plt.title('const')

plt.figure(figsize=(10, 2))
plot_cn_profile(adata, 'type_hTERT', value_layer_name='smoothed', rawy=True, s=1, chromosome=chromosome)
# # scgenome.pl.plot_cn_profile(adata, 'type_hTERT', value_layer_name='betas_min', rawy=True, s=1, chromosome=chromosome)
# # scgenome.pl.plot_cn_profile(adata, 'type_hTERT', value_layer_name='betas_max', rawy=True, s=1, chromosome=chromosome)
plt.title('type_hTERT')

plt.figure(figsize=(10, 2))
scgenome.pl.plot_cn_profile(adata, 'type_hTERT', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'type_hTERT', value_layer_name='betas_min', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'type_hTERT', value_layer_name='betas_max', rawy=True, s=1, chromosome=chromosome)
plt.title('type_hTERT')

plt.figure(figsize=(10, 2))
scgenome.pl.plot_cn_profile(adata, 'type_HGSOC', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'type_HGSOC', value_layer_name='betas_min', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'type_HGSOC', value_layer_name='betas_max', rawy=True, s=1, chromosome=chromosome)
plt.title('type_HGSOC')

plt.figure(figsize=(10, 2))
scgenome.pl.plot_cn_profile(adata, 'type_TNBC', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'type_TNBC', value_layer_name='betas_min', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'type_TNBC', value_layer_name='betas_max', rawy=True, s=1, chromosome=chromosome)
plt.title('type_TNBC')

plt.figure(figsize=(10, 2))
scgenome.pl.plot_cn_profile(adata, 'signature_FBI', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'signature_FBI', value_layer_name='betas_min', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'signature_FBI', value_layer_name='betas_max', rawy=True, s=1, chromosome=chromosome)
plt.title('signature_FBI')

plt.figure(figsize=(10, 2))
scgenome.pl.plot_cn_profile(adata, 'signature_HRD', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'signature_HRD', value_layer_name='betas_min', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'signature_HRD', value_layer_name='betas_max', rawy=True, s=1, chromosome=chromosome)
plt.title('signature_HRD')

plt.figure(figsize=(10, 2))
scgenome.pl.plot_cn_profile(adata, 'wgd', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'wgd', value_layer_name='betas_min', rawy=True, s=1, chromosome=chromosome)
scgenome.pl.plot_cn_profile(adata, 'wgd', value_layer_name='betas_max', rawy=True, s=1, chromosome=chromosome)
plt.title('wgd')

plt.figure(figsize=(10, 2))
scgenome.pl.plot_cn_profile(adata, 'signature_FBI', value_layer_name='pvalue', rawy=True, s=1, chromosome=chromosome)
plt.title('wgd')

plt.figure(figsize=(10, 2))
scgenome.pl.plot_cn_profile(adata, 'signature_HRD', value_layer_name='pvalue', rawy=True, s=1, chromosome=chromosome)
plt.title('wgd')

```

```python

import pandas as pd

df = pd.read_csv(
    '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/spectrum/SPECTRUM-OV-081/s_phase_cells_with_scRT_filtered.csv.gz',
    low_memory=False, usecols=['cell_id', 'PERT_phase', 'clone_id'])
df

```

```python

```
