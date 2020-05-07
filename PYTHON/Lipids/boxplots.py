import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np
current_palette = sns.color_palette()


#%% remove rt info from lipid ids
def reformat_lipids(lipids):
    _lipids = {}
    __lipids = []
    for l in lipids:
        _l = (l.split(' ')[0]).split('_')[0]
        if(_l in _lipids.keys()):
            _lipids[_l]=_lipids[_l]+1
            __lipids.append(_l+"_{0}".format(_lipids[_l]+1))
        else:
            _lipids[_l]=1
            __lipids.append(_l)
    return __lipids

# %%
lipidData = pd.read_csv('../../../Data/Lipids/lipid_data.csv')
lipidData.columns = reformat_lipids(lipidData.columns)
lipidData.info()

sampleGroup = lipidData.group.astype('int')
phenoType = lipidData.pheno.astype('int')
family = lipidData.fam.astype('int')

pvalues = pd.read_csv('../../../Data/Lipids/pvalues_lipids_tot.csv')
pvalues.index =  reformat_lipids(pvalues.lipid)

# %%
sampleSel = (sampleGroup == 0) | (sampleGroup == 1) | (sampleGroup == 2)
features_ov = pvalues[(pvalues.ptot * 1160) < 0.05].index.values
features_ald = pvalues[pvalues.pval < 0.05].index.values

features_tot = np.unique(list(features_ov)+list(features_ald))
features_both = np.intersect1d(features_ov,features_ald)

ordered_data = pvalues.loc[features_tot].sort_values(['ptot','pval']).copy()
bpData = pd.concat([lipidData.loc[sampleSel, ordered_data.index], lipidData.loc[sampleSel,lipidData.columns[-3:]]], axis=1)

pval_ov = pvalues.loc[features_tot].ptot
pval_ald = pvalues.loc[features_tot].pval

# %% 38 features in total create grid of 6x7

nrcols = 6
nrrows = 7

f, axes = plt.subplots(nrows=nrrows, ncols=nrcols,figsize=(nrcols*2,nrrows*2))

for i in range(0, nrrows):
    for j in range(0, nrcols):
        if i * nrcols + j < len(features_tot):
            axs = axes[i, j]
            outlierprops = dict(markerfacecolor='0.75', markersize=2, linestyle='none')

            dd = sns.boxplot(x="pheno", y=bpData.iloc[:, i * nrcols + j], data=bpData, ax=axs, order=[2, 0, 1],linewidth=0.5,flierprops=outlierprops)
            dd.xaxis.set_ticklabels(['CTRL', 'non-CALD', 'CALD'], fontsize=8,rotation=0)
            lbs = axs.yaxis.get_ticklabels()

            ylbl = dd.yaxis.get_label().get_text()

            pvalov = pval_ov[ylbl]
            pvalald = pval_ald[ylbl]

            lbtxt = ''
            if ylbl in features_ov:
                lbtxt = lbtxt + '^:{0:1.3} '.format(pvalov)
            if ylbl in features_ald:
                lbtxt = lbtxt + '*:{0:1.3}'.format(pvalald)

            axs.set_title('{0}\n{1}'.format(ylbl,lbtxt), fontdict={'fontsize': 8, 'fontweight': 'light'})

            axs.set_ylabel(axs.get_ylabel(), fontsize=10)
            axs.yaxis.set_label_text('')

plt.setp(axes, xlabel='')

# remove remaining axes from plot
for i in range(2,6):
     f.delaxes(axes[6, i])
plt.tight_layout()
plt.text(0.5, 0.02, '^: significant ALD/CTRL adjusted, *: significant CALD/non-CALD unadjusted', fontsize=8, transform=plt.gcf().transFigure)
plt.show()
