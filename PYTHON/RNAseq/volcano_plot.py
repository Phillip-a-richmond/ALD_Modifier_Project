import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.lines as mlines

current_palette = sns.color_palette()

# %% load the data
data = pd.read_csv('../../../Data/RNA/tmm_norm_counts_5_2_2020.csv')
pvalues = pd.read_csv('../../../Data/RNA/pvalues_all_23_01_2019.csv')
# rename first column
pvalues = pvalues.rename(columns={pvalues.columns[0]: 'id'})

#%% select the data
allfam = pvalues.iloc[:, [0, 1, 2, -5]].copy()
allfam['logp'] = -np.log10(allfam.all_families_p)
allfam['grp'] = (allfam.logp > (-np.log10(0.05))) & (abs(allfam.all_families_logFC) > 1)
# %% create volcano plot
plt.close('all')

f = plt.figure(figsize=(8, 7))

p = sns.scatterplot(x='all_families_logFC', y='logp', hue='grp', data=allfam, legend=False)
ylim = p.axes.get_ylim()

p.axes.set_xlim([-4.5, 4.5])
xlim = p.axes.get_xlim()

plt.plot([-1, -1], ylim, linewidth=1, dashes=[6, 2], color='lightgrey')
plt.plot([1, 1], ylim, linewidth=1, dashes=[6, 2], color='lightgrey')
plt.plot(xlim, [-np.log10(0.05), -np.log10(0.05)], color='lightgrey', dashes=[6, 1])

sns.despine()
plt.rc('text', usetex=False)

# table for annotation
tbl_1 = allfam.loc[(allfam.logp > (-np.log10(0.05))) & (abs(allfam.all_families_logFC) > 1)]

for i in range(0, tbl_1.shape[0]):
    if not list(tbl_1.hgnc_symbol.isnull())[i]:
        txt = tbl_1.iloc[i].hgnc_symbol
    else:
        txt = tbl_1.iloc[i].id
# uncomment next line to display labels
        #p.text(tbl_1.iloc[i]['all_families_logFC'],tbl_1.iloc[i]['logp']*1.01,txt,fontsize=6)

plt.text(0.85, 0.02, '*: unadjusted', fontsize=8, transform=plt.gcf().transFigure)

params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)

plt.xlabel('log2(FoldChange)')
plt.ylabel('$-log10(pvalue)^*$')

plt.show()
