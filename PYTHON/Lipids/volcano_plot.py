
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import numpy as np
current_palette = sns.color_palette()

# %%
lipidData = pd.read_csv('../../../Data/Lipids/lipid_data.csv')
lipidData.info()

sampleGroup = lipidData.group.astype('int')
phenoType = lipidData.pheno.astype('int')
family = lipidData.fam.astype('int')

pvalues = pd.read_csv('../../../Data/Lipids/pvalues_lipids_tot.csv')
pvalues.index = pvalues.lipidid

# %% add descriptors to pvalue table

pvc = pvalues.copy()

for i in pvc.index:
    j = i.split('(')[1].split(')')[0].split(':')[0]
    if not j[0].isdigit():
        if j[0] == 't':
            tp = 3
        if j[0] == 'd':
            tp = 2
        lgt = abs(int(j[1:]))
    else:
        tp = 1
        lgt = abs(int(j))
    db = i.split('(')[1].split(')')[0].split(':')[1]
    tpl = i.split('(')[0]

    pvc.loc[i, 'lipidclass'] = tpl  # np.where(ulip==tpl)[0][0]
    pvc.loc[i, 'tp'] = tp
    pvc.loc[i,'clf']=1

    if lgt < 20:
        pvc.loc[i, 'len'] = '<20'
        pvc.loc[i, 'clf'] = 1
    if (lgt >= 20) & (lgt < 40):
        pvc.loc[i, 'len'] = '>=20&<40'
        pvc.loc[i, 'clf'] = 2
    if (lgt >= 40) & (lgt < 60):
        pvc.loc[i, 'len'] = '>=40&<60'
        pvc.loc[i, 'clf'] = 3
    if (lgt >= 60) & (lgt < 80):
        pvc.loc[i, 'len'] = '>=60&<80'
        pvc.loc[i, 'clf'] = 4
    if (lgt >= 80):
        pvc.loc[i, 'len'] = '>=80'
        pvc.loc[i, 'clf'] = 5

    pvc.loc[i, 'dbs'] = db

pvc.rename(columns={'len': 'chain length', 'lipidclass': 'lipid class'}, inplace=True)
# %% make volcano plot of CALD vs non-CALD

plt.close('all')
f = plt.figure(figsize=(8, 7))
p = sns.scatterplot(x='log2fc', y='log10p', size='chain length',alpha=0.7,sizes=(100,10) , hue='lipid class', data=pvc)
ylim = p.axes.get_ylim()
p.axes.set_xlim([-1.5,1.5])
xlim = p.axes.get_xlim()

plt.plot([-1, -1], ylim, linewidth=1, dashes=[6, 2], color='lightgrey')
plt.plot([1, 1], ylim, linewidth=1, dashes=[6, 2], color='lightgrey')
plt.plot([-1.3,1.3], [-np.log10(0.05), -np.log10(0.05)], color='lightgrey', dashes=[6, 1])

sns.despine(trim=False)

plt.rc('text', usetex=False)
box = p.get_position()
p.set_position([box.x0, box.y0, box.width * 0.85, box.height])  # resize position
leg = p.legend()

# Put a legend to the right side
p.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1, frameon=False)

# add some labels
tbl = pvalues.loc[(pvalues.log2fc>=0) & (pvalues.log10p>(-np.log10(0.025)))]
tbl_1 = pvc.loc[(pvc.log10p > (-np.log10(0.05)))]
for i in range(0, tbl_1.shape[0]):
    txt = tbl_1.index[i]
    txt = txt.split('_')[0]
    txt = txt.split('[')[0]
    p.text(tbl_1.iloc[i]['log2fc'], tbl_1.iloc[i]['log10p'] * 1.01, txt, fontsize=6)

pvalues.sort_values('log2fc', ascending=True)[0:10]
pvalues.loc[pvalues.grp == 2].sort_values('log2fc', ascending=False)[0:10]

params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)

plt.xlabel('log2(FoldChange)')
plt.ylabel('$-log10(pvalue)^*$')

plt.text(0.85, 0.02, '*: unadjusted', fontsize=8, transform=plt.gcf().transFigure)

plt.tight_layout()
plt.show()

#%% create ALD vs CTRL volcano plot

plt.close('all')
f = plt.figure(figsize=(8, 7))

p = sns.scatterplot(x=-pvc['log2fc_overall'], y='log10t', alpha=0.7,sizes=(100,10),size='chain length',hue='lipid class', data=pvc)
ylim = p.axes.get_ylim()
xlim = p.axes.get_xlim()
plt.plot([-1, -1], ylim, linewidth=1, dashes=[6, 2], color='lightgrey')
plt.plot([1, 1], ylim, linewidth=1, dashes=[6, 2], color='lightgrey')
plt.plot(xlim, [-np.log10(0.05 / 1160), -np.log10(0.05 / 1160)], color='lightgrey', dashes=[6, 1])

sns.despine(trim=False)

box = p.get_position()
p.set_position([box.x0, box.y0, box.width * 0.85, box.height])  # resize position
leg = p.legend()
# Put a legend to the right side
p.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1, frameon=False)

# add some text labels

tbl_2 = pvalues.loc[(pvalues.log10t > -np.log10(0.05 / 1160))]
for i in range(0, tbl_2.shape[0]):
    txt = tbl_2.index[i]
    txt = txt.split('_')[0]
    txt = txt.split('[')[0]
    p.text(-tbl_2.iloc[i]['log2fc_overall'], tbl_2.iloc[i]['log10t'] * 1.01, txt, fontsize=6)

plt.xlabel('log2(FoldChange)')
plt.ylabel('$-log10(pvalue)^*$')

plt.text(0.85, 0.02, '*: unadjusted', fontsize=8, transform=plt.gcf().transFigure)

plt.tight_layout()
plt.show()