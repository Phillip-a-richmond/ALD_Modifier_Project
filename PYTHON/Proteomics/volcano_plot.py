
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

current_palette = sns.color_palette()

#%%
pvalues = pd.read_csv('../../../Data/Proteomics/pvalues.csv')

#%%
sns.set_style("whitegrid", {'axes.grid' : False})

plt.close('all')

f = plt.figure(figsize=(8, 7))
# make the plot
p = sns.scatterplot(x='log2fc', y='log10p', hue='grp', data=pvalues, palette=current_palette[0:3], legend=False)

# format the plot
ylim = p.axes.get_ylim()
p.axes.set_xlim([-3,3])
xlim = p.axes.get_xlim()
plt.plot([-1, -1], ylim, linewidth=1, dashes=[6, 2], color='lightgrey')
plt.plot([1, 1], ylim, linewidth=1, dashes=[6, 2], color='lightgrey')
plt.plot(xlim, [-np.log10(0.05), -np.log10(0.05)], color='lightgrey', dashes=[6, 1])

sns.despine()
plt.rc('text', usetex=False)

box = p.get_position()
p.set_position([box.x0, box.y0, box.width * 0.85, box.height])  # resize position

# add text to 'significant' labels
tbl_1 = pvalues.loc[(pvalues.log10p > (-np.log10(0.05)))]
for i in range(0, tbl_1.shape[0]):
    txt = tbl_1.geneid.iloc[i]
    if str(txt) != 'nan':
        txt = txt.split(';')[0]
        txt = txt.split('[')[0]
    else:
        txt = tbl_1.protacc.iloc[i]
        txt = txt.split(';')[0]
        txt = txt.split('[')[0]

    p.text(tbl_1.iloc[i]['log2fc'], tbl_1.iloc[i]['log10p'] * 1.01, txt, fontsize=6)

pvalues.sort_values('log2fc', ascending=True)[0:10]
pvalues.loc[pvalues.grp == 2].sort_values('log2fc', ascending=False)[0:10]

params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)

plt.xlabel('log2(FoldChange)')
plt.ylabel('$-log10(pvalue)^*$')

plt.text(0.85, 0.02, '*: unadjusted', fontsize=8, transform=plt.gcf().transFigure)

plt.show()