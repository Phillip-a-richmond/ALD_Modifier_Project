

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines

current_palette = sns.color_palette()

# %%
proteomicData = pd.read_csv('../../../Data/Proteomics/proteomic_data.csv',index_col=0)
pvalues = pd.read_csv('../../../Data/Proteomics/pvalues.csv', index_col=0)

# %% select those proteomes with low p-value (<0.05)
sampleSel = pvalues.loc[(pvalues.log10p>(-np.log10(0.05)))]
bpdata = proteomicData.loc[:,sampleSel.index.to_list()+['pheno','fam']].copy()
bpdata = bpdata.sort_values(['fam','pheno'])

#%%
# define grid of 4x5
totnrplots = bpdata.shape[1]
dims = [4,5]
f, axes = plt.subplots(nrows=dims[0], ncols=dims[1], figsize=(dims[0]*3.5, dims[1]*2))

for i in range(0, dims[0]):
    for j in range(0, dims[1]):
        if i * dims[1] + j <= (totnrplots - 1):
            axs = axes[i, j]

            yd = bpdata.iloc[:, i * dims[0] + j].copy()

            sns.boxplot(x="pheno", y=yd, data=bpdata, ax=axs, order=['AMN', 'CALD'])
            plt.setp(axs.artists, edgecolor='k', facecolor='w')
            plt.setp(axs.lines, color='k')

            # if there are missings then a line cannot be drawn
            for _r in range(6):
                _ids = [_r*2,_r*2+1]
                if sum(np.isnan(yd[_ids]))==0:
                    sns.lineplot(x=bpdata.pheno[_ids],y=yd[_ids],ax=axs,marker='o')

            axs.xaxis.set_major_formatter(ticker.EngFormatter())

            yorig = axs.yaxis.get_label().get_text()
            ylabel = pvalues.loc[yorig].geneid

            if str(ylabel) != 'nan':
                ylabel = ylabel.split(';')[0]
                ylabel = ylabel.split('[')[0]
            else:
                ylabel = pvalues.loc[yorig].protacc
                ylabel = ylabel.split(';')[0]
                ylabel = ylabel.split('[')[0]

            axs.set_ylabel(ylabel)

alines = []
labels = []
# define legend entries
for f in range(1, 7):
    aline = mlines.Line2D([], [], color=current_palette[f - 1], label='fam {0}'.format(f))
    alines.append(aline)
    labels.append(aline.get_label())

plt.figlegend(handles=alines, labels=labels, loc='lower right')
xt = plt.setp(axes, xticks=[0, 1], xticklabels=['non-CALD', 'CALD'])
plt.setp(axes, xlabel='')

# remove remaining axes from 4th row
for i in range(1,5):
    plt.delaxes(axes[3,i])

# display plot
plt.tight_layout()
plt.show()
