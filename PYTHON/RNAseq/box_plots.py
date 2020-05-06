# %%
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

# %%
# copy the results
famResults = pvalues.columns[1:-5]
# copy the annotation data
annotData = pvalues.columns[-5:]


# %% copy data for with all families included
allfam = pvalues.iloc[:, [0, 1, 2, -5]].copy()
allfam['logp'] = -np.log10(allfam.all_families_p)
allfam['grp'] = (allfam.logp > (-np.log10(0.05))) & (abs(allfam.all_families_logFC) > 1)



# %%
dnorm = data.copy()
# rescale to 1e6 reads
dnorm.iloc[:, 1:] = dnorm.iloc[:, 1:] / dnorm.iloc[:, 1:].sum() * 1e6
dnorm = dnorm.rename(columns={dnorm.columns[0]: 'id'})
dm = dnorm.merge(allfam,left_on='id',right_on='id')

#%% create table with p <0.05 AND logFC>1 to use for boxplots
tbl_1 = dm.loc[(dm.logp > (-np.log10(0.05))) & (abs(dm.all_families_logFC) > 1)].sort_values('logp',ascending=False).copy()
tbl_1.index = tbl_1.id

genes = tbl_1.hgnc_symbol.copy()

tbl_1.drop(columns=['id'],inplace=True)
bpdata = tbl_1.iloc[:,:-5].transpose().sort_index().copy()
bpdata['pheno'] = [1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0]
bpdata['fam'] = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6, 3]
#
bpdata.pheno = bpdata.pheno.astype('category')
bpdata.fam = bpdata.fam.astype('category')

# %% create the boxplots

totnrplots = bpdata.shape[1] - 2

dims = [int(np.ceil(np.sqrt(totnrplots))), int(np.ceil(np.sqrt(totnrplots)))]
if dims[0] * (dims[1] - 1) >= totnrplots:
    dims[0] = dims[0] - 1

entrezs = bpdata.columns[:-2]
f, axes = plt.subplots(nrows=dims[0], ncols=dims[1], figsize=(dims[0]*2.5, dims[1]*2))

_id = 0

for i in range(0, dims[0]):
    for j in range(0, dims[1]):
        if i * dims[1] + j <= (totnrplots - 1):
            axs = axes[i, j]

            sns.boxplot(x="pheno", y=bpdata.iloc[:, _id], data=bpdata, ax=axs, order=[0, 1],color='white')
            plt.setp(axs.artists, edgecolor='k', facecolor='w')
            plt.setp(axs.lines, color='k')

            pp = sns.pointplot(x="pheno", y=bpdata.iloc[:, _id], hue='fam', data=bpdata,
                               jitter=0, dodge=True, linewidth=0.05, width=0.05, ax=axs, alpha=0.1,scale=0.75)

            ylabel = axs.yaxis.get_label().get_text()
            ylabel = ((ylabel.split('_')[0]).split('#')[0]).strip()

            if str(genes[_id]) != 'nan':
                ylabel = genes[_id]
            else:
                ylabel = genes.index[_id]

            axs.set_ylabel(ylabel)

            _id += 1

            lgd = axs.get_legend()
            lgd.remove()


# define and plot legends
alines = []
labels = []
for f in range(1, 7):
    aline = mlines.Line2D([], [], color=current_palette[f - 1], label='fam {0}'.format(f))
    alines.append(aline)
    labels.append(aline.get_label())

plt.figlegend(handles=alines, labels=labels, loc='best', bbox_to_anchor=(0.9,0.12))
xt = plt.setp(axes, xticks=[0, 1], xticklabels=['non-CALD', 'CALD'])
plt.setp(axes, xlabel='')

# remove remaining axes still shown
for i in range(3,8):
    plt.delaxes(axes[6, i])

plt.tight_layout()

plt.show()