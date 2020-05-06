# %%
import pandas as pd

import numpy as np
import seaborn as sns
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests

from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl
# plt.style.use('seaborn-white')
import matplotlib.lines as mlines

# from mpl_toolkits.mplot3d import axes3d, Axes3D
current_palette = sns.color_palette()

# %% read the data

# df_pvalues = pd.read_csv('DEC2019_UPDATED_RESULTS\\ALD_DMRs_Annotated_Dec2019.csv')
pvalues = pd.read_csv('../../../Data/Methylation/ALD_Limma_Final_Dec2019_CHR.csv')
pvalues['log10p'] = -np.log10(pvalues.Nominal_P)
pvalues['absdb'] = abs(pvalues.Delta_Beta)
pvalues['grp'] = 0
pvalues.loc[(pvalues.log10p > (-np.log10(0.0005))), 'grp'] = 3
pvalues.loc[(pvalues.Delta_Beta < -0.05) & (pvalues.log10p > (-np.log10(0.0005))), 'grp'] = 1
pvalues.loc[(pvalues.Delta_Beta > 0.05) & (pvalues.log10p > (-np.log10(0.0005))), 'grp'] = 2

pvalues = pvalues.sort_values(['CHR', 'Coordinate'])

betas = pd.read_csv('../../../Data/Methylation/ALD_Deconvoluted_Betas_Dec2019.csv')
betas.rename(columns={'Unnamed: 0': 'cpg'}, inplace=True)
colnames = betas.cpg
betas = betas.drop('cpg', axis=1)
betas = betas.transpose().sort_index()
betas.columns = colnames

# add hardcoded phenotypes etc

betas['pheno'] = ['CALD', 'AMN', 'CALD', 'AMN', 'CALD', 'AMN', 'CALD', 'AMN', 'CALD', 'CALD', 'AMN', 'AMN']
betas.pheno = betas.pheno.astype('category')
betas['fam'] = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6, 3]
betas.fam = betas.fam.astype('category')

# read the data for all DMRs
all_dmrs = pd.read_csv('../../../Data/Methylation/ALD_DMRs_Annotated_Dec2019_0.10DB.csv')
all_dmrs['abs_db'] = abs(all_dmrs.DB)

# %% find DMR regions
# first sort all DMRs according to their FDR and abs(delta beta)
mysel = all_dmrs.sort_values(['FDR', 'abs_db'], ascending=[True, False]).copy()

# find subset of single loci that are in the list of significant DMRs
subset = pvalues.CpG.isin(mysel.cpg)
# merge both sets
tbl_2 = pvalues.merge(mysel, right_on='cpg', left_on='CpG')

# copy the pvalues data frame to use for the volcano plot
df_pval_cp = pvalues.copy()
df_pval_cp['in_dmr'] = 0
df_pval_cp.loc[subset, 'in_dmr'] = 1
df_pval_cp.loc[subset & (df_pval_cp.grp == 0), 'grp'] = 5

# find those single loci that are significant
tbl_cpg = pvalues.loc[((pvalues.log10p > (-np.log10(0.0005))) & (abs(pvalues.Delta_Beta) > 0.05))]
# merge with dmr to find those points to annotate later on
tbl_both = tbl_cpg.merge(all_dmrs, left_on='CpG', right_on='cpg')
# find points to annotate... top 15
tbl_annot = pd.concat([tbl_both, tbl_2.sort_values('abs_db', ascending=False)[0:15]]).drop_duplicates()

#%% define function to remove nested lists
from collections import Iterable


def single_list(list, ignore_types=(str)):
    for item in list:
        if isinstance(item, Iterable) and not isinstance(item, ignore_types):
            yield from single_list(item, ignore_types=(str))
        else:
            yield item
# %% plot DMRs from coordinates to cpgs to betas
plt.close('all')

uranges = tbl_annot.coord.unique()
for r in uranges:
    print(r)

    chrom = int(r.split(':')[0].strip('chr'))
    minr = int(r.split(':')[1].split('-')[0])
    maxr = int(r.split(':')[1].split('-')[1])

    _selection = (pvalues.CHR == chrom) & (pvalues.Coordinate >= minr) & (pvalues.Coordinate <= maxr)
    if sum(_selection)==0:
        continue

    rng = pvalues.loc[_selection].copy()

    if rng.shape[0] < 2:
        continue

    posx = rng.Coordinate.values
    _betas = betas[rng.CpG].copy()
    betasT = pd.DataFrame(_betas.unstack().copy())
    betasT.rename(columns={0: 'beta'}, inplace=True)

    nrrows = rng.shape[0]
    xpos = [j for j in single_list([[posx[i]] * 12 for i in range(nrrows)])]

    betasT['X'] = xpos
    betasT['pheno'] = ['CALD', 'non-CALD', 'CALD', 'non-CALD', 'CALD', 'non-CALD', 'CALD', 'non-CALD', 'CALD', 'CALD', 'non-CALD',
                       'non-CALD'] * nrrows
    betasT['fam'] = ['1', '1', '2', '2', '3', '4', '4', '5', '5', '6', '6', '3'] * nrrows
    betasT.pheno = betasT.pheno.astype('category')
    betasT.fam = betasT.fam.astype('category')

    betasT = betasT.sort_values(['X'])

    if rng.UCSC_RefGene_Name.any():
        txt = rng.UCSC_RefGene_Name.any()
        txt = txt.split(';')[0]
    else:
        txt = r

   f, ax = plt.subplots(figsize=(6, 6))

    sns.scatterplot(x='X', y='beta', style='fam', palette=[current_palette[f] for f in range(6)], data=betasT,
                    hue='fam', legend=False, ax=ax)
    sns.lineplot(x='X', y='beta', style='pheno', palette=[current_palette[f] for f in range(2)], hue='pheno',
                 data=betasT, legend='full', ax=ax, ci=95)

    handles, labels = ax.axes.get_legend_handles_labels()
    ax.set(xlabel='relative chromosomal coordinate (chr{0})'.format(chrom))
    plt.title(txt)
    
    sns.despine(offset=10, trim=True)

    plt.show()