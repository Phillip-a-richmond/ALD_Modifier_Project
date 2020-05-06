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

# %% create volcono plot with significant DMRs
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

# %% create the volcano plot with DMRs highlighted
plt.close('all')
f = plt.figure(figsize=(5, 5))

kwargs  =   {'rasterized':True}

# the CpGs with low beta and high pvalue are not drawn individually
sns.scatterplot(x='Delta_Beta', y='log10p', style='in_dmr', alpha=0.01, data=df_pval_cp.loc[df_pval_cp.grp == 0, :],palette=current_palette[0],legend=False,**kwargs)
sns.scatterplot(x='Delta_Beta', y='log10p', style='in_dmr', alpha=1, data=df_pval_cp.loc[df_pval_cp.grp == 1, :],palette=current_palette[4],legend=False)
sns.scatterplot(x='Delta_Beta', y='log10p', style='in_dmr', alpha=1, data=df_pval_cp.loc[df_pval_cp.grp == 2, :],palette=current_palette[2],legend=False)
sns.scatterplot(x='Delta_Beta', y='log10p', size=5, marker='x', alpha=1,data=df_pval_cp.loc[df_pval_cp.grp == 5, :], palette=current_palette[1],legend=False)
p = sns.scatterplot(x='Delta_Beta', y='log10p', style='in_dmr', alpha=0.01, data=df_pval_cp.loc[df_pval_cp.grp == 3, :],palette=current_palette[0],legend=False)

mrk_y = plt.scatter([],[],color='red',marker='x',label='Yes')
mrk_n = plt.scatter([],[],color='black',marker='o',label='No')
plt.legend(handles=[mrk_n,mrk_y],frameon=False,title='In DMR')

ylim = p.axes.get_ylim()
p.axes.set_xlim([-0.2, 0.2])
xlim = p.axes.get_xlim()

plt.plot([-.05, -.05], ylim, linewidth=1, dashes=[6, 2], color='grey')
plt.plot([.05, .05], ylim, linewidth=1, dashes=[6, 2], color='grey')
plt.plot(xlim, [-np.log10(0.0005), -np.log10(0.0005)], color='grey', dashes=[6, 1])

sns.despine()
plt.rc('text', usetex=False)
plt.text(0.85, 0.02, '*: unadjusted', fontsize=8, transform=plt.gcf().transFigure)

params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)

plt.xlabel('ß change', fontsize=10)
plt.ylabel('$-log10(pvalue)^*$',fontsize=10)

plt.show()