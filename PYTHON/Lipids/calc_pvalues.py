import pandas as pd
import seaborn as sns
import numpy as np
from statsmodels.formula.api import ols
current_palette = sns.color_palette()
from scipy import stats
# %%
lipidData = pd.read_csv('../../../Data/Lipids/lipid_data.csv')
lipidData.index = lipidData.iloc[:,0]
lipidData.drop(["Unnamed: 0"],axis=1,inplace=True)
lipidData.info()

sampleGroup = lipidData.group.astype('int')
phenoType = lipidData.pheno.astype('int')
family = lipidData.fam.astype('int')

#%%
isoLabels = lipidData.columns[:-3]
sampleSel = np.where((sampleGroup == 0) | (sampleGroup == 1) | (sampleGroup==2))[0]
workData = lipidData.loc[sampleSel].copy()

# %% determine pvalues for ALD vs CTRL and CALD vs non-CALD, and calculate FC
pvalues = pd.DataFrame(index=isoLabels, columns=['pval', 'ptot', 'ols'])

for l in isoLabels:
    wd = workData.loc[:, [l, 'group', 'fam']].copy()
    wd.rename(columns={l: 'thevar'}, inplace=True)
    wd.loc[wd.group == 0, 'group'] = 1

    model1 = stats.ttest_ind(wd.loc[wd.group == 1, 'thevar'], wd.loc[wd.group == 2, 'thevar'], equal_var=False)
    pvalues.loc[l, 'ptot'] = model1[1]

    meangroups = pd.concat([np.log2(wd.thevar), wd.group], axis=1).groupby('group').mean().loc[[1, 2]]
    pvalues.loc[l, 'log2fc_overall'] = meangroups.thevar.diff()[2]


    wd = workData.loc[(workData.group == 1) | (workData.group == 0), [l, 'group', 'fam']].copy()
    wd.rename(columns={l: 'thevar'}, inplace=True)
    wd.group = wd.group.astype('int').astype('category')
    wd.fam = wd.fam.astype('int').astype('category')

    # the model between AMN and CALD .. pval == ols
    model = stats.ttest_rel(wd.loc[wd.group == 0, 'thevar'], wd.loc[wd.group == 1, 'thevar'])
    pvalues.loc[l, 'pval'] = model[1]

    model = ols("thevar ~ group + fam", wd).fit()
    pvalues.loc[l, 'ols'] = model.pvalues[1]

    meangroups = pd.concat([np.log2(wd.thevar), wd.group], axis=1).groupby('group').mean().loc[[0, 1]]
    pvalues.loc[l, 'log2fc'] = meangroups.thevar.diff()[1]

pvalues['log10p'] = -np.log10(pvalues.pval.astype('float'))  # log10 cALD vs AMN
pvalues['log10t'] = -np.log10(pvalues.ptot.astype('float'))  # log10 ALD vs CTRL uncorrected
# %% save to any file
# pvalues.to_csv('')