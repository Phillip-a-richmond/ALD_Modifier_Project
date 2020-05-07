# %%
import pandas as pd
import numpy as np
from statsmodels.formula.api import ols
import warnings
# %% read raw data
df = pd.read_csv('../../../Data/Proteomics/Report_Precursor_Proteins.txt', sep='\t', na_values='Filtered')
# rename columns
for i in df.columns[3:]:
    df.rename(columns={i: i.split('.')[0].split('_')[1]}, inplace=True)

#%% drop id columns and transpose data
dfT = df.drop(['PG.Genes', 'PG.ProteinAccessions', 'PG.ProteinDescriptions'], axis=1).T
# rename columns to var0, var1 etc
for col in dfT.columns:
    dfT[col] = dfT[col].astype('float64')
    dfT.rename(columns={col: "var{0}".format(col)}, inplace=True)

# hard code order of phenotypes here
dfT['pheno'] = ['CALD', 'AMN', 'CALD', 'AMN', 'CALD', 'AMN', 'CALD', 'AMN', 'CALD', 'CALD', 'AMN', 'AMN']
dfT.pheno = dfT.pheno.astype('category')
# hard code family order
dfT['fam'] = [1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6, 3]
dfT.fam = dfT.fam.astype('category')

# save to file
dfT.to_csv('../../../Data/Proteomics/proteomic_data.csv')

#%% calculate the p-values and FC
warnings.simplefilter("ignore")
pvalues = pd.DataFrame(index=[v for v in dfT.columns if 'var' in v], columns=['pval'])
for col in dfT.columns:
    if "var" in col:
        workData = dfT.loc[:,[col,'pheno','fam']].copy()

        model = ols("{0} ~ pheno + fam".format(col), workData).fit()
        pvalues.loc[col, 'pval'] = model.pvalues['pheno[T.CALD]']
        meangroups = pd.concat([np.log2(dfT[col]), dfT.pheno], axis=1).groupby('pheno').mean().loc[['AMN', 'CALD']]
        pvalues.loc[col, 'log2fc'] = meangroups[col].diff()[1]

        pvalues.loc[col, 'geneid'] = df.loc[int(col.strip('var')), 'PG.Genes']
        pvalues.loc[col, 'protacc'] = df.loc[int(col.strip('var')), 'PG.ProteinAccessions']

# add some extra descriptors
pvalues['log10p'] = -np.log10(pvalues.pval.astype('float'))
pvalues['nrnan'] = 0
pvalues['nrnan'] = dfT.isnull().sum()[:-2]

pvalues['grp'] = 0
pvalues.loc[(pvalues.log10p > (-np.log10(0.05))), 'grp'] = 3
pvalues.loc[(pvalues.log2fc < -1) & (pvalues.log10p > (-np.log10(0.05))), 'grp'] = 1
pvalues.loc[(pvalues.log2fc > 1) & (pvalues.log10p > (-np.log10(0.05))), 'grp'] = 2

warnings.simplefilter("default")
# %%  save pvalues if necessary
pvalues.to_csv('../../../Data/Proteomics/pvalues.csv')

