import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

current_palette = sns.color_palette()

# %% load the data
lipidData = pd.read_csv('data\lipid_data.csv')
lipidData.info()


# %%
sampleGroup = lipidData.group.astype('int')
phenoType = lipidData.pheno.astype('int')
family = lipidData.fam.astype('int')

# %% perform pca with control samples

sampleSel = (sampleGroup == 0) | (sampleGroup == 1) | (sampleGroup == 2)
# define pca data as everything except last 3 columns
pcaData = lipidData.loc[sampleSel,lipidData.columns[1:-3]]

# mean center data and autoscale the data
X = StandardScaler(with_std=True).fit_transform(pcaData)

# calculate first 2 components
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(X)
dfPCA = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'], index=pcaData.index)

# add sample information
dfPCAplot = pd.concat([dfPCA, pd.concat([sampleGroup[sampleSel], phenoType[sampleSel],family[sampleSel]],axis=1)], axis=1)

# %% create the PCA plot

plt.close("all")
fig = plt.figure(figsize=[8, 8], dpi=300)
p = sns.scatterplot(x=-dfPCAplot['PC1'], y=dfPCAplot['PC2'], s=100, hue=dfPCAplot.group, style=dfPCAplot.group,
                    palette=current_palette[0:3])

p.set_xlabel('PC1 ({0:2.1f}%)'.format(pca.explained_variance_ratio_[0] * 100))
p.set_ylabel('PC2 ({0:2.1f}%)'.format(pca.explained_variance_ratio_[1] * 100))

l = p.axes.get_legend()
new_labels = ['Phenotype', 'non-CALD', 'CALD', 'Control', 'Batch Control']
for t, l in zip(l.texts, new_labels): t.set_text(l)

sns.despine(offset=1, trim=False)

ylim = p.axes.get_ylim()
xlim = p.axes.get_xlim()
plt.plot([0, 0], ylim, linewidth=1, dashes=[6, 2], color='lightgrey')
plt.plot(xlim, [0, 0], linewidth=1, dashes=[6, 2], color='lightgrey')

plt.show()
