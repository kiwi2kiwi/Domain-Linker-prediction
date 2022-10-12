import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import umap

protein_file = open("protein_list.pkl","rb")
prots = pickle.load(protein_file)
sample_n = 13119#13119
alpha=0.2
size=4

linkers_long = prots.loc[prots['is_domain']==0]
linkers = linkers_long.sample(n=sample_n, random_state=1)
domains = prots.loc[prots['is_domain']==1]
balanced_domains = domains.sample(n=sample_n, random_state=1)
balanced_proteins = pd.concat([linkers,balanced_domains])
X = balanced_proteins.loc[:, balanced_proteins.columns != 'is_domain']
Y = balanced_proteins.loc[:,'is_domain']

reducer = umap.UMAP(random_state=1)
embedding = reducer.fit_transform(X)

tsne_domains=embedding[sample_n:,:]
tsne_linkers=embedding[:sample_n,:]

plt.scatter(tsne_domains[:, 0], tsne_domains[:, 1], c=["r"], label="domain", alpha=alpha, s=size)
plt.scatter(tsne_linkers[:, 0], tsne_linkers[:, 1], c=["b"], label="linker", alpha=alpha, s=size)
plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP projection of the linker/domain dataset', fontsize=24)
plt.suptitle("samples: " + str(sample_n * 2))
plt.show()

print("end")

