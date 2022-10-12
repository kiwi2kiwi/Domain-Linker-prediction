import pickle
from sklearn.manifold import TSNE
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
#import seaborn as sns
#from time import time

protein_file = open("protein_list.pkl","rb")
prots = pickle.load(protein_file)
sample_n = 17722#17722

linkers_long = prots.loc[prots['is_domain']==0]
linkers = linkers_long.sample(n=sample_n, random_state=1)
domains = prots.loc[prots['is_domain']==1]
balanced_domains = domains.sample(n=sample_n, random_state=1)
balanced_proteins = pd.concat([linkers,balanced_domains])
X = balanced_proteins.loc[:, balanced_proteins.columns != 'is_domain']
Y = balanced_proteins.loc[:,'is_domain']

perplex = round((sample_n*2)*0.01)
perplex = 50
learn_rate = 200
learn_rate = round((sample_n*2)/12)
components=2
m = TSNE(n_components=components,learning_rate=learn_rate,perplexity=perplex, random_state=1, init="pca").fit_transform(X)
print(m.shape)

plt.rcParams['font.size'] = 15
target_ids = [0,1]
colors = ["r","b"]
#zettel = np.concatenate((np.full((sample_n),"linker"),np.full((sample_n),"domain")))
colors = np.concatenate((np.full((sample_n),"r"),np.full((sample_n),"b")))
colors_linkers = np.full((sample_n),"r")
colors_domains = np.full((sample_n),"b")
alpha=0.2
size=4
tsne_domains=m[sample_n:,:]
tsne_linkers=m[:sample_n,:]

print(m[:200,0].shape)
fig, ax = plt.subplots()
if components == 1:
    plt=ax.scatter(tsne_linkers, tsne_linkers, c=colors_linkers, label="linker", alpha=alpha, s=size)
    ax.scatter(tsne_domains, tsne_domains, c=colors_domains, label="domain", alpha=alpha, s=size)
if components == 2:
    plt=ax.scatter(tsne_linkers[:,0], tsne_linkers[:,1], c=colors_linkers, label="linker", alpha=alpha, s=size)
    ax.scatter(tsne_domains[:,0], tsne_domains[:,1], c=colors_domains, label="domain", alpha=alpha, s=size)
if components == 3:
    plt=ax.axes(projection='3d')
    ax.scatter(tsne_linkers[:,0],tsne_linkers[:,1],tsne_linkers[:,2], c=colors_linkers, label="linker", alpha=alpha, s=size)
    ax.scatter(tsne_domains[:,0],tsne_domains[:,1],tsne_domains[:,2], c=colors_domains, label="domain", alpha=alpha, s=size)

fig.suptitle("samples: " + str(sample_n * 2) + "\nperplexity: " + str(perplex) + "\nlearning rate: " + str(learn_rate) + "\nn_components: " + str(components))
ax.legend()
fig.show()

print("ende")
#nparray = np.array([[4,5],[6,7]])
#bools = [1,0]
#bools = np.array(bools)
#first_table=np.concatenate((bools[:,None],nparray),axis=1)
#
#nparray2 = np.array([[7,8],[8,9]])
#bools2 = [1,0]
#bools2 = np.array(bools2)
#second_table=np.concatenate((bools[:,None],nparray2),axis=1)
##print("second_table\n"+str(second_table))
#
#third_table = np.concatenate((first_table,second_table),axis=0)
#test_df = pd.DataFrame(third_table,columns=[1,2,3])
#
#
#print("third_table\n"+str(third_table))
#
#
#print(np.concatenate((np.array([1]), np.array([1,2,3]))))
#
##dictionary = open("stripped_dictionary.pkl","rb")
##stripped_dictionary = pickle.load(dictionary)
##dict_frame = pd.DataFrame(stripped_dictionary)
##print(dict_frame.shape)
##del stripped_dictionary
##print(dict_frame.shape)
#array = np.array([1,2,3,4])
#array = pd.DataFrame(pd.Series(array))
#array = array.T
##bool=[True]
##bool_df=pd.DataFrame(data=[bool],columns=["is_domain"])
##test=pd.concat([bool_df,array],axis=1)
#a=[1,2,3,4]
#a=np.char.mod('%d', a)
##print(a)
##print(np.char.mod('%d', a))
#bool=np.array([True])
#columns=["is_domain"]
#columns=columns+a
#bool_df=pd.DataFrame(data=[bool,array.T],columns=columns)
##test=pd.concat([bool_df,array],axis=1) np.array([1,2,3,4])
##test=pd.DataFrame(data=[bool,array])
#print(bool_df)

#linker_file = open("linker_list.pkl","rb")
#linkers = pickle.load(linker_file)
#linker_frame = pd.DataFrame(linkers)
#print(linker_frame.shape)
#
#domain_file = open("domain_list.pkl","rb")
#domains = pickle.load(domain_file)
#domain_frame = pd.DataFrame(domains)
#print(domain_frame.shape)


