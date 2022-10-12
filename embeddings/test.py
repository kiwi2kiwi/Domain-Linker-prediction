import time
import h5py
import numpy as np
import pickle
import pandas as pd

class Protein:
   def __init__(self, name, domains: list, regions: list):
       self.name = name
       self.domains = [1]
       self.regions = []
       self.domains[0]=domains
       # diese regionenliste enthaelt nur domains
       for region in regions:
           self.regions.append(region)
       self.linkers = []
       self.domain_counter = 0
       self.longlinker = False
       self.longerthan40linker = False
       self.shortlinker = False
       self.domainresidues = 0
       self.linkerresidues = 0
       self.sequence = ""
       self.missing_in_dssp = False
       self.incutoff = False

   def add_domain(self,domain):
       self.domains.append(domain)

   def add_region(self,region):
       self.regions.append(region)

class Domain:
   def __init__(self, name, regions: list):
       self.name = name
       self.regions = []
       self.regions.append(regions)

   def add_region(self,region):
       self.regions.append(region)

class Region:
   def __init__(self, start, end):
       self.start = start
       self.end = end
       self.length = end - start + 1
       self.sequence = ""

a = np.array([[[1,2]],[[3,4]],[[5,6]],[[7,8]]])
b = np.array([[1,2],[3,4],[5,6],[7,8]])
print(a)
print(b[0,1])
c = [[[]]]
#print(a[0,0,0])
#print(a[0,0])
#print(a[0])
#print(a)
c[0][0] = [0,1]
#print(np.append(c,[[[9,10,11]]],axis = 0))



protein_file = open("protein_list.pkl","rb")
vectors = pickle.load(protein_file)

dictionary = open("../automated_protein_statistics/40_cutoff_then_0_20.pkl","rb")
prot_dict = pickle.load(dictionary)

filename = "protT5_xl_u50_v3_domains_cutoff_40_redundancy_0.20.h5"
f = h5py.File(filename,"r")
total=len(list(f.keys()))

vectors=vectors.to_numpy()
passed = 0
counter = 0
protein_clustered = []#np.array([[[]]])
started = False
print("proteins : " + str(len(f)))
for key in list(f.keys()):
    print("protein length: " + str(len(list(f[key]))))
    start = time.time()
    length = len(prot_dict[key].sequence)
    if not started:
        protein_clustered=[vectors[passed:passed+length]]
        started = True
    else:
        protein_clustered.append(vectors[passed:passed+length])
    #    protein_clustered = np.concatenate((protein_clustered[counter],vectors[passed,passed+length]),axis=0)
    passed = passed+length
    counter+=1
    end=time.time()
    print("time taken: " + str(end-start))
    print(str(counter)+"/"+str(total))
#    print(protein_clustered.shape)

protein_names = open("protein_names.pkl","wb")
pickle.dump(list(f.keys()),protein_names)
protein_names.close()

per_protein = open("per_protein_vectors.pkl","wb")
pickle.dump(protein_clustered,per_protein)
per_protein.close()
