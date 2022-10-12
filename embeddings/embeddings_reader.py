import pickle
import h5py
import pandas as pd
import numpy as np
import time
from pathlib import Path

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


class fasta_protein:
   def __init__(self, id):
       self.id = id
       self.range = range
       self.seq = ""

   def add_seq(self, seq):
       self.seq=self.seq + seq

np.random.seed(29)
test1 = np.random.randint(0, 2, size=(3,5))
test2 = np.random.randint(0, 2, size=(2,5))
test=[]
test.append(test1)
test.append(test2)
test_df = pd.DataFrame(test)
print(test)

dictionary = open("../automated_protein_statistics/40_cutoff_then_0_20.pkl","rb")
prot_dict = pickle.load(dictionary)

filename = "protT5_xl_u50_v3_domains_cutoff_40_redundancy_0.20.h5"
f = h5py.File(filename,"r")
print("Keys: %s" % f.keys())

counter=0
total=len(list(f.keys()))

#erstellt eine liste die angibt ob der residue domain oder linker ist
def classification_output(key):
    isdomain = np.full((len(prot_dict[key].sequence)), True)
    for linker in prot_dict[key].linkers:
        for i in np.arange(linker.start, linker.end + 1):
            isdomain[i - 1] = False
    return isdomain

columns=np.char.mod('%d', np.insert(np.arange(0,1024),0,1)).tolist()
columns[0]="is_domain"
linker_dataframe=pd.DataFrame(columns=columns)
domain_dataframe =pd.DataFrame(columns=columns)
proteins=pd.DataFrame(columns=columns)

protein_array =[]
started=False
print("proteins : " + str(len(f)))
f.keys()
namelist = np.array(list(f.keys()))
#pickle.dump(namelist,open("protein_name_list.pkl","wb"))
WD = Path(__file__).resolve().parents[1]
np.save(file=WD / 'embeddings' / f'protein_name_list.npy', arr=namelist, allow_pickle=False)

for key in list(f.keys()):
    print("protein length: "+str(len(list(f[key]))))
    start=time.time()
    isdomain=classification_output(key)

    is_domain_array=np.array(isdomain)[:,None]
    to_append = np.concatenate((is_domain_array,list(f[key])),axis=1)
    #    print(to_append.shape)
#    if started==False:
#        protein_array=to_append
#        started=True
#    else:
#        print(protein_array.shape)
#        print(to_append.shape)
#    print(len(protein_array)," ",len(protein_array[0])," ",len(protein_array[0][0]))
#    print(len(to_append)," ",len(to_append[0])," ",len(to_append[0][0]))
    protein_array.append(to_append)
    #    print(protein_array.shape)
    counter+=1
    end=time.time()
    print("time taken: " + str(end-start))
    print(str(counter)+"/"+str(total))

#domain_dataframe= domain_dataframe.iloc[:len(linker_dataframe.index), ]

protein_embeddings_list = open("protein_emb_list.pkl","wb")
pickle.dump(protein_array,protein_embeddings_list)
protein_embeddings_list.close()

#protein_df = pd.DataFrame(protein_array, columns=columns)
#
#protein_file = open("protein_list.pkl","wb")
#pickle.dump(protein_df, protein_file)
#protein_file.close()

#linker_file = open("linker_list.pkl","wb")
#pickle.dump(linker_dataframe, linker_file)
#linker_file.close()
#
#
#domain_file = open("domain_list.pkl","wb")
#pickle.dump(domain_dataframe, domain_file)
#domain_file.close()



print("end")


