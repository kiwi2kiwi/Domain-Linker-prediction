import pickle

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

dictionary = open("40_cutoff_then_0_20.pkl","rb")
dict = pickle.load(dictionary)

fasta = open("cutoff_40_redundancy_0.20.fasta","w")

for key in dict:
    protein = dict[key]
    print(">"+key+"\n")
    fasta.write(">"+key+"\n")
    print(protein.sequence+"\n")
    fasta.write(protein.sequence+"\n")

fasta.close()



