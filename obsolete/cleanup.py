import csv
import optparse
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
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

from Bio.PDB import PDBList
'''Selecting structures from PDB'''
pdbl = PDBList()

#this programm only plots and doesn't write


dictionary = open("C:\\Users\\yanni\\Desktop\\Bachelor\\clean_Protein_dictionary.pkl","rb")
prot_dict = pickle.load(dictionary)
pdb_ids = prot_dict.keys()

plt.rcParams['font.size'] = 15

#TODO wie lange ist domain und linker im schnitt
def protein_linker(domain_lengths,linker_lengths,protein_lengths):
    domain_lengths_sorted = sorted(domain_lengths)
    linker_lengths_sorted = sorted(linker_lengths)
    domain_lengths_capped = domain_lengths_sorted[
                            :len(domain_lengths_sorted) - (round(0.05 * len(domain_lengths_sorted)))]
    linker_lengths_capped = linker_lengths_sorted[
                             :len(linker_lengths_sorted) - (round(0.05 * len(linker_lengths_sorted)))]
    f, (ax) = plt.subplots(1)
    # domain lengths
    pro = sns.kdeplot(data=protein_lengths, label="Proteins", ax=ax, color="red")
    # domain lengths
    dom = sns.kdeplot(data=domain_lengths, label="Domains", ax=ax, color="blue")
    # dom = sns.kdeplot(data = pd.DataFrame({"domain_lengths_capped": domain_lengths_capped},columns=["domain_lengths_capped"]), x="domain_lengths_capped")
    # linker lengths
    lin = sns.kdeplot(data=linker_lengths, label="Linkers", ax=ax, color="orange")
    # pro = sns.kdeplot(data = pd.DataFrame({"protein_lengths_capped": protein_lengths_capped},columns=["protein_lengths_capped"]), x="protein_lengths_capped")
    pro.set(xticks=[10, 25, 50, 100, 250, 500])
    # lin.set(xticks=[10, 25, 50, 100, 250, 500])
    ax.legend()
    plt.xlim([0, 500])
    mean_p = np.percentile(protein_lengths,[50])
    mean_d = np.percentile(domain_lengths,[50])
    mean_l = np.percentile(linker_lengths,[50])
    textstr = "median protein length: " + str(round(mean_p[0])) + "\n" + \
              "median domain length: " + str(round(mean_d[0])) + "\n" + \
              "median linker length: " + str(round(mean_l[0]))
    props = dict(boxstyle='round', facecolor='wheat', alpha=1)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=15, verticalalignment='top', bbox=props)
    plt.grid()
    plt.title("Comparison of Domain, Protein and Linker lengths for the clean data")
    plt.show()

def com_cumul(domain_lengths, linker_lengths):
    f, (ax) = plt.subplots(1)
    dom_cumul = sns.ecdfplot(data=domain_lengths, label="Domains", ax=ax)
    # dom = sns.kdeplot(data = pd.DataFrame({"domain_lengths_capped": domain_lengths_capped},columns=["domain_lengths_capped"]), x="domain_lengths_capped")
    # protein lengths
    lin_cumul = sns.ecdfplot(data=linker_lengths, label="Linkers", ax=ax, color="orange", complementary=True)
    # pro = sns.kdeplot(data = pd.DataFrame({"protein_lengths_capped": protein_lengths_capped},columns=["protein_lengths_capped"]), x="protein_lengths_capped")
    dom_cumul.set(xticks=[10, 25, 50, 100, 250, 500])
    ax.legend()
    ps_d = np.percentile(domain_lengths, [5, 95, 50])
    ps_l = np.percentile(linker_lengths, [5, 95, 50])
    print("5%-tile cut-off for domains: ", ps_d[0])
    print("95%-tile cut-off for domains: ", ps_d[1])
    len_ten = len([1 for i in domain_lengths if i < 10])
    len_fif = len([1 for i in domain_lengths if i < 50])
    len_four = len([1 for i in domain_lengths if i < 40])
    print(str(len_ten), " of domains <10: ", (len_ten / len(domain_lengths)))
    print(str(len_four), " of domains <40: ", (len_four / len(domain_lengths)))
    print("mean domain length: ", round(np.mean(domain_lengths)))
    print("Standard deviation of domain length: ", np.std(domain_lengths))
    print("5%-tile cut-off for linkers: ", ps_l[0])
    print("95%-tile cut-off for linkers: ", ps_l[1])
    print("mean linker length: ", round(np.mean(linker_lengths)))
    print("Standard deviation of linker length: ", np.std(linker_lengths))
    _, ymax = ax.get_ybound()
    ax.axvline(ps_d[0], label="5%", color="lightblue", linestyle="dashed", linewidth=2)
    ax.axvline(ps_d[1], label="95%", color="blue", linestyle="dashed", linewidth=2)
    ax.axvline(ps_d[2], label="median", color="darkblue", linestyle="dashed", linewidth=2)
    ax.axvline(ps_l[0], label="5%", color="pink", linestyle="dashed", linewidth=2)
    ax.axvline(ps_l[1], label="95%", color="red", linestyle="dashed", linewidth=2)
    ax.axvline(ps_l[2], label="median", color="darkred", linestyle="dashed", linewidth=2)
    plt.legend()
    plt.xlim(0, 500)
    plt.title("comparison of cumulative distributions of domain lengths and linker lengths", fontsize=15)
    textstr = "5%-tile cut-off for domain length: " + str(round(ps_d[0])) + "\n" + \
              "95%-tile cut-off for domain length: " + str(round(ps_d[1])) + "\n" + \
              "5%-tile cut-off for linker length: " + str(round(ps_l[0])) + "\n" + \
              "95%-tile cut-off for linker length: " + str(round(ps_l[1])) + "\n" + \
              str(len_ten) + " of domains <10: " + str(round((100 * (len_ten / len(domain_lengths))), 2)) + "%" + "\n" + \
              str(len_four) + " of domains <40: " + str(round((100 * (len_four / len(domain_lengths))), 2)) + "%"
    props = dict(boxstyle='round', facecolor='wheat', alpha=1)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=15, verticalalignment='top', bbox=props)
    plt.show()

domain_lengths_clean=[]
linker_lengths_clean=[]
protein_lengths_clean=[]
for prot in prot_dict:
    p1 = prot_dict[prot]
    proteinstart = p1.regions[0].start
    proteinend = p1.regions[-1].end
    proteinlaenge = proteinend-proteinstart+1
    protein_lengths_clean.append(proteinlaenge)
    for link in p1.linkers:
        linker_lengths_clean.append(link.length)
    for doma in p1.domains:
        domainresidues = 0
        for domregion in doma.regions:
            domainresidues+=(domregion[0]).length
        domain_lengths_clean.append(domainresidues)

#com_cumul(domain_lengths_clean, linker_lengths_clean)
#protein_linker(domain_lengths_clean,linker_lengths_clean,protein_lengths_clean)


domain_residue_list=[]
linker_residue_list=[]
for key in prot_dict:
    protein = prot_dict[key]
    domain_residue_list.append(protein.domainresidues)
    linker_residue_list.append(protein.linkerresidues)


def lin_dom_res():
    linker_domain_residues = pd.DataFrame(data = {"domain_residues":domain_residue_list,"linker_residues":linker_residue_list}, columns=["domain_residues","linker_residues"])
    t = sns.JointGrid(data=linker_domain_residues,x="domain_residues",y="linker_residues", ylim=[0,400], xlim=[0,1500])
    t.plot_joint(sns.scatterplot,alpha=.2)
    t.plot_marginals(sns.kdeplot)
    ax = sns.regplot(data=linker_domain_residues,x="domain_residues",y="linker_residues", scatter=False, color="orange", ax=t.ax_joint)
    #ldr = sns.scatterplot(data = linker_domain_residues, x = "domain_residues", y = "linker_residues")
    #plt.title("Clean data: proteins with linkers:\n domain linker residue relation")
    plt.suptitle("clean data: if a protein has >=1 linker(s), the total number of linker/domain residues are compared in this plot")
    mdr = np.percentile(domain_residue_list,[50])
    mlr = np.percentile(linker_residue_list,[50])
    print("mean domain residues: ",mdr[0])
    print("mean linker residues: ",mlr[0])
    corr= linker_domain_residues.corr("pearson")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    textstr = "Pearson correlation: "+ str(round(corr.iloc[1, 0], 2))+"\n"+"mean domain residues: "+str(round(mdr[0]))+"\n"+"mean linker residues: "+str(round(mlr[0]))
    t.ax_joint.text(0.05, 0.95, textstr, transform=t.ax_joint.transAxes, fontsize=15, verticalalignment='top', bbox=props)
    plt.ylim([0,150])
    plt.show()





