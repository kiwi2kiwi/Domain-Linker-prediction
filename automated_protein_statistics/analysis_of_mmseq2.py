import csv
import optparse
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from Bio import SeqIO
import re

parser = optparse.OptionParser()
parser.add_option('-d', '--dictionary',
   action="store", dest="dictionary",
   help="protein dictionary", default="clean_pdb_Protein_dictionary.pkl")
parser.add_option('-r', '--raw_dictionary',
   action="store", dest="raw_dictionary",
   help="raw protein dictionary", default="raw_Protein_dictionary.pkl")
parser.add_option('-i', '--redundancy_reduced_file',
   action="store", dest="input",
   help="redundancy reduced file", default="domain_PIDE20_rep_seq.fasta")
parser.add_option('-c', '--cutoff',
   action="store", dest="cutoff",
   help="cutoff linker length", default="40")
options, args = parser.parse_args()


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

def linkershorterthan(protein,cutoff):
    if protein.linkers:
        for linker in protein.linkers:
            if linker.length>int(cutoff):
                return False
    else:
        return False
    return True

dictionary = open(options.dictionary,"rb")
prot_dict = pickle.load(dictionary)

raw_dictionary = open(options.raw_dictionary,"rb")
raw_prot_dict = pickle.load(raw_dictionary)

WD = Path(__file__).resolve().parents[1]
print(WD)
til = str(WD / 'Prediction' / "train_input_list_clean.pkl")
input_file = "C:\\Users\\yanni\\Desktop\\Git_bachelors_thesis\\automated_protein_statistics\\cutoff_40_redundancy_0.20.fasta"
#options.input
print("now using: "+input_file)
stream = open(input_file,"r")
set = set()
id = ""
for line in stream:
    if line.startswith(">"):
        id = line.split()
        id = id[0].replace("_","").replace(">","")
        set.add(id)


#prot_dict_mmseq = {}
#for prot_id in set:
#    if prot_id in prot_dict:
#        protein_to_fetch = prot_dict[prot_id]
#        if "X" not in protein_to_fetch.sequence:
#            prot_dict_mmseq[prot_id]= protein_to_fetch
#        else:
#            print(prot_id + " has invalid residue X")

pipeline = True

disqualified_proteins = 0
disqualified_domain_residues = 0
disqualified_linker_residues = 0
disqualified_domain_count = 0
disqualified_linker_count = 0

reduced_domain_residues = 0
reduced_linker_residues = 0
reduced_domain_count = 0
reduced_linker_count = 0

cutoff_domain_residues = 0
cutoff_linker_residues = 0
cutoff_domain_count = 0
cutoff_linker_count = 0

redundancy_prot_dict_mmseq = {}
cutoff_prot_dict_mmseq = {}
missing_sequences=0
for prot_id in set:
    if prot_id in raw_prot_dict:
        protein_to_fetch = raw_prot_dict[prot_id]
        if pipeline:
            if protein_to_fetch.linkers:
                disqualified_proteins += 1
                for domain in protein_to_fetch.domains:
                    disqualified_domain_count += 1
                    for region in domain.regions[0]:
                        disqualified_domain_residues += region.length
                for linker in protein_to_fetch.linkers:
                    disqualified_linker_count += 1
                    disqualified_linker_residues += linker.length
            for domain in protein_to_fetch.domains:
                reduced_domain_count+=1
                for region in domain.regions[0]:
                    reduced_domain_residues += region.length
            for linker in protein_to_fetch.linkers:
                reduced_linker_count+=1
                reduced_linker_residues+=linker.length
            if linkershorterthan(protein_to_fetch,options.cutoff):
                for domain in protein_to_fetch.domains:
                    cutoff_domain_count+=1
                    for region in domain.regions[0]:
                        cutoff_domain_residues += region.length
                for linker in protein_to_fetch.linkers:
                    cutoff_linker_count+=1
                    cutoff_linker_residues+=linker.length
        if protein_to_fetch.sequence == "":
            missing_sequences+=1
        if "X" not in protein_to_fetch.sequence:
            redundancy_prot_dict_mmseq[prot_id]= protein_to_fetch
            if linkershorterthan(protein_to_fetch,options.cutoff):
                cutoff_prot_dict_mmseq[prot_id]=protein_to_fetch
        else:
            print(prot_id + " has invalid residue X")


#m = re.search("\.(\d+)\_", options.input)
#min_seq_id = str("0."+m.group(1))
#print("after "+str(min_seq_id)+" redundancy reduction there are:")
print("reduced_domain_residues: "+str(reduced_domain_residues))
print("reduced_linker_residues: "+str(reduced_linker_residues))
print("reduced_domain_count: "+str(reduced_domain_count))
print("reduced_linker_count: "+str(reduced_linker_count))
print("reduced_protein_count: "+str(len(redundancy_prot_dict_mmseq)))

if pipeline:
    print("after redundancy reduction and eliminating linkerless proteins there are:")
    print("disqualified_domain_residues: "+str(disqualified_domain_residues))
    print("disqualified_linker_residues: "+str(disqualified_linker_residues))
    print("disqualified_domain_count: "+str(disqualified_domain_count))
    print("disqualified_linker_count: "+str(disqualified_linker_count))
    print("disqualified_proteins: "+str(disqualified_proteins))

    print("after redundancy reduction, eliminating linkerless proteins and "+str(options.cutoff)+" cutoff there are:")
    print("cutoff_domain_residues: "+str(cutoff_domain_residues))
    print("cutoff_linker_residues: "+str(cutoff_linker_residues))
    print("cutoff_domain_count: "+str(cutoff_domain_count))
    print("cutoff_linker_count: "+str(cutoff_linker_count))
    print("cutoff_protein_count: "+str(len(cutoff_prot_dict_mmseq)))

dictionary = open("40_cutoff_then_0_20.pkl","wb")
pickle.dump(cutoff_prot_dict_mmseq, dictionary)
dictionary.close()



plt.rcParams['font.size'] = 10

def protein_linker(domain_lengths,linker_lengths,protein_lengths,title):
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
    plt.title(title)
    plt.show()

def com_cumul(domain_lengths, linker_lengths,title):
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
    plt.title(title, fontsize=15)
    textstr = "5%-tile cut-off for domain length: " + str(round(ps_d[0])) + "\n" + \
              "95%-tile cut-off for domain length: " + str(round(ps_d[1])) + "\n" + \
              "5%-tile cut-off for linker length: " + str(round(ps_l[0])) + "\n" + \
              "95%-tile cut-off for linker length: " + str(round(ps_l[1])) + "\n" + \
              str(len_ten) + " of domains <10: " + str(round((100 * (len_ten / len(domain_lengths))), 2)) + "%" + "\n" + \
              str(len_four) + " of domains <40: " + str(round((100 * (len_four / len(domain_lengths))), 2)) + "%"
    props = dict(boxstyle='round', facecolor='wheat', alpha=1)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=15, verticalalignment='top', bbox=props)
    plt.show()

redundancy_domain_lengths=[]
redundancy_linker_lengths=[]
redundancy_protein_lengths=[]
for prot in redundancy_prot_dict_mmseq:
    p1 = redundancy_prot_dict_mmseq[prot]
    proteinstart = p1.regions[0].start
    proteinend = p1.regions[-1].end
    proteinlaenge = proteinend-proteinstart+1
    redundancy_protein_lengths.append(proteinlaenge)
    for link in p1.linkers:
        redundancy_linker_lengths.append(link.length)
    for doma in p1.domains:
        domainresidues = 0
        for domregion in doma.regions:
            domainresidues+=(domregion[0]).length
        redundancy_domain_lengths.append(domainresidues)

cutoff_domain_lengths=[]
cutoff_linker_lengths=[]
cutoff_protein_lengths=[]
for prot in cutoff_prot_dict_mmseq:
    p1 = cutoff_prot_dict_mmseq[prot]
    proteinstart = p1.regions[0].start
    proteinend = p1.regions[-1].end
    proteinlaenge = proteinend-proteinstart+1
    cutoff_protein_lengths.append(proteinlaenge)
    for link in p1.linkers:
        cutoff_linker_lengths.append(link.length)
    for doma in p1.domains:
        domainresidues = 0
        for domregion in doma.regions:
            domainresidues+=(domregion[0]).length
        cutoff_domain_lengths.append(domainresidues)

#protein_linker(redundancy_domain_lengths,redundancy_linker_lengths,redundancy_protein_lengths,"Comparison of Domain, Protein and Linker lengths on the mmseq2 data")
com_cumul(redundancy_domain_lengths,redundancy_linker_lengths,"comparison of cumulative distributions of domain lengths and linker lengths on the mmseq2 data")
#protein_linker(cutoff_domain_lengths,cutoff_linker_lengths,cutoff_protein_lengths,"Comparison of Domain, Protein and Linker lengths on the mmseq2 + cutoff data")
com_cumul(cutoff_domain_lengths,cutoff_linker_lengths,"comparison of cumulative distributions of domain lengths and linker lengths on the mmseq2 + cutoff data")


def linkertodomain_ratio():
    domain_residue_list = []
    linker_residue_list = []
    for key in cutoff_prot_dict_mmseq:
        protein = cutoff_prot_dict_mmseq[key]
        domain_residue_list.append(protein.domainresidues)
        linker_residue_list.append(protein.linkerresidues)

    linker_domain_residues = pd.DataFrame(
        data={"domain_residues": domain_residue_list, "linker_residues": linker_residue_list},
        columns=["domain_residues", "linker_residues"])
    t = sns.JointGrid(data=linker_domain_residues, x="domain_residues", y="linker_residues", ylim=[0, 400],
                      xlim=[0, 1500])
    t.plot_joint(sns.scatterplot, alpha=.2)
    t.plot_marginals(sns.kdeplot)
    ax = sns.regplot(data=linker_domain_residues, x="domain_residues", y="linker_residues", scatter=False,
                     color="orange", ax=t.ax_joint)
    # ldr = sns.scatterplot(data = linker_domain_residues, x = "domain_residues", y = "linker_residues")
    # plt.title("Clean data: proteins with linkers:\n domain linker residue relation")
    plt.suptitle("mmseq2 data:the total number of linker/domain residues are compared in this plot")
    mdr = np.percentile(domain_residue_list, [50])
    mlr = np.percentile(linker_residue_list, [50])
    print("mean domain residues: ", mdr[0])
    print("mean linker residues: ", mlr[0])
    corr = linker_domain_residues.corr("pearson")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    textstr = "Pearson correlation: " + str(round(corr.iloc[1, 0], 2)) + "\n" + "mean domain residues: " + str(
        round(mdr[0])) + "\n" + "mean linker residues: " + str(round(mlr[0]))
    t.ax_joint.text(0.05, 0.95, textstr, transform=t.ax_joint.transAxes, fontsize=15, verticalalignment='top',
                    bbox=props)
    plt.ylim([0, 150])
    plt.show()

#linkertodomain_ratio()

def analyze_text(text, letter='e'):
    n = len([x for x in text if x.isalpha()])
    freq = text.count(letter)
    percent = round(float(freq) / n * 100)
    text = "The text contains {} alphabetic characters, of which {} ({}%) are '{}'.".format(n, freq, percent, letter)
#    print(text)
    return percent

cutoff_domain_residues=""
cutoff_linker_residues=""
for prot in cutoff_prot_dict_mmseq:
    protein = cutoff_prot_dict_mmseq[prot]
    for domain in protein.domains:
        for region in domain.regions[0]:
            cutoff_domain_residues += region.sequence
    for linker in protein.linkers:
        if linker.length < 6 and linker.length > 1:
            cutoff_linker_residues += linker.sequence

print("total domain residue count: "+str(len(cutoff_domain_residues)))
print("total linker residue count: "+str(len(cutoff_linker_residues)))

#print("linker residues left: " + str(len(cutoff_linker_residues)))
#print("domain residues left: " + str(len(cutoff_domain_residues)))
#print("protein residues left: " + str((len(cutoff_linker_residues)+len(cutoff_domain_residues))))
#print(str(len(prot_dict_mmseq)))

#statistics = open("residue_statistics.txt","a")
#statistics.write("linker\t" + str(len(cutoff_linker_residues)) +"\tdomain\t" + str(len(cutoff_domain_residues)) +"\tprotein\t" + str((len(cutoff_linker_residues)+len(cutoff_domain_residues))) +"\tmin_seq_id\t" + str(min_seq_id) +"\tcutoff\t" + str(options.cutoff) +"\tlinker_count\t" + str(len(redundancy_linker_lengths)) + "\tdomain_count\t" + str(len(redundancy_domain_lengths)) + "\tprotein_count\t" + str(len(prot_dict_mmseq)))
#statistics.close()
#print("protein statistics done")

def split(word):
    return [char for char in word]

dom_AA_dict = {}
lin_AA_dict = {}
AAs = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
for AA in AAs:
    dom_AA_dict[AA] = analyze_text(cutoff_domain_residues,AA)
    lin_AA_dict[AA] = analyze_text(cutoff_linker_residues,AA)

dom_AA = pd.DataFrame(data = {"A": ["A",dom_AA_dict["A"],"domain"],
                              "R": ["R",dom_AA_dict["R"],"domain"],
                              "N": ["N",dom_AA_dict["N"],"domain"],
                              "D": ["D",dom_AA_dict["D"],"domain"],
                              "C": ["C",dom_AA_dict["C"],"domain"],
                              "E": ["E",dom_AA_dict["E"],"domain"],
                              "Q": ["Q",dom_AA_dict["Q"],"domain"],
                              "G": ["G",dom_AA_dict["G"],"domain"],
                              "H": ["H",dom_AA_dict["H"],"domain"],
                              "I": ["I",dom_AA_dict["I"],"domain"],
                              "L": ["L",dom_AA_dict["L"],"domain"],
                              "K": ["K",dom_AA_dict["K"],"domain"],
                              "M": ["M",dom_AA_dict["M"],"domain"],
                              "F": ["F",dom_AA_dict["F"],"domain"],
                              "P": ["P",dom_AA_dict["P"],"domain"],
                              "S": ["S",dom_AA_dict["S"],"domain"],
                              "T": ["T",dom_AA_dict["T"],"domain"],
                              "W": ["W",dom_AA_dict["W"],"domain"],
                              "Y": ["Y",dom_AA_dict["Y"],"domain"],
                              "V": ["V",dom_AA_dict["V"],"domain"]})
lin_AA = pd.DataFrame(data = {"A": ["A",lin_AA_dict["A"],"linker"],
                              "R": ["R",lin_AA_dict["R"],"linker"],
                              "N": ["N",lin_AA_dict["N"],"linker"],
                              "D": ["D",lin_AA_dict["D"],"linker"],
                              "C": ["C",lin_AA_dict["C"],"linker"],
                              "E": ["E",lin_AA_dict["E"],"linker"],
                              "Q": ["Q",lin_AA_dict["Q"],"linker"],
                              "G": ["G",lin_AA_dict["G"],"linker"],
                              "H": ["H",lin_AA_dict["H"],"linker"],
                              "I": ["I",lin_AA_dict["I"],"linker"],
                              "L": ["L",lin_AA_dict["L"],"linker"],
                              "K": ["K",lin_AA_dict["K"],"linker"],
                              "M": ["M",lin_AA_dict["M"],"linker"],
                              "F": ["F",lin_AA_dict["F"],"linker"],
                              "P": ["P",lin_AA_dict["P"],"linker"],
                              "S": ["S",lin_AA_dict["S"],"linker"],
                              "T": ["T",lin_AA_dict["T"],"linker"],
                              "W": ["W",lin_AA_dict["W"],"linker"],
                              "Y": ["Y",lin_AA_dict["Y"],"linker"],
                              "V": ["V",lin_AA_dict["V"],"linker"]})
def show_percentages_on_bars(plot):
    for p in plot.patches:
        percentage = p.get_height()
        x = p.get_x() + p.get_width()-.2
        y = p.get_height()
        plot.annotate(int(percentage), (x, y),ha='center')

def plot_residue_comparison(dom_AA,lin_AA):
    fig, ax = plt.subplots(1)
    dom_AA=dom_AA.T
    lin_AA=lin_AA.T
    combined = pd.concat([dom_AA,lin_AA],axis=0)
    combined.columns=["residue","percentage","region"]
    ax=sns.barplot(data = combined, y="percentage", x="residue", hue="region")
    #iterk = pd.DataFrame(combined["percentage"])
    #for index, row in iterk.iterrows():
    #    ax.text(row.name, row.percentage, row.percentage, color='black', ha="center")
    show_percentages_on_bars(ax)
    plt.title("percentage residue composition of domain and linker \nat linker length cutoff: 5 no 1")
    #plt.suptitle("total domain residue count: "+str(len(cutoff_domain_residues))+"\ntotal linker residue count: "+str(len(cutoff_linker_residues)))
    plt.legend()
    fig.savefig("anom_5_no_1.png", dpi=300)
    plt.show()

plot_residue_comparison(dom_AA,lin_AA)

print("end")
