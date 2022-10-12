import csv
import optparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

parser = optparse.OptionParser()
parser.add_option('-q', '--query',
   action="store", dest="path",
   help="path", default="cath-domain-boundaries-seqreschopping.txt")
parser.add_option('-c', '--cutoff',
   action="store", dest="cutoff",
   help="cutoff for linker length", default="40")
parser.add_option('-o', '--output',
   action="store", dest="output",
   help="clean protein dicitonary output", default="cutoff_linkerproteins_Protein_dictionary.pkl")
parser.add_option('-f', '--fasta',
   action="store", dest="fasta",
   help="give out fasta off all proteins prior to cutoff", default="false")
options, args = parser.parse_args()

#diese datei sortiert die linker ab einer bestimmten länge raus und speichert alle proteine in ein lexikon
print("diese datei sortiert die linker ab einer bestimmten länge raus und speichert alle proteine in ein lexikon")
# -q C:\\Users\\yanni\\Downloads\\cath-domain-boundaries-seqreschopping.txt
# philip ording - 99 variationen eines beweises

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
       self.longerthancutoff = False
       self.shortlinker = False
       self.domainresidues = 0
       self.linkerresidues = 0
       self.incutoff = False
       self.sequence=""

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

protein_lexikon={}
def readfile(filename):
   with open(filename) as bf:
       rd = csv.reader(bf, delimiter="\t", quotechar='"')
       for row in rd:
           #zeile einlesen
           protein_gruppe = row[0]
           proteinname = protein_gruppe[0:5]
           domain_to_be_split = row[1]
           domain_for_regions = domain_to_be_split.split(",")
           #wenn protein schon existiert dann hinzufuegen
           if proteinname in protein_lexikon:
               regionlist = list()
               #die regionen in die domain einfuegen
               for region in domain_for_regions:
                   region_prepare = region.split("-")
                   reigonen_hinzufuegen = Region(int(region_prepare[0]), int(region_prepare[1]))
                   regionlist.append(reigonen_hinzufuegen)
               domaene_hinzufuegen = Domain(protein_gruppe[0:7], regionlist)
               placeholder = protein_lexikon[proteinname]
               placeholder.add_domain(domaene_hinzufuegen)
               for region in regionlist:
                   placeholder.regions.append(region)
               placeholder.domain_counter +=1

               protein_lexikon[proteinname] = placeholder
           else:
               regionlist = list()
               # die regionen in die domain einfuegen
               for region in domain_for_regions:
                   region_prepare = region.split("-")
                   reigonen_hinzufuegen = Region(int(region_prepare[0]), int(region_prepare[1]))
                   regionlist.append(reigonen_hinzufuegen)

               domaene_hinzufuegen = Domain(protein_gruppe[0:7], regionlist)
               protein_hinzufuegen = Protein(proteinname,domaene_hinzufuegen,regionlist)
               protein_hinzufuegen.domain_counter =1
               protein_lexikon[proteinname] = protein_hinzufuegen


readfile(options.path)

def write_fasta(fasta):
    raw_protein_fasta = open(fasta,"w")
    for key in protein_lexikon:
        protein = protein_lexikon[key]
        line = ">"+str(key)
        raw_protein_fasta.write(line+"\n")
        raw_protein_fasta.write("placeholder" + "\n")
    raw_protein_fasta.close()

# ########### for test reasons not executed
#if options.fasta != "false":
#    write_fasta(options.fasta)



#TODO   Methode schreiben um regions zu sortieren in der region liste des Proteins
def sort_regions():
   #durch jedes Protein iterieren
   for unsortiertes_protein in protein_lexikon.values():
       #regionen-liste sortieren
       unsortiertes_protein.regions.sort(key=lambda x:x.start, reverse=False)

def relative_position(protein, position):
    #start and end are the start and end of the first and last domain
   protein_regions_start = protein.regions[0]
   protein_regions_end = 0
   for region in protein.regions:
       if region.end > protein_regions_end:
           protein_regions_end = region.end
   protein_start_position = protein_regions_start.start
   protein_end_position = protein_regions_end
   protein_end = protein_end_position-protein_start_position
   position_rel = position-protein_start_position
   relative_position_result = position_rel/protein_end
   return relative_position_result


sort_regions()

domaenen_counter = []
linkercollection = []
linkercollection_relative = []
linker_lengths = []
starting_delay = []
linker_domain_ratio = []
ldr_l = []
ldr_d = []
linker_counter = []
protein_lengths = []
domain_lengths = []
domain_lengths_of_all_proteins = []

#TODO   Methode um linker-regionen zu berechnen
def linker_calculation():
   #durch jedes Protein iterieren
   for key in protein_lexikon:
       protein = protein_lexikon[key]
       protein_lengths.append(protein.regions[-1].end-protein.regions[0].start+1)
       for domaene in protein.domains:
           domaenenlaenge = 0
           for region in domaene.regions:
               domaenenlaenge += region[0].end - region[0].start + 1
           domain_lengths_of_all_proteins.append(domaenenlaenge)

       domaenen_counter.append(protein.domain_counter)


       if protein.regions[0].start > 1:
           starting_delay.append(protein.regions[0].start)
#           if protein.regions[0].start > 1000:
#               print("Protein: " + key + " starts at: " + str(protein.regions[0].start))
       iterators = len(protein.regions)
       if iterators >1:
           iterators -=1
           iterators = range(iterators)
           for i in iterators:
               this_region = protein.regions[i]
               i += 1
               next_region = protein.regions[i]
               nextstart = next_region.start-1
               thisend = this_region.end
               linkerlength = nextstart-thisend
               if linkerlength > 0:
                   if linkerlength < 200:
                       linker_lengths.append(linkerlength)
#                   else:
#                       print(str(linkerlength) + " " + protein.name)
                   if linkerlength > int(options.cutoff):
                       protein.longlinker = True
                   if linkerlength > 40:
                       protein.longerthan40linker = True
                   if linkerlength > int(options.cutoff):
                       protein.longerthancutoff = True
                   if linkerlength < 41:
                       protein.shortlinker = True
                   linker = Region(thisend+1,nextstart)
#                    print("protein " + key + " has a linker at: " + str(linker.start) + ", " + str(linker.end))

                   linkercollection.extend(list(range(linker.start, linker.end+1)))
                   protein.linkers.append(linker)
                   protein.regions.append(linker)
                   # linker positin in the protein is normed to the known structure, not the number position
                   relative_linker_start = relative_position(protein,linker.start)
                   relative_linker_end = relative_position(protein, linker.end)
                   percentage_linker_start = round((relative_linker_start*100))
                   percentage_linker_end = round((relative_linker_end * 100))
#                   if percentage_linker_start ==1 or percentage_linker_end == 100:
#                       print(protein.name)
                   linkercollection_relative.extend(list(range(percentage_linker_start,percentage_linker_end)))
       if len(protein.linkers) != 0:
           linker_domain_ratio.append(len(protein.linkers)/protein.domain_counter)
           ldr_l.append(len(protein.linkers))
           ldr_d.append(protein.domain_counter)
           linker_counter.append(len(protein.linkers))
           for domaene in protein.domains:
               domaenenlaenge = 0
               for region in domaene.regions:
                   domaenenlaenge += region[0].end - region[0].start + 1
               domain_lengths.append(domaenenlaenge)


linker_calculation()
sort_regions()

proteins_with_linkers = {}
# How many proteins have linkers that are longer than 30?
longlinkerproteins = 0
missing_linker_proteins = 0
domain_residues_before_cutoff = 0
linker_residues = 0
domain_residue_list = []
linker_residue_list = []
for key in protein_lexikon:
    protein = protein_lexikon[key]
    if protein.longlinker:
        longlinkerproteins+=1
    if not protein.linkers:
        missing_linker_proteins += 1
    if protein.linkers:
        for dom in protein.regions:
            protein.domainresidues += dom.length
        domain_residue_list.append(protein.domainresidues)
        for lin in protein.linkers:
            protein.linkerresidues += lin.length
        linker_residue_list.append(protein.linkerresidues)
        proteins_with_linkers[key] = protein

cutoff_dictionary = {}
for key in proteins_with_linkers:
    protein = proteins_with_linkers[key]
    if not protein.longlinker:
        cutoff_dictionary[key]=protein




#print("Proteins with linker lengths above 30: " + str(longlinkerproteins))
#print("Proteins with no linker: " + str(missing_linker_proteins))

domaenen = set(domaenen_counter)
for i in domaenen:
   summe = sum(value==i for value in domaenen_counter)
#   print("Proteine mit domaenenanzahl " + str(i) + ":" + str(summe))
   prozent = summe/len(domaenen_counter)
#   print(str(prozent*100)[:6] + "%")

domaenen_counter
linkercollection
starting_delay
linker_domain_ratio



# Data cleanup
#clean_protein_lexikon = {}
#proteins_with_small_linkers=0
#proteins_with_only_small_linkers=0
#for key in protein_lexikon:
#    protein = protein_lexikon[key]
#    if protein.shortlinker and len(protein.linkers)>0:
#        proteins_with_small_linkers+=1
#    if not protein.longerthancutoff and len(protein.linkers)>0:
#        protein.incutoff = True
#        clean_protein_lexikon[key]=protein
#        proteins_with_only_small_linkers+=1

#print("proteins_with_only_small_linkers: " + str(proteins_with_only_small_linkers))

qualified_proteins = 0
qualified_domain_residues = 0
qualified_linker_residues = 0
qualified_domain_count = 0
qualified_linker_count = 0

cutoff_proteins = 0
cutoff_domain_residues = 0
cutoff_linker_residues = 0
cutoff_domain_count = 0
cutoff_linker_count = 0

qualified_proteins = len(proteins_with_linkers)
cutoff_proteins = len(cutoff_dictionary)

for key in proteins_with_linkers:
    protein = proteins_with_linkers[key]
    for domain in protein.domains:
        qualified_domain_count += 1
        for region in domain.regions[0]:
            qualified_domain_residues += region.length
    for linker in protein.linkers:
        qualified_linker_count += 1
        qualified_linker_residues += linker.length

for key in cutoff_dictionary:
    protein = cutoff_dictionary[key]
    for domain in protein.domains:
        cutoff_domain_count += 1
        for region in domain.regions[0]:
            cutoff_domain_residues += region.length
    for linker in protein.linkers:
        cutoff_linker_count += 1
        cutoff_linker_residues += linker.length


print("after eliminating linkerless proteins there are:")
print("qualified_domain_residues: " + str(qualified_domain_residues))
print("qualified_linker_residues: " + str(qualified_linker_residues))
print("qualified_domain_count: " + str(qualified_domain_count))
print("qualified_linker_count: " + str(qualified_linker_count))
print("qualified_proteins: " + str(qualified_proteins))

print("after eliminating linkerless proteins and "+str(options.cutoff)+" cutoff there are:")
print("cutoff_domain_residues: "+str(cutoff_domain_residues))
print("cutoff_linker_residues: "+str(cutoff_linker_residues))
print("cutoff_domain_count: "+str(cutoff_domain_count))
print("cutoff_linker_count: "+str(cutoff_linker_count))
print("cutoff_protein_count: "+str(cutoff_proteins))

def write_to_files():
#    to_michael = open("C:\\Users\\yanni\\Desktop\\Bachelor\\to_michael.txt","w")
#    to_michael.write("Protein id\tstart\tend\n")
#    for key in cutoff_dictionary:
#        protein = cutoff_dictionary[key]
#        start = (protein.regions[0]).start
#        end = (protein.regions[-1]).end
#        write=str(key)+"\t"+str(start)+"\t"+str(end)+"\n"
#        to_michael.write(write)
#    to_michael.close()
    dictionary = open(options.output,"wb")
    pickle.dump(cutoff_dictionary, dictionary)
    dictionary.close()
    dictionary2 = open("raw_protein_dictionary.pkl","wb")
    pickle.dump(protein_lexikon, dictionary2)
    dictionary2.close()

# ########### for test reasons not executed
# write_to_files()
# print("protein_reader finished")

original_anzahl_an_domaenen = set()
original_anzahl_an_linkern = set()
for key in protein_lexikon:
    protein = protein_lexikon[key]
    for domaene in protein.domains:
        original_anzahl_an_domaenen.add(domaene)
    for linker in protein.linkers:
        original_anzahl_an_linkern.add(linker)

#print("sum of domain length: ", sum([Domain.length for Domain in original_anzahl_an_domaenen]))
#print("sum of linker length: ", sum([Linker.length for Linker in original_anzahl_an_linkern]))

#statistics = open("region_statistics.txt","a")
#statistics.write("linker\t"+str(len(linker_residues))+"\tdomain\t"+str(len(domain_residues))+"\tprotein\t"+str((len(linker_residues)+len(domain_residues)))+"\tmin_seq_id\t"+options.min_seq_id+"\tcutoff\t"+options.cutoff)
#original_statistics= open("")

def show_percentages_on_bars(plot):
    for p in plot.patches:
        percentage = p.get_height()
        x = p.get_x() + p.get_width()
        y = p.get_height()
        plot.annotate(percentage, (x, y),ha='center')


def unique_combinations(x,y):
    ctn = []
    vert_clean = []
    hori_clean = []
    ctn_clean = []
    ctn_codes = []
    for idx, v in enumerate(x):
        ctnv = str(v) + " " + str(y[idx])
        ctn.append(ctnv)
    for idx, ctnv in enumerate(ctn):
        if ctnv not in ctn_codes:
            ctn_codes.append(ctnv)
            ctn_clean.append(1)
            vert_clean.append(int((ctnv.split())[0]))
            hori_clean.append(int((ctnv.split())[1]))
        else:
            ctn_clean[ctn_codes.index(ctnv)] += 1
    return vert_clean,hori_clean,ctn_clean

print("check")

def plots():
    # shows how many domains proteins have histogram
    domaenen_plot = sns.histplot(domaenen_counter, edgecolor="black", bins=(np.arange(1, 21) - 0.5),
                                 log_scale=(False, True))
    cntr = 0
    for p in domaenen_plot.patches:
        if cntr < 5:
            domaenen_plot.annotate("%.0f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()),
                                   ha='center', va='center', fontsize=7, color='black', xytext=(7, 5),
                                   textcoords='offset points')
        else:
            domaenen_plot.annotate("%.0f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()),
                                   ha='center', va='center', fontsize=7, color='black', xytext=(0, 5),
                                   textcoords='offset points')
        cntr += 1
    # print(str(domaenen_plot[1][i]) + " " + str(domaenen_plot[0][i]))
    # plt.style.use('fivethirtyeight')
    plt.title("domain counter")
    plt.xlabel("y proteins have x domains")
    plt.tight_layout()
    plt.xticks(range(1, 21))
    plt.show()
    plt.figure(2)




    # DATA VISUALISATION
    plt.rcParams['font.size'] = 15

    f,(ax)= plt.subplots(1)
    dom_cumul = sns.ecdfplot(data = domain_lengths, label="Domains", ax = ax)
    #dom = sns.kdeplot(data = pd.DataFrame({"domain_lengths_capped": domain_lengths_capped},columns=["domain_lengths_capped"]), x="domain_lengths_capped")
    # protein lengths
    lin_cumul = sns.ecdfplot(data = linker_lengths, label="Linkers", ax = ax, color = "orange", complementary = True)
    #pro = sns.kdeplot(data = pd.DataFrame({"protein_lengths_capped": protein_lengths_capped},columns=["protein_lengths_capped"]), x="protein_lengths_capped")
    dom_cumul.set(xticks=[10,25,50,100,250,500])
    ax.legend()
    ps_d = np.percentile(domain_lengths,[5,95,50])
    ps_l = np.percentile(linker_lengths,[5,95,50])
    print("5%-tile cut-off for domains: ", ps_d[0])
    print("95%-tile cut-off for domains: ", ps_d[1])
    len_ten = len([1 for i in domain_lengths if i < 10])
    len_fif = len([1 for i in domain_lengths if i < 50])
    len_four = len([1 for i in domain_lengths if i < 40])
    print(str(len_ten)," of domains <10: ",(len_ten/len(domain_lengths)))
    print(str(len_four)," of domains <40: ",(len_four/len(domain_lengths)))
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
    plt.xlim(0,500)
    plt.title("comparison of cumulative distributions of domain lengths and linker lengths",fontsize=15)
    textstr="5%-tile cut-off for domain length: "+str(ps_d[0])+"\n"+\
            "95%-tile cut-off for domain length: "+str(ps_d[1])+"\n"+\
            "5%-tile cut-off for linker length: "+str(ps_l[0])+"\n"+\
            "95%-tile cut-off for linker length: "+str(ps_l[1])+"\n"+\
            str(len_ten)+" of domains <10: "+str(round((100*(len_ten/len(domain_lengths))),2))+"%"+"\n"+\
            str(len_four)+" of domains <40: "+str(round((100*(len_four/len(domain_lengths))),2))+"%"
    props = dict(boxstyle='round', facecolor='wheat', alpha=1)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=15, verticalalignment='top', bbox=props)
    plt.show()


    #todo only domain lengths from domains that are in proteins with linkers
    #todo find out how many domains are shorter than 30&40
    print("total amount of domains: " + str(len(domain_lengths)))
    print("total amount of linkers: " + str(len(linker_lengths)))
    print("domains shorter than 41: " + str(sum(i <41 for i in domain_lengths)))
    print("domains shorter than 31: " + str(sum(i <31 for i in domain_lengths)))
    print("linkers longer than 40: " + str(sum(i >40 for i in linker_lengths)))
    print("linkers longer than 30: " + str(sum(i >30 for i in linker_lengths)))

    #length correlation between domains and linkers
    f,(ax)= plt.subplots(1)
    dom_hist = sns.histplot(data = domain_lengths, label="Domains", ax = ax, binwidth=5)
    lin_hist = sns.histplot(data = linker_lengths, label="Linkers", ax = ax, binwidth=5, color = "orange")
    dom_hist.set(xticks=[10,25,50,100,250,500])
    ax.legend()
    _, ymax = ax.get_ybound()
    ps_d = np.percentile(domain_lengths,[5,95,50])
    ps_l = np.percentile(linker_lengths,[5,95,50])
    ax.axvline(ps_d[0], label="5%", color="lightblue", linestyle="dashed", linewidth=2)
    ax.axvline(ps_d[1], label="95%", color="blue", linestyle="dashed", linewidth=2)
    label = "domain length mean " + str(round(np.mean(domain_lengths)))
    ax.axvline(round(np.mean(domain_lengths)), label=label, color="black", linestyle="-", linewidth=2)
    ax.axvline(ps_l[0], label="5%", color="pink", linestyle="dashed", linewidth=2)
    ax.axvline(ps_l[1], label="95%", color="red", linestyle="dashed", linewidth=2)
    label = "linker length mean " + str(round(np.mean(linker_lengths)))
    ax.axvline(round(np.mean(linker_lengths)), label=label, color="grey", linestyle="-", linewidth=2)
    plt.legend()
    plt.xlim(0,500)
    plt.title("Length correlation between domains and linkers")
    plt.show()



    domain_lengths_sorted = sorted(domain_lengths)
    protein_lengths_sorted = sorted(protein_lengths)
    domain_lengths_capped = domain_lengths_sorted[:len(domain_lengths_sorted)-(round(0.05*len(domain_lengths_sorted)))]
    protein_lengths_capped = protein_lengths_sorted[:len(protein_lengths_sorted)-(round(0.05*len(protein_lengths_sorted)))]
    f,(ax)= plt.subplots(1)
    # domain lengths
    dom = sns.kdeplot(data = domain_lengths, label="Domains", ax = ax, color="blue")
    #dom = sns.kdeplot(data = pd.DataFrame({"domain_lengths_capped": domain_lengths_capped},columns=["domain_lengths_capped"]), x="domain_lengths_capped")
    # protein lengths
    pro = sns.kdeplot(data = protein_lengths, label="Proteins", ax = ax, color="red")
    #pro = sns.kdeplot(data = pd.DataFrame({"protein_lengths_capped": protein_lengths_capped},columns=["protein_lengths_capped"]), x="protein_lengths_capped")
    # protein lengths
    pro = sns.kdeplot(data=linker_lengths, label="Linkers", ax=ax, color="orange")
    dom.set(xticks=[10,25,50,100,250,500])
    pro.set(xticks=[10,25,50,100,250,500])
    plt.xlim([0,500])
    ax.legend()
    mean_p = np.percentile(protein_lengths,[50])
    mean_d = np.percentile(domain_lengths,[50])
    mean_l = np.percentile(linker_lengths,[50])
    textstr = "median protein length: " + str(round(mean_p[0])) + "\n" + \
              "median domain length: " + str(round(mean_d[0])) + "\n" + \
              "median linker length: " + str(round(mean_l[0]))
    props = dict(boxstyle='round', facecolor='wheat', alpha=1)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=15, verticalalignment='top', bbox=props)
    plt.grid()
    plt.title("Comparison of Domain, Protein and Linker lengths")
    plt.show()



    # relative linker position
    #f,(ax)= plt.subplots(1)
    fig, ax = plt.subplots(1)
    hist1 = sns.histplot(data=pd.DataFrame({"relative linker positions": linkercollection_relative},
                                           columns=["relative linker positions"]), x="relative linker positions",
                         stat="frequency", bins=np.arange(0, 101, 1), kde=True)
    total = float(len(linkercollection_relative))
    # _, ymax = ax.get_ybound()
    # test = (ax.containers[0])
    # ax.axhline(10000, label="12.5%", color="lightblue", linestyle="dashed", linewidth=2)
    plt.subplots_adjust(bottom=0.15, left=0.14)
    plt.title("relative linker position in the protein")
    fig.savefig("relative.png", dpi=300)
    plt.show()



    #f,(ax)=plt.subplots(1)
    linker_domain_residues = pd.DataFrame(data = {"domain_residues":domain_residue_list,"linker_residues":linker_residue_list}, columns=["domain_residues","linker_residues"])
    t = sns.JointGrid(data=linker_domain_residues,x="domain_residues",y="linker_residues", ylim=[0,400], xlim=[0,1500])
    t.plot_joint(sns.scatterplot,alpha=.2)
    t.plot_marginals(sns.kdeplot)
    ax = sns.regplot(data=linker_domain_residues,x="domain_residues",y="linker_residues", scatter=False, color="orange", ax=t.ax_joint)
    #ldr = sns.scatterplot(data = linker_domain_residues, x = "domain_residues", y = "linker_residues")
    plt.title("proteins with linkers:\n domain linker residue relation")
    plt.suptitle("if a protein has >=1 linker(s), the total number of linker/domain residues are compared in this plot")
    mdr = np.percentile(domain_residue_list,[50])
    mlr = np.percentile(linker_residue_list,[50])
    print("mean domain residues: ",mdr[0])
    print("mean linker residues: ",mlr[0])
    corr= linker_domain_residues.corr("pearson")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    textstr = "Pearson correlation: "+ str(round(corr.iloc[1, 0], 2))+"\n"+"mean domain residues: "+str(round(mdr[0]))+"\n"+"mean linker residues: "+str(round(mlr[0]))
    t.ax_joint.text(0.05, 0.95, textstr, transform=t.ax_joint.transAxes, fontsize=15, verticalalignment='top', bbox=props)
    plt.show()

    dp_data = pd.DataFrame(data={"number of domain(s)": domaenen_counter, "protein lengths": protein_lengths}, columns=["number of domain(s)", "protein lengths"])
    t = sns.JointGrid(data=dp_data, x="protein lengths", y="number of domain(s)",xlim=[0,1500])
    t.plot_joint(sns.scatterplot)
    sns.kdeplot(data=dp_data, x="protein lengths", ax=t.ax_marg_x, fill=True)
    sns.histplot(data=dp_data, y="number of domain(s)", binwidth=1, ax=t.ax_marg_y)
    sns.regplot(data=dp_data, x="protein lengths", y="number of domain(s)", order=3, ax = t.ax_joint, color="r", scatter=False)
    corr= dp_data.corr("pearson")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    textstr = "Pearson correlation: "+ str(round(corr.iloc[1, 0], 2))+"\n"+"mean protein length: "+str(round((np.percentile(protein_lengths,[50]))[0]))
    t.ax_joint.text(0.05, 0.95, textstr, transform=t.ax_joint.transAxes, fontsize=15, verticalalignment='top', bbox=props)
    plt.ylim([0,10])
    #plt.title("number of domains to protein length correlation")
    plt.show()

    #linkerperdomain = sns.scatterplot(data=dp_data, x="protein lenghths", y="number of domain(s)")
    #regression_line = sns.regplot(data=dp_data, x="protein lenghths", y="number of domain(s)")

    vert_clean,hori_clean,ctn_clean = unique_combinations(ldr_l, ldr_d)
    lpd_data = pd.DataFrame(data = {"number of linker(s)":vert_clean,"number of domain(s)":hori_clean,"counts":ctn_clean},columns=["number of linker(s)","number of domain(s)","counts"])
    cmap = sns.cubehelix_palette(rot=-1, as_cmap=True)
    g = sns.scatterplot(
        data=lpd_data,
        y="number of linker(s)", x="number of domain(s)",
        size="counts", hue="counts",
        palette=cmap, sizes=(30,300)
    )
    regression_line = sns.regplot(scatter=False, data=pd.DataFrame(data = {"number of linker(s)":ldr_l, "number of domain(s)":ldr_d}, columns = ["number of linker(s)","number of domain(s)"]), y = "number of linker(s)", x="number of domain(s)", ci=None)
    plt.title("analysis of proteins with linkers")
    print("analysis of proteins with linkers")
    pd.DataFrame(data = {"number of linker(s)":ldr_l, "number of domain(s)":ldr_d}, columns = ["number of linker(s)","number of domain(s)"]).corr("pearson")
    #g.set(xscale="linear", yscale="linear")
    #g.ax.xaxis.grid(linewidth=.25)
    #g.ax.yaxis.grid(linewidth=.25)
    plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,20])
    plt.legend(prop={'size': 7})
    plt.xlim(0,10.5)
    plt.show()





    plt.figure(3)
    relative_linker_plot = plt.hist(linkercollection_relative, edgecolor="black", bins=100, log=False)
    for i in range(100):
       plt.text(relative_linker_plot[1][i],relative_linker_plot[0][i],str(round(relative_linker_plot[0][i])))
    plt.style.use('fivethirtyeight')
    plt.title("position relative linker distribution")
    plt.xlabel("what relative positions are linkers")
    plt.show()
    plt.figure(4)
    linker_length_plot = plt.hist(linker_lengths, edgecolor="black", bins=100, log=True)
    for i in range(100):
       plt.text(linker_length_plot[1][i],linker_length_plot[0][i],str(round(linker_length_plot[0][i])))
    plt.style.use('fivethirtyeight')
    plt.title("what lengths are the linkers")
    plt.xlabel("linker length")
    plt.show()
    plt.figure(6)
    linker_counter_plot = plt.hist(linker_counter, edgecolor="black", bins=100, log=True)
    for i in range(100):
       plt.text(linker_counter_plot[1][i],linker_counter_plot[0][i],str(round(linker_counter_plot[0][i])))
    plt.style.use('fivethirtyeight')
    plt.title("proteins with linkers")
    plt.xlabel("amount of linkers")
    plt.show()
    plt.figure(7)

    fig, ax = plt.subplots(1)
    # domaenen_plot = plt.hist(starting_delay, edgecolor="black", bins=70, log=True)
    delay_plot = sns.histplot(starting_delay, edgecolor="black", bins=75, log_scale=(False, True))
    cntr = 0
    for p in delay_plot.patches:
        if cntr < 7:
            # (p.get_x() + p.get_width() / 2., p.get_height())
            delay_plot.annotate("%.0f" % p.get_height(), ((p.get_x() + p.get_width() / 2.) + 0.2, p.get_height()),
                                   xytext=(p.get_x() + 0.3, 30), ha='center', va='center', fontsize=10, color='black',
                                   textcoords='offset points',
                                   arrowprops=dict(facecolor='black', shrink=0.05, headwidth=5, width=1))
        if cntr > 13:
            delay_plot.annotate("%.0f" % p.get_height(), ((p.get_x() + p.get_width() / 2.), p.get_height()),
                                   xytext=(0, 5), ha='center', va='center', fontsize=10, color='black',
                                   textcoords='offset points')
        cntr += 1
    # plt.title("protein ")
    plt.xlabel("start index of first domain")
    plt.ylabel("protein count")
    plt.tight_layout()
    # plt.xticks(range(1, 21))
    fig.tight_layout()
    plt.subplots_adjust(bottom=0.15, left=0.14)
    fig.savefig("index.png", dpi=300)
    plt.show()





