import csv
import optparse
import re
parser = optparse.OptionParser()
parser.add_option('-q', '--query',
    action="store", dest="path",
    help="path", default="spam")
options, args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as plt
protein_domain_counter_lexikon={}
protein_length_lexikon={}
def readfile(filename):
    filecounter=0
    linecounter = 0
    with open(filename) as bf:
        rd = csv.reader(bf, delimiter="\t", quotechar='"')
        lastproteinint= 1
        lastprotein=""
        commas = 0
        for row in rd:
            linecounter += 1
            protein_gruppe = row[0]
            protein = protein_gruppe[0:5]
            commas = commas + row[1].count(",")
            if (protein in protein_domain_counter_lexikon):
                protein_domain_counter_lexikon.update({protein: protein_domain_counter_lexikon[protein] + 1})
                lastprotein=protein
                lastproteinint=protein_domain_counter_lexikon[protein]
            else:
                protein_domain_counter_lexikon[protein]= 1
                if lastproteinint != 1:
                    print(lastprotein + " " + str(lastproteinint) + " " + str(commas))
                lastprotein = protein
                lastproteinint = 1
                commas = 0

            #hier wird die länge der proteine in ein lexikon gepackt
            length = re.search(r'\d*\Z', row[1]).group(0)
#            print(length)
            if (protein in protein_length_lexikon):
                if protein_length_lexikon[protein] < length:
                    protein_length_lexikon.update({protein:length})
            else:
                protein_length_lexikon[protein]=length
    filecounter += 1


readfile(options.path)


liste = protein_domain_counter_lexikon.values()
domaenen = set(liste)
for i in domaenen:
    if i > 14:
        print()
    summe = sum(value==i for value in protein_domain_counter_lexikon.values())
    print("Proteine mit domaenenanzahl " + str(i) + ":" + str(summe))
    prozent = summe/len(protein_domain_counter_lexikon)
    print(str(prozent*100)[:6] + "%")

plt.hist(liste,100)#bins=len(set(liste)))
plt.title("variations of proteins")
plt.show()
plt.close()

liste = protein_length_lexikon.values()
plt.hist(liste,bins=len(set(liste)))
plt.title("length of proteins")
plt.show()

print ('Query string:', options.path)
