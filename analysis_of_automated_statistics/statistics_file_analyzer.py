import csv
import optparse
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import style
import re


parser = optparse.OptionParser()
parser.add_option('-s', '--statistics',
   action="store", dest="statistics",
   help="statistics file from analysis_of_mmseq2", default="..\\automated_protein_statistics\\residue_statistics_demo.txt")
parser.add_option('-t', '--type',
   action="store", dest="type",
   help="what should the z axis depict? options are:protein, linker and domain", default="linker")
options, args = parser.parse_args()


#statistics = open("residue_statistics.txt","a")
#statistics.write("linker\t"+str(len(linker_residues))+"\tdomain\t"+str(len(domain_residues))+"\tprotein\t"+str((len(linker_residues)+len(domain_residues)))+"\tmin_seq_id\t"+str(options.min_seq_id)+"\tcutoff\t"+str(options.cutoff)+"\tlinker_count\t"+str(linker_lengths_clean)+"\tdomain_count\t"+str(domain_lengths_clean)+"\tprotein_count\t"+str(len(prot_dict_mmseq)))

statistics = open(options.statistics,"r")

stats = pd.DataFrame(columns=["linker residues","domain residues","protein residues","min_seq_id","cutoff","linker count","domain count","protein count"])
for line in statistics:
    line = line.split()
    zeile = pd.Series([float(line[1]),float(line[3]),float(line[5]),float(line[7]),float(line[9]),float(line[11]),float(line[13]),float(line[15])],index=["linker residues","domain residues","protein residues","min_seq_id","cutoff","linker count","domain count","protein count"])
    stats = stats.append(zeile,ignore_index=True)


from matplotlib import cm
def show_in_3d(x_ax,y_ax,z_ax):
    X = stats[x_ax]
    Y = stats[y_ax]
    Z = stats[z_ax]
    # Normalize the colors based on Z value
    norm = plt.Normalize(Z.min(), Z.max())
    colors = cm.jet(norm(Z))
    ax = plt.axes(projection='3d', proj_type="ortho")
    ax.scatter(X, Y, Z, c=Z, cmap='BrBG', linewidth=1)
    plt.rcParams["figure.figsize"] = 12.8, 9.6
    plt.title(z_ax)
    plt.xlabel(x_ax)
    plt.ylabel(y_ax)
    ax.set_zlabel(z_ax)
    plt.show()

show_in_3d("cutoff", "min_seq_id", str(options.type)+" residues")
show_in_3d("cutoff", "min_seq_id", str(options.type)+" count")
print("3D visualization done")
