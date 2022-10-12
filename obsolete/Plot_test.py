import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from itertools import chain, combinations
from collections import Counter

regions = []
regions.append("test1")
regions.append("test2")


p=sns.load_dataset("planets")
vert = [1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 5, 5, 7, 7, 4, 3, 3, 7, 6, 5, 3, 7, 9, 3, 2, 6, 6, 6, 4, 3]
hori = [6, 3, 5, 5, 3, 1, 2, 3, 2, 3, 4, 4, 7, 7, 1, 5, 9, 2, 7, 7, 4, 3, 1, 3, 1, 2, 6, 6, 6, 1]
comb = {}
ctn = []
vert_clean = []
hori_clean = []
ctn_clean = []
ctn_codes = []
for idx,v in enumerate(vert):
    ctnv = str(v) +" "+ str(hori[idx])
    ctn.append(ctnv)
for idx, ctnv in enumerate(ctn):
    if ctnv not in ctn_codes:
        ctn_codes.append(ctnv)
        ctn_clean.append(1)
        vert_clean.append((ctnv.split())[0])
        hori_clean.append((ctnv.split())[1])
    else:
        ctn_clean[ctn_codes.index(ctnv)]+=1

linkers = ['1', '1', '2', '2', '3', '1', '2', '2', '1', '1', '4', '1', '2', '3', '3', '5', '2', '6', '5', '3', '5', '3', '3', '6', '4', '2', '1', '4', '1', '2', '1', '1', '4', '3', '3', '5', '2', '2', '4', '2', '5', '7', '4', '1', '1', '4', '3', '6', '5', '3', '2', '7', '1', '6', '4', '5', '4']
domains = ['2', '3', '2', '4', '4', '4', '3', '7', '6', '5', '4', '8', '6', '5', '3', '5', '5', '4', '3', '2', '6', '6', '7', '3', '3', '8', '9', '5', '1', '9', '10', '7', '2', '13', '10', '10', '1', '11', '6', '13', '4', '7', '10', '13', '11', '7', '8', '5', '12', '1', '20', '8', '20', '8', '8', '8', '14']
counts = [8437, 3547, 811, 515, 210, 1398, 1154, 26, 360, 328, 159, 37, 60, 91, 194, 7, 275, 18, 1259, 100, 17, 12, 1, 76, 201, 13, 17, 3, 388, 6, 16, 34, 5, 2, 2, 1, 24, 2, 17, 22, 12, 1, 6, 2, 2, 2, 5, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1]
#planets = pd.DataFrame(data = {"vert":vert_clean,"hori":hori_clean,"ctn":ctn_clean},columns=["vert","hori","ctn"])
planets = pd.DataFrame(data = {"vert":vert,"hori":hori},columns=["vert","hori"])
corr = planets.corr("pearson")
print(round(corr.iloc[1,0],2))
t = sns.JointGrid(data=planets, x="vert", y="hori", xlim=[0,10])
t.plot_joint(sns.scatterplot, s=100, alpha=.5)
sns.kdeplot(data=planets, x="vert", ax = t.ax_marg_x)
sns.histplot(data=planets, y="hori", ax = t.ax_marg_y, binwidth=1)
sns.regplot(data=planets, x="vert", y="hori", ax = t.ax_joint)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr = "Pearson correlation: "+ str(round(corr.iloc[1, 0], 2))
t.ax_joint.text(0.05, 0.95, textstr, transform=t.ax_joint.transAxes, fontsize=10, verticalalignment='top', bbox=props)
plt.show()





cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
g = sns.relplot(
    data=planets,
    x="vert", y="hori", hue="ctn",
    palette=cmap, kind="scatter")
g.set(xscale="linear", yscale="linear")
g.ax.xaxis.grid(True, "minor", linewidth=.25)
g.ax.yaxis.grid(True, "minor", linewidth=.25)
g.despine(left=True, bottom=True)
plt.gca().invert_yaxis()
regp = pd.DataFrame(data = {"vert":vert,"hori":hori},columns=["vert","hori"])
sns.regplot(
    data = regp,
    x="vert", y="hori",
    order=3
)
plt.show()

penguins = sns.load_dataset("penguins")

f, (ax) = plt.subplots()
label_test1 = pd.DataFrame(data = {"one":[1,1,1,1,3]})
label_test2 = pd.DataFrame(data = {"two":[2,2,4,4,4]})
lt1 = sns.ecdfplot(data = label_test1, ax=ax, complementary=True, color="orange")
lt2 = sns.ecdfplot(data = label_test2, ax=ax)
ax.legend(["one","two"])
ps_d = np.percentile(label_test1,[5,95])
ps_l = np.percentile(label_test2,[5,95])
_, ymax = ax.get_ybound()
#heights_d = lt1[0][np.searchsorted(lt1[1], ps_d, side="left")-1]
#heights_l = lt2[0][np.searchsorted(lt2[1], ps_l, side="left")-1]
ax.axvline(ps_d[0], label="5%", color="lightblue", linestyle="dashed", linewidth=2)
ax.axvline(ps_d[1], label="95%", color="blue", linestyle="dashed", linewidth=2)
ax.axvline(ps_l[0], label="5%", color="pink", linestyle="dashed", linewidth=2)
ax.axvline(ps_l[1], label="95%", color="red", linestyle="dashed", linewidth=2)
plt.legend()
plt.show()


# scatterplot
bill_depth = penguins.loc[:,"bill_depth_mm"]
bill_length = penguins.loc[:,"bill_length_mm"]
species = penguins.loc[:,"species"]

# create a dataframe with lists inserted as columns, put in a dictionary
dataframe = pd.DataFrame(data = {"bill_depth_mm":bill_depth,"bill_length_mm":bill_length,"species":species},columns = ["bill_depth_mm","bill_length_mm","species"])
# create a dataframe with lists inserted as rows, put in lists
dataframe2 = pd.DataFrame(data = [bill_depth,bill_length,species], index = ["bill_depth_mm","bill_length_mm","species"])
# seaborn only accepts dataframes in which it can search for column names
penguins_plot = sns.scatterplot(data = dataframe, x = "bill_depth_mm", y = "bill_length_mm", hue = "species")
penguins_plot.set(xticks=[15,17,19,21,23,25,27,29])
plt.show()

# to display multiple classes in multiple plots, leave out multiple and put col = hue
penguins_plot2 = sns.displot(data=penguins, x="bill_depth_mm", hue="island", col="island", kind="kde")
plt.show()

plt.figure(0)
penguins_plot3 = sns.kdeplot(data=penguins, x="bill_depth_mm", hue="island", multiple="stack")
plt.figure(1)
# funktioniert genauso, hat aber keine schöne box aussenrum
penguins_plot4 = sns.displot(data=penguins, x="bill_depth_mm", hue="island", multiple="stack", kind="kde")

# x = lable of the y axis line you want to use in the plot
# hue = lable of the y axis line you want to color in the plot
# data = the dataframe you want to use
plt.show()




# relative frequency plot
a = [-0.126,1,9,72.3,-44.2489,87.44]
bins = np.arange(-180,181,20)
hist, edges = np.histogram(a, bins)
freq = hist/float(hist.sum())
plt.bar(bins[:-1], freq, width=20, align="edge", ec="k" )
plt.show()
