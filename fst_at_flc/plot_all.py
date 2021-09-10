#!/usr/bin/python3

import pandas as pd
import sys, os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Need an input file
if len(sys.argv) < 2:
    sys.exit("Need table as input")

# -----------------------------------------------------------------------------
#     Input
# -----------------------------------------------------------------------------

# Read the file
df = pd.read_csv( sys.argv[1], sep='\t' )
colnames = list(df.columns)
print("Found", len(df), "rows and", len(colnames), "cols")

# Prepare our needed columns
if not colnames[0] == "Chr" or not colnames[1] == "Pos":
    sys.exit("Expect columns Chr and Pos")
colnames.remove("Chr")
colnames.remove("Pos")
print("Columns:", colnames)

# Output file name based on input file
fn = os.path.splitext(sys.argv[1])[0]
# fn = os.path.splitext(os.path.basename(sys.argv[1]))[0]

# -----------------------------------------------------------------------------
#     Plot Properties
# -----------------------------------------------------------------------------

# Number of chromosomes to process.
# We expect the column "Chr" to contain these as strings/integers
num_chrom = 1

# Prepare plot with subplot for each chromosome.
fig, axes = plt.subplots(1,num_chrom, figsize=(num_chrom*8, 8))
plt.margins(0,0)

# If we have min max from the cmd line, use it
# Get min max if available from command line
if len(sys.argv) > 2:
    miny = float(sys.argv[2])
    maxy = float(sys.argv[3])
    print("Using miny", miny, "maxy", maxy)
else:
    miny = min([ np.nanmin(df[c].values) for c in colnames])
    maxy = max([ np.nanmax(df[c].values) for c in colnames])
    diff = maxy-miny
    miny = miny-diff/10.0
    maxy = maxy+diff/10.0
    print("Using miny", round(miny, 5), "maxy", round(maxy, 5))
plt.setp(axes, ylim=(miny, maxy))

# Remove all but the first y-axis labels and ticks
for ax in range(1,num_chrom):
    axes[ax].set_yticks([])
#plt.tight_layout()

# Remove all x margins within the subplot boxes
#for ax in range(0,num_chrom):
#    axes[ax].margins(x=0)

# -----------------------------------------------------------------------------
#     Colors and Lines
# -----------------------------------------------------------------------------

# Solid lines in separate colors for each replicate (hard coded for this dataset).
# We assume that the column name order is the one that we want!
# Then let's give each replicate its own color.
if len(colnames) == 12:
    # Use nice colors, with three consecutive samples (replicates) in the same color
    color_selection = [ "C0", "C1", "C2", "C4" ]
    pal = {colnames[key]: color_selection[int(key/3)] for key in range(12)}
    # pal = {
    #     "S1": "C0", "S2": "C0", "S3": "C0",
    #     "S4": "C1", "S5": "C1", "S6": "C1",
    #     "S7": "C2", "S8": "C2", "S9": "C2",
    #     "S10": "C4", "S11": "C4", "S12": "C4"
    # }
elif len(colnames) == 4:
    color_selection = [ "C0", "C1", "C2", "C3" ]
    pal = {colnames[key]: color_selection[key] for key in range(4)}
	# pal = {
	# 	"T0": "C0",
	# 	"T1": "C1",
	# 	"T2": "C2",
	# 	"T3": "C3"
	# }
elif len(colnames) == 3:
    color_selection = [ "C0", "C1", "C2" ]
    pal = {colnames[key]: color_selection[key] for key in range(3)}
    # pal = {
    #     "R1": "C0",
    #     "R2": "C1",
    #     "R3": "C2"
    # }
else:
    sys.exit("Invalid number of columns for this dataset")

linestyles = [""] * len(colnames)
linealpha=0.8

# -----------------------------------------------------------------------------
#     Plot the Data
# -----------------------------------------------------------------------------

# Hard coded chromosome list for now
#for chrom in range(1, num_chrom+1):
for chrom in range(5, 6):
    print("Chrom", str(chrom))

    # Select all data for that chromosome
    chr_data = df.loc[df['Chr'] == chrom][["Pos"] + colnames]
    data_melt = pd.melt(chr_data, ["Pos"])
    #data_melt = pd.melt(chr_data, ["Pos"]).head()
    #print(data_melt)

    # Make the plot in the subfigure of the current chromosome
    sns_plot = sns.scatterplot(
        x="Pos", y="value", hue="variable", data=data_melt, 
        palette=pal, alpha=linealpha
    )
    
axes.axvline(3173382)
axes.axvline(3179448) 
    
# -----------------------------------------------------------------------------
#     Save the Figure
# -----------------------------------------------------------------------------

# Remove all but the last legend. Have to do this at the end, because matplotlib
# otherwise has not created them yet.
for ax in range(0,num_chrom-1):
    axes[ax].legend().remove()

# Apparently, we cannot change the legend alpha and legend position at the same time,
# as changing one will reset the other... what a shitshow
# leg=plt.legend()
# for l in leg.get_lines():
#     l.set_alpha(linealpha)

# Move the last remaining legend out of the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# Save and finish the plot
#plt.show()
plt.savefig(fn + ".png", bbox_inches = 'tight')
plt.savefig(fn + ".svg", bbox_inches = 'tight')
plt.close()

# -----------------------------------------------------------------------------
#     Scatter Plots
# -----------------------------------------------------------------------------

if len(colnames) != 3:
    sys.exit()

# Prepare plot with subplot for each chromosome.
fig, axes = plt.subplots(1,3, figsize=(3*8, 8))
plt.margins(0,0)
plt.setp(axes, xlim=(miny, maxy), ylim=(miny, maxy))

def do_scatplot( col1, col2, figidx ):
    # compute correlation
    dfcorr = df[[ colnames[col1], colnames[col2] ]]
    corr = dfcorr.corr().iloc[0][1]
    corrtext = "r  = " + str(round(corr,3)) + "\nr2 = " + str(round(corr*corr,3))

    sns_plot = sns.scatterplot(
        x=colnames[col1], y=colnames[col2], data=df, ax=axes[figidx], color="#666666", alpha=0.6, marker='o', s=0.7
    )
    
    sns_plot.text(0.80, 0.02, corrtext, horizontalalignment='left', size='medium', color='black', weight='semibold')
    sns_plot.plot([miny, maxy], [miny, maxy], linewidth=1)

do_scatplot( 0, 1, 0 )
do_scatplot( 0, 2, 1 )
do_scatplot( 1, 2, 2 )

#plt.show()
plt.savefig(fn + "_scatter.png", bbox_inches = 'tight')
plt.close()
