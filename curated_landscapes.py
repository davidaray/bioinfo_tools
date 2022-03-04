# RM-like Landscapes
# Nicole S. Paulat
# 1 February 2022
# nicole.paulat@ttu.edu

# Input: 
# bedfiles. One bed file per taxon named as TaxonAbbreviation_rm.bed
#     genomesizefile. Tab delimited with TaxonName, GenomeSize, MutationRate, and TaxonAbbreviation. 
#     For example, Ixodes_scapularis	2226883318	3.39917E-09	iSca. May include multiple 
#     lines for each taxon being evaluated.
#Output:
# CSV files of TE subfamilies of interest (columns) genome contribution (bp) by age (rows)
# PNG landscape plot of TE class accumulation by proxy age (divergence, x-axis) and genome 
#     proportion (y-axis)



#############################################
############## Import modules ###############
#############################################
import argparse
import os
import sys
import re
import math
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from cycler import cycler
plt.switch_backend('Agg')
from pylab import savefig
import seaborn as sns

################# Arguments ##################

def get_args():
    parser = argparse.ArgumentParser(description="Will generate pie chart and landscape plot along with multiple table files when given repeatmasker.bed files and a genomesize file.")
    parser.add_argument("-d", "--divergence", required=True, type = int, help="Please enter your maximum divergence value. Default = 50.", default = 50)
    parser.add_argument("-g", "--genomesize", required=True, type = str, help="Enter path to genomesizefile. Tab delimited with TaxonName, GenomeSize, MutationRate, and TaxonAbbreviation. For example, Ixodes_scapularis	2226883318	3.39917E-09	iSca. May include multiple lines for each taxon being evaluated.")
#    parser.add_argument("-m", "--minimum100bp", type = str, help="Only count hits that are at least 100 bp long? y or n? Default = y. Optional.", default = 'y')
    parser.add_argument("-m", "--minimum100bp", action='store_false', help="If entered, count hits that are at under 100 bp long. Otherwise, omit them. Default (not entered) is to omit.")

    args = parser.parse_args()
    DIVERGENCE = args.divergence
    SIZEFILE = args.genomesize
    MIN100 = args.minimum100bp

    if MIN100 is True:
        print('minimum100bp not entered. Default = y. Not counting hits < 100 bp.')
    else:
        print('minimum100bp flag entered. Counting hits < 100 bp.')
#    if MIN100 == 'y':
#        MIN100 == 'y'
#        print('minimum100bp entered as y. Not counting hits < 100 bp.')
#    elif MIN100 == 'n':
#        MIN100 == 'n'
#        print('minimum100bp input is n. Counting hits < 100 bp.')
#    else: 
#        MIN100 == 'y'
#        print('minimum100bp intput is not y or n. Value set to y. Not counting hits < 100 bp.')
    
    return DIVERGENCE, SIZEFILE, MIN100

################# Plotting ###################

def plot(SIZEFILE, DIVERGENCE, MIN100):
    print('Reading in genome sizes file, ' + SIZEFILE + '.')
    GENOMESIZES = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu', 'taxon_abbrev'], index_col=0)
    GENOMESIZES = GENOMESIZES.squeeze()
    GENOMESIZESFRAME = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu', 'taxon_abbrev'])
    GENOMESIZESFRAME = GENOMESIZESFRAME.squeeze()
    TAXA = GENOMESIZESFRAME['taxon'].tolist()
    TAXA.reverse()
    IDS = GENOMESIZESFRAME['taxon_abbrev'].tolist()
    IDS.reverse()
    TAXALENGTH = len(TAXA)
    print('Will plot divergences and accumulation for ' + str(TAXALENGTH) + ' species.')
    NAMES = [x.replace("_", " ") for x in TAXA]

    ALLTAXA = pd.DataFrame()
    ALLTAXA2 = pd.DataFrame()
    ALLTAXA3 = pd.DataFrame()

    TETYPES = ['SINE', 'LINE', 'LTR', 'RC', 'DNA', 'Unknown', 'Satellite', 'Simple_repeat']

    for TAXON, ID in zip(TAXA, IDS):
        # Read in the filename for a given taxon
        BEDFILE = ID + '_rm.bed'
        print('Reading in ' + BEDFILE)
        COLS = ['Scaffold', 'Start', 'End', 'TE', 'Size', 'Strand', 'Class', 'Family', 'Div', 'id']
        FILE = pd.read_csv(BEDFILE, sep="\t", names=COLS)
        
        GENSIZE = GENOMESIZES.genomesize[TAXON]

        # Filter to only keep elements that are more than 100 bp in length
        FILE.drop(['Scaffold', 'Start', 'End', 'Strand', 'id'], axis=1, inplace=True)
        FILE=FILE[FILE['Class'].isin(TETYPES)]
        FILE['Size'] = pd.to_numeric(FILE['Size'])
        FILE['Div'] = pd.to_numeric(FILE['Div'])
        if MIN100 == "y":
            FILE=FILE[FILE.Size >= 100]
        FILE['Proportion']= FILE['Size']/GENSIZE
        #FILE['Age'] = FILE['Age'].map(lambda x: round(x, -7))
        #AGE_DF = FILE.groupby(["Class", "Age"])[['Proportion']].sum().reset_index()
        FILE['Div'] = FILE['Div'].map(lambda x: round(x, 0))
        DIV_DF = FILE.groupby(["Class", "Div"])[['Proportion']].sum().reset_index()
        CLASS_DF = FILE.groupby(["Class"])[['Proportion']].sum().reset_index()
        # Filter to only keep TEs that occur at least 100 times
        #FILE = (FILE.groupby("TE").filter(lambda x : len(x)>100))
        # Create new dataframe with only taxon id and age
        #FILE = FILE[['Taxon', 'Div']]
        #FILE = FILE[['Taxon', 'Class', 'Age']]
        #print(AGE_DF.head())
        #AGE_DF['Taxon'] = NAME
        #DIV_DF['Taxon'] = NAME
        #CLASS_DF['Taxon'] = NAME
        DIV_DF['Taxon'] = ID
        CLASS_DF['Taxon'] = ID
        # Append to growing dataframe
        #ALLTAXA1 = pd.concat([ALLTAXA, AGE_DF], ignore_index=True)
        ALLTAXA2 = pd.concat([ALLTAXA2, DIV_DF], ignore_index=True)
        ALLTAXA3 = pd.concat([ALLTAXA3, CLASS_DF], ignore_index=True)
        # Optional: filter to just young elements (ex. 10my)
        #ALLTAXA1=ALLTAXA1[ALLTAXA1.Age <= AGE]
        ALLTAXA2=ALLTAXA2[ALLTAXA2.Div <= DIVERGENCE]

        # Save to file
        #ALLTAXA1.to_csv(TAXON + '_curated_repeats_' + str(AGE) + 'my_lineplotframe.txt', index=False)
        ALLTAXA2.to_csv(TAXON + '_curated_repeats_' + str(DIVERGENCE) + '_div_lineplotframe.txt', index=False)
        ALLTAXA3.to_csv(TAXON + '_curated_repeats_' + str(DIVERGENCE) + '_proportions_lineplotframe.txt', index=False)

        ##Make landscape plots per species
        sns.set_style("ticks")
        plt.rcParams['font.family'] = "sans-serif"
        plt.rcParams['text.usetex'] = False
        hue_colors = {'SINE': '#ffb3e6', 'LINE': 'mediumturquoise', 'LTR': 'lightgreen', 'RC': 'deeppink', 'DNA': '#ffcc99', 'Unknown':'gray', 'Satellite':'darkgray', 'Simple_repeat':'black'}
        g = sns.relplot(x = 'Div', y = 'Proportion', hue = 'Class', hue_order = TETYPES, palette = hue_colors, col = "Taxon", height = 3, kind = "line", estimator = None, data = ALLTAXA2)
        #g.fig.subplots_adjust(top = 0.8)
        #g.set_titles('${col_name}$')
        plt.title(label='Ixodes scapularis', fontstyle='italic')
        #x_ticks = np.arange(0, 150000001, 10000000)
        x_ticks = np.arange(0, 51, 10)
        plt.xticks(x_ticks)
        #plt.ticklabel_format(axis='x', style='sci', scilimits=(6,6), useMathText=True)
        #plt.xlim(0, 150000000)
        plt.xlim(0, DIVERGENCE)
        g.set_ylabels('Genome Proportion')
        #g.set_xlabels('Age (My)')
        g.set_xlabels('Divergence')
        g.savefig(TAXON + '_curated_rep_div' + str(DIVERGENCE) + '_landscape' + '.png', dpi=300)
        g.savefig(TAXON + '_curated_rep_div' + str(DIVERGENCE) + '_landscape' + '.svg')

        ##### Generate piechart of repeat genome proportions #####
        data = ALLTAXA3.copy()
        cols = 4
        d = data.groupby(["Taxon"])[["Proportion"]].sum().reset_index()
        d["Class"] = "Unmasked"
        d["Proportion"] = 1 - d["Proportion"]
        d = d[["Class", "Proportion", "Taxon"]]
        data = pd.concat([data, d], ignore_index=True)
        data = data.sort_values(by=['Taxon', 'Class'], ascending=False)
        ALLTAXA3.to_csv(TAXON + '_curated_repeats_proportions_lineplotframe.txt', index=False)
        p = data.groupby("Taxon")
        rows = int(np.ceil(len(p)/cols))

        LABELS = 'Unmasked', 'Unknown', 'DNA', 'RC', 'LTR', 'LINE', 'SINE', 'Satellite', 'Simple Repeat'
        COLORS = ['white', 'silver', '#ffcc99', 'deeppink', 'lightgreen', 'mediumturquoise', '#ffb3e6', 'grey', 'darkgrey']
        #LABELS = 'DNA', 'LINE', 'LTR', 'RC', 'SINE', 'Unknown', 'Unmasked'
        #COLORS = ['#ffcc99', 'mediumturquoise', 'lightgreen', '#ffb3e6', 'deeppink', 'silver', 'white']

        fig, axes = plt.subplots(ncols=cols, nrows=rows)
        for (c, grp), ax in zip(p, axes.flat):
            #ax.pie(grp.Proportion, labels=grp.Class)
            ax.pie(grp.Proportion, labels=grp.Class, colors=COLORS, autopct='%1.1f%%', shadow=True, startangle=140, textprops={'fontsize': 6})
            plt.axis('equal')
            ax.set_title(c, fontdict = {'fontsize' : 8})

        if len(p) < cols*rows:
            for ax in axes.flatten()[len(p):]:
                ax.axis("off")

        plt.tight_layout()

        fig.savefig(TAXON + '_curated_repeats_piechart' + '.png', dpi=300)
        fig.savefig(TAXON + '_curated_repeats_piechart' + '.svg')


########## MAIN function ###########
def main():	
##Get input arguments
    DIVERGENCE, SIZEFILE, MIN100 = get_args()

##Start the plotting
    plot(SIZEFILE, DIVERGENCE, MIN100)

if __name__ =="__main__":main()

