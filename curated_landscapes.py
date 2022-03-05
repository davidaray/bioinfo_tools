# RM-like Landscapes
# Nicole S. Paulat and David Ray
# 5 February 2022
# nicole.paulat@ttu.edu, david.4.ray@gmail.com

# Command line input: 
# divergence - Maximum divergence to be plotted. Default = 50.
# bedfiles - One bed file per taxon named as TaxonAbbreviation_rm.bed
# genomesizefile - Tab delimited with TaxonName, GenomeSize, MutationRate, and TaxonAbbreviation. For 
#     example, Ixodes_scapularis	2226883318	3.39917E-09	iSca. May include multiple lines for 
#     each taxon being evaluated.
# minimum100 bp - flag (-m) to indicate if hits <100 bp should be included in calculations. Optional.
#     Default = False (do not include hits <100 bp).
#
# Other input:
# bedfiles - One bed file per taxon named as TaxonAbbreviation_rm.bed
#
# Output:
# CSV files of TE subfamilies of interest (columns) genome contribution (bp) by age (rows)
# PNG and SVG landscape plots of TE class accumulation by proxy age (divergence, x-axis) and genome 
#     proportion (y-axis).
# PNG and SVG pie charts of genome proportions of TE classes.

#Usage: python curated_landscapes.py -d 50 -g sizefile.txt [-m]


###### CURRENT ISSUES: Currently unable to properly add taxon names to combined landscape plots. 
# Lines 270-294. Unsure how to fix.


############## Import modules ###############
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
    parser.add_argument("-m", "--minimum100bp", action='store_false', help="If entered, count hits that are at under 100 bp long. Otherwise, omit them. Default (not entered) is to omit.")

    args = parser.parse_args()
    DIVERGENCE = args.divergence
    SIZEFILE = args.genomesize
    MIN100 = args.minimum100bp

    if MIN100 is True:
        print('minimum100bp not entered. Default = y. Not counting hits < 100 bp.')
    else:
        print('minimum100bp flag entered. Counting hits < 100 bp.')
    
    return DIVERGENCE, SIZEFILE, MIN100

################# Gather necessary data from input ###################

def generate_df(SIZEFILE, DIVERGENCE, MIN100):
    print('Reading in genome sizes file, ' + SIZEFILE + '.')
    # Read in size file and assign headers to columns
    GENOMESIZES = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu', 'taxon_abbrev'], index_col=0)
    # Read in size file and assign headers to columns
    GENOMESIZESFRAME = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu', 'taxon_abbrev'])
    #Create lists of taxa, IDs, genome sizes, and mutation rates.
    TAXA = GENOMESIZESFRAME['taxon'].tolist()
    TAXA.reverse()
    IDS = GENOMESIZESFRAME['taxon_abbrev'].tolist()
    IDS.reverse()
    GSIZES = GENOMESIZESFRAME['genomesize'].tolist()
    GSIZES.reverse()
    MU = GENOMESIZESFRAME['mu'].tolist()
    MU.reverse()

    #Count how many taxa to process.
    TAXALENGTH = len(TAXA)
    print('Will plot divergences and accumulation for ' + str(TAXALENGTH) + ' species.')
#    NAMES = [x.replace("_", " ") for x in TAXA]

    #Initiate dataframes for downstream
    ALLTAXADIV = pd.DataFrame()
    ALLTAXACLASS = pd.DataFrame()

    #Define TE types to consider
    TETYPES = ['SINE', 'LINE', 'LTR', 'RC', 'DNA', 'Unknown', 'Satellite', 'Simple_repeat']

    # For each group of paired species names and IDs
    for TAXON, ID in zip(TAXA, IDS):
        # Read in the filename for a given taxon
        BEDFILE = ID + '_rm.bed'
        print('Reading in ' + BEDFILE)
        print('TAXON = ' + TAXON)
        print('ID = ' + ID)
#       Assign column names and read in the file
        COLS = ['Scaffold', 'Start', 'End', 'TE', 'Size', 'Strand', 'Class', 'Family', 'Div', 'id']
        FILE = pd.read_csv(BEDFILE, sep="\t", names=COLS)
        
        #Get genome size from the genomesizes file
        GENSIZE = GENOMESIZES.genomesize[TAXON]
        print(TAXON + ' genome size = ' + str(GENSIZE) + ' bp.')

        #Drop unnecessary columns in dataframe.
#        print("printing original dataframe")
#        print(FILE)
        FILE.drop(['Scaffold', 'Start', 'End', 'Strand', 'id'], axis=1, inplace=True)
#        print("printing dropped dataframe")
#        print(FILE)
        #Only allow the Classes in 'TETYPES'
        FILE=FILE[FILE['Class'].isin(TETYPES)]
#        print("printing Class dataframe")
#        print(FILE)
        #Convert values in Size column to numeric
        FILE['Size'] = pd.to_numeric(FILE['Size'])
#        print('printing sizes dataframe')
#        print(FILE)
        #Convert values in Div column to numeric
        FILE['Div'] = pd.to_numeric(FILE['Div'])
#        print('printing numeric div dataframe')
#        print(FILE)
        #Eliminate hits < 100 bp if minimum100bp flag was set to y
        if MIN100 is True:
            FILE=FILE[FILE.Size >= 100]
        #Add Proportion column
        print('Proportion = each value in Size column dividied by ' + str(GENSIZE)) 
        FILE['Proportion']= FILE['Size']/GENSIZE
#        print('printing proportion dataframe')
#        print(FILE)
        #Round divergence column (why am I doing this?)
        FILE['Div'] = FILE['Div'].map(lambda x: round(x, 0))
#        print('printing divergence dataframe')
#        print(FILE)
        #Create Divergence dataframe from FILE. Has proportions by divergence bin.
        DIV_DF = FILE.groupby(["Class", "Div"])[['Proportion']].sum().reset_index()
#        print('printing DIV_DF')
#        print(DIV_DF)
        #Create Class dataframe from FILE. Has proportions by Class.
        CLASS_DF = FILE.groupby(["Class"])[['Proportion']].sum().reset_index()
#        print('printing CLASS_DF')
#        print(CLASS_DF)
        #add TAXON column to DIV_DF
        DIV_DF['Taxon'] = TAXON
        DIV_DF.to_csv(TAXON + '_curated_repeats_' + str(DIVERGENCE) + '_div_lineplotframe.txt', sep='\t', index=False)
#        print('printing DIV_DF')
#        print(DIV_DF)
        #add TAXON column to CLASS_DF
        CLASS_DF['Taxon'] = TAXON
        CLASS_DF.to_csv(TAXON + '_curated_repeats_' + str(DIVERGENCE) + '_proportionsframe.txt', sep='\t', index=False)
#        print('printing CLASS_DF')
#        print(CLASS_DF)
        # Append to growing dataframe
        #ALLTAXA1 = pd.concat([ALLTAXA, AGE_DF], ignore_index=True)
        ALLTAXADIV = pd.concat([ALLTAXADIV, DIV_DF], ignore_index=True)
#        print('printing ALLTAXADIV')
#        print(ALLTAXADIV)
        ALLTAXACLASS = pd.concat([ALLTAXACLASS, CLASS_DF], ignore_index=True)
#        print('printing ALLTAXACLASS')
#        print(ALLTAXACLASS)
        #Filter to just elements less than the maximum divergence
        ALLTAXADIV=ALLTAXADIV[ALLTAXADIV.Div <= DIVERGENCE]
        
#    print('printing ALLTAXADIV')
#    print(ALLTAXADIV)
#    print('printing ALLTAXACLASS')
#    print(ALLTAXACLASS)
        
    # Save all dataframes to appropriate files
    ALLTAXADIV.to_csv('curated_repeats_' + str(DIVERGENCE) + '_div_lineplotframe.txt', sep='\t', index=False)
    ALLTAXACLASS.to_csv('curated_repeats_' + str(DIVERGENCE) + '_proportionsframe.txt', sep='\t', index=False)
    
################# Generate the various plots ###################      
def plot(SIZEFILE, DIVERGENCE):
    #Read in all the same information as in the generate_df function
    GENOMESIZES = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu', 'taxon_abbrev'], index_col=0)
    GENOMESIZESFRAME = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu', 'taxon_abbrev'])
    TAXA = GENOMESIZESFRAME['taxon'].tolist()
    TAXA.reverse()
    IDS = GENOMESIZESFRAME['taxon_abbrev'].tolist()
    IDS.reverse()
    GSIZES = GENOMESIZESFRAME['genomesize'].tolist()
    GSIZES.reverse()
    MU = GENOMESIZESFRAME['mu'].tolist()
    MU.reverse()
    TAXALENGTH = len(TAXA)
    TETYPES = ['SINE', 'LINE', 'LTR', 'RC', 'DNA', 'Unknown', 'Satellite', 'Simple_repeat']
    
    #For every taxon, read in the saved files as new dataframes
    for TAXON, ID in zip(TAXA, IDS):
        DIVFILE = TAXON + '_curated_repeats_' + str(DIVERGENCE) + '_div_lineplotframe.txt'
        TAXADIV = pd.read_csv(DIVFILE, sep='\t', header=0)
        TAXADIV['Taxon'] = TAXADIV['Taxon'].str.replace('_', ' ')
        CLASSFILE = TAXON + '_curated_repeats_' + str(DIVERGENCE) + '_proportionsframe.txt'
        TAXACLASS = pd.read_csv(CLASSFILE, sep='\t', header=0)
        TAXACLASS['Taxon'] = TAXACLASS['Taxon'].str.replace('_', ' ')
#        print(TAXADIV)
#        print(TAXACLASS)

        ##Make landscape plots per species
        #Set the parameters that control the general style of the plots
        sns.set_style("ticks")
        #Each time Matplotlib loads, it defines a runtime configuration (rc) containing the default styles for every plot element you create. We'll start by saving a copy of the current rcParams dictionary, so we can easily reset these changes in the current session
        plt.rcParams['font.family'] = "sans-serif"
        plt.rcParams['text.usetex'] = False
        #Set the colors for each TE type
        hue_colors = {'SINE': '#ffb3e6', 'LINE': 'mediumturquoise', 'LTR': 'lightgreen', 'RC': 'deeppink', 'DNA': '#ffcc99', 'Unknown':'gray', 'Satellite':'darkgray', 'Simple_repeat':'black'}
        #The sns.relplot function provides access to several different axes-level functions that show the relationship between two variables with semantic mappings of subsets.
        g = sns.relplot(x = 'Div', y = 'Proportion', hue = 'Class', hue_order = TETYPES, palette = hue_colors, col = "Taxon", height = 3, kind = "line", estimator = None, data = TAXADIV)
        #Remove the _ from the taxon name before adding the title.
        TAXONNAME = TAXON.replace("_", " ")
        #Set the title for the individual plots. 
        plt.title(label=TAXONNAME, fontstyle='italic')
        #np.arange = Return evenly spaced values within a given interval.
        x_ticks = np.arange(0, 51, 10)
        #Get or set the current tick locations and labels of the x-axis
        plt.xticks(x_ticks)
        #Get or set the x limits of the current axesplt.xlim(0, DIVERGENCE)
        #Set the y label
        g.set_ylabels('Genome Proportion')
        #Set the x label
        g.set_xlabels('Divergence')
        #Save the files as png and svg.
        g.savefig(TAXON + '_curated_rep_div' + str(DIVERGENCE) + '_landscape' + '.png', dpi=300)
        g.savefig(TAXON + '_curated_rep_div' + str(DIVERGENCE) + '_landscape' + '.svg')
        
        ##Generate piechart of repeat genome proportions per species
        #Import the data from TAXACLASS dataframe
        DATA = TAXACLASS.copy()
        #There are four columns in the dataframe
        COLS = 4
        #Group DataFrame using a mapper or by a Series of columns. In this case, group each Taxon and assign the sum of the proportions for the Taxon to each group. Essentially, summarize the data by taxon.
        D = DATA.groupby(["Taxon"])[["Proportion"]].sum().reset_index()
        #Add a new class called Unmasked calculated by subtracting the total of all proportions from 1.
        D["Class"] = "Unmasked"
        D["Proportion"] = 1 - D["Proportion"]
        D = D[["Class", "Proportion", "Taxon"]]
        #Concatenate the Unmasked dataframe with the existing dataframe
        DATA = pd.concat([DATA, D], ignore_index=True)
        #Sort the dataframe from largest to smallest proportion
        DATA = DATA.sort_values(by=['Taxon', 'Class'], ascending=False)
        #Save the dataframe to a file.
        TAXACLASS.to_csv(TAXON + '_curated_repeats_proportions_lineplotframe.txt', index=False)
        
        #I'm not quite sure what this next bit does
        P = DATA.groupby("Taxon")        
        #Determine ROWS by dividing the number of species by the number of columns in TAXACLASS 
#        print('length of P = ' + str(len(P)))
        ROWS = int(np.ceil(len(P)/COLS))
#        print('ROWS = ' + str(ROWS))

        #Set labels and colors
        LABELS = 'Unmasked', 'Unknown', 'DNA', 'RC', 'LTR', 'LINE', 'SINE', 'Satellite', 'Simple Repeat'
        COLORS = ['white', 'silver', '#ffcc99', 'deeppink', 'lightgreen', 'mediumturquoise', '#ffb3e6', 'grey', 'darkgrey']

        fig, axes = plt.subplots(ncols=COLS, nrows=ROWS)
        for (c, grp), ax in zip(P, axes.flat):
            ax.pie(grp.Proportion, labels=grp.Class, colors=COLORS, autopct='%1.1f%%', shadow=True, startangle=140, textprops={'fontsize': 6})
            plt.axis('equal')
            ax.set_title(c, fontdict = {'fontsize' : 8})

        if len(P) < COLS*ROWS:
            for ax in axes.flatten()[len(P):]:
                ax.axis("off")

        plt.tight_layout()

        fig.savefig(TAXON + '_curated_repeats_piechart' + '.png', dpi=300)
        fig.savefig(TAXON + '_curated_repeats_piechart' + '.svg')

        
    ##Make landscape plots all species
    ALLDIVFILE = 'curated_repeats_' + str(DIVERGENCE) + '_div_lineplotframe.txt'
    ALLTAXADIV = pd.read_csv(ALLDIVFILE, sep='\t', header=0)
    ALLTAXADIV['Taxon'] = ALLTAXADIV['Taxon'].str.replace('_', ' ')
    sns.set_style("ticks")
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['text.usetex'] = False
    hue_colors = {'SINE': '#ffb3e6', 'LINE': 'mediumturquoise', 'LTR': 'lightgreen', 'RC': 'deeppink', 'DNA': '#ffcc99', 'Unknown':'gray', 'Satellite':'darkgray', 'Simple_repeat':'black'}
    g = sns.relplot(x = 'Div', y = 'Proportion', hue = 'Class', hue_order = TETYPES, palette = hue_colors, col = "Taxon", height = 3, kind = "line", estimator = None, data = ALLTAXADIV)
    #### See https://datavizpyr.com/customize-titles-in-multi-panel-plots-with-seaborn/
    #g.fig.subplots_adjust(top = 0.8)
    #g.set_titles('${col_name}$')
#    TAXONNAME = TAXON.replace("_", " ")
    plt.title(label=TAXON, fontstyle='italic')
    #x_ticks = np.arange(0, 150000001, 10000000)
    x_ticks = np.arange(0, 51, 10)
    plt.xticks(x_ticks)
    #plt.ticklabel_format(axis='x', style='sci', scilimits=(6,6), useMathText=True)
    #plt.xlim(0, 150000000)
    plt.xlim(0, DIVERGENCE)
    g.set_ylabels('Genome Proportion')
    #g.set_xlabels('Age (My)')
    g.set_xlabels('Divergence')
    g.savefig('curated_rep_div' + str(DIVERGENCE) + '_landscape' + '.png', dpi=300)
    g.savefig('curated_rep_div' + str(DIVERGENCE) + '_landscape' + '.svg')

    ##### Generate piechart of repeat genome proportions #####
    ALLCLASSFILE = 'curated_repeats_' + str(DIVERGENCE) + '_proportionsframe.txt'
    ALLTAXACLASS = pd.read_csv(ALLCLASSFILE, sep='\t', header=0)
    ALLTAXACLASS['Taxon'] = ALLTAXACLASS['Taxon'].str.replace('_', ' ')
    DATA = ALLTAXACLASS.copy()
    COLS = 4
    D = DATA.groupby(["Taxon"])[["Proportion"]].sum().reset_index()
    D["Class"] = "Unmasked"
    D["Proportion"] = 1 - D["Proportion"]
    D = D[["Class", "Proportion", "Taxon"]]
    DATA = pd.concat([DATA, D], ignore_index=True)
    DATA = DATA.sort_values(by=['Taxon', 'Class'], ascending=False)
    ALLTAXACLASS.to_csv('curated_repeats_proportions_lineplotframe.txt', index=False)
    P = DATA.groupby("Taxon")
    ROWS = int(np.ceil(len(P)/COLS))

    LABELS = 'Unmasked', 'Unknown', 'DNA', 'RC', 'LTR', 'LINE', 'SINE', 'Satellite', 'Simple Repeat'
    COLORS = ['white', 'silver', '#ffcc99', 'deeppink', 'lightgreen', 'mediumturquoise', '#ffb3e6', 'grey', 'darkgrey']
    #LABELS = 'DNA', 'LINE', 'LTR', 'RC', 'SINE', 'Unknown', 'Unmasked'
    #COLORS = ['#ffcc99', 'mediumturquoise', 'lightgreen', '#ffb3e6', 'deeppink', 'silver', 'white']

    fig, axes = plt.subplots(ncols=COLS, nrows=ROWS)
    for (c, grp), ax in zip(P, axes.flat):
        #ax.pie(grp.Proportion, labels=grp.Class)
        ax.pie(grp.Proportion, labels=grp.Class, colors=COLORS, autopct='%1.1f%%', shadow=True, startangle=140, textprops={'fontsize': 6})
        plt.axis('equal')
        ax.set_title(c, fontdict = {'fontsize' : 8})

    if len(P) < COLS*ROWS:
        for ax in axes.flatten()[len(P):]:
            ax.axis("off")

    plt.tight_layout()

    fig.savefig('curated_repeats_piechart' + '.png', dpi=300)
    fig.savefig('curated_repeats_piechart' + '.svg')


########## MAIN function ###########
def main():	
##Get input arguments
    DIVERGENCE, SIZEFILE, MIN100 = get_args()

##Generate data from bedfiles
    generate_df(SIZEFILE, DIVERGENCE, MIN100)
    
##Generate the plots
    plot(SIZEFILE, DIVERGENCE)
    
if __name__ =="__main__":main()

