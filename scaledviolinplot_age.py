# Violin Plots
# Jenna R. Grimshaw
# Oct 11, 2019
# Modified by Jenny Korstian 4/3/2020
# Modified by David Ray 5/21/2021

import pandas as pd
import os
import glob
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
from pylab import savefig

# 1. Loop through all rm bed files
# 2. Concatenate all into a single file
# 3. Create violin plots

#######
## This script assumes you want to generate a violin plot from the 'class' files generated
## by cat_data_props.py followed by filtering via filter_beds.py.
#######

##### Housekeeping #####
## Decide if you want to do this on TE class or TE family
CATEGORY = 'class'
#CATEGORY = 'family'

## Update paths here if necessary
PROCESSEDBEDSPATH = '/lustre/scratch/daray/npaulat_beds/' + CATEGORY + '_files'
WORKPATH = '/lustre/scratch/daray/npaulat_beds/testing'

## Decide what class of TE you want to plot
TETYPE = 'DNA'
#TETYPE = 'RC'
## Decide whether you want 'all' rows of bed files or some subcategory filtered by filter_beds.py
#RANGE = 'all'
RANGE = '50my'

## Provide a file that has a list of the names of each species, one per line. This should be available as the size file provided to 'filter_beds.py'.
SIZEFILE = WORKPATH + '/genome_sizes_short.txt'

# Initialize a dataframe
ALLTAXA = pd.DataFrame()

# Read in the genomesizes to get a taxon list.
print('Reading in genome sizes file.')
GENOMESIZES = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu'], squeeze=True, index_col=0)
GENOMESIZESFRAME = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu'], squeeze=True)
TAXA = GENOMESIZESFRAME['taxon'].tolist()
TAXALENGTH = len(TAXA)
print('Will plot ages and accumulation for ' + str(TAXALENGTH) + ' species.')

## Check file names. As written this limits to only DNAs and 50my files
for TAXON in TAXA:
    # Read in the filename for a givan taxon
    print('Reading in ' + TAXON + '_' + TETYPE + '_' + RANGE + '_' + CATEGORY +'_processed_beds.txt')
    FILE = pd.read_csv(PROCESSEDBEDSPATH + '/' + TAXON + '_' + TETYPE + '_' + RANGE + '_' + CATEGORY +'_processed_beds.txt', sep="\t")
    # Rename the columns to make downstream processing a bit easier
    FILE = FILE.rename(columns={TAXON + '_TE' : 'TE', TAXON + '_size' : 'Size', TAXON + '_prop' : 'Prop', TAXON + '_class' : 'Class', TAXON + '_family' : 'Family', TAXON + '_div' : 'Div', TAXON + '_age' : 'Age'})
    # Add column with species ID
    FILE['Taxon'] = TAXON
    # Filter to only keep elements that are more than 100 bp in length
    FILE=FILE[FILE.Size >= 100]
    # Filter to only keep TEs that occur at least 100 times
    FILE = (FILE.groupby("TE").filter(lambda x : len(x)>100))	
    # Create new dataframe with only taxon id and age
    FILE = FILE[['Taxon', "Age"]]
    print(FILE.head())
    # Append to growing dataframe
    ALLTAXA = pd.concat([ALLTAXA, FILE], ignore_index=True)
ALLTAXA.to_csv(WORKPATH + '/violinframe.txt', index=False)
# Optional: filter to just young elements (ex. 10my)
ALLTAXA=ALLTAXA[ALLTAXA.Age <= 50000000]
FIG, ax = plt.subplots()
sns.violinplot('Taxon', 'Age', data=ALLTAXA, scale="count", ax=ax, order = TAXA)
ax.set_title(TETYPE + 's')
ax.yaxis.grid(True)
ax.set_xlabel('Taxon')
ax.set_ylabel('Age')
FIG.savefig(WORKPATH + '/scaled' + TETYPE + 's_violin' + '.png')
