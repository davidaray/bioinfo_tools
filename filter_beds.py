import pandas as pd
import os
import itertools
import argparse

####MAIN function
def main():

##Use the get_args function
    SIZEFILE, PREFIX, DIV, AGE = get_args()
    print('Input file = ' + SIZEFILE)

#Check for optional flags
    FLAGS = [PREFIX, AGE, DIV]
    if not FLAGS:
        print('No filtering based on age or divergence. Prefix is "all".\n')
        PREFIX = 'all'
    else:
        if not PREFIX:
            print('--prefix flag not used. Prefix is "all".\n')
            PREFIX = 'all'
        else:
            print('Custom filtering used.\n')
        if DIV:
            print('--divergence flag used.\n')
            print('Maximum divergence is ' + str(DIV) + '.\n')
        if AGE:
            print('--age flag used.\n')
            print('Maximum age is ' + str(AGE) + '.\n')

#Opens genomesizes.txt as a dictionary. Values taken from RepeatMasker output files 
    print('Reading in genome sizes file.')
    GENOMESIZES = pd.read_table('test_sizes_mrates.txt', sep='\t', names=['taxon', 'genomesize', 'mu'], squeeze=True, index_col=0)
    GENOMESIZESFRAME = pd.read_table('test_sizes_mrates.txt', sep='\t', names=['taxon', 'genomesize', 'mu'], squeeze=True)

#Get taxon list from genome sizes file.
    TAXA = GENOMESIZESFRAME['taxon'].tolist()

    TAXALENGTH = len(TAXA)

    for TAXON in TAXA:
        print('Reading in rm.bed for ' + TAXON + '. Writing out ' + TAXON + '_' + PREFIX +'_processed.bed.')
        TAXONMU = GENOMESIZES.at[TAXON, 'mu']
        print('Mutation rate for ' + TAXON + ' = '+ str(TAXONMU))
        TAXONRMBED = pd.read_table(TAXON + '_rm.bed', sep='\t', index_col=False, names=[TAXON + '_chromosome', 'start', 'stop', 'TE', 'size', 'orientation','class', 'family', 'div', 'ID'])
        TAXONPROCESSEDBED = TAXONRMBED[['TE', 'size', 'class', 'family', 'div']].copy()
        TAXONPROCESSEDBED['age'] = TAXONPROCESSEDBED['div'].div(TAXONMU)
        TAXONPROCESSEDBED['age'] = TAXONPROCESSEDBED['age'].div(100)
        TAXONPROCESSEDBED = TAXONPROCESSEDBED.round({'age': 1})
        if AGE:
            TAXONPROCESSEDBED = TAXONPROCESSEDBED[TAXONPROCESSEDBED['age'] <= AGE]
            TAXONPROCESSEDBED.to_csv(TAXON + '_' + PREFIX +'_filtered.bed', sep='\t', index=False)
        if DIV:
            TAXONPROCESSEDBED = TAXONPROCESSEDBED[TAXONPROCESSEDBED['div'] <= DIV]
            TAXONPROCESSEDBED.to_csv(TAXON + '_' + PREFIX +'_filtered.bed', sep='\t', index=False)
        else:
            TAXONPROCESSEDBED.to_csv(TAXON + '_' + PREFIX +'_unfiltered.bed', sep='\t', index=False)

##Get arguments function
def get_args():
	parser = argparse.ArgumentParser(description="If provided with a genomesize file consisting of three columns: genome id, genome size (bp), and mutation rate, this script will filter an rm.bed file output from RM2bed.py based on age or divergence. If no age or divergence filtering is proposed, it will simply add an age column.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#	parser.add_argument('-i', '--input', type=str, help='Name of the rm.bed file to be parsed.', required=True)
	parser.add_argument('-g', '--sizefile', type=str, help='File containing two corresponding columns of taxon abbreviations and genome sizes in bp.', required=True)
	parser.add_argument('-d', '--divergence', type=int, help='Maximum divergence allowed. 0-100. Optional')
	parser.add_argument('-p', '--prefix', type=str, help='Prefix to put after taxon id. I use this when separating young and old elements. Should be descriptive. For example if filtering to only 50my or younger, use "50my". Input name is aJam_rm.bed, output name is aJam_50my_rm.bed.')
	parser.add_argument('-a', '--age', type=int, help='Maximum age allowed. Requires file with muation rate values. This is optional.')
#	parser.add_argument('-n', '--minhitnum', type=int, help='Minimum number of hits in file before being created. Only implemented if --split option is invoked. Optional.')
#	parser.add_argument('-d', '--maxdiverge', type=float, help='Maximum divergence allowed in output file.')
#	parser.add_argument('-dmin', '--mindiverge', type=float, help='Minimum divergence allowed in output file.')

	args = parser.parse_args()
#	INPUT = args.input
	SIZEFILE = args.sizefile
	PREFIX = args.prefix
	DIV = args.divergence
	AGE = args.age
#	HITS = args.minhitnum
#	MAX = args.maxdiverge
#	MINDIV = args.mindiverge

	return SIZEFILE, PREFIX, DIV, AGE

if __name__ =="__main__":main()
	
