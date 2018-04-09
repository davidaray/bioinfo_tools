#USAGE: RM2bed.py -i <input> -m <minsize> -s <sort criterion> -p <prefix(optional) -sp <split file by>
##Will process a .align.gz or .align file from RepeatMasker 
##to generate single or multiple sorted output 
##files in .bed format.
##VERY IMPORTANT CAVEAT - RepeatMasker output MUST have #Class/Family format.
##If the second term is missing, the output files will not have 

#To do
#1. convert to python
#2. add input options
#2a. detect align.gz vs .align
#2b. detect minimum size of hit
#2c. detect output type for sorting, i.e. by size, by name, by family, etc.

from Bio import SeqIO
import argparse
import os
import re
import gzip
import subprocess
import fileinput
import pandas as pd
import sys


####MAIN function
def main():

##Use the get_args function
	ALIGN, MINSIZE, PREFIX, CRITERION, SPLIT, HITS = get_args()
	print('Input file = ' + ALIGN)
	
##If no prefix, get prefix from filename using the first part of the filename.
	if PREFIX is None:
		PREFIX = re.split("[_.]", ALIGN)[0]
		print('No prefix provided. Will use ' + PREFIX)
	else:
		print('Prefix = ' + PREFIX)
	if CRITERION is None:
		pass
	else:
		print('Sort criterion is ' + CRITERION)

##Determine file type using extension and create a tmp.file to use.
	if ALIGN.lower().endswith('gz'):
		print('Using .gz file.')
		TMP = PREFIX + '.tmp.file'
		with gzip.open(ALIGN, 'rt') as IN:
			with open(TMP, 'w') as OUT:
				for LINE in IN:
					OUT.write(LINE)
	elif ALIGN.lower().endswith('align'):
		print('Using .align file')
		TMP = PREFIX + '.tmp.file'
		print('TMP = ' + TMP)
		IN = open(ALIGN, 'r')
		with open(TMP, 'w') as OUT:
			for LINE in IN:
				OUT.write(LINE)
		IN.close()
	else:
		print('Input file must be either .align.gz or .align.')
		exit()

##Call subprocess, calcDivergenceFromAlign.Ray.pl, using perl. Keep the new file for troubleshooting.
	print('Calculating divergences and saving to div.out file.')
	subprocess.check_call('perl /lustre/work/daray/software/RepeatMasker/util/calcDivergenceFromAlign.Ray.pl -noCpG {} >{}'.format(TMP, PREFIX + '_div.out'), shell=True)

##Delete tmp file.
	print('Removing tmp file.')
	os.remove(TMP)

##Grep to remove lines with simple and satellite repeats. 
##Also replace '/' and '#' with tabs and delete 'kimura='.
	print('Removing Simple and Satellite repeats.')
	print('Replacing hashtags and slashes with tabs.')
	print('Removing kimura=')
	BADWORDS = ['Simple','Satellite']
	DIV_OUT_LIST = []
	with open(PREFIX + '_div.out') as DIV_OUT_FILE:
		for LINE in DIV_OUT_FILE:
			LINE = re.sub(' +','\t', LINE)
			LINE = LINE.lstrip()
			LINE = LINE.replace('/', '\t').replace('kimura=', '').replace('#', '\t')		
			DIV_OUT_LIST.append(LINE)
		if not any(WORD in LINE for WORD in BADWORDS):
			DIV_OUT_LIST.append(LINE)

##Create a tmp file from the new list for use as an dataframe in pandas
##For some reason, there is an extra copy of the last line in the list. the 'for' line gets rid of it.
	NEWFILE = open(PREFIX + '.tmp.bed', 'w')
	for LINE in DIV_OUT_LIST[0:len(DIV_OUT_LIST)-1]:
#		print(LINE)
		NEWFILE.write(LINE)
	NEWFILE.close()

##Create a pandas dataframe from the tmp file, adding headers as you do so.
	print('Reading in dataframe.')
	OUT_ARRAY = pd.read_table(PREFIX + '.tmp.bed', sep='\t', names=['A', 'B', 'C', 'D', 'chrom', 'start', 'stop', 'E', 'strand', 'name', 'class', 'family', 'F', 'G', 'H', 'I', 'J', 'diverge'])
#	OUT_ARRAY.to_csv(PREFIX + '_import_tmp.bed', sep='\t', header=False, index=False)

##Delete the tmp file
#	os.remove(PREFIX + '.tmp.bed')
	
##Calculate size of insertion and add column to end of lines
	print('Adding size column.')
	OUT_ARRAY['size'] = OUT_ARRAY['stop'].subtract(OUT_ARRAY['start'] -1)
#	OUT_ARRAY.to_csv(PREFIX + '_addsize_tmp.bed', sep='\t', header=True, index=True)

##Rearrange the columns as you want
	print('Rearranging columns.')
	OUT_ARRAY = OUT_ARRAY[['chrom', 'start', 'stop', 'name', 'size', 'strand', 'class', 'family', 'diverge']]
#	OUT_ARRAY.to_csv(PREFIX + '_arrange_tmp.bed', sep='\t', header=True, index=True)

##Replace 'C' with '-' in the 'strand' column. 
	print('Replacing C with -.')
	OUT_ARRAY.strand.replace('C', '-', inplace=True)
#	OUT_ARRAY.to_csv(PREFIX + '_replaceC_tmp.bed', sep='\t', header=True, index=True)

##Filter by minsize
	print('Filtering by minsize.')
	OUT_ARRAY = OUT_ARRAY[OUT_ARRAY['size'] >= MINSIZE]
#	OUT_ARRAY.to_csv(PREFIX + '_minsize_tmp.bed', sep='\t', header=True, index=True)

##Sort main output if asked.	
	if CRITERION is not None:
		if CRITERION in ['name', 'family', 'class']:
			OUT_ARRAY = OUT_ARRAY.sort_values([CRITERION])
		elif CRITERION in ['size']:
			OUT_ARRAY = OUT_ARRAY.sort_values([CRITERION], ascending = [0])
		elif CRITERION in ['diverge']:
			OUT_ARRAY = OUT_ARRAY.sort_values([CRITERION], ascending = [1])
		else:
			print('Choices are size, name, family, class, or diverge. Not sorting.')	

##Write the dataframe to a file
	OUT_ARRAY.to_csv(PREFIX + '_rm.bed', sep='\t', header=False, index=False)

##Split into files if asked. Also check to see if there is a minumum hit number and act accordingly.
	if SPLIT is not None:
		print('Split files by ' + SPLIT)
		if SPLIT in ['name', 'family', 'class']:
			CLUSTERED = OUT_ARRAY.sort_values([SPLIT])
			for SPLITVALUE in CLUSTERED[SPLIT].unique():
				if HITS is not None:
					print('Will only output files with at least ' + HITS + ' hits.')
					CLUSTEREDW = CLUSTERED[CLUSTERED[SPLIT]==SPLITVALUE]
					COUNT_ROW = CLUSTEREDW.shape[0]
					if COUNT_ROW >= HITS:
						CLUSTEREDW.to_csv(PREFIX + '_' + SPLITVALUE + '_rm.bed', sep='\t', header=False, index=False)
				else:	
					CLUSTEREDW = CLUSTERED[CLUSTERED[SPLIT]==SPLITVALUE]
					CLUSTEREDW.to_csv(PREFIX + '_' + SPLITVALUE + '_rm.bed', sep='\t', header=False, index=False)
		else:
			print('Splitting options are by name, family, and class.')			

##Get arguments function
def get_args():
	parser = argparse.ArgumentParser(description="Will process a .align or .align.gz file from RepeatMasker to generate multiple sorted output files in .bed format", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', type=str, help='Name of the .align or .align.gz file to be parsed', required=True)
	parser.add_argument('-m', '--minsize', type=int, help='Minimum size of hit to include in sorted file', default = 100)
	parser.add_argument('-s', '--sortcriterion', type=str, help='Sort criterion, i.e. size, name, family, class, size, or divergence (diverge), etc.')
	parser.add_argument('-p', '--prefix', type=str, help='Prefix to use for output file - default is first field of input filename')
	parser.add_argument('-sp', '--split', type=str, help='Split into files based on name, family, class? This is optional.')
	parser.add_argument('-n', '--minhitnum', type=int, help='Minimum number of hits in file before being created. Only implemented if --split option is invoked. Optional.')

	args = parser.parse_args()
	ALIGN = args.input
	MINSIZE = args.minsize
	PREFIX = args.prefix
	CRITERION = args.sortcriterion
	SPLIT = args.split
	HITS = args.minhitnum

	return ALIGN, MINSIZE, PREFIX, CRITERION, SPLIT, HITS

if __name__ =="__main__":main()
		
#gzip $ABBREV".align" &

#gzip $ABBREV"_div.out" &

#gzip $ABBREV"_wCpG_div.out" &


