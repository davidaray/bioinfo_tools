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
	ALIGN, MINSIZE, PREFIX, AGE, CRITERION, SPLIT, HITS, MAX, MINDIV, RERUN = get_args()
	print('Input file = ' + ALIGN)
	
##If no prefix, get prefix from filename using the first part of the filename.
	if PREFIX is None and AGE is None:
		CUTPREFIX = re.split("[_.]", ALIGN)[0]
		print('No prefix provided. No age provided. Will use ' + CUTPREFIX + ' as prefix.')
	else:
		CUTPREFIX = re.split("[_.]", ALIGN)[0]
		print('Age is ' + AGE + '.')
		PREFIXAGE = CUTPREFIX + '_' + AGE
		print ('Prefix = ' + PREFIXAGE)
	if CRITERION is None:
		pass
	else:
		print('Sort criterion is ' + CRITERION)

##If this is a rerun, skip some of these steps
##Determine file type using extension and create a tmp.file to use.
	if RERUN == 'n':
		if ALIGN.lower().endswith('gz'):
			print('Using .gz file.')
			TMP = CUTPREFIX + '.tmp.file'
			with gzip.open(ALIGN, 'rt') as IN:
				with open(TMP, 'w') as OUT:
					for LINE in IN:
						OUT.write(LINE)
		elif ALIGN.lower().endswith('align'):
			print('Using .align file')
			TMP = CUTPREFIX + '.tmp.file'
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
		print('Calculating divergences and saving to tmp_div.out file.')
		subprocess.check_call('perl /lustre/work/daray/software/RepeatMasker/util/calcDivergenceFromAlign.Ray.pl -noCpG {} >{}'.format(TMP, CUTPREFIX + '_tmp_div.out'), shell=True)

##Delete tmp file.
		print('Removing tmp file.')
		os.remove(TMP)
		
##Simple repeats are a problem with .out files. They sometimes lack the Class/Family structure. Fix this here.
		INFILE = CUTPREFIX + '_tmp_div.out'
		OUTFILE1 = open(CUTPREFIX + '_div1.out', 'w')
		INFILE1 = CUTPREFIX + '_div1.out'
		OUTFILE2 = open(CUTPREFIX + '_div2.out', 'w')
		INFILE2 = CUTPREFIX + '_div2.out'
		OUTFILE3 = open(CUTPREFIX + '_div3.out', 'w')
		INFILE3 = CUTPREFIX + '_div3.out'
		OUTFILE4 = open(CUTPREFIX + '_div4.out', 'w')
		INFILE4 = CUTPREFIX + '_div4.out'
		OUTFILE5 = open(CUTPREFIX + '_div5.out', 'w')
		INFILE5 = CUTPREFIX + '_div5.out'
		OUTFILE6 = open(CUTPREFIX + '_div.out', 'w')
		subprocess.call(['sed', 's/#Simple_repeat/#Simple_repeat\/Simple_repeat/g', INFILE], stdout=OUTFILE1)
		subprocess.call(['sed', 's/#Simple_repeat\/Simple_repeat\/Simple_repeat/#Simple_repeat\/Simple_repeat/g', INFILE1], stdout=OUTFILE2)
		subprocess.call(['sed', 's/#Simple_repeat\/Simple_repeat\/Satellite/#Simple_repeat\/Satellite/g', INFILE2], stdout=OUTFILE3) 
		subprocess.call(['sed', 's/#Simple_repeat\/Simple_repeat\/Low_complexity/#Simple_repeat\/Low_complexity/g', INFILE3], stdout=OUTFILE4) 
		subprocess.call(['sed', 's/#rRNA/#rRNA\/rRNA/g', INFILE4], stdout=OUTFILE5) 
		subprocess.call(['sed', 's/#rRNA\/rRNA\/rRNA/#rRNA\/rRNA/g', INFILE5], stdout=OUTFILE6) 
		OUTFILE1.close()
		OUTFILE2.close()
		OUTFILE3.close()
		OUTFILE4.close()
		OUTFILE5.close()
		os.remove(CUTPREFIX + '_div1.out')
		os.remove(CUTPREFIX + '_div2.out')
		os.remove(CUTPREFIX + '_div3.out')
		os.remove(CUTPREFIX + '_div4.out')
		os.remove(CUTPREFIX + '_div5.out')
	else:
		pass

##Also replace '/' and '#' with tabs and delete 'kimura='.
	print('Replacing hashtags and slashes with tabs.')
	print('Removing kimura=')
	DIV_OUT_LIST = []
	with open(CUTPREFIX + '_div.out') as DIV_OUT_FILE:
		for LINE in DIV_OUT_FILE:
			LINE = re.sub(' +','\t', LINE)
			LINE = LINE.lstrip()
			LINE = LINE.replace('/', '\t').replace('kimura=', '').replace('#', '\t')		
			DIV_OUT_LIST.append(LINE)

##Create a tmp file from the new list for use as an dataframe in pandas
##For some reason, there is an extra copy of the last line in the list. the 'for' line gets rid of it.
	NEWFILE = open(CUTPREFIX + '.tmp.bed', 'w')
	for LINE in DIV_OUT_LIST[0:len(DIV_OUT_LIST)-1]:
#		print(LINE)
		NEWFILE.write(LINE)
	NEWFILE.close()

##Create a pandas dataframe from the tmp file, adding headers as you do so.
	print('Reading in dataframe.')
	OUT_ARRAY = pd.read_table(CUTPREFIX + '.tmp.bed', sep='\t', names=['A', 'B', 'C', 'D', 'chrom', 'start', 'stop', 'E', 'strand', 'name', 'class', 'family', 'F', 'G', 'H', 'I', 'J', 'diverge'])
#	OUT_ARRAY.to_csv(PREFIX + '_import_tmp.bed', sep='\t', header=False, index=False)

##Delete the tmp file
	os.remove(CUTPREFIX + '.tmp.bed')
	
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

##Announce max hits limit			
	if HITS is not None:
		print('Will only output files with at least ' + str(HITS) + ' hits.')
	
##Apply max divergence if requested
	if MAX is not None:
		print('Will only output hits with divergences less than ' + str(MAX) + '.')
		OUT_ARRAY = OUT_ARRAY[OUT_ARRAY['diverge'] <= MAX]

##Apply min divergence if requested
	if MINDIV is not None:
		print('Will only output hits with divergences greater than ' + str(MINDIV) + '.')
		OUT_ARRAY = OUT_ARRAY[OUT_ARRAY['diverge'] >= MINDIV]
			
##Write the dataframe to a file
	if PREFIX is None and AGE is None:
		print('Writing main output to ' + CUTPREFIX + '_rm.bed.')
		OUT_ARRAY.to_csv(CUTPREFIX + '_rm.bed', sep='\t', header=False, index=False)
	##Split into files if asked. Also check to see if there is a minumum hit number and act accordingly.
	##If no age is provided	
		if SPLIT is not None: 
			print('Split files by ' + SPLIT)
			if SPLIT in ['name', 'family', 'class']:
				CLUSTERED = OUT_ARRAY.sort_values([SPLIT])
				for SPLITVALUE in CLUSTERED[SPLIT].unique():
					if HITS is not None:
						CLUSTEREDW = CLUSTERED[CLUSTERED[SPLIT]==SPLITVALUE]
						COUNT_ROW = CLUSTEREDW.shape[0]
						if COUNT_ROW >= HITS:
							CLUSTEREDW.to_csv(CUTPREFIX + '_' + SPLITVALUE + '_rm.bed', sep='\t', header=False, index=False)
					else:	
						CLUSTEREDW = CLUSTERED[CLUSTERED[SPLIT]==SPLITVALUE]
						CLUSTEREDW.to_csv(CUTPREFIX + '_' + SPLITVALUE + '_rm.bed', sep='\t', header=False, index=False)
			else:
				print('Splitting options are by name, family, and class.')			
	##if an age is provided....
	else: 
		print('Writing main output to ' + PREFIXAGE + '_rm.bed.')
		OUT_ARRAY.to_csv(PREFIXAGE + '_rm.bed', sep='\t', header=False, index=False)
	##Split into files if asked. Also check to see if there is a minumum hit number and act accordingly.
		if SPLIT is not None:
			print('Split files by ' + SPLIT)
			if SPLIT in ['name', 'family', 'class']:
				CLUSTERED = OUT_ARRAY.sort_values([SPLIT])
				for SPLITVALUE in CLUSTERED[SPLIT].unique():
					if HITS is not None:
						CLUSTEREDW = CLUSTERED[CLUSTERED[SPLIT]==SPLITVALUE]
						COUNT_ROW = CLUSTEREDW.shape[0]
						if COUNT_ROW >= HITS:
							CLUSTEREDW.to_csv(PREFIXAGE + '_' + SPLITVALUE + '_rm.bed', sep='\t', header=False, index=False)
					else:	
						CLUSTEREDW = CLUSTERED[CLUSTERED[SPLIT]==SPLITVALUE]
						CLUSTEREDW.to_csv(PREFIXAGE + '_' + SPLITVALUE + '_rm.bed', sep='\t', header=False, index=False)
			else:
				print('Splitting options are by name, family, and class.')			


##Get arguments function
def get_args():
	parser = argparse.ArgumentParser(description="Will process a .align or .align.gz file from RepeatMasker to generate multiple sorted output files in .bed format", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', type=str, help='Name of the .align or .align.gz file to be parsed', required=True)
	parser.add_argument('-m', '--minsize', type=int, help='Minimum size of hit to include in sorted file', default = 100)
	parser.add_argument('-s', '--sortcriterion', type=str, help='Sort criterion, i.e. size, name, family, class, size, or divergence (diverge), etc.')
	parser.add_argument('-p', '--prefix', type=str, help='Prefix to use for output file - default is first field of input filename')
	parser.add_argument('-a', '--age', type=str, help='Age designation, for example young, old, etc. depending on divergence values allowed. Optional.')
	parser.add_argument('-sp', '--split', type=str, help='Split into files based on name, family, class? This is optional.')
	parser.add_argument('-n', '--minhitnum', type=int, help='Minimum number of hits in file before being created. Only implemented if --split option is invoked. Optional.')
	parser.add_argument('-d', '--maxdiverge', type=float, help='Maximum divergence allowed in output file.')
	parser.add_argument('-dmin', '--mindiverge', type=float, help='Minimum divergence allowed in output file.')
	parser.add_argument('-r', '--rerun', type=str, help='Is this a rerun? Does a div.out file already exist? y or n.', default = 'n')

	args = parser.parse_args()
	ALIGN = args.input
	MINSIZE = args.minsize
	PREFIX = args.prefix
	AGE = args.age
	CRITERION = args.sortcriterion
	SPLIT = args.split
	HITS = args.minhitnum
	MAX = args.maxdiverge
	MINDIV = args.mindiverge
	RERUN = args.rerun

	return ALIGN, MINSIZE, PREFIX, AGE, CRITERION, SPLIT, HITS, MAX, MINDIV, RERUN

if __name__ =="__main__":main()
		
#gzip $ABBREV".align" &

#gzip $ABBREV"_div.out" &

#gzip $ABBREV"_wCpG_div.out" &
