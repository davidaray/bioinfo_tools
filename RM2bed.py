#USAGE: RM2bed.py -i <input> -m <minsize> -s <sortcriterion> -p <prefix(optional)>
#Will process a .align.gz or .align file from RepeatMasker to generate multiple sorted output 
#files in .bed format

#To do
#1. convert to python
#2. add input options
#2a. detect align.gz vs .align
#2b. detect minimum size of hit
#2c. detect output type for sorting, i.e. by size, by name, by family, etc.

from Bio import SeqIO
import argparse
import os
impore re
import gzip
import subprocess
import fileinput
import pandas as pd

##Get arguments function
def get_args():
	parser = argparse.ArgumentParser(description="Will process a .align or .align.gz file from RepeatMasker to generate multiple sorted output files in .bed format", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', type=str, help='Name of the .align or .align.gz file to be parsed', required=True)
	parser.add_argument('-m', '--minsize', type=int, help='Minimum size of hit to include in sorted file', default = 100)
	parser.add_argument('-s', '--sortcriterion', type=str, help='Sort criterion, i.e. size, name, family, class, etc.', required=True, default=size)
	parser.add_argument('-p', '--prefix', type=str, help='Prefix to use for output file - default is first field of input filename')

	args = parser.parse_args()
	ALIGN = args.input
	MINSIZE = args.minsize
	PREFIX = args.prefix
	CRITERION = args.sortcriterion
	
	return ALIGN, MINSIZE, PREFIX, CRITERION

##use the get_args function
ALIGN, MINSIZE, PREFIX, CRITERION = get_args()

##If no prefix, get prefix from filename using the first part of the filename.
if PREFIX = ""
        PREFIX = re.split("[_.]", INPUT)[0]

##Determine file type using extension.
if ALIGN.lower().endswith(gz):
	TMP = PREFIX + '.tmp.file'
	IN = gzip.open(ALIGN, 'rb')
	with open(TMP, 'wb') as OUT:
		OUT.write(IN.read())
	IN.close()
	OUT.close()
elif ALIGN.lower().endswith(align):
	TMP = PREFIX + '.tmp.file'
	IN = open(ALIGN, 'r')
	with open(TMP, 'wb') as OUT:
		OUT.write(IN.read())
	IN.close()
	OUT.close()
else:
	print('Input file must be either .align.gz or .align.')
	exit()

##Call subprocess, calcDivergenceFromAlign.Ray.pl, using perl. Keep the new file for troubleshooting.
subprocess.check_call('perl /lustre/work/daray/software/RepeatMasker/util/calcDivergenceFromAlign.Ray.pl -noCpG {} >{}'.format(TMP, PREFIX + '_div.out'), shell=True)

##Call subprocess, RM.addsize.div.pl, to add size column	
subprocess.check_call('perl /lustre/work/daray/software/RM.addsize.div.pl {}'.format(PREFIX + '_div.out', PREFIX + '_div_size.out'), shell = TRUE)
 
#grep to remove lines with simple and satellite repeats 
BADWORDS = ['Simple','Satellite']
SIZE_OUT_LIST = []
with open(PREFIX + '_div_size.out') as SIZE_OUT_FILE:
        for LINE in SIZE_OUT_FILE:
			if not any(BADWORD in LINE for BADWORD in BADWORDS):
				SIZE_OUT_LIST.append(LINE)

#sed to remove extra '\', "#", and 'kimura=' and replace the first two with tabs.			
for LINE in SIZE_OUT_LIST:
	LINE.replace('\\', '\t')
	LINE.replace('kimura=', '')
	LINE.replace('#', '\t')

##Create a tmp file from the new list for use as an dataframe in pandas
NEWFILE = open('tmp.bed', 'w')
for LINE in SIZE_OUT_LIST:
        NEWFILE.write(LINE)
NEWFILE.close()

##Create a pandas dataframe from the tmp file, adding headers as you do so.
OUT_ARRAY = pd.read_table('tmp.bed', names=['A', 'B', 'C', 'D', 'chrom', 'start', 'stop', 'size', 'E', 'strand', 'name', 'class', 'family', 'F', 'G', 'H', 'I', 'div'])

##Delete the tmp file
os.remove('tmp.bed')

##Rearrange the headers as you want
RM_OUT_ARRAY = RM_OUT_ARRAY[['chrom', 'start', 'stop', 'name', 'size', 'strand' 'class', 'family', 'div']]
#print(RM_OUT_ARRAY)

##Replace 'C' with '-' in the 'strand' column. 
RM_OUT_ARRAY.strand = RM_OUT_ARRAY.strand.replace({'C':'-'})

##Write the dataframe to a file
RM_OUT_ARRAY.to_csv(PREFIX + '_rm.bed', sep='\t', header=False, index=False)

#gzip $ABBREV".align" &

#gzip $ABBREV"_div.out" &

#gzip $ABBREV"_wCpG_div.out" &


