import argparse
import pandas as pd
import re
import os
import subprocess

####MAIN function
def main():

##Use the get_args function
	BED = get_args()
	BASEPREFIX = os.path.basename(BED)
	PREFIX = re.split("[.]", BASEPREFIX)[0]
	ANIMAL = re.split("[_rm]", BASEPREFIX)[0]
	print(ANIMAL)

	IN_BED = pd.read_table(BED, sep='\t', names=['chrom', 'start', 'stop', 'name', 'size', 'strand', 'class', 'family', 'diverge'])

	CLASSES = ['Unknown', 'DNA', 'SINE', 'LINE', 'LTR', 'RC']
	
	def DO_COUNT(CLASSTOCOUNT):
		CLASS_IN_BED = IN_BED['class'] == CLASSTOCOUNT
		IN_BED_CLASS = IN_BED[CLASS_IN_BED]
		IN_UNIQUE = IN_BED_CLASS['name'].unique()
		
		COUNTS = []
		for UNIQUE in IN_UNIQUE:
			COUNT = len(IN_BED_CLASS[IN_BED_CLASS['name'] == UNIQUE])
			COUNTS.append(COUNT)
	
		OUTFRAME = pd.DataFrame({'NAME': IN_UNIQUE, 'COUNT': COUNTS})
		OUTFRAME.sort_values('COUNT', ascending=False, inplace=True)
		OUTFRAME.to_csv(PREFIX + '_' + CLASSTOCOUNT + '_counts.txt', sep='\t', index=False)
		
		DIR = r'/lustre/scratch/daray/200mammals/SINEandLINE/'
		if os.path.isdir(DIR + '/' + PREFIX + '/'):
#			print(PREFIX + ' exists.')
			OUTFRAME1000 = OUTFRAME[OUTFRAME['COUNT'] >= 1000]
			OUTFRAME1000.to_csv(DIR + '/' + PREFIX + '/' + ANIMAL + '_' + CLASSTOCOUNT + '_outframe1000.txt', sep='\t', index=False)
#			print(OUTFRAME1000)
			TELIST = OUTFRAME1000['NAME'].tolist()
#			print(TELIST)
			for TENAME in TELIST:
				subprocess.check_call('cp /lustre/scratch/daray/200mammals/analyses/{}/extract_align/muscle/{}*.muscle.fas /lustre/scratch/daray/200mammals/SINEandLINE/{}/' .format (ANIMAL, TENAME + '_', PREFIX), shell=True) 
		else:
#			print(PREFIX + ' does not exist.')
			os.mkdir(DIR + '/' + PREFIX + '/')
			OUTFRAME1000 = OUTFRAME[OUTFRAME['COUNT'] >= 1000]
			OUTFRAME1000.to_csv(DIR + '/' + PREFIX + '/' + ANIMAL + '_' + CLASSTOCOUNT + '_outframe1000.txt', sep='\t', index=False)
#			print(OUTFRAME1000)
			TELIST = OUTFRAME1000['NAME'].tolist()
#			print(TELIST)
			for TENAME in TELIST:
				subprocess.check_call('cp /lustre/scratch/daray/200mammals/analyses/{}/extract_align/muscle/{}*.muscle.fas /lustre/scratch/daray/200mammals/SINEandLINE/{}/' .format (ANIMAL, TENAME + '_', PREFIX), shell=True) 
	
	for CLASS in CLASSES:
		DO_COUNT(CLASS)

##Get arguments function
def get_args():
	parser = argparse.ArgumentParser(description="Will process a .bed file, identify putative families categorized as 'Unknown' and count the number of instances in the genome of that family.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-b', '--bed', type=str, help='Name of the .bed file to be parsed', required=True)
#	parser.add_argument('-m', '--minsize', type=int, help='Minimum size of hit to include in sorted file', default = 100)
#	parser.add_argument('-s', '--sortcriterion', type=str, help='Sort criterion, i.e. size, name, family, class, size, or divergence (diverge), etc.')
#	parser.add_argument('-p', '--prefix', type=str, help='Prefix to use for output file - default is first field of input filename')
#	parser.add_argument('-sp', '--split', type=str, help='Split into files based on name, family, class? This is optional.')
#	parser.add_argument('-n', '--minhitnum', type=int, help='Minimum number of hits in file before being created. Only implemented if --split option is invoked. Optional.')
#	parser.add_argument('-d', '--maxdiverge', type=float, help='Maximum divergence allowed in output file.')
#	parser.add_argument('-dmin', '--mindiverge', type=float, help='Minimum divergence allowed in output file.')

	args = parser.parse_args()
	BED = args.bed
#	MINSIZE = args.minsize
#	PREFIX = args.prefix
#	CRITERION = args.sortcriterion
#	SPLIT = args.split
#	HITS = args.minhitnum
#	MAX = args.maxdiverge
#	MINDIV = args.mindiverge

	return BED

if __name__ =="__main__":main()
		
#gzip $ABBREV".align" &

#gzip $ABBREV"_div.out" &

#gzip $ABBREV"_wCpG_div.out" &
