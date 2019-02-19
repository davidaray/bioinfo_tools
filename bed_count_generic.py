import argparse
import pandas as pd
import re

####MAIN function
def main():

##Use the get_args function
	BED = get_args()
	PREFIX = re.split("[.]", BED)[0]

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
		OUTFRAME.sort_values('COUNT', ascending=False)
		OUTFRAME.to_csv(PREFIX + '_' + CLASSTOCOUNT + '_counts.csv', sep='\t', index=False)
	
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
