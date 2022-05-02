from Bio import SeqIO
import argparse

def main():
##Use the get_args function
	FAMILIES, LIST, OUT = get_args()

	with open(LIST, 'r') as FILE:
		TELIST = FILE.read().splitlines()
	
	with open(OUT, 'w') as OUTPUT:
		for SEQRECORD in SeqIO.parse(FAMILIES, 'fasta'):
			ID = SEQRECORD.id.split('#')[0]
			if ID not in TELIST:
				SeqIO.write(SEQRECORD, OUTPUT, 'fasta')

def get_args():
	parser = argparse.ArgumentParser(description="Will pull sequences from a modified RepeatModeler output families.fa file using a text file processed using our 200 mammals qsub.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-l', '--list', type=str, help='Your list of TEs to extract. Typically, $SUBNAME_young_TEs.txt.', required=True)
	parser.add_argument('-f', '--families', type=str, help='RepeatModeler families file in .fasta format. Typically, $SUBNAME-families_100.fa.', required=True)
	parser.add_argument('-o', '--output', type=str, help='The name of your output file.', required=True)
	args = parser.parse_args()
	FAMILIES = args.families
	LIST = args.list
	OUT = args.output
	
	return FAMILIES, LIST, OUT

if __name__ =="__main__":main()
