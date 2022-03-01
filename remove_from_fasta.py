from Bio import SeqIO
import argparse

def main():
##Use the get_args function
	FASTA, LIST, OUT = get_args()

	with open(LIST, 'r') as FILE:
		HEADERLIST = FILE.read().splitlines()
	
	with open(OUT, 'w') as OUTPUT:
		for SEQRECORD in SeqIO.parse(FASTA, 'fasta'):
			ID = SEQRECORD.id.split('#')[0]
			if ID not in HEADERLIST:
				SeqIO.write(SEQRECORD, OUTPUT, 'fasta')

def get_args():
	parser = argparse.ArgumentParser(description="Will remove sequences from a fasta file using a provided list of headers.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-l', '--list', type=str, help='Your list of TEs to remove. Header names should trimmed of #Class/Familyy', required=True)
	parser.add_argument('-f', '--fasta', type=str, help='Target fasta file in .fasta format. ', required=True)
	parser.add_argument('-o', '--output', type=str, help='The name of your output file.', required=True)
	args = parser.parse_args()
	FASTA = args.fasta
	LIST = args.list
	OUT = args.output
	
	return FASTA, LIST, OUT

if __name__ =="__main__":main()