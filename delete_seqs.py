
from Bio import SeqIO
import argparse
import sys

def main():
##Use the get_args function
        FAMILIES, LIST, OUT = get_args()

        TELIST = set(LINE.strip() for LINE in open(LIST, 'r'))


        with open(OUT, 'w') as OUTPUT:
                for SEQRECORD in SeqIO.parse(FAMILIES, 'fasta'):
                        try:
                                TELIST.remove(SEQRECORD.name)
                        except KeyError:
                                SeqIO.write(SEQRECORD, OUTPUT, 'fasta')
                                continue
                if len(TELIST) != 0:
                        print(len(TELIST),'of the headers from list were not identified in the input fasta file.', file=sys.stderr)

def get_args():
        parser = argparse.ArgumentParser(description="Will delete sequences from a modified RepeatModeler output families.fa file using a text file list of headers. Headers must match exactly.", formatter_$
        parser.add_argument('-l', '--list', type=str, help='Your list of TEs to delete. Typically, $SUBNAME_young_TEs.txt.', required=True)
        parser.add_argument('-f', '--families', type=str, help='RepeatModeler families file in .fasta format. Typically, $SUBNAME-families_100.fa.', required=True)
        parser.add_argument('-o', '--output', type=str, help='The name of your output file.', required=True)
        args = parser.parse_args()
        FAMILIES = args.families
        LIST = args.list
        OUT = args.output

        return FAMILIES, LIST, OUT

if __name__ =="__main__":main()
