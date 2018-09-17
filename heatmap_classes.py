import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import seaborn as sns
from pylab import savefig
import argparse

####MAIN function
def main():

##Use the get_args function
	PREFIX, WIDTH, HEIGHT = get_args()
	print('File is = ' + PREFIX + '_all_taxa_classes_merged_cats.txt')
	print('Figure width is ' + str(WIDTH) + '.')
	print('Figure height is ' + str(HEIGHT) + '.')

	FIGUREFRAME = pd.read_table(PREFIX + '_all_taxa_classes_merged_cats.txt', sep='\t', index_col=0)
	plt.figure(figsize=(WIDTH, HEIGHT))
	FIGURE = sns.heatmap(FIGUREFRAME, annot=False, cmap='viridis', xticklabels=True, yticklabels=True)
	FIGURE.figure.savefig(PREFIX + '_heatmap_classes.png')


##Get arguments function
def get_args():
	parser = argparse.ArgumentParser(description="Will process an classs_merged_cats.txt file output from catdata_props script into a heatmap.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-a', '--age', type=str, help='Prefix to put after taxon id. I use this when separating young and old elements.')
	parser.add_argument('-w', '--width', type=float, help='Width of figure in inches.', required=True)
	parser.add_argument('-t', '--height', type=float, help='Height of figure in inches.', required=True)
#	parser.add_argument('-d', '--maxdiverge', type=float, help='Maximum divergence allowed in output file.')
#	parser.add_argument('-dmin', '--mindiverge', type=float, help='Minimum divergence allowed in output file.')

	args = parser.parse_args()
#	INPUT = args.input
	PREFIX = args.age
	WIDTH = args.width
	HEIGHT = args.height
#	SPLIT = args.split
#	HITS = args.minhitnum
#	MAX = args.maxdiverge
#	MINDIV = args.mindiverge

	return PREFIX, WIDTH, HEIGHT

if __name__ =="__main__":main()





	

