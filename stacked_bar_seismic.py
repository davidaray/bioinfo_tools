import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
from pylab import savefig
import numpy as np
import re
import argparse
from matplotlib.colors import ListedColormap


####MAIN function
def main():

##Use the get_args function
	INPUT, TWEAK, LEGEND, DATATYPE, BOXWIDTH, SUBPLOTHORIZONAL, SUBPLOTVERTICAL, SUBPLOTWIDTH, SUBPLOTHEIGHT, COLUMNS = get_args()
	print('Input file = ' + INPUT)

	PREFIX = re.split("[.]", INPUT)[0]
	
	if TWEAK == 'y':
		print('Going with a tweaked plot. Gonna need some input.')
#Import the input dataframe.
		INPUTFRAME = pd.read_table(INPUT, sep = '\t', index_col=0)
#Output from 'catdata' will need to be transposed. My 'catdata' script is pretty personalized.
		INPUTFRAME = INPUTFRAME.transpose()
#Create a figure object
		plt.figure()
		colors = plt.cm.seismic(np.linspace(0, 1, 24)) #see https://matplotlib.org/tutorials/colors/colormaps.html
		if LEGEND == 'n':
#Create a stacked bar plot without a legend embedded
			FIG = INPUTFRAME.plot.bar(stacked=True, legend=False, color=colors)
#Save the figure to a png
			FIG.figure.savefig(PREFIX + '_nolegend.png')
		else:
#Create a subplot for the legend (subplot(nrows, ncols, index, **kwargs) https://matplotlib.org/api/_as_gen/matplotlib.pyplot.subplot.html)
			plt.subplot(111)
#Create a stacked bar plot without a legend embedded
			FIG = INPUTFRAME.plot.bar(stacked=True, legend=False, color=colors)
#Get position data for the bar plot
			BOX = FIG.get_position()
#Narrow the width of the plot to allow for the external legend
			FIG.set_position([BOX.x0, BOX.y0, BOX.width * BOXWIDTH, BOX.height])
#Create the legend and add to the plot. 
			plt.legend(loc=0, bbox_to_anchor=(SUBPLOTHORIZONAL, SUBPLOTVERTICAL, SUBPLOTWIDTH, 	SUBPLOTHEIGHT), ncol=COLUMNS, fontsize='small', borderaxespad=0.)
#Save the figure to a png
			FIG.figure.savefig(PREFIX + '_' + str(SUBPLOTHORIZONAL) + '_' +  str(SUBPLOTVERTICAL) + '_' +  str(SUBPLOTWIDTH) + '_' +  str(SUBPLOTHEIGHT) + '.png')
		
	else:
		INPUTFRAME = pd.read_table(INPUT, sep = '\t', index_col=0)
		INPUTFRAME = INPUTFRAME.transpose()
		plt.figure()
		colors = plt.cm.seismic(np.linspace(0, 1, 24))
		if LEGEND == 'n':
			FIG = INPUTFRAME.plot.bar(stacked=True, legend=False, color=colors)
			FIG.figure.savefig(PREFIX + '_nolegend.png')
		else:
			plt.subplot(111) 
			FIG = INPUTFRAME.plot.bar(stacked=True, legend=False, color=colors)
			lgd = plt.legend(loc=2, bbox_to_anchor=(1.01, 1), ncol=COLUMNS, borderaxespad=0.)
			FIG.figure.savefig(PREFIX + 'seismic' + '.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

##Get arguments function
def get_args():
	parser = argparse.ArgumentParser(description="Will process a saved pandas dataframe of data from RM2bed.py and processed through any of several catdata python scripts", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', type=str, help='Name of the input file to be parsed', required=True)
	parser.add_argument('-t', '--tweak', type=str, help='y = personalized, need extra arguments. Anything else = standard plot.', required=True)
	parser.add_argument('-l', '--legend', type=str, help='n = do not include legend. Default = y.', default='y', required=True)
	parser.add_argument('-d', '--datatype', type=str, help='CLASS, FAMILY, SINE, LINE, LTR, etc., depending on file')
	parser.add_argument('-b', '--boxwidth', type=float, help='Fraction of the standard size plot to allow for plotting the legend outside of that box. 0.8 (80%) will usually work')
	parser.add_argument('-ho', '--subplothorizontal', type=float, help='Number to indicate where to horizontally place the legend. Default is 1.01, just to the right of the main plot.', default = 1.01)
	parser.add_argument('-v', '--subplotvertical', type=float, help='Number to indicate where to vertically place the legend. Default is 1, just at the top of the main plot.', default = 1)
	parser.add_argument('-w', '--subplotwidth', type=float, help='Number to indicate how wide the legend should be. Default is .1, enough for one column.', default = .1)
	parser.add_argument('-sh', '--subplotheight', type=float, help='Number to indicate how tall the legend should be. Usually done automatically. Default = .1', default = 1)
	parser.add_argument('-c', '--numcol', type=int, help='Number of columns in the legend. Default = 1', default = 1)
#	parser.add_argument('-co', '--colorpalette', type=str, help='Choose from here, https://matplotlib.org/tutorials/colors/colormaps.html', default = 'RdYlGn' )

	args = parser.parse_args()
	INPUT = args.input
	TWEAK = args.tweak
	LEGEND = args.legend
	DATATYPE = args.datatype
	BOXWIDTH = args.boxwidth
	SUBPLOTHORIZONAL = args.subplothorizontal
	SUBPLOTVERTICAL = args.subplotvertical
	SUBPLOTWIDTH = args.subplotwidth
	SUBPLOTHEIGHT = args.subplotheight
	COLUMNS = args.numcol
#	COLORS = args.colorpalette

	return INPUT, TWEAK, LEGEND, DATATYPE, BOXWIDTH, SUBPLOTHORIZONAL, SUBPLOTVERTICAL, SUBPLOTWIDTH, SUBPLOTHEIGHT, COLUMNS

if __name__ =="__main__":main()	
	