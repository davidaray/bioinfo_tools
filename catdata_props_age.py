import pandas as pd
import os
import itertools
import argparse

####MAIN function
def main():

##Use the get_args function
	SIZEFILE, PREFIX = get_args()
	print('Input file = ' + SIZEFILE)
	print('Prefix is = ' + PREFIX)

#Opens genomesizes.txt as a dictionary. Values taken from RepeatMasker output files 
	print('Reading in genome sizes file.')
	GENOMESIZES = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu'], squeeze=True, index_col=0)
	GENOMESIZESFRAME = pd.read_table(SIZEFILE, sep='\t', names=['taxon', 'genomesize', 'mu'], squeeze=True)

#Get taxon list from genome sizes file.
	TAXA = GENOMESIZESFRAME['taxon'].tolist()

	TAXALENGTH = len(TAXA)

#Initialize lists of names
	CLASSLIST = []
	FAMILYLIST = []
	TELIST = []

	for TAXON in TAXA:
		print('Reading in ' + str(TAXON) + '_' + str(PREFIX) + '_filtered.bed. Writing out as processed.bed file.')
		TAXONPROCESSEDBED = pd.read_table(TAXON + '_' + PREFIX + '_filtered.bed', sep='\t', index_col=False)
		TAXONPROCESSEDBED.columns = [TAXON + '_TE', TAXON + '_size', TAXON + '_class', TAXON + '_family', TAXON + '_div', TAXON + '_age']

		TAXONPROCESSEDBED.to_csv(TAXON + '_' + PREFIX +'_processed.bed', sep='\t', index=False)
	

#Generate lists of class names, family names, and TE names as well as files that correspond to the classes and families. Too many to save files for each individual TE.

	for TAXON in TAXA:
		print('Generating basic files for ' + str(TAXON) + '.')
		#Read in processed bed file that has all of the data from RepeatMasker. These files were	output by 
		#for TAXON in $TAXA; do cd /lustre/work/daray/heliconiine/all_TEs_analysis/hexplusbutterfly_repeatmasker/$TAXON; awk '{print $4 "\t" $5 "\t" $7 "\t" $8 "\t" $9}' $TAXON"_old_rm.bed" >../$TAXON"_old_processed.bed"; cd..; done
		#...old_rm.bed was generated using RM2bed.py, selecting for only old divergences
		TAXONFRAME = pd.read_table(TAXON + '_' + PREFIX + '_processed.bed', sep='\t', index_col=False)
		#calculate proportions and add them to the dataframe
		TAXONFRAME[TAXON + '_prop'] = TAXONFRAME[TAXON + '_size']/GENOMESIZES.at[TAXON, 'genomesize']	
		#rearrange the columns
		TAXONFRAME = TAXONFRAME[[TAXON + '_TE', TAXON + '_size', TAXON + '_prop', TAXON + '_class', TAXON + '_family', TAXON + '_div', TAXON + '_age']]
		#write the proportions dataframe to a file
		TAXONFRAME.to_csv(TAXON + '_' + PREFIX + '_processed_beds.txt', sep='\t', index=False)
		#Get a list of class names
		LIST = TAXONFRAME[TAXON + '_class'].tolist()
		CLASSLIST = CLASSLIST + LIST
		#Get a list of family names
		LIST = TAXONFRAME[TAXON + '_family'].tolist()
		FAMILYLIST = FAMILYLIST + LIST
		#Get a list of TE names
		LIST = TAXONFRAME[TAXON + '_TE'].tolist()
		TELIST = TELIST + LIST	
		#Sort the entire dataframe by class
		CLASSFRAME = TAXONFRAME.sort_values([TAXON + '_class'])

	CLASSLIST = set(CLASSLIST)
	FAMILYLIST = set(FAMILYLIST)
	TELIST = set(TELIST)	

#For every class name split the dataframe into a separate dataframe and save it
	if not os.path.exists('class_files'):
		os.makedirs('class_files')
	if not os.path.exists('family_files'):
		os.makedirs('family_files')
	for TAXON in TAXA:
		TAXONFRAME = pd.read_table(TAXON + '_' + PREFIX + '_processed_beds.txt', sep='\t', index_col=False)
		for SPLITVALUE in CLASSLIST:
			CLASSFRAMEOUT = TAXONFRAME[TAXONFRAME[TAXON + '_class']==SPLITVALUE]
			CLASSFRAMEOUT.to_csv('class_files/' + TAXON + '_' + SPLITVALUE + '_' + PREFIX + '_class_processed_beds.txt', sep='\t', index=False)
		#For every family name split the dataframe into a separate dataframe and save it
		for SPLITVALUE in FAMILYLIST:
			FAMILYFRAMEOUT = TAXONFRAME[TAXONFRAME[TAXON + '_family']==SPLITVALUE]
			FAMILYFRAMEOUT.to_csv('family_files/' + TAXON + '_' + SPLITVALUE + '_' + PREFIX + '_family_processed_beds.txt', sep='\t', index=False)

	print('Generating a combined dataframe for all data.')
	FIRSTTAXON = GENOMESIZESFRAME.at[0, 'taxon']
	SHORTTAXA=TAXA[1:]
	ALLCOMBINEDDATAFRAME = pd.read_table(FIRSTTAXON + '_' + PREFIX + '_processed_beds.txt', sep='\t', index_col=False, low_memory=False)
	for TAXON in SHORTTAXA:
		IMPORT = pd.read_table(TAXON + '_' + PREFIX + '_processed_beds.txt', sep='\t', index_col=False)
		ALLCOMBINEDDATAFRAME = pd.concat([ALLCOMBINEDDATAFRAME, IMPORT], axis=1, ignore_index=False)

	ALLCOMBINEDDATAFRAME.to_csv(PREFIX + '_all_taxa_all_hits_unmerged.txt', sep='\t', index=True)

	print('Generating a combined dataframe for classes.')
	#Sum proportions for each class and generate an output file
	COMBINEDCLASSFRAME = pd.DataFrame(0, index=range(len(CLASSLIST)), columns=range(TAXALENGTH))
	COMBINEDCLASSFRAME.index = CLASSLIST
	NEWNAMES = []
	for TAXON in TAXA:
		NEWNAMES.append(TAXON)   # deleted  + '_prop'
	COMBINEDCLASSFRAME.columns = NEWNAMES

	#Fill in the proportions for each class name
	for TAXON in TAXA:
		for CLASSNAME in CLASSLIST:
			#print('Working on ' + str(TAXON) + ' and ' + str(CLASSNAME))
			CLASSPROP = ALLCOMBINEDDATAFRAME.loc[ALLCOMBINEDDATAFRAME[TAXON + '_class'] == CLASSNAME, TAXON + '_prop'].sum()
			COMBINEDCLASSFRAME.loc[CLASSNAME, TAXON] = CLASSPROP  #deleted  + '_prop'

	COMBINEDCLASSFRAME.sort_index(inplace=True)		
	COMBINEDCLASSFRAME.to_csv(PREFIX + '_all_taxa_classes_merged_cats.txt', sep='\t', index=True)

	print('Generating a combined dataframe for all families.')
	#Sum proportions for each family and generate an output file
	COMBINEDFAMILYFRAME = pd.DataFrame(0, index=range(len(FAMILYLIST)), columns=range(TAXALENGTH))
	COMBINEDFAMILYFRAME.index = FAMILYLIST
	COMBINEDFAMILYFRAME.columns = NEWNAMES

	#Fill in the proportions for each family name
	for TAXON in TAXA:
		for FAMILYNAME in FAMILYLIST:
			#print('Working on ' + str(TAXON) + ' and ' + str(FAMILYNAME))
			FAMILYPROP = ALLCOMBINEDDATAFRAME.loc[ALLCOMBINEDDATAFRAME[TAXON + '_family'] == FAMILYNAME, TAXON + '_prop'].sum() 
			COMBINEDFAMILYFRAME.loc[FAMILYNAME, TAXON] = FAMILYPROP #deleted  + '_prop'

	COMBINEDFAMILYFRAME.sort_index(inplace=True)		
	COMBINEDFAMILYFRAME.to_csv(PREFIX + '_all_taxa_families_merged_cats.txt', sep='\t', index=True)

#	print('Generating a combined dataframe for all TEs.')
	#Sum proportions for each family and generate an output file
#	COMBINEDTEFRAME = pd.DataFrame(0, index=range(len(TELIST)), columns=range(TAXALENGTH))
#	COMBINEDTEFRAME.index = TELIST
#	COMBINEDTEFRAME.columns = NEWNAMES

	#Fill in the proportions for each family name
#	for TAXON in TAXA:
#		for TENAME in TELIST:
#			#print('Working on ' + str(TAXON) + ' and ' + str(TENAME))
#			TEPROP = ALLCOMBINEDDATAFRAME.loc[ALLCOMBINEDDATAFRAME[TAXON + '_TE'] == TENAME, TAXON + '_prop'].sum() 
#			COMBINEDTEFRAME.loc[TENAME, TAXON] = TEPROP

#	COMBINEDTEFRAME.sort_index(inplace=True)		
#	COMBINEDTEFRAME.to_csv(PREFIX + '_all_taxa_TEs_merged_cats.txt', sep='\t', index=True)

	#Get lists of each subset of families in each class 
	print('Generating lists of families within each class')
	for TAXON in TAXA:	
		DNANAME = pd.read_table('class_files/' + TAXON + '_DNA_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, names=[TAXON + '_TE', TAXON + '_size', TAXON + '_class', TAXON + '_family', TAXON + '_div', TAXON + '_age'])
		SINEFRAME = pd.read_table('class_files/' + TAXON + '_SINE_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, names=[TAXON + '_TE', TAXON + '_size', TAXON + '_class', TAXON + '_family', TAXON + '_div', TAXON + '_age'])
		LINEFRAME = pd.read_table('class_files/' + TAXON + '_LINE_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, names=[TAXON + '_TE', TAXON + '_size', TAXON + '_class', TAXON + '_family', TAXON + '_div', TAXON + '_age'])
		RCFRAME = pd.read_table('class_files/' + TAXON + '_RC_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, names=[TAXON + '_TE', TAXON + '_size', TAXON + '_class', TAXON + '_family', TAXON + '_div', TAXON + '_age'])
		LTRFRAME = pd.read_table('class_files/' + TAXON + '_LTR_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, names=[TAXON + '_TE', TAXON + '_size', TAXON + '_class', TAXON + '_family', TAXON + '_div', TAXON + '_age'])
	
		DNALIST = DNANAME[TAXON + '_family'].tolist()
		SINELIST = SINEFRAME[TAXON + '_family'].tolist()
		LINELIST = LINEFRAME[TAXON + '_family'].tolist()
		RCLIST = RCFRAME[TAXON + '_family'].tolist()
		LTRLIST = LTRFRAME[TAXON + '_family'].tolist()

	DNALIST = set(DNALIST)
	SINELIST = set(SINELIST)
	LINELIST = set(LINELIST)
	RCLIST = set(RCLIST)
	LTRLIST = set(LTRLIST)


	print('Generating combined class files.')
	#This section will create a dataframe with all lines from all classes. There will multiple lines for each TE with the details for that individual insertion.
	#For the next part, we need to load the first taxon in the list first and the rest in subsequent iterations. Hence the shorter taxon list.
	SHORTTAXA=TAXA[1:]

	print('SINEs')
	#Build SINE tables
	#Read in the first SINE bed file created above.
	SINEFRAME = pd.read_table('class_files/' + FIRSTTAXON + '_SINE_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, low_memory=False)

	#For every species, import the bed file and concatenate it with the first.
	for SINETAXON in SHORTTAXA:
		IMPORT = pd.read_table('class_files/' + SINETAXON + '_SINE_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, low_memory=False)
		SINEFRAME = pd.concat([SINEFRAME, IMPORT], axis=1, ignore_index=False)

	SINEFRAME.to_csv(PREFIX + '_all_taxa_SINE_processed_cats.txt', sep='\t', index=False)


	#This section will combine the data for individual TEs
	#Create a complete list of SINE names
	SINELIST = []
	for TAXON in TAXA:
		THISTAXONSINECOLUMN = TAXON + '_TE'
		THISSINELIST = SINEFRAME[THISTAXONSINECOLUMN].unique()
		for LINE in THISSINELIST:
			SINELIST.append(LINE)
	FINALSINELIST = set(SINELIST)

	#Create a blank dataframe to fill with dimensions needed and the names of the indices and columns
	COMBINEDSINEFRAME = pd.DataFrame(0, index=range(len(FINALSINELIST)), columns=range(TAXALENGTH))
	COMBINEDSINEFRAME.index = FINALSINELIST
	COMBINEDSINEFRAME.columns = NEWNAMES

	#Fill in the proportions for each TE name
	for TAXON in TAXA:
		for SINENAME in FINALSINELIST:
			#print('Working on ' + str(TAXON) + 'and ' + str(SINENAME))
			SINENAMEPROP = SINEFRAME.loc[SINEFRAME[TAXON + '_TE'] == SINENAME, TAXON + '_prop'].sum()
			COMBINEDSINEFRAME.loc[SINENAME, TAXON] = SINENAMEPROP

	COMBINEDSINEFRAME.sort_index(inplace=True)	
#	COMBINEDSINEFRAME.replace({'False': 0}, inplace=True)
	COMBINEDSINEFRAME.to_csv(PREFIX + '_all_taxa_SINE_merged_cats.txt', sep='\t', index=True)
		
	print('DNA transposons')
	#Build DNA tables
	#Repeat above for DNA transposons
	DNAFRAME = pd.read_table('class_files/' + FIRSTTAXON + '_DNA_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, low_memory=False)

	for DNATAXON in SHORTTAXA:
		IMPORT = pd.read_table('class_files/' + DNATAXON + '_DNA_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, low_memory=False)
		DNAFRAME = pd.concat([DNAFRAME, IMPORT], axis=1, ignore_index=False)

	DNAFRAME.to_csv(PREFIX + '_all_taxa_DNA_processed_cats.txt', sep='\t', index=False)

	DNALIST = []
	for TAXON in TAXA:
		THISTAXONDNACOLUMN = TAXON + '_TE'
		THISDNALIST = DNAFRAME[THISTAXONDNACOLUMN].unique()
		for LINE in THISDNALIST:
			DNALIST.append(LINE)
	FINALDNALIST = set(DNALIST)

	COMBINEDDNAFRAME = pd.DataFrame(0, index=range(len(FINALDNALIST)), columns=range(TAXALENGTH))
	COMBINEDDNAFRAME.index = FINALDNALIST
	COMBINEDDNAFRAME.columns = NEWNAMES

	for TAXON in TAXA:
		for DNANAME in FINALDNALIST:
			#print('Working on ' + str(TAXON) + ' and ' + str(DNANAME))
			DNANAMEPROP = DNAFRAME.loc[DNAFRAME[TAXON + '_TE'] == DNANAME, TAXON + '_prop'].sum()
			COMBINEDDNAFRAME.loc[DNANAME, TAXON] = DNANAMEPROP
		
	COMBINEDDNAFRAME.sort_index(inplace=True)
	COMBINEDDNAFRAME.to_csv(PREFIX + '_all_taxa_DNA_merged_cats.txt', sep='\t', index=True)

	print('LINEs')
	#Build LINE tables
	LINEFRAME = pd.read_table('class_files/' + FIRSTTAXON + '_LINE_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, low_memory=False)

	for LINETAXON in SHORTTAXA:
		IMPORT = pd.read_table('class_files/' + LINETAXON + '_LINE_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, low_memory=False)
		LINEFRAME = pd.concat([LINEFRAME, IMPORT], axis=1, ignore_index=False)

	LINEFRAME.to_csv(PREFIX + '_all_taxa_LINE_processed_cats.txt', sep='\t', index=False)

	LINELIST = []
	for TAXON in TAXA:
		THISTAXONLINECOLUMN = TAXON + '_TE'
		THISLINELIST = LINEFRAME[THISTAXONLINECOLUMN].unique()
		for LINE in THISLINELIST:
			LINELIST.append(LINE)
	FINALLINELIST = set(LINELIST)

	COMBINEDLINEFRAME = pd.DataFrame(0, index=range(len(FINALLINELIST)), columns=range(TAXALENGTH))
	COMBINEDLINEFRAME.index = FINALLINELIST
	COMBINEDLINEFRAME.columns = NEWNAMES

	for TAXON in TAXA:
		for LINENAME in FINALLINELIST:
			#print('Working on ' + str(TAXON) + ' and ' + str(LINENAME))
			LINENAMEPROP = LINEFRAME.loc[LINEFRAME[TAXON + '_TE'] == LINENAME, TAXON + '_prop'].sum()
			COMBINEDLINEFRAME.loc[LINENAME, TAXON] = LINENAMEPROP
		
	COMBINEDLINEFRAME.sort_index(inplace=True)
	COMBINEDLINEFRAME.to_csv(PREFIX + '_all_taxa_LINE_merged_cats.txt', sep='\t', index=True)

	print('LTR retrotransposons')
	#Build LTR tables
	LTRFRAME = pd.read_table('class_files/' + FIRSTTAXON + '_LTR_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, low_memory=False)

	for LTRTAXON in SHORTTAXA:
		IMPORT = pd.read_table('class_files/' + LTRTAXON + '_LTR_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, low_memory=False)
		LTRFRAME = pd.concat([LTRFRAME, IMPORT], axis=1, ignore_index=False)

	LTRFRAME.to_csv(PREFIX + '_all_taxa_LTR_processed_cats.txt', sep='\t', index=False)

	LTRLIST = []
	for TAXON in TAXA:
		THISTAXONLTRCOLUMN = TAXON + '_TE'
		THISLTRLIST = LTRFRAME[THISTAXONLTRCOLUMN].unique()
		for LTR in THISLTRLIST:
			LTRLIST.append(LTR)
	FINALLTRLIST = set(LTRLIST)

	COMBINEDLTRFRAME = pd.DataFrame(0, index=range(len(FINALLTRLIST)), columns=range(TAXALENGTH))
	COMBINEDLTRFRAME.index = FINALLTRLIST
	COMBINEDLTRFRAME.columns = NEWNAMES

	for TAXON in TAXA:
		for LTRNAME in FINALLTRLIST:
			#print('Working on ' + str(TAXON) + ' and ' + str(LTRNAME))
			LTRNAMEPROP = LTRFRAME.loc[LTRFRAME[TAXON + '_TE'] == LTRNAME, TAXON + '_prop'].sum()
			COMBINEDLTRFRAME.loc[LTRNAME, TAXON] = LTRNAMEPROP
		
	COMBINEDLTRFRAME.sort_index(inplace=True)
	COMBINEDLTRFRAME.to_csv(PREFIX + '_all_taxa_LTR_merged_cats.txt', sep='\t', index=True)

	print('Rolling circle transposons')
	#Build RC tables
	RCFRAME = pd.read_table('class_files/' + FIRSTTAXON + '_RC_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False, low_memory=False)

	for RCTAXON in SHORTTAXA:
		IMPORT = pd.read_table('class_files/' + RCTAXON + '_RC_' + PREFIX + '_class_processed_beds.txt', sep='\t', index_col=False)
		RCFRAME = pd.concat([RCFRAME, IMPORT], axis=1, ignore_index=False)

	RCFRAME.to_csv(PREFIX + '_all_taxa_RC_processed_cats.txt', sep='\t', index=False)

	RCLIST = []
	for TAXON in TAXA:
		THISTAXONRCCOLUMN = TAXON + '_TE'
		THISRCLIST = RCFRAME[THISTAXONRCCOLUMN].unique()
		for RC in THISRCLIST:
			RCLIST.append(RC)
	FINALRCLIST = set(RCLIST)

	COMBINEDRCFRAME = pd.DataFrame(0, index=range(len(FINALRCLIST)), columns=range(TAXALENGTH))
	COMBINEDRCFRAME.index = FINALRCLIST
	COMBINEDRCFRAME.columns = NEWNAMES

	for TAXON in TAXA:
		for RCNAME in FINALRCLIST:
			#print('Working on ' + str(TAXON) + ' and ' + str(RCNAME))
			RCNAMEPROP = RCFRAME.loc[RCFRAME[TAXON + '_TE'] == RCNAME, TAXON + '_prop'].sum()
			COMBINEDRCFRAME.loc[RCNAME, TAXON] = RCNAMEPROP
		
	COMBINEDRCFRAME.sort_index(inplace=True)
#	COMBINEDRCFRAME.replace({'False': 0}, inplace=True)
	COMBINEDRCFRAME.to_csv(PREFIX + '_all_taxa_RC_merged_cats.txt', sep='\t', index=True)


##Get arguments function
def get_args():
	parser = argparse.ArgumentParser(description="Will process a filtered.bed file output from filter_beds.py into various subfiles for downstream processing.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#	parser.add_argument('-i', '--input', type=str, help='Name of the rm.bed file to be parsed.', required=True)
	parser.add_argument('-g', '--sizefile', type=str, help='File containing two corresponding columns of taxon abbreviations and genome sizes in bp.', required=True)
#	parser.add_argument('-s', '--sortcriterion', type=str, help='Sort criterion, i.e. size, name, family, class, size, or divergence (diverge), etc.')
	parser.add_argument('-p', '--prefix', type=str, help='Prefix to put after taxon id. I use this when separating young and old elements.')
#	parser.add_argument('-sp', '--split', type=str, help='Split into files based on name, family, class? This is optional.')
#	parser.add_argument('-n', '--minhitnum', type=int, help='Minimum number of hits in file before being created. Only implemented if --split option is invoked. Optional.')
#	parser.add_argument('-d', '--maxdiverge', type=float, help='Maximum divergence allowed in output file.')
#	parser.add_argument('-dmin', '--mindiverge', type=float, help='Minimum divergence allowed in output file.')

	args = parser.parse_args()
#	INPUT = args.input
	SIZEFILE = args.sizefile
	PREFIX = args.prefix
#	CRITERION = args.sortcriterion
#	SPLIT = args.split
#	HITS = args.minhitnum
#	MAX = args.maxdiverge
#	MINDIV = args.mindiverge

	return SIZEFILE, PREFIX

if __name__ =="__main__":main()
	
