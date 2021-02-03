#!/bin/bash 
#SBATCH --job-name=<NAME>_RM2bed
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --mem-per-cpu=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-user david.a.ray@ttu.edu

## Usage:
## You will need to replace <NAME> with the genome assembly identifier as appropriate.
## DIR may need modifying.
## sed "s|DIRPATH|path to your directory|g" template_RM_rm2bed.sh ><new file name>
## You will need to include the path to your genome assembly, GENOMEPATH.
## Genome assembly should be named <NAME>.fa
## If running lots of genomes, here's a quick way to do it using a 'list.txt' file that has all of your taxon names in it.
## cat list.txt | while read i; do sed "s/<NAME>/$i/g" <new file name> >${i}_RM2bed.sh; done
## Then, submit the scripts using a similar loop.

#This script uses a python program. You need to make sure python is loaded and ready to go by activating conda.
#Your conda installation must have biopython, pandas, and pyfaidx installed.
. ~/conda/etc/profile.d/conda.sh
conda activate

#makes variables used in this script and in the other
GENOME=<NAME>
RUNTYPE=${GENOME}_RM
DIR=DIRPATH/$RUNTYPE

#Creates a working directory and then goes into it to make all the files
#mkdir /lustre/scratch/daray/zoonomia_rm/$RUNTYPE
cd $DIR

#Runs a python script RM2bed to generate one complete .bed file and several subfiles subdivided by TE class. Merges overlapping hits based using 'lower_divergence' criterion. 
[ ! -f ${GENOME}_rm.bed ] python  ~/gitrepositories/bioinfo_tools/RM2bed.py -d . -sp class -p ${GENOME} -o lower_divergence ${GENOME}.fa.align.gz

#Calculate genome size by counting all characters in .masked.fa.gz file that are not on a line begining with >. Assign the value to GENOMESIZE. 
echo "calculating genome size"
GENOMESIZE=$(cat ${GENOME}.fa |  grep -v ">"  | wc | awk '{print $3-$1}')

echo "working on <NAME>"
#Examine selected .bed files of interest to us. Assign the total number of bases in each file to a variable. For example, DNABEDBP would represent the total number of bases in that assembly categorized as being derived from DNA transposons.
#UnspecifiedBEDBP=$(awk '{SUM += $5} END {print SUM}' ${GENOME}_Unspecified_rm.bed)
UnknownBEDBP=$(awk '{SUM += $5} END {print SUM}' ${GENOME}_Unknown_rm.bed)
SINEBEDBP=$(awk '{SUM += $5} END {print SUM}' ${GENOME}_SINE_rm.bed)
RetroposonBEDBP=$(awk '{SUM += $5} END {print SUM}' ${GENOME}_Retroposon_rm.bed)
RCBEDBP=$(awk '{SUM += $5} END {print SUM}' ${GENOME}_RC_rm.bed)
LTRBEDBP=$(awk '{SUM += $5} END {print SUM}' ${GENOME}_LTR_rm.bed)
LINEBEDBP=$(awk '{SUM += $5} END {print SUM}' ${GENOME}_LINE_rm.bed)
DNABEDBP=$(awk '{SUM += $5} END {print SUM}' ${GENOME}_DNA_rm.bed)
echo "all subgroups calculated"

#Sums all the groups
TEBEDBP=$(awk "BEGIN {print ($UnknownBEDBP + $SINEBEDBP + $RetroposonBEDBP + $RCBEDBP + $LTRBEDBP + $LINEBEDBP + $DNABEDBP)}")
echo "total calculated"

#Divide value of whateverBP from previous steps by value of GENOMESIZE. Assign to whateverPROPORTION.
TEPROPORTION=$(awk "BEGIN {print ($TEBEDBP/$GENOMESIZE)}")
LINEPROPORTION=$(awk "BEGIN {print ($LINEBEDBP/$GENOMESIZE)}")
SINEPROPORTION=$(awk "BEGIN {print ($SINEBEDBP/$GENOMESIZE)}")
UNKNOWNPROPORTION=$(awk "BEGIN {print ($UnknownBEDBP/$GENOMESIZE)}")
#UNSPECIFIEDPROPORTION=$(awk "BEGIN {print ($UnspecifiedBEDBP/$GENOMESIZE)}")
RCPROPORTION=$(awk "BEGIN {print ($RCBEDBP/$GENOMESIZE)}")
LTRPROPORTION=$(awk "BEGIN {print ($LTRBEDBP/$GENOMESIZE)}")
DNAPROPORTION=$(awk "BEGIN {print ($DNABEDBP/$GENOMESIZE)}")

#Create output file for all TEs
echo "Species TE_bp LINE_bp SINE_bp Genome_size TE_proportion LINE_proportion SINE_proportion LTR_proportion DNA_proportion RC_proportion Unknown_proportion" >te_table_total.txt 

#Prints in the total_te_table.txt the data previously found
echo <NAME> $TEBEDBP $LINEBEDBP $SINEBEDBP $GENOMESIZE $TEPROPORTION $LINEPROPORTION $SINEPROPORTION $LTRPROPORTION $DNAPROPORTION $RCPROPORTION $UNKNOWNPROPORTION >> te_table_total.txt 

#Does the same thing but with a divergence of all the TEs under .22
#YOUNGUnspecifiedBEDBP=$(awk '{ if ($9 <= 22) {print}}' <NAME>_Unspecified_rm.bed | awk '{SUM += $5} END {print SUM}' )
YOUNGUnknownBEDBP=$(awk '{ if ($9 <= 22) {print}}' ${GENOME}_Unknown_rm.bed | awk '{SUM += $5} END {print SUM}' )
YOUNGSINEBEDBP=$(awk '{ if ($9 <= 22) {print}}' ${GENOME}_SINE_rm.bed | awk '{SUM += $5} END {print SUM}' )
YOUNGRetroposonBEDBP=$(awk '{ if ($9 <= 22) {print}}' ${GENOME}_Retroposon_rm.bed | awk '{SUM += $5} END {print SUM}' )
YOUNGRCBEDBP=$(awk '{ if ($9 <= 22) {print}}' ${GENOME}_RC_rm.bed | awk '{SUM += $5} END {print SUM}' )
YOUNGLTRBEDBP=$(awk '{ if ($9 <= 22) {print}}' ${GENOME}_LTR_rm.bed | awk '{SUM += $5} END {print SUM}' )
YOUNGLINEBEDBP=$(awk '{ if ($9 <= 22) {print}}' ${GENOME}_LINE_rm.bed | awk '{SUM += $5} END {print SUM}' )
YOUNGDNABEDBP=$(awk '{ if ($9 <= 22) {print}}' ${GENOME}_DNA_rm.bed | awk '{SUM += $5} END {print SUM}' )
echo "all young subgroups calculated"

#Sums all the groups that was previously made
YOUNGTEBEDBP=$(awk "BEGIN {print ($YOUNGUnknownBEDBP + $YOUNGSINEBEDBP + $YOUNGRetroposonBEDBP + $YOUNGRCBEDBP + $YOUNGLTRBEDBP + $YOUNGLINEBEDBP + $YOUNGDNABEDBP)}")
echo "total young calculated"

YOUNGTEPROPORTION=$(awk "BEGIN {print ($YOUNGTEBEDBP/$GENOMESIZE)}")
YOUNGLINEPROPORTION=$(awk "BEGIN {print ($YOUNGLINEBEDBP/$GENOMESIZE)}")
YOUNGSINEPROPORTION=$(awk "BEGIN {print ($YOUNGSINEBEDBP/$GENOMESIZE)}")
YOUNGUNKNOWNPROPORTION=$(awk "BEGIN {print ($YOUNGUnknownBEDBP/$GENOMESIZE)}")
#YOUNGUNSPECIFIEDPROPORTION=$(awk "BEGIN {print ($YOUNGUnspecifiedBEDBP/$GENOMESIZE)}")
YOUNGRCPROPORTION=$(awk "BEGIN {print ($YOUNGRCBEDBP/$GENOMESIZE)}")
YOUNGLTRPROPORTION=$(awk "BEGIN {print ($YOUNGLTRBEDBP/$GENOMESIZE)}")
YOUNGDNAPROPORTION=$(awk "BEGIN {print ($YOUNGDNABEDBP/$GENOMESIZE)}")

#Create output file
echo "Species TE_bp Young_LINE_bp Young_SINE_bp Genome_size Young_TE_proportion Young_LINE_proportion Young_SINE_proportion Young_LTR_proportion Young_DNA_proportion Young_RC_proportion Young_Unknown_proportion" >te_table_young.txt 

#makes the young table by adding in the data
echo <NAME> $TEBEDBP $YOUNGLINEBEDBP $YOUNGSINEBEDBP $GENOMESIZE $YOUNGTEPROPORTION $YOUNGLINEPROPORTION $YOUNGSINEPROPORTION $YOUNGLTRPROPORTION $YOUNGDNAPROPORTION $YOUNGRCPROPORTION $YOUNGUNKNOWNPROPORTION >> te_table_young.txt



