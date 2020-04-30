#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N psmc
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe sm 36
#$ -P quanah


####USAGE and PREPWORK
# qsub psmc_pipeline.sh <R1 reads> <R2 reads> <Assemblies directory> <Assembly for mapping> <Prefix> <Working directory> 
#
#NOTE: do not use relative paths. Haven't figured out how to handle those yet.
##
# Input $1 <R1 reads> = fastq.gz containing first set of paired illumina reads for the target species
# Input $2 <R2 reads> = fastq.gz containing second set of paired illumina reads for the target species
# Input $3 <Assemblies directory> for mapping> = Path to the genome assemblies required. Ex. 'cPer' for the bat, Carollia perspicillata 
# Input $4 <Assembly for mapping> = Filename of the assembly to which the reads will be mapped. Should be saved in <Assemblies directory>. Ex. 'cPer.fa' or something similar
# Input $5 <Prefix> = The prefix you want to use for your output. Ex. mAus for Myotis austroriparius
# Input $6 <Working directory> = The main working directory for your analyses. Subdirectories 
# will be created within.
# 
# Example -
# cd <working directory> 
# qsub psmc_pipeline.sh R1_file.fq R2_file.fq assemblies cPer.fa cPer_psmc .
# translation = "qsub this script using R1_file.fq, R2_file.fq, with the 
#                cPer.fa genome assembly in the 'assemblies' 
#                directory, and using the prefix 'cPer_psmc' for all output files."
#

R1=$1 #Read1 files
R2=$2 #Read2 files
ASSEMBLIES=$3
GENOME=$4 #Genome fasta to which you are mapping with simple name.
MAP2=$(echo $GENOME | awk -F'[.]' '{print $1}') #The genome to which you're mapping, truncated name.
PREFIX=$5 #The species you're analyzing
WORKDIR=$6/$PREFIX/$MAP2
#Source genome, MyoMyo.fa is renamed from 
#/lustre/scratch/daray/bat1k_TE_analyses/maskerruns/assemblies/mMyoMyo_m19_AffsNnoesSC.p1.fa

module load intel bowtie2 bcftools/1.9 samtools 
module swap intel gnu/5.4.0 
module load gnuplot
PSMC=/lustre/work/daray/software/psmc
PSMCUTILS=/lustre/work/daray/software/psmc/utils
PICARD=/lustre/work/daray/software/picardtools/picard-tools-2.5.0
mkdir -p $WORKDIR
cd $WORKDIR

# Generate a consensus sequence from the mapped reads with SAMtools’s mpileup or 
# bcftools command (samtools mpileup is deprecated.)
# In the PSMC tutorial they use samtools but is deprecated,
# it is important to notice that at the tutorial the argument for the reference is -uf, 
# however in bcftools I did not find this argument, only -f.
# -d sets and minimum read depth and -D sets the maximum.
# It is recommended to set -d to a third of the average depth and -D to twice.

#If index does not exist, make it usig bowtie2.
[ ! -f $ASSEMBLIES/$MAP2".1.bt2" ] && bowtie2-build --threads 36 $ASSEMBLIES/$MAP2.fa $ASSEMBLIES/$MAP2

#Map the reads to the selected assembly and generate a bam file.
[ ! -f $PREFIX".bam" ] && bowtie2 --threads 36 -x $ASSEMBLIES/$MAP2 -p 36 -1 $R1 -2 $R2 \
	| samtools view -Sb >$PREFIX".bam" /dev/fd/0
	
#Sort the bam file
[ ! -f $PREFIX".sort.bam" ] && java -jar $PICARD/picard.jar SortSam I=$PREFIX.bam O=$PREFIX".sort.bam" SORT_ORDER=coordinate

#create a pileup in and map to the genome
[ ! -f $PREFIX"2"$MAP2".fq.gz" ] && bcftools mpileup --threads 36 -C 50 -f $ASSEMBLIES/$MAP2".fa" $PREFIX".sort.bam" | bcftools call -c - | vcfutils.pl vcf2fq -d 7 -D 50 | gzip > $PREFIX"2"$MAP2".fq.gz"

# Convert this fastq file to the input format for PSMC
[ ! -f $PREFIX"2"$MAP2".psmcfa" ] && $PSMCUTILS/fq2psmcfa $PREFIX"2"$MAP2".fq.gz" > $PREFIX"2"$MAP2".psmcfa"

# Split sequences to perform bootstrapping
[ ! -f $PREFIX"2"$MAP2".split.psmcfa" ] &&  $PSMCUTILS/splitfa $PREFIX"2"$MAP2".psmcfa" > $PREFIX"2"$MAP2".split.psmcfa"

# Run PSMC, using the default options
[ ! -f $PREFIX"2"$MAP2".psmcfa" ] && $PSMC/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $PREFIX"2"$MAP2".psmc" $PREFIX"2"$MAP2".psmcfa" 

#seq 100 | xargs -i echo $PSMC/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.psmc $PREFIX"2"$MAP2".split.psmcfa" | sh

[ ! -f $PREFIX"2"$MAP2".combined.psmc" ] && seq 100 | parallel -I% --max-args 1 $PSMC/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-%.psmc $PREFIX"2"$MAP2".split.psmcfa" | sh

[ ! -f $PREFIX"2"$MAP2".combined.psmc" ] && cat $PREFIX"2"$MAP2".psmc" round-*.psmc > $PREFIX"2"$MAP2".combined.psmc"

rm round-*.psmc

# Generate PSMC plot, using the per-generation mutation rate -u and the generation time in years -g.
$PSMCUTILS/psmc_plot.pl -u 2.2e-9 -g 4 -p  $PREFIX"2"$MAP2".combined.g4.psmc" $PREFIX"2"$MAP2".combined.psmc"

$PSMCUTILS/psmc_plot.pl -u 2.2e-9 -g 1.2 -p $PREFIX"2"$MAP2".combined.g1.2.psmc" $PREFIX"2"$MAP2".combined.psmc" 

