#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N mitobim
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q raycc,ray512cc
#$ -pe fill 1
#$ -P communitycluster

#This file will process raw Illumina data using Trimmomatic.  This will be followed by mapping to a reference mitochondrial genome using MIRA4 to create a new mitochondrial genome assembly.

BASEDIR=/lustre/scratch/daray/mitobim
WORKDIR=$BASEDIR/output

mkdir $WORKDIR
cd $WORKDIR

REFGENOME=mBra_mito_full.fa  	# your reference genome for the assembly
REF_HOME=$BASEDIR/reference	#the location of your reference genome
REF=$(basename $REFGENOME _full.fa)

RAW_READS_HOME=/lustre/scratch/daray/uce.partition1/fastq   #the location of your raw data
mkdir $BASEDIR/data_raw
UNZIPPED_RAW_HOME=$BASEDIR/data_raw
mkdir $BASEDIR/support_files
SUPPORT_FILES=$BASEDIR/support_files	#where the support files like the adapter sequences will be located.

######
#set up alias' for major programs
######
BWA_HOME=/lustre/work/apps/bwa-0.6.2
SAMTOOLS_HOME=/lustre/work/apps/samtools-1.2
SAMTOOLS1_8_HOME=/lustre/work/apps/samtools-0.1.18
PICARD_HOME=/lustre/work/apps/picard-tools-1.91
BCFTOOLS_HOME=/lustre/work/apps/samtools-0.1.18/bcftools
RAY_SOFTWARE=/lustre/work/daray/software
TRIM_HOME=/lustre/work/apps/Trimmomatic-0.27
FASTX_HOME=/lustre/work/apps/fastx_toolkit-0.0.14/bin
VCFTOOLS_HOME=/lustre/work/daray/software/vcftools_0.1.12b/bin
MIRA_HOME=/lustre/work/apps/mira


for RAW_READ_FILE in $RAW_READS_HOME/*R1_001.fastq.gz
do
	SAMPLE_ID=$(basename $RAW_READ_FILE _L001_R1_001.fastq.gz)
#Unzip the raw reads into the processed_reads directory
	gunzip -c $RAW_READS_HOME/$SAMPLE_ID"_L001_R1_001.fastq.gz" >$UNZIPPED_RAW_HOME/$SAMPLE_ID"_R1.fastq"
	gunzip -c $RAW_READS_HOME/$SAMPLE_ID"_L001_R2_001.fastq.gz" >$UNZIPPED_RAW_HOME/$SAMPLE_ID"_R2.fastq"

mkdir $WORKDIR/$SAMPLE_ID
cd $WORKDIR/$SAMPLE_ID


#======================
#MIRA4 assembly 
#Create manifest.config for MIRA
echo -e "\n#manifest file for basic mapping assembly with illumina data using MIRA 4\n\nproject = initial-mapping-of-"$SAMPLE_ID"-to-mBra-mt\n\njob=genome,mapping,accurate\n\nparameters = -NW:mrnl=0 -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no\n\nreadgroup\nis_reference\ndata = $REF_HOME/$REFGENOME\nstrain = mBra-mt-genome\n\nreadgroup = reads\ndata = "$UNZIPPED_RAW_HOME/$SAMPLE_ID"_R1.fastq" $UNZIPPED_RAW_HOME/$SAMPLE_ID"_R2.fastq\ntechnology = solexa\nstrain = "$SAMPLE_ID"\n" > $SUPPORT_FILES/$SAMPLE_ID"_manifest.conf"

#Run MIRA
$MIRA_HOME/bin/mira $SUPPORT_FILES/$SAMPLE_ID"_manifest.conf"

#======================
#MITObim assembly
#Bait and iteratively map to the reference genome using MITObim
perl $RAY_SOFTWARE/MITObim_1.8.pl \
	-start 1 \
	-end 10 \
	-sample $SAMPLE_ID \
	-ref mBra-mt-genome \
	-readpool $UNZIPPED_RAW_HOME/$SAMPLE_ID"_R1.fastq" $UNZIPPED_RAW_HOME/$SAMPLE_ID"_R2.fastq" \
	-maf $WORKDIR/$SAMPLE_ID/"initial-mapping-of-"$SAMPLE_ID"-to-mBra-mt_assembly"/"initial-mapping-of-"$SAMPLE_ID"-to-mBra-mt_d_results"/"initial-mapping-of-"$SAMPLE_ID"-to-mBra-mt_out.maf" \
	&> $SAMPLE_ID".log"

cd ..
done

#-bash-4.1$ /PATH/TO/MITObim.pl -start 1 -end 10 -sample testpool -ref Salpinus_mt_genome -readpool reads.fastq -maf initial-mapping-testpool-to-Salpinus-mt_assembly/initial-mapping-testpool-to-Salpinus-mt_d_results/initial-mapping-testpool-to-Salpinus-mt_out.maf &> log








