#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N TH-1
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q Yoda
#$ -pe fill 10 
#$ -P communitycluster

BASEDIR=/lustre/scratch/daray/Ray_low_cov_work/transcript1
WORKDIR=$BASEDIR/bin

cd $WORKDIR

THREADS=9  

ALN_HOME=$BASEDIR/alignments	#Where your initial tophat alignments will end up
ASSEMBLIES_HOME=$BASEDIR/assemblies	#where your cufflinks assemblies will go
MERGED=$BASEDIR/merged
RAW_READS_HOME=$BASEDIR/reads   #the location of your raw data
PROCESSED_READS_HOME=$BASEDIR/processed_reads
INDEX=$BASEDIR/index_files
GENOMES_HOME=$BASEDIR/genomes	#the location of your processed data
SUPPORT_FILES=$BASEDIR/support_files	#where the support files like the adapter sequences will be located.

######
#set up alias' for major programs
######
BWA_HOME=/lustre/work/apps/bwa-0.7.12
SAMTOOLS_HOME=/lustre/work/apps/samtools-1.2
SAMTOOLS1_8_HOME=/lustre/work/apps/samtools-0.1.18
PICARD_HOME=/lustre/work/apps/picard-tools-1.91
BCFTOOLS_HOME=/lustre/work/apps/samtools-0.1.18/bcftools
RAY_SOFTWARE=/lustre/work/daray/software
TRIM_HOME=/lustre/work/apps/Trimmomatic-0.27
FASTX_HOME=/lustre/work/apps/fastx_toolkit-0.0.14/bin
VCFTOOLS_HOME=/lustre/work/daray/software/vcftools_0.1.12b/bin
BEDTOOLS_HOME=/lustre/work/apps/bedtools-2.17.0/bin
TOPHAT_HOME=/lustre/work/apps/tophat-2.1.0.Linux_x86_64
CUFFLINKS_HOME=/lustre/work/apps/cufflinks-2.2.1/bin
BOWTIE2_HOME=/lustre/work/apps/bowtie2-2.0.5

################################################################################
# 1 Build index of genome for mapping reads using Bowtie2
#~~~~~~~~~~~
#string $BOWTIE2_HOME/bowtie2-build	\
#	-f	\
#string	$GENOMES_HOME/M8132_mem.fa	\
#string	$INDEX/M8132_mem

#echo "M8132_index_finished" |  mailx -s "M8132_index_finished" david.4.ray@gmail.com	

################################################################################
# 2 Map RNA-Seq reads to genome with TopHat
#~~~~~~~~~~~

ABBREV=TK186200-T

$TOPHAT_HOME/tophat2 \
	-p $THREADS	\
	-o $ALN_HOME	\
	$INDEX/M8132_mem	\
	$PROCESSED_READS_HOME/$ABBREV"_R1_paired.fastq",$PROCESSED_READS_HOME/$ABBREV"_RX_unpaired.fastq"	\
	$PROCESSED_READS_HOME/$ABBREV"_R1_paired.fastq"

echo "M8132-TK186200-T_mapping_finished" |  mailx -s "M8132-TK186200-T_mapping_finished" david.4.ray@gmail.com	

################################################################################
# 3 Assemble transcripts for each sample
#~~~~~~~~~~~
#string $CUFFLINKS_HOME/cufflinks	\
#	-p $THREADS	\
#	-o $ASSEMBLIES_HOME	\
#string	$ALN_HOME/

#echo "M8132-TK186200-T_assembly_finished" |  mailx -s "M8132-TK186200-T_assembly_finished" david.4.ray@gmail.com	

################################################################################
# 4 Create list of assemblies
#~~~~~~~~~~~
#cd $ASSEMBLIES_HOME
#ls -d -1 $PWD/** > $SUPPORT_FILES/M8132_mem_assemblies.txt
#cd $WORKDIR

################################################################################
# 4 Create a single merged transcriptome annotation using cuffmerge
#~~~~~~~~~~~
#string $CUFFLINKS_HOME/cufflinks	\
#	-s $GENOMES_HOME/M8132_mem.fa	\
#	-o $MERGED	\
#string	$SUPPORT_FILES/M8132_mem_assemblies.txt
	
#echo "M8132-TK186200-T_merge_finished" |  mailx -s "M8132-TK186200-T_merge_finished" david.4.ray@gmail.com	

#sleep 5

