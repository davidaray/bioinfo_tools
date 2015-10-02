#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N myoYum_STAR
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q R2D2,Chewie
#$ -pe fill 9 
#$ -P communitycluster

BASEDIR=/lustre/scratch/daray/Ray_low_cov_work/star-1
WORKDIR=$BASEDIR/output

cd $WORKDIR

THREADS=8  

ALN_HOME=$BASEDIR/alignments	#Where your initial tophat alignments will end up
ASSEMBLIES_HOME=$BASEDIR/assemblies	#where your cufflinks assemblies will go
MERGED=$BASEDIR/merged
RAW_READS_HOME=$BASEDIR/reads   #the location of your raw data
PROCESSED_READS_HOME=$BASEDIR/processed_reads
DRAFTS_HOME=$BASEDIR/drafts #the location of the genome drafts
GENOME_HOME=$BASEDIR/genome	#the location of your genome indexes
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
STAR_HOME=/lustre/work/apps/STAR-2.4/bin


################################################################################
# 1 Generate index for genome with STAR
#~~~~~~~~~~~
GENOME=TK186200
DRAFT=myoYum
RNASEQ_TAXON=myoYum

$STAR_HOME/STAR \
	--runThreadN $THREADS	\
	--runMode genomeGenerate	\
	--genomeDir $GENOME_HOME	\
	--genomeFastaFiles $DRAFTS_HOME/$GENOME"_mem.fa"

echo $DRAFT"_index_finished" |  mailx -s $DRAFT"_index_finished" david.4.ray@gmail.com	

################################################################################
# 2 Map RNA-Seq reads to genome with STAR
#~~~~~~~~~~~

SAMPLE1=TK186200-L
SAMPLE2=TK186200-T
SAMPLE3=TK186220-L
SAMPLE4=TK186220-T

$STAR_HOME/STAR \
        --runThreadN $THREADS   \
        --genomeDir $GENOME_HOME        \
        --readFilesIn \
		$PROCESSED_READS_HOME/$SAMPLE1"_R1_paired.fastq.gz",$PROCESSED_READS_HOME/$SAMPLE2_R1_paired.fastq.gz",$PROCESSED_READS_HOME/$SAMPLE3_R1_paired.fastq",$PROCESSED_READS_HOME/$SAMPLE4_R1_paired.fastq" $PROCESSED_READS_HOME/$SAMPLE1"_R2_paired.fastq.gz",$PROCESSED_READS_HOME/$SAMPLE2_R2_paired.fastq.gz",$PROCESSED_READS_HOME/$SAMPLE2_R3_paired.fastq",$PROCESSED_READS_HOME/$SAMPLE4_R2_paired.fastq"	\
        --readFilesCommand zcat \
        --quantMode TranscriptomeSAM    \
        --outFileNamePrefix $DRAFT"_v_"$RNASEQ_TAXON	\
		--outSAMtype BAM SortedByCoordinate	\
		--outBAMsortingThreadN $THREADS	\

echo $DRAFT"_v_"$RNASEQ_TAXON"_mapping_finished" |  mailx -s $DRAFT"_v_"$RNASEQ_TAXON"_mapping_finished" david.4.ray@gmail.com	

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

