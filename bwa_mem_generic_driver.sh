#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N massem1
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q Yoda
#$ -pe fill 35 
#$ -P communitycluster

#This file will process raw Illumina data using Trimmomatic.  This will be followed by mapping to a reference genome to create a new genome assembly.

#Depending on the size of the genome to be assembled, you may run into problems with memory. I had no trouble using Yoda but occassionally was unable to complete this using Chewie.  In those cases, it failed at the "sort BAM file" step.  I eventually began the process using lots of processors on Chewie and then, when it failed, would restart at the sorting step on Yoda.  

BASEDIR=/lustre/scratch/daray/Ray_low_cov_work/memassemblies1
WORKDIR=$BASEDIR/output

cd $WORKDIR

THREADS=34  #Line 9 sets this up to run 20 processors.  If you want fewer, make sure to change that line as well as this one.

#for RAW_READ_FILE in $BASEDIR/data_links/*.fastq.gz;
#do gunzip -c $BASEDIR/data_links/$RAW_READ_FILE >$RAW_READ_FILE.fastq &;	\
#done 
#wait

REFGENOME=myoLuc2.fa  	# your reference genome for the assembly
REF_HOME=$BASEDIR/Mluc.reference	#the location of your reference genome

RAW_READS_HOME=$BASEDIR/data_links   #the location of your raw data
PROCESSED_READS_HOME=$BASEDIR/data_processed	#the location of your processed data
QUALITY_INFO=$PROCESSED_READS_HOME/QC_files	#Where the quality stats will be saved
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


######
#Set up insert size.  This will be specific to the insert size for your particular taxon's library.
######
insSize=1000  

######
#make sure your genome file has no blank lines  - ALREADY DONE, NOT BEING USED HERE
######
#sed '/^$/d' $REF_HOME/$refgenome >tempGenome
#cp tempGenome $REF_HOME/$refgenome".clean"

#echo "spaces" |  mailx -s "spaces" david.4.ray@gmail.com

################################################################################
# Map reads to genome with BWA
#~~~~~~~~~~~

#[1a] Use bwa to index the genome.  May not be needed in some cases.
$BWA_HOME/bwa index \
	-a bwtsw \
	$REF_HOME/$REFGENOME


#####
# Preprocessing and QC steps - may not be necessary, depending on the job
#####	
for RAW_READ_FILE in $BASEDIR/data_links/*_R1.fastq
do ABBREV=$(basename $RAW_READ_FILE _R1.fastq) #This will be the name you use to process your files.  You will need to change the RP1, RP2, RU1, and RU2 slots just below here accordingly.

#These files will be generated by trimmomatic
RP1=$ABBREV"_R1_paired.fastq"
RP2=$ABBREV"_R2_paired.fastq"
RU1=$ABBREV"_R1_unpaired.fastq"
RU2=$ABBREV"_R2_unpaired.fastq"


###########
#Prepare the reads and get quality stats
###########
#[1] Use trimmomatic to generate paired and unpaired reads files
java -jar $TRIM_HOME/trimmomatic-0.27.jar \
	PE \
	-threads $THREADS \
	-phred33 \
	$RAW_READS_HOME/$ABBREV"_R1.fastq" \
	$RAW_READS_HOME/$ABBREV"_R2.fastq" \
	$PROCESSED_READS_HOME/$RP1 \
	$PROCESSED_READS_HOME/$RU1 \
	$PROCESSED_READS_HOME/$RP2 \
	$PROCESSED_READS_HOME/$RU2 \
	ILLUMINACLIP:$SUPPORT_FILES/TruSeq4-PE.fa:2:30:10 \
	LEADING:20 \
	TRAILING:20 \
	SLIDINGWINDOW:4:20 \
	MINLEN:33

	#[2] Generate quality stats 
	for FASTQ in $PROCESSED_READS_HOME/$ABBREV*.fastq
		do 	PROCESSED_FILE=$(basename $FASTQ .fastq)

		$FASTX_HOME/fastx_quality_stats 					\
			-Q33		 						\
			-o $QUALITY_INFO/$PROCESSED_FILE".stats" 	\
			-i $PROCESSED_READS_HOME/$PROCESSED_FILE".fastq"

		$FASTX_HOME/fastx_nucleotide_distribution_graph.sh 			\
			-i $QUALITY_INFO/$PROCESSED_FILE".stats"		\
			-o $QUALITY_INFO/$PROCESSED_FILE"_NUC.png"	\
			-t $QUALITY_INFO/$PROCESSED_FILE"_clipped"		

		$FASTX_HOME/fastq_quality_boxplot_graph.sh 				\
			-i $QUALITY_INFO/$PROCESSED_FILE".stats" 	\
			-o $QUALITY_INFO/$PROCESSED_FILE"_BOX.png"	\
			-t $QUALITY_INFO/$PROCESSED_FILE"_clipped"

	done
done

echo "_qc_finished" |  mailx -s "_qc_finished" david.4.ray@gmail.com

sleep 5

################################################################################
# Map reads to genome with BWA
#~~~~~~~~~~~

for RAW_READ_FILE in $BASEDIR/data_links/*_R1.fastq
	do ABBREV=$(basename $RAW_READ_FILE _R1.fastq)
echo "NOW WORKING "$ABBREV

#These files were generated by trimmomatic
RP1=$ABBREV"_R1_paired.fastq"
RP2=$ABBREV"_R2_paired.fastq"
RU1=$ABBREV"_R1_unpaired.fastq"
RU2=$ABBREV"_R2_unpaired.fastq"
		
#===================
# [1b] Map the reads to the genome
$BWA_HOME/bwa mem 			\
	-M	\
	-R "@RG\tID:$ABBREV\tSM:$ABBREV\tLB:lc_paired\tPL:illumina"	\
	$REF_HOME/$REFGENOME			\
    $PROCESSED_READS_HOME/$RP1	\
    $PROCESSED_READS_HOME/$RP2	\
	-t $THREADS	\
	> $ABBREV"_aln_pe.sam" 	
		
						
echo $ABBREV"_R1&2_map" |  mailx -s $ABBREV"_R1&2_map" david.4.ray@gmail.com	

	$BWA_HOME/bwa mem 			\
	-M	\
	-R "@RG\tID:$ABBREV\tSM:$ABBREV\tLB:lc_paired\tPL:illumina"     \
	$REF_HOME/$REFGENOME			\
       	$PROCESSED_READS_HOME/$ABBREV"_RX_cat.fastq"	\
	-t $THREADS	\
	> $ABBREV"_aln_se.sam" 	

echo $ABBREV"_RX_map" |  mailx -s $ABBREV"_RX_map" david.4.ray@gmail.com	
		
		
#===================
# [2] calculate the max insertion size
	maxInsSize=insSize*2  
	
#===================
# [3] use SAMtools to create bam files of the mapped reads

#convert paired sam file to bam
	$SAMTOOLS_HOME/samtools view 			\
		-Sb 					\
		-o $ABBREV"_aln_pe.bam" 	\
		$ABBREV"_aln_pe.sam"		

#convert unpaired sam file to bam
	$SAMTOOLS_HOME/samtools view 			\
		-Sb 					\
		-o $ABBREV"_aln_se.bam" 	\
		$ABBREV"_aln_se.sam" 		

#merge the paired and unpaired mapped reads to a single bam file
	$SAMTOOLS_HOME/samtools merge	\
		$ABBREV"_merge.bam" \
		$ABBREV"_aln_pe.bam"	\
		$ABBREV"_aln_se.bam"	

#Not sure what this does - sort bam file?
	$SAMTOOLS_HOME/samtools view 			\
		-F 4 					\
		-q 20 					\
		-b						\
		-o $ABBREV"_R3.bam" 	\
		$ABBREV"_merge.bam" 

#### for samtools v1.2	This is where queues other than Yoda sometimes run into trouble. Restart here if the run dies on the sorting step.	
	$SAMTOOLS_HOME/samtools sort 			\
		-O bam	\
		-o $ABBREV"_R3_sorted.bam"	\
		-T $ABBREV"_R3_sorted" \
		-@ $THREADS 					\
		$ABBREV"_R3.bam" 		

echo $ABBREV"_RX_sort" |  mailx -s $ABBREV"_RX_sort" david.4.ray@gmail.com	
		
#===================
# [4] remove sequencing duplicates from the sorted bam file w/ PICARD	
	java 						\
        -Xmx24g 				\
		-Djava.io.tmpdir=tmp 			\
		-jar $PICARD_HOME/MarkDuplicates.jar 	\
        	I=$ABBREV"_R3_sorted.bam" 	\
       		O=$ABBREV"_R3_noDup.bam" 	\
        	M=$ABBREV"_R3_dupMetric.out" 	\
        	REMOVE_DUPLICATES=true 			\
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 	\
		VALIDATION_STRINGENCY=SILENT 		\
		ASSUME_SORTED=TRUE 			\
		TMP_DIR=tmp

#==================
#[5] create pileup from noDup.bam

$SAMTOOLS1_8_HOME/samtools mpileup \
	-C50 \
	-f $REF_HOME/$REFGENOME \
	$ABBREV"_R3_noDup.bam" \
	>$ABBREV"_mPileUp_0_1_18.vcf"
	
--echo $ABBREV"_pileup_finished" |  mailx -s $ABBREV"_pileup_finished" david.4.ray@gmail.com

	
#=======================	
#[6] generate fasta consensus from pileup

perl $RAY_SOFTWARE/pileup2fasta_v1-4.pl \
	-i $ABBREV"_mPileUp_0_1_18.vcf" \
	-o $ABBREV".fa"	\
	-g $ABBREV".gff"	\
	-b 8	\
	-s		\
	-V		

	
echo $ABBREV"_fasta_finished" |  mailx -s $ABBREV"_fasta_finished" david.4.ray@gmail.com

	
$BEDTOOLS_HOME/bedtools genomecov \
	-ibam $ABBREV"_R3_noDup.bam" \
	-g $ABBREV".fa" \
	| grep genome >$ABBREV"_genome_cov.txt" 
	
echo $ABBREV"_assembly_finished" |  mailx -s $ABBREV"_assembly_finished" david.4.ray@gmail.com


done
sleep 5

