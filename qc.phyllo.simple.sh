#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N QC
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q R2D2
#$ -pe fill 5 
#$ -P communitycluster



cd /lustre/scratch/daray/ray_phyllo/HiSeq/Sample_46394

zcat 46394_R1.fastq.gz | /lustre/work/apps/fastx_toolkit-0.0.14/bin/fastq_quality_filter -Q33 -q20 -p 50 -o 46394_R1_QC.fastq

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastx_quality_stats 					\
	-Q33		 						\
	-o 46394_R1_QC.stats 	\
	-i 46394_R1_QC.fastq

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastx_nucleotide_distribution_graph.sh 			\
	-i 46394_R1_QC.stats		\
	-o 46394_R1_QC_NUC.png	\
	-t 46394_R1_QC		

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastq_quality_boxplot_graph.sh 				\
	-i 46394_R1_QC.stats	\
	-o 46394_R1_QC_BOX.png	\
	-t 46394_R1_QC		

zcat 46394_R2.fastq.gz | /lustre/work/apps/fastx_toolkit-0.0.14/bin/fastq_quality_filter -Q33 -q20 -p 50 -o 46394_R2_QC.fastq

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastx_quality_stats 					\
	-Q33		 						\
	-o 46394_R2_QC.stats 	\
	-i 46394_R2_QC.fastq

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastx_nucleotide_distribution_graph.sh 			\
	-i 46394_R2_QC.stats		\
	-o 46394_R2_QC_NUC.png	\
	-t 46394_R2_QC		

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastq_quality_boxplot_graph.sh 				\
	-i 46394_R2_QC.stats	\
	-o 46394_R2_QC_BOX.png	\
	-t 46394_R2_QC

cd /lustre/scratch/daray/ray_phyllo/HiSeq/Sample_46395

zcat 46395_R1.fastq.gz | /lustre/work/apps/fastx_toolkit-0.0.14/bin/fastq_quality_filter -Q33 -q20 -p 50 -o 46395_R1_QC.fastq

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastx_quality_stats 					\
	-Q33		 						\
	-o 46395_R1_QC.stats 	\
	-i 46395_R1_QC.fastq

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastx_nucleotide_distribution_graph.sh 			\
	-i 46395_R1_QC.stats		\
	-o 46395_R1_QC_NUC.png	\
	-t 46395_R1_QC		

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastq_quality_boxplot_graph.sh 				\
	-i 46395_R1_QC.stats	\
	-o 46395_R1_QC_BOX.png	\
	-t 46395_R1_QC		
	
zcat 46395_R2.fastq.gz | /lustre/work/apps/fastx_toolkit-0.0.14/bin/fastq_quality_filter -Q33 -q20 -p 50 -o 46395_R2_QC.fastq

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastx_quality_stats 					\
	-Q33		 						\
	-o 46395_R2_QC.stats 	\
	-i 46395_R2_QC.fastq

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastx_nucleotide_distribution_graph.sh 			\
	-i 46395_R2_QC.stats		\
	-o 46395_R2_QC_NUC.png	\
	-t 46395_R2_QC		

/lustre/work/apps/fastx_toolkit-0.0.14/bin/fastq_quality_boxplot_graph.sh 				\
	-i 46395_R2_QC.stats	\
	-o 46395_R2_QC_BOX.png	\
	-t 46395_R2_QC	




