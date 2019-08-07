#!/bin/bash
#PBS -N ALD-MEI
#PBS -V
#PBS -o /mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/$PBS_JOBNAME.$PBS_JOBID.o
#PBS -e /mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/$PBS_JOBNAME.$PBS_JOBID.e
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
## Set the total memory for the job
#PBS -l mem=30G
## Set the max walltime for the job
#PBS -l walltime=500:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=8


# Check to load environment
source /opt/tools/hpcenv.sh

WORKING_DIR=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/
ANALYSIS_DIR=${WORKING_DIR}MELT/
MELT_DIR=/opt/tools/MELTv2.1.5/
MEI_LIST=${MELT_DIR}/me_refs/1KGP_Hg19/mei_list.txt
GENOME_FASTA=/mnt/causes-vnx1/GENOMES/hg19/hg19_bwa.fa
GENE_ANNO=/opt/tools/MELTv2.1.5/add_bed_files/1KGP_Hg19/hg19.genes.bed
BAM1=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD010_BWAmem_dupremoved_realigned.sorted.bam
BAM2=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD011_BWAmem_dupremoved_realigned.sorted.bam
BAM3=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD026_BWAmem_dupremoved_realigned.sorted.bam
BAM4=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD027_BWAmem_dupremoved_realigned.sorted.bam
BAM5=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD036_BWAmem_dupremoved_realigned.sorted.bam
BAM6=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD041_BWAmem_dupremoved_realigned.sorted.bam
BAM7=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD042_BWAmem_dupremoved_realigned.sorted.bam
BAM8=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD048_BWAmem_dupremoved_realigned.sorted.bam
BAM9=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD049_BWAmem_dupremoved_realigned.sorted.bam
BAM10=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD058_BWAmem_dupremoved_realigned.sorted.bam
BAM11=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD059_BWAmem_dupremoved_realigned.sorted.bam
BAM12=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/ALD065_BWAmem_dupremoved_realigned.sorted.bam

## Step 1 - Preprocess
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM1 \
#	-h $GENOME_FASTA
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM2 \
#	-h $GENOME_FASTA
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM3 \
#	-h $GENOME_FASTA
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM4 \
#        -h $GENOME_FASTA
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM5 \
#        -h $GENOME_FASTA
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM6 \
#        -h $GENOME_FASTA
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM7 \
#        -h $GENOME_FASTA
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM8 \
#        -h $GENOME_FASTA
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM9 \
#        -h $GENOME_FASTA
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM10 \
#        -h $GENOME_FASTA
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM11 \
#        -h $GENOME_FASTA
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Preprocess \
#        -bamfile $BAM12 \
#        -h $GENOME_FASTA
#
#
## Step 2 - Individual Analysis
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#	-w $ANALYSIS_DIR \
#	-t $MEI_LIST \
#	-c 30 \
#	-bamfile $BAM1 \
#	-h $GENOME_FASTA
#
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#        -w $ANALYSIS_DIR \
#        -t $MEI_LIST \
#        -c 30 \
#        -bamfile $BAM2 \
#        -h $GENOME_FASTA
#
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#        -w $ANALYSIS_DIR \
#        -t $MEI_LIST \
#        -c 30 \
#        -bamfile $BAM3 \
#        -h $GENOME_FASTA
#
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#        -w $ANALYSIS_DIR \
#        -t $MEI_LIST \
#        -c 30 \
#        -bamfile $BAM4 \
#        -h $GENOME_FASTA
#
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#        -w $ANALYSIS_DIR \
#        -t $MEI_LIST \
#        -c 30 \
#        -bamfile $BAM5 \
#        -h $GENOME_FASTA
#
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#        -w $ANALYSIS_DIR \
#        -t $MEI_LIST \
#        -c 30 \
#        -bamfile $BAM6 \
#        -h $GENOME_FASTA
#
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#        -w $ANALYSIS_DIR \
#        -t $MEI_LIST \
#        -c 30 \
#        -bamfile $BAM7 \
#        -h $GENOME_FASTA
#
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#        -w $ANALYSIS_DIR \
#        -t $MEI_LIST \
#        -c 30 \
#        -bamfile $BAM8 \
#        -h $GENOME_FASTA
#
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#        -w $ANALYSIS_DIR \
#        -t $MEI_LIST \
#        -c 30 \
#        -bamfile $BAM9 \
#        -h $GENOME_FASTA
#
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#        -w $ANALYSIS_DIR \
#        -t $MEI_LIST \
#        -c 30 \
#        -bamfile $BAM10 \
#        -h $GENOME_FASTA
#
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#        -w $ANALYSIS_DIR \
#        -t $MEI_LIST \
#        -c 30 \
#        -bamfile $BAM11 \
#        -h $GENOME_FASTA
#
#java -Xmx6G -jar ${MELT_DIR}MELT.jar IndivAnalysis \
#        -w $ANALYSIS_DIR \
#        -t $MEI_LIST \
#        -c 30 \
#        -bamfile $BAM12 \
#        -h $GENOME_FASTA
#
#
## Step 3 Group Discovery
#java -Xmx4G -jar ${MELT_DIR}MELT.jar GroupAnalysis \
#	-discoverydir $ANALYSIS_DIR \
#	-w $ANALYSIS_DIR \
#	-t $MEI_LIST \
#	-h $GENOME_FASTA \
#	-n $GENE_ANNO 
#	
#
## Step 4 Genotype
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#	-bamfile $BAM1 \
#	-t $MEI_LIST \
#	-h $GENOME_FASTA \
#	-w $ANALYSIS_DIR \
#	-p $ANALYSIS_DIR 
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#        -bamfile $BAM2 \
#        -t $MEI_LIST \
#        -h $GENOME_FASTA \
#        -w $ANALYSIS_DIR \
#        -p $ANALYSIS_DIR 
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#        -bamfile $BAM3 \
#        -t $MEI_LIST \
#        -h $GENOME_FASTA \
#        -w $ANALYSIS_DIR \
#        -p $ANALYSIS_DIR 
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#        -bamfile $BAM4 \
#        -t $MEI_LIST \
#        -h $GENOME_FASTA \
#        -w $ANALYSIS_DIR \
#        -p $ANALYSIS_DIR
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#        -bamfile $BAM5 \
#        -t $MEI_LIST \
#        -h $GENOME_FASTA \
#        -w $ANALYSIS_DIR \
#        -p $ANALYSIS_DIR
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#        -bamfile $BAM6 \
#        -t $MEI_LIST \
#        -h $GENOME_FASTA \
#        -w $ANALYSIS_DIR \
#        -p $ANALYSIS_DIR
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#        -bamfile $BAM7 \
#        -t $MEI_LIST \
#        -h $GENOME_FASTA \
#        -w $ANALYSIS_DIR \
#        -p $ANALYSIS_DIR
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#        -bamfile $BAM8 \
#        -t $MEI_LIST \
#        -h $GENOME_FASTA \
#        -w $ANALYSIS_DIR \
#        -p $ANALYSIS_DIR
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#        -bamfile $BAM9 \
#        -t $MEI_LIST \
#        -h $GENOME_FASTA \
#        -w $ANALYSIS_DIR \
#        -p $ANALYSIS_DIR
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#        -bamfile $BAM10 \
#        -t $MEI_LIST \
#        -h $GENOME_FASTA \
#        -w $ANALYSIS_DIR \
#        -p $ANALYSIS_DIR
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#        -bamfile $BAM11 \
#        -t $MEI_LIST \
#        -h $GENOME_FASTA \
#        -w $ANALYSIS_DIR \
#        -p $ANALYSIS_DIR
#
#java -Xmx2G -jar ${MELT_DIR}MELT.jar Genotype \
#        -bamfile $BAM12 \
#        -t $MEI_LIST \
#        -h $GENOME_FASTA \
#        -w $ANALYSIS_DIR \
#        -p $ANALYSIS_DIR

# Make VCF
java -Xmx2G -jar ${MELT_DIR}MELT.jar MakeVCF \
	-genotypingdir $ANALYSIS_DIR \
	-h $GENOME_FASTA \
	-t $MEI_LIST \
	-w $ANALYSIS_DIR \
	-p $ANALYSIS_DIR \
	-o $ANALYSIS_DIR


