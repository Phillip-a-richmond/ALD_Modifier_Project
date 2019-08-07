#!/bin/bash
#PBS -N Population_CNV_Annotation
#PBS -V
#PBS -o /mnt/causes-vnx1/CAUSES/CNVS/Annotation.o
#PBS -e /mnt/causes-vnx1/CAUSES/CNVS/Annotation.e
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
## Set the total memory for the job
#PBS -l mem=24gb
## Set the max walltime for the job
#PBS -l walltime=200:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=6
NSLOTS=$PBS_NUM_PPN
umask 0002
source /opt/tools/hpcenv.sh

export ANNOTSV=/mnt/causes-vnx1/DATABASES/AnnotSV_2.1/
source /opt/tools/hpcenv.sh
source activate SVcalling


WORKDIR=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD_Smoove/
cd $WORKDIR
INPUT_VCF=${WORKDIR}/ALD-smoove.genotyped.vcf
OUTPUT_DIR=$WORKDIR
OUTPUT_FILE=ALD_SMOOVE_AnnotSV

$ANNOTSV/bin/AnnotSV \
	-SVinputFile $INPUT_VCF \
	-SVinputInfo 0 \
	-SVminSize 100 \
	-genomeBuild GRCh37 \
	-bedtools /opt/tools/bedtools/bin/bedtools \
	-outputDir $OUTPUT_DIR \
	-outputFile ${OUTPUT_FILE}_split \
	-typeOfAnnotation split \
	-reciprocal yes \
	-overlap 70

$ANNOTSV/bin/AnnotSV \
	-SVinputFile $INPUT_VCF \
	-SVinputInfo 0 \
	-SVminSize 100 \
	-genomeBuild GRCh37 \
	-bedtools /opt/tools/bedtools/bin/bedtools \
	-outputDir $OUTPUT_DIR \
	-outputFile ${OUTPUT_FILE}_full \
	-typeOfAnnotation full \
	-reciprocal yes \
	-overlap 70

