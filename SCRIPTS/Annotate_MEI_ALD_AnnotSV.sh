#!/bin/bash
#PBS -N MEI_Annotation_ALD
#PBS -V
#PBS -o /mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/MELT/Annotation.o
#PBS -2 /mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/MELT/Annotation.e
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

export ANNOTSV=/mnt/causes-vnx1/DATABASES/AnnotSV_2.2/
source /opt/tools/hpcenv.sh
source activate SVcalling


WORKDIR=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/MEI/MELT/
cd $WORKDIR

# ALU
INPUT_VCF=${WORKDIR}/ALU.final_comp.vcf
OUTPUT_DIR=$WORKDIR
OUTPUT_FILE=ALD_MEI_ALU_AnnotSV

$ANNOTSV/bin/AnnotSV/AnnotSV.tcl \
	-SVinputFile $INPUT_VCF \
	-SVinputInfo 0 \
	-SVminSize 100 \
	-genomeBuild GRCh37 \
	-bedtools /opt/tools/bedtools/bin/bedtools \
	-outputDir $OUTPUT_DIR \
	-outputFile ${OUTPUT_FILE}_split \
	-typeOfAnnotation split \

$ANNOTSV/bin/AnnotSV/AnnotSV.tcl\
	-SVinputFile $INPUT_VCF \
	-SVinputInfo 0 \
	-SVminSize 100 \
	-genomeBuild GRCh37 \
	-bedtools /opt/tools/bedtools/bin/bedtools \
	-outputDir $OUTPUT_DIR \
	-outputFile ${OUTPUT_FILE}_full \
	-typeOfAnnotation full \


# SVA
INPUT_VCF=${WORKDIR}/SVA.final_comp.vcf
OUTPUT_DIR=$WORKDIR
OUTPUT_FILE=ALD_MEI_SVA_AnnotSV

$ANNOTSV/bin/AnnotSV/AnnotSV.tcl\
        -SVinputFile $INPUT_VCF \
        -SVinputInfo 0 \
        -genomeBuild GRCh37 \
        -bedtools /opt/tools/bedtools/bin/bedtools \
        -outputDir $OUTPUT_DIR \
        -outputFile ${OUTPUT_FILE}_split \
        -typeOfAnnotation split \

$ANNOTSV/bin/AnnotSV/AnnotSV.tcl\
        -SVinputFile $INPUT_VCF \
        -SVinputInfo 0 \
        -genomeBuild GRCh37 \
        -bedtools /opt/tools/bedtools/bin/bedtools \
        -outputDir $OUTPUT_DIR \
        -outputFile ${OUTPUT_FILE}_full \
        -typeOfAnnotation full \

#LINE1

INPUT_VCF=${WORKDIR}/LINE1.final_comp.vcf
OUTPUT_DIR=$WORKDIR
OUTPUT_FILE=ALD_MEI_LINE1_AnnotSV

$ANNOTSV/bin/AnnotSV/AnnotSV.tcl\
        -SVinputFile $INPUT_VCF \
        -SVinputInfo 0 \
        -genomeBuild GRCh37 \
        -bedtools /opt/tools/bedtools/bin/bedtools \
        -outputDir $OUTPUT_DIR \
        -outputFile ${OUTPUT_FILE}_split \
        -typeOfAnnotation split \

$ANNOTSV/bin/AnnotSV/AnnotSV.tcl\
        -SVinputFile $INPUT_VCF \
        -SVinputInfo 0 \
        -genomeBuild GRCh37 \
        -bedtools /opt/tools/bedtools/bin/bedtools \
        -outputDir $OUTPUT_DIR \
        -outputFile ${OUTPUT_FILE}_full \
        -typeOfAnnotation full \







