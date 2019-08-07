#!/bin/bash
#PBS -N STRCalling_ALD
#PBS -V
#PBS -o /mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/$PBS_JOBNAME.$PBS_JOBID.o
#PBS -e /mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/$PBS_JOBNAME.$PBS_JOBID.e
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
## Set the total memory for the job
#PBS -l mem=30G
## Set the max walltime for the job
#PBS -l walltime=500:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=8
# The job array portion. Since I have 10 files, I'll run from 0-9 (indexes of the array)
#PBS -t 0-11%12

source /opt/tools/hpcenv.sh

# This script has been designed to run analysis on the entire set of ALD samples. By keeping the same structure in the header
# you are able to run some analysis on all the BAM files

# Define some variables relevant to your analysis
WORKDIR=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/
REF='/mnt/causes-vnx1/GENOMES/hg19/FASTA/hg19.fa'
CHROM='/mnt/causes-vnx1/GENOMES/hg19/FASTA/'
GENOME_FASTA='/mnt/causes-vnx1/GENOMES/hg19/FASTA/hg19.fa'

# Set the ALD directories into a list
ALD_DIRS=(${WORKDIR}ALD*)
echo ${ALD_DIRS[@]}

# These are relevant directories
BAMDIR=${ALD_DIRS[$MOAB_JOBARRAYINDEX]}/
WORKDIR=${ALD_DIRS[$MOAB_JOBARRAYINDEX]}/

# Change into the working directory
cd $WORKDIR

# Extract the BAM and SAMPLE identifiers
IFS='/' read -a array <<< "${ALD_DIRS[$MOAB_JOBARRAYINDEX]}"
SAMPLE=${array[-1]}
echo $SAMPLE
BAM=${SAMPLE}_BWAmem_dupremoved_realigned.sorted.bam
ls -lahtr $BAM


# Run ExpansionHunterDN
/mnt/causes-vnx1/PIPELINES/ExpansionHunter_DeNovo/ExpansionHunterDenovo-v0.6.2 --bam $BAMDIR$BAM \
        --reference $GENOME_FASTA \
        --output $WORKDIR${SAMPLE}_ExpansionHunterDenovo.json



