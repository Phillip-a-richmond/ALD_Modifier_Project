#!/bin/bash
#PBS -N SVCalling_ALD
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
#PBS -t 10-12%2

source /opt/tools/hpcenv.sh

# This script has been designed to run analysis on the entire set of ALD samples. By keeping the same structure in the header
# you are able to run some analysis on all the BAM files

# Define some variables relevant to your analysis
WORKING_DIR=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/
REF='/mnt/causes-vnx1/GENOMES/hg19/FASTA/hg19.fa'
CHROM='/mnt/causes-vnx1/GENOMES/hg19/FASTA/'
GENOME_FASTA='/mnt/causes-vnx1/GENOMES/hg19/FASTA/hg19.fa'

# Set the ALD directories into a list
ALD_DIRS=(${WORKING_DIR}ALD*)
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

#Index bam file

/opt/tools/samtools-1.2/samtools index $BAM

#Extract discordant reads
/opt/tools/samtools-1.2/samtools view \
	-b \
	-F 1294 \
$BAMDIR/$BAM > $WORKDIR/${SAMPLE}.discordants.unsorted.bam

#Extract split reads
/opt/tools/samtools-1.2/samtools view -h $BAMDIR/$BAM \
	| /opt/tools/lumpy-0.2.11/scripts/extractSplitReads_BwaMem -i stdin \
	| /opt/tools/samtools-1.2/samtools view -Sb - > $WORKDIR/${SAMPLE}.splitters.unsorted.bam

##Sort alignments
/opt/tools/samtools-1.2/samtools sort $WORKDIR/${SAMPLE}.discordants.unsorted.bam $WORKDIR/${SAMPLE}.discordants
/opt/tools/samtools-1.2/samtools sort $WORKDIR/${SAMPLE}.splitters.unsorted.bam $WORKDIR/${SAMPLE}.splitters 

#Obtain insert size metrics
/opt/tools/samtools-1.2/samtools view $BAMDIR/$BAM \
	| tail -n+100000 \
	| /opt/tools/lumpy-0.2.11/scripts/pairend_distro.py \
	-r 150 \
	-X 4 \
	-N 10000 \
	-o $WORKDIR/${SAMPLE}.histo > $WORKDIR/${SAMPLE}.insertsize

INSMEAN=$(awk '{print $1}' $WORKDIR/${SAMPLE}.insertsize | cut -d ':' -f2)
INSSTD=$(awk '{print $2}' $WORKDIR/${SAMPLE}.insertsize | cut -d ':' -f2)

#Lumpy
/opt/tools/lumpy-0.2.11/lumpy \
-mw 4 \
-tt 0 \
-pe id:$WORKDIR/$SAMPLE,bam_file:$WORKDIR/${SAMPLE}.discordants.bam,histo_file:$WORKDIR/${SAMPLE}.histo,mean:$INSMEAN,stdev:$INSSTD,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,min_mapping_threshold:20 \
-sr id:$WORKDIR/$SAMPLE,bam_file:$WORKDIR/${SAMPLE}.splitters.bam,back_distance:20,weight:1,min_mapping_threshold:20 > $WORKDIR/${SAMPLE}.LUMPY.vcf

##CNVnator
/opt/tools/CNVnator/src/cnvnator -root $WORKDIR/${SAMPLE}.root -genome $REF -tree $BAMDIR/$BAM
/opt/tools/CNVnator/src/cnvnator -root $WORKDIR/${SAMPLE}.root -genome $REF -his 100 -d $CHROM
/opt/tools/CNVnator/src/cnvnator -root $WORKDIR/${SAMPLE}.root -genome $REF -stat 100
/opt/tools/CNVnator/src/cnvnator -root $WORKDIR/${SAMPLE}.root -genome $REF -partition 100
/opt/tools/CNVnator/src/cnvnator -root $WORKDIR/${SAMPLE}.root -genome $REF -call 100 > $WORKDIR/${SAMPLE}_CNVcall_100

##Pindel
python /opt/tools/pindel-0.2.5b6/pindel_config.py $BAMDIR/$BAM $INSMEAN $WORKDIR/$SAMPLE
/opt/tools/pindel-0.2.5b6/pindel \
        -f $REF \
        -i $BAMDIR/${BAM}_config.txt \
        -c ALL \
        -o $WORKDIR/$SAMPLE \
        --number_of_threads 12 \
        -M 4 \
        -N \
        -r false \
        -t false \
        -I false \
        -x 2 \
        -J /mnt/causes-vnx1/DATABASES/RLCRs_no_repeatMaster.bed


# Re-run HaplotypeCaller
/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T HaplotypeCaller \
       -nct 4 --never_trim_vcf_format_field \
       --genotyping_mode DISCOVERY \
       --standard_min_confidence_threshold_for_calling 10 \
 --standard_min_confidence_threshold_for_emitting 10 \
       --min_mapping_quality_score 0 \
       --min_base_quality_score 10 \
       --minReadsPerAlignmentStart 5 \
 --minPruning 2 \
       --pcr_indel_model NONE \
       --dbsnp /opt/tools/GATK-3.5-0-g36282e4/resources/dbsnp_138.b37.excluding_sites_after_129.vcf \
       -R $GENOME_FASTA \
       -I $BAMDIR$BAM \
       -o $BAMDIR$SAMPLE'_dupremoved_realigned_HaplotypeCaller.vcf'



# CNV Calling with ERDS
perl /opt/tools/erds1.1/erds_pipeline.pl \
        -b $BAMDIR/$BAM \
        -v $BAMDIR/${SAMPLE}_BWAmem_dupremoved_realigned_HaplotypeCaller.vcf \
        -o $WORKDIR/${SAMPLE}.ERDS/ \
        -r $GENOME_FASTA

##metaSV
run_metasv.py --bam $BAMDIR/$BAM --reference $REF \
        --lumpy_vcf $WORKDIR/${SAMPLE}.LUMPY.vcf \
        --cnvnator_native $WORKDIR/${SAMPLE}_CNVcall_100 \
        --pindel_native $WORKDIR/${SAMPLE}_D $WORKDIR/${SAMPLE}_SI \
        --filter_gaps \
        --sample $WORKDIR/$SAMPLE --num_threads 12 \
        --outdir $WORKDIR/${SAMPLE}Out \
        --workdir $WORKDIR/${SAMPLE}Work \
        --spades /opt/tools/SPAdes-3.10.1/bin/spades.py \
        --age /opt/tools/AGE/age_align \
        --svs_to_assemble INV \
        --minsvlen 50

