#!/bin/bash
#PBS -N ALD065_new_PrimaryPipeline
#PBS -V
#PBS -o /mnt/causes-data04/PROCESS/ALD_GENOME/ALD065/ALD065new.o
#PBS -e /mnt/causes-data04/PROCESS/ALD_GENOME/ALD065/ALD065new.e
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
## Set the total memory for the job
#PBS -l mem=30
## Set the max walltime for the job
#PBS -l walltime=200:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=8
NSLOTS=$PBS_NUM_PPN
umask 0002

SAMPLE_ID='ALD065_BWAmem'
WORKING_DIR='/mnt/causes-data04/PROCESS/ALD_GENOME/ALD065/'
FASTQR1='/mnt/causes-data04/RAW/ALD_GENOME/1707KHX-0054_hdd1/ALD065/ALD065_R1.fastq.gz'
FASTQR2='/mnt/causes-data04/RAW/ALD_GENOME/1707KHX-0054_hdd1/ALD065/ALD065_R2.fastq.gz'
BOWTIE2_INDEX='/mnt/causes-data01/data/GENOMES/hg19/hg19'
BWA_INDEX='/mnt/causes-data01/data/GENOMES/hg19/hg19_bwa'
GENOME_FASTA='/mnt/causes-data01/data/GENOMES/hg19/FASTA/hg19.fa'
CHROM='/mnt/causes-data01/data/GENOMES/hg19/FASTA/'
GENOMEFILE=/mnt/causes-data01/data/GENOMES/hg19/hg19_bwa.genome


 echo "Primary Analysis Started"
date

#FastQC
#/opt/tools/FastQC/fastqc --extract $FASTQR1 $FASTQR2 -o $WORKING_DIR
#
##Map with BWA
#/opt/tools/bwa-0.7.12/bwa mem $BWA_INDEX -t $NSLOTS -R "@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID\tPL:illumina" -M $FASTQR1 $FASTQR2 > $WORKING_DIR$SAMPLE_ID'.sam'
#
##Convert to binary, sort, and index
#/opt/tools/samtools-1.2/samtools view -@ $NSLOTS -u -bS $WORKING_DIR$SAMPLE_ID'.sam' | /opt/tools/samtools-1.2/samtools sort -@ $NSLOTS -m 3G - $WORKING_DIR$SAMPLE_ID'.sorted'
#/opt/tools/samtools-1.2/samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'
#
##Remove Duplicates
#TMPDIR=$WORKING_DIR'picardtmp/'
#mkdir $TMPDIR
#/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/picard-tools-1.139/picard.jar MarkDuplicates I=$WORKING_DIR$SAMPLE_ID'.sorted.bam' O=$WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam' REMOVE_DUPLICATES=false TMP_DIR=$TMPDIR M=$WORKING_DIR$SAMPLE_ID'_DuplicateResults.txt'
#/opt/tools/samtools-1.2/samtools index $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam'
#
##Realign
#/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt $NSLOTS -R $GENOME_FASTA -minReads 5 -I $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam' -o $WORKING_DIR$SAMPLE_ID'_indelsites.intervals' 
#/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T IndelRealigner -model USE_READS -R $GENOME_FASTA -targetIntervals $WORKING_DIR$SAMPLE_ID'_indelsites.intervals' -I $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam' -o $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam'
#/opt/tools/samtools-1.2/samtools index $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam'
#
##Clean Up
#ValidateOutput=`cat $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam_PicardSAMValidate.txt'`
#rm $WORKING_DIR$SAMPLE_ID'.sam'
#rm $WORKING_DIR$SAMPLE_ID'.sorted.bam'
#rm $WORKING_DIR$SAMPLE_ID'_dupremoved.sorted.bam'
#zip $WORKING_DIR$SAMPLE_ID'_R1_chastitypassed_fastqc.zip' $WORKING_DIR$SAMPLE_ID'_R1_chastitypassed_fastqc'
#rm $WORKING_DIR$SAMPLE_ID'_R1_chastitypassed_fastqc'
#zip $WORKING_DIR$SAMPLE_ID'_R2_chastitypassed_fastqc.zip' $WORKING_DIR$SAMPLE_ID'_R2_chastitypassed_fastqc'
#rm $WORKING_DIR$SAMPLE_ID'_R2_chastitypassed_fastqc'
#
##Read Metrics
#/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/picard-tools-1.139/picard.jar CollectMultipleMetrics INPUT=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' OUTPUT=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted' REFERENCE_SEQUENCE=$GENOME_FASTA
#
##SNP Calling HaplotypeCaller
#/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T HaplotypeCaller \
#	-nct 4 --never_trim_vcf_format_field \
#	--genotyping_mode DISCOVERY \
#	--standard_min_confidence_threshold_for_calling 10 \
# --standard_min_confidence_threshold_for_emitting 10 \
#	--min_mapping_quality_score 0 \
#	--min_base_quality_score 10 \
#	--minReadsPerAlignmentStart 5 \
# --minPruning 2 \
#	--pcr_indel_model NONE \
#	--dbsnp /opt/tools/GATK-3.5-0-g36282e4/resources/dbsnp_138.b37.excluding_sites_after_129.vcf \
#	-R $GENOME_FASTA \
#	-I $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' \
#	-o $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_HaplotypeCaller.vcf'

#SNP Calling HaplotypeCaller GVCFmode
/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar \
-T HaplotypeCaller -nct 4 --emitRefConfidence GVCF -R $GENOME_FASTA -I $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' -o $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_HaplotypeCaller.g.vcf'

exit

#Picard Validate BAM
/opt/tools/jdk1.7.0_79/bin/java -jar  /opt/tools/picard-tools-1.139/picard.jar ValidateSamFile M=VERBOSE MO=10000000 I=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' O=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam_PicardSAMValidate.txt'

#GATK Depth Of Coverage
/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R $GENOME_FASTA -I $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' -o $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_Coverage'

 echo "Primary Analysis Finished"
date
