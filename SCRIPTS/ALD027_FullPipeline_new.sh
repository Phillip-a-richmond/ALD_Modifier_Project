#!/bin/bash
#PBS -N ALD027_new_PrimaryPipeline
#PBS -V
#PBS -o /mnt/causes-data04/PROCESS/ALD_GENOME/ALD027/ALD027new.o
#PBS -e /mnt/causes-data04/PROCESS/ALD_GENOME/ALD027/ALD027new.e
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

SAMPLE_ID='ALD027_BWAmem'
WORKING_DIR='/mnt/causes-data04/PROCESS/ALD_GENOME/ALD027/'
FASTQR1='/mnt/causes-data04/RAW/ALD_GENOME/1707KHX-0054_hdd1/ALD027/ALD027_R1.fastq.gz'
FASTQR2='/mnt/causes-data04/RAW/ALD_GENOME/1707KHX-0054_hdd1/ALD027/ALD027_R2.fastq.gz'
BOWTIE2_INDEX='/mnt/causes-data01/data/GENOMES/hg19/hg19'
BWA_INDEX='/mnt/causes-data01/data/GENOMES/hg19/hg19_bwa'
GENOME_FASTA='/mnt/causes-data01/data/GENOMES/hg19/FASTA/hg19.fa'
CHROM='/mnt/causes-data01/data/GENOMES/hg19/FASTA/'
GENOMEFILE=/mnt/causes-data01/data/GENOMES/hg19/hg19_bwa.genome


# echo "Primary Analysis Started"
#date
#
##FastQC
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
#	--genotypeing_mode DISCOVERY \
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
#
##SNP Calling HaplotypeCaller GVCFmode
#/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar \
#-T HaplotypeCaller -nct 4 --emitRefConfidence GVCF -R $GENOME_FASTA -I $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' -o $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_HaplotypeCaller.g.vcf'
#
##Picard Validate BAM
#/opt/tools/jdk1.7.0_79/bin/java -jar  /opt/tools/picard-tools-1.139/picard.jar ValidateSamFile M=VERBOSE MO=10000000 I=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' O=$WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam_PicardSAMValidate.txt'
#
##GATK Depth Of Coverage
#/opt/tools/jdk1.7.0_79/bin/java -jar /opt/tools/GATK-3.4-46/GenomeAnalysisTK.jar -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R $GENOME_FASTA -I $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned.sorted.bam' -o $WORKING_DIR$SAMPLE_ID'_dupremoved_realigned_Coverage'
#
# echo "Primary Analysis Finished"
#date
#
#echo "CNV Analysis Started"
#date
##Running CNVNATOR Windowsize 100
#
##Defining Variables
#BAM=${SAMPLE_ID}_dupremoved_realigned.sorted.bam 
#WIN=100
#/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR${BAM}.root -genome $GENOME_FASTA -tree $WORKING_DIR$BAM 
#/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -his $WIN -d $CHROM 
#/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -stat $WIN 
#/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -partition $WIN 
#/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -call $WIN > $WORKING_DIR/${SAMPLE_ID}_CNVnatorCall_$WIN 
## Preparing for use within Lumpy
#/opt/tools/lumpy-0.2.11/scripts/cnvanator_to_bedpes.py -c $WORKING_DIR${SAMPLE_ID}_CNVnatorCall_$WIN -b 600 --del_o $WORKING_DIR${SAMPLE_ID}.del.$WIN.bedpe --dup_o $WORKING_DIR${SAMPLE_ID}.dup.$WIN.bedpe 
#/opt/tools/lumpy-0.2.11/scripts/bedpe_sort.py -b $WORKING_DIR${SAMPLE_ID}.del.$WIN.bedpe -g $GENOMEFILE > $WORKING_DIR${SAMPLE_ID}.del.$WIN.sorted.bedpe 
#/opt/tools/lumpy-0.2.11/scripts/bedpe_sort.py -b $WORKING_DIR${SAMPLE_ID}.dup.$WIN.bedpe -g $GENOMEFILE > $WORKING_DIR${SAMPLE_ID}.dup.$WIN.sorted.bedpe 
##Running CNVNATOR Windowsize 1000
#
#WIN=1000
#/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -his $WIN -d $CHROM 
#/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -stat $WIN 
#/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -partition $WIN 
#/opt/tools/CNVnator/src/cnvnator -root $WORKING_DIR/${BAM}.root -genome $GENOME_FASTA -call $WIN > $WORKING_DIR/${SAMPLE_ID}_CNVnatorCall_$WIN 
#
#
#Lumpy CNV Caller

##PreProcess
BAM_FILE=$SAMPLE_ID'_dupremoved_realigned'
/opt/tools/samtools-1.2//samtools view -@ $NSLOTS -b -F 1294 $WORKING_DIR$BAM_FILE'.sorted.bam' > $WORKING_DIR$BAM_FILE'_discordants.bam'
/opt/tools/samtools-1.2/samtools view -@ $NSLOTS -h $WORKING_DIR$BAM_FILE'.sorted.bam' | python /opt/tools/lumpy/extractSplitReads_BwaMem.py -i stdin | /opt/tools/samtools-1.2/samtools view -@ $NSLOTS -Sb - > $WORKING_DIR$BAM_FILE'_splitters.bam'

##Sort the splitters and discordants
/opt/tools/samtools-1.2/samtools sort -@ $NSLOTS $WORKING_DIR$BAM_FILE'_discordants.bam' $WORKING_DIR$BAM_FILE'_discordants.sorted'
/opt/tools/samtools-1.2/samtools sort -@ $NSLOTS $WORKING_DIR$BAM_FILE'_splitters.bam' $WORKING_DIR$BAM_FILE'_splitters.sorted'

##Generate Empirical insert size stats
/opt/tools/samtools-1.2/samtools view -r $SAMPLE_ID $WORKING_DIR$BAM_FILE'.sorted.bam' | tail -n+1000000 | python /opt/tools/lumpy/pairend_distro.py -r 150 -X 4 -N 1000000 -o $WORKING_DIR$BAM_FILE'.histo' > $WORKING_DIR$BAM_FILE'.insertStats'
MEAN=`cat $WORKING_DIR$BAM_FILE'.insertStats' | sed -E 's/\s+/,/' | cut -d, -f1 | sed -E 's/mean://' | xargs printf "%.0f"`
STDEV=`cat $WORKING_DIR$BAM_FILE'.insertStats' | sed -E 's/\s+/,/' | cut -d, -f2 | sed -E 's/stdev://' | xargs printf "%.0f"`
echo "$STDEV"
echo "$MEAN"

##Run LUMPY
/opt/tools/lumpy/lumpy -mw 4 -tt 1.0 \
-bedpe bedpe_file:$WORKING_DIR${SAMPLE_ID}.dup.100.sorted.bedpe,weight:3,id:DUP \
-bedpe bedpe_file:$WORKING_DIR${SAMPLE_ID}.del.100.sorted.bedpe,weight:3,id:DEL \
 -pe  bam_file:$WORKING_DIR$BAM_FILE'.sorted.bam',id:$SAMPLE_ID,histo_file:$WORKING_DIR$BAM_FILE'.histo',read_length:150,min_non_overlap:100,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:40,mean:$MEAN,stdev:$STDEV \
 -sr id:$SAMPLE_ID,bam_file:$WORKING_DIR$BAM_FILE'_splitters.sorted.bam',back_distance:10,weight:1,min_mapping_threshold:40 > $WORKING_DIR$SAMPLE_ID'_Lumpy.vcf'

#Genotype,Filter,Annotate
##Run SVtyper
/opt/tools/svtyper-0.1.0/svtyper -B $WORKING_DIR$BAM_FILE'.sorted.bam' -i $WORKING_DIR$SAMPLE_ID'_Lumpy.vcf' -o $WORKING_DIR$SAMPLE_ID'_Lumpy_genotyped.vcf'

##Filter genotyped Lumpy (Custom script)
python /opt/tools/VariantAnnotation/FilterLumpyGenotyped.py $WORKING_DIR$SAMPLE_ID'_Lumpy_genotyped.vcf'\
 $WORKING_DIR$SAMPLE_ID'_Lumpy_genotyped_filtered.vcf' $WORKING_DIR$SAMPLE_ID'_Lumpy_DelDup.bed' \
 $WORKING_DIR$SAMPLE_ID'_Lumpy_DelDup.vcf' 0 3000000 

##Convert to ANNOVAR AVinput
perl /opt/tools/annovar/convert2annovar.pl -format vcf4 -allsample --comment --includeinfo -withfreq $WORKING_DIR$SAMPLE_ID'_Lumpy_genotyped_filtered.vcf' > $WORKING_DIR$SAMPLE_ID'_Lumpy.avinput'

##Run ANNOVAR DGV
/mnt/causes-data01/data/Databases/annovar/new_table_annovar.pl --otherinfo --buildver hg19 --protocol refgene,dgvMerged_commafix --operation g,r --argument '','-minqueryfrac 0.8 --colsWanted 2&3&4&10&17&18&19&21' \
$WORKING_DIR$SAMPLE_ID'_Lumpy.avinput' /mnt/causes-data01/data/Databases/annovar/humandb -out $WORKING_DIR$SAMPLE_ID'_Lumpy_annovar' 
python /opt/tools/VariantAnnotation/PruneAnnotatedCNVs.py  $WORKING_DIR$SAMPLE_ID'_Lumpy_annovar.hg19_multianno.txt' $WORKING_DIR$SAMPLE_ID'_Lumpy_annovar.hg19_multianno.filtered.txt' 
python /opt/tools/VariantAnnotation/AddSummaryToAnnovar_CNV.py $WORKING_DIR$SAMPLE_ID'_Lumpy_annovar.hg19_multianno.filtered.txt' $WORKING_DIR$SAMPLE_ID'_Lumpy_annovar.hg19_multianno.filtered.omimAndSummary.txt' 
# Running Pindel 

#Defining Variables
BAM=${SAMPLE_ID}_dupremoved_realigned.sorted.bam 

##Generate Empirical insert size stats
/opt/tools/samtools-1.2/samtools view -r $SAMPLE_ID $WORKING_DIR$BAM | tail -n+1000000 | python /opt/tools/lumpy/pairend_distro.py -r 150 -X 4 -N 1000000 -o $WORKING_DIR${BAM}.histo > $WORKING_DIR${BAM}.insertStats
MEAN=`cat $WORKING_DIR${BAM}.insertStats | sed -E 's/\s+/,/' | cut -d, -f1 | sed -E 's/mean://' | xargs printf "%.0f"`
python /mnt/causes-data01/data/PipelineControl/processingpipeline/pindel_config.py $WORKING_DIR$BAM $MEAN $WORKING_DIR$SAMPLE_ID 
/opt/tools/pindel-0.2.5b6/pindel --number_of_threads $NSLOTS \
-f $GENOME_FASTA \
-i $WORKING_DIR${BAM}_config.txt \
-c ALL -o $WORKING_DIR$SAMPLE_ID \
-M 4 \
-N -x 3 \
-I false -t false -r false \
-J /mnt/causes-data01/data/Databases/hg19_centromeres_telomeres.bed  


# MetaSV 

BAM=${SAMPLE_ID}_dupremoved_realigned.sorted.bam 
run_metasv.py --bam $WORKING_DIR/$BAM --reference $GENOME_FASTA \
--enable_per_tools_output \
--pindel_native $WORKING_DIR/${SAMPLE_ID}_D $WORKING_DIR/${SAMPLE_ID}_TD \
--lumpy_vcf $WORKING_DIR${SAMPLE_ID}_Lumpy.vcf \
--cnvnator_native ${SAMPLE_ID}_CNVnatorCall_1000 \
--filter_gaps --sample $WORKING_DIR$SAMPLE_ID --num_threads $NSLOTS \
--outdir $WORKING_DIR${SAMPLE_ID}MetaSVout \
--workdir $WORKING_DIR${SAMPLE_ID}MetaSVwork \
--spades /opt/tools/SPAdes-3.10.1/bin/spades.py \
--age /opt/tools/AGE/age_align \

echo "CNV Analysis Finished"
date
