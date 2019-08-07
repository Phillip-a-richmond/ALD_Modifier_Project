#!/bin/bash
#PBS -N ALD_Bam2Gemini
#PBS -V
#PBS -o /mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/DeepVariantv0.8/ALD.o
#PBS -e /mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/DeepVariantv0.8/ALD.e
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
## Set the total memory for the job
#PBS -l mem=60G
## Set the max walltime for the job
#PBS -l walltime=240:00:00
## Set the total number of processors for the job
#PBS -l nodes=1:ppn=16
NSLOTS=$PBS_NUM_PPN
umask 0002
source /opt/tools/hpcenv.sh

# Define Tool paths. If they are in your path, simply change these full filepaths to only be the final command
# For example: Change BCFTOOLS=/opt/tools/bcftools-1.8/bin/bcftools to be BCFTOOLS=bcftools if it's in your path 

ANNOTVARDIR=/mnt/causes-vnx1/PIPELINES/AnnotateVariants/
SNPEFFJAR=/mnt/causes-vnx1/PIPELINES/SNPEff/snpEff/snpEff.jar
BCFTOOLS=/opt/tools/bcftools-1.8/bin/bcftools
VCFANNO=/opt/tools/vcfanno/vcfanno
VCF2DB=/opt/tools/vcf2db/vcf2db.py
GATKJAR=/opt/tools/GATK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
JAVA=/opt/tools/jdk1.8.0_92/bin/java
SNPEFFJAVA=/opt/tools/jdk1.8.0_92/bin/java
BGZIP=/opt/tools/tabix/bgzip
TABIX=/opt/tools/tabix/tabix
VT=/opt/tools/vt/vt

# Define variables for what you are working with

FAMILY_ID='ALD'
WORKING_DIR='/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/DeepVariantv0.8/'
GENOME_FASTA='/mnt/causes-vnx1/GENOMES/hg19/FASTA/hg19.fa'
PED_FILE=$WORKING_DIR/ALD.ped
TMPDIR=${WORKING_DIR}tmpdir/
mkdir $TMPDIR
GEMINIDB=$WORKING_DIR${FAMILY_ID}.db
VCF=$WORKING_DIR${FAMILY_ID}_DeepVariant_MergedVCFTools.vcf.gz
NORMVCF=$WORKING_DIR${FAMILY_ID}.deepvariant.merged.norm.vcf.gz
NORMFILTERVCF=$WORKING_DIR${FAMILY_ID}.deepvariant.merged.norm.filter.vcf.gz
ANNOVCF=$WORKING_DIR${FAMILY_ID}.deepvariant.merged.norm.vcfanno.vcf.gz 
SAMPLE1_VCF=ALD010_DeepVariant0.8.vcf
SAMPLE2_VCF=ALD011_DeepVariant0.8.vcf
SAMPLE3_VCF=ALD026_DeepVariant0.8.vcf
SAMPLE4_VCF=ALD027_DeepVariant0.8.vcf
SAMPLE5_VCF=ALD036_DeepVariant0.8.vcf
SAMPLE6_VCF=ALD041_DeepVariant0.8.vcf
SAMPLE7_VCF=ALD042_DeepVariant0.8.vcf
SAMPLE8_VCF=ALD048_DeepVariant0.8.vcf
SAMPLE9_VCF=ALD049_DeepVariant0.8.vcf
SAMPLE10_VCF=ALD058_DeepVariant0.8.vcf
SAMPLE11_VCF=ALD059_DeepVariant0.8.vcf
SAMPLE12_VCF=ALD065_DeepVariant0.8.vcf

# Merge VCFs

#$JAVA -Djava.io.tmpdir=$TMPDIR -jar $GATKJAR -T CombineVariants \
#-nt $NSLOTS \
#-R $GENOME_FASTA \
#--variant $WORKING_DIR${SAMPLE1_VCF} \
#--variant $WORKING_DIR${SAMPLE2_VCF} \
#--variant $WORKING_DIR${SAMPLE3_VCF} \
#--variant $WORKING_DIR${SAMPLE4_VCF} \
#--variant $WORKING_DIR${SAMPLE5_VCF} \
#--variant $WORKING_DIR${SAMPLE6_VCF} \
#--variant $WORKING_DIR${SAMPLE7_VCF} \
#--variant $WORKING_DIR${SAMPLE8_VCF} \
#--variant $WORKING_DIR${SAMPLE9_VCF} \
#--variant $WORKING_DIR${SAMPLE10_VCF} \
#--variant $WORKING_DIR${SAMPLE11_VCF} \
#--variant $WORKING_DIR${SAMPLE12_VCF} \
#-o $WORKING_DIR${FAMILY_ID}.deepvariant.merged.vcf 
#
#Get Rid of non-chr chromosomes

#  Normalize merged VCF, annotate with SNPeff

zless $VCF \
	| sed 's/ID=AD,Number=./ID=AD,Number=R/' \
	| $VT decompose -s - \
	| $VT normalize -r $GENOME_FASTA - \
	| $SNPEFFJAVA -Xmx10g -jar $SNPEFFJAR -dataDir /mnt/causes-vnx1/PIPELINES/SNPEff/snpEff/data/ GRCh37.87 \
	| $BGZIP -c > $NORMVCF 
$TABIX -p vcf $NORMVCF

# Filter Merged, normalized VCF

$BCFTOOLS filter \
	 --include 'FORMAT/AD[*:1]>=7 && FORMAT/DP[*] < 600' \
	 -m + \
	 -s + \
	 -O z \
	 --output $NORMFILTERVCF \
	 $NORMVCF 

$TABIX $NORMFILTERVCF \


# VCFAnno - Turn your VCF file into an annotated VCF file
$VCFANNO -lua $ANNOTVARDIR/VCFAnno/custom.lua \
-p $NSLOTS -base-path /mnt/causes-vnx1/DATABASES/ \
$ANNOTVARDIR/VCFAnno/VCFAnno_Config_20190321_GAC.toml \
$NORMFILTERVCF > $ANNOVCF 


# VCF2DB - Turn your annotated VCF file into a GEMINI DB

python $VCF2DB \
--expand gt_quals --expand gt_depths --expand gt_alt_depths --expand gt_ref_depths --expand gt_types \
 --a-ok InHouseDB_AC  --a-ok in_segdup --a-ok AF --a-ok AC --a-ok AN --a-ok MLEAC --a-ok MLEAF \
 --a-ok cpg_island --a-ok common_pathogenic --a-ok cse-hiseq --a-ok DS --a-ok ConfidentRegion \
--a-ok gnomad_exome_ac_global --a-ok gnomad_exome_ac_popmax --a-ok gnomad_exome_an_global --a-ok gnomad_exome_an_popmax --a-ok gnomad_exome_hom_controls --a-ok gnomad_exome_hom_global \
--a-ok gnomad_exome_hom_popmax --a-ok gnomad_exome_popmax --a-ok gnomad_genome_ac_global --a-ok gnomad_genome_ac_popmax --a-ok gnomad_genome_an_global --a-ok gnomad_genome_an_popmax \
--a-ok gnomad_genome_hom_controls --a-ok gnomad_genome_hom_global --a-ok gnomad_genome_hom_popmax --a-ok gnomad_genome_popmax \
$ANNOVCF $PED_FILE $GEMINIDB 
# Create Query Script within working directory 
rm ALD_GeminiQueryScript.header
rm ALD_GeminiQueryScript.sh
echo "WORKING_DIR=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/DeepVariantv0.8/" >> $WORKING_DIR/ALD_GeminiQueryScript.header
echo "GEMINIDB=ALD.db" >> $WORKING_DIR/ALD_GeminiQueryScript.header
echo "FAMILY_ID=ALD" >> $WORKING_DIR/ALD_GeminiQueryScript.header
echo "TableAnnotator=/mnt/causes-vnx1/PIPELINES/AnnotateVariants//TableAnnotators/GeminiTable2CVL.py" >> $WORKING_DIR/ALD_GeminiQueryScript.header
cat $WORKING_DIR/ALD_GeminiQueryScript.header /mnt/causes-vnx1/PIPELINES/AnnotateVariants//GeminiQueryScripts/GeminiQueries.sh > $WORKING_DIR/ALD_GeminiQueryScript.sh
#Run Query Script you just generated
sh $WORKING_DIR/ALD_GeminiQueryScript.sh

