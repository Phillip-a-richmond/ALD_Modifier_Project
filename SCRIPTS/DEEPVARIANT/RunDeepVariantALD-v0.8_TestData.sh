source /opt/tools/hpcenv.sh
#source activate DeepVariant

BIN_VERSION="0.8.0"

NSLOTS=4
HOME=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/
OUTPUT_DIR=${HOME}/DeepVariantv0.8/
GENOME_DIR=/mnt/causes-vnx1/GENOMES/hg19/
REF=hg19_bwa.fa
INPUT_DIR=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/
BAM=NA12878_S1.chr20.10_10p1mb.bam
SAMPLE_ID=NA12878_test
LOGDIR=${HOME}/logs
OUTPUT_VCF=${SAMPLE_ID}_DeepVariant0.8.vcf
OUTPUT_GVCF=${SAMPLE_ID}_DeepVariant0.8.g.vcf
N_SHARDS="4"

mkdir -p $OUTPUT_DIR

docker run \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}":"/output" \
  -v "${GENOME_DIR}":"/reference" \
  gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --regions "chr20:10,000,000-10,010,000" \
  --model_type=WGS \
  --ref="/reference/${REF}" \
  --reads="/input/${BAM}" \
  --output_vcf=/output/$OUTPUT_VCF \
  --output_gvcf=/output/$OUTPUT_GVCF \
  --num_shards=${N_SHARDS}



