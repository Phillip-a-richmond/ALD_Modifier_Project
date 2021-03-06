source /opt/tools/hpcenv.sh
source activate DeepVariant


BIN_VERSION="0.7.2"
MODEL_VERSION="0.7.2"

MODEL_NAME="DeepVariant-inception_v3-${MODEL_VERSION}+data-wgs_standard"
MODEL_HTTP_DIR="https://storage.googleapis.com/deepvariant/models/DeepVariant/${MODEL_VERSION}/${MODEL_NAME}"
DATA_HTTP_DIR="https://storage.googleapis.com/deepvariant/quickstart-testdata"

HOME=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/
OUTPUT_DIR=${HOME}/DeepVariantVCFs/
REF='/mnt/causes-vnx1/GENOMES/hg19/hg19_bwa.fa'
MODEL="/mnt/causes-vnx1/PIPELINES/DeepVariant/${MODEL_NAME}/model.ckpt"
BAM=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/STRs/ALD010_BWAmem_dupremoved_realigned.sorted.bam
SAMPLEID=ALD010
# Make Examples
docker run \
  -v /mnt:/mnt \
  gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/make_examples \
  --mode calling   \
  --ref "${REF}"   \
  --reads "${BAM}" \
  --examples "${OUTPUT_DIR}/examples.tfrecord.gz"

LOGDIR=${HOME}/logs

# Call Variants
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/call_variants_output.tfrecord.gz"

docker run \
  -v ${HOME}:${HOME} \
  gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/call_variants \
 --outfile "${CALL_VARIANTS_OUTPUT}" \
 --examples "${OUTPUT_DIR}/examples.tfrecord.gz" \
 --checkpoint "${MODEL}"

# Postprocess variants
FINAL_OUTPUT_VCF="${OUTPUT_DIR}/${SAMPLE_ID}_DeepVariant.vcf.gz"

docker run \
  -v ${HOME}:${HOME} \
  gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/postprocess_variants \
  --ref "${REF}" \
  --infile "${CALL_VARIANTS_OUTPUT}" \
  --outfile "${FINAL_OUTPUT_VCF}"







