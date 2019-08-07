# Check to load environment
source /opt/tools/hpcenv.sh

WORKING_DIR=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/ALD010/
BAM=ALD010_BWAmem_dupremoved_realigned.sorted.bam 
MELT_DIR=/opt/tools/MELTv2.1.5/
MEI_LIST=${MELT_DIR}/me_refs/1KGP_Hg19/mei_list.txt

java -jar ${MELT_DIR}MELT.jar Single \
	-a \
	-b hs37d5/NC_007605 \
	-c 8 \
	-h /mnt/causes-vnx1/GENOMES/hg19/FASTA/hg19.fa \
	-bamfile $WORKING_DIR$BAM \
	-n ${MELT_DIR}add_bed_files/1KGP_Hg19/hg19.genes.bed \
	-t $MEI_LIST \
	-w $WORKING_DIR
