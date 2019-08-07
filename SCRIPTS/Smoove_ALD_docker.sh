# Running SMOOVE through the docker on hpc01
BAM1=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD010_BWAmem_dupremoved_realigned.sorted.bam
BAM2=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD011_BWAmem_dupremoved_realigned.sorted.bam
BAM3=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD026_BWAmem_dupremoved_realigned.sorted.bam
BAM4=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD027_BWAmem_dupremoved_realigned.sorted.bam
BAM5=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD036_BWAmem_dupremoved_realigned.sorted.bam
BAM6=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD041_BWAmem_dupremoved_realigned.sorted.bam
BAM7=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD042_BWAmem_dupremoved_realigned.sorted.bam
BAM8=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD048_BWAmem_dupremoved_realigned.sorted.bam
BAM9=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD049_BWAmem_dupremoved_realigned.sorted.bam
BAM10=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD058_BWAmem_dupremoved_realigned.sorted.bam
BAM11=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD059_BWAmem_dupremoved_realigned.sorted.bam
BAM12=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/ALD065_BWAmem_dupremoved_realigned.sorted.bam
NAME=ALD
OUTDIR=/mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/SVs/${NAME}_Smoove
PROCS=8
GENOME_FASTA=/mnt/causes-vnx1/GENOMES/hg19/hg19_bwa.fa

mkdir $OUTDIR


docker run -v /mnt:/mnt brentp/smoove smoove call \
	--name $NAME -p $PROCS -x \
	--outdir $OUTDIR \
	--genotype --duphold -f $GENOME_FASTA \
	$BAM1 $BAM2 $BAM3 $BAM4 $BAM5 $BAM6 $BAM7 $BAM8 $BAM9 $BAM10 $BAM11 $BAM12

