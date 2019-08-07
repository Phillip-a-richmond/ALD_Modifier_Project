python /mnt/causes-vnx1/PIPELINES/AnnotateVariants/PipelineScripts/Bam2Gemini.py \
	-d /mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/DeepVariantv0.8/ \
	-P /mnt/causes-vnx2/TIDE/PROCESS/ALD_GENOME/VCFs/ALD.ped \
	-G hg19 \
	-E prichmond@cmmt.ubc.ca \
	-T Genome \
	-F ALD \
	-v VCF \
	-V ALD010_DeepVariant0.8.vcf,ALD011_DeepVariant0.8.vcf,ALD026_DeepVariant0.8.vcf,ALD027_DeepVariant0.8.vcf,ALD036_DeepVariant0.8.vcf,ALD041_DeepVariant0.8.vcf,ALD042_DeepVariant0.8.vcf,ALD048_DeepVariant0.8.vcf,ALD049_DeepVariant0.8.vcf,ALD058_DeepVariant0.8.vcf,ALD059_DeepVariant0.8.vcf,ALD065_DeepVariant0.8.vcf \
	-A /mnt/causes-vnx1/PIPELINES/AnnotateVariants/ \
	-D /mnt/causes-vnx1/DATABASES/ \
	-m 90G \
	-p 24 \


