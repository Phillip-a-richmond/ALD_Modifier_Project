#for vcf in $(ls *DeepVariant0.8.vcf.gz)
#do
#	echo $vcf
#	#bgzip $vcf
#	tabix $vcf
#done

bcftools merge

