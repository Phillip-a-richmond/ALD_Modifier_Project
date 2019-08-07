for VCF in $(ls *vcf)
do
	echo $VCF
	/opt/tools/IGVTools/igvtools index $VCF
done

