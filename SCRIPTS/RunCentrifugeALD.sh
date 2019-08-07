# Step 1 - Get Centrifuge binary
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/downloads/centrifuge-1.0.3-beta-Linux_x86_64.zip
unzip centrifuge-1.0.3-beta-Linux_x86_64.zip

# Step 2 - Get centrifuge library for bacteria + viruses + human
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p_compressed+h+v.tar.gz
tar -xzvf p_compressed+h+v.tar.gz


# Step 3 - Run centrifuge for all ALD RNA-seq Fastqs
declare -a SampleList=('ALD10' 'ALD11' 'ALD26' 'ALD27' 'ALD41' 'ALD42' 'ALD48' 'ALD49' 'ALD58' 'ALD59' 'ALD36' 'ALD65')

for sample in "${SampleList[@]}"
do
	echo $sample

	./centrifuge -p 16 -x p_compressed+h+v -1 /mnt/causes-vnx2/TIDE/RAW/ALD_GENOME/${sample}.dedupped_1.fastq.gz -2 /mnt/causes-vnx2/TIDE/RAW/ALD_GENOME/${sample}.dedupped_2.fastq.gz --report-file ${sample}_report.tsv > ${sample}_centrifugeresults
done

exit

