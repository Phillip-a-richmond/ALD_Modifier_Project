Filename=SummaryFile_ALD_Centrifuge.tsv
rm $Filename
for file in *tsv
do
	#echo $file
	IFS='_' read -a array <<< "$file"
	sampleID=${array[0]}
	echo "$sampleID:	Viral Analysis With Centrifuge" >> $Filename
	cat $file >> $Filename
	echo ""	>> $Filename
done

