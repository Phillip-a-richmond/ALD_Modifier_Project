# Phillip Richmond
# August 16th, 2018
# This is the workflow for generating the gene comparisons between discordant siblings

# Assumed start position: You have run the script: ALD_GeminiQueries.sh
# This will give you, from the gemini database, 4 files per family:
# Family1_HIGHMED_DominantProtective.tsv
# Family1_HIGHMED_RecessiveProtective.tsv
# Family1_HIGHMED_DominantDamaging.tsv
# Family1_HIGHMED_RecessiveDamaging.tsv


# 1) Now that you have your *tsv files, extract SVids with this command:
for file in *tsv
do
		cut -f1 $file > $file.SVid
done

# 2) Run intersections with intervene
# You need to have intervene installed for this, which you can do if you have miniconda installed
# Just run: pip install intervene
# Then make sure you have bedtools: conda install -c bioconda bedtools
# Then, run these commands:

# To generate the overlaps, and the upset diagrams
# ALL
intervene upset --ninter 22 --showzero --order degree --figsize 20 8 -i Family*DominantDamaging.tsv.SVid --save-overlaps --type list  --names Family1,Family2,Family3,Family4,Family5,Family6 --project All_VarLevel_SVid_DominantDamaging -o All_VarLevel_SVid_DominantDamaging_INTERVENE_UPSET
intervene upset --ninter 22 --showzero --order degree --figsize 20 8 -i Family*RecessiveDamaging.tsv.SVid --save-overlaps --type list --names Family1,Family2,Family3,Family4,Family5,Family6 --project All_VarLevel_SVid_RecessiveDamaging -o All_VarLevel_SVid_RecessiveDamaging_INTERVENE_UPSET
intervene upset --ninter 22 --showzero --order degree --figsize 20 8 -i Family*DominantProtective.tsv.SVid --save-overlaps --type list --names Family1,Family2,Family3,Family4,Family5,Family6 --project All_VarLevel_SVid_DominantProtective -o All_VarLevel_SVid_DominantProtective_INTERVENE_UPSET
intervene upset --ninter 22 --showzero --order degree --figsize 20 8 -i Family*RecessiveProtective.tsv.SVid --save-overlaps --type list --names Family1,Family2,Family3,Family4,Family5,Family6 --project All_VarLevel_SVid_RecessiveProtective -o All_VarLevel_SVid_RecessiveProtective_INTERVENE_UPSET

