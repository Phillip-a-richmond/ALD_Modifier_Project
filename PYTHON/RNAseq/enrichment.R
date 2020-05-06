rm(list=ls(all=TRUE))

# use rstudioapi library
library(rstudioapi) 
library(readr)


# use edgeR library from bioconductor

#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")


library(edgeR)

# use homo sapiens db
library(org.Hs.eg.db)

# set path
current_path <- getActiveDocumentContext()$path 
# The next line set the working directory to the relevant one:
setwd(dirname(current_path ))


# load local annotation function
source("annotateD.R")

# load the original count data

count_data <- read_delim("../../../Data/RNA/count_data.txt","\t", escape_double = FALSE, trim_ws = TRUE)
lengths = count_data[,2]
dc = as.data.frame(count_data)
rownames(dc)<-dc[,1]
# remove genenames
dc = dc[,-c(1,2)]

# hard code the user ids, family labels and phenotypes
t <- c(10,11,26,27,36,41,42,48,49,58,59,65,
       1,0,1,0,1,0,1,0,1,1,0,0,
       1,1,2,2,3,4,4,5,5,6,6,3)

# reformat the sample info

sampleInfo=matrix(t,nrow=12,ncol=3,byrow=FALSE)
colnames(sampleInfo) <- c("ID", "CALD","Family")


# order sampleInfo by family, then CALD status
sampleInfo = sampleInfo[order(sampleInfo[,3],sampleInfo[,2]),]
# remove ensembleids and read length labels
samples = colnames(count_data)[-c(1,2)]
# sampleID contains order in count_data and dc table 
sampleID = as.numeric(gsub("ALD","",gsub(".bam","",samples)))

# sampleID[idmatch] == sampleInfo[,1]
idmatch = match(sampleInfo[,1],sampleID)

# copy dc table in sorted order
df_counts_match = dc[,idmatch]

# create the factors for the model matrix
fam = factor(sampleInfo[,3])
pheno = factor(sampleInfo[,2])
group = pheno

# show the table
table(pheno,fam)

# define the model design
design = model.matrix(~ fam + pheno)


# create the DGEList objects
y = DGEList(counts=df_counts_match,group=group,genes=lengths)

# find proper mapping ids and symbols for go/kegg analysis later...
y$genes$Symbol=mapIds(org.Hs.eg.db,rownames(df_counts_match),keytype = "ENSEMBL",column="SYMBOL")
y$genes$Entrez=mapIds(org.Hs.eg.db,rownames(df_counts_match),keytype = "ENSEMBL",column="ENTREZID")
y$genes <- rename(y$genes, c("Entrez"="entrezgene", "Symbol"="hgnc_symbol"))
y <- annotateD(y, y$genes$hgnc_symbol, geneFilter = "hgnc_symbol")


# calculate the p-values for all the families

ywork <-y
famwork = factor(fam)
phenowork = factor(pheno)
designwork = model.matrix(~ famwork + phenowork)

nrfam = 6
keep<- rowSums(cpm(ywork)[,phenowork==1]>0.5)== nrfam | rowSums(cpm(ywork)[,phenowork==0]>0.5)==nrfam
ywork<- ywork[keep, , keep.lib.sizes=FALSE]

# calculate the normalization factors
ywork_norm <- calcNormFactors(ywork) 
# estimate the dispersion paramters
ywork_norm <- estimateDisp(ywork_norm,designwork,robust=TRUE)

# quasi likelihood fitting
ywork_fit <- glmQLFit(ywork_norm,designwork,robust=TRUE) 
ywork_qlf2 <- glmQLFTest(ywork_fit,coef=7) 

# only select those that have at least a pvalue < 0.05

id_sig_0.05 = ywork_qlf2$table$PValue<0.05
mygeneswork = ywork_qlf2[id_sig_0.05,]
si=sort(mygeneswork$table$PValue,index.return=TRUE)
myvalswork = mygeneswork[si$ix,]$table
subsetwork = ywork_fit[id_sig_0.05,]
mygenes_sorted = subsetwork[si$ix,]


goterms = goana(de=mygenes_sorted$genes$entrezgene,Species="Hs")
topGO(goterms,ontology="BP")
topGO(goterms,ontology="MF")
topGO(goterms,ontology="CC")

go_data = topGO(goterms,ontology=c("BP","MF","CC"))
#wite.csv(go_data,'go_results.csv')


keggterms = kegga(de=mygenes_sorted$genes$entrezgene_id,Species="Hs")
kegg_data = topKEGG(keggterms)
print(kegg_data)
#write.csv(kegg_data,'kegg_results.csv')

