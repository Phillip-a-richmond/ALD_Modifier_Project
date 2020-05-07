library(edgeR)
library(readr)
#'
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


# create the DGEList objects
ywork = DGEList(counts=df_counts_match,group=group,genes=lengths)

# create the design
famwork = factor(fam)
phenowork = factor(pheno)
designwork = model.matrix(~ famwork + phenowork)

nrfam = dim(ywork$counts)[2]/2
keep<- rowSums(cpm(ywork)[,phenowork==1]>0.5)== nrfam | rowSums(cpm(ywork)[,phenowork==0]>0.5)==nrfam
dgList <- ywork[keep, , keep.lib.sizes=FALSE]
# estimate normalization factors
dgList <- calcNormFactors(dgList)
# estimate dispersion
dgList <- estimateDisp(dgList,designwork,robust=TRUE)

# estimate common and tagwise dispersion
dgList <- estimateCommonDisp(dgList)
dgList <- estimateTagwiseDisp(dgList)

# multiply the pseudocounts by the norm factors to obtain the TMM counts
norm_counts.table <- t(t(dgList$pseudo.counts)*(dgList$samples$norm.factors))
# write.csv(norm_counts.table,'tmm_norm_counts_5_2_2020.csv')
# use cpm to calcuate the cpm counts
cpm_counts.table<-cpm(dgList$counts)
# write.csv(cpm_counts.table,'cpm_norm_counts_5_2_2020.csv')
