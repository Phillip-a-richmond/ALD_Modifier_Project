annotateD <- function(D, geneSet, species="human", geneFilter="hgnc_symbol", archive="www.ensembl.org",...){

  # AJ - 12092016
  
  # Load necessarry packages
  library(biomaRt)
  require(plyr)
  
  #### Annotation ####
  # This is the point to add the gene info (Entrez ID, Ensembl ID, description)
  #  geneSet <- rownames(D$counts)
  # or:
  #  geneSet <- sapply(strsplit(rownames(D$counts),split=":"),'[',1)
  #  transcripts <- sapply(strsplit(rownames(D$counts),split=":"),'[',3)
  #  rownames(D$counts) <- paste(sapply(rownames(total),function(x) strsplit(x,":")[[1]][1]),sapply(rownames(total),function(x) strsplit(x,":")[[1]][2]),sep=":")
  #  geneSet <- sapply(strsplit(geneSet,split="_dup"),'[',1)
  # The latter has been used in combination with Lifescope (i.e. SOLiD) pipeline. The 'transcript' column is dropped here...
  
  # https://support.bioconductor.org/p/74322/
  # Biomart doesn't provide the service any longer...
  # "www.ensembl.org" reflects the latest version
  # NOTE that this means the latest GRCh, GRCm, hg or mm as well and that this might be of importance when
  # using coordinates etc. Here we only use identifiers, so we should be save..
  #
  
  
  
  if (species == "human"){
    ensembl  <- useEnsembl(biomart = "ensembl", 
                            dataset = "hsapiens_gene_ensembl", 
                            mirror = "useast")
    #ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host=archive)
    # Write out the versions used...
    write.table(listMarts(ensembl),paste0("biomaRt.",species,".versions"), sep="\t", quote=FALSE)
    symbol = "hgnc_symbol"
  } else if (species == "mouse"){
    ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host=archive)
    # Write out the versions used...
    write.table(listMarts(ensembl),paste0("biomaRt.",species,".versions"), sep="\t", quote=FALSE)
    symbol = "mgi_symbol"
  } else {
    warning("\tEither choose \"human\" or \"mouse\", other species are not (yet?) implemented...")
  }
  
  # Get annotation
  annot.biomart = getBM(attributes = c(symbol, "entrezgene_id", "ensembl_gene_id","description"), filters = geneFilter, values = geneSet,mart = ensembl, uniqueRows=FALSE)

  saveRDS(annot.biomart,file="annot.biomart.rds")
  rm(ensembl)
  annot.biomart = readRDS("annot.biomart.rds")
  
  ## remove duplicated combinations of GeneSymbol, Entrez Gene ID and Ensembl Gene ID
  annot.biomart = subset(annot.biomart,!duplicated(annot.biomart[,c(symbol, "entrezgene_id", "ensembl_gene_id")]))
  
  ## concatenate probes mapping to multiple Entrez Gene IDs
  df = aggregate(annot.biomart,list(annot.biomart[,geneFilter]),function(x) paste(unlist(x),collapse="//"))
  
  # If entrezgene etc. is a concatenation of two (or more) the same ids, uniguefy and save, otherwise just leave it ..
  df$entrezgene_id = unlist(lapply(df$entrezgene_id,function(x){a = unique(unlist(strsplit(x,"//")));ifelse(length(a) == 1,a,x)}))
  df$ensembl_gene_id = unlist(lapply(df$ensembl_gene_id,function(x){a = unique(unlist(strsplit(x,"//")));ifelse(length(a) == 1,a,x)}))
  df$description = unlist(lapply(df$description,function(x){a = unique(unlist(strsplit(x,"//")));ifelse(length(a) == 1,a,x)}))
  df[,symbol] = unlist(lapply(df[,symbol],function(x){a = unique(unlist(strsplit(x,"//")));ifelse(length(a) == 1,a,paste(a,collapse="//"))}))

  idx <- which(colnames(df)==geneFilter)
  df <- df[,-c(idx)]
  colnames(df)[1] <- geneFilter
  
  genes <- as.data.frame(rownames(D$counts))
  genes <- cbind(genes,geneSet)
  colnames(genes) <- c("original",geneFilter)
  
  # Make sure to use 'join' from plyr package
  genesListNew <- plyr::join(genes,df,by=geneFilter,type="left")
  
  # Replace the old genenames with the new annotation information
  D$genes <- genesListNew
  
  # Clean
  rm(annot.biomart,genes,geneSet,genesListNew,df)
  
  # Return object
  return(D)
}
