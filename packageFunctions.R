library(dplyr)
library(DT)
library(data.table)
# Clear workspace
rm(list = ls())
# some helper functions
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
# end of helper functions



#' createMAE
createMAE.local <- function (exp, 
                       met, 
                       colData, 
                       sampleMap,
                       linearize.exp = FALSE,
                       filter.probes = NULL,
                       met.na.cut = 0.2,
                       filter.genes = NULL,
                       met.platform = "450K",
                       genome = NULL,
                       save = TRUE,
                       save.filename,
                       TCGA = FALSE) {
  
  if(missing(genome)) stop("Please specify the genome (hg38, hg19)")
  
  # Check if input are path to rda files
  if(is.character(exp)) exp <- get(load(exp))
  if(is.character(met)) met <- get(load(met))
  
  suppressMessages({
    
    if(!missing(colData)) { 
      if(is.character(colData)) { 
        colData <- as.data.frame(read_tsv(colData))
      }
      if (!"primary" %in% colnames(colData)) stop("No primary column in colData input")
      rownames(colData) <- colData$primary
    }
    if(!missing(sampleMap)) { 
      if(is.character(sampleMap)) sampleMap <- read_tsv(sampleMap)
      if (!all(c("assay","colname","primary") %in% colnames(sampleMap))) 
        stop("All assay, primary and colname columns should be in sampleMap input")
    }
  })  
  
  # Expression data must have the ensembl_gene_id (Ensemble ID) and external_gene_name (Gene Symbol)
  required.cols <- c("external_gene_name", "ensembl_gene_id")
  # If my input is a data frame we will need to add metadata information for the ELMER analysis steps
  if(!is(exp,"RangedSummarizedExperiment")){
    exp <- makeSummarizedExperimentFromGeneMatrix(exp, genome)
  }
  # Add this here ?
  if(linearize.exp) assay(exp) <- log2(assay(exp) + 1)
  
  
  if(!is(met,"RangedSummarizedExperiment")){
    met <- makeSummarizedExperimentFromDNAMethylation(met, genome, met.platform)
  }
  met <- met[rowMeans(is.na(assay(met))) < met.na.cut, ]
  
  # Select the regions from DNA methylation that overlaps enhancer.
  if(!is.null(filter.probes)){
    if(is.character(filter.probes)){
      filter.probes <- get(load(filter.probes))
    }
  } 
  if(!is.null(filter.probes) & !is.null(met)){
    met <- met[rownames(met) %in% names(filter.probes),]
  }
  if(!is.null(filter.genes) & !is.null(exp)){
    exp <- exp[rownames(exp) %in% names(filter.genes),]
  } 
  
  # We will need to check if the fields that we need exists.
  # Otherwise we will need to create them
  if(!is(exp,"RangedSummarizedExperiment")){
    required.cols <- required.cols[!required.cols %in% colnames(values(exp))]
    if(length(required.cols) > 0) {
      gene.info <- get.GRCh(genome)
      colnames(gene.info)[grep("external_gene", colnames(gene.info))] <- "external_gene_name"
      if(all(grepl("ENSG",rownames(exp)))) {
        extra <- as.data.frame(gene.info[match(rownames(exp),gene.info$ensembl_gene_id),required.cols])
        colnames(extra) <- required.cols
        values(exp) <- cbind(values(exp),extra)
      } else {
        stop("Please the gene expression matrix should receive ENSEMBLE IDs")
      }
    }
  } 
  if(TCGA){
    message("Checking if samples have both DNA methylation and Gene expression and if they are in the same order...")
    # If it is not TCGA we will assure the sample has both DNA methylation and gene expression
    ID <- intersect(substr(colnames(met),1,16), substr(colnames(exp),1,16))
    
    # Get only samples with both DNA methylation and Gene expression
    met <- met[,match(ID,substr(colnames(met),1,16))]
    exp <- exp[,match(ID,substr(colnames(exp),1,16))]
    stopifnot(all(substr(colnames(exp),1,16) == substr(colnames(met),1,16)))
    stopifnot(ncol(exp) == ncol(met))
    
    # Get clinical information
    if(missing(colData)) {
      colData <- TCGAbiolinks::colDataPrepare(colnames(met))
      # This will keep the same strategy the old ELMER version used:
      # Every type of tumor samples (starts with T) will be set to tumor and
      # every type of normal samples   (starts with N) will be set to normal 
      # See : https://github.com/lijingya/ELMER/blob/3e050462aa41c8f542530ccddc8fa607207faf88/R/Small.R#L8-L48
      colData$TN <- NA
      colData[grep("^N",colData$shortLetterCode),"TN"] <- "Normal" 
      colData[grep("^T",colData$shortLetterCode),"TN"] <- "Tumor" 
      
      colData$barcode <- NULL
      colData <- colData[!duplicated(colData),]      
      rownames(colData) <- colData$sample
    } 
    if(missing(sampleMap)) {
      sampleMap <- DataFrame(assay = c(rep("DNA methylation", length(colnames(met))), rep("Gene expression", length(colnames(exp)))),
                             primary = substr(c(colnames(met),colnames(exp)),1,16),
                             colname = c(colnames(met),colnames(exp)))
    }
    
    message("Creating MultiAssayExperiment")
    mae <- MultiAssayExperiment(experiments=list("DNA methylation" = met,
                                                 "Gene expression" = exp),
                                colData = colData,   
                                sampleMap = sampleMap,
                                metadata = list(TCGA= TRUE, genome = genome, met.platform = met.platform ))
  } else {
    
    if(missing(colData)){
      message <- paste("Please set colData argument. A data frame with samples", 
                       "information. All rownames should be colnames of DNA",
                       "methylation and gene expression. An example is showed",
                       "in MultiAssayExperiment documentation",
                       "(access it with ?MultiAssayExperiment)")
      stop(message)
    }
    
    if(missing(sampleMap)){
      # Check that we have the same number of samples
      message("Removing samples not found in both DNA methylation and gene expression (we are considering the names of the gene expression and DNA methylation columns to be the same) ")
      ID <- intersect(colnames(met), colnames(exp))
      met <- met[,match(ID,colnames(met))]
      exp <- exp[,match(ID,colnames(exp))]
      
      if(!all(colnames(exp) == colnames(met))) 
        stop("Error DNA methylation matrix and gene expression matrix are not in the same order")
      
      colData <- colData[match(ID,rownames(colData)),,drop = FALSE]
      sampleMap <- DataFrame(assay= c(rep("DNA methylation", length(colnames(met))), 
                                      rep("Gene expression", length(colnames(exp)))),
                             primary = c(colnames(met),colnames(exp)),
                             colname=c(colnames(met),colnames(exp)))
      mae <- MultiAssayExperiment(experiments=list("DNA methylation" = met,
                                                   "Gene expression" = exp),
                                  colData = colData,
                                  sampleMap = sampleMap,
                                  metadata = list(TCGA=FALSE, genome = genome, met.platform = met.platform ))
    } else {
      # Check that we have the same number of samples
      if(!all(c("primary","colname") %in% colnames(sampleMap))) 
        stop("sampleMap should have the following columns: primary (sample ID) and colname(DNA methylation and gene expression sample [same as the colnames of the matrix])")
      #if(!any(rownames(colData) %in% sampleMap$primary))
      #  stop("colData row names should be mapped to sampleMap primary column ")
      # Find which samples are DNA methylation and gene expression
      sampleMap.met <- sampleMap[sampleMap$assay %in% "DNA methylation",,drop = FALSE]
      sampleMap.exp <- sampleMap[sampleMap$assay %in% "Gene expression",,drop = FALSE]
      
      # Which ones have both DNA methylation and gene expression ?
      commun.samples <- intersect(sampleMap.met$primary,sampleMap.exp$primary)
      
      # Remove the one that does not have both data
      sampleMap.met <- sampleMap.met[match(sampleMap.met$primary,commun.samples),,drop = FALSE]
      sampleMap.exp <- sampleMap.exp[match(sampleMap.exp$primary,commun.samples),,drop = FALSE]
      
      # Ordering samples to be matched
      met <- met[,sampleMap.met$colname,drop = FALSE]
      exp <- exp[,sampleMap.exp$colname,drop = FALSE]
      
      if(!all(sampleMap.met$primary == sampleMap.exp$primary)) 
        stop("Error DNA methylation matrix and gene expression matrix are not in the same order")
      
      colData <- colData[match(commun.samples,colData$primary),,drop = FALSE]
      sampleMap <- DataFrame(assay= c(rep("DNA methylation", length(colnames(met))), 
                                      rep("Gene expression", length(colnames(exp)))),
                             primary = commun.samples,
                             colname=c(colnames(met),colnames(exp)))
      mae <- MultiAssayExperiment(experiments=list("DNA methylation" = met,
                                                   "Gene expression" = exp),
                                  colData = colData,
                                  sampleMap = sampleMap,
                                  metadata = list(TCGA=FALSE, genome = genome, met.platform = met.platform ))
    }
  }
  if(save) {
    if(missing(save.filename)) save.filename <- paste0("mae_",genome,"_",met.platform,".rda")
    save(mae, file = save.filename,compress = "xz")
    message("MAE saved as: ", save.filename)
  }
  return(mae)
}

# From: 5.calculate similarities

# convertEnsembleID2EntrezID <- function(ensemblgenes){
#   
#   exists_index <- sapply(ensemblgenes, function(x) exists(x, org.Hs.egENSEMBL2EG))
#   existed_genes <- ensemblgenes[exists_index]
#   ensemble2eg <- as.list(org.Hs.egENSEMBL2EG)
#   existed_genes_entrezIDs <- ensemble2eg[existed_genes]
#   return(existed_genes_entrezIDs)
# }

calculateGeneSimilaritiesBasedOnImmuneEscapeGenes <- function(geneList,geneFilePath){
  tmp.path<-getwd()
  setwd(geneFilePath) #setwd(immuneGenesFilePath)
  # Convert symbols to entrez ids
  idMapping <- as.list(org.Hs.egALIAS2EG)
  idMapping <- idMapping[!is.na(idMapping)]
  protoImmuneSymbols.old <- c('IL10', 'TGFB1', 'PTGER2','PTGER4', 'PTGER3', 'PTGER1', 'VEGFA',
                          'CTLA4', 'IDO1','IDO2', 'ICAM1', 'CD86', 'CD80', 'CAMP', 'FOXP3',
                          'IL2', 'IL12B','IL12A' , 'IFNG', 'SELL', 'CXCR1', 'CXCR2', 'ARG1',
                          'HLA-E', 'HLA-G', 'FASLG', 'CD274', 'TNFRSF6B', 'TNFSF15','TNFSF14')
  protoImmuneSymbols.file<-read.csv("all_human_immune-escape_searched.csv",stringsAsFactors = F)
  protoImmuneSymbols<-protoImmuneSymbols.file$Gene
  protoImmuneSymbols<-unique(c(protoImmuneSymbols,protoImmuneSymbols.old))
  print(protoImmuneSymbols)

  #idvals <- idMapping[protoImmuneSymbols.old]
  idvals <- idMapping[protoImmuneSymbols]
  protoImmune <- unlist(idvals, use.names=FALSE)
  geneIDVector <- rapply(geneList,c)
  simDF <- getGeneSimPrototypes(geneIDVector, prototypes = protoImmune, similarity = "max",
                                similarityTerm = "JiangConrath", pca = FALSE, # pca = TRUE
                                normalization = TRUE, verbose = TRUE)
  # Change Nan to 0
  is.nan.data.frame <- function(x) {
    do.call(cbind, lapply(x, is.nan)) 
  }
  simil <- simDF$similarity
  simil[is.nan(simil)] <- 0
  setwd(tmp.path) # back to default path before go...
  return(simil)
}

# From 6: Clustering
convertEnsembleID2EntrezID4DEGenes <- function(DEGenes){
  DEGenes  <- DEGenes %>% rownames_to_column('gene')
  ensemblegenes <- DEGenes$gene
  exists_index <- sapply(ensemblegenes, function(x) exists(x, org.Hs.egENSEMBL2EG))
  existed_genes <- ensemblegenes[exists_index]
  DEGenes$exists <- exists_index
  DEGenes <- DEGenes[DEGenes$exists == TRUE,]
  ensemble2eg <- as.list(org.Hs.egENSEMBL2EG)
  DEGenes$gene <- ensemble2eg[existed_genes]
  return(DEGenes)
}

clusterGenes <- function(savePath,tumorType,conditionType,similarityMatrix){
  tmp.path<-getwd()
  setwd(savePath)
  if(sum(dim(similarityMatrix))>4){
    h=hclust(as.dist(1-similarityMatrix),"average")
    # check that is there a folder to save the png files if not create one
    if (file.exists(tumorType)){
      cat("=========> ",tumorType, "folder exists...","\n")
    } else {
      dir.create(file.path(clusterResultsPath, tumorType))
      cat("=========> ",tumorType, "folder created successfully...","\n")
    }
    png(paste0(savePath,tumorType,"/",conditionType,".meth.",tumorType,".cluster.png"))
    suppressWarnings(plot(h))
    dev.off()
    groups <- cutree(h, h=0.3) # 0.36 # 0.3 PPAR
    if(length(unique(groups))>=7){
      png(paste0(savePath,tumorType,"/",conditionType,".meth.",tumorType,".silhouette.png"))
      suppressWarnings(plot(silhouette(groups, dist=as.dist(similarityMatrix), main = "Silhouette")))
      dev.off()
      clusters <- c()
      for (i in 1:max(groups)){
        if(length(names(groups[groups == i])) >= 7){
          clusters[i] <- list(names(groups[groups == i]))    
        }
      }
      clusters <- Filter(Negate(is.null), clusters) # returns logical(0) if filter ends with nothing
      setwd(tmp.path)
      return(clusters)
    }else{
      cat("=========> ",tumorType, "silhoutte plot can't be plotted not enough groups...","\n")
      setwd(tmp.path)
      return(NA)
    }
  }else{
    cat("=========> ",tumorType, "has not enough similar genes so iteration skipped for this type...","\n")
    setwd(tmp.path)
    return(NA)
  }
}

# From 7: Pathfinder

