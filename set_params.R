# Clear workspace
rm(list = ls())
library(TCGAbiolinks)
library(SummarizedExperiment)
library(miceadds)
#library(ELMER)
library(MultiAssayExperiment)
library(GOSim)
library("org.Hs.eg.db") #bioclite("org.Hs.eg.db")
library(SPIA)
library(cluster) # for silhouette
library(dplyr)
library(tibble)
#library(pathfindR)

library(annotate)

# job 0. Path Assignments
initPath<-"/home/burak/Desktop/new_meth_jobs2"
coreFuncPath<-"/home/burak/Desktop/new_meth_jobs2/core functions/pathFinder"
setwd(initPath)
source("packageFunctions.R")
source("ensembleIDtoGeneID.R") #convertEnsembleID2GeneID(my.genes,"hypo","blca")
setwd(coreFuncPath)
source("core_functions.R") #run_pathfindR.local
source("active_snw_functions.R") #active_snw_search.local
source("enrichment_functions.R") #enrichment_analyses.local
source("visualization_functions.R") #visualize_terms.local
setwd(initPath)

metFilesPath<-"/home/burak/Desktop/meth_pan_cancer/rdaMET_files"
expFilesPath<-"/home/burak/Desktop/meth_pan_cancer/rdaEXP_files"
maeFilePath<-"/home/burak/Desktop/meth_pan_cancer/outputs/rdaMAE_files/"
localFilesPath<-"/home/burak/Desktop/new_meth_jobs2/csv_files"
oldRunFilesPath<-"/home/burak/Desktop/meth/results"
sigDiffResultsPath<-"/home/burak/Desktop/meth_pan_cancer/outputs/results/sig_diff_meth/"
dEGresultsPath<-"/home/burak/Desktop/meth_pan_cancer/outputs/results/degs/"
immuneGenesFilePath<-"/home/burak/Desktop/meth_pan_cancer/"
clusterResultsPath<-"/home/burak/Desktop/meth_pan_cancer/outputs/results/clusters/"
pathFinderResultsPath<-"/home/burak/Desktop/meth_pan_cancer/outputs/results/pathfinders/"
rootResultsPath<-"/home/burak/Desktop/meth_pan_cancer/outputs/results/"

# Num 1:
convertToMAE<-function(cancerDF.csv,inputPath.met,inputPath.exp,outPath,localPath,oldRunsPath){
  # fnc job 0. change the wd to the local path
  setwd(localPath)
  outOfJob.files<-list.files(oldRunsPath)
  `%notin%` <- Negate(`%in%`) # r does not have default 'notin' operator so we just make one
  # fnc job 1. read the df and pull the cancer names
  tempDF<-read.csv(cancerDF.csv, header = T, stringsAsFactors = F)
  # fnc job 2. pull and shape the necessary data from the df
  for (cancer.type in tempDF$Project){
    filename.prefix<-tolower(sub(".*-", "", cancer.type))
    # fnc job 3. filter out previously worked cancer types
      if(filename.prefix %notin% outOfJob.files){
        
        print(filename.prefix)
        
        # fnc 4. convert to mae job definitions
        #load methylation data
        setwd(inputPath.met)
        #met.file<-load.Rdata(filename=paste0(filename.prefix, "DNAmet.rda"), paste0(filename.prefix, "DNAmet.rda"))
        met.file <- miceadds::load.Rdata2(filename=paste0(filename.prefix, "DNAmet.rda"))

        # na.omit
        met.file <- subset(met.file, subset = (rowSums(is.na(assay(met.file))) == 0))

        #load gene expression data
        setwd(inputPath.exp)
        #exp.file<-load.Rdata(filename=paste0(filename.prefix, "Exp.rda"), paste0(filename.prefix, "Exp.rda"))
        exp.file <- miceadds::load.Rdata2(filename=paste0(filename.prefix, "Exp.rda"))
        
        distal.probes <- get.feature.probe(genome = "hg38", met.platform = "450K")

        mae <- createMAE.local(exp = exp.file, 
                         met = met.file,
                         save = TRUE,
                         linearize.exp = TRUE,
                         filter.probes = distal.probes,
                         save.filename = paste0(outPath,"mae_",filename.prefix,".rda"),
                         met.platform = "450K",
                         genome = "hg38",
                         TCGA = TRUE)
      }
  }
  return(NA)
}
# Run Step: 1 (this is MAE file construction step)
#convertToMAE("tcgaProjects.csv",metFilesPath,expFilesPath,maeFilePath,localFilesPath,oldRunFilesPath)

# Num 2: Identifying differentially methylated probes 
diffMethAnalyzer<-function(cancerType,mae.file,maeFilePath,outFilePath,group_1,group_2,direction){
  setwd(outFilePath)
  sig.diff <- get.diff.meth(data = mae.file, 
                            group.col = "definition",
                            group1 =  group_1, # "Primary solid Tumor"
                            group2 = group_2, # "Solid Tissue Normal"
                            minSubgroupFrac = 0.2, # if supervised mode set to 1
                            sig.dif = 0.3,
                            diff.dir = direction, # "hypo", # Search for hypomethylated probes in group 1
                            cores = 1, 
                            dir.out = paste0(outFilePath,direction,"/",cancerType), 
                            pvalue = 0.01)
  
  head(sig.diff)  %>% datatable(options = list(scrollX = TRUE))
  # get.diff.meth automatically save output files. 
  # - getMethdiff.hypo.probes.csv contains statistics for all the probes.
  # - getMethdiff.hypo.probes.significant.csv contains only the significant probes which
  # is the same with sig.diff
  # - a volcano plot with the diff mean and significance levels
  dir(path = paste0(outFilePath,direction,"/",cancerType), pattern = "getMethdiff")  
  
  group1 <-  group_1
  group2 <- group_2
  out <- readr::read_csv(dir(path = paste0(outFilePath,direction,"/",cancerType), pattern = paste0("getMethdiff.",direction,".probes.csv"),full.names = TRUE))
  TCGAbiolinks:::TCGAVisualize_volcano(x = as.data.frame(out)[,grep("Minus",colnames(out),value = T)],
                                       y = out$adjust.p, 
                                       title =  paste0("Volcano plot - Probes ",
                                                       direction,"methylated in ", group1, " vs ", group2,"\n"),
                                       filename = NULL,
                                       label =  c("Not Significant",
                                                  paste0("Hypermethylated in ",group1),
                                                  paste0("Hypomethylated in ",group1)),
                                       ylab =  expression(paste(-Log[10],
                                                                " (FDR corrected P-values) [one tailed test]")),
                                       xlab =  expression(paste(
                                         "DNA Methylation difference (",beta,"-values)")
                                       ),
                                       x.cut = 0.3, 
                                       y.cut = 0.01)
  # job finished now back to old path
  setwd(maeFilePath)
  
}

createDirsFordiffMethAnalysis<-function(defaultPath,cancerType,direction,permuFolder="permu"){
  if (file.exists(paste0(defaultPath,direction,"/",cancerType))){
    print(paste0(cancerType," dir already exists..."))
  } else {
    dir.create(paste0(defaultPath,direction,"/",cancerType))
    print(paste0(cancerType," dir succesfully created..."))
  }
  if (file.exists(paste0(defaultPath,direction,"/",cancerType,"/",permuFolder))){
    print("dir already exists...")
  } else {
    dir.create(paste0(defaultPath,direction,"/",cancerType,"/",permuFolder))
    print("dir succesfully created...")
  }
  return(NA)
}

createDirsForDEGAnalysis<-function(defaultPath,cancerType){
  if (file.exists(paste0(defaultPath,cancerType,"/"))){
    print(paste0(cancerType," dir already exists..."))
    outPath<-paste0(defaultPath,cancerType,"/")
  } else {
    dir.create(paste0(defaultPath,cancerType,"/"))
    print(paste0(cancerType," dir succesfully created..."))
    outPath<-paste0(defaultPath,cancerType,"/")
  }
  return(outPath)
}

diffMethAnalyzer.caller<-function(maeFilePath,directionMeth,jobType=c("diff.meth","near.gene")){
  setwd(maeFilePath)
  maeFiles<-list.files(".")
  #maeFiles<-"mae_blca.rda" # this line is for debugging...
  for (maeFile in maeFiles){
    maeWord<-gsub(".*[_]([^.]+)[.].*", "\\1", maeFile) # https://stackoverflow.com/questions/23518325/how-to-extract-substring-between-patterns-and-in-r
    mae.file<-miceadds::load.Rdata2(filename=maeFile)
    print(maeWord)
    #
    if(length(unique(mae.file$definition))==1){
      checkGroups<-unique(mae.file$definition)
    }else{
      checkGroups<-names(sort(table(mae.file$definition),decreasing = T))
      #checkGroups <- names(checkGroups)
    }
    #
    if(length(checkGroups)>1){
      createDirsFordiffMethAnalysis(sigDiffResultsPath,maeWord,directionMeth,"permu")
      if(length(checkGroups)>2){ # if more than two we use only first 2 so we filter the mae file
        group.names<-checkGroups[1:2]
        print(group.names)
        met.rowNums<-which(mae.file@colData$definition %in% group.names)
        print(length(met.rowNums))
        exp.rowNums<-met.rowNums+nrow(mae.file@colData)
        sampleSet.rowNums<-c(met.rowNums,exp.rowNums)
        # now we have indices now is the time for filter...
        filteredMet<-mae.file@colData[met.rowNums,]
        filteredSamples<-mae.file@sampleMap[sampleSet.rowNums,]
        filteredDNA_methylation<-mae.file@ExperimentList$`DNA methylation`[,met.rowNums]
        filteredGene_expression<-mae.file@ExperimentList$`Gene expression`[,met.rowNums]
        # now attach filtered values to the existing mae file...
        mae.file@colData<-filteredMet
        print(nrow(mae.file@colData))
        mae.file@sampleMap<-filteredSamples
        mae.file@ExperimentList$`DNA methylation`<-filteredDNA_methylation
        mae.file@ExperimentList$`Gene expression`<-filteredGene_expression
        # now call the diff analysis function
        print("filtered one is running...")
        if(jobType=="diff.meth"){
          diffMethAnalyzer(maeWord,mae.file,maeFilePath,sigDiffResultsPath,group.names[1],group.names[2],directionMeth)
        }else if(jobType=="near.gene"){
          probeGenePairsFinder4CancerTypes(maeWord,mae.file,maeFile.Path,directionMeth,sigDiffResultsPath,group.names[1],group.names[2])
        }
      }else{
        # no filter justcall the diff analysis function
        print("not filtered one is running...")
        group.names<-unique(mae.file$definition)
        if(jobType=="diff.meth"){
          diffMethAnalyzer(maeWord,mae.file,maeFilePath,sigDiffResultsPath,group.names[1],group.names[2],directionMeth)
        }else if(jobType=="near.gene"){
          probeGenePairsFinder4CancerTypes(maeWord,mae.file,maeFile.Path,directionMeth,sigDiffResultsPath,group.names[1],group.names[2])
        }
      }
    }else{
      print(paste0(maeWord," has only one group so analyze for this type is skipped..."))
    }
  }# end for

  return(NA)
  
}
# Run Step: 2 (this is meth analysis step)
#diffMethAnalyzer.caller(maeFilePath,"hypo","diff.meth")

## ANALİZ SONUÇLARI İÇİN DEĞİŞİK PARAMETRE GİRİŞLERİ AŞAĞIDAN İTİBAREN YAPILACAK

# Num 3: Identifying putative probe-gene pairs
probeGenePairsFinder4CancerTypes<-function(cancerType,mae.file,maeFile.Path,directionMeth,outFilePath,group_1,group_2){
  setwd(outFilePath)
  sig.diff <- read.csv(paste0(sigDiffResultsPath,directionMeth,"/",cancerType,"/","getMethdiff.",directionMeth,".probes.significant.csv"))
  nearGenes <- GetNearGenes(data = mae.file, 
                            probes = sig.diff$probe, 
                            numFlankingGenes = 40) # 10 upstream and 10 dowstream genes (default is 20)

  HypoOrHyper.pair <- get.pair(data = mae.file,
                               group.col = "definition",
                               group1 =  group_1, # "Primary solid Tumor"
                               group2 = group_2, # "Solid Tissue Normal",
                               nearGenes = nearGenes,
                               mode = "unsupervised",
                               permu.dir = paste0(outFilePath,directionMeth,"/",cancerType,"/","permu"), #"result/permu"
                               permu.size = 1000, # Please set to 100000 to get significant results (default is 100)
                               raw.pvalue = 0.05, # 0.05
                               Pe = 0.05, # Please set to 0.001 to get significant results (default is 0.01)
                               filter.probes = TRUE, # See preAssociationProbeFiltering function
                               filter.percentage = 0.05,
                               filter.portion = 0.3,
                               dir.out = paste0(outFilePath,directionMeth,"/",cancerType),
                               cores = 1, # default is 1
                               label = directionMeth) # "hypo") 
  HypoOrHyper.pair %>% datatable(options = list(scrollX = TRUE))
  saveRDS(HypoOrHyper.pair, file = paste0(directionMeth,"/",cancerType,"/",cancerType,"_",directionMeth,"_meth.rds"))
  # get.pair automatically save output files. 
  # getPair.hypo.all.pairs.statistic.csv contains statistics for all the probe-gene pairs.
  # getPair.hypo.pairs.significant.csv contains only the significant probes which is 
  # same with Hypo.pair.
  dir(path = paste0(outFilePath,directionMeth,"/",cancerType), pattern = "getPair") 
  # job finished now back to old path
  setwd(maeFilePath)
    #
}
# Run Step: 3 (this is the near gene step)

# runmethods.neargene<-c("hypo","hyper")
# for (runType in runmethods.neargene){
#   diffMethAnalyzer.caller(maeFilePath,runType,"near.gene")
# }

# Num 4: Perform DEA (Differential expression analysis) to identify differentially expressed genes (DEGs)
DiffExpressionAnalyzer<-function(initPath,outPath){
  setwd(initPath)
  expFiles<-list.files(".")
  #expFiles<-"blcaExp.rda" # this line is for debugging...
  for (expFile in expFiles){
    expWord<-sub("Exp.rda", "", expFile)
    exp.file<-miceadds::load.Rdata2(filename=expFile)
    print(expWord)
    conditionTypes<-c() # just init for hold the data later
    if(length(unique(exp.file$shortLetterCode))==1){
      checkGroups<-unique(exp.file$shortLetterCode)
    }else{
      checkGroups<-names(sort(table(exp.file$shortLetterCode),decreasing = T))
    }
    if(length(checkGroups)>1){
      # create the save dir if not exists
      setwd(createDirsForDEGAnalysis(outPath,expWord))
      # Normalization common steps...
      dataPrep <- TCGAanalyze_Preprocessing(object = exp.file, 
                                            cor.cut = 0.6)
      dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                            geneInfo = geneInfoHT,
                                            method = "gcContent")
      dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                        method = "quantile", 
                                        qnt.cut =  0.25)
      # Normalization additional common steps..
      # UseRaw_afterFilter: Keep raw counts after filtering
      dataPrep_raw <- UseRaw_afterFilter(dataPrep, dataFilt)
      # save DEG
      saveRDS(dataPrep_raw, file = paste0(expWord,"_dataPrep_raw.rds")) #like acc_DEG.rds
      # TCGA_MolecularSubtype: Query subtypes for cancer data:
      data.subtypes<-TCGA_MolecularSubtype(colnames(dataFilt))
      saveRDS(data.subtypes, file = paste0(expWord,"_data.subtypes.rds")) #like acc_DEG.rds
      if(length(checkGroups)>2){ # if more than two we use only first 2 so we filter the mae file
        group.names<-checkGroups[1:2]
        print(group.names)
        exp.rowNums<-which(exp.file@colData$shortLetterCode %in% group.names)
        print(length(exp.rowNums))
        # now we have indices now is the time for filter...
        filteredExp<-exp.file@colData[exp.rowNums,]
        # now attach filtered values to the existing mae file...
        exp.file@colData<-filteredExp
        print(nrow(exp.file@colData))
        # now call the diff analysis function
        print("filtered one is running...")
        # Normalization seperate steps..
        samples1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = group.names[1])
        samples2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = group.names[2])
        # short job before run the DEGS we here make the condition types dynamic
        my.vector<-exp.file$shortLetterCode
        names(my.vector)<-exp.file$definition
        for (group.name in group.names){ # for exact match we use loop
          val<-unique(names(my.vector)[my.vector%in%group.name])
          conditionTypes<-c(conditionTypes,val)
        }
        #conditionTypes<-stringr::word(unique(names(my.vector)[my.vector%in%group.names]),-1)
        #conditionTypes<-unique(names(my.vector)[my.vector%in%group.names]) # this is before loop used code 
        print(paste0(group.names,"<---------->",conditionTypes))
        # run DEG
        dataDEGs <- TCGAanalyze_DEA(
          batch.factors = c("Plate"),
          mat1 = dataFilt[, samples1],
          mat2 = dataFilt[, samples2],
          Cond1type = conditionTypes[1],
          Cond2type = conditionTypes[2]
        )
        # save DEG
        saveRDS(dataDEGs, file = paste0(expWord,"_DEG.rds")) #like acc_DEG.rds
      }else{
        # no filter justcall the diff analysis function
        print("not filtered one is running...")
        group.names<-unique(exp.file$shortLetterCode)
        
        # Normalization seperate steps..
        samples1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = group.names[1])
        samples2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = group.names[2])
        # short job before run the DEGS we here make the condition types dynamic
        my.vector<-exp.file$shortLetterCode
        names(my.vector)<-exp.file$definition
        for (group.name in group.names){ # for exact match we use loop
          val<-unique(names(my.vector)[my.vector%in%group.name])
          conditionTypes<-c(conditionTypes,val)
        }
        #conditionTypes<-stringr::word(unique(names(my.vector)[my.vector%in%group.names]),-1)
        #conditionTypes<-unique(names(my.vector)[my.vector%in%group.names]) # this is before loop used code 
        print(paste0(group.names,"<---------->",conditionTypes))
        # run DEG
        dataDEGs <- TCGAanalyze_DEA(
          batch.factors = c("Plate"),
          mat1 = dataFilt[, samples1],
          mat2 = dataFilt[, samples2],
          Cond1type = conditionTypes[1],
          Cond2type = conditionTypes[2]
        )
        # save DEG
        saveRDS(dataDEGs, file = paste0(expWord,"_DEG.rds")) #like acc_DEG.rds
      }
    }else{
      print(paste0(expWord," has only one group so analyze for this type is skipped..."))
    }
    # return back to init path in order to pull next exp file...
    setwd(initPath)
    }
    #
  
}
# Run Step: 4 (DEG analysis step)
#DiffExpressionAnalyzer(expFilesPath,dEGresultsPath)

# Num 5: Calculate similarities based on immune escape genes
calculateGeneSimilaritiesBasedOnImmuneEscapeGenes4CancerTypes<-function(initPath,methType=c("hypo","hyper")){
  setwd(initPath) # this is sig fiff meth folder we want to start one folder up
  setwd("..") # main results path
  for(run.type in methType){
    all.file.names<-list.files(paste0("sig_diff_meth/",run.type,"/"))
    processed.file.names<-sub("\\/.*", "", list.files(paste0("sig_diff_meth/",run.type,"/"),pattern = "similarities",recursive = T))
    file.names<-setdiff(all.file.names, processed.file.names)
    file.names<-"ucec" # this line is for debugging...
    for(file.name in file.names){
      print(paste0(run.type,"<------>",file.name," job is running..."))
      # meth_genes covers both hypo and hyper
      DEGenes <- readRDS(paste0("degs/",file.name,"/",file.name,"_DEG.rds")) # blca_DEG.rds  DEGenes <- readRDS("results/luad/luad_DEG.rds")

      # Read from rds or csv doesn't matter contains the same data...
      #meth_genes<-readRDS(paste0("sig_diff_meth/",run.type,"/",file.name,"/",file.name,"_",run.type,"_meth.rds")) # blca_hypo_meth.rds)

      if(run.type=="hypo"){
        overExpressedGenes  <- DEGenes %>% rownames_to_column('gene') %>%  filter(FDR<0.05 & logFC >= 1.5)
        meth_genes<-read.csv(file=paste0("sig_diff_meth/",run.type,"/",file.name,"/getPair.",run.type,".pairs.significant.csv"))
        meth.over.exp.genes <- intersect(meth_genes$GeneID,overExpressedGenes$gene) 
        if(length(meth.over.exp.genes)>1){ #  Gene list should contain at least 2 elements!
          meth.over.exp.genes <- convertEnsembleID2GeneID(meth.over.exp.genes,run.type,file.name)
          saveRDS(meth.over.exp.genes, file = paste0("sig_diff_meth/",run.type,"/",file.name,"/",run.type,".meth.over.exp.genes.rds")) #same as converted_geneIDs.csv
          print(paste0('converted enseml to genes are:',meth.over.exp.genes))
          hypo.meth.similarities <- calculateGeneSimilaritiesBasedOnImmuneEscapeGenes(meth.over.exp.genes,immuneGenesFilePath)
          saveRDS(hypo.meth.similarities, file = paste0(initPath,run.type,"/",file.name,"/",run.type,".meth.similarities.rds")) #hypo.meth.similarities.rds
          cat(" ===> ",run.type," <---------> ",file.name," similarity file saved successfully...","\n")
        }else{
          cat(" ===> there is no or <2 ensemble gene to convert so iteration skipped for ",file.name," cancer type...","\n")
        }

      }else if(run.type=="hyper"){
        underExpressedGenes  <- DEGenes %>% rownames_to_column('gene') %>%  filter(FDR<0.05 & logFC <= -1.5)
        meth_genes<-read.csv(file=paste0("sig_diff_meth/",run.type,"/",file.name,"/getPair.",run.type,".pairs.significant.csv"))
        meth.under.exp.genes  <- intersect(meth_genes$GeneID,underExpressedGenes$gene)
        if(length(meth.under.exp.genes)>1){ #  Gene list should contain at least 2 elements!
          meth.under.exp.genes <- convertEnsembleID2GeneID(meth.under.exp.genes,run.type,file.name)
          saveRDS(meth.under.exp.genes, file = paste0("sig_diff_meth/",run.type,"/",file.name,"/",run.type,".meth.under.exp.genes.rds")) #hyper.meth.over.exp.genes.rds
          hyper.meth.similarities <- calculateGeneSimilaritiesBasedOnImmuneEscapeGenes(meth.under.exp.genes,immuneGenesFilePath)
          saveRDS(hyper.meth.similarities, file = paste0(initPath,run.type,"/",file.name,"/",run.type,".meth.similarities.rds")) #hyper.meth.similarities.rds
          cat(" ===> ",run.type," <---------> ",file.name," similarity file saved successfully...","\n")
        }else{
          cat(" ===> there is no or <2 ensemble gene to convert so iteration skipped for ",file.name," cancer type...","\n")
        }

      }
    }
  }
return(NA)
  
}
#calculateGeneSimilaritiesBasedOnImmuneEscapeGenes4CancerTypes(sigDiffResultsPath)

# Num 6: Clustering genes
clusterGenes4TumorType <- function(initPath,outPath,methType=c("hypo","hyper")){
  setwd(initPath)
  for(run.type in methType){
    #file.names<-list.files(paste0(initPath,run.type))
    file.names<-sub("\\/.*", "", list.files(paste0(run.type,"/"),pattern = "similarities",recursive = T))
    #file.names<-"blca"
    for(file.name in file.names){
      hypo_hyper.meth.similarities <- readRDS(paste0(run.type,"/",file.name,"/",run.type,".meth.similarities.rds"))
      # check that similartiy results are valid
      if(sum(dim(hypo_hyper.meth.similarities))>0){
        clustersHypo_Hyper <- clusterGenes(outPath,file.name,run.type,hypo_hyper.meth.similarities)
        if(!length(clustersHypo_Hyper)==F){ # to avoid logical(0) !length used <-> ol code is: !is.na(clustersHypo_Hyper)
          saveRDS(clustersHypo_Hyper, file = paste0(initPath,run.type,"/",file.name,"/clusters",firstup(run.type),".rds")) # like: clustersHyper.rds
        }
      }else{
        cat(" ===> no similarity found for", run.type," <---------> ",file.name," so this job analysis skipped...","\n")
      }

    } # enf for file.names
  } # end for run.types

}
#clusterGenes4TumorType(sigDiffResultsPath,clusterResultsPath)

# Num 7: Pathfinder Analysis on Enrichment Genes

pathfindeR_analysis <- function(outPath,cancerType,runType,clusters){ # output:pathFinderResultsPath
  geneID.fromCluster<-unlist(clusters)
  DEG_with_Entrez_IDs<-getEnsembleIDfromGeneIDThenFilterForPathFinder(geneID.fromCluster,cancerType)
  print(DEG_with_Entrez_IDs)
  # main part is here...
  for (analyzeType in c("Reactome","KEGG", "BioCarta","GO-BP", "GO-CC", "GO-MF")){
    clusterResults  <- c()
    for (i in 1:length(clusters)){
      success = FALSE
      # cluster <- DEG <- subset(DEG <- with <- Entrez <- IDs,gene %in% unlist(clustersHypo[1]))
      #------------------------
      cluster_DEG <- subset(DEG_with_Entrez_IDs,GeneID %in% unlist(clusters[i]))
      cluster_DEG$geneSymbol <- getSYMBOL(c(unlist(cluster_DEG$GeneID)), data='org.Hs.eg')
      sigGenes  <- data.frame(as.character(cluster_DEG$geneSymbol),
                              cluster_DEG$logFC,
                              cluster_DEG$FDR)
      #------------------------

      colnames(sigGenes) <- c("Gene.symbol","logFC","adj.P.Val")
      sigGenes  <- sigGenes[complete.cases(sigGenes), ]
      for(pinDB in c("GeneMania", "IntAct", "KEGG","Biogrid")){
        out = paste0(outPath,cancerType,"/",runType,"/",analyzeType,"/",pinDB)
        tryCatch({
          print(paste0("Analyze starting for:",analyzeType,"on",pinDB))
          tmpOutput <- run_pathfindR.local(sigGenes, gene_sets =analyzeType,pin_name_path = pinDB,
                                     output_dir=out,silent_option = TRUE, adj_method="BH",p_val_threshold=0.1,
                                     n_processes = 10)
        }, error = function(e) {
          print(paste('error:', e))
        })
      }
    }
  }
  # end of the main part
  
}
#"hypo",
enrichWithPathFinder4TumorType <- function(initPath,outPath,conditionType=c("hypo")){# input: rootResultsPath output: pathFinderResultsPath
  setwd(initPath)
  for(run.type in conditionType){
    file.names<-sub("\\/.*", "", list.files(paste0("sig_diff_meth/",run.type,"/"),pattern = paste0("clusters",firstup(run.type),".rds"),recursive = T))
    file.names<-"ucec"
    if(length(file.names)>0){ # run if is there a file to run for...
      for(file.name in file.names){
        cat(" ===> now trying to run for ", run.type," <---------> ",file.name," pathfinding jobs...","\n")
        clustersHypoOrHyper <- readRDS(paste0("sig_diff_meth/",run.type,"/",file.name,"/clusters",firstup(run.type),".rds"))
        pathFindeR_results <- pathfindeR_analysis(outPath,file.name,run.type,clustersHypoOrHyper)
      }
    }else{
      cat(" ===> none cancer type found for ", run.type," pathfinding jobs so whole process skipped...","\n")
    }
  }# end for run type
}

enrichWithPathFinder4TumorType(rootResultsPath,pathFinderResultsPath)

