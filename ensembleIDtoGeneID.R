############# connection hepler functions #############################
connect2Net <- function(){
  print("now opening the connection...")
  chromeverVec <- unlist(binman::list_versions("chromedriver"))
  geckodriverVec <- unlist(binman::list_versions("geckodriver"))
  rD <- RSelenium::rsDriver(browser=c("firefox"), chromever=chromeverVec[1], geckover="latest", port = 4446L) # rD <- rsDriver()
  remDr <- rD[["client"]]                     #geckodriverVec[3]
  return(list(connVar1 = remDr, connVar2 = rD))
}
closeTheConnection <- function(remoteDr,rDrv){
  print("now closing the connection...")
  remoteDr$close() # for avoid server signals port is already in use
  rDrv$server$stop()
  rm(remoteDr, rDrv) # for avoid server signals port is already in use
  gc() # for avoid server signals port is already in use
  return(NA)
}
############# end of connection helper functions #######################

############# job helper functions #####################################
jobFileDiagnostics.write <- function(file.name,jobOutDF){ # works under main results path : "/media/burak/Seagate Backup Plus Drive/meth_pan_cancer/outputs/results/
  # job 1. check if file exists in the disk
  if(file.exists(file.name)){
    print("there it is...")
    write.table(jobOutDF,file=paste0(file.name), append = T, sep=',', row.names=F, col.names=F) 
  }else{
    print("no job file found...")
    write.csv(jobOutDF,paste0(file.name),row.names = F)
  }
}

jobFileDiagnostics.read <- function(file.name){ # works under main results path : "/media/burak/Seagate Backup Plus Drive/meth_pan_cancer/outputs/results/
  # job 1. check if file exists in the disk
  if(file.exists(file.name)){
    print("there it is...")
    tempDF<-read.csv(paste0(file.name), header = T, stringsAsFactors = F)
    return(tempDF)
  }else{
    print("no job file found...")
    return(NA)
  }
}

jobFileDiagnostics.matchSame <- function(ens.ids,file.name){
  tempDF<-jobFileDiagnostics.read(file.name)
  if(!is.na(tempDF)){
    if(length(which(tempDF$EnsembleID %in% ens.ids))>0){
      matchedEnsembls <- with(tempDF, EnsembleID[which(EnsembleID %in% ens.ids)]) # pull matching ensembl ones
      matchedGeneIDs <- with(tempDF, GeneID[which(EnsembleID %in% ens.ids)]) # pull matching gene ones
      return(list(matchedVar1 = matchedEnsembls, matchedVar2 = matchedGeneIDs))
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }

}
############# end of job helper functions ##############################

convertEnsembleID2GeneID<-function(ens.ids,run.type,cancer.type,file.name="converted_geneIDs.csv"){ # setwd(sigDiffResultsPath)
  # run path arrangements
  tmp.path<-getwd()
  setwd(sigDiffResultsPath)
  # variable declarations
  Gene.ids<-c()
  Ens.idholder<-c()
  Pull.type<-c()
  # job0. Rselenium connection arrengement steps
  connectionVals<-connect2Net()
  remDr<-connectionVals$connVar1
  rD<-connectionVals$connVar2
  # job1. check for matching existing ids in order to save time
  matchedVals<-jobFileDiagnostics.matchSame(ens.ids,file.name)
  if(!is.na(matchedVals)){
    `%notin%` <- Negate(`%in%`) # r does not have default 'notin' operator so we just make one
    ens.ids<-ens.ids[ens.ids %notin% matchedVals$matchedVar1]
    #print(paste0('new matched ids:',ens.ids))
    Ens.idholder <- matchedVals$matchedVar1
    Gene.ids <- matchedVals$matchedVar2
    Pull.type <- rep("Synonym", each=length(Ens.idholder))
  }
  for (ens.id in ens.ids){
    pageTitleInfo <- "dummyControlWord"
    appURL <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=",ens.id) # like: https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000109501
    remDr$navigate(appURL)
    # control for page to be loaded...
    while(pageTitleInfo=="dummyControlWord"){ # we check that if list empty so the web elem doesn't show up yet
      pageTitleInfo <- tryCatch({unlist(remDr$getTitle())},
                                error = function(e){"dummyControlWord"})}
    # control that given ens id has really a gene id...
    checkElem.script <- "return document.getElementsByClassName('geneid').length > 0"
    checkElem.result <- unlist(remDr$executeScript(checkElem.script))  # returns empty list first like: list()
    if(checkElem.result){
      # get the gene ID
      checkElem.script <- "return document.querySelector('span.geneid').innerText"
      checkElem.result <- remDr$executeScript(checkElem.script)  # returns empty list first like: list()
      # convert pulled data to the proper format 
      checkElem.result<-unlist(checkElem.result) # Gene ID: 7466, updated on 13-Dec-2020
      gene.id<-trimws(stringr::str_extract(checkElem.result, ".\\d+(?=,)")) # now we have just 7466
      # update the value holders
      Gene.ids<-c(Gene.ids,gene.id)
      Ens.idholder<-c(Ens.idholder,ens.id)
      Pull.type<-c(Pull.type,"Original")
      cat("# ",which(ens.ids==ens.id),"/",length(ens.ids)," ===> ENSEMBL ID: ", ens.id," successfully converted to the <---------> GENE ID: ",gene.id,"\n")
    }else{
      cat("# ",which(ens.ids==ens.id),"/",length(ens.ids)," ===> ENSEMBL ID: ", ens.id," has no gene information!","\n")
    }
  }
  # # arrange the results to df
  tempDF<-data.frame(Gene.ids,Ens.idholder,Pull.type,cancer.type,run.type)
  colnames(tempDF)<-c("GeneID","EnsembleID","PullType","CancerType","RunType")
  # write.csv(tempDF,paste0("sig_diff_meth/",run.type,"/",cancer.type,"/converted_geneIDs.csv"),row.names = F)
  jobFileDiagnostics.write(file.name,tempDF)

  # Close the connection
  closeTheConnection(remoteDr=remDr,rDrv=rD)
  # Return back to default working directory
  setwd(tmp.path)
  return(as.list(Gene.ids)) # character type gives error so we convert it to list
}

convertEnsembleID2GeneID4DEGs<-function(DEG_file,run.type,cancer.type){ # setwd(sigDiffResultsPath)
  # job 0. convert DEG file to proper format
  DEGenes  <- DEG_file %>% rownames_to_column('gene')
  ens.ids <- DEGenes$gene[1:10]
  # job 1. init the necessary variables
  Gene.ids<-c()
  Ens.ids<-c()
  # job 2. Rselenium connection arrengement steps
  connectionVals<-connect2Net()
  remDr<-connectionVals$connVar1
  rD<-connectionVals$connVar2
  for (ens.id in ens.ids){
    pageTitleInfo <- "dummyControlWord"
    appURL <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=",ens.id) # like: https://www.ncbi.nlm.nih.gov/gene/?term=ENSG00000109501
    remDr$navigate(appURL)
    # control for page to be loaded...
    while(pageTitleInfo=="dummyControlWord"){ # we check that if list empty so the web elem doesn't show up yet
      pageTitleInfo <- tryCatch({unlist(remDr$getTitle())},
                                error = function(e){"dummyControlWord"})}
    # control that given ens id has really a gene id...
    checkElem.script <- "return document.getElementsByClassName('geneid').length > 0"
    checkElem.result <- unlist(remDr$executeScript(checkElem.script))  # returns empty list first like: list()
    if(checkElem.result){
      # get the gene ID
      checkElem.script <- "return document.querySelector('span.geneid').innerText"
      checkElem.result <- remDr$executeScript(checkElem.script)  # returns empty list first like: list()
      # convert pulled data to the proper format 
      checkElem.result<-unlist(checkElem.result) # Gene ID: 7466, updated on 13-Dec-2020
      gene.id<-trimws(stringr::str_extract(checkElem.result, ".\\d+(?=,)")) # now we have just 7466
      # update the value holders
      Gene.ids<-c(Gene.ids,gene.id)
      Ens.ids<-c(Ens.ids,ens.id)
      cat("# ",which(ens.ids==ens.id),"/",length(ens.ids)," ===> ENSEMBL ID: ", ens.id," successfully converted to the <---------> GENE ID: ",gene.id,"\n")
    }else{
      cat("# ",which(ens.ids==ens.id),"/",length(ens.ids)," ===> ENSEMBL ID: ", ens.id," has no gene information!","\n")
    }
  }
  # # arrange the results to df
  tempDF<-data.frame(Gene.ids,Ens.ids)
  colnames(tempDF)<-c("GeneID","EnsembleID")
  write.csv(tempDF,paste0("sig_diff_meth/",run.type,"/",cancer.type,"/from_DEG_converted_geneIDs.csv"),row.names = F)
  
  # Close the connection
  closeTheConnection(remoteDr=remDr,rDrv=rD)
  return(tempDF) # character type gives error so we convert it to list
}

getEnsembleIDfromGeneIDThenFilterForPathFinder<-function(gene.ids,cancerType){
  # job 1. init the necessary variables
  Ens.ids<-c()
  # job 2. Rselenium connection arrengement steps
  connectionVals<-connect2Net()
  remDr<-connectionVals$connVar1
  rD<-connectionVals$connVar2
  # job 3. convert Gene IDs to ensembl IDs
  for (gene.id in gene.ids){
    pageTitleInfo <- "dummyControlWord"
    appURL <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=",gene.id) # like: https://www.ncbi.nlm.nih.gov/gene/?term=7466
    remDr$navigate(appURL)
    # control for page to be loaded...
    while(pageTitleInfo=="dummyControlWord"){ # we check that if list empty so the web elem doesn't show up yet
      pageTitleInfo <- tryCatch({unlist(remDr$getTitle())},
                                error = function(e){"dummyControlWord"})}
    # control that given ens id has really a gene id...
    webElem<-remDr$findElement(using = "partial link text", value = "Ensembl")
    webElem.result<-unlist(webElem$getElementAttribute("text")) # result is like: "Ensembl:ENSG00000101605"
    ensembl.id<-sub('.*:', '', webElem.result) # we pull ENSG00000101605
    Ens.ids<-c(Ens.ids,ensembl.id)

  }
  # job 4.search matching ids in the corresponding DEG file by filtering
  DEG_file <- readRDS(paste0("degs/",cancerType,"/",cancerType,"_DEG.rds"))
  beforeNum<-nrow(DEG_file)
  DEGenes  <- DEG_file %>% rownames_to_column('EnsemblID') %>% filter(EnsemblID %in% Ens.ids) %>% mutate(GeneID = gene.ids)
  afterNum<-nrow(DEGenes)
  cat(" =======> ",beforeNum," Ensembl ids succesfully converted anf filtered to ",afterNum,"\n")

  # Close the connection
  closeTheConnection(remoteDr=remDr,rDrv=rD)
  return(DEGenes) # character type gives error so we convert it to list
}
