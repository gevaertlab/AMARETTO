#' AMARETTO_Download
#'
#' Downloading TCGA dataset for AMARETTO analysis
#' @param CancerSite TCGA cancer code for data download
#' @param TargetDirectory Directory path to download data
#' @param downloadData
#'
#' @return
#' @import RCurl
#' @import curatedTCGAData
#' @import limma
#' @importFrom RCurl getURL
#' @importFrom limma strsplit2
#' @export
#'
#' @examples
#' TargetDirectory <- "./Downloads/" # path to data download directory
#' CancerSite <- "LIHC"
#' DataSetDirectories <- AMARETTO_Download(CancerSite,TargetDirectory)
AMARETTO_Download <- function(CancerSite,TargetDirectory,downloadData=TRUE) {
  cat('Downloading MA, CNV and MET data for:',CancerSite,'\n')
  Cancers=c('BLCA','BRCA','LUAD','LUSC','COADREAD','HNSC','KIRC','GBM','OV','LAML','UCEC','COAD','READ')
  if (!(CancerSite %in% Cancers)) {
    cat('This TCGA cancer site/type was not tested, continue at your own risk.\n')
  }
  dir.create(TargetDirectory,showWarnings=FALSE)
  TCGA_acronym_uppercase=toupper(CancerSite)
  assays <- c("RNASeq2GeneNorm")
  MAEO <- curatedTCGAData::curatedTCGAData(CancerSite, assays, FALSE)
  saveRDS(MAEO, file=paste0(TargetDirectory, CancerSite, "_RNASeq_MAEO.rds"))
  dataType='analyses'
  dataFileTag='CopyNumber_Gistic2.Level_4'
  cat('Searching CNV data for:',CancerSite,'\n')
  CNVdirectory=get_firehoseData(downloadData,saveDir=TargetDirectory,TCGA_acronym_uppercase=TCGA_acronym_uppercase,dataType=dataType,dataFileTag=dataFileTag)
  setwd("..")
  return(list(MAdirectory=TargetDirectory,CNVdirectory=CNVdirectory))
}

#' get_firehoseData
#'
#' Downloading TCGA dataset via firehose
#' @param downloadData
#' @param saveDir
#' @param TCGA_acronym_uppercase
#' @param dataType
#' @param dataFileTag
#' @param FFPE
#' @param fileType
#' @param gdacURL
#' @param untarUngzip
#' @param printDisease_abbr
#'
#' @return
#' @keywords internal
#' @examples
#' data_download <- get_firehoseData(downloadData=TRUE,saveDir = "./",TCGA_acronym_uppercase = "LUAD",dataType="stddata",dataFileTag = "mRNAseq_Preprocess.Level_3", FFPE=FALSE,fileType= "tar.gz",  gdacURL= "http://gdac.broadinstitute.org/runs/",untarUngzip=TRUE,printDisease_abbr=FALSE)

get_firehoseData <- function(downloadData=TRUE,saveDir = "./",TCGA_acronym_uppercase = "LUAD",dataType="stddata",dataFileTag = "mRNAseq_Preprocess.Level_3",
                             FFPE=FALSE,fileType= "tar.gz",  gdacURL= "http://gdac.broadinstitute.org/runs/",untarUngzip=TRUE,printDisease_abbr=FALSE){

  # Cases Shipped by BCR  # Cases with Data*  Date Last Updated (mm/dd/yy)
  cancers <- c("Acute Myeloid Leukemia [LAML] \n","Adrenocortical carcinoma [ACC]	\n",
               "Bladder Urothelial Carcinoma [BLCA] \n",	"Brain Lower Grade Glioma [LGG] \n",
               "Breast invasive carcinoma [BRCA] \n","Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] \n",
               "Cholangiocarcinoma [CHOL] \n",	"Colon adenocarcinoma [COAD] \n",	"Esophageal carcinoma [ESCA] \n",
               "Glioblastoma multiforme [GBM] \n",	"Head and Neck squamous cell carcinoma [HNSC]	\n",
               "Kidney Chromophobe [KICH]	\n","Kidney renal clear cell carcinoma [KIRC]	\n",
               "Kidney renal papillary cell carcinoma [KIRP]	\n","Liver hepatocellular carcinoma [LIHC]	\n",
               "Lung adenocarcinoma [LUAD]	\n", "Lung squamous cell carcinoma [LUSC] \n",
               "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma [DLBC]	\n","Mesothelioma [MESO] \n",
               "Ovarian serous cystadenocarcinoma [OV]	\n","Pancreatic adenocarcinoma [PAAD]	\n",
               "Pheochromocytoma and Paraganglioma [PCPG] \n","Prostate adenocarcinoma [PRAD] \n",
               "Rectum adenocarcinoma [READ]	\n","Sarcoma [SARC]	\n","Skin Cutaneous Melanoma [SKCM]	\n",
               "Stomach adenocarcinoma [STAD] \n","Testicular Germ Cell Tumors [TGCT] \n","Thymoma [THYM] \n",
               "Thyroid carcinoma [THCA]	\n","Uterine Carcinosarcoma [UCS]	 \n",
               "Uterine Corpus Endometrial Carcinoma [UCEC]	\n","Uveal Melanoma [UVM] \n");

  cancers_acronyms <-c("LAML","ACC","BLCA","LGG","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","LIHC","LUAD",
                       "LUSC","DLBC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THYM","THCA","UCS","UCEC","UVM")

  if(printDisease_abbr){
    return(cat("here are the possible TCGA database disease acronyms. \nRe-run this function with printDisease_abbr=FALSE to then run an actual query.\n\n",cancers));
  }
  if (TCGA_acronym_uppercase %in% cancers_acronyms){
    gdacURL_orig <- gdacURL
    urlData <- getURL(gdacURL)
    urlData <- strsplit2(urlData,paste(dataType,"__",sep=""))
    urlData <- urlData[,2:dim(urlData)[2]]
    urlData <- strsplit2(urlData,"/")
    urlData <- urlData[,1]
    urlData <- as.POSIXct(strptime(urlData, "%Y_%m_%d"))
    dateData <- as.Date(as.character(urlData[which(!is.na(urlData))]))
    lastDate <- dateData[match( summary(dateData)[which(names(summary(dateData))=="Max.")], dateData)]
    lastDate <- gsub("-","_",as.character(lastDate))
    lastDateCompress <- gsub("_","",lastDate)
    gdacURL <- paste(gdacURL,dataType,"__",lastDate,"/data/",TCGA_acronym_uppercase,"/",lastDateCompress,"/",sep="")
    urlData <- getURL(gdacURL)
    urlData <- strsplit2(urlData,"href=\\\"")
    while(length(grep("was not found",urlData))>0) {
      cat(file.path("\tNOTE: the TCGA run dated ",lastDate, "for ", dataType," for disease ",TCGA_acronym_uppercase," isn't available 	for download yet.\n"))
      cat("\tTaking the run dated just before this one.\n")
      dateData <-  dateData[-which(dateData==(summary(dateData)[which(names(summary(dateData))=="Max.")]))]
      lastDate <- dateData[match( summary(dateData)[which(names(summary(dateData))=="Max.")], dateData)]
      lastDate <- gsub("-","_",as.character(lastDate))
      lastDateCompress <- gsub("_","",lastDate)
      gdacURL <- paste(gdacURL_orig,dataType,"__",lastDate,"/data/",TCGA_acronym_uppercase,"/",lastDateCompress,"/",sep="")
      urlData <- getURL(gdacURL)
      urlData <- strsplit2(urlData,"href=\\\"")
      if(length(dateData)<=1){
        break
      }
    }
    if(length(grep("was not found",urlData))>0){
      #this disease may not even be in the analyses directory yet.
      stop( paste0("\nNot returning any viable url data paths after searching by date for disease ",TCGA_acronym_uppercase,". No data was downloaded.\n"))
    }
    if (FFPE) {
      urlData <- urlData[grep("FFPE",urlData)]
      if(length(urlData)==0){
        stop("\nNo FFPE data found for this query. Try FFPE=FALSE.\n")
      }
    } else {
      if(length(grep("FFPE",urlData))>0){
        urlData <- urlData[-grep("FFPE",urlData)]
      }
      if(length(urlData)==0){
        stop("\nNo non-FFPE data found for this query. Try FFPE=TRUE.\n")
      }
    }
    fileName <- urlData[grep(dataFileTag,urlData)]
    if(length(fileName)==0){
      warnMessage <- paste0("\nNot returning any viable url data paths after searching by date for disease ",TCGA_acronym_uppercase," 	for data type ",dataFileTag ,".No data was downloaded.\n")
      warning(warnMessage)
      return(NA)
    }
    fileName <- strsplit2(fileName,"tar.gz")[1,1]
    fileName <- paste(fileName,fileType,sep="")
    gdacURL <- paste(gdacURL,fileName,sep="")
    print(fileName)
    saveDir <- paste(saveDir,"gdac_",lastDateCompress,'/',sep="")
    if(downloadData){
      cat("\tDownloading",dataFileTag,"data, version:",lastDate,"\n")
      cat("\tThis may take 10-60 minutes depending on the size of the data set.\n")
      dir.create(saveDir,showWarnings=FALSE)
      setwd(saveDir)
      download.file(gdacURL,fileName,quiet=TRUE,mode="wb")
      if(fileType=="tar.gz" && untarUngzip) {
        cat("\tUnpacking data.\n")
        tarfile=fileName
        untar(tarfile)
        fileToRemove <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
        file.remove(paste0(fileToRemove))
      } else if(untarUngzip) {
        warning("File expansion/opening only built in for tar.gz files at the moment.\n")
      }
      finalDir <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
      finalDir <- strsplit2(finalDir,fileType)
      finalDir <- substr(finalDir,start=0,stop=(nchar(finalDir)-1))
      finalDir <- paste0(saveDir,finalDir)
      cat("\tFinished downloading",dataFileTag,"data to",finalDir,"\n")
    } else {
      cat("download data url is :\n ",gdacURL,'\n')
      finalDir <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
      finalDir <- strsplit2(finalDir,fileType)
      finalDir <- substr(finalDir,start=0,stop=(nchar(finalDir)-1))
      finalDir <- paste0(saveDir,finalDir)
    }
    DownloadedFile=paste0(finalDir,'/')
    return(DownloadedFile)
  } else{
    return(cat(paste0("No data correspond to cancer ",TCGA_acronym_uppercase,"\n")))
  }
}


#' AMARETTO_DownloadResults
#' Retrieve a download of all the data linked with the run (including heatmaps)
#' @param AMARETTOinit AMARETTO initialize output
#' @param AMARETTOresults AMARETTO results output
#' @param data_address Directory to save data folder
#' @param Heatmaps Output heatmaps as pdf
#' @import doParallel
#' @import tidyverse
#' @return
#' @export
#'
#' @examples
#' AMARETTO_DownloadResults(AMARETTOinit,AMARETTOresults,"./")

AMARETTO_ExportResults <-function(AMARETTOinit,AMARETTOresults,data_address,Heatmaps=TRUE){

  if (!dir.exists(data_address)){
    stop("Output directory is not existing.")
  }

  #add a date stamp to the output directory
  output_dir<-paste0("AMARETTOresults_",gsub("-|:","",gsub(" ","_",Sys.time())))
  dir.create(file.path(data_address,output_dir))
  
  #parallelize the heatmap production
  cluster <- makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
  registerDoParallel(cluster,cores=NrCores)

  if(Heatmaps==TRUE){
    foreach (ModuleNr = 1:NrModules, .packages = c('AMARETTO')) %dopar% {
      pdf(file=paste(data_address,output_dir,"/Module_",as.character(ModuleNr),".pdf",sep=""))
      AMARETTO_VisualizeModule(AMARETTOinit, AMARETTOresults=AMARETTOresults, CNV_matrix, MET_matrix, ModuleNr=ModuleNr)
      dev.off()
    }
  }

  stopCluster(cluster)
  #save rdata files for AMARETTO_Run and AMARETTO_Initialize output
  save(AMARETTOresults, file=paste0(data_address,output_dir,"/amarettoResults.RData"))
  save(AMARETTOinit, file=paste0(data_address,output_dir,"/amarettoInit.RData"))
  
  #save some tables that might be useful for further analysis
  write_gct(AMARETTOresults$ModuleData,paste0(data_address,output_dir,'/ModuleData_amaretto.gct'))
  write_gct(AMARETTOresults$ModuleMembership,paste0(data_address,output_dir,'/ModuleMembership_amaretto.gct'))
  write_gct(AMARETTOresults$RegulatoryProgramData,paste0(data_address,output_dir,'/RegulatoryProgramData_amaretto.gct'))
  write_gct(AMARETTOresults$RegulatoryPrograms,paste0(data_address,output_dir,'/RegulatoryPrograms_amaretto.gct'))
  readr::write_tsv(as.data.frame(AMARETTOresults$AllGenes),paste0(data_address,output_dir,'/AllGenes_amaretto.tsv'))
  readr::write_tsv(as.data.frame(AMARETTOresults$AllRegulators),paste0(data_address,output_dir,'/AllRegulators_amaretto.tsv'))
  readr::write_tsv(as.data.frame(AMARETTOresults$NrModules),paste0(data_address,output_dir,'/NrModules_amaretto.tsv'))

  #zip the file
  zip(zipfile = paste0(data_address,output_dir),files=paste0(data_address,output_dir))
}


#' write_gct
#' Write a gct file with the proper header.
#' @param data_in
#' @param file_address
#'
#' @return
#' @keywords internal
#' @examples
#' write_gct(AMARETTOresults$RegulatoryPrograms,'RegulatoryPrograms.gct')
write_gct<-function(data_in,file_address){
  header_gct<-paste0('#1.2\n',nrow(data_in),'\t',ncol(data_in))
  data_in<-rownames_to_column(as.data.frame(data_in),"Name") %>% dplyr::mutate(Description=Name) %>% dplyr::select(Name,Description,everything())
  write(header_gct,file=file_address,append = FALSE)
  readr::write_tsv(data_in,file_address,append = TRUE,col_names = TRUE)
}
