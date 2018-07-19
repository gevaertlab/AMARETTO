
##############################################################################################################################
# dataType can be stddata or analyses
# dataFileTag options for stddata are: mRNAseq_Preprocess.Level_3, 
# dataFileTag options for analyses are: CopyNumber_Gistic2.Level_4, Pathway_Paradigm_RNASeq_And_Copy_Number.Level_4, 
#                   Pathway_Paradigm_mRNA_And_Copy_Number.Level_4 and MutSigNozzleReport2.0.Level_4
# TCGA_acronym_uppercase: BRCA, LUAD, COADREAD, etc. 
##############################################################################################################################

##############################################################################################################################
#code for testing this script
##############################################################################################################################
# gdacURL= "http://gdac.broadinstitute.org/runs/"
# TCGA_acronym_uppercase="COADREAD"
# fileType= "tar.gz"
# downloadData=TRUE
# dataType="stddata"
# saveDir = "./"
# tmpData=get_firehoseData(downloadData,saveDir,TCGA_acronym_uppercase)
###############################################################################################################################

Download_CancerSite <- function(CancerSite,TargetDirectory,downloadData=TRUE) {    
	cat('Downloading MA, CNV and MET data for:',CancerSite,'\n')
	Cancers=c('BLCA','BRCA','LUAD','LUSC','COADREAD','HNSC','KIRC','GBM','OV','LAML','UCEC','COAD','READ')
	if (!(CancerSite %in% Cancers)) {
	    cat('This TCGA cancer site/type was not tested, continue at your own risk.\n')
	}
	
	# Creating the top level directory where all data will be stored.
	#command <- paste0("mkdir -p ",TargetDirectory)
	#system(command)
	dir.create(TargetDirectory,showWarnings=FALSE)
	
	# Settings
	TCGA_acronym_uppercase=toupper(CancerSite)
	
	# get RNA seq data (GBM does not have much RNAseq data.)
	dataType='stddata'	
	dataFileTag='mRNAseq_Preprocess.Level_3'	 

	#special case for GBM and OV, not enough RNAseq data, so using the microarray data instead
	if (CancerSite=="GBM") { 	             
		dataFileTag=c('Merge_transcriptome__agilentg4502a_07_1__unc_edu__Level_3__unc_lowess_normalization_gene_level__data','Merge_transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data')        	         
	} else if(CancerSite=="OV") {	               
		dataFileTag='Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data'        
	}	
	cat('Searching MA data for:',CancerSite,"\n")
	if (length(dataFileTag)==1) {	  
			MAdirectory=get_firehoseData(downloadData,saveDir=TargetDirectory,TCGA_acronym_uppercase=TCGA_acronym_uppercase,dataFileTag=dataFileTag)    	
		} else {	    
			MAdirectory=c()	  
		for (i in 1:length(dataFileTag)) {
			MAdirectory=c(MAdirectory,get_firehoseData(downloadData,saveDir=TargetDirectory,TCGA_acronym_uppercase=TCGA_acronym_uppercase,dataFileTag=dataFileTag[i]))	 
		}        
	}
	
	# get CNV GISTIC data.
	dataType='analyses'
	dataFileTag='CopyNumber_Gistic2.Level_4'
	cat('Searching CNV data for:',CancerSite,'\n')
	CNVdirectory=get_firehoseData(downloadData,saveDir=TargetDirectory,TCGA_acronym_uppercase=TCGA_acronym_uppercase,dataType=dataType,dataFileTag=dataFileTag)	
	return(list(MAdirectory=MAdirectory,CNVdirectory=CNVdirectory))
}



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
	
		#first column is junk
		urlData <- urlData[,2:dim(urlData)[2]]
		urlData <- strsplit2(urlData,"/")
		urlData <- urlData[,1]
		#see the strptime codes here: http://stat.ethz.ch/R-manual/R-devel/library/base/html/strptime.html
		#I am transferring this text to dates: that way, we can programmatically find the latest one.
		#the POSIXct class in R is the date-time class
		urlData <- as.POSIXct(strptime(urlData, "%Y_%m_%d"))
		#want to remove time zone: do as.Date
		dateData <- as.Date(as.character(urlData[which(!is.na(urlData))]))
		lastDate <- dateData[match( summary(dateData)[which(names(summary(dateData))=="Max.")], dateData)]
	
		#sub back _ symbols
		lastDate <- gsub("-","_",as.character(lastDate))
	
		lastDateCompress <- gsub("_","",lastDate)
		#need last "/" for it to find this page.
	
		gdacURL <- paste(gdacURL,dataType,"__",lastDate,"/data/",TCGA_acronym_uppercase,"/",lastDateCompress,"/",sep="")
	
		#now get full dataset name. we just have the http link to the last page we select this dataset from right now.
		urlData <- getURL(gdacURL)
	
	
		#regular expressions: need \ to have R recognize any " or \ that's actually in our text
		urlData <- strsplit2(urlData,"href=\\\"")
	
		while(length(grep("was not found",urlData))>0) {
			cat(file.path("\tNOTE: the TCGA run dated ",lastDate, "for ", dataType," for disease ",TCGA_acronym_uppercase," isn't available 	for download yet.\n"))
  	        cat("\tTaking the run dated just before this one.\n")
			dateData <-  dateData[-which(dateData==(summary(dateData)[which(names(summary(dateData))=="Max.")]))]
			lastDate <- dateData[match( summary(dateData)[which(names(summary(dateData))=="Max.")], dateData)]
			#sub back _ symbols
			lastDate <- gsub("-","_",as.character(lastDate))
			lastDateCompress <- gsub("_","",lastDate)
			#need last "/" for it to find this page.
			gdacURL <- paste(gdacURL_orig,dataType,"__",lastDate,"/data/",TCGA_acronym_uppercase,"/",lastDateCompress,"/",sep="")
		
			#now get full dataset name. we just have the http link to the last page we select this dataset from right now.
			urlData <- getURL(gdacURL)
			#regular expressions: need \ to have R recognize any " or \ that's actually in our text
			urlData <- strsplit2(urlData,"href=\\\"")
		
			#did we reach the end of the dates - ie only one left? leave the loop then.
			if(length(dateData)<=1){		
				break	
			}  
		} 

		#cat("Using data from date ",lastDate,"\n")
		#should be OK now!
		if(length(grep("was not found",urlData))>0){  
			#this disease may not even be in the analyses directory yet.
			stop( paste0("\nNot returning any viable url data paths after searching by date for disease ",TCGA_acronym_uppercase,". No data was downloaded.\n"))
		}
	
		#remove any FFPE datasets, or only keep those depending on user inputs.
		if (FFPE) { 
			urlData <- urlData[grep("FFPE",urlData)]	
			if(length(urlData)==0){		
				stop("\nNo FFPE data found for this query. Try FFPE=FALSE.\n")		
			}	  
		} else {	
			#we DON'T want FFPE data.
			#but if no FFPE data to begin with: don't subset on this.
			if(length(grep("FFPE",urlData))>0){		
				urlData <- urlData[-grep("FFPE",urlData)]		
			}
			if(length(urlData)==0){		
				stop("\nNo non-FFPE data found for this query. Try FFPE=TRUE.\n")		
			}
		}
		#now get full dataset name.
		fileName <- urlData[grep(dataFileTag,urlData)]
	
		if(length(fileName)==0){	  
			warnMessage <- paste0("\nNot returning any viable url data paths after searching by date for disease ",TCGA_acronym_uppercase," 	for data type ",dataFileTag ,".No data was downloaded.\n")
			warning(warnMessage)
			return(NA)	  
		}
		#some redundancy..but that' OK because we'll add back on the unique tar.gz file tag.
		#first file is one we want - not md5 file.
		fileName <- strsplit2(fileName,"tar.gz")[1,1]
		fileName <- paste(fileName,fileType,sep="")
	
		#final download url
		gdacURL <- paste(gdacURL,fileName,sep="")
	
		# Beed the savedir when we don't download !!!!!!!!
		saveDir <- paste(saveDir,"gdac_",lastDateCompress,'/',sep="")
	
		if(downloadData){		
			cat("\tDownloading",dataFileTag,"data, version:",lastDate,"\n")				
			cat("\tThis may take 10-60 minutes depending on the size of the data set.\n")
		
			# create dirs
			dir.create(saveDir,showWarnings=FALSE)
			
			# download file		
			setwd(saveDir)				
			download.file(gdacURL,fileName,quiet=TRUE,mode="wb")
			
			#this assumes a tar.gz file.
			if(fileType=="tar.gz" && untarUngzip) {					
				cat("\tUnpacking data.\n")
				tarfile=paste0(saveDir,fileName)
				untar(tarfile)
				  
				#remove tarred file
				fileToRemove <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
				file.remove(paste0(saveDir,fileToRemove))
		
			} else if(untarUngzip) {		
				warning("File expansion/opening only built in for tar.gz files at the moment.\n")		
			}
		
			finalDir <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
			finalDir <- strsplit2(finalDir,fileType)		
			#must remove LAST period (ie laster character) only. 
			finalDir <- substr(finalDir,start=0,stop=(nchar(finalDir)-1))
			finalDir <- paste0(saveDir,finalDir)		
			cat("\tFinished downloading",dataFileTag,"data to",finalDir,"\n")
	  
		} else {
			#just spit out the command you need
			cat("download data url is :\n ",gdacURL,'\n')
			finalDir <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
			finalDir <- strsplit2(finalDir,fileType)
		
			#must remove LAST period (ie laster character) only. 
			finalDir <- substr(finalDir,start=0,stop=(nchar(finalDir)-1))
			finalDir <- paste0(saveDir,finalDir)	
		}
    	DownloadedFile=paste0(finalDir,'/')
    	return(DownloadedFile)
    } else{
    	return(cat(paste0("No data correspond to cancer ",TCGA_acronym_uppercase,"\n")))
    }
}


