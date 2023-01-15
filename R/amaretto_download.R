#' AMARETTO_Download
#'
#' Downloading TCGA dataset for AMARETTO analysis
#' @param CancerSite TCGA cancer code for data download
#' @param TargetDirectory Directory path to download data
#' @return result
#' @importFrom curatedTCGAData curatedTCGAData
#' @importFrom httr GET stop_for_status
#' @importFrom limma strsplit2
#' @importFrom BiocFileCache BiocFileCache bfcadd bfcquery
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr  everything  mutate  select
#' @importFrom foreach foreach
#' @import grDevices
#' @importFrom parallel makeCluster stopCluster
#' @importFrom readr  write_tsv
#' @importFrom tibble rownames_to_column
#' @importFrom utils  untar zip
#' @export
#' @examples
#' TargetDirectory <- file.path(getwd(),"Downloads/");dir.create(TargetDirectory)
#' CancerSite <- 'CHOL'
#' DataSetDirectories <- AMARETTO_Download(CancerSite,TargetDirectory = TargetDirectory)
AMARETTO_Download <- function(CancerSite = "CHOL", 
    TargetDirectory = TargetDirectory) {
    ori.dir <- getwd()
    message("Downloading Gene Expression and Copy Number Variation data for: ", 
        CancerSite, "\n")
    Cancers = c("BLCA", "BRCA", "LUAD", "LUSC", "COADREAD", 
        "HNSC", "KIRC", "GBM", "OV", "LAML", "UCEC", 
        "COAD", "READ")
    if (!(CancerSite %in% Cancers)) {
        message("This TCGA cancer site/type was not tested, continue at your own risk.\n")
    }
    if (!file.exists(TargetDirectory)) 
        dir.create(TargetDirectory, showWarnings = FALSE)
    TCGA_acronym_uppercase = toupper(CancerSite)
    assays <- c("RNASeq2GeneNorm")
    MAEO <- suppressMessages(curatedTCGAData::curatedTCGAData(
        CancerSite, assays, version = "1.1.38", dry.run = FALSE
    ))
    
    saveRDS(MAEO, file = paste0(TargetDirectory, "/", CancerSite, "_RNASeq_MAEO.rds"))
    
    dataType = "analyses"
    dataFileTag = "CopyNumber_Gistic2.Level_4"
    message("Searching CNV data for:", CancerSite, 
        "\n")
    CNVdirectory = get_firehoseData(TargetDirectory = TargetDirectory, 
        TCGA_acronym_uppercase = TCGA_acronym_uppercase, 
        dataType = dataType, dataFileTag = dataFileTag)
    
    on.exit(setwd(ori.dir))
    return(list(CancerSite = CancerSite, MAdirectory = TargetDirectory, 
        CNVdirectory = CNVdirectory))
    
}

#' get_firehoseData
#'
#' Downloading TCGA dataset via firehose
#' @return result
#' @keywords internal
get_firehoseData <- function(TargetDirectory = "./", 
    TCGA_acronym_uppercase = "LUAD", dataType = "stddata", 
    dataFileTag = "mRNAseq_Preprocess.Level_3", FFPE = FALSE, 
    fileType = "tar.gz", gdacURL = "https://gdac.broadinstitute.org/runs/", 
    untarUngzip = TRUE, printDisease_abbr = FALSE) {
    # Cases Shipped by BCR # Cases with Data* Date Last
    # Updated (mm/dd/yy)
    ori.dir <- getwd()
    cancers <- c("Acute Myeloid Leukemia [LAML] \n", 
        "Adrenocortical carcinoma [ACC]\t\n", "Bladder Urothelial Carcinoma [BLCA] \n", 
        "Brain Lower Grade Glioma [LGG] \n", "Breast invasive carcinoma [BRCA] \n", 
        "Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] \n", 
        "Cholangiocarcinoma [CHOL] \n", "Colon adenocarcinoma [COAD] \n", 
        "Esophageal carcinoma [ESCA] \n", "Glioblastoma multiforme [GBM] \n", 
        "Head and Neck squamous cell carcinoma [HNSC]\t\n", 
        "Kidney Chromophobe [KICH]\t\n", "Kidney renal clear cell carcinoma [KIRC]\t\n", 
        "Kidney renal papillary cell carcinoma [KIRP]\t\n", 
        "Liver hepatocellular carcinoma [LIHC]\t\n", 
        "Lung adenocarcinoma [LUAD]\t\n", "Lung squamous cell carcinoma [LUSC] \n", 
        "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma [DLBC]\t\n", 
        "Mesothelioma [MESO] \n", "Ovarian serous cystadenocarcinoma [OV]\t\n", 
        "Pancreatic adenocarcinoma [PAAD]\t\n", "Pheochromocytoma and Paraganglioma [PCPG] \n", 
        "Prostate adenocarcinoma [PRAD] \n", "Rectum adenocarcinoma [READ]\t\n", 
        "Sarcoma [SARC]\t\n", "Skin Cutaneous Melanoma [SKCM]\t\n", 
        "Stomach adenocarcinoma [STAD] \n", "Testicular Germ Cell Tumors [TGCT] \n", 
        "Thymoma [THYM] \n", "Thyroid carcinoma [THCA]\t\n", 
        "Uterine Carcinosarcoma [UCS]\t \n", "Uterine Corpus Endometrial Carcinoma [UCEC]\t\n", 
        "Uveal Melanoma [UVM] \n")
    
    cancers_acronyms <- c("LAML", "ACC", "BLCA", "LGG", 
        "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
        "HNSC", "KICH", "KIRC", "LIHC", "LUAD", "LUSC", 
        "DLBC", "MESO", "OV", "PAAD", "PCPG", "PRAD", 
        "READ", "SARC", "SKCM", "STAD", "TGCT", "THYM", 
        "THCA", "UCS", "UCEC", "UVM")
    
    if (printDisease_abbr) {
        message(cat("Here are the possible TCGA database disease acronyms. \nRe-run this function with printDisease_abbr=FALSE to then run an actual query.\n\n", 
            cancers))
    }
    if (TCGA_acronym_uppercase %in% cancers_acronyms) {
        
        gdacURL_orig <- gdacURL
        message("\t Printing gdacURL.\n")
        print(gdacURL)
        urlData <- web.lnk <- httr::GET(gdacURL)
        urlData <- limma::strsplit2(urlData, paste(dataType, 
            "__", sep = ""))
        urlData <- urlData[, 2:dim(urlData)[2]]
        urlData <- limma::strsplit2(urlData, "/")
        urlData <- urlData[, 1]
        urlData <- as.POSIXct(strptime(urlData, "%Y_%m_%d"))
        message("\t Printing urlData \n")
        print(urlData)
        dateData <- as.Date(as.character(urlData[which(!is.na(urlData))]))
        print(dateData)
        lastDate <- dateData[match(summary(dateData)[which(names(summary(dateData)) == 
            "Max.")], dateData)]
        lastDate <- gsub("-", "_", as.character(lastDate))
        lastDateCompress <- gsub("_", "", lastDate)
        gdacURL <- paste(gdacURL, dataType, "__", lastDate, 
            "/data/", TCGA_acronym_uppercase, "/", 
            lastDateCompress, "/", sep = "")
        
        urlData <- web.lnk <- httr::GET(gdacURL)
        urlData <- limma::strsplit2(urlData, "href=\\\"")
        while (length(grep("was not found", urlData)) > 
            0) {
            message(paste0("\tNOTE: the TCGA run dated ", 
                lastDate, " for ", TCGA_acronym_uppercase, 
                " isn't available for download yet. \n"))
            message("\tTaking the run dated just before this one.\n")
            dateData <- dateData[-which(dateData == 
                (summary(dateData)[which(names(summary(dateData)) == 
                  "Max.")]))]
            lastDate <- dateData[match(summary(dateData)[which(names(summary(dateData)) == 
                "Max.")], dateData)]
            lastDate <- gsub("-", "_", as.character(lastDate))
            lastDateCompress <- gsub("_", "", lastDate)
            gdacURL <- paste(gdacURL_orig, dataType, 
                "__", lastDate, "/data/", TCGA_acronym_uppercase, 
                "/", lastDateCompress, "/", sep = "")
            urlData <- web.lnk <- httr::GET(gdacURL)
            
            urlData <- limma::strsplit2(urlData, "href=\\\"")
            if (length(dateData) <= 1) {
                break
            }
        }
        httr::stop_for_status(web.lnk, task = "FALIED to download input TCGA data type")
        if (FFPE) {
            urlData <- urlData[grep("FFPE", urlData)]
            if (length(urlData) == 0) {
                stop("\nNo FFPE data found for this query. Try FFPE=FALSE.\n")
            }
        } else {
            if (length(grep("FFPE", urlData)) > 0) {
                urlData <- urlData[-grep("FFPE", urlData)]
            }
            if (length(urlData) == 0) {
                stop("\nNo non-FFPE data found for this query. Try FFPE=TRUE.\n")
            }
        }
        fileName <- urlData[grep(dataFileTag, urlData)]
        if (length(fileName) == 0) {
            warnMessage <- paste0("\nNot returning any viable url data paths after searching by date for disease ", 
                TCGA_acronym_uppercase, " \tfor data type ", 
                dataFileTag, ".No data was downloaded.\n")
            warning(warnMessage)
            return(NA)
        }
        fileName <- limma::strsplit2(fileName, "tar.gz")[1, 
            1]
        fileName <- paste(fileName, fileType, sep = "")
        gdacURL <- paste(gdacURL, fileName, sep = "")
        
        cancer_url <- computeGisticURL(url = gdacURL)
        cache_target <- cacheResource(TargetDirectory=TargetDirectory,resource = cancer_url)
        utils::untar(cache_target$rpath, exdir = TargetDirectory)
        DownloadedFile <- list.dirs(TargetDirectory, 
            full.names = TRUE)[grep(TCGA_acronym_uppercase, list.dirs(TargetDirectory, 
            full.names = TRUE))]
        DownloadedFile <- paste0(DownloadedFile, "/")
        return(DownloadedFile)
    }
    on.exit(setwd(ori.dir))
}


#' AMARETTO_ExportResults
#'
#' Retrieve a download of all the data linked with the run (including heatmaps)
#' @param AMARETTOinit AMARETTO initialize output
#' @param AMARETTOresults AMARETTO results output
#' @param data_address Directory to save data folder
#' @param Heatmaps Output heatmaps as pdf
#' @param CNV_matrix CNV_matrix
#' @param MET_matrix MET_matrix
#' @return result
#' @examples
#' data('ProcessedDataLIHC')
#' TargetDirectory <- file.path(getwd(),"Downloads/");dir.create(TargetDirectory)
#' AMARETTOinit <- AMARETTO_Initialize(ProcessedData = ProcessedDataLIHC,
#'                                     NrModules = 2, VarPercentage = 50)
#' 
#' AMARETTOresults <- AMARETTO_Run(AMARETTOinit)
#' AMARETTO_ExportResults(AMARETTOinit,AMARETTOresults,TargetDirectory,Heatmaps = FALSE)
#' @export
AMARETTO_ExportResults <- function(AMARETTOinit, AMARETTOresults, 
    data_address, Heatmaps = TRUE, CNV_matrix = NULL, 
    MET_matrix = NULL) {
  `%dopar%` <- foreach::`%dopar%`
    if (!dir.exists(data_address)) {
        stop("Output directory is not existing.")
    }
    
    # add a date stamp to the output directory
    output_dir <- paste0("AMARETTOresults_", gsub("-|:", 
        "", gsub(" ", "_", Sys.time())))
    dir.create(file.path(data_address, output_dir))
    
    NrCores <- AMARETTOinit$NrCores
    NrModules <- AMARETTOresults$NrModules
    ModuleNr <- NULL
    
    # parallelize the heatmap production
#    cluster <- parallel::makeCluster(c(rep("localhost", 
#        NrCores)), type = "SOCK")
#    doParallel::registerDoParallel(cluster, cores = NrCores)
    
    if (Heatmaps == TRUE) {
        foreach::foreach(ModuleNr = 1:NrModules, .packages = c("AMARETTO")) %dopar%
            {
                pdf(file = file.path(data_address, 
                  output_dir, paste0("Module_", as.character(ModuleNr), 
                    ".pdf")))
                AMARETTO_VisualizeModule(AMARETTOinit, 
                  AMARETTOresults, CNV_matrix, MET_matrix, 
                  ModuleNr = ModuleNr)
                dev.off()
            }
    }
    
#    parallel::stopCluster(cluster)
    # save rdata files for AMARETTO_Run and
    # AMARETTO_Initialize output
    save(AMARETTOresults, file = file.path(data_address, 
        output_dir, "/amarettoResults.RData"))
    save(AMARETTOinit, file = file.path(data_address, 
        output_dir, "/amarettoInit.RData"))
    
    # save some tables that might be useful for further
    # analysis
    write_gct(AMARETTOresults$ModuleData, file.path(data_address, 
        output_dir, "/ModuleData_amaretto.gct"))
    write_gct(AMARETTOresults$ModuleMembership, file.path(data_address, 
        output_dir, "/ModuleMembership_amaretto.gct"))
    write_gct(AMARETTOresults$RegulatoryProgramData, 
        file.path(data_address, output_dir, "/RegulatoryProgramData_amaretto.gct"))
    write_gct(AMARETTOresults$RegulatoryPrograms, file.path(data_address, 
        output_dir, "/RegulatoryPrograms_amaretto.gct"))
    readr::write_tsv(as.data.frame(AMARETTOresults$AllGenes), 
        file.path(data_address, output_dir, "/AllGenes_amaretto.tsv"))
    readr::write_tsv(as.data.frame(AMARETTOresults$AllRegulators), 
        file.path(data_address, output_dir, "/AllRegulators_amaretto.tsv"))
    readr::write_tsv(as.data.frame(AMARETTOresults$NrModules), 
        file.path(data_address, output_dir, "/NrModules_amaretto.tsv"))
    
    # zip the file
    utils::zip(zipfile = file.path(data_address, output_dir), 
        files = file.path(data_address, output_dir))
}


#' write_gct
#' 
#' @return result
#' @keywords internal
write_gct <- function(data_in, file_address) {
    Name <- Description <- NULL
    header_gct <- paste0("#1.2\n", nrow(data_in), "\t", 
        ncol(data_in))
    data_in <- tibble::rownames_to_column(as.data.frame(data_in), 
        "Name") %>% dplyr::mutate(Description = Name) %>% 
        dplyr::select(Name, Description, dplyr::everything())
    write(header_gct, file = file_address, append = FALSE)
    readr::write_tsv(data_in, file_address, append = TRUE, 
        col_names = TRUE)
}


#' computeGisticURL
#' 
#' @return result
#' @keywords internal
computeGisticURL <- function(url = NULL, acronym = "CHOL") {
    if (!is.null(url)) 
        return(url)
    sprintf("http://gdac.broadinstitute.org/runs/analyses__2016_01_28/data/%s/20160128/gdac.broadinstitute.org_%s-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz", 
        acronym, acronym)
}


#' cacheResource
#' 
#' @return result
#' @keywords internal
cacheResource <- function(TargetDirectory=TargetDirectory, 
    resource = resource) {
    cache = BiocFileCache::BiocFileCache(TargetDirectory)
    chk = bfcquery(cache, resource)
    if (nrow(chk) == 0) {
        message("downloading ", resource)
        BiocFileCache::bfcadd(cache, resource)
        return(bfcquery(cache, resource))
    }
    chk
}

