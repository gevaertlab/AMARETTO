#' CreateRegulatorData
#'
#' Determine potential regulator genes.
#' @return result
#' @keywords internal
CreateRegulatorData <- function(MA_matrix = MA_matrix, CNV_matrix = NULL, MET_matrix = NULL, Driver_list = NULL, PvalueThreshold = 0.001, RsquareThreshold = 0.1, method = "union") {
    if (is.null(Driver_list)) 
        DriversList <- NULL
    if (is.null(CNV_matrix)) 
        CNV_matrix <- matrix(0, nrow = 0, ncol = 0)
    if (nrow(CNV_matrix) > 1) {
        GeneVariances = rowVars(CNV_matrix)
        CNV_matrix = CNV_matrix[GeneVariances >= 1e-04,]
    }
    if (nrow(CNV_matrix) > 0) {
        CNV_matrix = FindTranscriptionallyPredictive_CNV(MA_matrix, CNV_matrix, PvalueThreshold = PvalueThreshold, RsquareThreshold = RsquareThreshold)
    }
    cat("\tFound", length(rownames(CNV_matrix)), "CNV driver genes.\n")
    if (is.null(MET_matrix)) 
        MET_matrix <- matrix(0, nrow = 0, ncol = 0)
    if (nrow(MET_matrix) > 1) {
        METcounts = apply(MET_matrix, 1, function(x) length(unique(x)))
        MET_matrix = MET_matrix[METcounts == 2, ]
        GeneVariances = rowVars(MET_matrix)
        MET_matrix = MET_matrix[GeneVariances >= 1e-04,]
    }
    MET_drivers <- intersect(rownames(MET_matrix),rownames(MA_matrix))
    cat("\tFound", length(MET_drivers), "MethylMix driver genes.\n")
    if (!is.null(Driver_list)) {
        DriversList <- Driver_list
        DriversList <- intersect(DriversList, rownames(MA_matrix))
        cat("\tFound", length(DriversList), "driver genes from the input list.\n")
    }
    if (is.null(Driver_list)) 
        Drivers <- unique(c(rownames(CNV_matrix), MET_drivers,DriversList))
    dataset_drivers <- unique(c(rownames(CNV_matrix),MET_drivers))
    if (!is.null(Driver_list) & method == "union") 
        Drivers <- union(dataset_drivers, DriversList)
    if (!is.null(Driver_list) & method == "intersect") 
        Drivers <- intersect(dataset_drivers, DriversList)
    cat("\tFound a total of", length(Drivers), "unique drivers with your selected method.\n")
    if (length(Drivers) == 0) {
        cat("AMARETTO doesn't find any driver genes.")
        return("No driver")
    } else {
        if (length(Drivers) == 1) {
            RegulatorData_temp <- matrix(0, 1, ncol(MA_matrix))
            colnames(RegulatorData_temp) <- colnames(MA_matrix)
            rownames(RegulatorData_temp) <- Drivers
            RegulatorData_temp[1, ] <- MA_matrix[Drivers,]
            RegulatorData <- RegulatorData_temp
        } else {
            RegulatorData = MA_matrix[Drivers, ]
        }
        RegulatorData = t(scale(t(RegulatorData)))
        MET_aberrations <- matrix(0, ncol = 3, nrow = length(MET_drivers))
        rownames(MET_aberrations) <- MET_drivers
        colnames(MET_aberrations)<-c("Hyper-methylated","Hypo-methylated","No_change")
        if (length(MET_drivers) > 0) {
            MET_aberrations[, "Hyper-methylated"] <- rowSums(MET_matrix[MET_drivers,] > 0)/ncol(MET_matrix)
            MET_aberrations[, "Hypo-methylated"] <- rowSums(MET_matrix[MET_drivers,] < 0)/ncol(MET_matrix)
            MET_aberrations[, "No_change"] <- rowSums(MET_matrix[MET_drivers,] == 0)/ncol(MET_matrix)
        }
        CNV_alterations <- matrix(0, ncol = 3, nrow = nrow(CNV_matrix))
        rownames(CNV_alterations) <- rownames(CNV_matrix)
        colnames(CNV_alterations)<-c("Amplification","Deletion","No_change")
        if (nrow(CNV_alterations) > 0) {
            CNV_alterations[, "Amplification"] <- rowSums(CNV_matrix > 0)/ncol(CNV_matrix)
            CNV_alterations[, "Deletion"] <- rowSums(CNV_matrix < 0)/ncol(CNV_matrix)
            CNV_alterations[, "No_change"] <- rowSums(CNV_matrix == 0)/ncol(CNV_matrix)
        }
        driverList_alterations <- matrix(1, ncol = 1, nrow = length(DriversList))
        rownames(driverList_alterations) <- DriversList
        Alterations <- matrix(0, nrow = length(Drivers), ncol = 3)
        colnames(Alterations) <- c("CNV", "MET", "Driver List")
        rownames(Alterations) <- Drivers
        Alterations[which(Drivers %in% rownames(CNV_matrix)), 1] <- rep(1, length(which(Drivers %in% rownames(CNV_matrix))))
        Alterations[which(Drivers %in% rownames(MET_matrix)), 2] <- rep(1, length(which(Drivers %in% rownames(MET_matrix))))
        Alterations[which(Drivers %in% rownames(driverList_alterations)), 3] <- rep(1, length(which(Drivers %in% rownames(driverList_alterations))))
        Alterations <- list(MET = MET_aberrations, 
                            CNV = CNV_alterations,
                            Driver_list = driverList_alterations, 
                            Summary = Alterations)
        return(list(RegulatorData = RegulatorData, Alterations = Alterations))
    }
}

#' FindTranscriptionallyPredictive_CNV
#'
#' Function to identify which genes CNV significantly predict expression of that gene.
#' @return result
#' @keywords internal
FindTranscriptionallyPredictive_CNV <- function(MA_matrix, 
    CNV_matrix, PvalueThreshold = 0.001, RsquareThreshold = 0.1) {
    OverlapGenes = Reduce(intersect, list(rownames(MA_matrix), 
        rownames(CNV_matrix)))
    OverlapSamples = Reduce(intersect, list(colnames(MA_matrix), 
        colnames(CNV_matrix)))
    if (length(OverlapGenes) == 1) {
        CNV_TCGA_temp = MA_matrix_temp <- matrix(0, 
            length(OverlapGenes), length(OverlapSamples))
        rownames(CNV_TCGA_temp) = rownames(MA_matrix_temp) <- OverlapGenes
        colnames(CNV_TCGA_temp) = colnames(MA_matrix_temp) <- OverlapSamples
        CNV_TCGA_temp[1, ] <- CNV_matrix[OverlapGenes, 
            OverlapSamples]
        MA_matrix_temp[1, ] <- MA_matrix[OverlapGenes, 
            OverlapSamples]
        CNV_matrix = CNV_TCGA_temp
        MA_matrix = MA_matrix_temp
    } else {
        CNV_matrix = CNV_matrix[OverlapGenes, OverlapSamples]
        MA_matrix = MA_matrix[OverlapGenes, OverlapSamples]
    }
    if (length(OverlapGenes) > 0 && length(OverlapSamples) > 
        0) {
        CNVdrivers = c()
        for (i in 1:length(rownames(CNV_matrix))) {
            res = lm(MA_matrix[i, ] ~ CNV_matrix[i, 
                ])
            res.summary = summary(res)
            if (res$coefficients[2] > 0 & res.summary$coefficients[2, 
                4] < PvalueThreshold & res.summary$r.squared > 
                RsquareThreshold) {
                CNVdrivers = c(CNVdrivers, rownames(CNV_matrix)[i])
            }
        }
        if (length(CNVdrivers) == 1) {
            CNV_matrix_temp <- matrix(0, 1, ncol(CNV_matrix))
            colnames(CNV_matrix_temp) = colnames(CNV_matrix)
            rownames(CNV_matrix_temp) = CNVdrivers
            CNV_matrix_temp[1, ] <- CNV_matrix[CNVdrivers, 
                ]
            CNV_matrix = CNV_matrix_temp
        } else {
            CNV_matrix = CNV_matrix[CNVdrivers, ]
        }
    }
    return(CNV_matrix = CNV_matrix)
}

#' geneFiltering
#'
#' Function to filter gene expression matrix
#' @return result
#' @keywords internal
geneFiltering <- function(Type, MAdata, Percentage) {
    switch(Type, Variance = {
        GeneVariances = rowVars(MAdata)
        tmpResult = sort(GeneVariances, decreasing = TRUE)
        SortedGenes = names(tmpResult)
        tmpNrGenes = round(length(rownames(MAdata)) * 
            Percentage/100)
        MAdata_Filtered = MAdata[SortedGenes[1:tmpNrGenes], 
            ]
    }, MAD = {
        GeneVariances = rowMads(MAdata)
        names(GeneVariances) = rownames(MAdata)
        tmpResult = sort(GeneVariances, decreasing = TRUE)
        SortedGenes = names(tmpResult)
        tmpNrGenes = round(length(rownames(MAdata)) * 
            Percentage/100)
        MAdata_Filtered = MAdata[SortedGenes[1:tmpNrGenes], 
            ]
    })
    return(MAdata_Filtered)
}

#' read_gct
#' 
#' Function to turn a .gct data files into a matrix format
#'
#' @param file_address 
#'
#' @importFrom utils read.delim
#' @return result
#' @export
#' @examples
#' data_matrix<-read_gct(file_address="")
read_gct <- function(file_address){
  if(file_address==""){
    print("No gct file address is provided.")
    return(NULL)
  }
  else{
    data_fr <-read.delim(file_address, skip=2, sep="\t", header=TRUE, row.names=1)
    data_fr <- as.matrix(subset(data_fr, select=-c(Description)))
    return(data_fr)
  }
}

#' printf
#'
#' Wrapper function for C-style formatted output.
#' @return result
#' @keywords internal
printf <- function(...) {
    cat(sprintf(...))
}
