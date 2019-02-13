#' CreateRegulatorData
#'
#' Determine potential regulator genes.
#' Function to identify which genes' CNV and/or methylation statuses significantly predict expression of that gene. Input is output from Preprocess_CancerSite(). This function is used by AMARETTO_Initialize().
#' @param MA_matrix Processed gene expression data matrix ($MA_matrix element from Preprocess_CancerSite() list output).
#' @param CNV_matrix Processed CNV matrix ($CNV_matrix element from Preprocess_CancerSite() list output).
#' @param MET_matrix Processed methylation matrix ($MET_matrix element from Preprocess_CancerSite() list output).
#' @param Driver_list Custom list of driver genes.
#' @param PvalueThreshold Threshold used to find relevant driver genes with CNV alterations: maximal p-value.
#' @param RsquareThreshold Threshold used to find relevant driver genes with CNV alterations: minimal R-square value between CNV and gene expression data.
#' @param method Perform union or intersection of the driver genes evaluated from the input data matrices and custom driver gene list provided.
#' @return result
#' @export
#' @keywords internal
#' @examples
#' data('ProcessedDataLIHC')
#' CreateRegulatorData(MA_matrix = ProcessedDataLIHC$MA_matrix,
#'                     CNV_matrix = ProcessedDataLIHC$CNV_matrix,
#'                      MET_matrix= ProcessedDataLIHC$MET_matrix,
#'                      PvalueThreshold = 0.001, RsquareThreshold = 0.1)
#'
CreateRegulatorData <- function(MA_matrix = MA_matrix, 
    CNV_matrix = NULL, MET_matrix = NULL, Driver_list = NULL, 
    PvalueThreshold = 0.001, RsquareThreshold = 0.1, 
    method = "union") {
    if (is.null(Driver_list)) 
        DriversList <- NULL
    if (is.null(CNV_matrix)) 
        CNV_matrix <- matrix(0, nrow = 0, ncol = 0)
    if (nrow(CNV_matrix) > 1) {
        GeneVariances = rowVars(CNV_matrix)
        CNV_matrix = CNV_matrix[GeneVariances >= 1e-04, 
            ]
    }
    if (nrow(CNV_matrix) > 0) {
        CNV_matrix = FindTranscriptionallyPredictive_CNV(MA_matrix, 
            CNV_matrix, PvalueThreshold = PvalueThreshold, 
            RsquareThreshold = RsquareThreshold)
    }
    cat("\tFound", length(rownames(CNV_matrix)), "CNV driver genes.\n")
    if (is.null(MET_matrix)) 
        MET_matrix <- matrix(0, nrow = 0, ncol = 0)
    if (nrow(MET_matrix) > 1) {
        METcounts = apply(MET_matrix, 1, function(x) length(unique(x)))
        MET_matrix = MET_matrix[METcounts == 2, ]
        GeneVariances = rowVars(MET_matrix)
        MET_matrix = MET_matrix[GeneVariances >= 1e-04, 
            ]
    }
    MET_drivers <- intersect(rownames(MET_matrix), 
        rownames(MA_matrix))
    cat("\tFound", length(MET_drivers), "MethylMix driver genes.\n")
    if (!is.null(Driver_list)) {
        DriversList <- Driver_list
        DriversList <- intersect(DriversList, rownames(MA_matrix))
        cat("\tFound", length(DriversList), "driver genes from the input list.\n")
    }
    if (is.null(Driver_list)) 
        Drivers <- unique(c(rownames(CNV_matrix), MET_drivers, 
            DriversList))
    dataset_drivers <- unique(c(rownames(CNV_matrix), 
        MET_drivers))
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
            RegulatorData_temp[1, ] <- MA_matrix[Drivers, 
                ]
            RegulatorData <- RegulatorData_temp
        } else {
            RegulatorData = MA_matrix[Drivers, ]
        }
        RegulatorData = t(scale(t(RegulatorData)))
        MET_aberrations <- matrix(0, ncol = 3, nrow = length(MET_drivers))
        colnames(MET_aberrations) <- c("Hypo-methylated", 
            "No_change", "Hyper-methylated")
        rownames(MET_aberrations) <- MET_drivers
        if (length(MET_drivers) > 0) {
            MET_aberrations[, 1] <- rowSums(MET_matrix[MET_drivers, 
                ] > 0)/ncol(MET_matrix)
            MET_aberrations[, 2] <- rowSums(MET_matrix[MET_drivers, 
                ] < 0)/ncol(MET_matrix)
            MET_aberrations[, 3] <- rowSums(MET_matrix[MET_drivers, 
                ] == 0)/ncol(MET_matrix)
        }
        CNV_alterations <- matrix(0, ncol = 3, nrow = nrow(CNV_matrix))
        colnames(CNV_alterations) <- c("Amplification", 
            "No_change", "Deletion")
        rownames(CNV_alterations) <- rownames(CNV_matrix)
        if (nrow(CNV_alterations) > 0) {
            CNV_alterations[, 1] <- rowSums(CNV_matrix > 
                0)/ncol(CNV_matrix)
            CNV_alterations[, 2] <- rowSums(CNV_matrix < 
                0)/ncol(CNV_matrix)
            CNV_alterations[, 3] <- rowSums(CNV_matrix == 
                0)/ncol(CNV_matrix)
        }
        driverList_alterations <- matrix(1, ncol = 1, 
            nrow = length(DriversList))
        rownames(driverList_alterations) <- DriversList
        Alterations <- matrix(0, nrow = length(Drivers), 
            ncol = 3)
        colnames(Alterations) <- c("CNV", "MET", "Driver List")
        rownames(Alterations) <- Drivers
        Alterations[which(Drivers %in% rownames(CNV_matrix)), 
            1] <- rep(1, length(which(Drivers %in% 
            rownames(CNV_matrix))))
        Alterations[which(Drivers %in% rownames(MET_matrix)), 
            2] <- rep(1, length(which(Drivers %in% 
            rownames(MET_matrix))))
        Alterations[which(Drivers %in% rownames(driverList_alterations)), 
            3] <- rep(1, length(which(Drivers %in% 
            rownames(driverList_alterations))))
        Alterations <- list(MET = MET_aberrations, 
            CNV = CNV_alterations, Driver_list = driverList_alterations, 
            Summary = Alterations)
        return(list(RegulatorData = RegulatorData, 
            Alterations = Alterations))
    }
}

#' FindTranscriptionallyPredictive_CNV
#'
#' Function to identify which genes CNV significantly predict expression of that gene.
#' @param MA_matrix Processed gene expression data matrix ($MA_matrix element from Preprocess_CancerSite() list output).
#' @param CNV_matrix Processed CNV matrix ($CNV_matrix element from Preprocess_CancerSite() list output).
#' @param PvalueThreshold Threshold used to find relevant driver genes with CNV alterations: maximal p-value.
#' @param RsquareThreshold Threshold used to find relevant driver genes with CNV alterations: minimal R-square value between CNV and gene expression data.
#'
#' @return result
#' @keywords internal
#'
#' @examples
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
#' @param Type
#' @param MAdata Processed gene expression data matrix ($MA_matrix element from Preprocess_CancerSite() list output).
#' @param Percentage Minimum percentage by variance for filtering of genes; for example, 75\% would indicate that the CreateRegulatorData() function only analyses genes that have a variance above the 75th percentile across all samples.
#'
#' @return result
#' @keywords internal
#'
#' @examples
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

#' printf
#'
#' Wrapper function for C-style formatted output.
#'
#' @param ...
#' @return result
#' @keywords internal
#' @examples
printf <- function(...) {
    cat(sprintf(...))
}
