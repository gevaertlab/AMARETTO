#' AMARETTO_LarsenBased
#' 
#' @import MultiAssayExperiment
#' @import graphics
#' @return result
#' @keywords internal
AMARETTO_LarsenBased <- function(Data, Clusters, RegulatorData, Parameters, NrCores) {
    registerDoParallel(cores = NrCores)
    ptm1 <- proc.time()
    RegulatorData_rownames = rownames(RegulatorData)
    Data_rownames = rownames(Data)
    AutoRegulation = Parameters$AutoRegulation
    RegulatorSign = array(0, length(RegulatorData_rownames))
    Lambda = Parameters$Lambda2
    OneRunStop = Parameters$OneRunStop
    if (AutoRegulation == 1) {
        cat("\tAutoregulation is turned ON.\n")
    } else if (AutoRegulation == 2) {
        cat("\tAutoregulation is turned ON.\n")
    } else {
        cat("\tAutoregulation is turned OFF.\n")
    }
    NrReassignGenes = length(Data_rownames)
    NrReassignGenes_history <- NrReassignGenes
    error_history<-list()
    index=1
    while (NrReassignGenes > 0.01 * length(Data_rownames)) {
        ptm <- proc.time()
        switch(Parameters$Mode, larsen = {
            regulatoryPrograms <- AMARETTO_LearnRegulatoryProgramsLarsen(Data, Clusters, RegulatorData, RegulatorSign, Lambda, AutoRegulation, alpha = Parameters$alpha, pmax = Parameters$pmax)
        })
        error_history[[index]]<-regulatoryPrograms$error
        index<-index+1
        ptm <- proc.time() - ptm
        printf("Elapsed time is %f seconds\n", ptm[3])
        NrClusters = length(unique(Clusters))
        sum = 0
        for (i in 1:NrClusters) {
            sum = sum + Matrix::nnzero(regulatoryPrograms$Beta[i,])
        }
        avg = sum/NrClusters
        printf("Average nr of regulators per module: %f \n", avg)
        PreviousClusters = Clusters
        if (OneRunStop == 1) {
            break
        }
        ptm <- proc.time()
        ReassignGenesToClusters <- AMARETTO_ReassignGenesToClusters(Data, RegulatorData, regulatoryPrograms$Beta, Clusters, AutoRegulation)
        ptm <- proc.time() - ptm
        printf("Elapsed time is %f seconds\n", ptm[3])
        NrReassignGenes = ReassignGenesToClusters$NrReassignGenes
        Clusters = ReassignGenesToClusters$Clusters
        printf("Nr of reassignments is: %i \n", NrReassignGenes)
        NrReassignGenes_history <- c(NrReassignGenes_history, NrReassignGenes)
    }
    ptm1 <- proc.time() - ptm1
    printf("Elapsed time is %f seconds\n", ptm1[3])
    ModuleMembership = as.matrix(PreviousClusters)
    rownames(ModuleMembership) = rownames(Data)
    colnames(ModuleMembership) = c("ModuleNr")
    result <- list(NrModules = length(unique(Clusters)), RegulatoryPrograms = regulatoryPrograms$Beta, AllRegulators = rownames(RegulatorData),
                   AllGenes = rownames(Data), ModuleMembership = ModuleMembership, AutoRegulationReport = regulatoryPrograms$AutoRegulationReport,
                   run_history = list(NrReassignGenes_history = NrReassignGenes_history, error_history = error_history))
    return(result)
}

#' AMARETTO_LearnRegulatoryProgramsLarsen
#'
#' @return result
#' @keywords internal
AMARETTO_LearnRegulatoryProgramsLarsen <- function(Data, Clusters, RegulatorData, RegulatorSign, Lambda, AutoRegulation, alpha, pmax) {
    `%dopar%` <- foreach::`%dopar%`
    RegulatorData_rownames = rownames(RegulatorData)
    Data_rownames = rownames(Data)
    trace = 0
    NrFolds = 10
    NrClusters = length(unique(Clusters))
    NrGenes = nrow(Data)
    NrSamples = ncol(Data)
    NrInterpolateSteps = 100
    if (AutoRegulation >= 1) {}
    else if (AutoRegulation == 0) {
        BetaSpecial = list(NrClusters, 1)
        RegulatorPositions = list(NrClusters, 1)
    }
    y_all = mat.or.vec(NrClusters, NrSamples)
    ClusterIDs = unique(Clusters)
    ClusterIDs = sort(ClusterIDs, decreasing = FALSE)
    cnt <- 1:NrClusters
    ptm1 <- proc.time()
    BetaY_all <- foreach(i = 1:NrClusters, .combine = cbind, .init = list(list(), list(), list()), .packages = "glmnet") %dopar% {
            if (length(which(Clusters == ClusterIDs[i])) > 1) {
                y = apply((Data[which(Clusters == ClusterIDs[i]),]), 2, mean)
            } else {
                y = Data[which(Clusters == ClusterIDs[i]),]
            }
            CurrentClusterPositions = which(Clusters %in% ClusterIDs[i])
            nrGenesInClusters = length(CurrentClusterPositions)
            if (AutoRegulation >= 1) {
                X = RegulatorData
            } else if (AutoRegulation == 0) {
                X = RegulatorData[setdiff(RegulatorData_rownames, Data_rownames[CurrentClusterPositions]),]
            }
            suppressWM = function(...) suppressWarnings(suppressMessages(...))
            fit = suppressWM(cv.glmnet(t(X), y, alpha = alpha, pmax = pmax, lambda = Lambda_Sequence(t(X), y)))
            nonZeroLambdas <- fit$lambda[which(fit$nzero > 0)]
            nonZeroCVMs <- fit$cvm[which(fit$nzero > 0)]
            if (length(which(nonZeroCVMs == min(nonZeroCVMs, na.rm = TRUE))) == 0) {
                warnMessage <- paste0("\nOn cluster ", i, " there were no cv.glm results that gave non-zero coefficients.")
                message(warnMessage)
            }
            bestNonZeroLambda <- nonZeroLambdas[which(nonZeroCVMs == min(nonZeroCVMs, na.rm = TRUE))]
            b_o = coef(fit, s = bestNonZeroLambda)
            b_opt <- c(b_o[2:length(b_o)])
            if (AutoRegulation == 2) {
                CurrentUsedRegulators = RegulatorData_rownames[which(b_opt != 0, arr.ind = TRUE)]
                CurrentClusterMembers = Data_rownames[CurrentClusterPositions]
                nrIterations = 0
                while (length(CurrentClusterMembers[CurrentClusterMembers %in% CurrentUsedRegulators]) != 0) {
                  CurrentClusterMembers = setdiff(CurrentClusterMembers, CurrentUsedRegulators)
                  nrCurrentClusterMembers = length(CurrentClusterMembers)
                  if (nrCurrentClusterMembers > 0) {
                    names = Data_rownames %in% CurrentClusterMembers
                    if (length(which(names == TRUE)) > 1) {
                      y = apply((Data[names, ]), 2, mean)
                    } else {
                      y = Data[names, ]
                    }
                    fit = suppressWM(cv.glmnet(t(X), y, alpha = alpha, pmax = pmax, lambda = Lambda_Sequence(t(X), y)))
                    nonZeroLambdas <- fit$lambda[which(fit$nzero > 0)]
                    nonZeroCVMs <- fit$cvm[which(fit$nzero > 0)]
                    if (length(which(nonZeroCVMs == min(nonZeroCVMs, na.rm = TRUE))) == 0) {
                      warnMessage <- paste0("\nOn cluster ", i, " there were no cv.glm results that gave non-zero coefficients during the Autoregulation step.")
                      message(warnMessage)
                    }
                    bestNonZeroLambda <- nonZeroLambdas[which(nonZeroCVMs == min(nonZeroCVMs, na.rm = TRUE))]
                    new_b_o = coef(fit, s = bestNonZeroLambda)
                    new_b_opt <- c(new_b_o[2:length(b_o)])
                    CurrentUsedRegulators = RegulatorData_rownames[which(new_b_opt != 0)]
                    nrIterations = nrIterations + 1
                    b_opt = new_b_opt
                  } else {
                    b_opt = rep(0, length(RegulatorData_rownames))
                  }
                }
                Report <- c(length(CurrentClusterPositions), length(CurrentClusterMembers), nrIterations)
            }
            if (sum(RegulatorSign[which(RegulatorSign != 0)]) > 0) {
                RegulatorCheck = RegulatorSign * t(b_opt)
                WrongRegulators = which(RegulatorCheck < 0)
                if (length(WrongRegulators) == 0) {
                  b_opt[WrongRegulators] = 0
                }
            }
            if (AutoRegulation >= 1) {
            } else {
                BetaSpecial[i] = b_opt
                RegulatorPositions[i] = (RegulatorData_rownames %in% setdiff(RegulatorData_rownames, Data_rownames[CurrentClusterPositions]))
            }
            list(b_opt, y, Report)
        }
    if (AutoRegulation == 0) {
        for (i in 1:NrClusters) {
            Beta[i, RegulatorPositions[i]] = BetaSpecial[i]
        }
    }
    tmpPos = NrClusters + 1
    Beta <- do.call(cbind, BetaY_all[1, 2:tmpPos])
    Beta = t(Beta)
    colnames(Beta) = RegulatorData_rownames
    rownames(Beta) = gsub("result.", "Module_", rownames(Beta))
    y_all <- do.call(cbind, BetaY_all[2, 2:tmpPos])
    y_all = t(y_all)
    rownames(y_all) = gsub("result.", "Module_", rownames(y_all))
    AutoRegulationReport <- do.call(cbind, BetaY_all[3, 2:tmpPos])
    AutoRegulationReport = t(AutoRegulationReport)
    rownames(AutoRegulationReport) = gsub("result.", "Module_", rownames(AutoRegulationReport))
    error = y_all - (Beta %*% RegulatorData)
    result <- list(Beta = Beta, error = error, AutoRegulationReport = AutoRegulationReport)
    return(result)
}

#' AMARETTO_ReassignGenesToClusters
#'
#' @return result
#' @importFrom Matrix nnzero
#' @keywords internal
AMARETTO_ReassignGenesToClusters <- function(Data, RegulatorData, Beta, Clusters, AutoRegulation) {
    `%dopar%` <- foreach::`%dopar%`
    RegulatorData_rownames = rownames(RegulatorData)
    Data_rownames = rownames(Data)
    NrGenes = nrow(Data)
    NrSamples = ncol(Data)
    NrReassignGenes = 0
    X = RegulatorData
    X1 = data.matrix(X)
    ModuleVectors = Beta %*% X1
    GeneNames = rownames(Data)
    ptm1 <- proc.time()
    i <- NULL
    nc <- foreach(i = 1:NrGenes, .combine = c) %dopar%
        {
            OldModule = Clusters[i]
            CurrentGeneVector = Data[i, , drop = FALSE]
            Correlations = cor(t(CurrentGeneVector), t(ModuleVectors))
            corr = data.matrix(Correlations, rownames.force = NA)
            MaxCorrelation = max(corr, na.rm = TRUE)
            MaxPosition = which(signif(corr, digits = 7) == signif(MaxCorrelation, digits = 7))
            MaxPosition = MaxPosition[1]
            if (AutoRegulation > 0) {
                if (MaxPosition != OldModule) {
                  NrReassignGenes = NrReassignGenes + 1
                }
                NewClusters = MaxPosition
            } else {
                if (nnzero(rownames(RegulatorData_rownames) %in% GeneNames[i]) != 0) {
                  if (nnzero(which(which(GeneNames %in% rownames(RegulatorData_rownames)) %in% i) %in% which(Beta[MaxPosition,] != 0)) != 0) {
                    if (MaxPosition != OldModule) {
                      NrReassignGenes = NrReassignGenes + 1
                    }
                    NewClusters = MaxPosition
                  } else {
                    NewClusters = OldModule
                  }
                } else {
                  if (MaxPosition != OldModule) {
                    NrReassignGenes = NrReassignGenes + 1
                  }
                  NewClusters = MaxPosition
                }
            }
        }
    ptm1 <- proc.time() - ptm1
    NrReassignGenes = length(which(nc != Clusters))
    result <- list(NrReassignGenes = NrReassignGenes, Clusters = nc)
    return(result)
}

#' AMARETTO_CreateModuleData
#'
#' @param AMARETTOinit List output from AMARETTO_Initialize().
#' @param AMARETTOresults List output from AMARETTO_Run()
#'
#' @return result
#' @export
#' @examples
#' data('ProcessedDataLIHC')
#' AMARETTOinit <- AMARETTO_Initialize(ProcessedData = ProcessedDataLIHC,
#'                                     NrModules = 2, VarPercentage = 50)
#' AMARETTOresults <- AMARETTO_Run(AMARETTOinit)
#' AMARETTO_MD <- AMARETTO_CreateModuleData(AMARETTOinit, AMARETTOresults)
AMARETTO_CreateModuleData <- function(AMARETTOinit, AMARETTOresults) {
    ModuleData = matrix(0, AMARETTOresults$NrModules, length(colnames(AMARETTOinit$MA_matrix_Var)))
    rownames(ModuleData) = rownames(AMARETTOresults$AutoRegulationReport)
    colnames(ModuleData) = colnames(AMARETTOinit$MA_matrix_Var)
    for (ModuleNr in 1:AMARETTOresults$NrModules) {
        currentModuleData = AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership[, 1] == ModuleNr, ]
        if (length(which(AMARETTOresults$ModuleMembership[, 1] == ModuleNr)) > 1) {
            ModuleData[ModuleNr, ] = colMeans(currentModuleData)
        } else {
            ModuleData[ModuleNr, ] = currentModuleData
        }
    }
    return(ModuleData)
}

#' AMARETTO_CreateRegulatorPrograms
#'
#' @param AMARETTOinit  List output from AMARETTO_Initialize().
#' @param AMARETTOresults List output from AMARETTO_Run()
#'
#' @return result
#' @export
#' @examples
#' data('ProcessedDataLIHC')
#' AMARETTOinit <- AMARETTO_Initialize(ProcessedData = ProcessedDataLIHC,
#'                                     NrModules = 2, VarPercentage = 50)
#' AMARETTOresults <- AMARETTO_Run(AMARETTOinit)
#' AMARETTO_RP <- AMARETTO_CreateRegulatorPrograms(AMARETTOinit,AMARETTOresults)
AMARETTO_CreateRegulatorPrograms <- function(AMARETTOinit, AMARETTOresults) {
    RegulatorProgramData = matrix(0, AMARETTOresults$NrModules, length(colnames(AMARETTOinit$MA_matrix_Var)))
    rownames(RegulatorProgramData) = rownames(AMARETTOresults$AutoRegulationReport)
    colnames(RegulatorProgramData) = colnames(AMARETTOinit$MA_matrix_Var)
    RegulatorNames = rownames(AMARETTOinit$RegulatorData)
    for (ModuleNr in 1:AMARETTOresults$NrModules) {
        currentRegulators = RegulatorNames[which(AMARETTOresults$RegulatoryPrograms[ModuleNr, ] != 0)]
        weights = AMARETTOresults$RegulatoryPrograms[ModuleNr, currentRegulators]
        RegulatorData = AMARETTOinit$RegulatorData[currentRegulators, ]
        RegulatorProgramData[ModuleNr, ] = weights %*% RegulatorData
    }
    return(RegulatorProgramData)
}


#' Lambda_Sequence
#'
#' @return result
#' @keywords internal
Lambda_Sequence <- function(sx, sy) {
    n <- nrow(sx)
    lambda_max <- max(abs(colSums(sx * sy)))/n
    epsilon <- 1e-04
    K <- 100
    lambdaseq <- round(exp(seq(log(lambda_max), log(lambda_max * epsilon), length.out = K)), digits = 10)
    return(lambdaseq)
}
