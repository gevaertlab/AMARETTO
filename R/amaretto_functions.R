#' AMARETTO_Initialize
#' Code used to initialize the seed clusters for an AMARETTO run. 
#' Requires processed gene expression (rna-seq or microarray), CNV (usually from a GISTIC run), and methylation (from MethylMix, provided in this package) data. 
#' Uses the function CreateRegulatorData() and results are fed into the function AMARETTO_Run().
#' @param MA_matrix Expression matrix, with genes in rows and samples in columns.
#' @param CNV_matrix CNV matrix, with genes in rows and samples in columns.
#' @param MET_matrix Methylation matrix, with genes in rows and samples in columns.
#' @param Driver_list Custom list of driver genes to be considered in analysis
#' @param NrModules How many gene co-expression modules should AMARETTO search for? Usually around 100 is acceptable, given the large number of possible driver-passenger gene combinations.
#' @param VarPercentage Minimum percentage by variance for filtering of genes; for example, 75\% would indicate that the CreateRegulatorData() function only analyses genes that have a variance above the 75th percentile across all samples.
#' @param PvalueThreshold Threshold used to find relevant driver genes with CNV alterations: maximal p-value.
#' @param RsquareThreshold Threshold used to find relevant driver genes with CNV alterations: minimal R-square value between CNV and gene expression data.
#' @param pmax "pmax" variable for glmnet function from glmnet package; the maximum number of variables aver to be nonzero. Should not be changed by user unless she/he fully understands the AMARETTO algorithm and how its parameters choices affect model output.
#' @param NrCores A numeric variable indicating the number of computer/server cores to use for paralellelization. Default is 1, i.e. no parallelization. Please check your computer or server's computing capacities before increasing this number.  Parallelization is done via the RParallel package. Mac vs. Windows environments may behave differently when using parallelization. 
#' @param method Perform union or intersection of the driver genes evaluated from the input data matrices and custom driver gene list provided.
#' @export
#' @import matrixStats 
#' @importFrom matrixStats rowVars
#' @importFrom matrixStats rowMads
#' @examples 
#' 
#' AMARETTOinit <- AMARETTO_Initialize(MA_matrix=MA_matrix, CNV_matrix= CNV_matrix,MET_matrix= MET_matrix, Driver_list = NULL, NrModules= Nr, VarPercentage= Var, 
#'                                     PvalueThreshold = 0.001, RsquareThreshold = 0.1, pmax = 10, 
#'                                     NrCores = 1, OneRunStop = 0)
#'                                     
#'  In absense of either of data matrix and providing a list of driver genes, perform a union/intersection of driver genes evaluated from input data matrices and custom drivers.                                    
#'  AMARETTOinit <- AMARETTO_Initialize(MA_matrix=MA_matrix, CNV_matrix= NULL,MET_matrix= MET_matrix, Driver_list = custom_drivers, NrModules= Nr, VarPercentage= Var, 
#'                                     PvalueThreshold = 0.001, RsquareThreshold = 0.1, pmax = 10, 
#'                                     NrCores = 1, OneRunStop = 0 , method= "union")                                   
AMARETTO_Initialize <- function(MA_matrix=MA_matrix,CNV_matrix=NULL,MET_matrix=NULL, Driver_list = NULL,NrModules,
                                VarPercentage,PvalueThreshold=0.001,RsquareThreshold=0.1,pmax=10,NrCores=1,OneRunStop=0, method= "union"){
  if(is.null(MET_matrix) & is.null(CNV_matrix) & is.null(Driver_list)) {
    stop("Please select the correct input data types")
  }
  if(is.null(CNV_matrix)  & is.null(Driver_list))  
  {
    if (ncol(MA_matrix)<2 || ncol(MET_matrix)==1){
      stop("AMARETTO cannot be run with less than two samples.\n")
    }
  }
  if(is.null(MET_matrix)  & is.null(Driver_list)) {
    if (ncol(MA_matrix)<2 || ncol(CNV_matrix)==1 ){
      stop("AMARETTO cannot be run with less than two samples.\n")
    }
  }
  
  if(!is.null(MA_matrix) & !is.null(CNV_matrix) & !is.null(MET_matrix)  & !is.null(Driver_list)) {
    if (ncol(MA_matrix)<2 || ncol(CNV_matrix)==1 || ncol(MET_matrix)==1){
      stop("AMARETTO cannot be run with less than two samples.\n")
    }
  }
  AutoRegulation=2; Lambda2=0.0001; alpha=1-1e-06
  Parameters <- list(AutoRegulation=AutoRegulation,OneRunStop=OneRunStop,Lambda2=Lambda2,Mode='larsen',pmax=pmax,alpha=alpha)
  RegulatorInfo=CreateRegulatorData(MA_matrix=MA_matrix,CNV_matrix=CNV_matrix,MET_matrix=MET_matrix,Driver_list=Driver_list,PvalueThreshold=PvalueThreshold,RsquareThreshold=RsquareThreshold,method=method)
  if (length(RegulatorInfo)>1){
    RegulatorData=RegulatorInfo$RegulatorData
    Alterations = RegulatorInfo$Alterations
    MA_matrix_Var=geneFiltering('MAD',MA_matrix,VarPercentage)
    MA_matrix_Var=t(scale(t(MA_matrix_Var)))
    if (NrModules>=nrow(MA_matrix_Var)){
      stop(paste0("The number of modules is too large compared to the number of genes. Choose a number of modules smaller than ",nrow(MA_matrix_Var),".\n"))
    }
    KmeansResults=kmeans(MA_matrix_Var,NrModules,iter.max=100)
    Clusters=as.numeric(KmeansResults$cluster)
    names(Clusters) <- rownames(MA_matrix_Var)
    return(list(MA_matrix_Var=MA_matrix_Var,RegulatorData=RegulatorData,RegulatorAlterations=Alterations,ModuleMembership=Clusters,Parameters=Parameters,NrCores=NrCores))
  }
}


#' AMARETTO_Run
#' Function to run AMARETTO, a statistical algorithm to identify cancer drivers by integrating a variety of omics data from cancer and normal tissue.
#' @param AMARETTOinit List output from AMARETTO_Initialize().
#'
#' @return
#' @export
#' @import doParallel
#' @import foreach
#' @import glmnet
#' @import parallel
#' @importFrom doParallel registerDoParallel
#' 
#' @examples AMARETTOresults<-AMARETTO_Run(AMARETTOinit)

AMARETTO_Run <- function(AMARETTOinit) {
  if (length(AMARETTOinit)==0){
    cat('For cancer ',CancerSite,' no drivers were find during the initialization step of AMARETTO')
  } else{
    if (nrow(AMARETTOinit$RegulatorData)==1){
      cat('For cancer ',CancerSite,' only one driver is detected. AMARETTO cannot be run with less than two drivers.\n')
    } else {
      cat('Running AMARETTO on',length(rownames(AMARETTOinit$MA_matrix_Var)),'genes and',length(colnames(AMARETTOinit$MA_matrix_Var)),'samples.\n')
      cat('\tStopping if less then',0.01*length(rownames(AMARETTOinit$MA_matrix_Var)),'genes reassigned.\n')
      result=AMARETTO_LarsenBased(AMARETTOinit$MA_matrix_Var,AMARETTOinit$ModuleMembership,AMARETTOinit$RegulatorData,AMARETTOinit$Parameters,AMARETTOinit$NrCores)
      result$ModuleData=AMARETTO_CreateModuleData(AMARETTOinit,result)
      result$RegulatoryProgramData=AMARETTO_CreateRegulatorPrograms(AMARETTOinit,result)
      return(result)
    }
  }
}

#' AMARETTO_EvaluateTestSet
#' Code to evaluate AMARETTO on a new gene expression test set. Uses output from AMARETTO_Run() and CreateRegulatorData().
#' @param AMARETTOresults AMARETTO output from AMARETTO_Run().
#' @param MA_Data_TestSet Gene expression matrix from a test set (that was not used in AMARETTO_Run()).
#' @param RegulatorData_TestSet Test regulator data from CreateRegulatorData().
#'
#' @return
#' @export
#'
#' @examples AMARETTOtestReport<-AMARETTO_EvaluateTestSet(AMARETTOresults=AMARETTOresults, MA_Data_TestSet=AMARETTOinit$MA_matrix_Var, RegulatorData_TestSet=AMARETTOinit$RegulatorData)
AMARETTO_EvaluateTestSet <- function(AMARETTOresults,MA_Data_TestSet,RegulatorData_TestSet) {
  nrSamples = ncol(MA_Data_TestSet)
  RegulatorNames=rownames(RegulatorData_TestSet)
  stats = mat.or.vec(AMARETTOresults$NrModules,9)
  Rsquare = mat.or.vec(AMARETTOresults$NrModules,1)
  RsquareAjusted = mat.or.vec(AMARETTOresults$NrModules,1)
  modules <- list()
  for (i in 1:AMARETTOresults$NrModules){
    currentRegulators = RegulatorNames[which(AMARETTOresults$RegulatoryPrograms[i,] != 0)]
    nrPresentRegulators = sum((rownames(RegulatorData_TestSet)  %in% currentRegulators))
    currentPresentRegulators = (currentRegulators %in% rownames(RegulatorData_TestSet))
    stats[i,1] = nrPresentRegulators
    stats[i,2] = length(currentRegulators)
    currentClusterGenes = AMARETTOresults$AllGenes[which(AMARETTOresults$ModuleMembership[,1] == i)]
    nrPresentClusterGenes = sum((rownames(MA_Data_TestSet) %in% currentClusterGenes))
    stats[i,3] = nrPresentClusterGenes
    stats[i,4] = length(currentClusterGenes)
    stats[i,5] = stats[i,3] / stats[i,4] * 100
    currentWeights = AMARETTOresults$RegulatoryPrograms[i,which(AMARETTOresults$RegulatoryPrograms[i,] != 0)]
    totalWeights = sum(abs(currentWeights))
    presentWeights = currentWeights[currentPresentRegulators]
    presentRegulators = currentRegulators[currentPresentRegulators]
    totalWeightsPercentage = sum(abs(presentWeights)) / totalWeights * 100
    modules[[i]] = currentClusterGenes[currentClusterGenes %in% rownames(MA_Data_TestSet)]
    if (nrPresentRegulators > 0) {
      predictions = (t(RegulatorData_TestSet[presentRegulators,,drop=FALSE])) %*% (presentWeights) # need to make sure that the first argument remains a matrix.
      predictions = data.matrix(predictions)
      if (length(modules[[i]]) !=0) {
        if (length(currentClusterGenes)>1){
          outcome = apply(MA_Data_TestSet[currentClusterGenes,],2,mean)
        } else {
          outcome = MA_Data_TestSet[currentClusterGenes,]
        }
        residuals = predictions-outcome
        SStot = sum((outcome-mean(outcome))^2)
        SSres = sum((predictions-outcome)^2)
        Rsquare[i] = 1 - (SSres / SStot)
        RsquareAjusted[i] = Rsquare[i] - (nrPresentRegulators/(nrSamples - 1 - nrPresentRegulators))*(1-Rsquare[i])
        MSE = (1/nrSamples) * sum((predictions-outcome)^2)
        
        stats[i,6] = totalWeightsPercentage
        stats[i,7] = Rsquare[i]
        stats[i,8] = RsquareAjusted[i]
        stats[i,9] = MSE
      } else {
        stats[i,6] = totalWeightsPercentage
        stats[i,7] = 0
        stats[i,8] = 0
        stats[i,9] = 0
      }
    } else {
      stats[i,6] = 0
      stats[i,7] = 0
      stats[i,8] = 0
      stats[i,9] = 0
    }
  }
  dimnames(stats) <- list(rownames(stats, do.NULL = FALSE, prefix = "Module_"),
                          c("nrPresReg" ,"nrTotalReg", "nrPresGen", "nrTotGen", "percPresGen",
                            "percWeightPresent", "Rsquare", "RsquareAdjusted","MSE"))
  return(stats)
}