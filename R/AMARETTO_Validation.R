#' AMARETTO_Validation
#'
#'Requires processed gene expressiosn (rna-seq or microarray), CNV (usually from a GISTIC run), and methylation (from MethylMix, provided in this package) data.
#' Uses the functions AMARETTO_Initialize() and AMARETTO_Run() and AMARETTO_EvaluateTestSet()
#' MA_matrix Processed gene expression data matrix ($MA_matrix element from Preprocess_CancerSite() list output).
#' CNV_matrix Processed CNV matrix ($CNV_matrix element from Preprocess_CancerSite() list output).
#' MET_matrix Processed methylation matrix ($MET_matrix element from Preprocess_CancerSite() list output).
#' NrModules How many gene co-expression modules should AMARETTO search for? Usually around 100 is acceptable, given the large number of possible driver-passenger gene combinations.
#' VarPercentage Minimum percentage by variance for filtering of genes; for example, 75\% would indicate that the CreateRegulatorData() function only analyses genes that have a variance above the 75th percentile across all samples.
#' split The ratio of the data set to be used for training, example 0.8 means 80% of the data set will be used for taining and 20% for validation, default=0.8



AMARETTO_Validation <- function(MA_matrix, CNV_matrix, MET_matrix, split=0.8, NrModules,
                                VarPercentage, NrCores=10,  PvalueThreshold = 0.001,
                                RsquareThreshold = 0.1,
                                method = "union", Driver_list = NULL)
{
  library(caret)
  
  
  trainIndex <- createDataPartition(MA_matrix[1,], p = split, 
                                    list = FALSE, 
                                    times = 1)
  
  sample_names<-colnames(MA_matrix)[trainIndex]
  not_sample_names<- setdiff(colnames(MA_matrix),sample_names)
  
  
  MA_matrix_train <- MA_matrix[,sample_names]
  
  MA_matrix_val<- MA_matrix[,not_sample_names]
  
  
  CNV_matrix_train <- CNV_matrix[,intersect(sample_names,colnames(CNV_matrix))]
  
  CNV_matrix_val<- CNV_matrix[, intersect(not_sample_names,colnames(CNV_matrix))]
  
  
  
  MET_matrix_train <- MET_matrix[,intersect(sample_names,colnames(MET_matrix))]
  
  MET_matrix_val<- MET_matrix[,intersect(not_sample_names,colnames(MET_matrix))]
  
  
  
  ####
  
  
  AMARETTOinit <-
    AMARETTO::AMARETTO_Initialize(
      MA_matrix = MA_matrix_train,
      CNV_matrix = CNV_matrix_train,
      MET_matrix = MET_matrix_train,
      NrModules = NrModules,
      VarPercentage = VarPercentage,  NrCores=10
    )
  
  
  
  ########################################################
  # Running AMARETTO
  ########################################################
  AMARETTOresults <- AMARETTO::AMARETTO_Run(AMARETTOinit)
  
  ########################################################
  # Evaluate AMARETTO Results
  ########################################################
  
  
  
  
  
  RegulatorInfo = AMARETTO::CreateRegulatorData(
    MA_matrix = MA_matrix_val,
    CNV_matrix = CNV_matrix_val,
    MET_matrix = MET_matrix_val,
    Driver_list = Driver_list,
    PvalueThreshold = PvalueThreshold,
    RsquareThreshold = RsquareThreshold,
    method = method
  )
  
  
  
  AMARETTO_valReport <- AMARETTO_Evaluate_Validation(AMARETTOresults,
                                                        MA_matrix_val,
                                                        RegulatorInfo$RegulatorData)
  
  AMARETTO_trainReport <- AMARETTO_Evaluate_Validation(AMARETTOresults,
                                                       AMARETTOinit$MA_matrix_Var,
                                                       AMARETTOinit$RegulatorData)
  
  
  return(list(AMARETTO_valReport=AMARETTO_valReport,AMARETTO_trainReport=AMARETTO_trainReport) )
} 




#' AMARETTO_Evaluate_Validation
#'
#' Code to evaluate AMARETTO on a gene expression validation set. Uses output from AMARETTO_Run() 
#' AMARETTOresults AMARETTO output from AMARETTO_Run().
#'MA_Data_TestSet Gene expression matrix from a validation set (that was not used in AMARETTO_Run()).
#' RegulatorData_TestSet Test regulator data with  sameregulators created by AMARETTO_Run().
#'
#' 


AMARETTO_Evaluate_Validation <- function(AMARETTOresults = AMARETTOresults,MA_Data_TestSet =MA_Data_TestSet,RegulatorData_TestSet=RegulatorData_TestSet) {
  
  MA_Data_TestSet<-t(scale(t(MA_Data_TestSet)))
  #RegulatorData_TestSet<-t(scale(t(RegulatorData_TestSet)))
  
  nrSamples = ncol(MA_Data_TestSet)
  RegulatorNames=rownames(RegulatorData_TestSet)
  stats = mat.or.vec(AMARETTOresults$NrModules,10)
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
          outcome_m= MA_Data_TestSet[currentClusterGenes,]
          #PCA 
          PCA_output=prcomp(outcome_m, retx = TRUE, center = TRUE, scale. = FALSE,  tol = NULL, rank. = NULL)
          eigvalues=PCA_output$sdev^2
          perVarianceExplained=100*eigvalues[1]/sum(eigvalues[2:length(eigvalues)])
          
        } else {
          outcome = MA_Data_TestSet[currentClusterGenes,]
          perVarianceExplained=NULL
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
        stats[i,10] = perVarianceExplained
        
      } else {
        stats[i,6] = totalWeightsPercentage
        stats[i,7] = 0
        stats[i,8] = 0
        stats[i,9] = 0
        stats[i,10] = 0
      }
    } else {
      stats[i,6] = 0
      stats[i,7] = 0
      stats[i,8] = 0
      stats[i,9] = 0
      stats[i,10] = 0
    }
  }
  dimnames(stats) <- list(rownames(stats, do.NULL = FALSE, prefix = "Module_"),
                          c("nrPresReg" ,"nrTotalReg", "nrPresGen", "nrTotGen", "percPresGen",
                            "percWeightPresent", "Rsquare", "RsquareAdjusted","MSE", "perVarianceExplained"))
  return(stats)
}

