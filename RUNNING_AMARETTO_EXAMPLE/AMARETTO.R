#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
##  AMARETTO: Regulatory network inference and driver gene evaluation using
##            integrative multi-omics analysis and penalized regression
##  For details on the implementation visit
##  https://github.com/gevaertlab/AMARETTO
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
##########################################################
###################### Example R Script ##################
##########################################################
#-----------------------------------------------------------------------------------------
# 1. Installing AMARETTO and loading the package:
#-----------------------------------------------------------------------------------------
install.packages("devtools")
library(devtools)
devtools::install_github("gevaertlab/AMARETTO",ref="master")
library(AMARETTO)

#-----------------------------------------------------------------------------------------
# 2. Dowloading TCGA input data for analysis
#-----------------------------------------------------------------------------------------
#TargetDirectory <- "./Downloads/" # path to data download directory
#CancerSite <- "READ"
#DataSetDirectories <- AMARETTO_Download(CancerSite,TargetDirectory)

#-----------------------------------------------------------------------------------------
# 3. Preprocessing the downloaded TCGA data
#-----------------------------------------------------------------------------------------
#ProcessedData <- AMARETTO_Preprocess(CancerSite,DataSetDirectories)

#data(MethylStates) #MethylMix preprocessed data for CancerSite
#met <- MethylStates[CancerSite] #MethylMix preprocessed data for CancerSite

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Loading preprocessed TCGA-LAML example dataset:
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
data("ProcessedDataLIHC")

#-----------------------------------------------------------------------------------------
# 4. Running AMARETTO
#-----------------------------------------------------------------------------------------
AMARETTOinit <- AMARETTO_Initialize(MA_matrix = ProcessedDataLIHC$MA_matrix,
                                    CNV_matrix = ProcessedDataLIHC$CNV_matrix,
                                    MET_matrix = ProcessedDataLIHC$MET_matrix,
                                    NrModules = 20, VarPercentage = 50)


AMARETTOresults <- AMARETTO_Run(AMARETTOinit)

AMARETTOtestReport <- AMARETTO_EvaluateTestSet(AMARETTOresults,AMARETTOinit$MA_matrix_Var,AMARETTOinit$RegulatorData)

#-----------------------------------------------------------------------------------------
# 4. Visualize AMARETTO modules
#-----------------------------------------------------------------------------------------
ModuleNr <- 1 #define the module number to visualize

AMARETTO_VisualizeModule(AMARETTOinit,AMARETTOresults,ProcessedDataLIHC$CNV,ProcessedDataLIHC$MET,ModuleNr)

#-----------------------------------------------------------------------------------------
# 5. Get HTML report for AMARETTO modules
#-----------------------------------------------------------------------------------------

AMARETTO_HTMLreport(AMARETTOinit,AMARETTOresults,CNV_matrix=ProcessedDataLIHC$CNV,MET_matrix = ProcessedDataLIHC$MET,VarPercentage=10,hyper_geo_test_bool=FALSE,output_address='./')
