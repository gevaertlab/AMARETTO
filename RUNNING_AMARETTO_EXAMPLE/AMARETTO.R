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
data(ProcessedDataLAML)

#-----------------------------------------------------------------------------------------
# 4. Running AMARETTO 
#-----------------------------------------------------------------------------------------
AMARETTOinit <- AMARETTO_Initialize(MA_matrix = ProcessedDataLAML$MA_TCGA,
                                    CNV_matrix = ProcessedDataLAML$CNV_TCGA,
                                    MET_matrix = ProcessedDataLAML$MET_TCGA,
                                    NrModules = 10, VarPercentage = 5)


AMARETTOresults <- AMARETTO_Run(AMARETTOinit)

AMARETTOtestReport <- AMARETTO_EvaluateTestSet(AMARETTOresults,AMARETTOinit$MA_TCGA_Var,AMARETTOinit$RegulatorData)

#-----------------------------------------------------------------------------------------
# 4. Visualize AMARETTO modules
#-----------------------------------------------------------------------------------------
ModuleNr <- 1 #define the module number to visualize

AMARETTO_VisualizeModule(AMARETTOinit,AMARETTOresults,ProcessedDataLAML$CNV,ProcessedDataLAML$MET,ModuleNr)
