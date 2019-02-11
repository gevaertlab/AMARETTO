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
install.packages("BiocManager")
BiocManager::install("gevaertlab/AMARETTO")
library(AMARETTO)

resdir <- "AMARETTO_Results";if(!file.exists(resdir))	dir.create(resdir)
setwd(resdir)
#-----------------------------------------------------------------------------------------
# 2. Dowloading TCGA input data for analysis
#-----------------------------------------------------------------------------------------
#TargetDirectory <- file.path(getwd(),"Downloads/") #Absolute path to data download directory
#CancerSite <- "LIHC"
#DataSetDirectories <- AMARETTO_Download(CancerSite,TargetDirectory)

#-----------------------------------------------------------------------------------------
# 3. Preprocessing the downloaded TCGA data
#-----------------------------------------------------------------------------------------
#ProcessedData <- AMARETTO_Preprocess(CancerSite,DataSetDirectories)

#load(data/MethylStates_TCGA.rda) #MethylMix preprocessed data for CancerSite
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

AMARETTOtestReport <- AMARETTO_EvaluateTestSet(AMARETTOresults = AMARETTOresults,
                                               MA_Data_TestSet = AMARETTOinit$MA_matrix_Var,
                                               RegulatorData_TestSet = AMARETTOinit$RegulatorData)

#-----------------------------------------------------------------------------------------
# 4. Visualize AMARETTO modules
#-----------------------------------------------------------------------------------------
ModuleNr <- 1 #define the module number to visualize

AMARETTO_VisualizeModule(AMARETTOinit = AMARETTOinit,AMARETTOresults = AMARETTOresults,
                         CNV_matrix = ProcessedDataLIHC$CNV,MET_matrix = ProcessedDataLIHC$MET,
                         ModuleNr = ModuleNr)

#-----------------------------------------------------------------------------------------
# 5. Get HTML report for AMARETTO modules
#-----------------------------------------------------------------------------------------

AMARETTO_HTMLreport(AMARETTOinit,AMARETTOresults,CNV_matrix=ProcessedDataLIHC$CNV,MET_matrix = ProcessedDataLIHC$MET,VarPercentage=10,hyper_geo_test_bool=FALSE,output_address='./')
