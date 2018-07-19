\name{AMARETTO-package}
\alias{AMARETTO}
\docType{package}
\title{
AMARETTO
}
\description{
	Functions to analyze gene expression, CNV and methylation data via AMARETTO, 
	a gene module network method. Package provides tools to download and analyze 
	publically available TCGA cancer datasets.
	}
\details{
\tabular{ll}{
Package: \tab AMARETTO \cr
Type: \tab Package \cr
Version: \tab 0.99.1 \cr
Date: \tab 2016-06-06\cr
License: \tab GPL (>=2) \cr
}
}
\author{
Magali Champion, Stanford University \cr
Katie Planey, Stanford University \cr
Olivier Gevaert, Stanford University \cr
}
\references{
	Gevaert O, Villalobos V, Sikic BI, Plevritis SK. Identification of ovarian 
	cancer driver genes by using module network integration of multi-omics data. 
	Interface Focus. 5(2) (2013)
	}
\examples{
########################################################
# Loading the data
########################################################
data(ProcessedDataLAML)
NrModules <- 5
VarPercentage <- 2
dim(ProcessedDataLAML$MA_TCGA)
dim(ProcessedDataLAML$CNV_TCGA)
dim(ProcessedDataLAML$MET_TCGA)
AMARETTOinit<-AMARETTO_Initialize(MA_matrix=ProcessedDataLAML$MA_TCGA,
	CNV_matrix=ProcessedDataLAML$CNV_TCGA,MET_matrix=ProcessedDataLAML$MET_TCGA, 
	NrModules=NrModules,VarPercentage=VarPercentage)
  
########################################################
# Running AMARETTO
########################################################
AMARETTOresults<-AMARETTO_Run(AMARETTOinit)

########################################################
# Evaluate AMARETTO Results
########################################################    
AMARETTOtestReport<-AMARETTO_EvaluateTestSet(AMARETTOresults,
	AMARETTOinit$MA_matrix_Var,AMARETTOinit$RegulatorData)

########################################################
# Visualize a certain module
########################################################   
ModuleNr=2
AMARETTO_VisualizeModule(AMARETTOinit,AMARETTOresults=AMARETTOresults,
 CNV_matrix=ProcessedDataLAML$CNV_TCGA,MET_matrix=ProcessedDataLAML$MET_TCGA,
 ModuleNr  =ModuleNr) 
}