\name{AMARETTO_EvaluateTestSet}
\alias{AMARETTO_EvaluateTestSet}
\title{
AMARETTO_EvaluateTestSet}
\description{
Code to evaluate AMARETTO on a new gene expression test set. Uses output from AMARETTO_Run() and CreateRegulatorData().
}
\usage{
AMARETTO_EvaluateTestSet(AMARETTOresults, MA_Data_TestSet, RegulatorData_TestSet)
}
\arguments{
  \item{AMARETTOresults}{
AMARETTO output from AMARETTO_Run().
}
  \item{MA_Data_TestSet}{
		Gene expression matrix from a test set (that was not used in 
		AMARETTO_Run()).
}
  \item{RegulatorData_TestSet}{
		Test regulator data from CreateRegulatorData().
}
}
\details{
	Calculates the R squared values between the driver and passenger genes on 
	the new gene expression matrix.
}
\value{
	stats, a matrix with the following column variables: "nrPresReg" (number of 
	present regulators),"nrTotalReg" (number of total regulators), "nrPresGen" 
	(number of present genes), "nrTotGen" (number of total genes), "percPresGen" 
	(percent of present genes), "percWeightPresent" (percent of regulation 
	accounted for by model weight), "Rsquare" (R squared metric), 
	"RsquareAdjusted" (adjusted R squared metric).
}

\author{
	Olivier Gevaert <ogevaert@stanford.edu>
}
\seealso{
	\code{\link{CreateRegulatorData}}, \code{\link{AMARETTO_Run}}
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
# In the absence of a true test set, test this on the training set
AMARETTOtestReport<-AMARETTO_EvaluateTestSet(AMARETTOresults,
	AMARETTOinit$MA_matrix_Var,AMARETTOinit$RegulatorData)
}