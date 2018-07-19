\name{AMARETTO_Run}
\alias{AMARETTO_Run}

\title{
Run AMARETTO
}
\description{
  Function to run AMARETTO, a statistical algorithm to identify cancer drivers by integrating a variety of omics data from cancer and normal tissue.
}
\usage{
  AMARETTO_Run(AMARETTOinit)
}
\arguments{
  \item{AMARETTOinit}{
  List output from AMARETTO_Initialize().
}
}
\value{
  A list with the following elements:
	\item{NrModules}{Number of final gene clusters (modules).}
	\item{RegulatoryPrograms}{Beta values from regulatory programs.}
	\item{AllRegulators}{Names of regulator genes.}
	\item{AllGenes}{Names of passenger genes.}
	\item{ModuleMembership}{Gene module (cluster) assignments for genes.}	
	\item{AutoRegulationReport}{AutoRegulation report for each module, which 
	includes information such as the R squared value to estimate how well the 
	driver genes control the passenger genes' expresion levels.}
	\item{ModuleData}{Module data based on the average expression of the genes in the module.}
	\item{RegulatoryProgramData}{Predicted module data, based on AMARETTO learning program.}
}
\references{
	Gevaert O, Villalobos V, Sikic BI, Plevritis SK. Identification of ovarian 
	cancer driver genes by using module network integration of multi-omics data. 
	Interface Focus. 5(2) (2013)}
\author{
	Olivier Gevaert <ogevaert@stanford.edu>
}

\seealso{
	\code{\link{AMARETTO_Initialize}}, \code{\link{AMARETTO_EvaluateTestSet}}
}

\examples{
########################################################
# Loading the data
########################################################
data(ProcessedDataLAML)
NrModules <- 5
VarPercentage <- 2
AMARETTOinit<-AMARETTO_Initialize(MA_matrix=ProcessedDataLAML$MA_TCGA,
	CNV_matrix=ProcessedDataLAML$CNV_TCGA,MET_matrix=ProcessedDataLAML$MET_TCGA, 
	NrModules=NrModules,VarPercentage=VarPercentage)
  
########################################################
# Running AMARETTO
########################################################
AMARETTOresults<-AMARETTO_Run(AMARETTOinit)
}