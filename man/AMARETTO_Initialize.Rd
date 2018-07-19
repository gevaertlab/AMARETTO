\name{AMARETTO_Initialize}
\alias{AMARETTO_Initialize}

\title{
AMARETTO_Initialize
}
\description{
  Code used to initialize the seed clusters for an AMARETTO run. Requires processed gene expression (rna-seq or microarray), CNV (usually from a GISTIC run), and methylation (from MethylMix, provided in this package) data. Uses the function CreateRegulatorData() and results are fed into the function AMARETTO_Run().
}
\usage{
AMARETTO_Initialize(MA_matrix, CNV_matrix, MET_matrix, NrModules, VarPercentage, PvalueThreshold = 0.001, RsquareThreshold = 0.1, pmax = 10, NrCores = 1, OneRunStop = 0)
}

\arguments{
  \item{MA_matrix}{
  Expression matrix, with genes in rows and samples in columns.
}
  \item{CNV_matrix}{
  CNV matrix, with genes in rows and samples in columns.
}
  \item{MET_matrix}{
  Methylation matrix, with genes in rows and samples in columns.
}
  \item{NrModules}{
  How many gene co-expression modules should AMARETTO search for? Usually around 100 is acceptable, given the large number of possible driver-passenger gene combinations.
}
  \item{VarPercentage}{
  Minimum percentage by variance for filtering of genes; for example, 75\% would indicate that the CreateRegulatorData() function only analyses genes that have a variance above the 75th percentile across all samples.
}
  \item{PvalueThreshold}{
Threshold used to find relevant driver genes with CNV alterations: maximal p-value.
}
  \item{RsquareThreshold}{
Threshold used to find relevant driver genes with CNV alterations: minimal R-square value between CNV and gene expression data.
}
  \item{pmax}{
  "pmax" variable for glmnet function from glmnet package; the maximum number of variables aver to be nonzero. Should not be changed by user unless she/he fully understands the AMARETTO algorithm and how its parameters choices affect model output.
}
  \item{NrCores}{
		A numeric variable indicating the number of computer/server cores to use 
		for paralellelization. Default is 1, i.e. no parallelization. Please 
		check your computer or server's computing capacities before increasing 
		this number.  Parallelization is done via the RParallel package. Mac 
		vs. Windows environments may behave differently when using 
		parallelization.
}
  \item{OneRunStop}{
  Binary variable. Must be set to 1 if the user wants to run one iteration only (could be useful for large comparisons).
}
}

\value{
	A list with the following elements:
		\item{MA_matrix_Var}{Gene expression data matrix filtered by VarPercentage.}
  	\item{RegulatorData}{RegulatorData output from CreateRegulatorData().}
    \item{RegulatorAlterations}{A list containing information about the copy number and methylation status of each gene: MET, matrix with the percentage of patients with hyper/hypo or no methylation for each gene; CNV, matrix with the percentage of patients with amplification/deletion or normal CNV status for each gene; Summary, a matrix with the status of each driver gene (CNV alterations or methylation aberrations).}
	  \item{ModuleMembership}{Initial "seed" gene expression cluster assignments 
		(clustered by genes, not patients.) There will be NrModules gene 
		clusters.}
	  \item{Parameters}{A list containing several of the inputs to the function, 
		to feed into AMARETTO_Run(): AutoRegulation,OneRunStop,Lambda2,Mode,
		pmax,alpha)}
	  \item{NrCores}{Number of parallelization cores chosen.}
}
\references{
	Gevaert O, Villalobos V, Sikic BI, Plevritis SK. Identification of ovarian 
	cancer driver genes by using module network integration of multi-omics data. 
	Interface Focus. 5(2) (2013)}
\author{
	Olivier Gevaert <ogevaert@stanford.edu>
}
\note{
	\code{\link{Preprocess_CancerSite}},\code{\link{AMARETTO_Run}},
	\code{\link{CreateRegulatorData}}
	}


\examples{
########################################################
# Loading the data
########################################################
data(ProcessedDataLAML)
NrModules <- 100
VarPercentage <- 25
dim(ProcessedDataLAML$MA_TCGA)
dim(ProcessedDataLAML$CNV_TCGA)
dim(ProcessedDataLAML$MET_TCGA)
AMARETTOinit <- AMARETTO_Initialize(MA_matrix=ProcessedDataLAML$MA_TCGA,
	CNV_matrix=ProcessedDataLAML$CNV_TCGA,MET_matrix=ProcessedDataLAML$MET_TCGA, 
	NrModules=NrModules,VarPercentage=VarPercentage)
}