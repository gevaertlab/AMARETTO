\name{CreateRegulatorData}
\alias{CreateRegulatorData}

\title{
Determine potential regulator genes
}
\description{
	Function to identify which genes' CNV and/or methylation statuses 
	significantly predict expression of that gene. Input is output from 
	Preprocess_CancerSite(). This function is used by AMARETTO_Initialize().
	}
\usage{
CreateRegulatorData(MA_matrix, CNV_matrix, MET_matrix, PvalueThreshold = 0.001, RsquareThreshold = 0.1)
}

\arguments{
  \item{MA_matrix}{
		Processed gene expression data matrix ($MA_matrix element from 
		Preprocess_CancerSite() list output).
}
  \item{CNV_matrix}{
		Processed CNV matrix ($CNV_matrix element from Preprocess_CancerSite() 
		list output).
}
  \item{MET_matrix}{
		Processed methylation matrix ($MET_matrix element from 
		Preprocess_CancerSite() list output).
		}
  \item{PvalueThreshold}{
Threshold used to find relevant driver genes with CNV alterations: maximal p-value.
}
  \item{RsquareThreshold}{
Threshold used to find relevant driver genes with CNV alterations: minimal R-square value between CNV and gene expression data.
}
}
\details{
	This function identifies genes whose methylation status significantly 
	predicts negative expression levels and/or whose CNV status significantly 
	predicts positive expression levels of the same gene. These are potential 
	regulator genes.
}
\value{
	A list with the following elements:
	\item{RegulatorData}{Matrix of gene expression values of hypothesized 
	regulator genes.}
    \item{RegulatorAlterations}{A list containing information about the copy number and methylation status of each gene: MET, matrix with the percentage of patients with hyper/hypo or no methylation for each gene; CNV, matrix with the percentage of patients with amplification/deletion or normal CNV status for each gene; Summary, a matrix with the status of each driver gene (CNV alterations or methylation aberrations).}
    }
\references{
	Gevaert O, Villalobos V, Sikic BI, Plevritis SK. Identification of ovarian 
	cancer driver genes by using module network integration of multi-omics data. 
	Interface Focus. 5(2) (2013)
}
\author{
	Olivier Gevaert <ogevaert@stanford.edu>
}
\note{
	While there remain biological exceptions, generally, methylation values at a 
	promotor site induce a decrease in gene expression, while CNV amplifications 
	induce an increase in gene expression. This function models these 
	assumptions.
}

\seealso{
	\code{\link{Preprocess_CancerSite}},\code{\link{AMARETTO_Initialize}} 
}
\examples{
########################################################
# Identify potential regulator genes
######################################################## 
data(ProcessedDataLAML)

RegulatorInfo=CreateRegulatorData(MA_matrix=ProcessedDataLAML$MA_TCGA,
	CNV_matrix=ProcessedDataLAML$CNV_TCGA,MET_matrix=ProcessedDataLAML$MET_TCGA)
}