\name{AMARETTO_VisualizeModule}
\alias{AMARETTO_VisualizeModule}
\title{
Visualize AMARETTO module}
\description{
Plot a heatmap of the regulator and passenger gene expression levels for a certain AMARETTO module. Methylation and CNV profiles of driver genes can be included as well as sample annotations.
}
\usage{
AMARETTO_VisualizeModule(AMARETTOinit, AMARETTOresults, CNV_matrix, MET_matrix, ModuleNr, SAMPLE_annotation, ID, order_samples)
}

\arguments{
  \item{AMARETTOinit}{
		Output from AMARETTO_Initialize().
}
  \item{AMARETTOresults}{
		Output from AMARETTO_Run().
 }
  \item{CNV_matrix}{
  CNV matrix, with genes in rows and samples in columns. When absent, CNV profiles are not added.
}
  \item{MET_matrix}{
  Methylation matrix, with genes in rows and samples in columns. When absent, methylation profiles are not added.
}
  \item{ModuleNr}{
		Module number. Selecting a module with a high R squared is a good rule 
		of thumb, as this represents a module that may be tightly regulated.
}
 \item{SAMPLE_annotation}{
  A matrix or dataframe with sample annotations (samples as rows) is optional to add.
}
 \item{ID}{
 The name of the column with the IDs matching those in the data matrixes.
 }
 \item{order_samples}{
 Samples can be ordered on mean expression (standard) or clustering (use ="clust") of the module.
 }
}
\value{
	Does not return output, but plots heatmap to screen.
}

\author{
	Olivier Gevaert <ogevaert@stanford.edu>
}

\seealso{
	\code{\link{AMARETTO_Initialize}}, \code{\link{AMARETTO_Run}}
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
# Visualize a certain module
########################################################   
ModuleNr=2
AMARETTO_VisualizeModule(AMARETTOinit,AMARETTOresults=AMARETTOresults,
 CNV_matrix=ProcessedDataLAML$CNV_TCGA,MET_matrix=ProcessedDataLAML$MET_TCGA,
 ModuleNr=ModuleNr) 
}