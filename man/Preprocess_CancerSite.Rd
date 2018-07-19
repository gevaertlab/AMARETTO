\name{Preprocess_CancerSite}
\alias{Preprocess_CancerSite}

\title{
	Preprocess CancerSite expression, CNV and methylation data
}
\description{
	A function to process gene expression, CNV and methylation data (from called 
	Dwonload_CancerSite) for a specific TCGA cancer into a format readily 
	accepted by the AMARETTO algorithm.
}
\usage{
	Preprocess_CancerSite(CancerSite, DataSetDirectories)
}

\arguments{
  \item{CancerSite}{
		TCGA cancer acronym (function analyzes one cancer at a time).
}
  \item{DataSetDirectories}{
		Director(ies) where data is stored; outputted by TCGAdownload(). Most 
		cancers will have just one data directory, but sometimes, say for GBM 
		and microarrays, two platforms were used.  The user does not need to 
		specifically know this, as the Download_CancerSite function returns the 
		correct number and names of the directories needed.
}
}
\details{
	Calls functions Preprocess_MAdata, TCGA_Load_GISTICdata, and 
	FindTranscriptionallyPredictive_CNV
	}
\value{
	A list of gene expression, CNV and methlation matrices to be fed into 
	AMARETTO.
	\item{MA_TCGA}{Processed gene expression data matrix}
	\item{CNV_TCGA}{Processed CNV data matrix}
	\item{MET_TCGA}{Processed methylation data matrix}
}
\note{
	GISTIC is the algorithm used to analyze TCGA CNV data found at 
	http://gdac.broadinstitute.org/
}
\seealso{
	\code{\link{Download_CancerSite}},\code{\link{Save_CancerSite}},
}

\examples{
#Takes a while to download an individual cancer data set. Estimated time 
#is between 10-60 minutes depending on data set size (READ is smallest, 
#BRCA is largest).
#
#To run this set the TargetDirectory to your path and uncomment the 
#Download_CancerSite and Preprocess_CancerSite function calls. READ is a good 
#example to start as it takes about 5 minutes to download.
CancerSite='READ'
TargetDirectory="./"
    		
########################################################
# Download TCGA expression and CNV data
########################################################       
#cat('Loading data for',CancerSite,"\n")
#DataSetDirectories <- Download_CancerSite(CancerSite,TargetDirectory)			

########################################################
# Prepare/process expression, CNV and methylation data
######################################################## 
#ProcessedData <- Preprocess_CancerSite(CancerSite,DataSetDirectories)
			
}