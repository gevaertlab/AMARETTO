\name{Save_CancerSite}
\alias{Save_CancerSite}

\title{
	Save preprocessed cancer data
}
\description{
	Saves all three final preprocessed gene expression, CNV, and methylation 
	datasets.
}
\usage{
	Save_CancerSite(CancerSite, TargetDirectory, DataSetDirectories, 
	ProcessedData)
}

\arguments{
  \item{CancerSite}{
		Character string denoting specific TCGA disease.
}
  \item{TargetDirectory}{
		Character string specifying directory in which to save processed data.
}
  \item{DataSetDirectories}{
		List of CNV and gene expression dataset directories; output from 
		Download_CancerSite.
}
  \item{ProcessedData}{
		List of processed data as outputted from Preprocess_CancerSite
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	\item{SaveFile}{Character string of the full path to the final saved RData 
		object.}
}
\seealso{
	\code{\link{Download_CancerSite}},\code{\link{Preprocess_CancerSite}}
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
			
########################################################
# save data (if desired)
########################################################     
#Save_CancerSite(CancerSite,TargetDirectory,DataSetDirectories,ProcessedData)
}