\name{Download_CancerSite}
\alias{Download_CancerSite}

\title{
Download TCGA cancer gene expression, CNV and methylation data.
}
\description{
	Function to download the gene expression, and CNV (GISTIC) data required for 
	AMARETTO for one TCGA cancer (CancerSite.) Requires writing of files to a 
	TargetDirectory, of which the user must have write access.
}

\arguments{
	\item{CancerSite}{
		A string vector specifying the CancerSite, i.e. TCGA cancer. For 
		example, "BRCA" is breast cancer. Use the function get_firehoseData() to 
		display a list of all possible cancers string names.  Names are case-
		sensitive (use caps.)
	}
	\item{TargetDirectory}{
		String vector containing directory to write files to.  User must have 
		write access to this directory.
	}
	\item{downloadData}{ 
		Boolean setting if the data has to be downloaded. Setting 
		to false can be used when re-running a previously downloaded data set 
		and only using the returned data set directories. 
	}
}

\value{
	Returns a list with the two data types:
	\item{MAdirectory}{Directory location of downloaded gene expression data.}
	\item{CNVdirectory}{Directory location of downloaded CNV (GISTIC) data.}
}
\details{
	Uses the function get_firehoseData to download to directly download 
	processed TCGA data from GDAC.
}
\author{
	Olivier Gevaert <ogevaert@stanford.edu>
}
\note{
	Methylation data is also required to run AMARETTO, but this data is taken 
	directly from the AMARETTO package data object MethylStates, and 
	is loaded in Preprocess_CancerSite and processed alongside the expression 
	and CNV data.
}

\seealso{
	\code{\link{get_firehoseData}}, \code{\link{Preprocess_CancerSite}}
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

}