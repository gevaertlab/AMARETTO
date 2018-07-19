\name{BatchData}
\alias{BatchData}
\docType{data}
\title{
	Batch (site collection) data for TCGA cancers.
}
\description{
	This dataset provides the batch IDs for each sample in 12 TCGA cancers. This 
	information is used internally in the TCGA data processing to perform Combat 
	batch correction. 
}
\usage{BatchData}
\format{
	A data frame with 3 variables for each sample in all 12 TCGA cancers: 
	ArrayName, Sample Name (TCGA sample ID), and Batch.
}
\value{
	No value returned, as this is a dataset.
}
\details{
	The data was provided by the main TCGA website and is publically available, 
	but in standard text document formats, in the Aliquot Report.  Processed 
	gdac/Firehose data does not correct for the potential batch effects. 
	AMARETTO first checks to see if a batch effect exists, and if it does, 
	performs Combat batch correction.  The type of array is also included as a 
	covariate in batch correction.
}
\source{
	Aliquot report: 
	https://wiki.nci.nih.gov/display/TCGA/Experiment+Aliquot+Report
	Another potentially helpful file that contains information in addition to 
	standard TCGA data downloads-AliquotBreakdown Report: 
	https://wiki.nci.nih.gov/display/TCGA/Aliquot+ID+Breakdown+Report
}

%\examples{
%data(BatchData)
%}
\keyword{datasets}