\name{MethylStates}
\alias{MethylStates}
\docType{data}
\title{
	Methylation states list data object
}
\description{
	Methylation data used in the AMARETTO algorithm for 25 TCGA cancers.
}
\usage{MethylStates}
	\format{
	A list, named MethylStates, with the methylation states for each of the 25 
	TCGA cancers. The columns are the methylation features (clusters of 
	methylation probes) that were found to be significantly hyper or hypo-
	methylated in cancer vs. normal tissue for a certain TCGA disease.  The 
	method used to derive these features and significance is MethylMix 
	(available via the MethylMix package in Biocondctor.)
}
\value{
	No value returned, as this is a dataset.
}
\details{
	The method used to derive these features and significance is MethylMix 
	(available via the MethylMix package in Biocondctor.) Actual downloading and 
	processing of raw methylation data is computationally intensive, and thus we 
	provide 	this data directly to the user. This data is used in the 
	AMARETTO algorithm to discover potential cancer drivers.
}
\source{
	Methylation data was processed and analyzed by the Gevaert Lab at Stanford 
	University.
}
\references{
	Gevaert O. MethylMix: an R package for identifying DNA methylation-driven 
	genes. Bioinformatics. (Jan 2015)
}

\keyword{datasets}