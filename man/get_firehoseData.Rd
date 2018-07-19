\name{get_firehoseData}
\alias{get_firehoseData}

\title{
	get_firehoseData
}
\description{
	A method to directly download processed TCGA data from GDAC. Traditionally, 
	a user would have to make queries by first downloading the firehose command 
	line tool, but this function allows you to download any data type directly 
	in R. It always searched for the most recent available dataset by date. Used 
	by Download_CancerSite().
}
\usage{
get_firehoseData(downloadData = TRUE, saveDir = "./", TCGA_acronym_uppercase = "LUAD", dataType = "stddata", dataFileTag = "mRNAseq_Preprocess.Level_3", FFPE = FALSE, fileType = "tar.gz", gdacURL = "http://gdac.broadinstitute.org/runs/", untarUngzip = TRUE, printDisease_abbr = FALSE)
}

\arguments{
  \item{downloadData}{
		A boolean. To downloadData, set to TRUE (default.)
}
  \item{saveDir}{
		A boolean. If downloadData=TRUE, where the data will be downloaded. The 
		user must have system permissions to download data here, otherwise the 
		function will throw a permissions error.
		}
  \item{TCGA_acronym_uppercase}{
		The TCGA disease character string, such as "BRCA" for breast cancer. 
		Options can be viewed by running with function first with 
		printDisease_abbr set to TRUE.
		}
  \item{dataType}{
		A character string specifying either "stddata" or "analyses". An example 
		of a subsequent dataFileTag options for stddata is 
		mRNAseq_Preprocess.Level_3. For analyses, examples are 
		CopyNumber_Gistic2.Level_4, 
		Pathway_Paradigm_RNASeq_And_Copy_Number.Level_4,  
		Pathway_Paradigm_mRNA_And_Copy_Number.Level_4 and 
		MutSigNozzleReport2.0.Level_4. Stddata and analyses are separated into 
		two different file paths at http://gdac.broadinstitute.org/runs/, hence 
		the need for the user to understand if she/he is simply downloaded 
		processed data like an expression matrix ("stddata") or output from an 
		analytical pipeline like GISTIC ("analyses").
}
  \item{dataFileTag}{
		A character string specifying what the "grep" or search for among all of 
		the possible datasets on http://gdac.broadinstitute.org/runs/. The more 
		specific, the better, although it is not necessary to include a tag with 
		the specific TCGA acronym. For example, the default dataFileTag = 
		"mRNAseq_Preprocess.Level_3" is appropriate for any TCGA cancer when a 
		user wants rna-seq expression data that is already pre-processed and 
		read to be analyzed.
}
  \item{FFPE}{
		If FFPE = FALSE, any datasets that are FFPE will not be used. The 
		default is FFPE=FALSE because FFPE samples tend to be of lower quality.  
		If FFPE=TRUE, only FFPE datasets will be returned (if they exist).
}
  \item{fileType}{
		There are several data types; the most common a user will want to use is 
		the default, "tar.gz". But one can look at 
		http://gdac.broadinstitute.org/runs/ to see that there are, for example, 
		md5 files also. This is simply a keyword to grep for at the end of 
		files.
}
  \item{gdacURL}{
		Most likely won't be changed by the user, but this is the default gdac 
		URL where the firehose processed data is stored. It's used to build the 
		more detailed URL queries as this function travels through the site 
		pages.
		}
  \item{untarUngzip}{
		Boolean. Should the downloaded file also be untarred/ungzipped? 
		Currently only available for tar.gz files.
		}
  \item{printDisease_abbr}{
		A boolean. If TRUE, the function only returns all possible TCGA 
		acronyms. The user can run this function, see which acronym fits the 
		desired disease for analysis, and then run the function again with with 
		this variable set to FALSE and the correct acronym to download the data. 
		Default is FALSE, as this allows for data downloads.
		}
}

\note{
	This function can be used to download a variety of files. To see which 
	dataFileTag should be used for your desired dataset type, click on one of 
	the data dump dates at http://gdac.broadinstitute.org/runs/ and peruse the 
	different datasets available.
}

\seealso{
	\code{\link{Download_CancerSite}}
}

\examples{
#look up all of the possible TCGA cancers/diseases we can query:
get_firehoseData(printDisease_abbr=TRUE)

#don't download data - just confirm that data for your dataFileTag exists
#and see what the most recent date for this data is
get_firehoseData(downloadData=FALSE,TCGA_acronym_uppercase = "READ",
	dataType="stddata")

	\dontrun{
		\donttest{
#now actually download the data:
#This can take between 10-60 minutes depending on data set size.
dataDownloadDirectory <- get_firehoseData(downloadData=FALSE,
	TCGA_acronym_uppercase = "READ",dataType="stddata")

		}
	}
}