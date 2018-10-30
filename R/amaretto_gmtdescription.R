#' GMT Gene Set Description
#' Function to visualize the gene modules
#' @param filename GMT file to read descriptions from

GeneSetDescription<-function(filename){
  gmtLines<-strsplit(readLines(filename),"\t")
  gmtLines_description <- lapply(gmtLines, function(x) {
    c(x[[1]],x[[2]],length(x)-2)
  })
  gmtLines_description<-data.frame(matrix(unlist(gmtLines_description),byrow=T,ncol=3),stringsAsFactors=FALSE)
  rownames(gmtLines_description)<-NULL
  colnames(gmtLines_description)<-c("GeneSet","Description","NumberGenes")
  gmtLines_description$NumberGenes<-as.numeric(gmtLines_description$NumberGenes)
  return(gmtLines_description)
}