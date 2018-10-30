#' readGMT
#' function to read gmt files
#' @param filename gmt filename

readGMT<-function(filename){
  gmtLines<-strsplit(readLines(filename),"\t")
  gmtLines_genes <- lapply(gmtLines, tail, -2)
  names(gmtLines_genes) <- sapply(gmtLines, head, 1)
  return(gmtLines_genes)  
}