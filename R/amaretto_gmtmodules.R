#' AMARETTO_GmtFunctions
#' Output a GMT file for the modules, read gmt, read description from gmt
#' @param AMARETTOinit List output from AMARETTO_Initialize().
#' @param AMARETTOresults List output from AMARETTO_Run().
#' @param filename GMT file to read descriptions from


GmtFromModules <- function(AMARETTOinit,AMARETTOresults){
  
  ModuleMembership<-rownames_to_column(as.data.frame(AMARETTOresults$ModuleMembership),"GeneNames")
  NrModules<-AMARETTOresults$NrModules
  ModuleMembership<-ModuleMembership %>% arrange(GeneNames)
  
  ModuleMembers_list<-split(ModuleMembership$GeneNames,ModuleMembership$ModuleNr)
  names(ModuleMembers_list)<-paste0("Module_",names(ModuleMembers_list))
  
  gmt_file="./Modules_targets_only.gmt"
  write.table(sapply(names(ModuleMembers_list),function(x) paste(x,paste(ModuleMembers_list[[x]],collapse="\t"),sep="\t")),gmt_file,quote = FALSE,row.names = TRUE,col.names = FALSE,sep='\t')
}

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

readGMT<-function(filename){
  gmtLines<-strsplit(readLines(filename),"\t")
  gmtLines_genes <- lapply(gmtLines, tail, -2)
  names(gmtLines_genes) <- sapply(gmtLines, head, 1)
  return(gmtLines_genes)  
}