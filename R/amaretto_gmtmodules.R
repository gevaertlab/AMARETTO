#' AMARETTO_GmtFromModules
#' Output a GMT file for the modules
#' @param AMARETTOinit List output from AMARETTO_Initialize().
#' @param AMARETTOresults List output from AMARETTO_Run().

GmtFromModules <- function(AMARETTOinit,AMARETTOresults){
  
  ModuleMembership<-rownames_to_column(as.data.frame(AMARETTOresults$ModuleMembership),"GeneNames")
  NrModules<-AMARETTOresults$NrModules
  ModuleMembership<-ModuleMembership %>% arrange(GeneNames)
  
  ModuleMembers_list<-split(ModuleMembership$GeneNames,ModuleMembership$ModuleNr)
  names(ModuleMembers_list)<-paste0("Module_",names(ModuleMembers_list))
  
  gmt_file="./Modules_targets_only.gmt"
  write.table(sapply(names(ModuleMembers_list),function(x) paste(x,paste(ModuleMembers_list[[x]],collapse="\t"),sep="\t")),gmt_file,quote = FALSE,row.names = TRUE,col.names = FALSE,sep='\t')
}