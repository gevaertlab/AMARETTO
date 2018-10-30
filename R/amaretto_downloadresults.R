#' AMARETTO_DownloadData
#' Retrieve a download of all the data linked with the run (including heatmaps)
#' @param AMARETTOinit AMARETTO initialize output
#' @param AMARETTOresults AMARETTO results output
#' @param data_address Directory to save data folder
#' @param Heatmaps Output heatmaps as pdf
#' @param data_in A dataframe that needs to be converted to a gct file
#' @param file_address File adress of the gct file
#' @import doParallel
#' @import tidyverse
#' @return
#' @export
#'
#' @examples

AMARETTO_DownloadResults <-function(AMARETTOinit,AMARETTOresults,data_address,Heatmaps=TRUE){
  
  cluster <- makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
  registerDoParallel(cluster,cores=NrCores)
  
  if (!dir.exists(data_address)){
    stop("Output directory is not existing.")
  }
  if (!endsWith(data_address,"/")){
    output_address=paste0(data_address,"/")
  }
  
  output_dir<-paste0("AMARETTOresults_",gsub("-|:","",gsub(" ","_",Sys.time())))
  dir.create(paste0(data_address,output_dir))
  
  if(Heatmaps==TRUE){
    foreach (ModuleNr = 1:NrModules, .packages = c('AMARETTO')) %dopar% {
      pdf(file=paste(data_address,output_dir,"/Module_",as.character(ModuleNr),".pdf",sep=""))
      AMARETTO_VisualizeModule(AMARETTOinit, AMARETTOresults=AMARETTOresults, CNV_matrix, MET_matrix, ModuleNr=ModuleNr) 
      dev.off()
    }
  }
  
  stopCluster(cluster)
  
  save(AMARETTOresults, file=paste0(data_address,output_dir,"/amarettoResults.RData"))
  write_gct(AMARETTOresults$ModuleData,paste0(data_address,output_dir,'/ModuleData_amaretto.gct'))
  write_gct(AMARETTOresults$ModuleMembership,paste0(data_address,output_dir,'/ModuleMembership_amaretto.gct'))
  write_gct(AMARETTOresults$RegulatoryProgramData,paste0(data_address,output_dir,'/RegulatoryProgramData_amaretto.gct'))
  write_gct(AMARETTOresults$RegulatoryPrograms,paste0(data_address,output_dir,'/RegulatoryPrograms_amaretto.gct'))
  readr::write_tsv(as.data.frame(AMARETTOresults$AllGenes),paste0(data_address,output_dir,'/AllGenes_amaretto.tsv'))
  readr::write_tsv(as.data.frame(AMARETTOresults$AllRegulators),paste0(data_address,output_dir,'/AllRegulators_amaretto.tsv'))
  readr::write_tsv(as.data.frame(AMARETTOresults$NrModules),paste0(data_address,output_dir,'/NrModules_amaretto.tsv'))
  
  zip(zipfile = paste0(data_address,output_dir),files=paste0(data_address,output_dir))
}

write_gct<-function(data_in,file_address){
  header_gct<-paste0('#1.2\n',nrow(data_in),'\t',ncol(data_in))
  data_in<-rownames_to_column(as.data.frame(data_in),"Name") %>% dplyr::mutate(Description=Name) %>% dplyr::select(Name,Description,everything())
  write(header_gct,file=file_address,append = FALSE)
  readr::write_tsv(data_in,file_address,append = TRUE,col_names = TRUE)
}