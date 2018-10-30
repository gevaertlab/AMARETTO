#' writegct
#' write gct files from data
#' @param data_in A dataframe that needs to be converted to a gct file
#' @param file_address File adress of the gct file
#' @import tidyverse

write_gct<-function(data_in,file_address){
  header_gct<-paste0('#1.2\n',nrow(data_in),'\t',ncol(data_in))
  data_in<-rownames_to_column(as.data.frame(data_in),"Name") %>% dplyr::mutate(Description=Name) %>% dplyr::select(Name,Description,everything())
  write(header_gct,file=file_address,append = FALSE)
  readr::write_tsv(data_in,file_address,append = TRUE,col_names = TRUE)
}