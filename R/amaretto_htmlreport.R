#' AMARETTO_HTMLreport
#'
#' Retrieve an interactive html report, including gene set enrichment analysis if asked for.
#'
#' @param AMARETTOinit AMARETTO initialize output
#' @param AMARETTOresults AMARETTO results output
#' @param ProcessedData List of processed input data
#' @param SAMPLE_annotation SAMPLE annotation will be added to heatmap
#' @param ID ID column of the SAMPLE annotation data frame
#' @param hyper_geo_reference Either GMT file address for genesets or computed GSEA dataframe using HyperGeoEnrichmentTest()
#' @param output_address Output directory for the html files.
#' @param show_row_names if True, sample names will appear in the heatmap
#' @param driverGSEA if TRUE, module drivers will also be included in the hypergeometric test.
#' @param phenotype_association_table Optional, Phenotype Association table.
#'
#' @import dplyr
#' @importFrom doParallel registerDoParallel
#' @importFrom DT datatable formatRound formatSignif  formatStyle styleColorBar styleInterval
#' @importFrom reshape2 melt
#' @importFrom dplyr arrange group_by left_join mutate select summarise rename filter case_when
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom knitr knit_meta
#' @importFrom utils write.table
#' @importFrom tibble rownames_to_column
#' @importFrom stats p.adjust  phyper
#' @importFrom rmarkdown render
#' @return result
#' @export
#' @examples
#'\dontrun{
#' data('ProcessedDataLIHC')
#' AMARETTOinit <- AMARETTO_Initialize(ProcessedData = ProcessedDataLIHC,
#'                                     NrModules = 2, VarPercentage = 50)
#'
#' AMARETTOresults <- AMARETTO_Run(AMARETTOinit)
#'
#' AMARETTO_HTMLreport(AMARETTOinit= AMARETTOinit,AMARETTOresults= AMARETTOresults,
#'                     ProcessedData = ProcessedDataLIHC,
#'                     hyper_geo_test_bool=FALSE,
#'                     output_address='./')
#'}
AMARETTO_HTMLreport <- function(AMARETTOinit,
                                AMARETTOresults,
                                ProcessedData,
                                show_row_names = FALSE,
                                SAMPLE_annotation = NULL,
                                ID = NULL,
                                hyper_geo_reference = NULL,
                                genetic_pert_hyper_geo_reference = NULL,
                                chem_pert_hyper_geo_reference = NULL,
                                output_address = './',
                                driverGSEA = TRUE,
                                phenotype_association_table = NULL){
  
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  CNV_matrix <- ProcessedData[[2]]
  MET_matrix <- ProcessedData[[3]]
  NrModules <- AMARETTOresults$NrModules
  VarPercentage <- AMARETTOinit$Parameters$VarPercentage
  
  #set number of cores and check
  NrCores <- AMARETTOinit$NrCores
  MaxCores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
  
  options('DT.warn.size'=FALSE)
  
  if(MaxCores < NrCores){
  stop(paste0("The number of cores that is asked for (",NrCores,"), is more than what's avalaible. Changes can be made on AMARETTOinit$NrCores."))
  }
  
  #check directory
  if (!dir.exists(output_address)){
    stop("Output directory is not existing.")
  }
  

  report_address <- file.path(output_address)
  dir.create(paste0(report_address, "/AMARETTOhtmls/modules"), recursive = TRUE, showWarnings = FALSE)
  cat("The output folder structure is created.\n")
  
  #==============================================================================================================
  hyper_geo_test_bool<-TRUE
  if(is.null(hyper_geo_reference)){
    hyper_geo_test_bool<-FALSE
  }else if (is.data.frame(hyper_geo_reference)){
    output_hgt<-hyper_geo_reference
  }else if(is.character(hyper_geo_reference)&file.exists(hyper_geo_reference[1])){
    output_hgt <-HyperGeoEnrichmentTest(AMARETTOinit, AMARETTOresults, hyper_geo_reference, driverGSEA, NrCores)
  }else {
    stop("The hyper_geo_reference is not properly provided. It should be either an address to an existing .gmt file or a hyper-geo-test dataframe table\n")
  }
  #======================================
  genetic_pert_hyper_geo_test_bool<-TRUE
  if(is.null(genetic_pert_hyper_geo_reference)){
    genetic_pert_hyper_geo_test_bool<-FALSE
  }else if (is.data.frame(genetic_pert_hyper_geo_reference)){
    genetic_pert_output_hgt<-genetic_pert_hyper_geo_reference
  }else if(is.character(genetic_pert_hyper_geo_reference)&file.exists(genetic_pert_hyper_geo_reference[1])){
    genetic_pert_output_hgt <-HyperGeoEnrichmentTest(AMARETTOinit, AMARETTOresults, genetic_pert_hyper_geo_reference, driverGSEA, NrCores)
  }else {
    stop("The genetic_pert_hyper_geo_reference is not properly provided. It should be either an address to an existing .gmt file or a hyper-geo-test dataframe table\n")
  }
  #======================================
  chem_pert_hyper_geo_test_bool<-TRUE
  if(is.null(genetic_pert_hyper_geo_reference)){
    chem_pert_hyper_geo_test_bool<-FALSE
  }else if (is.data.frame(chem_pert_hyper_geo_reference)){
    chem_pert_output_hgt<-chem_pert_hyper_geo_reference
  }else if(is.character(chem_pert_hyper_geo_reference)&file.exists(chem_pert_hyper_geo_reference[1])){
    chem_pert_output_hgt <-HyperGeoEnrichmentTest(AMARETTOinit, AMARETTOresults, chem_pert_hyper_geo_reference, driverGSEA, NrCores)
  }else{
    stop("The chem_pert_hyper_geo_reference is not properly provided. It should be either an address to an existing .gmt file or a hyper-geo-test dataframe table\n")
  }
  #==============================================================================================================
  
  #Parallelizing
  cluster <- parallel::makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
  doParallel::registerDoParallel(cluster,cores=NrCores)

  full_path <- normalizePath(report_address)
  unlink(paste0(full_path,"/*"))

  ModuleOverviewTable <- NULL
  
  yml_file <- paste0(full_path,"/AMARETTOhtmls/modules/_site.yml")
  file.copy(system.file("templates/module_templates/_site.yml",package="AMARETTO"),yml_file)
  
  ModuleOverviewTable<-foreach (ModuleNr = 1:NrModules, .packages = c('AMARETTO','tidyverse','DT','rmarkdown')) %dopar% {
  #for(ModuleNr in 1:NrModules){
    #get heatmap
    
    print(paste0("ModuleNr = ",ModuleNr))
    heatmap_module <- AMARETTO_VisualizeModule(AMARETTOinit, AMARETTOresults, ProcessedData, show_row_names = show_row_names, SAMPLE_annotation=SAMPLE_annotation, ID=ID, ModuleNr=ModuleNr,printHM = FALSE)
    print("The Heatmap is visualised.")
    
    # create datables that are supplied to the RMarkdown file
    ModuleRegulators <- AMARETTOresults$RegulatoryPrograms[ModuleNr,which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
    print("ModuleRegulators are defined.")
    filename_table <- paste0("regulators_module",ModuleNr)
    buttons_list <- list(list(extend ='csv',filename=filename_table), list(extend ='excel',filename=filename_table), list(extend = 'pdf', pageSize = 'A4', orientation = 'landscape',filename=filename_table),list(extend ='print'), list(extend ='colvis'))
    column_markup <- list(list(width = '200px', className = 'dt-head-center', targets = "_all"), list(className = 'text-left', targets = "_all"))
    dt_regulators <- DT::datatable(tibble::rownames_to_column(as.data.frame(ModuleRegulators),"RegulatorIDs") %>% 
                                     dplyr::rename(Weights="ModuleRegulators") %>% 
                                     dplyr::mutate(Weights=signif(Weights, digits = 3)) %>% 
                                     dplyr::mutate(RegulatorIDs=paste0('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=',RegulatorIDs,'">',RegulatorIDs,'</a>')) %>% 
                                     dplyr::arrange(-Weights),
                                class = 'display', filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE, 
                                options = list(columnDefs = column_markup, pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), 
                                               keys = TRUE, dom = 'Blfrtip', buttons = buttons_list),
                               colnames = c("Driver Gene", "Weight"), escape = 'Weight') %>% 
                      DT::formatStyle('Weights',color = DT::styleInterval(0, c('darkblue', 'darkred')))
    print("Data tabel of regulators is created.")
    
    filename_table <- paste0("targets_module",ModuleNr)
    buttons_list <- list(list(extend ='csv',filename=filename_table), list(extend ='excel',filename=filename_table), list(extend = 'pdf', pageSize = 'A4', orientation = 'landscape',filename=filename_table),list(extend ='print'), list(extend ='colvis'))
    dt_targets <- DT::datatable(as.data.frame(AMARETTOresults$ModuleMembership) %>% 
                                  tibble::rownames_to_column("TargetIDs") %>% 
                                  dplyr::arrange(TargetIDs) %>% 
                                  dplyr::rename(moduleNr=ModuleNr) %>% 
                                  dplyr::filter(moduleNr==ModuleNr) %>% 
                                  dplyr::select(-moduleNr) %>% 
                                  dplyr::mutate(TargetIDs=paste0('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=',TargetIDs,'">',TargetIDs,'</a>')),
                                class = 'display', filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE, 
                                options = list(columnDefs = column_markup, pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), 
                                               keys = TRUE, dom = 'Blfrtip', buttons = buttons_list),
                                colnames = c("Target Gene"),escape = FALSE)
    print("Data tabel of targets is created.")
    #=========================================================================================================
    # create GSEA output table, taking into account the resource of the GMT file (eg. MSIGDB)
    if (hyper_geo_test_bool){
      dt_genesets_list<-create_hgt_datatable(output_hgt=output_hgt, module_table=TRUE, ModuleNr = ModuleNr)
      dt_genesets<-dt_genesets_list$dt_genesets
      ngenesets <- dt_genesets_list$ngenesets

    } else {
      dt_genesets <- "Genesets were not analysed as they were not provided."
      ngenesets<-"NA"
      #ngenesets <- "NA"
    }
    print("Data Table for GSEA results is created.")
    #=========================================================
    # create GSEA output table, taking into account the resource of the GMT file (eg. MSIGDB)
    if (genetic_pert_hyper_geo_test_bool){
      dt_genesets_genetic_pert<-create_hgt_datatable(output_hgt=genetic_pert_output_hgt, module_table=TRUE, ModuleNr = ModuleNr)
      dt_genesets_genetic_pert<-dt_genesets_genetic_pert$dt_genesets
    } else {
      dt_genesets_genetic_pert <- "Genesets were not analysed as they were not provided."
      #ngenesets <- "NA"
    }
    print("Data Table for GSEA results is created.")
    #=========================================================
    # create GSEA output table, taking into account the resource of the GMT file (eg. MSIGDB)
    if (chem_pert_hyper_geo_test_bool){
      dt_genesets_chem_pert<-create_hgt_datatable(output_hgt=chem_pert_output_hgt, module_table=TRUE, ModuleNr = ModuleNr)
      dt_genesets_chem_pert<-dt_genesets_chem_pert$dt_genesets
    } else {
      dt_genesets_chem_pert <- "Genesets were not analysed as they were not provided."
      #ngenesets <- "NA"
    }
    print("Data Table for GSEA results is created.")
    #=========================================================
    #created datatable for phenotype associations
    if (!is.null(phenotype_association_table)){
      filename_table <- paste0("phenotypes_module",ModuleNr)
      buttons_list <- list(list(extend ='csv',filename=filename_table), list(extend ='excel',filename=filename_table), list(extend = 'pdf', pageSize = 'A4', orientation = 'landscape',filename=filename_table),list(extend ='print'), list(extend ='colvis'))
      dt_phenotype_association <- DT::datatable(phenotype_association_table %>% 
                                                                dplyr::filter(ModuleNr==paste0("Module ",!!ModuleNr)) %>% 
                                                                dplyr::mutate(p.value = signif(p.value, digits = 3), q.value = signif(q.value, digits = 3)) %>% 
                                                                dplyr::arrange(q.value) %>%
                                                                dplyr::select(-ModuleNr), class='display', filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE, 
                                                              options = list(pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',buttons = buttons_list),
                                                              colnames=c("Phenotype","Statistics Test","P-value","FDR Q-value","Descriptive Statistics"),escape = FALSE) %>% 
                                                DT::formatSignif(c('p.value','q.value'), 2)
    } else{
      dt_phenotype_association <- "Phenotype association resuls were not provided."
    }
    print("The datatable with phenotype association results is created.")
    
    #copy the template file, needed when parallelized
    modulemd <- paste0(full_path,"/AMARETTOhtmls/modules/module",ModuleNr,".rmd")
    file.copy(system.file("templates/module_templates/TemplateReportModule.Rmd",package="AMARETTO"),modulemd)
    
    print("The copy of the template file is created.")
    # output_format<-system.file("templates/module_templates/TemplateReportModule.Rmd",package="AMARETTO")
    knitr::knit_meta(class=NULL, clean = TRUE)
    rmarkdown::render(modulemd, 
                      output_file = paste0("module",ModuleNr,".html"),
                      params = list(
                      report_address = report_address,
                      ModuleNr = ModuleNr,
                      heatmap_module = heatmap_module,
                      dt_regulators = dt_regulators,
                      dt_targets = dt_targets,
                      dt_phenotype_association = dt_phenotype_association,
                      dt_genesets = dt_genesets,
                      dt_genesets_genetic_pert = dt_genesets_genetic_pert,
                      dt_genesets_chem_pert = dt_genesets_chem_pert), knit_meta=knitr::knit_meta(class=NULL, clean = TRUE),quiet = TRUE)
    print("Rmarkdown created the module html page.")
    
    #remove rmd copy of template
    file.remove(modulemd)
    #file.remove(paste0(full_path,"/AMARETTOhtmls/modules/module",ModuleNr,"_files"))
    print("file removed successfully :) Done!")
    #ModuleOverviewTable<-rbind(ModuleOverviewTable,c(ModuleNr,length(which(AMARETTOresults$ModuleMembership==ModuleNr)),length(ModuleRegulators),ngenesets))
    while (!is.null(dev.list()))  dev.off()
    return(c(ModuleNr, length(which(AMARETTOresults$ModuleMembership==ModuleNr)), length(ModuleRegulators),ngenesets))
    
    # },error=function(e){message(paste("an error occured for Module", ModuleNr))})
  }
  
  suppressWarnings(suppressMessages(file.remove(paste0(full_path,"/AMARETTOhtmls/modules/_site.yml"))))
  file_remove<-suppressWarnings(suppressMessages(file.remove(paste0(full_path,"/AMARETTOhtmls/modules/module",c(1:NrModules),"_files"))))
  parallel::stopCluster(cluster)
  
  cat("All module htmls are created.\n")
  ModuleOverviewTable <- data.frame(matrix(unlist(ModuleOverviewTable), byrow=TRUE, ncol=4), stringsAsFactors=FALSE)
  colnames(ModuleOverviewTable)<-c("ModuleNr","NrTarGenes","NrRegGenes","SignGS")
  
  if (!is.null(CNV_matrix)){
    nCNV = ncol(CNV_matrix)
  } else {nCNV = NA}
  if (!is.null(MET_matrix)){
    nMET = ncol(MET_matrix)
  } else {nMET = NA}
  
  nExp = ncol(AMARETTOresults$RegulatoryProgramData)
  nGenes = length(AMARETTOresults$AllGenes)
  nMod = AMARETTOresults$NrModules
  
  options('DT.warn.size'=FALSE) # avoid showing datatable size-related warnings.
  
  filename_table <- "overview_modules"
  buttons_list <- list(list(extend ='csv',filename=filename_table), list(extend ='excel',filename=filename_table), list(extend = 'pdf', pageSize = 'A4', orientation = 'landscape',filename=filename_table),list(extend ='print'), list(extend ='colvis'))
  
  dt_overview<-DT::datatable(ModuleOverviewTable %>% 
                               dplyr::mutate(ModuleNr=paste0('<a href="./modules/module',ModuleNr,'.html">Module ',ModuleNr,'</a>')), 
                             class = 'display', filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE, colnames =c("Module","# Target Genes", "# Driver Genes", "# Gene Sets"),
                             options = list(pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100, 200), keys = TRUE, dom = 'Blfrtip',buttons = buttons_list,columnDefs = list(list(className = 'dt-head-center', targets = "_all"),list(className = 'text-left', targets = "_all"))),
                             escape = FALSE)
  
  all_targets<-tibble::rownames_to_column(data.frame(AMARETTOresults$ModuleMembership),"Genes") %>% 
    dplyr::rename(Module="ModuleNr") %>%
    dplyr::mutate(value=0) %>% 
    dplyr::mutate(Type="Target") %>% 
    dplyr::select(Genes,Module,value,Type)
  
  all_regulators <- reshape2::melt(tibble::rownames_to_column(as.data.frame(AMARETTOresults$RegulatoryPrograms),"Module"),id.vars = "Module") %>% 
    dplyr::filter(value!=0) %>% dplyr::mutate(Module=sub("Module_","",Module),Type="Driver") %>% 
    dplyr::rename(Genes='variable') %>%
    dplyr::select(Genes,Module,value,Type)
  
  

  all_genes <- rbind(all_targets,all_regulators) %>% 
    dplyr::arrange(Genes) %>% 
    dplyr::mutate(Genes=paste0('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=',Genes,'">',Genes,'</a>')) %>% 
    dplyr::mutate(Module=paste0('<a href="./modules/module',Module,'.html">Module ',Module,'</a>'))
  
  all_genes <- all_genes %>%
    dplyr::mutate(Color=dplyr::case_when(
          is.na(as.numeric(value))~"",
          as.numeric(value)>0~"darkred",
          as.numeric(value)<0~"darkblue",
          TRUE~"darkgreen")) %>%
    dplyr::mutate(Type=paste0('<font color=',Color,'>',Type,'</font>')) %>% 
    dplyr::select(-Color,-value)
  
  filename_table <- "genes_to_modules"
  buttons_list <- list(list(extend ='csv',filename=filename_table), list(extend ='excel',filename=filename_table), list(extend = 'pdf', pageSize = 'A4', orientation = 'landscape',filename=filename_table),list(extend ='print'), list(extend ='colvis'))
  
  dt_genes <- DT::datatable(all_genes, 
                          class = 'display', filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE,colnames =c("Gene","Module","Gene Type"),
                          options = list(deferRender=TRUE,columnDefs = list(list(className = 'dt-head-center', targets = "_all"), list(className = 'text-left', targets = "_all")), pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',buttons = buttons_list),
                          escape = FALSE)
  #=================================================================================
  if (hyper_geo_test_bool){
    dt_genesetsall <- create_hgt_datatable(output_hgt = output_hgt, module_table = FALSE)
  } else {
    dt_genesetsall <- data.frame(Hyper_Geometric_Test="Genesets were not analysed as they were not provided.")
  }
  #=============================
  if (genetic_pert_hyper_geo_test_bool){
    dt_genesetsall_genetic_pert <- create_hgt_datatable(output_hgt = genetic_pert_output_hgt, module_table = FALSE)
  } else {
    dt_genesetsall_genetic_pert <- data.frame(Hyper_Geometric_Test="Genesets were not analysed as they were not provided.")
  }
  #=============================
  if (chem_pert_hyper_geo_test_bool){
    dt_genesetsall_chem_pert <- create_hgt_datatable(output_hgt = chem_pert_output_hgt, module_table = FALSE)
  } else {
    dt_genesetsall_chem_pert <- data.frame(Hyper_Geometric_Test="Genesets were not analysed as they were not provided.")
  }
  #=============================
  #created phenotype table for index page
  if (!is.null(phenotype_association_table)){
    filename_table <- "phenotypes_all_modules"
    buttons_list <- list(list(extend ='csv',filename=filename_table), list(extend ='excel',filename=filename_table), list(extend = 'pdf', pageSize = 'A4', orientation = 'landscape',filename=filename_table),list(extend ='print'), list(extend ='colvis'))
    
    dt_phenotype_association_all <- DT::datatable(phenotype_association_table %>% 
                                                    dplyr::mutate(p.value=signif(p.value, digits = 3), q.value=signif(q.value, digits = 3)) %>% 
                                                    dplyr::mutate(ModuleNr=paste0('<a href="./modules/module',gsub("Module ","",ModuleNr),'.html">',ModuleNr,'</a>'))%>%
                                                    dplyr::arrange(q.value), class='display',filter = 'top', extensions = c('Buttons','KeyTable'),rownames = FALSE,
                                               options = list(pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',buttons = buttons_list, columnDefs = list(list(className = 'dt-head-center', targets = "_all"),list(className = 'text-left', targets = "_all"))),colnames=c("Module","Phenotype","Statistics Test","P-value","FDR Q-value","Descriptive Statistics"),
                                               escape = FALSE) %>% DT::formatSignif(c('p.value','q.value'),2)
  }
  else{
    dt_phenotype_association_all <- data.frame(Phenotype_Association="Phenotype association resuls were not provided.")
  }
  
  #Render index page
  rmarkdown::render(system.file("templates/TemplateIndexPage.Rmd",package="AMARETTO"), output_dir=paste0(full_path,"/AMARETTOhtmls/"),output_file= "index.html", params = list(
    nExp = nExp,
    nCNV = nCNV,
    nMET = nMET,
    nGenes = nGenes,
    VarPercentage = VarPercentage,
    nMod = nMod,
    dt_overview = dt_overview),quiet = TRUE)
  
  rmarkdown::render(system.file("templates/TemplateIndexPage_Overview.Rmd",package="AMARETTO"), output_dir=paste0(full_path,"/AMARETTOhtmls/"),output_file= "index_Overview.html", params = list(
    dt_overview = dt_overview),quiet = TRUE)
  
  rmarkdown::render(system.file("templates/TemplateIndexPage_AllGenes.Rmd",package="AMARETTO"), output_dir=paste0(full_path,"/AMARETTOhtmls/"),output_file= "index_AllGenes.html", params = list(
    dt_genes = dt_genes),quiet = TRUE)
  
  rmarkdown::render(system.file("templates/TemplateIndexPage_GenesetsEnrichment.Rmd",package="AMARETTO"), output_dir=paste0(full_path,"/AMARETTOhtmls/"),output_file= "index_GenesetsEnrichment.html", params = list(
    dt_gensesetsall = dt_gensesetsall),quiet = TRUE)
  
  rmarkdown::render(system.file("templates/TemplateIndexPage_GenesetsEnrichment_gp.Rmd",package="AMARETTO"), output_dir=paste0(full_path,"/AMARETTOhtmls/"),output_file= "index_GenesetsEnrichment_gp.html", params = list(
    dt_genesetsall_genetic_pert = dt_genesetsall_genetic_pert),quiet = TRUE)
  
  rmarkdown::render(system.file("templates/TemplateIndexPage_GenesetsEnrichment_cp.Rmd",package="AMARETTO"), output_dir=paste0(full_path,"/AMARETTOhtmls/"),output_file= "index_GenesetsEnrichment_cp.html", params = list(
    dt_genesetsall_chem_pert = dt_genesetsall_chem_pert),quiet = TRUE)
  
  rmarkdown::render(system.file("templates/TemplateIndexPage_PhenoAssociation.Rmd",package="AMARETTO"), output_dir=paste0(full_path,"/AMARETTOhtmls/"),output_file= "index_PhenoAssociation.html", params = list(
    dt_phenotype_association_all = dt_phenotype_association_all),quiet = TRUE)
  
  
  dir.create(paste0(report_address, "/AMARETTOhtmls/Report_data"), recursive = TRUE, showWarnings = FALSE)
  
  report_data <- list(nExp = nExp,
                      nCNV = nCNV,
                      nMET = nMET,
                      nGenes = nGenes,
                      VarPercentage = VarPercentage,
                      nMod = nMod,
                      ModuleOverviewTable = ModuleOverviewTable,
                      all_genes = all_genes,
                      dt_genesetsall = dt_genesetsall,
                      dt_genesetsall_genetic_pert = dt_genesetsall_genetic_pert,
                      dt_genesetsall_chem_pert = dt_genesetsall_chem_pert,
                      dt_phenotype_association_all = dt_phenotype_association_all,
                      AMARETTOinit = AMARETTOinit,
                      AMARETTOresults = AMARETTOresults)
  
  saveRDS(report_data, file = paste0(report_address, "/AMARETTOhtmls/Report_data/AMARETTOreport_data.rds"))
  cat("The full report is created and ready to use.\n")
  return(report_data)
}

#' Hyper Geometric Geneset Enrichement Test
#'
#' Calculates the p-values for unranked gene set enrichment based on two gmt files as input and the hyper geometric test.
#' @return result
#' @param gmtfile The gmt file with reference gene set.
#' @param testgmtfile The gmt file with gene sets to test. In our case, the gmt file of the modules.
#' @param NrCores Number of cores used for parallelization.
#' @param ref.numb.genes The total number of genes teste, standard equal to 45 956 (MSIGDB standard).
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @keywords internal
HyperGTestGeneEnrichment<-function(gmtfile,testgmtfile,NrCores,ref.numb.genes=45956){
  
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  test.gmt<-readGMT(testgmtfile) # our gmt_file_output_from Amaretto
  gmt.path<-readGMT(gmtfile)  # the hallmarks_and_co2...

  ###########################  Parallelizing :
  cluster <- parallel::makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
  doParallel::registerDoParallel(cluster,cores=NrCores)
  
  resultloop<-foreach(j=1:length(test.gmt), .combine='rbind') %do% {
    foreach(i=1:length(gmt.path),.combine='rbind') %dopar% {
      l<-length(gmt.path[[i]])
      k<-sum(gmt.path[[i]] %in% test.gmt[[j]])
      m<-ref.numb.genes
      n<-length(test.gmt[[j]])
      p1<-stats::phyper(k-1,l,m-l,n,lower.tail=FALSE)
    
      if (k>0){
        overlapping.genes<-gmt.path[[i]][gmt.path[[i]] %in% test.gmt[[j]]]
        overlapping.genes<-paste(overlapping.genes,collapse = ', ')
        # resultloop<-rbind(resultloop,c(Geneset=names(gmt.path[i]),Testset=names(test.gmt[j]),p_value=p1,n_Overlapping=k,Overlapping_genes=overlapping.genes))
        c(Geneset=names(gmt.path[i]),Testset=names(test.gmt[j]),Geneset_length=l,p_value=p1,n_Overlapping=k,Overlapping_genes=overlapping.genes)
      }
    }
  }

  parallel::stopCluster(cluster)
  resultloop<-as.data.frame(resultloop,stringsAsFactors=FALSE)
  resultloop$p_value<-as.numeric(resultloop$p_value)
  resultloop$n_Overlapping<-as.numeric((resultloop$n_Overlapping))
  resultloop$Geneset_length<-as.numeric(resultloop$Geneset_length)
  resultloop[,"padj"]<-stats::p.adjust(resultloop[,"p_value"],method='BH')
  return(resultloop)
}

#' GmtFromModules
#' @return result
#'
#' @param driverGSEA if TRUE , module driver genes will also be added to module target genes for GSEA.
#' @param AMARETTOresults List output from AMARETTO_Run().
#'
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom dplyr arrange mutate select rename  filter 
#' @importFrom utils write.table
#' @keywords internal
GmtFromModules <- function(AMARETTOresults,driverGSEA){
  ModuleMembership <- tibble::rownames_to_column(as.data.frame(AMARETTOresults$ModuleMembership),"GeneNames")
  if(driverGSEA){
    all_regulators <- reshape2::melt(tibble::rownames_to_column(as.data.frame(AMARETTOresults$RegulatoryPrograms),"Module"), id.vars = "Module") %>%
      dplyr::filter(value!=0) %>% dplyr::select(variable, Module) %>% 
      dplyr::mutate(Module = sub("Module_", "", Module)) %>% 
      dplyr::rename(GeneNames = "variable")%>% 
      dplyr::rename(ModuleNr = "Module")
    ModuleMembership <- rbind(ModuleMembership, all_regulators)
  }
  NrModules <- AMARETTOresults$NrModules
  ModuleMembership <- ModuleMembership %>% dplyr::arrange(GeneNames)

  ModuleMembers_list <- split(ModuleMembership$GeneNames,ModuleMembership$ModuleNr)
  names(ModuleMembers_list) <- paste0("Module_",names(ModuleMembers_list))

  gmt_file="./Modules_genes.gmt"
  utils::write.table(sapply(names(ModuleMembers_list), function(x) paste(x,paste(ModuleMembers_list[[x]],collapse="\t"),sep="\t")), gmt_file, quote = FALSE, row.names = TRUE, col.names = FALSE,sep='\t')
}


#' readGMT
#'
#' @param filename
#'
#' @return result
#' @keywords internal
readGMT <- function(filename){
  gmtLines <- strsplit(readLines(filename),"\t")
  gmtLines_genes <- lapply(gmtLines, tail, -2)
  names(gmtLines_genes) <- sapply(gmtLines, head, 1)
  return(gmtLines_genes)
}

#' Title plot_run_history
#'
#' @param AMARETTOinit 
#' @param AMARETTOResults 
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom stats sd
#' @return plot
#' @export
#'
#' @examples  plot_run_history(AMARETTOinit,AMARETTOResults)
plot_run_history <- function(AMARETTOinit,AMARETTOResults){
  means <- unlist(lapply(AMARETTOResults$run_history$error_history, mean))
  stds <- unlist(lapply(AMARETTOResults$run_history$error_history, sd))
  iterationNr <- c(1:length(means))
  NrReassignGenes <- AMARETTOResults$run_history$NrReassignGenes_history[-1]
  threshold <- AMARETTOinit$Parameters$convergence_cutoff*nrow(AMARETTOinit$MA_matrix_Var)
  TotGenesNr <- nrow(AMARETTOinit$MA_matrix_Var)
  
  df <- data.frame(iterationNr = iterationNr,
                 means = means,
                 stds = stds,
                 NrReassignGenes = NrReassignGenes,
                 threshold = threshold,
                 TotGenesNr = TotGenesNr,
                 stringsAsFactors = FALSE)
  p1 <- ggplot2::qplot(x = iterationNr, y = means, data = df) + 
    ggplot2::geom_errorbar(ggplot2::aes(x=iterationNr, ymin=means-stds, ymax=means+stds),data=df,width=0.25) + 
    ggplot2::xlab("Iteration Number") + ggplot2::ylab("Mean Square Error") + 
    ggplot2::geom_line() + 
    ggplot2::geom_point()
  p2 <- ggplot2::qplot(x = iterationNr, y = NrReassignGenes) + 
    ggplot2::geom_hline(yintercept = TotGenesNr, linetype="dashed", color = "blue") + 
    ggplot2::geom_hline(yintercept = threshold, linetype="dashed", color = "red") + 
    ggplot2::xlab("Iteration Number") + 
    ggplot2::ylab("Target Gene Reassignments Number") + 
    ggplot2::geom_line() + 
    ggplot2::geom_point() + 
    ggplot2::scale_y_continuous(trans='log2') 
  
  gridExtra::grid.arrange(p1, p2, nrow = 2)
}

#' Title HyperGeoEnrichmentTest
#'
#' @param AMARETTOresults AMARETTO results output
#' @param hyper_geo_reference GMT file with gene sets to compare with.
#' @param driverGSEA if TRUE, module drivers will also be included in the hypergeometric test.
#' @param NrCores Number of cores for parallel processing. 
#'
#' @return Hyper-Geometric Enrichment Test table
#' @export
#'
#' @examples HyperGeoEnrichmentTest(AMARETTOresults=NULL, hyper_geo_reference, driverGSEA=TRUE, MSIGDB=TRUE, NrCores=4)
HyperGeoEnrichmentTest<-function(AMARETTOresults, hyper_geo_reference, driverGSEA=TRUE, NrCores=4){
  output_hgt_all<-NULL
  for(i in 1:length(hyper_geo_reference)){
    if (is.null(AMARETTOresults)){
      return(1)
    }
    GmtFromModules(AMARETTOresults, driverGSEA)
    output_hgt <- HyperGTestGeneEnrichment(hyper_geo_reference[i], "./Modules_genes.gmt", NrCores)
    utils::data(MsigdbMapping)
    MsigdbMapping<-MsigdbMapping%>%dplyr::mutate(url=paste0('<a href="http://software.broadinstitute.org/gsea/msigdb/cards/',geneset,'">',gsub("_"," ",geneset),'</a>'))
    output_hgt<-output_hgt%>%dplyr::left_join(MsigdbMapping,by=c("Geneset"="geneset"))%>%
      dplyr::mutate(description=ifelse(is.na(description),Geneset,description))%>%
      dplyr::mutate(Geneset=ifelse(is.na(url),Geneset,url))%>%dplyr::rename("Description"="description")%>%dplyr::select(-url)
    cat("The hyper geometric test results are calculated.\n")
    output_hgt_all<-rbind(output_hgt_all,output_hgt)
  }
  return(output_hgt_all)
}

#' Title create_hgt_datatable
#'
#' @param output_hgt GSEA test dataframe from HyperGeoEnrichmentTest function.
#' @param module_table If TRUE, makes the ModuleNr datatable, If FALSE, makes the index page datatable.
#' @param ModuleNr The module number. 
#'
#' @return result
#' @examples 
create_hgt_datatable<-function(output_hgt, module_table, ModuleNr = 1){
  
  if (module_table){
  ##################################################################
    # filter results from module from all datatable with all GSEA results
    output_hgt_filter <- output_hgt %>% dplyr::filter(Testset==paste0("Module_",as.character(ModuleNr))) %>% dplyr::arrange(padj)
    output_hgt_filter <- output_hgt_filter %>% dplyr::mutate(overlap_perc=n_Overlapping/Geneset_length) %>% 
      mutate(overlap_perc=signif(overlap_perc, digits = 3)) %>% 
      dplyr::select(Geneset,Description,Geneset_length, n_Overlapping, Overlapping_genes, overlap_perc, p_value,padj) %>% 
      arrange(padj) %>% 
      mutate(Geneset_length=as.integer(Geneset_length), n_Overlapping=as.integer(n_Overlapping))
    
    filename_table <- paste0("gsea_module",ModuleNr)
    buttons_list <- list(list(extend ='csv',filename=filename_table), list(extend ='excel',filename=filename_table), list(extend = 'pdf', pageSize = 'A4', orientation = 'landscape',filename=filename_table),list(extend ='print'), list(extend ='colvis'))
    #create interactive tables
    options('DT.warn.size'=FALSE)
    
    dt_genesets <- DT::datatable(output_hgt_filter, 
                                 #dplyr::mutate(Geneset=paste0('<a href="http://software.broadinstitute.org/gsea/msigdb/cards/',Geneset,'">',gsub("_"," ",Geneset),'</a>')),
                                 class = 'display', filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE,
                                 options = list(pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',buttons = buttons_list,columnDefs = list(list(className = 'dt-head-center', targets = "_all"),list(className = 'text-left', targets = "_all"))),
                                 colnames=c("Gene Set Name", "Gene Set Description", "# Genes in Gene Set", "# Genes in Overlap", "Genes in Overlap", "% Genes in overlap", "P-value", "FDR Q-value"), escape = FALSE) %>%
      DT::formatSignif(c('p_value','padj','overlap_perc'),2) %>% 
      DT::formatStyle('overlap_perc', background = DT::styleColorBar(c(0,1), 'lightblue'), backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat', backgroundPosition = 'center') %>% 
      DT::formatStyle(columns = c(5), fontSize = '60%')
    
    ngenesets <- nrow(output_hgt_filter %>% dplyr::filter(padj<0.05))
    
    return(list(dt_genesets=dt_genesets,ngenesets=ngenesets))
  } 
  ##################################################################
  else{
    genesetsall<-output_hgt %>% 
      dplyr::mutate(Testset=paste0('<a href="./modules/module',sub("Module_","",Testset),'.html">',paste0(Testset,paste0(rep("&nbsp",14),collapse = "")),'</a>')) %>% 
      dplyr::mutate(Modules=gsub("_","&nbsp",Testset))%>%dplyr::mutate(overlap_perc=n_Overlapping/Geneset_length) %>%
      dplyr::mutate(overlap_perc=signif(overlap_perc, digits = 3))
    
    genesetsall<-genesetsall %>%
      select(Modules, Geneset, Description, Geneset_length, n_Overlapping, Overlapping_genes, overlap_perc, p_value, padj) %>% 
      dplyr::arrange(padj) %>% 
      dplyr::filter(n_Overlapping>2) %>%
      dplyr::mutate(Geneset_length=as.integer(Geneset_length),n_Overlapping=as.integer(n_Overlapping))

    genesetsall<-as.matrix(genesetsall)
    
    filename_table <- "gsea_all_modules"
    buttons_list <- list(list(extend ='csv',filename=filename_table), list(extend ='excel',filename=filename_table), list(extend = 'pdf', pageSize = 'A4', orientation = 'landscape',filename=filename_table),list(extend ='print'), list(extend ='colvis'))
    
    dt_genesetsall<-DT::datatable(genesetsall,class = 'display',filter = 'top', extensions = c('Buttons'), rownames = FALSE,
                                  options = list(deferRender=TRUE,paging =TRUE, pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',buttons = buttons_list,columnDefs = list(list(className = 'dt-head-center', targets = "_all"),list(className = 'text-left', targets = "_all"))),
                                  colnames=c("Module","Gene Set Name","Gene Set Description","# Genes in Gene Set","# Genes in Overlap","Genes in Overlap","% Genes in overlap","P-value","FDR Q-value"),
                                  escape = FALSE) %>%
      DT::formatSignif(c('p_value','padj','overlap_perc'),2) %>% 
      DT::formatStyle('overlap_perc',background = DT::styleColorBar(c(0,1), 'lightblue'),backgroundSize = '98% 88%',backgroundRepeat = 'no-repeat', backgroundPosition = 'center') %>%
      DT::formatStyle(columns = c(6), fontSize = '60%')
    return(dt_genesetsall)
    }
  

}

#' Title driver_genes_summary
#' Provide summary of all driver genes in that 
#'
#' @param AMARETTOresults 
#' @param weight_threshold 
#' @param plot_wordcloud 
#'
#' @return all_regulators_grouped
#' @export
#'
#' @examples 
driver_genes_summary<-function(AMARETTOresults,weight_threshold=0.001, plot_wordcloud=FALSE){
  all_regulators <- reshape2::melt(tibble::rownames_to_column(as.data.frame(AMARETTOresults$RegulatoryPrograms),"Module"),id.vars = "Module") %>% 
    dplyr::filter(value!=0) %>% dplyr::mutate(Module=sub("Module_","",Module),Type="Driver") %>% 
    dplyr::rename(Genes='variable') %>%
    dplyr::select(Genes,Module,value,Type)
  
  all_regulators_grouped<-all_regulators%>%filter(abs(value)>weight_threshold)%>%group_by(Genes)%>%summarise(Modules=paste(Module,collapse = ","),
                                                                       Weights=paste(value,collapse=","))
  all_regulators_grouped$frequency<-unlist(lapply(strsplit(all_regulators_grouped$Modules,","),length))
  all_regulators_grouped<-all_regulators_grouped%>%arrange(-frequency)
  word_cloud<-NA
  if(plot_wordcloud){
    word_cloud<-wordcloud::wordcloud(as.character(all_regulators_grouped$Genes),
                          as.numeric(all_regulators_grouped$frequency),
                          min.freq=1,
                          max.words=Inf,
                          random.order=FALSE,
                          rot.per=.15,
                          colors=brewer.pal(10,"Dark2"))
  }
  return(list(all_regulators=all_regulators,
              all_regulators_grouped=all_regulators_grouped,
              word_cloud=word_cloud))
  
}
