#' AMARETTO_HTMLreport_serial
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
AMARETTO_HTMLreport_serial <- function(AMARETTOinit,
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

  ModuleOverviewTable <- data.frame()
  
  yml_file <- paste0(full_path,"/AMARETTOhtmls/modules/_site.yml")
  file.copy(system.file("templates/module_templates/_site.yml",package="AMARETTO"),yml_file)
  
  #ModuleOverviewTable<-foreach (ModuleNr = 1:NrModules, .packages = c('AMARETTO','tidyverse','DT','rmarkdown')) %dopar% {
  for(ModuleNr in 1:NrModules){
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
    #return(c(ModuleNr, length(which(AMARETTOresults$ModuleMembership==ModuleNr)), length(ModuleRegulators),ngenesets))
    
    # },error=function(e){message(paste("an error occured for Module", ModuleNr))})
    
    ModuleOverviewTable<-rbind(ModuleOverviewTable,
                               data.frame(ModuleNr,
                                          length(which(AMARETTOresults$ModuleMembership==ModuleNr)),
                                          length(ModuleRegulators),
                                          ngenesets))
    }
  
  
  suppressWarnings(suppressMessages(file.remove(paste0(full_path,"/AMARETTOhtmls/modules/_site.yml"))))
  file_remove<-suppressWarnings(suppressMessages(file.remove(paste0(full_path,"/AMARETTOhtmls/modules/module",c(1:NrModules),"_files"))))
  parallel::stopCluster(cluster)
  
  cat("All module htmls are created.\n")
  #ModuleOverviewTable <- data.frame(matrix(unlist(ModuleOverviewTable), byrow=TRUE, ncol=4), stringsAsFactors=FALSE)
  colnames(ModuleOverviewTable)<-c("ModuleNr","NrTarGenes","NrRegGenes","SignGS")
  print(ModuleOverviewTable)
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
