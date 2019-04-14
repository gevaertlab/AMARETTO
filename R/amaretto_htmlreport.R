#' AMARETTO_HTMLreport
#'
#' Retrieve an interactive html report, including gene set enrichment analysis if asked for.
#'
#' @param AMARETTOinit AMARETTO initialize output
#' @param AMARETTOresults AMARETTO results output
#' @param ProcessedData List of processed input data
#' @param SAMPLE_annotation SAMPLE annotation will be added to heatmap
#' @param ID ID column of the SAMPLE annotation data frame
#' @param hyper_geo_test_bool Boolean if a hyper geometric test needs to be performed. If TRUE provide a GMT file in the hyper_geo_reference parameter.
#' @param hyper_geo_reference GMT file with gene sets to compare with.
#' @param output_address Output directory for the html files.
#' @param show_row_names if True, sample names will appear in the heatmap
#' @param driverGSEA if TRUE, module drivers will also be included in the hypergeometric test.
#' @param MSIGDB TRUE if gene sets were retrieved from MSIGDB. Links will be created in the report.
#' @param hyper_geo_filters filter values for hyper-geometric test, including the max number of reference geneset,max p-value, max q-value.
#'
#' @import dplyr
#' @importFrom doParallel registerDoParallel
#' @importFrom DT datatable formatRound formatSignif  formatStyle styleColorBar styleInterval
#' @importFrom reshape2 melt
#' @importFrom dplyr arrange group_by left_join mutate select summarise  rename  filter 
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom knitr knit_meta
#' @importFrom utils  write.table
#' @importFrom tibble rownames_to_column
#' @importFrom stats  p.adjust  phyper
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
                                hyper_geo_test_bool = FALSE,
                                hyper_geo_reference = NULL,
                                hyper_geo_filters = list(maxGeneLength = 1000, minPval= 0.05, minQval= 0.1),
                                output_address = './',
                                MSIGDB = TRUE,
                                driverGSEA = TRUE){
  
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  CNV_matrix <- ProcessedData[[2]]
  MET_matrix <- ProcessedData[[3]]
  NrModules<-AMARETTOresults$NrModules
  VarPercentage<-AMARETTOinit$Parameters$VarPercentage
  NrCores<-AMARETTOinit$NrCores
  if (!dir.exists(output_address)){
    stop("Output directory is not existing.")
  }
  if (hyper_geo_test_bool==TRUE){
    if (!file.exists(hyper_geo_reference)){
    stop("GMT for hyper geometric test is not existing.\n")
   }
  }
  report_address <- file.path(output_address)
  dir.create(paste0(report_address,"/AMARETTOhtmls/modules"),recursive = TRUE,showWarnings = FALSE)
  cat("The output folder structure is created.\n")
  if (hyper_geo_test_bool){
    GmtFromModules(AMARETTOinit,AMARETTOresults,driverGSEA)
    output_hgt<-HyperGTestGeneEnrichment(hyper_geo_reference, "./Modules_genes.gmt",NrCores)
    if(!is.null(hyper_geo_filters)){
      output_hgt<-output_hgt%>%filter(as.numeric(Geneset_length) < hyper_geo_filters$maxGeneLength)%>%filter(as.numeric(p_value) < hyper_geo_filters$minPval)%>%filter(as.numeric(padj) < hyper_geo_filters$minQval)
    }
    GeneSetDescriptions<-GeneSetDescription(hyper_geo_reference,MSIGDB)
  }
  cat("The hyper geometric test results are calculated.\n")
  cluster <- parallel::makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
  doParallel::registerDoParallel(cluster,cores=NrCores)

  full_path<-normalizePath(report_address)
  ModuleOverviewTable<-NULL
  ModuleOverviewTable<-foreach (ModuleNr = 1:NrModules, .packages = c('AMARETTO','tidyverse','DT','rmarkdown')) %dopar% {
  #for(ModuleNr in 1:NrModules){

    print(paste0("ModuleNr = ",ModuleNr))
    heatmap_module<-AMARETTO_VisualizeModule(AMARETTOinit, AMARETTOresults, ProcessedData, show_row_names = show_row_names, SAMPLE_annotation=SAMPLE_annotation, ID=ID, ModuleNr=ModuleNr)
    print("visualization is done")
    ModuleRegulators <- AMARETTOresults$RegulatoryPrograms[ModuleNr,which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
    print("ModuleRegulators is done")
    dt_regulators<-DT::datatable(tibble::rownames_to_column(as.data.frame(ModuleRegulators),"RegulatorIDs") %>% dplyr::rename(Weights="ModuleRegulators") %>% dplyr::mutate(RegulatorIDs=paste0('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=',RegulatorIDs,'">',RegulatorIDs,'</a>'))%>%dplyr::arrange(Weights),
                             class = 'display',filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE, options = list(
                               columnDefs = list(list(width = '200px', targets = "_all")),pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',
                               buttons = c('csv', 'excel', 'pdf', 'colvis')),escape = 'Weights') %>% DT::formatRound('Weights',3) %>% DT::formatStyle('Weights',color = DT::styleInterval(0, c('darkblue', 'darkred')))
    print("dt_regulators is done")
    dt_targets<-DT::datatable(as.data.frame(AMARETTOresults$ModuleMembership) %>% tibble::rownames_to_column("TargetIDs")%>% dplyr::arrange(TargetIDs) %>% dplyr::rename(moduleNr=ModuleNr) %>% dplyr::filter(moduleNr==ModuleNr) %>% dplyr::select(-moduleNr) %>% dplyr::mutate(TargetIDs=paste0('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=',TargetIDs,'">',TargetIDs,'</a>')),
                             class = 'display', filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE, options = list(
                               columnDefs = list(list(width = '200px', targets = "_all")),pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',
                               buttons = c('csv', 'excel', 'pdf', 'colvis')),escape = FALSE)
    print("dt_targets is done")
    if (hyper_geo_test_bool){
      output_hgt_filter<-output_hgt %>% dplyr::filter(Testset==paste0("Module_",as.character(ModuleNr))) %>% dplyr::arrange(padj)
      output_hgt_filter<-dplyr::left_join(output_hgt_filter,GeneSetDescriptions,by=c("Geneset"="GeneSet")) %>% dplyr::mutate(overlap_perc=n_Overlapping/NumberGenes) %>% dplyr::select(Geneset,Description,n_Overlapping,Overlapping_genes,overlap_perc,p_value,padj)
      if (MSIGDB==TRUE){
        dt_genesets<-DT::datatable(output_hgt_filter %>% dplyr::mutate(Geneset=paste0('<a href="http://software.broadinstitute.org/gsea/msigdb/cards/',Geneset,'">',gsub("_"," ",Geneset),'</a>')),class = 'display', filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE,
                                   options = list(pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',buttons = c('csv', 'excel', 'pdf', 'colvis')),
          colnames=c("Gene Set Name","Description","# Genes in Overlap","Overlapping Genes","Percent of GeneSet overlapping","p-value","FDR q-value"),escape = FALSE) %>%
          DT::formatSignif(c('p_value','padj','overlap_perc'),2) %>% DT::formatStyle('overlap_perc',background = DT::styleColorBar(c(0,1), 'lightblue'),backgroundSize = '98% 88%',backgroundRepeat = 'no-repeat', backgroundPosition = 'center')
      } 
      else{
        dt_genesets<-DT::datatable(output_hgt_filter,class = 'display', filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE,options = list(
          pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',
          buttons = c('csv', 'excel', 'pdf', 'colvis'))) %>% DT::formatSignif(c('p_value','padj','overlap_perc'),2)
      }
      ngenesets<-nrow(output_hgt_filter %>% dplyr::filter(padj<0.05))
    } else {
      dt_genesets<-"Genesets were not analysed as they were not provided."
      ngenesets<-"NA"
    }
    print("hypergeotest is done")
    # if (!is.null(phenotype_association_table)){
    #   module_number<-ModuleNr
    #   module_phenotype_association_datatable<-datatable(phenotype_association_table%>%filter(ModuleNr==paste0("Module ",module_number))%>%
    #                                                       select(-ModuleNr),class='display',extensions = 'Buttons',rownames = FALSE,options = list(
    #                                                         pageLength = 10,dom = 'Bfrtip',buttons = c('csv', 'excel', 'pdf'),colnames=c("Phenotypes","Statistical Test","p-value","FDR q-value","Descriptive Statistics"),escape = FALSE))%>%
    #                                                         formatSignif(c('p.value','q.value'),2)
    # }
    # else{
    #   module_phenotype_association_datatable<-"Phenotype association resuls were not provided."
    # }
    # print("phenotype is done!")
    modulemd<-paste0(full_path,"/AMARETTOhtmls/modules/module",ModuleNr,".rmd")
    file.copy(system.file("templates/TemplateReportModule.Rmd",package="AMARETTO"),modulemd)
    print("file.copy is done!")
    knitr::knit_meta(class=NULL, clean = TRUE)
    rmarkdown::render(modulemd,output_file = paste0("module",ModuleNr,".html"), params = list(
      report_address = report_address,
      ModuleNr = ModuleNr,
      heatmap_module = heatmap_module,
      dt_regulators = dt_regulators,
      dt_targets = dt_targets,
      #module_phenotype_association_datatable=module_phenotype_association_datatable,
      dt_genesets = dt_genesets),quiet = TRUE)
    print("rmarkdown is done and module html is created :)")
    file.remove(modulemd)
    print("file removed successfully :) Done!")
    #ModuleOverviewTable<-rbind(ModuleOverviewTable,c(ModuleNr,length(which(AMARETTOresults$ModuleMembership==ModuleNr)),length(ModuleRegulators),ngenesets))
    return(c(ModuleNr,length(which(AMARETTOresults$ModuleMembership==ModuleNr)),length(ModuleRegulators),ngenesets))
    dev.off()
    # },error=function(e){message(paste("an error occured for Module", ModuleNr))})
  }


  parallel::stopCluster(cluster)
  cat("The module htmls are finished.\n")
  ModuleOverviewTable<-data.frame(matrix(unlist(ModuleOverviewTable),byrow=TRUE,ncol=4),stringsAsFactors=FALSE)
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
  dt_overview<-DT::datatable(ModuleOverviewTable %>% dplyr::mutate(ModuleNr=paste0('<a href="./modules/module',ModuleNr,'.html">Module ',ModuleNr,'</a>')),class = 'display',filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE,options = list(
    pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',buttons = c('csv', 'excel', 'pdf', 'colvis')),escape = FALSE)
  
  all_targets<-tibble::rownames_to_column(data.frame(AMARETTOresults$ModuleMembership),"Genes") %>% dplyr::rename(Module="ModuleNr") %>% dplyr::mutate(Type="Target")
  all_regulators<-reshape2::melt(tibble::rownames_to_column(as.data.frame(AMARETTOresults$RegulatoryPrograms),"Module"),id.vars = "Module") %>% dplyr::filter(value>0) %>% dplyr::select(variable,Module) %>% dplyr::mutate(Module=sub("Module_","",Module),Type="Regulator") %>% dplyr::rename(Genes='variable')
  all_genes<-rbind(all_targets,all_regulators) %>% dplyr::arrange(Genes) %>% dplyr::mutate(Genes=paste0('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=',Genes,'">',Genes,'</a>')) %>% dplyr::mutate(Module=paste0('<a href="./modules/module',Module,'.html">Module ',Module,'</a>'))
  
  dt_genes<-DT::datatable(all_genes, class = 'display',filter = 'top',extensions = c('Buttons','KeyTable'), rownames = FALSE,options = list(
    pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',buttons = c('csv', 'excel', 'pdf', 'colvis')),escape = FALSE)
  if (hyper_geo_test_bool){
    genesetsall<-dplyr::left_join(output_hgt %>% dplyr::filter(padj<0.05 & n_Overlapping>1) %>% dplyr::group_by(Geneset) %>% dplyr::mutate(Testset=paste0('<a href="./modules/module',sub("Module_","",Testset),'.html">',Testset,'</a>')) %>% dplyr::summarise(Modules=paste(Testset,collapse=", ")),GeneSetDescriptions,by=c("Geneset"="GeneSet")) %>% dplyr::mutate(Modules=gsub("_"," ",Modules))
    if (MSIGDB==TRUE){
      genesetsall<-dplyr::mutate(genesetsall,Geneset=paste0('<a href="http://software.broadinstitute.org/gsea/msigdb/cards/',Geneset,'">',gsub("_"," ",Geneset),'</a>'))
    }
    dt_genesetsall<-DT::datatable(genesetsall,class = 'display',filter = 'top', extensions = c('Buttons','KeyTable'), rownames = FALSE,
                              options = list(pageLength = 10, lengthMenu = c(5, 10, 20, 50, 100), keys = TRUE, dom = 'Blfrtip',buttons = c('csv', 'excel', 'pdf', 'colvis')),escape = FALSE)
  }else{
    dt_genesetsall<-"Genesets were not analysed as they were not provided."
  }
  # if (!is.null(phenotype_association_table)){
  #   phenotype_association_datatable<-datatable(phenotype_association_table%>% mutate(ModuleNr=paste0('<a href="./modules/module',gsub("Module ","",ModuleNr),'.html">',ModuleNr,'</a>')),class='display',extensions = 'Buttons',rownames = FALSE,
  #                                              options = list(pageLength = 10,dom = 'Bfrtip',buttons = c('csv', 'excel', 'pdf'),colnames=c("ModuleNr","Phenotypes","Statistical Test","p-value","FDR q-value","Descriptive Statistics")),escape = FALSE)%>%formatSignif(c('p.value','q.value'),2)
  # }
  # else{
  #   phenotype_association_datatable<-"Phenotype association resuls were not provided."
  # }
  rmarkdown::render(system.file("templates/TemplateIndexPage.Rmd",package="AMARETTO"), output_dir=paste0(full_path,"/AMARETTOhtmls/"),output_file= "index.html", params = list(
    nExp = nExp,
    nCNV = nCNV,
    nMET = nMET,
    nGenes = nGenes,
    VarPercentage = VarPercentage,
    nMod = nMod,
    dt_overview = dt_overview,
    dt_genes=dt_genes,
    #phenotype_association_datatable=phenotype_association_datatable,
    dt_genesetsall = dt_gensesetsall),quiet = TRUE)
  cat("The report is ready to use\n")
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
  
  #resultloop<-c()
  resultloop<-foreach(j=1:length(test.gmt), .combine='rbind') %do% {
    #print(j)
    foreach(i=1:length(gmt.path),.combine='rbind') %dopar% {
      #print(i)
  # for(j in 1:length(test.gmt)){
  #   print(paste0("test_gmt = ",j))
  #   for(i in 1:length(gmt.path)){
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
  resultloop[,"padj"]<-stats::p.adjust(resultloop[,"p_value"],method='BH')
  return(resultloop)
}

#' GmtFromModules
#' @return result
#'
#' @param AMARETTOinit List output from AMARETTO_Initialize().
#' @param driverGSEA if TRUE , module driver genes will also be added to module target genes for GSEA.
#' @param AMARETTOresults List output from AMARETTO_Run().
#'
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom dplyr arrange mutate select rename  filter 
#' @importFrom utils write.table
#' @keywords internal
GmtFromModules <- function(AMARETTOinit,AMARETTOresults,driverGSEA){
  ModuleMembership<-tibble::rownames_to_column(as.data.frame(AMARETTOresults$ModuleMembership),"GeneNames")
  if(driverGSEA){
    all_regulators <-reshape2::melt(tibble::rownames_to_column(as.data.frame(AMARETTOresults$RegulatoryPrograms),"Module"), id.vars = "Module") %>%
      dplyr::filter(value > 0) %>% dplyr::select(variable, Module) %>% dplyr::mutate(Module = sub("Module_", "", Module)) %>% dplyr::rename(GeneNames = "variable")%>% dplyr::rename(ModuleNr = "Module")
    ModuleMembership<-rbind(ModuleMembership,all_regulators)
  }
  NrModules<-AMARETTOresults$NrModules
  ModuleMembership<-ModuleMembership %>% dplyr::arrange(GeneNames)

  ModuleMembers_list<-split(ModuleMembership$GeneNames,ModuleMembership$ModuleNr)
  names(ModuleMembers_list)<-paste0("Module_",names(ModuleMembers_list))

  gmt_file="./Modules_genes.gmt"
  utils::write.table(sapply(names(ModuleMembers_list),function(x) paste(x,paste(ModuleMembers_list[[x]],collapse="\t"),sep="\t")),gmt_file,quote = FALSE,row.names = TRUE,col.names = FALSE,sep='\t')
}

#' GeneSetDescription
#'
#' @param filename The name of the gmt file.
#' @param MSIGDB 
#'
#' @importFrom utils data
#' @return result
#' @keywords internal
GeneSetDescription<-function(filename,MSIGDB){
  utils::data(MsigdbMapping)
  gmtLines<-strsplit(readLines(filename),"\t")
  gmtLines_description <- lapply(gmtLines, function(x) {
    c(x[[1]],x[[2]],length(x)-2)
  })
  gmtLines_description<-data.frame(matrix(unlist(gmtLines_description),byrow=TRUE,ncol=3),stringsAsFactors=FALSE)
  rownames(gmtLines_description)<-NULL
  colnames(gmtLines_description)<-c("GeneSet","Description","NumberGenes")
  gmtLines_description$NumberGenes<-as.numeric(gmtLines_description$NumberGenes)
  if(MSIGDB){
  gmtLines_description$Description<-sapply(gmtLines_description$GeneSet, function(x) {
    index<-which(MsigdbMapping$geneset==x)
    ifelse(length(index)!=0, MsigdbMapping$description[index],gmtLines_description$Description[which(gmtLines_description$GeneSet==x)]) 
  })}
  return(gmtLines_description)
}

#' readGMT
#'
#' @param filename
#'
#' @return result
#' @keywords internal

readGMT<-function(filename){
  gmtLines<-strsplit(readLines(filename),"\t")
  gmtLines_genes <- lapply(gmtLines, tail, -2)
  names(gmtLines_genes) <- sapply(gmtLines, head, 1)
  return(gmtLines_genes)
}
