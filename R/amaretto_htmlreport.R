#' AMARETTO_HTMLreport
#'
#' Retrieve an interactive html report, including gene set enrichment analysis if asked for.
#' @param AMARETTOinit AMARETTO initialize output
#' @param AMARETTOresults AMARETTO results output
#' @param ProcessedData List of processed input data
#' @param SAMPLE_annotation SAMPLE annotation will be added to heatmap
#' @param ID ID column of the SAMPLE annotation data frame
#' @param VarPercentage Original Var Percentage used
#' @param hyper_geo_test_bool Boolean if a hyper geometric test needs to be performed. If TRUE provide a GMT file in the hyper_geo_reference parameter.
#' @param hyper_geo_reference GMT file with gene sets to compare with.
#' @param output_address Output directory for the html files.
#' @param MSIGDB TRUE if gene sets were retrieved from MSIGDB. Links will be created in the report.
#' @param GMTURL TRUE if second column of gmt contains URLs to gene set, FALSE if it contains a description
#' @importFrom doParallel registerDoParallel
#' @importFrom DT datatable formatRound formatSignif  formatStyle styleColorBar styleInterval
#' @importFrom reshape2 melt
#' @importFrom dplyr arrange  group_by  left_join mutate  select  summarise
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster stopCluster
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
#'                     VarPercentage=10,hyper_geo_test_bool=FALSE,
#'                     output_address='./')
#'}
AMARETTO_HTMLreport <- function(AMARETTOinit, AMARETTOresults, 
    ProcessedData, SAMPLE_annotation = NULL, ID = NULL, 
    VarPercentage, hyper_geo_test_bool = FALSE, hyper_geo_reference = NULL, 
    output_address = "./", MSIGDB = FALSE, GMTURL = FALSE) {
    
    `%dopar%` <- foreach::`%dopar%`
    CNV_matrix <- ProcessedData[[2]]
    MET_matrix <- ProcessedData[[3]]
    Description <- Overlapping_genes <- overlap_perc <- NULL
    p_value <- value <- variable <- Genes <- Modules <- dt_gensesetsall <- NULL
    dt_gensesetsall  
    NrModules <- AMARETTOresults$NrModules
    NrCores <- AMARETTOinit$NrCores
    if (!dir.exists(output_address)) {
        stop("Output directory is not existing.")
    }
    if (hyper_geo_test_bool == TRUE) {
        if (!file.exists(hyper_geo_reference)) {
            stop("GMT for hyper geometric test is not existing.\n")
        }
    }
    report_address = file.path(output_address, "report_html/")
    dir.create(paste0(report_address, "htmls/modules"), 
        recursive = TRUE, showWarnings = FALSE)
    cat("The output folder structure is created.\n")
    
    if (hyper_geo_test_bool) {
        GmtFromModules(AMARETTOinit, AMARETTOresults)
        output_hgt <- HyperGTestGeneEnrichment(hyper_geo_reference, 
            "./Modules_targets_only.gmt", NrCores)
        GeneSetDescriptions <- GeneSetDescription(hyper_geo_reference)
    }
    cat("The hyper geometric test results are calculated.\n")
    
#    cluster <- parallel::makeCluster(c(rep("localhost", 
#        NrCores)), type = "SOCK")
#    registerDoParallel(cluster, cores = NrCores)
    
    full_path <- normalizePath(report_address)
    ModuleNr <- NULL
    ModuleOverviewTable <- foreach::foreach(ModuleNr = 1:NrModules, 
        .packages = c("AMARETTO", "tidyverse", "DT", 
            "rmarkdown")) %dopar% {
        heatmap_module <- invisible(AMARETTO_VisualizeModule(AMARETTOinit, 
            AMARETTOresults, CNV_matrix, MET_matrix, 
            SAMPLE_annotation = SAMPLE_annotation, 
            ID = ID, ModuleNr = ModuleNr))
        ModuleRegulators <- AMARETTOresults$RegulatoryPrograms[ModuleNr, 
            which(AMARETTOresults$RegulatoryPrograms[ModuleNr, 
                ] != 0)]
        RegulatorIDs<-NULL
        dt_regulators <- DT::datatable(tibble::rownames_to_column(as.data.frame(ModuleRegulators), 
            "RegulatorIDs") %>% dplyr::rename(Weights = "ModuleRegulators") %>% 
            dplyr::mutate(RegulatorIDs = paste0("<a href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene=", 
                RegulatorIDs, "\">", RegulatorIDs, 
                "</a>")), class = "display", extensions = "Buttons", 
            rownames = FALSE, options = list(columnDefs = list(list(width = "200px", 
                targets = "_all")), pageLength = 10, 
                dom = "Bfrtip", buttons = c("csv", 
                  "excel", "pdf")), escape = "Weights") %>% 
            DT::formatRound("Weights", 3) %>% DT::formatStyle("Weights", 
            color = DT::styleInterval(0, c("darkblue", 
                "darkred")))
        
        if (hyper_geo_test_bool) {
            Testset<-NULL
            n_Overlapping<-NULL
            NumberGenes<-NULL
            padj<-NULL
            Geneset<-NULL
            output_hgt_filter <- output_hgt %>% dplyr::filter(Testset == 
                paste0("Module_", as.character(ModuleNr))) %>% 
                dplyr::arrange(padj)
            output_hgt_filter <- dplyr::left_join(output_hgt_filter, 
                GeneSetDescriptions, by = c(Geneset = "GeneSet")) %>% 
                dplyr::mutate(overlap_perc = n_Overlapping/NumberGenes) %>% 
                dplyr::select(Geneset, Description, 
                  n_Overlapping, Overlapping_genes, 
                  overlap_perc, p_value, padj)
            if (MSIGDB == TRUE & GMTURL == FALSE) {
                dt_genesets <- DT::datatable(output_hgt_filter %>% 
                  dplyr::mutate(Geneset = paste0("<a href=\"http://software.broadinstitute.org/gsea/msigdb/cards/", 
                    Geneset, "\">", gsub("_", " ", 
                      Geneset), "</a>")), class = "display", 
                  extensions = "Buttons", rownames = FALSE, 
                  options = list(pageLength = 10, dom = "Bfrtip", 
                    buttons = c("csv", "excel", "pdf")), 
                  colnames = c("Gene Set Name", "Description", 
                    "# Genes in Overlap", "Overlapping Genes", 
                    "Percent of GeneSet overlapping", 
                    "p-value", "FDR q-value"), escape = FALSE) %>% 
                  DT::formatSignif(c("p_value", "padj", 
                    "overlap_perc"), 2) %>% DT::formatStyle("overlap_perc", 
                  background = DT::styleColorBar(c(0, 
                    1), "lightblue"), backgroundSize = "98% 88%", 
                  backgroundRepeat = "no-repeat", backgroundPosition = "center")
            } else if (MSIGDB == TRUE & GMTURL == TRUE) {
                dt_genesets <- DT::datatable(output_hgt_filter %>% 
                  dplyr::select(-Description) %>% dplyr::mutate(Geneset = paste0("<a href=\"http://software.broadinstitute.org/gsea/msigdb/cards/", 
                  Geneset, "\">", gsub("_", " ", Geneset), 
                  "</a>")), class = "display", extensions = "Buttons", 
                  rownames = FALSE, options = list(pageLength = 10, 
                    dom = "Bfrtip", buttons = c("csv", 
                      "excel", "pdf")), colnames = c("Gene Set Name", 
                    "# Genes in Overlap", "Overlapping Genes", 
                    "Percent of GeneSet overlapping", 
                    "p-value", "FDR q-value"), escape = FALSE) %>% 
                  DT::formatSignif(c("p_value", "padj", 
                    "overlap_perc"), 2) %>% DT::formatStyle("overlap_perc", 
                  background = DT::styleColorBar(c(0, 
                    1), "lightblue"), backgroundSize = "98% 88%", 
                  backgroundRepeat = "no-repeat", backgroundPosition = "center")
            } else if (MSIGDB == FALSE & GMTURL == TRUE) {
                dt_genesets <- DT::datatable(output_hgt_filter %>% 
                  dplyr::mutate(Geneset = paste0("<a href=\"", 
                    Description, "\">", gsub("_", " ", 
                      Geneset), "</a>")) %>% dplyr::select(-Description), 
                  class = "display", extensions = "Buttons", 
                  rownames = FALSE, options = list(pageLength = 10, 
                    dom = "Bfrtip", buttons = c("csv", 
                      "excel", "pdf")), colnames = c("Gene Set Name", 
                    "Description", "# Genes in Overlap", 
                    "Overlapping Genes", "Percent of GeneSet overlapping", 
                    "p-value", "FDR q-value"), escape = FALSE) %>% 
                  DT::formatSignif(c("p_value", "padj", 
                    "overlap_perc"), 2) %>% DT::formatStyle("overlap_perc", 
                  background = DT::styleColorBar(c(0, 
                    1), "lightblue"), backgroundSize = "98% 88%", 
                  backgroundRepeat = "no-repeat", backgroundPosition = "center")
            } else {
                dt_genesets <- DT::datatable(output_hgt_filter, 
                  class = "display", extensions = "Buttons", 
                  rownames = FALSE, options = list(autoWidth = TRUE, 
                    pageLength = 10, dom = "Bfrtip", 
                    buttons = c("csv", "excel", "pdf"))) %>% 
                  DT::formatSignif(c("p_value", "padj", 
                    "overlap_perc"), 2)
            }
            ngenesets <- nrow(output_hgt_filter %>% 
                dplyr::filter(padj < 0.05))
        } else {
            dt_genesets <- "Genesets were not analysed as they were not provided."
            ngenesets <- "NA"
        }
        modulemd <- paste0(full_path, "/htmls/modules/module", 
            ModuleNr, ".rmd")
        file.copy(system.file("templates/TemplateReportModule.Rmd", 
            package = "AMARETTO"), modulemd)
        rmarkdown::render(modulemd, output_file = paste0("module", 
            ModuleNr, ".html"), params = list(report_address = report_address, 
            ModuleNr = ModuleNr, heatmap_module = heatmap_module, 
            dt_regulators = dt_regulators, dt_genesets = dt_genesets), 
            quiet = TRUE)
        file.remove(modulemd)
        return(c(ModuleNr, length(which(AMARETTOresults$ModuleMembership == 
            ModuleNr)), length(ModuleRegulators), ngenesets))
    }
    
#    parallel::stopCluster(cluster)
    cat("The module htmls are finished.\n")
    ModuleOverviewTable <- data.frame(matrix(unlist(ModuleOverviewTable), 
        byrow = TRUE, ncol = 4), stringsAsFactors = FALSE)
    colnames(ModuleOverviewTable) <- c("ModuleNr", 
        "NrTarGenes", "NrRegGenes", "SignGS")
    if (!is.null(CNV_matrix)) {
        nCNV = ncol(CNV_matrix)
    } else {
        nCNV = NA
    }
    if (!is.null(MET_matrix)) {
        nMET = ncol(MET_matrix)
    } else {
        nMET = NA
    }
    nExp = ncol(AMARETTOresults$RegulatoryProgramData)
    nGenes = length(AMARETTOresults$AllGenes)
    nMod = AMARETTOresults$NrModules
    
    dt_overview <- DT::datatable(ModuleOverviewTable %>% 
        dplyr::mutate(ModuleNr = paste0("<a href=\"./modules/module", 
            ModuleNr, ".html\">Module ", ModuleNr, 
            "</a>")), class = "display", extensions = "Buttons", 
        rownames = FALSE, options = list(pageLength = 10, 
            dom = "Bfrtip", buttons = c("csv", "excel", 
                "pdf")), escape = FALSE)
    
    all_targets <- tibble::rownames_to_column(data.frame(AMARETTOinit$ModuleMembership), 
        "Genes") %>% dplyr::rename(Module = "AMARETTOinit.ModuleMembership") %>% 
        dplyr::mutate(Type = "Target")
    all_regulators <- reshape2::melt(tibble::rownames_to_column(as.data.frame(AMARETTOresults$RegulatoryPrograms), 
        "Module"), id.vars = "Module") %>% dplyr::filter(value > 
        0) %>% dplyr::select(variable, Module) %>% 
        dplyr::mutate(Module = sub("Module_", "", Module), 
            Type = "Regulator") %>% dplyr::rename(Genes = "variable")
    all_genes <- rbind(all_targets, all_regulators) %>% 
        dplyr::arrange(Genes) %>% dplyr::mutate(Module = paste0("<a href=\"./modules/module", 
        Module, ".html\">Module ", Module, "</a>"))
    dt_genes <- DT::datatable(all_genes, class = "display", 
        extensions = "Buttons", rownames = FALSE, options = list(pageLength = 10, 
            dom = "Bfrtip", buttons = c("csv", "excel", 
                "pdf")), escape = FALSE)
    if (hyper_geo_test_bool) {
        dt_genesetsall <- DT::datatable(dplyr::left_join(output_hgt %>% 
            dplyr::filter(padj < 0.05 & n_Overlapping > 
                1) %>% dplyr::group_by(Geneset) %>% 
            dplyr::mutate(Testset = paste0("<a href=\"./modules/module", 
                sub("Module_", "", Testset), ".html\">", 
                Testset, "</a>")) %>% dplyr::summarise(Modules = paste(Testset, 
            collapse = ", ")), GeneSetDescriptions, 
            by = c(Geneset = "GeneSet")) %>% dplyr::mutate(Geneset = gsub("_", 
            " ", Geneset), Modules = gsub("_", " ", 
            Modules)), class = "display", extensions = "Buttons", 
            rownames = FALSE, options = list(pageLength = 10, 
                dom = "Bfrtip", buttons = c("csv", 
                  "excel", "pdf")), escape = FALSE)
    } else {
        dt_genesetsall <- "Genesets were not analysed as they were not provided."
    }
    
    rmarkdown::render(system.file("templates/TemplateIndexPage.Rmd", 
        package = "AMARETTO"), output_dir = paste0(full_path, 
        "/htmls/"), output_file = "index.html", params = list(nExp = nExp, 
        nCNV = nCNV, nMET = nMET, nGenes = nGenes, 
        VarPercentage = VarPercentage, nMod = nMod, 
        dt_overview = dt_overview, dt_genes = dt_genes, 
        dt_genesetsall = dt_gensesetsall), quiet = TRUE)
    cat("The report is ready to use\n")
}

#' Hyper Geometric Geneset Enrichement Test
#'
#' Calculates the p-values for unranked gene set enrichment based on two gmt files as input and the hyper geometric test.
#' @return result
#' @import doParallel
#' @keywords internal
HyperGTestGeneEnrichment <- function(gmtfile, testgmtfile, 
    NrCores, ref.numb.genes = 45956) {
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`
    test.gmt <- readGMT(testgmtfile)  # our gmt_file_output_from Amaretto
    gmt.path <- readGMT(gmtfile)  # the hallmarks_and_co2...
    i <- j <- NULL
    ########################### Parallelizing :
#    cluster <- parallel::makeCluster(c(rep("localhost", 
#        NrCores)), type = "SOCK")
#    registerDoParallel(cluster, cores = NrCores)
    
    resultloop <- foreach::foreach(j = 1:length(test.gmt), 
        .combine = "rbind") %do% {
        # print(j)
        foreach::foreach(i = 1:length(gmt.path), .combine = "rbind") %dopar%
            {
                # print(i)
                l <- length(gmt.path[[i]])
                k <- sum(gmt.path[[i]] %in% test.gmt[[j]])
                m <- ref.numb.genes
                n <- length(test.gmt[[j]])
                p1 <- stats::phyper(k - 1, l, m - l, 
                  n, lower.tail = FALSE)
                
                if (k > 0) {
                  overlapping.genes <- gmt.path[[i]][gmt.path[[i]] %in% 
                    test.gmt[[j]]]
                  overlapping.genes <- paste(overlapping.genes, 
                    collapse = ", ")
                  c(Geneset = names(gmt.path[i]), Testset = names(test.gmt[j]), 
                    p_value = p1, n_Overlapping = k, 
                    Overlapping_genes = overlapping.genes)
                }
            }
    }
    
#    parallel::stopCluster(cluster)
    resultloop <- as.data.frame(resultloop, stringsAsFactors = FALSE)
    resultloop$p_value <- as.numeric(resultloop$p_value)
    resultloop$n_Overlapping <- as.numeric((resultloop$n_Overlapping))
    resultloop[, "padj"] <- stats::p.adjust(resultloop[, 
        "p_value"], method = "BH")
    return(resultloop)
}

#' GmtFromModules
#' @return result
#' @keywords internal
GmtFromModules <- function(AMARETTOinit, AMARETTOresults) {
    ModuleMembership <- tibble::rownames_to_column(as.data.frame(AMARETTOresults$ModuleMembership), 
        "GeneNames")
    GeneNames <- NULL
    NrModules <- AMARETTOresults$NrModules
    ModuleMembership <- ModuleMembership %>% dplyr::arrange(GeneNames)
    
    ModuleMembers_list <- split(ModuleMembership$GeneNames, 
        ModuleMembership$ModuleNr)
    names(ModuleMembers_list) <- paste0("Module_", 
        names(ModuleMembers_list))
    
    gmt_file = "./Modules_targets_only.gmt"
    utils::write.table(sapply(names(ModuleMembers_list), 
        function(x) paste(x, paste(ModuleMembers_list[[x]], 
            collapse = "\t"), sep = "\t")), gmt_file, 
        quote = FALSE, row.names = TRUE, col.names = FALSE, 
        sep = "\t")
}

#' GeneSetDescription
#'
#' @return result
#' @keywords internal
GeneSetDescription <- function(filename) {
    gmtLines <- strsplit(readLines(filename), "\t")
    gmtLines_description <- lapply(gmtLines, function(x) {
        c(x[[1]], x[[2]], length(x) - 2)
    })
    gmtLines_description <- data.frame(matrix(unlist(gmtLines_description), 
        byrow = TRUE, ncol = 3), stringsAsFactors = FALSE)
    rownames(gmtLines_description) <- NULL
    colnames(gmtLines_description) <- c("GeneSet", 
        "Description", "NumberGenes")
    gmtLines_description$NumberGenes <- as.numeric(gmtLines_description$NumberGenes)
    return(gmtLines_description)
}

#' readGMT
#'
#' @return result
#' @keywords internal
#'
readGMT <- function(filename) {
    gmtLines <- strsplit(readLines(filename), "\t")
    gmtLines_genes <- lapply(gmtLines, utils::tail, -2)
    names(gmtLines_genes) <- sapply(gmtLines, utils::head, 
        1)
    return(gmtLines_genes)
}
