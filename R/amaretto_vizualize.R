#' AMARETTO_VisualizeModule
#'
#' Function to visualize the gene modules
#' @param AMARETTOinit List output from AMARETTO_Initialize().
#' @param AMARETTOresults List output from AMARETTO_Run().
#' @param ProcessedData List of processed input data
#' @param ModuleNr Module number to visualize
#' @param SAMPLE_annotation Matrix or Dataframe with sample annotation
#' @param ID Column used as sample name
#' @param order_samples Order samples in heatmap by mean or by clustering
#' @param printHM Boolean print heatmap directly
#' @importFrom circlize colorRamp2  rand_color
#' @importFrom grid gpar unit
#' @importFrom stats dist hclust
#' @import dplyr
#' @import grDevices
#' @import methods
#' @import ComplexHeatmap
#' @importFrom tibble column_to_rownames  rownames_to_column
#' @return result
#' @export
#'
#' @examples
#' data('ProcessedDataLIHC')
#' AMARETTOinit <- AMARETTO_Initialize(ProcessedData = ProcessedDataLIHC,
#'                                     NrModules = 2, VarPercentage = 50)
#'
#' AMARETTOresults <- AMARETTO_Run(AMARETTOinit)
#'
#' AMARETTO_VisualizeModule(AMARETTOinit = AMARETTOinit,AMARETTOresults = AMARETTOresults,
#'                          ProcessedData = ProcessedDataLIHC, ModuleNr = 1)
AMARETTO_VisualizeModule <- function(AMARETTOinit, 
    AMARETTOresults, ProcessedData, ModuleNr, SAMPLE_annotation = NULL, 
    ID = NULL, order_samples = NULL, printHM = FALSE) {
    CNV_matrix <- ProcessedData[[2]]
    MET_matrix <- ProcessedData[[3]]
    
    if (ModuleNr > AMARETTOresults$NrModules) {
        stop("\tCannot plot Module", ModuleNr, " since the total number of modules is", 
            AMARETTOresults$N, ".\n")
    }
    ModuleData <- as.data.frame(AMARETTOinit$MA_matrix_Var)[AMARETTOresults$ModuleMembership == 
        ModuleNr, ]
    ModuleRegulators <- AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr, 
        ] != 0)]
    RegulatorData <- as.data.frame(AMARETTOinit$RegulatorData)[ModuleRegulators, 
        ]
    ModuleGenes <- rownames(ModuleData)
    cat("Module", ModuleNr, "has", length(rownames(ModuleData)), 
        "genes and", length(ModuleRegulators), "regulators for", 
        length(colnames(ModuleData)), "samples.\n")
    Alterations <- tibble::rownames_to_column(as.data.frame(AMARETTOinit$RegulatorAlterations$Summary), 
        "HGNC_symbol") %>% dplyr::rename(DriverList = "Driver List") %>% 
        dplyr::filter(HGNC_symbol %in% ModuleRegulators)
    Alterations <- Alterations %>% dplyr::mutate(CNVMet_Alterations = dplyr::case_when(MET == 
        1 & CNV == 1 ~ "Methylation and copy number alterations", 
        CNV == 1 ~ "Copy number alterations", MET == 
            1 ~ "Methylation aberrations", MET == 0 & 
            CNV == 0 ~ "Not Altered"), DriversList_Alterations = dplyr::case_when(DriverList == 
        0 ~ "Driver not predefined", DriverList == 
        1 ~ "Driver predefined"))
    
    ha_drivers <- ComplexHeatmap::HeatmapAnnotation(df = tibble::column_to_rownames(Alterations, 
        "HGNC_symbol") %>% dplyr::select(CNVMet_Alterations, 
        DriversList_Alterations), col = list(CNVMet_Alterations = c(`Copy number alterations` = "#eca400", 
        `Methylation aberrations` = "#006992", `Methylation and copy number alterations` = "#d95d39", 
        `Not Altered` = "lightgray"), DriversList_Alterations = c(`Driver not predefined` = "lightgray", 
        `Driver predefined` = "#588B5B")), which = "column", 
        height = grid::unit(0.3, "cm"), name = "", 
        annotation_legend_param = list(title_gp = grid::gpar(fontsize = 8), 
            labels_gp = grid::gpar(fontsize = 6)))
    
    if (is.null(MET_matrix) && is.null(CNV_matrix)) {
        overlapping_samples <- colnames(ModuleData)
    } else if (is.null(MET_matrix)) {
        overlapping_samples <- Reduce(intersect, list(colnames(CNV_matrix), 
            colnames(ModuleData)))
    } else if (is.null(CNV_matrix)) {
        overlapping_samples <- Reduce(intersect, list(colnames(MET_matrix), 
            colnames(ModuleData)))
    } else {
        overlapping_samples <- Reduce(intersect, list(colnames(CNV_matrix), 
            colnames(MET_matrix), colnames(ModuleData)))
    }
    
    if (is.null(order_samples)) {
        overlapping_samples_clust <- overlapping_samples[order(colMeans(ModuleData[, 
            overlapping_samples]))]
    } else if (order_samples == "clust") {
        SampleClustering <- stats::hclust(stats::dist(t(ModuleData[, 
            overlapping_samples])), method = "complete", 
            members = NULL)
        overlapping_samples_clust <- overlapping_samples[SampleClustering$order]
    } else {
        print("ordering type not recognized, samples will be orderd based on mean expression of the module genes")
        overlapping_samples_clust <- overlapping_samples[order(colMeans(ModuleData[, 
            overlapping_samples]))]
    }
    
    ClustRegulatorData <- t(RegulatorData[, overlapping_samples_clust])
    ClustModuleData <- t(ModuleData[, overlapping_samples_clust])
    Regwidth <- ncol(ClustRegulatorData) * 0.5
    ha_Reg <- ComplexHeatmap::Heatmap(ClustRegulatorData, 
        name = "Gene Expression", column_title = "Regulator Genes\nExpression", 
        cluster_rows = FALSE, cluster_columns = TRUE, 
        show_column_dend = FALSE, show_column_names = TRUE, 
        show_row_names = FALSE, column_names_gp = grid::gpar(fontsize = 6), 
        top_annotation = ha_drivers, column_title_gp = grid::gpar(fontsize = 6, 
            fontface = "bold"), col = circlize::colorRamp2(c(-max(abs(ClustRegulatorData)), 
            0, max(abs(ClustRegulatorData))), c("darkblue", 
            "white", "darkred")), heatmap_legend_param = list(color_bar = "continuous", 
            legend_direction = "horizontal", title_gp = grid::gpar(fontsize = 8), 
            labels_gp = grid::gpar(fontsize = 6)), 
        width = grid::unit(Regwidth, "cm"))
    
    if (length(ClustModuleData) < 50) {
        fontsizeMo = 6
    } else if (length(ClustModuleData) < 200) {
        fontsizeMo = 4
    } else {
        fontsizeMo = 2
    }
    
    ha_Module <- ComplexHeatmap::Heatmap(ClustModuleData, 
        name = "", column_title = "Target Genes\nExpression", 
        cluster_rows = FALSE, cluster_columns = TRUE, 
        show_column_dend = FALSE, show_column_names = TRUE, 
        show_row_names = FALSE, column_names_gp = grid::gpar(fontsize = fontsizeMo), 
        show_heatmap_legend = FALSE, column_title_gp = grid::gpar(fontsize = 6, 
            fontface = "bold"), col = circlize::colorRamp2(c(-max(abs(ClustModuleData)), 
            0, max(abs(ClustModuleData))), c("darkblue", 
            "white", "darkred")), heatmap_legend_param = list(color_bar = "continuous", 
            legend_direction = "horizontal", title_gp = grid::gpar(fontsize = 8), 
            labels_gp = grid::gpar(fontsize = 6)))
    
    ha_list <- ha_Reg + ha_Module
    
    if (!is.null(MET_matrix)) {
        METreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$MET), 
            ModuleRegulators)
        print("MET regulators will be included when available")
        if (length(METreg) > 0) {
            MET_matrix = as.data.frame(MET_matrix)
            METData2 = METData = as.matrix(MET_matrix[unlist(Alterations %>% 
                dplyr::filter(MET == 1) %>% dplyr::select(HGNC_symbol)), 
                overlapping_samples_clust])
            METData2[which(METData > 0)] <- "Hyper-methylated"  # hyper
            METData2[which(METData < 0)] <- "Hypo-methylated"  # hypo
            METData2[which(METData == 0)] <- "Not altered"  # nothing
            METData2 <- t(METData2)
            Metwidth = ncol(METData2) * 0.5
            Met_col = structure(c("#006992", "#d95d39", 
                "white"), names = c("Hyper-methylated", 
                "Hypo-methylated", "Not altered"))
            ha_Met <- ComplexHeatmap::Heatmap(METData2, 
                name = "Methylation State", column_title = "Methylation State", 
                cluster_rows = FALSE, cluster_columns = TRUE, 
                show_column_dend = FALSE, show_column_names = TRUE, 
                show_row_names = FALSE, column_names_gp = grid::gpar(fontsize = 6), 
                show_heatmap_legend = TRUE, column_title_gp = grid::gpar(fontsize = 6, 
                  fontface = "bold"), col = Met_col, 
                width = grid::unit(Metwidth, "cm"), 
                heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8), 
                  labels_gp = grid::gpar(fontsize = 6)))
            ha_list <- ha_Met + ha_list
        }
    }
    
    if (!is.null(CNV_matrix)) {
        CNVreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$CNV), 
            ModuleRegulators)
        print("CNV regulators will be included when available")
        if (length(CNVreg) > 0) {
            CNV_matrix = as.data.frame(CNV_matrix)
            CNVData2 = CNVData = as.matrix(CNV_matrix[unlist(Alterations %>% 
                dplyr::filter(CNV == 1) %>% dplyr::select(HGNC_symbol)), 
                overlapping_samples_clust])
            CNVData2[which(CNVData >= 0.1)] <- "Amplified"  # amplification
            CNVData2[which(CNVData <= (-0.1))] <- "Deleted"  # deletion
            CNVData2[which(CNVData < 0.1 & CNVData > 
                (-0.1))] <- "Not altered"  # nothing
            CNVData2 <- t(CNVData2)
            CNVwidth = ncol(CNVData2) * 0.5
            CNV_col = structure(c("#006992", "#d95d39", 
                "white"), names = c("Deleted", "Amplified", 
                "Not altered"))
            ha_CNV <- ComplexHeatmap::Heatmap(CNVData2, 
                name = "CNV State", column_title = "CNV State", 
                cluster_rows = FALSE, cluster_columns = TRUE, 
                show_column_dend = FALSE, show_column_names = TRUE, 
                show_row_names = FALSE, column_names_gp = grid::gpar(fontsize = 6), 
                show_heatmap_legend = TRUE, column_title_gp = grid::gpar(fontsize = 6, 
                  fontface = "bold"), col = CNV_col, 
                width = grid::unit(CNVwidth, "cm"), 
                heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8), 
                  labels_gp = grid::gpar(fontsize = 6)))
            ha_list <- ha_CNV + ha_list
        }
    }
    
    if (!is.null(SAMPLE_annotation)) {
        if (ID %in% colnames(SAMPLE_annotation)) {
            SAMPLE_annotation_fil <- as.data.frame(SAMPLE_annotation) %>% 
                dplyr::filter(!!as.name(ID) %in% overlapping_samples_clust)
            suppressWarnings(SAMPLE_annotation_fil <- dplyr::left_join(as.data.frame(overlapping_samples_clust), 
                SAMPLE_annotation_fil, by = c(overlapping_samples_clust = ID)))
            SAMPLE_annotation_fil <- tibble::column_to_rownames(SAMPLE_annotation_fil, 
                "overlapping_samples_clust")
            cat(nrow(SAMPLE_annotation_fil), "samples have an annotation.\n")
            cat(ncol(SAMPLE_annotation_fil), "annotations are added")
            # define colors
            col <- c()
            for (sample_column in colnames(SAMPLE_annotation_fil)[colnames(SAMPLE_annotation_fil) != 
                ID]) {
                newcol <- circlize::rand_color(n = length(unique(SAMPLE_annotation_fil[, 
                  sample_column])), luminosity = "bright")
                names(newcol) <- unique(SAMPLE_annotation_fil[, 
                  sample_column])
                col <- c(col, newcol)
            }
            ha_anot <- ComplexHeatmap::Heatmap(SAMPLE_annotation_fil, 
                name = "Sample Annotation", column_title = "Sample\nAnnotation", 
                column_title_gp = grid::gpar(fontsize = 6, 
                  fontface = "bold"), col = col, show_row_names = FALSE, 
                width = grid::unit(4, "mm"), column_names_gp = grid::gpar(fontsize = 6), 
                heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8), 
                  labels_gp = grid::gpar(fontsize = 6)))
            ha_list <- ha_list + ha_anot
        } else {
            print("The ID is not identified as a column name in the annotation")
        }
    }
    if (printHM == TRUE) {
        ComplexHeatmap::draw(ha_list, heatmap_legend_side = "bottom", 
            annotation_legend_side = "bottom")
    } else {
        return(ha_list)
    }
}
