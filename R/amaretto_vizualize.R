#' AMARETTO_VisualizeModule
#'
#' Function to visualize the gene modules
#'
#' @param AMARETTOinit List output from AMARETTO_Initialize().
#' @param AMARETTOresults List output from AMARETTO_Run().
#' @param ProcessedData List of processed input data
#' @param ModuleNr Module number to visualize
#' @param SAMPLE_annotation Matrix or Dataframe with sample annotation
#' @param ID Column used as sample name
#' @param show_row_names 
#' @param order_samples Order samples in heatmap by mean or by clustering
#' @param printHM Boolean to print heatmap
#'
#' @importFrom circlize colorRamp2  rand_color
#' @importFrom grid gpar unit
#' @importFrom stats dist hclust
#' @importFrom dplyr  left_join mutate  select  summarise  rename  filter case_when
#' @import grDevices
#' @import methods
#' @importFrom ComplexHeatmap  HeatmapAnnotation Heatmap draw
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
AMARETTO_VisualizeModule <- function(AMARETTOinit,AMARETTOresults,ProcessedData,ModuleNr,show_row_names=FALSE, SAMPLE_annotation=NULL,ID=NULL,order_samples=NULL,printHM=TRUE) {
  CNV_matrix <- ProcessedData[[2]]
  MET_matrix <- ProcessedData[[3]]
  CNVMet_Alterations <- DriversList_Alterations <- MET <- HGNC_symbol <- CNV <- NULL
  if (ModuleNr>AMARETTOresults$NrModules){
    stop('\tCannot plot Module',ModuleNr,' since the total number of modules is',AMARETTOresults$N,'.\n')
  }
  ModuleData<-as.data.frame(AMARETTOinit$MA_matrix_Var)[AMARETTOresults$ModuleMembership==ModuleNr,]
  ModuleRegulators <- AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
  RegulatorData <- as.data.frame(AMARETTOinit$RegulatorData)[ModuleRegulators,]
  ModuleGenes <- rownames(ModuleData)
  cat('Module',ModuleNr,'has',length(rownames(ModuleData)),'genes and',length(ModuleRegulators),'regulators for',length(colnames(ModuleData)),'samples.\n')
  Alterations<- tibble::rownames_to_column(as.data.frame(AMARETTOinit$RegulatorAlterations$Summary),"HGNC_symbol") %>% dplyr::rename(DriverList="Driver List") %>% dplyr::filter(HGNC_symbol %in% ModuleRegulators)
  Alterations<- Alterations %>% dplyr::mutate(CNVMet_Alterations=case_when(MET==1 & CNV==1~"Methylation and copy number alterations",
                                                                           CNV==1~"Copy number alterations",
                                                                           MET==1~"Methylation aberrations",
                                                                           MET==0 & CNV==0 ~"Not Altered"),
                                              DriversList_Alterations=case_when(DriverList==0~"Driver not predefined",
                                                                                DriverList==1~"Driver predefined"))

  ha_drivers <- ComplexHeatmap::HeatmapAnnotation(df =tibble::column_to_rownames(Alterations,"HGNC_symbol") %>% dplyr::select(CNVMet_Alterations,DriversList_Alterations) , 
                                                  col = list(CNVMet_Alterations= c("Copy number alterations"="#eca400","Methylation aberrations"="#006992","Methylation and copy number alterations"="#d95d39","Not Altered"="lightgray"),
                                                             DriversList_Alterations=c("Driver not predefined"="lightgray","Driver predefined"="#588B5B")),
                                                  which = "column", height = grid::unit(0.3, "cm"),name="",show_annotation_name = FALSE,
                                  annotation_legend_param = list(title_gp = grid::gpar(fontsize = 8),labels_gp = grid::gpar(fontsize = 6)))

  overlapping_samples <- colnames(ModuleData)
  
  # Add NAs for Gene Expression samples not existing in CNV or MET data 

  if (!is.null(CNV_matrix)){
    non_CNV_sample_names<- overlapping_samples[!overlapping_samples%in%colnames(CNV_matrix)]
    print(non_CNV_sample_names)
    print(nrow(CNV_matrix))
    print(length(non_CNV_sample_names))
    CNV_extension_mat<- matrix(data=NA,nrow=nrow(CNV_matrix),ncol=length(non_CNV_sample_names))
    colnames(CNV_extension_mat)<-non_CNV_sample_names
    CNV_matrix<-cbind(CNV_matrix,CNV_extension_mat)
  }
  if (!is.null(MET_matrix)){
    non_MET_sample_names<- overlapping_samples[!overlapping_samples%in%colnames(MET_matrix)]
    MET_extension_mat<- matrix(data=NA,nrow=nrow(MET_matrix),ncol=length(non_MET_sample_names))
    colnames(MET_extension_mat)<-non_MET_sample_names
    MET_matrix<-cbind(MET_matrix,MET_extension_mat)
  }

  if(is.null(order_samples)){
    overlapping_samples_clust<-overlapping_samples[order(colMeans(ModuleData[,overlapping_samples]))]
  }else if(order_samples=="clust"){
    SampleClustering<-stats::hclust(stats::dist(t(ModuleData[,overlapping_samples])), method = "complete", members = NULL)
    overlapping_samples_clust<-overlapping_samples[SampleClustering$order]
  }else {
    print("ordering type not recognized, samples will be orderd based on mean expression of the module genes")
    overlapping_samples_clust<-overlapping_samples[order(colMeans(ModuleData[,overlapping_samples]))]
  }

  ClustRegulatorData <- t(RegulatorData[,overlapping_samples_clust])
  ClustModuleData <- t(ModuleData[,overlapping_samples_clust])
  
  if(length(ClustModuleData)<50){
    fontsizeMo=6
  } else if (length(ClustModuleData)<200){
    fontsizeMo=4
  } else {fontsizeMo=2}
  
  Regwidth <- ncol(ClustRegulatorData)*0.5
  ha_Reg <- Heatmap(ClustRegulatorData, name = "Gene Expression", column_title = "Regulator Genes\nExpression",cluster_rows=FALSE,cluster_columns=TRUE,show_column_dend=FALSE,show_column_names=TRUE,show_row_names=FALSE,column_names_gp = grid::gpar(fontsize = fontsizeMo),top_annotation = ha_drivers,
                    column_title_gp = grid::gpar(fontsize = 6, fontface = "bold"), col=circlize::colorRamp2(c(-max(abs(ClustRegulatorData)), 0, max(abs(ClustRegulatorData))), c("darkblue", "white", "darkred")),heatmap_legend_param = list(color_bar = "continuous",legend_direction = "horizontal",title_gp = grid::gpar(fontsize = 8),labels_gp = grid::gpar(fontsize = 6)), width = grid::unit(Regwidth, "cm"))


  ha_Module <- Heatmap(ClustModuleData, name = "", column_title = "Target Genes\nExpression",cluster_rows=FALSE,cluster_columns=TRUE,show_column_dend=FALSE,show_column_names=TRUE,show_row_names=show_row_names,column_names_gp = grid::gpar(fontsize = fontsizeMo),show_heatmap_legend = FALSE,
                       column_title_gp = grid::gpar(fontsize = 6, fontface = "bold"), col=circlize::colorRamp2(c(-max(abs(ClustModuleData)), 0, max(abs(ClustModuleData))), c("darkblue", "white", "darkred")),heatmap_legend_param = list(color_bar = "continuous",legend_direction = "horizontal",title_gp = grid::gpar(fontsize = 8),labels_gp = grid::gpar(fontsize = 6)))

  ha_list<- ha_Reg + ha_Module
  if (!is.null(MET_matrix)){
    METreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$MET),ModuleRegulators)
    print("MET regulators will be included when available")
    if (length(METreg)>0){
      MET_matrix = as.data.frame(MET_matrix)
      METData2 = METData = as.matrix(MET_matrix[unlist(Alterations %>% dplyr::filter(MET==1) %>% dplyr::select(HGNC_symbol)),overlapping_samples_clust])
      METData2[which(METData>0)] <- "Hyper-methylated"  # hyper
      METData2[which(METData<0)] <- "Hypo-methylated"  # hypo
      METData2[which(METData==0)] <- "Not altered" # nothing
      METData2<-t(METData2)
      Metwidth=ncol(METData2)*0.7
      Met_col=structure(c("#006992","#d95d39","white"),names=c("Hyper-methylated","Hypo-methylated","Not altered"))
      ha_Met <- Heatmap(METData2, name = "Methylation State", column_title = "Methylation", cluster_rows=FALSE,cluster_columns=TRUE,show_column_dend=FALSE,show_column_names=TRUE,show_row_names=FALSE,column_names_gp = grid::gpar(fontsize = 6),show_heatmap_legend = TRUE,
                        column_title_gp = gpar(fontsize = 6, fontface = "bold"), col = Met_col, width = grid::unit(Metwidth, "cm"),heatmap_legend_param = list(title_gp = gpar(fontsize = 8),labels_gp = grid::gpar(fontsize = 6)))
      ha_list<- ha_Met + ha_list
    }
  }

  if (!is.null(CNV_matrix)){
    CNVreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$CNV),ModuleRegulators)
    print("CNV regulators will be included when available")
    if (length(CNVreg)>0){
      CNV_matrix = as.data.frame(CNV_matrix)
      CNVData2 = CNVData = as.matrix(CNV_matrix[unlist(Alterations %>% dplyr::filter(CNV==1) %>% dplyr::select(HGNC_symbol)),overlapping_samples_clust])
      CNVData2[which(CNVData>=0.1)] <- "Amplified"  # amplification
      CNVData2[which(CNVData<=(-0.1))] <- "Deleted"  # deletion
      CNVData2[which(CNVData<0.1 & CNVData>(-0.1))] <- "Not altered" # nothing
      CNVData2<-t(CNVData2)
      CNVwidth=ncol(CNVData2)*0.7
      CNV_col=structure(c("#006992","#d95d39","white"),names=c("Deleted","Amplified","Not altered"))
      ha_CNV <- Heatmap(CNVData2, name = "CNV State", column_title = "CNV", cluster_rows=FALSE,cluster_columns=TRUE,show_column_dend=FALSE,show_column_names=TRUE,show_row_names=FALSE,column_names_gp = grid::gpar(fontsize = 6),show_heatmap_legend = TRUE,
                        column_title_gp = grid::gpar(fontsize = 6, fontface = "bold"),col = CNV_col,width = grid::unit(CNVwidth, "cm"),heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8),labels_gp = grid::gpar(fontsize = 6)))
      ha_list<-ha_CNV + ha_list
    }
  }

  if (!is.null(SAMPLE_annotation)){
    if (ID %in% colnames(SAMPLE_annotation)){
      SAMPLE_annotation_fil<-as.data.frame(SAMPLE_annotation) %>% dplyr::filter(!!as.name(ID) %in% overlapping_samples_clust)
      suppressWarnings(SAMPLE_annotation_fil<-dplyr::left_join(as.data.frame(overlapping_samples_clust),SAMPLE_annotation_fil,by=c("overlapping_samples_clust"=ID)))
      SAMPLE_annotation_fil<-tibble::column_to_rownames(SAMPLE_annotation_fil,"overlapping_samples_clust")
      cat(nrow(SAMPLE_annotation_fil),"samples have an annotation.\n")
      cat(ncol(SAMPLE_annotation_fil),"annotations are added")
      #define colors
      col<-c()
      # for (sample_column in colnames(SAMPLE_annotation_fil)[colnames(SAMPLE_annotation_fil) != ID]){
      #   newcol<-circlize::rand_color(n=length(unique(SAMPLE_annotation_fil[,sample_column])),luminosity = "bright")
      #   names(newcol)<-unique(SAMPLE_annotation_fil[,sample_column])
      #   col<-c(col,newcol)
      # }
      fsize<-6
      wsize=2
      if (ncol(SAMPLE_annotation_fil)>50){
        fsize<-5
        wsize=1
      }
      
      for (sample_column in colnames(SAMPLE_annotation_fil)[colnames(SAMPLE_annotation_fil) != ID]){
        # newcol<-circlize::rand_color(n=length(unique(SAMPLE_annotation_fil[,colnames(SAMPLE_annotation_fil)==sample_column])),luminosity = "bright")
        # names(newcol)<-unique(SAMPLE_annotation_fil[,colnames(SAMPLE_annotation_fil)==sample_column])
        annotation_data<-SAMPLE_annotation_fil[,colnames(SAMPLE_annotation_fil)==sample_column]
        unique_annotations<-as.vector(sort(unique(annotation_data),decreasing = FALSE))
        #unique_annotations <- unique_annotations[!is.na(unique_annotations)]
        if (length(unique_annotations)<6 & length(unique_annotations) > 1){
          
          annotation_data<-as.factor(annotation_data)
          if (length(unique_annotations) == 2){
            colors<- c("darkblue", "darkred")
          }else if (length(unique_annotations) == 3){
            colors<- c("darkblue", "darkgreen", "darkred")
          }else if (length(unique_annotations) == 4){
            colors<- c("darkblue", "darkgreen", "pink2", "darkred")
          }else {
            colors<- c("darkblue", "darkgreen", "pink2","yellow3", "darkred")
          }
          print("Hi")
          print(colors)
          print(unique_annotations)
          
          col_list<-structure(colors,names=unique_annotations)
          print(col_list)
          print("bye")
          # names(colors)<-unique_annotations
          # print(colors)
          # print(unique_annotations)
          # colormap<-circlize::colorRamp2(unique_annotations,colors)
          # print(colormap)
          
          # col_list = list(c("Copy number alterations"="#eca400","Methylation aberrations"="#006992","Methylation and copy number alterations"="#d95d39","Not Altered"="lightgray"),
          #            DriversList_Alterations=c("Driver not predefined"="lightgray","Driver predefined"="#588B5B"))
          # 
          #Met_col=structure(c("#006992","#d95d39","white"),names=c("Hyper-methylated","Hypo-methylated","Not altered"))
          # col_list <-list(colors)
          # names(col_list)<-sample_column
          # 
          # print(col_list)
          # col=col_list,
          ha_anot<-Heatmap(annotation_data, name=sample_column, column_title ="", column_title_gp = grid::gpar(fontsize = 5, fontface = "bold"),col=col_list, show_row_names=FALSE,width = unit(wsize, "mm"),
                           column_names_gp = gpar(fontsize = fsize),heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 6),labels_gp = grid::gpar(fontsize = 6),ncol = 1))
        }
        else{
          ha_anot<-Heatmap(annotation_data, name=sample_column, column_title ="", column_title_gp = grid::gpar(fontsize = 5, fontface = "bold"),  show_row_names=FALSE,width = unit(wsize, "mm"),
                           column_names_gp = gpar(fontsize = fsize),heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 6),labels_gp = grid::gpar(fontsize = 6),ncol = 1))
        }
        # ha_anot<-Heatmap(SAMPLE_annotation_fil[,colnames(SAMPLE_annotation_fil)==sample_column], name="Sample Annotation", column_title = "Sample\nAnnotation", column_title_gp = grid::gpar(fontsize = 5, fontface = "bold"), col=col[,k], show_row_names=FALSE,width = unit(wsize, "mm"),
        #                  column_names_gp = gpar(fontsize = fsize),heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 6),labels_gp = grid::gpar(fontsize = 6),ncol = 4))
        ha_list<-ha_list + ha_anot
      }
      
      
    } else {print("The ID is not identified as a column name in the annotation")}
  }
  if (printHM==TRUE){
  ComplexHeatmap::draw(ha_list,heatmap_legend_side = "bottom",annotation_legend_side="bottom")
    
  } else {
    return(ha_list)
  }
}
