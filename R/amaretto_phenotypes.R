#' @title AMARETTO_PhenoAssociation
#'
#' @param AMARETTOresults
#' @param annotation
#' @param idpatients
#' @param printplots
#' @param pdfname
#'
#' @return 
#' 
#' @importFrom dplyr arrange group_by inner_join mutate select summarise rename
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr drop_na
#' @import ggplot2
#' 
#' @return
#' @export
AMARETTO_PhenAssociation <- function(AMARETTOresults, annotation, idpatients, phenotypetests, printplots=FALSE, pdfname="phenotypes.pdf"){
  
  mean_expression_modules <- t(AMARETTOresults$ModuleData)
  mean_expression_modules <- tibble::rownames_to_column(as.data.frame(mean_expression_modules),"idpatients")
  ### test matching sample ids
  if (sum(annotation[,idpatients] %in% mean_expression_modules[,"idpatients"],na.rm=TRUE)==0){
    stop("No overlap between patients ids.")
  } else {
    print(paste0("Phenotypic association will be calculated on ",length(annotation[,idpatients] %in% mean_expression_modules[,"idpatients"]), " patients."))
  }
  annotation <- suppressMessages(dplyr::inner_join(annotation, mean_expression_modules %>% dplyr::rename(!!idpatients :="idpatients")))
  phenotypetests<-as.data.frame(phenotypetests)
  ### select phenotypes for which a ligid test is given
  result_df <- data.frame(ModuleNr=rep(paste0("Module_",1:AMARETTOresults$NrModules)))
  if (printplots==TRUE){pdf(pdfname)}
  
  for(parameter in phenotypetests[,"Phenotypes"]){
    sample_size <- nrow(annotation %>% dplyr::select(!!parameter) %>% drop_na())
    
    print(paste0("Phenotypic association is calculated for ",parameter, " on ",sample_size," patients."))
    test <- phenotypetests[which(phenotypetests$Phenotypes==parameter),"test"]
    
    if (!test %in% c("WILCOXONRANKSUMTEST","KRUSKALWALLISTEST","TTEST","ANOVATEST","PEARSONCORRTEST","SPEARMANCORRTEST","cathegorical","ordinal","continuous")){
      stop("There are no tests that match the options.")
    }
    
    if ((test == "cathegorical" && sample_size<30) || (test %in% c("WILCOXONRANKSUMTEST","KRUSKALWALLISTEST"))){
      annotation[,parameter] <- as.factor(annotation[,parameter])
      if ((nlevels(annotation[,parameter])==2) || test == "WILCOXONRANKSUMTEST"){
          print(paste0("A wilcox test is performed for ",parameter))
          colnames_result_df <- rep(paste0(parameter,c("_Wilcoxon_p","_Wilcoxon_padj","_Wilcoxon_95LI","_Wilcoxon_95HI")))
          result_df[,colnames_result_df] <- NA
          for(i in 1:AMARETTOresults$NrModules){
            moduleNr <- paste0("Module_",i)
            testresults <- wilcox.test(annotation[,moduleNr]~annotation[,parameter], conf.int = TRUE)
            
            result_df[i,colnames_result_df[1]]<-testresults$p.value
            result_df[i,colnames_result_df[3]]<-testresults$conf.int[1]
            result_df[i,colnames_result_df[4]]<-testresults$conf.int[2]
            
            if(printplots==TRUE){
              print(ggplot(annotation %>% tidyr::drop_na(!!parameter),aes(x=get(parameter), y=get(moduleNr), fill=get(parameter)))+
                geom_boxplot()+
                geom_jitter(color="lightgray")+
                theme_classic()+
                theme(legend.position = "none")+
                labs(x=parameter,y=moduleNr,caption=paste0("p=",round(testresults$p.value,4))))
            }
          }
         result_df[,colnames_result_df[2]] <- p.adjust(result_df[,colnames_result_df[1]],method = "BH")
      } else if ((nlevels(annotation[,parameter])>2) || test == "KRUSKALWALLISTEST"){
          print(paste0("A Kruskal-Wallis Rank sum test is performed for ",parameter))
          colnames_result_df <- rep(paste0(parameter,c("_KrusW_p","_KrusW_padj","_KrusW_stat")))
          result_df[,colnames_result_df] <- NA
          for(i in 1:AMARETTOresults$NrModules){
            moduleNr <- paste0("Module_",i)
              testresults<-kruskal.test(annotation[,i]~annotation[,parameter])
              result_df[i,colnames_result_df[1]]<-testresults$p.value
              result_df[i,colnames_result_df[3]]<-testresults$statistic
              if(printplots==TRUE){
                print(ggplot(annotation %>% drop_na(!!parameter),aes(x=get(parameter), y=get(moduleNr), fill=get(parameter)))+
                        geom_boxplot()+
                        geom_jitter(color="lightgray")+
                        theme_classic()+
                        theme(legend.position = "none")+
                        labs(x=parameter,y=moduleNr,caption=paste0("p=",round(testresults$p.value,4))))
              }
            }
          result_df[,colnames_result_df[2]] <- p.adjust(result_df[,colnames_result_df[1]],method = "BH")
      } else if (nlevels(annotation[,parameter])<2){
          stop(paste0(parameter, " has only one or no levels"))
        }
    } else if ((test == "cathegorical" && sample_size>=30) || (test %in% c("TTEST","ANOVATEST"))){
      annotation[,parameter]<-as.factor(annotation[,parameter])
      if ((nlevels(annotation[,parameter])==2) || test == "TTEST"){
        print(paste0("A t-test is performed for ",parameter))
        colnames_result_df<-rep(paste0(parameter,c("_Ttest_p","_Ttest_padj","_Ttest_95LI","_Ttest_95HI")))
        result_df[,colnames_result_df]<-NA
        for(i in 1:AMARETTOresults$NrModules){
          moduleNr <- paste0("Module_",i)
          testresults<-t.test(annotation[,moduleNr]~annotation[,parameter])
          #return results
          result_df[i,colnames_result_df[1]]<-testresults$p.value
          result_df[i,colnames_result_df[3]]<-testresults$conf.int[1]
          result_df[i,colnames_result_df[4]]<-testresults$conf.int[2]
          if(printplots==TRUE){
            print(ggplot(annotation %>% drop_na(!!parameter),aes(x=get(parameter), y=get(moduleNr), fill=get(parameter)))+
                    geom_boxplot()+
                    theme_classic()+
                    theme(legend.position = "none")+
                    labs(x=parameter,y=moduleNr,caption=paste0("p=",round(testresults$p.value,4))))
          }
        }
        result_df[,colnames_result_df[2]]<-p.adjust(result_df[,colnames_result_df[1]],method = "BH")
      } else if ((nlevels(annotation[,parameter])>2) || test =="ANOVATEST"){
        print(paste0("An ANOVA test is performed for ",parameter, ". The levels are: ",paste0(levels(annotation[,parameter]),collapse=", "),"."))
        colnames_result_df<-rep(paste0(parameter,c("_Anova_p","_Anova_padj")))
        result_df[,colnames_result_df]<-NA
        for(i in 1:AMARETTOresults$NrModules){
          moduleNr <- paste0("Module_",i)
          lmod<-lm(annotation[,moduleNr]~annotation[,parameter])
          testresults<-anova(lmod)
          #returns results
          result_df[i,colnames_result_df[1]]<-as.numeric(testresults[1,5])

          if(printplots==TRUE){
            print(ggplot(annotation %>% drop_na(!!parameter),aes(x=get(parameter), y=get(moduleNr), fill=get(parameter)))+
                    geom_boxplot()+
                    theme_classic()+
                    theme(legend.position = "none")+
                    labs(x=parameter,y=moduleNr,caption=paste0("p=",round(as.numeric(testresults[1,5])))))
          }
        }
        result_df[,colnames_result_df[2]]<-p.adjust(result_df[,colnames_result_df[1]],method = "BH")
      } else if (nlevels(annotation[,parameter])<2){
        stop(paste0(parameter, " has only one or no levels"))
      }
    } else if ((test %in% c("ordinal","SPEARMANCORRTEST")) || (test == "continuous" && sample_size<30)){
      if (test == "continuous" || test =="SPEARMANCORRTEST"){
        annotation[,parameter] <- as.numeric(annotation[,parameter])
        print(paste0("A Spearman Correlation is performed for ",parameter))
      } else {
        annotation[,parameter]<-as.factor(annotation[,parameter])
        print(paste0("A spearman correlation is calculated. The order of the factor levels is: ",paste0(levels(annotation[,parameter]),collapse=", ")))
        annotation[,parameter]<-as.numeric(annotation[,parameter])
      }
      colnames_result_df<-rep(paste0(parameter,c("_Spearman_p","_Spearman_padj","_Spearmanstatistic")))
      result_df[,colnames_result_df]<-NA
      for(i in 1:AMARETTOresults$NrModules){
        moduleNr <- paste0("Module_",i)
        testresults <- suppressWarnings(cor.test(annotation[,moduleNr],as.numeric(annotation[,parameter]), method="spearman", use="complete.obs"))
        
        result_df[i,colnames_result_df[1]]<-testresults$p.value
        result_df[i,colnames_result_df[3]]<-testresults$statistic
        if(printplots==TRUE){
          print(ggplot(annotation %>% drop_na(!!parameter),aes(x=get(parameter), y=get(moduleNr), fill=get(parameter)))+
                  geom_point()+
                  theme_classic()+
                  theme(legend.position = "none")+
                  labs(x=parameter,y=moduleNr,caption=paste0("p=",round(testresults$p.value,4))))
        }  
      }
      result_df[,colnames_result_df[2]]<-p.adjust(result_df[,colnames_result_df[1]],method = "BH")
    } else if ((test == "continuous" && sample_size>=30) || test == "PEARSONCORRTEST"){
      print(paste0("A Pearson Correlation is performed for ",parameter))
      colnames_result_df<-rep(paste0(parameter,c("_Pearson_p","_Pearson_padj","_Pearson95LI","_Pearson95HI")))
      result_df[,colnames_result_df]<-NA
      for(i in 1:AMARETTOresults$NrModules){
        moduleNr <- paste0("Module_",i)
        annotation[,parameter]<-as.numeric(annotation[,parameter])
        
        testresults <- suppressWarnings(cor.test(annotation[,moduleNr],as.numeric(annotation[,parameter]), method="pearson", use="complete.obs"))
        
        result_df[i,colnames_result_df[1]]<-testresults$p.value
        result_df[i,colnames_result_df[3]]<-testresults$conf.int[1]
        result_df[i,colnames_result_df[4]]<-testresults$conf.int[2]
        if(printplots==TRUE){
          print(ggplot(annotation %>% drop_na(!!parameter),aes(x=get(parameter), y=get(moduleNr), fill=get(parameter)))+
                  geom_point()+
                  geom_smooth(method="lm")+
                  theme_classic()+
                  theme(legend.position = "none")+
                  labs(x=parameter,y=moduleNr,caption=paste0("p=",round(testresults$p.value,4))))
        }
      }
      result_df[,colnames_result_df[2]]<-p.adjust(result_df[,colnames_result_df[1]],method = "BH")
    }
  }
  if (printplots==TRUE){dev.off()}
  
  return(result_df)
}

#' @title AMARETTO_phenotypesshiny
#'
#' @param AMARETTOresults
#' @param annotation
#' 
#' @return A set of HTMLs, giving caracteristics of the communities
#' 
#' @import tidyverse
#' @import shinythemes
#' @import shiny
#' @import DT
#' 
#' @return
#' @export
AMARETTO_phenotypesshiny<- function(AMARETTOresults, annotation, idpatients){
  
  
  mean_expression_modules <- t(AMARETTOresults$ModuleData)
  mean_expression_modules <- rownames_to_column(as.data.frame(mean_expression_modules),"idpatients")
  
  parameters <- colnames(annotation)
  parameters <- parameters[!parameters %in% idpatients]
  if (sum(annotation[,idpatients] %in% mean_expression_modules[,"idpatients"],na.rm=TRUE)==0){
    stop("No overlap between patients ids")
  } else {
    print(paste0("Phenotypic association will be calculated on ",length(annotation[,idpatients] %in% mean_expression_modules[,"idpatients"]), " patients."))
  }
  
  annotation <- suppressMessages(dplyr::inner_join(annotation, mean_expression_modules %>% dplyr::rename(!!idpatients :="idpatients")))
  
  modules <- paste0("Module ", seq(1:AMARETTOresults$NrModules))
  
  ui<-fluidPage(
    theme = shinytheme("paper"),
    titlePanel("Phenotype Analysis of the Modules"),
    sidebarPanel(
      selectInput(inputId = "Module", 
                  label="Which module do you want to analyse:", 
                  choices = modules, 
                  selected = modules[1]),
      
      selectInput(inputId = "Parameter", 
                  label="Which parameter do you want to analyse:", 
                  choices = parameters, 
                  selected = parameters[1]),
      
      selectInput(inputId = "TypeTest",
                  label = "What is the type of the parameter:",
                  choices = c("cathegorical","ordinal","continuous"),
                  selected = "cathegorical"),
      
      submitButton(text="Calculate")
      
    ),
    mainPanel(
      h3("Analysis Plot"),
      plotOutput(outputId = "AnalysisPlot",height = "600px"),
      br(),
      h3("Analysis Results"),
      verbatimTextOutput(outputId="StatisticsAnalysis"),
      dataTableOutput(outputId = "DTAnalysis")
    )
  )
  
  server<-function(input,output){
    
    testresults <- reactive({
      parameter<-input$Parameter
      test<-input$TypeTest
      moduleNr<-input$Module
      moduleNr <-sub(" ","_",moduleNr)
      sample_size<-nrow(annotation %>% dplyr::select(!!parameter) %>% drop_na())
      if (test == "cathegorical" && sample_size<30 ){
        annotation[,parameter]<-as.factor(annotation[,parameter])
        if (nlevels(annotation[,parameter])==2){
          #print(paste0("A wilcox test is performed for ",parameter))
          return(wilcox.test(annotation[,moduleNr]~annotation[,parameter], conf.int = TRUE))
        } else if (nlevels(annotation[,parameter])>2){
          #print(paste0("A Kruskal-Wallis Rank sum test is performed for ",parameter))
          return(kruskal.test(annotation[,i]~annotation[,parameter]))
        } else {
          return(paste0(parameter, " has less than two levels."))
        }
      } else if (test == "cathegorical" && sample_size>=30){
        annotation[,parameter]<-as.factor(annotation[,parameter])
        if (nlevels(annotation[,parameter])==2){
          #print(paste0("A t-test is performed for ",parameter))
          return(t.test(annotation[,moduleNr]~annotation[,parameter]))
        } else if (nlevels(annotation[,parameter])>2){
          #print(paste0("An ANOVA test is performed for ",parameter, ". The levels are: ",paste0(levels(annotation[,parameter]),collapse=", "),"."))
          lmod<-lm(annotation[,moduleNr]~annotation[,parameter])
          return(anova(lmod))
        } else {
          stop(paste0(parameter, " has only one level."))
        }
      } else if ( test=="ordinal" || (test == "continuous" && sample_size<30)){
        if (test == "continuous"){
          annotation[,parameter] <- as.numeric(annotation[,parameter])
          #print(paste0("A Spearman Correlation is performed for ",parameter))
        } else {
          annotation[,parameter]<-as.factor(annotation[,parameter])
          #print(paste0("A spearman correlation is calculated. The order of the factor levels is: ",paste0(levels(annotation[,parameter]),collapse=", ")))
          annotation[,parameter]<-as.numeric(annotation[,parameter])
        }
        return(suppressWarnings(cor.test(annotation[,moduleNr],as.numeric(annotation[,parameter]), method="spearman", use="complete.obs")))
      } else if ( test == "continuous" && sample_size>=30){
        #print(paste0("A Pearson Correlation is performed for ",parameter))
        return(suppressWarnings(cor.test(annotation[,moduleNr],as.numeric(annotation[,parameter]), method="pearson", use="complete.obs")))
      }
  })
    
    output$StatisticsAnalysis <- renderPrint({
      print(testresults())
    })
    
    output$AnalysisPlot<-renderPlot({
      moduleNr<-input$Module
      moduleNr <-sub(" ","_",moduleNr)
      parameter<-input$Parameter
      if(input$TypeTest!="continuous"){
      print(ggplot(annotation %>% drop_na(!!parameter),aes(x=get(parameter), y=get(moduleNr), fill=get(parameter)))+
              geom_boxplot()+
              geom_jitter(color="lightgray")+
              theme_classic()+
              theme(legend.position = "none")+
              labs(x=parameter,y=moduleNr))
      } else {
        print(ggplot(annotation %>% drop_na(!!parameter),aes(x=as.numeric(get(parameter)), y=get(moduleNr)))+
          geom_point(color="lightgray")+
          theme_classic()+
          theme(legend.position = "none")+
          labs(x=parameter,y=moduleNr))
      }
    })
    
    output$DTAnalysis <- renderDataTable({
      moduleNr<-input$Module
      moduleNr <-sub(" ","_",moduleNr)
      parameter<-input$Parameter
      if(input$TypeTest!="continuous"){
        data_summary<-annotation %>% dplyr::select(!!idpatients,!!parameter,!!moduleNr) %>% drop_na(!!parameter) %>% dplyr::group_by_at(vars(!!parameter)) %>% dplyr::summarise(mean=mean(get(!!moduleNr)),median=median(get(!!moduleNr)),nobservations=n(),sd=sd(get(!!moduleNr)))
        datatable_summary<-datatable(data_summary) %>% formatRound(c('mean','median','sd'),digits=3)
      } else {
        data_summary<-annotation %>% dplyr::select(!!idpatients,!!parameter,!!moduleNr) %>% drop_na(!!parameter) %>% dplyr::summarise(median=median(get(!!moduleNr)),sd=sd(get(!!moduleNr)),parameter_min=min(as.numeric(get(!!parameter)),na.rm=TRUE),parameter_max=max(as.numeric(get(!!parameter)),na.rm=TRUE),nobservations=n())
        datatable_summary<-datatable(data_summary) %>% formatRound(c('median','sd'),digits=3)
      }
      datatable_summary
    })
  }
  
  app<-shinyApp(ui=ui,server=server)
  runApp(app, port = 8888)
}
