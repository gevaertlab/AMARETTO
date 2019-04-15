#' @title AMARETTO_survival
#'
#' @param AMARETTOresults
#' @param SurvivalAnnotation
#' @param idpatients
#' @param idvitalstatus
#' @param iddaystodeath
#' @param idfollowupdays
#' @param daystoyears
#' @param rightcensored
#' @param printsurv
#'
#' @return A set of HTMLs, giving caracteristics of the communities
#' 
#' @import survminer
#' @import tidyverse
#' @import survival
#' @import gtools
#' @return
#' @export
AMARETTO_survival<- function(AMARETTOresults, SurvivalAnnotation, idpatients, idvitalstatus, iddaystodeath, idfollowupdays, daystoyears, rightcensored=FALSE, printHR=FALSE, printsurv=FALSE){
  
  mean_expression_modules <- t(AMARETTOresults$ModuleData)
  mean_expression_modules <- rownames_to_column(as.data.frame(mean_expression_modules),"idpatients")
  if (sum(SurvivalAnnotation[,idpatients] %in% mean_expression_modules[,"idpatients"],na.rm=TRUE)==0){
    stop("No overlap between patients ids")
  } else {
    print(paste0("Survival data will be calculated on ",length(SurvivalAnnotation[,idpatients] %in% mean_expression_modules[,"idpatients"]), " patients."))
  }
  
  SurvivalAnnotation <- suppressMessages(dplyr::inner_join(SurvivalAnnotation, mean_expression_modules %>% dplyr::rename(!!idpatients :="idpatients")))

  daysdeath <- as.numeric(SurvivalAnnotation[grep("dead", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE),iddaystodeath])
  dead <- SurvivalAnnotation[grep("dead", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE),idvitalstatus]
  daysfollowup <- as.numeric(SurvivalAnnotation[grep("alive", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE),idfollowupdays])
  alive <- SurvivalAnnotation[grep("alive", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE),idvitalstatus]
  time<- c(daysfollowup,daysdeath)
  status <- tolower(c(alive, dead))
  
  if (daystoyears == TRUE){
    time<-time/365
  }
  if (!rightcensored == FALSE){
    if (is.numeric(rightcensored)){
      status[time>rightcensored]<-"alive"
      time[time>rightcensored]<-rightcensored
    } else {
      stop("The time point to censor the data is incorrect.")
    }
  }
  
  status <- as.numeric(as.factor(status))
  c_surv<-Surv(time, status)

  #per module
  p_survival<-matrix(nrow=AMARETTOresults$NrModules, ncol=13)
  colnames(p_survival)<-c("coxregwaldtestp","coxregwaldtestpadj","coxreg_coef","coxreg_LL","coxreg_UL","logranktestp","logranktestpadj","medianOS1","medianOS2","Lower95CI1","Upper95CI1","Lower95CI2","Upper95CI2")
  for(i in 1:AMARETTOresults$NrModules){
    moduleNr <- paste0("Module_",i)
    mean_dead <- SurvivalAnnotation[grep("dead", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE), moduleNr]
    mean_alive <- SurvivalAnnotation[grep("alive", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE), moduleNr]
    mean_expression <- c(mean_alive,mean_dead)
    mean_expression_cut <- quantcut(mean_expression, q=2, na.rm=TRUE)
    datasurvfit <- as.data.frame(cbind(time,status,mean_expression,mean_expression_cut)) %>% drop_na()
    cox_reg <- coxph(Surv(time,status)~mean_expression,data=datasurvfit)
    if (printHR==TRUE){
      print(ggforest(cox_reg,data=datasurvfit))
    }
    p_survival[i,"coxregwaldtestp"] <- summary(cox_reg)$waldtest["pvalue"]
    p_survival[i,"coxreg_coef"] <- summary(cox_reg)$conf.int[1]
    p_survival[i,"coxreg_LL"] <- summary(cox_reg)$conf.int[3]
    p_survival[i,"coxreg_UL"] <- summary(cox_reg)$conf.int[4]
    
    
    #test
    c_meanexp <- survival::survfit(Surv(time,status) ~ factor(mean_expression_cut),data=datasurvfit)
    p_survival[i,"logranktestp"] <- survminer::surv_pvalue(c_meanexp,data=datasurvfit)$pval
    p_survival[i,"medianOS1"]<-summary(c_meanexp)$table[1,"median"]
    p_survival[i,"medianOS2"]<-summary(c_meanexp)$table[2,"median"]
    p_survival[i,"Lower95CI1"]<-summary(c_meanexp)$table[1,"0.95LCL"]
    p_survival[i,"Lower95CI2"]<-summary(c_meanexp)$table[2,"0.95LCL"]
    p_survival[i,"Upper95CI1"]<-summary(c_meanexp)$table[1,"0.95UCL"]
    p_survival[i,"Upper95CI2"]<-summary(c_meanexp)$table[2,"0.95UCL"]
    
    if (printsurv==TRUE){
      print(ggsurvplot(c_meanexp,data=datasurvfit,
                       pval = TRUE, conf.int = TRUE,
                       risk.table = TRUE, # Add risk table
                       risk.table.col = "strata", # Change risk table color by groups
                       linetype = "strata", # Change line type by groups
                       surv.median.line = NULL, # Specify median survival
                       font.tickslab = c(12, "plain", "black"),
                       ggtheme = theme_classic(), # Change ggplot2 theme
                       xlab = "Time in years",
                       title=moduleNr,
                       legend = "bottom",
                       legend.labs = c("Low Expression", "High Expression"),
                       risk.table.y.text = FALSE))
    }
  }
  p_survival[,"logranktestpadj"]<-p.adjust(p_survival[,"logranktestp"],method="BH")
  p_survival[,"coxregwaldtestpadj"]<-p.adjust(p_survival[,"coxregwaldtestp"],method="BH")
  return(p_survival)
}

#' @title AMARETTO_survivalshiny
#'
#' @param AMARETTOresults
#' @param SurvivalAnnotation
#' @param idpatients
#' @param idvitalstatus
#' @param iddaystodeath
#' @param idfollowupdays
#' @param daystoyears
#' @param rightcensored
#'
#' @return A set of HTMLs, giving caracteristics of the communities
#' 
#' @import survminer
#' @import tidyverse
#' @import survival
#' @import shinythemes
#' @import gtools
#' @import shiny
#' @return
#' @export
AMARETTO_survivalshiny<- function(AMARETTOresults, SurvivalAnnotation, idpatients, idvitalstatus, iddaystodeath, idfollowupdays){
  
  mean_expression_modules <- t(AMARETTOresults$ModuleData)
  mean_expression_modules <- rownames_to_column(as.data.frame(mean_expression_modules),"idpatients")
  if (sum(SurvivalAnnotation[,idpatients] %in% mean_expression_modules[,"idpatients"],na.rm=TRUE)==0){
    stop("No overlap between patients ids")
  } else {
    print(paste0("Survival data will be calculated on ",length(SurvivalAnnotation[,idpatients] %in% mean_expression_modules[,"idpatients"]), " patients."))
  }
  
  SurvivalAnnotation <- suppressMessages(dplyr::inner_join(SurvivalAnnotation, mean_expression_modules %>% dplyr::rename(!!idpatients :="idpatients")))
  
  daysdeath <- as.numeric(SurvivalAnnotation[grep("dead", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE),iddaystodeath])
  dead <- SurvivalAnnotation[grep("dead", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE),idvitalstatus]
  daysfollowup <- as.numeric(SurvivalAnnotation[grep("alive", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE),idfollowupdays])
  alive <- SurvivalAnnotation[grep("alive", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE),idvitalstatus]
  time<- c(daysfollowup,daysdeath)
  status <- tolower(c(alive, dead))
  
  modules <- paste0("Module ", seq(1:AMARETTOresults$NrModules))
  
  ui<-fluidPage(
    theme = shinytheme("paper"),
    titlePanel("Survival Analysis of the Modules"),
    sidebarPanel(
      selectInput(inputId = "Module", 
                  label="Which module do you want to analyse:", 
                  choices = modules, 
                  selected = modules[1]),
      
      checkboxInput(inputId = "daystoyears",
                    label = "Convert days to years",
                    value = FALSE),
      
      checkboxInput(inputId = "RCensor",
                    label = "Right censor the data",
                    value = FALSE),
      
      numericInput(inputId = "RCensorNumb",
                  label= "The data must be right censored at:",
                  value = 5),
      
      submitButton(text="Calculate")
      
    ),
    mainPanel(
      h3("Survival Plot"),
      plotOutput(outputId = "SurvivalPlot",height = "600px"),
      br(),
      h3("Hazard Ratio"),
      plotOutput(outputId = "HRPlot")
    )
  )
  
  server<-function(input,output){
    
    datasurvfit <- reactive({
      moduleNr <- input$Module
      moduleNr <-sub(" ","_",moduleNr)
      mean_dead <- SurvivalAnnotation[grep("dead", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE), moduleNr]
      mean_alive <- SurvivalAnnotation[grep("alive", SurvivalAnnotation[,idvitalstatus],ignore.case=TRUE), moduleNr]
      mean_expression <- c(mean_alive, mean_dead)
      mean_expression_cut <- quantcut(mean_expression, q=2, na.rm=TRUE)
      if(input$daystoyears == TRUE){
        time <- time/365
      }
      if (input$RCensor == TRUE){
        rightcensored<-input$RCensorNumb
        status[time>=rightcensored] <- "alive"
        time[time>=rightcensored] <- rightcensored
      }
      status <- as.numeric(as.factor(status))
      return(as.data.frame(cbind(time, status, mean_expression, mean_expression_cut)) %>% drop_na())
    })
    
    
    output$SurvivalPlot<-renderPlot({
      c_meanexp <- survival::survfit(Surv(time, status) ~ factor(mean_expression_cut), data=datasurvfit())
      ggsurvplot(c_meanexp, data=datasurvfit(),
                       pval = TRUE, conf.int = TRUE,
                       risk.table = TRUE, # Add risk table
                       risk.table.col = "strata", # Change risk table color by groups
                       linetype = "strata", # Change line type by groups
                       surv.median.line = NULL, # Specify median survival
                       font.tickslab = c(12, "plain", "black"),
                       ggtheme = theme_classic(), # Change ggplot2 theme
                       xlab = "Time in years",
                       legend = "bottom",
                       legend.labs = c("Low Expression", "High Expression"),
                       risk.table.y.text = FALSE)

    })
    
    output$HRPlot <- renderPlot({
      cox_reg <- coxph(Surv(time, status)~mean_expression_cut, data=datasurvfit())
      ggforest(cox_reg,data=datasurvfit(),fontsize = 1)
    })
  }
  
  app<-shinyApp(ui=ui,server=server)
  runApp(app, port = 8888)
}
