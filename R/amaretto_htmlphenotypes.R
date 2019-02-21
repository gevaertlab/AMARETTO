#' @title AMARETTO_HTMLCharacterisation
#'
#' @param AMARETTOresults
#' @param annotation
#' @param report_address
#' @param survivalanalysis
#' @param phenotypeanalysis
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
#'
AMARETTO_HTMLCharacterisation<-function(AMARETTOresults, annotation, report_address="./",survivalanalysis=FALSE, phenotypeanalysis=FALSE, idpatients=NULL, idvitalstatus=NULL, iddaystodeath=NULL, idfollowupdays=NULL, daystoyears=FALSE, rightcensored=FALSE, parameters=NULL, typelist=NULL){
  
  full_path<-normalizePath(report_address)
  
  mean_expression_modules <- t(AMARETTOresults$ModuleData)
  mean_expression_modules <- rownames_to_column(as.data.frame(mean_expression_modules),"idpatients")
  
  if (sum(annotation[,idpatients] %in% mean_expression_modules[,"idpatients"],na.rm=TRUE)==0){
    stop("No overlap between patients ids")
  } else {
    print(paste0("Survival data will be calculated on ",length(annotation[,idpatients] %in% mean_expression_modules[,"idpatients"]), " patients."))
  }
  
  annotation <- suppressMessages(inner_join(annotation, mean_expression_modules %>% dplyr::rename(!!idpatients :="idpatients")))
  
  datatable<-AMARETTO_PhenAssociation(AMARETTOresults, annotation=annotation, idpatients=idpatients, parameters=parameters, typelist=typelist, printplots=FALSE)
  if (survivalanalysis == TRUE){

    daysdeath <- as.numeric(annotation[grep("dead", annotation[,idvitalstatus],ignore.case=TRUE),iddaystodeath])
    dead <- annotation[grep("dead", annotation[,idvitalstatus],ignore.case=TRUE),idvitalstatus]
    daysfollowup <- as.numeric(annotation[grep("alive", annotation[,idvitalstatus],ignore.case=TRUE),idfollowupdays])
    alive <- annotation[grep("alive", annotation[,idvitalstatus],ignore.case=TRUE),idvitalstatus]
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

    colnames(p_survival) <- c("coxregwaldtestp","coxregwaldtestpadj","coxreg_coef","coxreg_LL","coxreg_UL","logranktestp","logranktestpadj","medianOS1","medianOS2","Lower95CI1","Upper95CI1","Lower95CI2","Upper95CI2")
    rownames(p_survival) <- paste0("Module ",1:AMARETTOresults$NrModules)
    
    cluster <- makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
    registerDoParallel(cluster,cores=NrCores)
    #AMARETTOresults$NrModules
    p_survival<-foreach (i = 1:AMARETTOresults$NrModules, .packages = c('tidyverse','rmarkdown')) %dopar% {
      moduleNr <- paste0("Module_",i)
      mean_dead <- annotation[grep("dead", annotation[,idvitalstatus],ignore.case=TRUE), moduleNr]
      mean_alive <- annotation[grep("alive", annotation[,idvitalstatus],ignore.case=TRUE), moduleNr]
      mean_expression <- c(mean_alive,mean_dead)
      mean_expression_cut <- quantcut(mean_expression, q=2, na.rm=TRUE)
      datasurvfit <- as.data.frame(cbind(time,status,mean_expression,mean_expression_cut)) %>% drop_na()
      cox_reg <- coxph(Surv(time,status)~mean_expression,data=datasurvfit)
      coxstats<- c(summary(cox_reg)$waldtest["pvalue"], summary(cox_reg)$conf.int[1], summary(cox_reg)$conf.int[3],summary(cox_reg)$conf.int[4])
      c_meanexp <- survival::survfit(Surv(time,status) ~ factor(mean_expression_cut),data=datasurvfit)
      kmstats <- c(survminer::surv_pvalue(c_meanexp,data=datasurvfit)$pval, summary(c_meanexp)$table[1,"median"], summary(c_meanexp)$table[2,"median"], summary(c_meanexp)$table[1,"0.95LCL"],
      summary(c_meanexp)$table[2,"0.95LCL"], summary(c_meanexp)$table[1,"0.95UCL"], summary(c_meanexp)$table[2,"0.95UCL"])
      
      modulemd<-paste0(full_path,"/AMARETTOhtmls/survival/modules/module",i,".rmd")
      file.copy(system.file("templates/TemplateSurvivalModule.Rmd",package="AMARETTO"),modulemd)
      rmarkdown::render(modulemd,output_file = paste0("module",i,".html"), params = list(
        cox_reg=cox_reg,
        c_meanexp=c_meanexp,
        datasurvfit=datasurvfit,
        i=i),quiet = TRUE)
      file.remove(modulemd)
      
      return(c(coxstats,kmstats))
    }
    stopCluster(cluster)
    cat("The survival module htmls are finished.\n")
    p_survival<-data.frame(matrix(unlist(p_survival),byrow=T,ncol=11),stringsAsFactors=FALSE)
    colnames(p_survival) <- c("coxregwaldtestp","coxreg_coef","coxreg_LL","coxreg_UL","logranktestp","medianOS1","medianOS2","Lower95CI1","Upper95CI1","Lower95CI2","Upper95CI2")
    rownames(p_survival) <- paste0("Module ",1:AMARETTOresults$NrModules)
    rownames(p_survival) <- paste0("Module ",1:10)
    p_survival[,"logranktestpadj"]<-p.adjust(p_survival[,"logranktestp"],method="BH")
    p_survival[,"coxregwaldtestpadj"]<-p.adjust(p_survival[,"coxregwaldtestp"],method="BH")
    
    rmarkdown::render(system.file("templates/TemplateSurvivalIndex.Rmd",package="AMARETTO"),output_file=paste0(full_path,"/AMARETTOhtmls/survival/survivalindex.html"), params = list(
      p_survival=p_survival))
  }
  
  if (phenotypeanalysis == TRUE){
    cluster <- makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
    registerDoParallel(cluster,cores=NrCores)
    foreach (i = 1:AMARETTOresults$NrModules, .packages = c('tidyverse','rmarkdown','kableExtra')) %dopar% {
      
      moduleNr <- paste0("Module_",i)
      
      modulemd<-paste0(full_path,"/AMARETTOhtmls/phenotypes/modules/module",i,".rmd")
      #file.copy(system.file("templates/TemplatePhenotypesModule.Rmd",package="AMARETTO"),modulemd)
      file.copy("TemplatePhenotypesModule.Rmd",modulemd)
      parsed = rmarkdown::render(modulemd,output_file = paste0("module",i,".html"), params = list(
        parameters=parameters,
        typelist=typelist,
        annotation=annotation,
        i=i),quiet = TRUE)
      file.remove(modulemd)
    }
    phenotypetable<-AMARETTO_PhenAssociation(AMARETTOresults, annotation = annotation, idpatients = "ID", parameters=parameters, typelist = typelist,printplots = FALSE)
    
    rmarkdown::render(system.file("templates/TemplatePhenotypesIndex.Rmd",package="AMARETTO"),output_file=paste0(full_path,"/AMARETTOhtmls/phenotypes/phenotypesindex.html"), params = list(
      p_survival=p_survival))
    
    stopCluster(cluster)
  }
}
