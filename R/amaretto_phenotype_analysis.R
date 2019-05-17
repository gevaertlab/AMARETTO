#' Title
#'
#' @param AMARETTOinit 
#' @param AMARETTOresults 
#' @param sample_annotations 
#' @param phenotype_tests 
#'
#' @return
#' @export
#'
#' @examples
AMARETTO_statistical_test<-function(AMARETTOinit,AMARETTOresults,sample_annotations,phenotype_tests){
  final_list<-list()
  index_current=1
  for (i in 1:length(phenotype_tests$Phenotypes)){
    if(phenotype_tests[i,2]=="COXPROPHAZARDTIMETOEVENT"){
      pheno1<-as.character(phenotype_tests[i,1])
      pheno2<-as.character(phenotype_tests[which(phenotype_tests[,2]=="COXPROPHAZARDRIGHTCENSORING"),1])
      #kk<-which(phenotype_tests[,2]=="COXPROPHAZARDRIGHTCENSORING")
      test_sample_annotations<-dplyr::select(sample_annotations,"Sample",pheno1,pheno2)
      result<-all_modules_cox_prop_hazard_test(AMARETTOinit,AMARETTOresults,test_sample_annotations)
      final_list[[index_current]]<-list(phenotype1=pheno1,phenotype2=pheno2,test_type="COXPROPHAZARD",test_result=result)
      index_current=index_current+1
    }
    else if(phenotype_tests[i,2]=="SPEARMANCORRTEST"){
      pheno<-as.character(phenotype_tests[i,1])
      test_sample_annotations<-dplyr::select(sample_annotations,"Sample",pheno)
      result<-all_modules_SPEARMANCOR_test(AMARETTOinit,AMARETTOresults,test_sample_annotations)
      final_list[[index_current]]<-list(phenotype=pheno,test_type="SPEARMANCORRTEST",test_result=result)
      index_current=index_current+1
    }
    else if(phenotype_tests[i,2]=="PEARSONCORRTEST"){
      pheno<-as.character(phenotype_tests[i,1])
      test_sample_annotations<-dplyr::select(sample_annotations,"Sample",pheno)
      result<-all_modules_PEARSONCOR_test(AMARETTOinit,AMARETTOresults,test_sample_annotations)
      final_list[[index_current]]<-list(phenotype=pheno,test_type="PEARSONCORRTEST",test_result=result)
      index_current=index_current+1
    }
    else if(phenotype_tests[i,2]=="WILCOXONRANKSUMTEST"){
      pheno<-as.character(phenotype_tests[i,1])
      test_sample_annotations<-dplyr::select(sample_annotations,"Sample",pheno)
      result<-all_module_WILCOXONRANKSUM_test(AMARETTOinit,AMARETTOresults,test_sample_annotations)
      final_list[[index_current]]<-list(phenotype=pheno,test_type="WILCOXONRANKSUMTEST",test_result=result)
      index_current=index_current+1
    }
    else if(phenotype_tests[i,2]=="WILCOXONRANKSUMTESTPAIRED"){
      pheno<-as.character(phenotype_tests[i,1])
      test_sample_annotations<-dplyr::select(sample_annotations,"Sample",pheno)
      result<-all_module_WILCOXONRANKSUMPAIRED_test(AMARETTOinit,AMARETTOresults,test_sample_annotations)
      final_list[[index_current]]<-list(phenotype=pheno,test_type="WILCOXONRANKSUMTESTPAIRED",test_result=result)
      index_current=index_current+1
    }
    else if(phenotype_tests[i,2]=="KRUSKALWALLISTEST"){
      pheno<-as.character(phenotype_tests[i,1])
      test_sample_annotations<-dplyr::select(sample_annotations,"Sample",pheno)
      result<-all_module_KRUSKALWALLIS_test(AMARETTOinit,AMARETTOresults,test_sample_annotations)
      final_list[[index_current]]<-list(phenotype=pheno,test_type="KRUSKALWALLISTEST",test_result=result)
      index_current=index_current+1
    }
    else if(phenotype_tests[i,2]=="TIMESERIESEDGETIME"){
      pheno1<-as.character(phenotype_tests[i,1])
      pheno2<-as.character(phenotype_tests[which(phenotype_tests[,2]=="TIMESERIESEDGECONDITION"),1])
      # kk<-which(phenotype_tests[,2]=="TIMESERIESEDGECONDITION")
      test_sample_annotations<-dplyr::select(sample_annotations,"Sample",pheno1,pheno2)
      result<-all_module_TIMESERIESEDGE_test(AMARETTOinit,AMARETTOresults,test_sample_annotations)
      final_list[[index_current]]<-list(phenotype1=pheno1,phenotype2=pheno2,test_type="TIMESERIESEDGETIME",test_result=result)
      index_current=index_current+1
    }
  }
  return(final_list)
}

#' Title
#'
#' @param results_list 
#' @import tibble
#' @import purrr
#' @return
#' @export
#'
#' @examples
AMARETTO_unite_results<-function(results_list){
  united_df<-NULL
  for (i in 1:length(results_list)){
    if (results_list[[i]]$test_type=="COXPROPHAZARD"){
      df<-results_list[[i]]$test_result
      Hazard_ratio<-unlist(purrr::map(strsplit(df$"HR (95% CI for HR)"," "),1))
      statistics_str<-gsub("\\(","[",unlist(purrr::map(strsplit(df$"HR (95% CI for HR)"," "),1)))
      statistics_str<-gsub("\\)","]",statistics_str)
      statistics_str<-gsub("-",",",statistics_str)
      df<-add_column(df,
                     Phenotypes=paste0(results_list[[i]]$phenotype1," (COXPROPHAZARDTIMETOEVENT), ",results_list[[i]]$phenotype2," (COXPROPHAZARDRIGHTCENSORING)"),
                     Statistical_Test="Survival Analysis: Cox proportional hazards regression (Wald test)")%>%
        mutate(ModuleNr=paste("Module",ModuleNr))%>%
        mutate(Descriptive_Statistics=paste0("Beta: ",df$beta,
                                                 ", Hazard Ratio: ",Hazard_ratio,
                                                 ", 95% CI: ",statistics_str,
                                                 ", Wald Statistic: ",df$wald.test))%>%
        select(ModuleNr,Phenotypes,Statistical_Test,p.value,q.value,Descriptive_Statistics)
    }
    else if (results_list[[i]]$test_type=="SPEARMANCORRTEST"){
      df<-results_list[[i]]$test_result
      df<-add_column(df,
                     Phenotypes=paste0(results_list[[i]]$phenotype," (SPEARMANCORRTEST)"),
                     Statistical_Test="Continuous or Ordinal Analysis: Spearman rank correlation analysis")%>%
        mutate(ModuleNr=paste("Module",ModuleNr))%>%
        add_column(Descriptive_Statistics=paste0("Correlation: ",unlist(df$estimate),
                                                  ", Statistic: ",unlist(df$statistic)))%>%
        select(ModuleNr,Phenotypes,Statistical_Test,p.value,q.value,Descriptive_Statistics)
      
      #%>%select(ModuleNr,p.value,q.value,estimate,statistic,parameter,null.value,alternative,method,data.name)
      
    }
    else if (results_list[[i]]$test_type=="PEARSONCORRTEST"){
      df<-results_list[[i]]$test_result
      df<-add_column(df,
                     Phenotypes=paste0(results_list[[i]]$phenotype," (PEARSONCORRTEST)"),
                     Statistical_Test="Continuous Analysis: Pearson linear correlation analysis")%>%
        mutate(ModuleNr=paste("Module",ModuleNr))%>%
        add_column(Descriptive_Statistics=paste0("Correlation: ",unlist(df$estimate),
                                                  ", 95% CI: ",df$conf.int,
                                                  ", Statistic: ",unlist(df$statistic)))%>%
        select(ModuleNr,Phenotypes,Statistical_Test,p.value,q.value,Descriptive_Statistics)
    }
    
    else if (results_list[[i]]$test_type=="KRUSKALWALLISTEST"){
      df<-results_list[[i]]$test_result
      df<-add_column(df,
                     Phenotypes=paste0(results_list[[i]]$phenotype," (KRUSKALWALLISTEST)"),
                     Statistical_Test="Nominal Multi-Class Analysis: Kruskal-Wallis test")%>%
        mutate(ModuleNr=paste("Module",ModuleNr))%>%
        add_column(Descriptive_Statistics=paste0(" Statistic: ",unlist(df$statistic)))%>%
        select(ModuleNr,Phenotypes,Statistical_Test,p.value,q.value,Descriptive_Statistics)
    }
    
    else if (results_list[[i]]$test_type=="WILCOXONRANKSUMTEST"){
      df<-results_list[[i]]$test_result
      statistics_str<-gsub("c","",as.character(df$conf.int))
      statistics_str<-gsub("\\(","[",statistics_str)
      statistics_str<-gsub("\\)","]",statistics_str)
      df<-add_column(df,
                     Phenotypes=paste0(results_list[[i]]$phenotype," (WILCOXONRANKSUMTEST)"),
                     Statistical_Test="Nominal Two-Class Analysis: Wilcoxon rank sum test")%>%
        mutate(ModuleNr=paste("Module",ModuleNr))%>%
        add_column(Descriptive_Statistics=paste0("Estimate: ",unlist(df$estimate),", 95% CI: ",statistics_str,", Statistics: ",unlist(df$statistic)))%>%
        select(ModuleNr,Phenotypes,Statistical_Test,p.value,q.value,Descriptive_Statistics)
    }
    
    else if (results_list[[i]]$test_type=="WILCOXONRANKSUMTESTPAIRED"){
      df<-results_list[[i]]$test_result
      df<-add_column(df,
                     Phenotypes=paste0(results_list[[i]]$phenotype," (WILCOXONRANKSUMTESTPAIRED)"),
                     Statistical_Test="Nominal Paired Two-Class Analysis: Wilcoxon rank sum test paired")%>%
        mutate(ModuleNr=paste("Module",ModuleNr))%>%
        add_column(Descriptive_Statistics=paste0("Statistic: ",unlist(df$statistic)))%>%
        select(ModuleNr,Phenotypes,Statistical_Test,p.value,q.value,Descriptive_Statistics)
    }
    
    else if (results_list[[i]]$test_type=="TIMESERIESEDGETIME"){
      df<-results_list[[i]]$test_result
      df<-add_column(df,
                     Phenotypes=paste0(results_list[[i]]$phenotype1," (TIMESERIESEDGECONDITION),",results_list[[i]]$phenotype2," (TIMESERIESEDGETIME)"),
                     Statistical_Test="Time Series Analysis: Edge time course analysis (Likelihood ratio test)")%>%
        mutate(ModuleNr=paste("Module",ModuleNr))%>%add_column(Descriptive_Statistics=" ")%>%arrange(q.value)
        select(ModuleNr,Phenotypes,Statistical_Test,p.value,q.value,Descriptive_Statistics)
    }
    united_df<-rbind(united_df,as.data.frame(df))
  }
  
  return(united_df)
}

#' Title create_feature_matrix
#'
#' @param AMARETTOinit 
#' @param AMARETTOresults 
#' @param sample_annotation_df 
#' @param module_number 
#' @param Sample 
#'
#' @return
#' @export
#'
#' @examples
create_feature_matrix<-function(AMARETTOinit,AMARETTOresults,sample_annotation_df,module_number,Sample='Sample'){
  gene_names<-names(AMARETTOresults$ModuleMembership[AMARETTOresults$ModuleMembership[,1]==module_number,])
  MA_matrix_module<-t(AMARETTOinit$MA_matrix_Var[gene_names,])
  nr_common_samples<-sum(sample_annotation_df$Sample%in%rownames(MA_matrix_module))
  print(paste0("There are common ",nr_common_samples, " samples between gene expression data and the phenotype data"))
  feature_matrix<-as.data.frame(MA_matrix_module)%>%rownames_to_column(var=Sample)%>%inner_join(sample_annotation_df,by=Sample)%>%mutate(mean_all = rowMeans(select(., gene_names), na.rm = TRUE))
  return(feature_matrix)
}

#' Title cox_prop_hazard_test
#'
#' @param xx 
#' @param time_to_event 
#' @param right_censoring 
#'
#' @return
#' @export
#'
#' @examples
cox_prop_hazard_test<-function(xx,time_to_event,right_censoring){
  df=data.frame(module_mean=xx,COXPROPHAZARDTIMETOEVENT=as.numeric(time_to_event),COXPROPHAZARDRIGHTCENSORING=as.numeric(right_censoring))
  coxph<-coxph(Surv(df$COXPROPHAZARDTIMETOEVENT,df$COXPROPHAZARDRIGHTCENSORING)~df$module_mean)
  coxphs<-summary(coxph)
  p.value<-signif(coxphs$waldtest["pvalue"], digits=5)
  wald.test<-signif(coxphs$waldtest["test"], digits=5)
  beta<-signif(coxphs$coefficients[1], digits=5)
  HR<-signif(coxphs$coefficients[2], digits=5)
  HR.confint.lower <- signif(coxphs$conf.int[,"lower .95"], 5)
  HR.confint.upper <- signif(coxphs$conf.int[,"upper .95"], 5)
  HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
  COXPROPHAZARD_res<-c(beta, HR, wald.test, p.value)
  names(COXPROPHAZARD_res)<-c("beta", "HR (95% CI for HR)", "wald.test", "p.value")
  COXPROPHAZARD_results <- t(as.data.frame(COXPROPHAZARD_res, check.names = FALSE,stringsAsFactors=FALSE))
  # COXPROPHAZARD_formulas <- sapply(c("xx"), function(x) as.formula(paste('Surv(COXPROPHAZARDTIMETOEVENT, COXPROPHAZARDRIGHTCENSORING)~', x)))
  # COXPROPHAZARD_models <- lapply(COXPROPHAZARD_formulas, function(x){coxph(x, data = df)})
  # COXPROPHAZARD_res <- lapply(COXPROPHAZARD_models,
  #                             function(x){ 
  #                               x <- summary(x)
  #                               p.value<-signif(x$wald["pvalue"], digits=5)
  #                               wald.test<-signif(x$wald["test"], digits=5)
  #                               beta<-signif(x$coef[1], digits=5);#coeficient beta
  #                               HR <-signif(x$coef[2], digits=5);#exp(beta)
  #                               HR.confint.lower <- signif(x$conf.int[,"lower .95"], 5)
  #                               HR.confint.upper <- signif(x$conf.int[,"upper .95"], 5)
  #                               HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
  #                               res<-c(beta, HR, wald.test, p.value)
  #                               names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", "p.value")
  #                               return(res)
  #                               #return(exp(cbind(coef(x),confint(x))))
  #                             })
  # COXPROPHAZARD_results <- t(as.data.frame(COXPROPHAZARD_res, check.names = FALSE,stringsAsFactors=FALSE))
  return(COXPROPHAZARD_results)
}

#' Title
#'
#' @param AMARETTOinit 
#' @param AMARETTOresults 
#' @param test_sample_annotations 
#'
#' @return
#' @export
#'
#' @examples
all_modules_cox_prop_hazard_test<-function(AMARETTOinit,AMARETTOresults,test_sample_annotations){
  df_result<-NULL
  for (module_number in 1:AMARETTOresults$NrModules){
    featureMatrix<-create_feature_matrix(AMARETTOinit,AMARETTOresults,sample_annotation_df=test_sample_annotations,module_number=module_number,Sample='Sample')
    time_to_event<-featureMatrix%>%dplyr::pull(colnames(test_sample_annotations)[2])  
    right_censoring<-featureMatrix%>%dplyr::pull(colnames(test_sample_annotations)[3])  
    result<-cox_prop_hazard_test(featureMatrix$mean_all,time_to_event,right_censoring)
    df_result<-rbind(df_result,result)
  }
  df_result<-as.data.frame(df_result,stringsAsFactors=FALSE)%>%rownames_to_column()
  df_result$p.value<-as.numeric(unlist(df_result$p.value))
  df_result$p.value<-signif(df_result$p.value, digits=5)
  q.value<-p.adjust(df_result$p.value,method = "BH")
  df_result<-add_column(df_result,q.value=q.value,ModuleNr=1:AMARETTOresults$NrModules)
  #result$ModuleNr=module_number
  return(df_result)
}


#' Title
#'
#' @param AMARETTOinit 
#' @param AMARETTOresults 
#' @param test_sample_annotations 
#'
#' @return
#' @export
#'
#' @examples
all_modules_SPEARMANCOR_test<-function(AMARETTOinit,AMARETTOresults,test_sample_annotations){
  df_result<-NULL
  for (module_number in 1:AMARETTOresults$NrModules){
    featureMatrix<-create_feature_matrix(AMARETTOinit,AMARETTOresults,sample_annotation_df=test_sample_annotations,module_number=module_number,Sample='Sample')
    spearmancorr_vector<-featureMatrix%>%dplyr::pull(colnames(test_sample_annotations)[2])  
    df=data.frame(module_average=featureMatrix$mean_all,SPEARMANCORRTEST=spearmancorr_vector)
    SPEARMANCORRTEST_results <- cor.test(df$module_average,df$SPEARMANCORRTEST, method = "spearman", alternative = "two.sided")
    df_result<-rbind(df_result,SPEARMANCORRTEST_results)
  }
  df_result<-as.data.frame(df_result,stringsAsFactors=FALSE)%>%rownames_to_column()
  df_result$p.value<-as.numeric(unlist(df_result$p.value))
  df_result$p.value<-signif(df_result$p.value, digits=5)
  q.value<- p.adjust(df_result$p.value,method = "BH")
  df_result<-add_column(df_result,q.value=q.value,ModuleNr=1:AMARETTOresults$NrModules)
  return(df_result)
}

#' Title
#'
#' @param AMARETTOinit 
#' @param AMARETTOresults 
#' @param test_sample_annotations 
#'
#' @return
#' @export
#'
#' @examples
all_modules_PEARSONCOR_test<-function(AMARETTOinit,AMARETTOresults,test_sample_annotations){
  df_result<-NULL
  for (module_number in 1:AMARETTOresults$NrModules){
    featureMatrix<-create_feature_matrix(AMARETTOinit,AMARETTOresults,sample_annotation_df=test_sample_annotations,module_number=module_number,Sample='Sample')
    pearsoncorr_vector<-featureMatrix%>%dplyr::pull(colnames(test_sample_annotations)[2])  
    df=data.frame(module_average=featureMatrix$mean_all,PEARSONCORRTEST=pearsoncorr_vector)
    PEARSONCORRTEST_results <- cor.test(df$module_average,df$PEARSONCORRTEST, method = "pearson",  conf.level = 0.95, alternative = "two.sided")
    df_result<-rbind(df_result,PEARSONCORRTEST_results)
  }
  df_result<-as.data.frame(df_result,stringsAsFactors=FALSE)%>%rownames_to_column()
  df_result$p.value<-as.numeric(unlist(df_result$p.value))
  df_result$p.value<-signif(df_result$p.value, digits=5)
  q.value<- p.adjust(df_result$p.value,method = "BH")
  df_result<-add_column(df_result,q.value=q.value,ModuleNr=1:AMARETTOresults$NrModules)
  return(df_result)
}


#' Title
#'
#' @param AMARETTOinit 
#' @param AMARETTOresults 
#' @param test_sample_annotations 
#'
#' @return
#' @export
#'
#' @examples
all_module_WILCOXONRANKSUM_test<-function(AMARETTOinit,AMARETTOresults,test_sample_annotations){
  df_result<-NULL
  for (module_number in 1:AMARETTOresults$NrModules){
    
    featureMatrix<-create_feature_matrix(AMARETTOinit,AMARETTOresults,sample_annotation_df=test_sample_annotations,module_number=module_number,Sample='Sample')
    wilcoxonsum_vector<-featureMatrix%>%dplyr::pull(colnames(test_sample_annotations)[2])  
    df<-data.frame(mean_all=featureMatrix$mean_all,wilcoxonsum_vector=as.factor(wilcoxonsum_vector),stringsAsFactors = FALSE)%>%drop_na()
    if(length(unique(df$wilcoxonsum_vector))!=2){
      stop(paste0("There is(are) ",length(unique(df$wilcoxonsum_vector)), " group(s) for ", colnames(test_sample_annotations)[2], ". Requiers two groups. The module Test is skiped."))
    }
    group1<-df$mean_all[df$wilcoxonsum_vector==unique(df$wilcoxonsum_vector)[1]]
    group2<-df$mean_all[df$wilcoxonsum_vector==unique(df$wilcoxonsum_vector)[2]]
    WILCOXONRANKSUMTEST_results <- wilcox.test(x=group1,y=group2, alternative = "two.sided", paired = FALSE, conf.int = TRUE, conf.level = 0.95)
    df_result<-rbind(df_result,WILCOXONRANKSUMTEST_results)
  }
  df_result<-as.data.frame(df_result,stringsAsFactors=FALSE)%>%rownames_to_column()
  df_result$p.value<-as.numeric(unlist(df_result$p.value))
  df_result$p.value<-signif(df_result$p.value, digits=5)
  q.value<- p.adjust(df_result$p.value,method = "BH")
  df_result<-add_column(df_result,q.value=q.value,ModuleNr=1:AMARETTOresults$NrModules)
  return(df_result)
}

#' Title
#'
#' @param AMARETTOinit 
#' @param AMARETTOresults 
#' @param test_sample_annotations 
#'
#' @return
#' @export
#'
#' @examples
all_module_WILCOXONRANKSUMPAIRED_test<-function(AMARETTOinit,AMARETTOresults,test_sample_annotations){
  df_result<-NULL
  for (module_number in 1:AMARETTOresults$NrModules){
    
    featureMatrix<-create_feature_matrix(AMARETTOinit,AMARETTOresults,sample_annotation_df=test_sample_annotations,module_number=module_number,Sample='Sample')
    wilcoxonsum_vector<-featureMatrix%>%dplyr::pull(colnames(test_sample_annotations)[2])  
    df<-data.frame(mean_all=featureMatrix$mean_all,wilcoxonsum_vector=as.factor(wilcoxonsum_vector),stringsAsFactors = FALSE)%>%drop_na()
    if(length(unique(df$wilcoxonsum_vector))!=2){
      stop(paste0("There is(are) ",length(unique(df$wilcoxonsum_vector)), " group(s) for ", colnames(test_sample_annotations)[2], ". Requiers two groups. The module Test is skiped."))
    }
    group1<-df$mean_all[df$wilcoxonsum_vector==df$wilcoxonsum_vector[[1]]]
    group2<-df$mean_all[df$wilcoxonsum_vector==df$wilcoxonsum_vector[[2]]]
    
    WILCOXONRANKSUMTESTPAIRED_results <- wilcox.test(x=group1,y=group2, alternative = "two.sided", paired = TRUE, exact = NULL, correct = TRUE, conf.int = FALSE, conf.level = 0.95)
    df_result<-rbind(df_result,WILCOXONRANKSUMTESTPAIRED_results)
  }
  df_result<-as.data.frame(df_result,stringsAsFactors=FALSE)%>%rownames_to_column()
  df_result$p.value<-as.numeric(unlist(df_result$p.value))
  df_result$p.value<-signif(df_result$p.value, digits=5)
  q.value<- p.adjust(df_result$p.value,method = "BH")
  df_result<-add_column(df_result,q.value=q.value,ModuleNr=1:AMARETTOresults$NrModules)
  return(df_result)
}

#' Title
#'
#' @param AMARETTOinit 
#' @param AMARETTOresults 
#' @param test_sample_annotations 
#'
#' @return
#' @export
#'
#' @examples
all_module_KRUSKALWALLIS_test<-function(AMARETTOinit,AMARETTOresults,test_sample_annotations){
  df_result<-NULL
  for (module_number in 1:AMARETTOresults$NrModules){
    featureMatrix<-create_feature_matrix(AMARETTOinit,AMARETTOresults,sample_annotation_df=test_sample_annotations,module_number=module_number,Sample='Sample')
    kruskalwallis_vector<-featureMatrix%>%dplyr::pull(colnames(test_sample_annotations)[2])  
    df=data.frame(module_average=featureMatrix$mean_all,KRUSKALWALLISTEST=as.factor(kruskalwallis_vector))
    KRUSKALWALLISTEST_results <- kruskal.test(x=df$module_average,g=df$KRUSKALWALLISTEST)
    df_result<-rbind(df_result,KRUSKALWALLISTEST_results)
  }
  df_result<-as.data.frame(df_result,stringsAsFactors=FALSE)%>%rownames_to_column()
  df_result$p.value<-as.numeric(unlist(df_result$p.value))
  df_result$p.value<-signif(df_result$p.value, digits=5)
  q.value<- p.adjust(df_result$p.value,method = "BH")
  df_result<-df_result%>%mutate(q.value=q.value)%>%mutate(ModuleNr=1:AMARETTOresults$NrModules)
  return(df_result)
}

#' Title all_module_TIMESERIESEDGE_test
#'
#' @param AMARETTOinit 
#' @param AMARETTOresults 
#' @param test_sample_annotations 
#' @import edge
#' @return
#' @export
#'
#' @examples
all_module_TIMESERIESEDGE_test<-function(AMARETTOinit,AMARETTOresults,test_sample_annotations){
  ModulesPhenotypes<-c()
  for (module_number in 1:AMARETTOresults$NrModules){
    featureMatrix<-create_feature_matrix(AMARETTOinit,AMARETTOresults,sample_annotation_df=test_sample_annotations,module_number=module_number,Sample='Sample')
    ModulesPhenotypes<-cbind(ModulesPhenotypes,featureMatrix$mean_all)
  }
  #View(ModulesPhenotypes)
  mods<-t(ModulesPhenotypes)
  cov.time<-featureMatrix%>%dplyr::pull(colnames(test_sample_annotations)[2])  
  cov.cond<-featureMatrix%>%dplyr::pull(colnames(test_sample_annotations)[3])
  length(cov.time)
  length(cov.cond)
  dim(mods)
  # cov.time<-TIMESERIESEDGETIME
  # cov.cond<-TIMESERIESEDGECONDITION
  # cov.replicate<-TIMESERIESEDGETIMEREPLICATE
  edge_obj <- edge::build_study(data = as.matrix(mods), adj.var = as.factor(cov.cond), tme = as.numeric(cov.time),sampling = "timecourse")
  #edge_obj <- build_study(data = as.matrix(mods), adj.var = cov.cond, tme = cov.time, ind = cov.replicate, sampling = "timecourse")
  #de_obj <- edge_obj 
  #full_matrix<-fullModel(de_obj)
  #null_matrix<-nullModel(de_obj)
  #ef_obj <- fit_models(de_obj, stat.type = "lrt")
  #alt_res <- resFull(ef_obj)
  #null_res <- resNull(ef_obj)
  #alt_fitted <- fitFull(ef_obj)
  #null_fitted <- fitNull(ef_obj)
  
  #de_lrt <- lrt(de_obj, nullDistn = "normal",pi0=1) 
  de_odp <- odp(edge_obj, bs.its = 50, verbose = FALSE,n.mods = 50)
  #summary(de_odp)
  sig_results <- qvalueObj(de_odp)
  # hist(sig_results)
  p.value<-as.numeric(unlist(sig_results$pvalues))
  p.value<-signif(p.value, digits=5)
  q.value<- p.adjust(p.value,method = "BH")
  TIMESERIESEDGE_results<-data.frame(ModuleNr=1:AMARETTOresults$NrModules,p.value=p.value,q.value=q.value,stringsAsFactors=FALSE)
  return(TIMESERIESEDGE_results)
}


