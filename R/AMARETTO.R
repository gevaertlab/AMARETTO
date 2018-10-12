AMARETTO_Initialize <- function(MA_matrix=MA_matrix,CNV_matrix=NULL,MET_matrix=NULL, Driver_list = NULL,NrModules,
                                VarPercentage,PvalueThreshold=0.001,RsquareThreshold=0.1,pmax=10,NrCores=1,OneRunStop=0, method= "union"){
  
  if(is.null(MET_matrix) & is.null(CNV_matrix) & is.null(Driver_list)) {
    stop("Please select the correct input data types")
  }
  if(is.null(CNV_matrix)  & is.null(Driver_list))  
  {
    if (ncol(MA_matrix)<2 || ncol(MET_matrix)==1){
      stop("AMARETTO cannot be run with less than two samples.\n")
    }
  }
  
  if(is.null(MET_matrix)  & is.null(Driver_list)) {
    if (ncol(MA_matrix)<2 || ncol(CNV_matrix)==1 ){
      stop("AMARETTO cannot be run with less than two samples.\n")
    }
  }
  
  if(!is.null(MA_matrix) & !is.null(CNV_matrix) & !is.null(MET_matrix)  & !is.null(Driver_list)) {
    if (ncol(MA_matrix)<2 || ncol(CNV_matrix)==1 || ncol(MET_matrix)==1){
      stop("AMARETTO cannot be run with less than two samples.\n")
    }
  }
  # Fixing default paramters
  AutoRegulation=2; Lambda2=0.0001; alpha=1-1e-06
  # Creating the parameters structure
  Parameters <- list(AutoRegulation=AutoRegulation,OneRunStop=OneRunStop,Lambda2=Lambda2,Mode='larsen',pmax=pmax,alpha=alpha)
  # Create cancer drivers
  RegulatorInfo=CreateRegulatorData(MA_matrix=MA_matrix,CNV_matrix=CNV_matrix,MET_matrix=MET_matrix,Driver_list=Driver_list,PvalueThreshold=PvalueThreshold,RsquareThreshold=RsquareThreshold,method=method)
  if (length(RegulatorInfo)>1){
    RegulatorData=RegulatorInfo$RegulatorData
    Alterations = RegulatorInfo$Alterations
    
    # Selecting gene expression data and Clustering
    MA_matrix_Var=geneFiltering('MAD',MA_matrix,VarPercentage)
    MA_matrix_Var=t(scale(t(MA_matrix_Var)))
    if (NrModules>=nrow(MA_matrix_Var)){
      stop(paste0("The number of modules is too large compared to the number of genes. Choose a number of modules smaller than ",nrow(MA_matrix_Var),".\n"))
    }
    KmeansResults=kmeans(MA_matrix_Var,NrModules,iter.max=100)
    Clusters=as.numeric(KmeansResults$cluster)
    names(Clusters) <- rownames(MA_matrix_Var)
    
    return(list(MA_matrix_Var=MA_matrix_Var,RegulatorData=RegulatorData,RegulatorAlterations=Alterations,ModuleMembership=Clusters,Parameters=Parameters,NrCores=NrCores))
  }
}


AMARETTO_Run <- function(AMARETTOinit) {
    if (length(AMARETTOinit)==0){
         cat('For cancer ',CancerSite,' no drivers were find during the initialization step of AMARETTO')
    } else{
        if (nrow(AMARETTOinit$RegulatorData)==1){
             cat('For cancer ',CancerSite,' only one driver is detected. AMARETTO cannot be run with less than two drivers.\n')
        } else {
             cat('Running AMARETTO on',length(rownames(AMARETTOinit$MA_matrix_Var)),'genes and',length(colnames(AMARETTOinit$MA_matrix_Var)),'samples.\n')
             cat('\tStopping if less then',0.01*length(rownames(AMARETTOinit$MA_matrix_Var)),'genes reassigned.\n')
             result=AMARETTO_LarsenBased(AMARETTOinit$MA_matrix_Var,AMARETTOinit$ModuleMembership,AMARETTOinit$RegulatorData,AMARETTOinit$Parameters,AMARETTOinit$NrCores)
    
             result$ModuleData=AMARETTO_CreateModuleData(AMARETTOinit,result)
             result$RegulatoryProgramData=AMARETTO_CreateRegulatorPrograms(AMARETTOinit,result)
    
             return(result)
        }
    }
}

AMARETTO_LarsenBased <- function(Data,Clusters,RegulatorData,Parameters,NrCores){
    # this will register nr of cores/threads, keep this here so the user can decide how many cores based on their hardware.
    registerDoParallel(cores=NrCores)
    ptm1 <- proc.time()

    RegulatorData_rownames=rownames(RegulatorData)
    Data_rownames=rownames(Data)

    AutoRegulation = Parameters$AutoRegulation
    RegulatorSign=array(0,length(RegulatorData_rownames))
    Lambda = Parameters$Lambda2
    OneRunStop = Parameters$OneRunStop
    if (AutoRegulation == 1){
        cat('\tAutoregulation is turned ON.\n')
    } else if (AutoRegulation == 2){
        cat('\tAutoregulation is turned ON.\n')
    } else {
        cat('\tAutoregulation is turned OFF.\n')
    }

    # main loop
    NrReassignGenes = length(Data_rownames)
    while (NrReassignGenes > 0.01*length(Data_rownames)){
  
        #STEP 1:  learning the regulatory program for each cluster
        ptm <- proc.time()
        switch(Parameters$Mode,
             larsen={
                regulatoryPrograms <- AMARETTO_LearnRegulatoryProgramsLarsen(Data,Clusters,RegulatorData,RegulatorSign,Lambda,AutoRegulation,alpha=Parameters$alpha,pmax=Parameters$pmax)
             }
        )
        ptm <- proc.time() - ptm
        printf("Elapsed time is %f seconds\n",ptm[3])
  
        NrClusters = length(unique(Clusters))
        sum = 0
        for(i in 1:NrClusters){
            sum = sum + Matrix::nnzero(regulatoryPrograms$Beta[i,] )
        }
        avg = sum / NrClusters
  
        printf("Average nr of regulators per module: %f \n",avg)
    
        PreviousClusters = Clusters # using the clusters where the regulatory program was trained and not the last clusters
        if (OneRunStop == 1){ break } 		# running only one iteration of optimization, useful in large comparisons
  
        #STEP 2: reassigning genes based on closed match to new regulatory programs
        ptm <- proc.time()
        ReassignGenesToClusters <- AMARETTO_ReassignGenesToClusters(Data,RegulatorData,regulatoryPrograms$Beta,Clusters,AutoRegulation)
        ptm <- proc.time() - ptm
        printf("Elapsed time is %f seconds\n",ptm[3])
  
        NrReassignGenes = ReassignGenesToClusters$NrReassignGenes
        Clusters = ReassignGenesToClusters$Clusters
        printf("Nr of reassignments is: %i \n",NrReassignGenes)
    }
    ptm1<- proc.time() - ptm1
    printf("Elapsed time is %f seconds\n",ptm1[3])
    
    # update results structure
    ModuleMembership=as.matrix(PreviousClusters)
    rownames(ModuleMembership)=rownames(Data)
    colnames(ModuleMembership)=c("ModuleNr")

    result <- list(NrModules = length(unique(Clusters)),RegulatoryPrograms = regulatoryPrograms$Beta,AllRegulators=rownames(RegulatorData),
        AllGenes = rownames(Data),ModuleMembership = ModuleMembership,AutoRegulationReport=regulatoryPrograms$AutoRegulationReport)
  
    return(result)
}

AMARETTO_LearnRegulatoryProgramsLarsen<-function(Data,Clusters,RegulatorData,RegulatorSign,Lambda,AutoRegulation,alpha,pmax){
  
    RegulatorData_rownames=rownames(RegulatorData)
    Data_rownames=rownames(Data)

    # stop has to be set because otherwise the algorithm continues until every
    # var is entered into the model
    #pmax = -10 # maximum nr of regulators that you want
    trace = 0
    NrFolds = 10
    NrClusters = length(unique(Clusters))
    NrGenes = nrow(Data)
    NrSamples = ncol(Data)
    NrInterpolateSteps = 100

    # autoregulation yes or no?
    if (AutoRegulation >= 1){
      #Beta = mat.or.vec((NrClusters),length(RegulatorData_rownames))
    } else if (AutoRegulation == 0) {
      BetaSpecial = list(NrClusters,1)
      RegulatorPositions = list(NrClusters,1)
    }
    #AutoRegulationReport = mat.or.vec(NrClusters,3)
    
    y_all = mat.or.vec(NrClusters,NrSamples)
    ClusterIDs = unique(Clusters)
    ClusterIDs = sort(ClusterIDs, decreasing = FALSE)
    cnt <- 1:NrClusters
    
    ptm1 <- proc.time()
    BetaY_all <- foreach(i=1:NrClusters,.combine=cbind,.init=list(list(),list(),list()),.packages = "glmnet") %dopar% {
        #for (i in 1:NrClusters){
      
        if (length(which(Clusters == ClusterIDs[i]))>1) {
            y = apply((Data[which(Clusters == ClusterIDs[i]),]),2,mean)
        } else {
            y = Data[which(Clusters == ClusterIDs[i]),]
        }
        CurrentClusterPositions = which(Clusters %in% ClusterIDs[i])
        nrGenesInClusters = length(CurrentClusterPositions)
      
        if (AutoRegulation >= 1){
            X = RegulatorData
        } else if (AutoRegulation == 0){
            X = RegulatorData[setdiff(RegulatorData_rownames,Data_rownames[CurrentClusterPositions]),]
        }
      
        fit = cv.glmnet(t(X), y,alpha = alpha, pmax = pmax)
      
        nonZeroLambdas <- fit$lambda[which(fit$nzero>0)]
        nonZeroCVMs <- fit$cvm[which(fit$nzero>0)]
      
        if(length(which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE)))==0){
        
            #for now: just print a warning, *although* this error WILL cause Amaretto to crash in a few steps.
            warnMessage <- paste0("\nOn cluster ",i," there were no cv.glm results that gave non-zero coefficients.")
            warning(warnMessage);
        
        }
      
        bestNonZeroLambda <- nonZeroLambdas[which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE))]
        b_o = coef(fit,s = bestNonZeroLambda)
        b_opt <- c(b_o[2:length(b_o)]) # removing the intercept.
      
        if (AutoRegulation == 2){ # autoregulation is allowed, if the regulator is stable also after removing it from the its regulated cluster.
        
            CurrentUsedRegulators = RegulatorData_rownames[which(b_opt!=0, arr.ind = T)]
            CurrentClusterMembers = Data_rownames[CurrentClusterPositions]
            nrIterations = 0 # 0 means no overlap initially
            while (length(CurrentClusterMembers[CurrentClusterMembers %in% CurrentUsedRegulators]) != 0){
                # keeping track of the removed Cluster Members:
                CurrentClusterMembers = setdiff(CurrentClusterMembers,CurrentUsedRegulators)  #problem here if the current cluster is empty
                nrCurrentClusterMembers = length(CurrentClusterMembers)
          
                if (nrCurrentClusterMembers > 0){
                    names = Data_rownames %in% CurrentClusterMembers
                    if (length(which(names==TRUE))>1){
                        y = apply((Data[names,]),2,mean) # only removing the used regulators from
                    } else {
                        y = Data[names,]
                    }
                    fit = cv.glmnet(t(X), y,alpha = alpha, pmax = pmax)
                    nonZeroLambdas <- fit$lambda[which(fit$nzero>0)]
                    nonZeroCVMs <- fit$cvm[which(fit$nzero>0)]
            
                    if(length(which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE)))==0){
              
                        #for now: just print a warning, *although* this error WILL cause Amaretto to crash in a few steps.
                        warnMessage <- paste0("\nOn cluster ",i," there were no cv.glm results that gave non-zero coefficients during the Autoregulation step.")
                        warning(warnMessage);
              
                    }
            
                    bestNonZeroLambda <- nonZeroLambdas[which(nonZeroCVMs==min(nonZeroCVMs,na.rm=TRUE))]
                    new_b_o = coef(fit,s = bestNonZeroLambda)
                    #what was used up until 09/09/2014 instead of bestNonZeroLambda:
                    #new_b_o = coef(fit,s = fit$lambda.1se)
                    new_b_opt <- c(new_b_o[2:length(b_o)])
            
                    CurrentUsedRegulators = RegulatorData_rownames[which(new_b_opt != 0)]
                    nrIterations = nrIterations + 1
                    b_opt = new_b_opt
                } else{
                    # no more cluster members left, the cluster is empty.
                    b_opt = rep(0,length(RegulatorData_rownames))
                }
            }
            Report <- c(length(CurrentClusterPositions),length(CurrentClusterMembers),nrIterations)
        
            #Report(1)=length(CurrentClusterPositions) #original clusters members
            #Report(2)=length(CurrentClusterMembers) #eventual nr of cluster members after removing members that were selected as regulator
            #Report(3)=nrIterations
            #AutoRegulationReport[i,]=t(Report)
            #AutoRegulationReport=Report
        }
  
        # need to do this after the autoregulation, otherwise autoregulation
        # can still add positive microRNA regulators
        if (sum(RegulatorSign[which(RegulatorSign != 0)]) > 0){ # there are limitations on the sign of the regulators
            RegulatorCheck = RegulatorSign * t(b_opt)
            WrongRegulators = which(RegulatorCheck < 0)
            if (length(WrongRegulators)  == 0){# just remove the wrong regulators
                b_opt[WrongRegulators] = 0
            }
        }
  
        if (AutoRegulation >= 1){
            #Beta[i,] = b_opt
        } else {
            BetaSpecial[i] = b_opt
            RegulatorPositions[i] = (RegulatorData_rownames %in% setdiff(RegulatorData_rownames,Data_rownames[CurrentClusterPositions])) # keeping track of the regulators' positions
        }
  
        #y_all[i,] = y
        list(b_opt,y,Report)
    }

    #ptm1<- proc.time() - ptm1
    #printf("Elapsed time is %f seconds\n",ptm1[3])
    if (AutoRegulation == 0){
        for (i in 1:NrClusters){
            Beta[i,RegulatorPositions[i]] = BetaSpecial[i]
        }
    }

    tmpPos=NrClusters+1

    Beta <- do.call(cbind, BetaY_all[1,2:tmpPos])
    Beta = t(Beta);
    colnames(Beta)=RegulatorData_rownames
    rownames(Beta)=gsub('result.','Module_',rownames(Beta))

    y_all<-do.call(cbind, BetaY_all[2,2:tmpPos])
    y_all = t(y_all);
    rownames(y_all)=gsub('result.','Module_',rownames(y_all))

    AutoRegulationReport<-do.call(cbind, BetaY_all[3,2:tmpPos])
    AutoRegulationReport = t(AutoRegulationReport)
    rownames(AutoRegulationReport)=gsub('result.','Module_',rownames(AutoRegulationReport))

    # calculating the error
    error = y_all - (Beta %*% RegulatorData)
    result <- list(Beta = Beta,error = error,AutoRegulationReport = AutoRegulationReport)
    return(result)
}

AMARETTO_ReassignGenesToClusters <- function(Data,RegulatorData,Beta,Clusters,AutoRegulation){
  
    RegulatorData_rownames=rownames(RegulatorData)
    Data_rownames=rownames(Data)

    NrGenes = nrow(Data)
    NrSamples = ncol(Data)
    NrReassignGenes = 0
    
    ##reassigning genes based on the Beta
    #getting the predictor data
    X = RegulatorData
    # creating the cluster "centroids"
    X1 = data.matrix(X)
    ModuleVectors = Beta %*% X1
    GeneNames = rownames(Data)
    
    #reassigning genes:
    #NewClusters = mat.or.vec(NrGenes,1)
    ptm1<- proc.time();
    nc <- foreach(i = 1:NrGenes, .combine = c) %dopar% {
        #for (i in 1:NrGenes){
        OldModule = Clusters[i]
        CurrentGeneVector = Data[i,,drop=FALSE]
        Correlations = cor(t(CurrentGeneVector),t(ModuleVectors))
        corr = data.matrix(Correlations,rownames.force = NA)
        #corr[is.na(corr)] <- -100000
        MaxCorrelation = max(corr,na.rm=TRUE)
        #MaxPosition = which(corr == max(corr,na.rm=TRUE))
        MaxPosition = which(signif(corr,digits=7) == signif(MaxCorrelation,digits=7))
        MaxPosition = MaxPosition[1] # this is new, to avoid two different reassignements
        
        if (AutoRegulation > 0){ #reassign unrestricted
            if (MaxPosition != OldModule){
                NrReassignGenes = NrReassignGenes + 1
            }
                #NewClusters[i] = MaxPosition
                NewClusters = MaxPosition
            }
        else {#only reassign if not a regulator THIS PART IS NOT VALIDATED !!! after update to RegulatorData instead of Regulators referred to as positions
            if (nnzero(rownames(RegulatorData_rownames) %in% GeneNames[i])!=0){
                 if (nnzero(which(which(GeneNames %in% rownames(RegulatorData_rownames)) %in% i) %in% which(Beta[MaxPosition,] != 0)) != 0){
                     if (MaxPosition != OldModule){
                         NrReassignGenes = NrReassignGenes + 1
                     }
                     #NewClusters[i] = MaxPosition
                     NewClusters = MaxPosition
                 } else {
                     #NewClusters[i] = OldModule
                     NewClusters = OldModule
                 }
            } else {
                if (MaxPosition != OldModule){
                    NrReassignGenes = NrReassignGenes + 1
                }
                #NewClusters[i] = MaxPosition
                NewClusters = MaxPosition
            }
        }
    }
    ptm1<- proc.time() - ptm1;
    NrReassignGenes = length(which(nc!=Clusters));
    result <- list(NrReassignGenes = NrReassignGenes,Clusters = nc)
    return(result)
}

AMARETTO_EvaluateTestSet <- function(AMARETTOresults,MA_Data_TestSet,RegulatorData_TestSet) {
    nrSamples = ncol(MA_Data_TestSet)
    RegulatorNames=rownames(RegulatorData_TestSet)
  
    #Iterating over the Modules
    stats = mat.or.vec(AMARETTOresults$NrModules,9)
    Rsquare = mat.or.vec(AMARETTOresults$NrModules,1)
    RsquareAjusted = mat.or.vec(AMARETTOresults$NrModules,1)
    modules <- list()
  
    for (i in 1:AMARETTOresults$NrModules){
        #check regulator presence
        currentRegulators = RegulatorNames[which(AMARETTOresults$RegulatoryPrograms[i,] != 0)]
        nrPresentRegulators = sum((rownames(RegulatorData_TestSet)  %in% currentRegulators))
        currentPresentRegulators = (currentRegulators %in% rownames(RegulatorData_TestSet))
        stats[i,1] = nrPresentRegulators
        stats[i,2] = length(currentRegulators)
  
        #checking the presence of the clusters
        currentClusterGenes = AMARETTOresults$AllGenes[which(AMARETTOresults$ModuleMembership[,1] == i)]
        nrPresentClusterGenes = sum((rownames(MA_Data_TestSet) %in% currentClusterGenes))
        stats[i,3] = nrPresentClusterGenes
        stats[i,4] = length(currentClusterGenes)
        stats[i,5] = stats[i,3] / stats[i,4] * 100
  
        #predict cluster expression in test set, always calculate but report
        #the totel percentage weight that is represented
        currentWeights = AMARETTOresults$RegulatoryPrograms[i,which(AMARETTOresults$RegulatoryPrograms[i,] != 0)]
        totalWeights = sum(abs(currentWeights))
        presentWeights = currentWeights[currentPresentRegulators]
        presentRegulators = currentRegulators[currentPresentRegulators]
        totalWeightsPercentage = sum(abs(presentWeights)) / totalWeights * 100
  
        modules[[i]] = currentClusterGenes[currentClusterGenes %in% rownames(MA_Data_TestSet)]
  
        # drop=FALSE, this solves the problem when you have only one regulator, so the previous version is not needed.
        if (nrPresentRegulators > 0) {
            predictions = (t(RegulatorData_TestSet[presentRegulators,,drop=FALSE])) %*% (presentWeights) # need to make sure that the first argument remains a matrix.
            predictions = data.matrix(predictions)
            if (length(modules[[i]]) !=0) {
                if (length(currentClusterGenes)>1){
                    outcome = apply(MA_Data_TestSet[currentClusterGenes,],2,mean)
                } else {
                     outcome = MA_Data_TestSet[currentClusterGenes,]
                }
                residuals = predictions-outcome
      
                # using explained variance as metric, since mean square error is not
                # enough, no baseline interpretation possible
                # SSreg=sumsqr(predictions-mean(outcome))
      
                SStot = sum((outcome-mean(outcome))^2)
                SSres = sum((predictions-outcome)^2)
        
                Rsquare[i] = 1 - (SSres / SStot)
                #Rsquare[i] = SSres
                RsquareAjusted[i] = Rsquare[i] - (nrPresentRegulators/(nrSamples - 1 - nrPresentRegulators))*(1-Rsquare[i])
                MSE = (1/nrSamples) * sum((predictions-outcome)^2)
      
                stats[i,6] = totalWeightsPercentage
                stats[i,7] = Rsquare[i]
                stats[i,8] = RsquareAjusted[i]
                stats[i,9] = MSE
            } else {
                stats[i,6] = totalWeightsPercentage
                stats[i,7] = 0
                stats[i,8] = 0
                stats[i,9] = 0
            }
        } else {
            stats[i,6] = 0
            stats[i,7] = 0
            stats[i,8] = 0
            stats[i,9] = 0
        }
    }
    #stats = list(stats,CellArrayOfNumToCellArrayofString(num2cell(1:AMARETTOresults.N)),{'nrPresReg' 'nrTotalReg' 'nrPresGen' 'nrTotGen' 'percPresGen' 'percWeightPresent' 'Rsquare' 'RsquareAdjusted'})

    dimnames(stats) <- list(rownames(stats, do.NULL = FALSE, prefix = "Module_"),
               c("nrPresReg" ,"nrTotalReg", "nrPresGen", "nrTotGen", "percPresGen",
               "percWeightPresent", "Rsquare", "RsquareAdjusted","MSE"))
    #res <- list(stats = stats,modules = modules)
  
    return(stats)
}

AMARETTO_VisualizeModule <- function(AMARETTOinit,AMARETTOresults,CNV_matrix=NULL,MET_matrix=NULL,ModuleNr,SAMPLE_annotation=NULL,ID=NULL,order_samples=NULL) {
  # getting the data
  if (ModuleNr>AMARETTOresults$NrModules){
    stop('\tCannot plot Module',ModuleNr,' since the total number of modules is',AMARETTOresults$N,'.\n')
  }
  ModuleData<-AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
  ModuleRegulators <- AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
  RegulatorData <- AMARETTOinit$RegulatorData[ModuleRegulators,]
  ModuleGenes <- rownames(ModuleData)
  cat('Module',ModuleNr,'has',length(rownames(ModuleData)),'genes and',length(ModuleRegulators),'regulators for',length(colnames(ModuleData)),'samples.\n')
  
  # create annotations
  Alterations<- rownames_to_column(as.data.frame(AMARETTOinit$RegulatorAlterations$Summary),"HGNC_symbol") %>% dplyr::rename(DriverList="Driver List") %>% dplyr::filter(HGNC_symbol %in% ModuleRegulators)
  Alterations<- Alterations %>% dplyr::mutate(CNVMet_Alterations=case_when(MET==1 & CNV==1~"Methylation and copy number alterations",
                                                                    CNV==1~"Copy number alterations",
                                                                    MET==1~"Methylation aberrations",
                                                                    MET==0 & CNV==0 ~"Not Altered"),
                                       DriversList_Alterations=case_when(DriverList==0~"Driver not predefined",
                                                                         DriverList==1~"Driver predefined"))
  
  ha_drivers <- HeatmapAnnotation(df = column_to_rownames(Alterations,"HGNC_symbol") %>% dplyr::select(CNVMet_Alterations,DriversList_Alterations), col = list(CNVMet_Alterations= c("Copy number alterations"="#eca400","Methylation aberrations"="#006992","Methylation and copy number alterations"="#d95d39","Not Altered"="lightgray"),DriversList_Alterations=c("Driver not predefined"="lightgray","Driver predefined"="#588B5B")),which = "column", height = unit(0.3, "cm"),name="",
                                  annotation_legend_param = list(title_gp = gpar(fontsize = 8),labels_gp = gpar(fontsize = 6)))
  
  if (is.null(MET_matrix) && is.null(CNV_matrix)){
    overlapping_samples<-colnames(ModuleData)
  } else if(is.null(MET_matrix)){
    overlapping_samples<-Reduce(intersect, list(colnames(CNV_matrix),colnames(ModuleData)))
  } else if(is.null(CNV_matrix)){
    overlapping_samples<-Reduce(intersect, list(colnames(MET_matrix),colnames(ModuleData)))
  } else {
    overlapping_samples<-Reduce(intersect, list(colnames(CNV_matrix),colnames(MET_matrix),colnames(ModuleData)))
  }
  
  #Clustering the samples based on target genes
  
  if(is.null(order_samples)){
    overlapping_samples_clust<-overlapping_samples[order(colMeans(ModuleData[,overlapping_samples]))]
  }else if(order_samples=="clust"){
    SampleClustering<-hclust(dist(t(ModuleData[,overlapping_samples])), method = "complete", members = NULL)
    overlapping_samples_clust<-overlapping_samples[SampleClustering$order]
  }else {
    print("ordering type not recognized, samples will be orderd based on mean expression of the module genes")
    overlapping_samples_clust<-overlapping_samples[order(colMeans(ModuleData[,overlapping_samples]))]
  }
  
  ClustRegulatorData <- t(RegulatorData[,overlapping_samples_clust])
  ClustModuleData <- t(ModuleData[,overlapping_samples_clust])
  Regwidth <- ncol(ClustRegulatorData)*0.5
  ha_Reg <- Heatmap(ClustRegulatorData, name = "Gene Expression", column_title = "Regulator Genes\nExpression",cluster_rows=FALSE,cluster_columns=TRUE,show_column_dend=FALSE,show_column_names=TRUE,show_row_names=FALSE,column_names_gp = gpar(fontsize = 6),top_annotation = ha_drivers,
                    column_title_gp = gpar(fontsize = 6, fontface = "bold"), col=colorRamp2(c(-max(abs(ClustRegulatorData)), 0, max(abs(ClustRegulatorData))), c("darkblue", "white", "darkred")),heatmap_legend_param = list(color_bar = "continuous",legend_direction = "horizontal",title_gp = gpar(fontsize = 8),labels_gp = gpar(fontsize = 6)), width = unit(Regwidth, "cm"))
  
  if(length(ClustModuleData)<50){
    fontsizeMo=6
  } else if (length(ClustModuleData)<200){
    fontsizeMo=4
  } else {fontsizeMo=2}
  
  ha_Module <- Heatmap(ClustModuleData, name = "", column_title = "Target Genes\nExpression",cluster_rows=FALSE,cluster_columns=TRUE,show_column_dend=FALSE,show_column_names=TRUE,show_row_names=FALSE,column_names_gp = gpar(fontsize = fontsizeMo),show_heatmap_legend = FALSE,
                       column_title_gp = gpar(fontsize = 6, fontface = "bold"), col=colorRamp2(c(-max(abs(ClustModuleData)), 0, max(abs(ClustModuleData))), c("darkblue", "white", "darkred")),heatmap_legend_param = list(color_bar = "continuous",legend_direction = "horizontal",title_gp = gpar(fontsize = 8),labels_gp = gpar(fontsize = 6)))
  
  ha_list<- ha_Reg + ha_Module
  
  if (!is.null(MET_matrix)){
    METreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$MET),ModuleRegulators)
    print("MET regulators will be included when available")
    if (length(METreg)>0){
      METData2 = METData = MET_matrix[unlist(Alterations %>% dplyr::filter(MET==1) %>% dplyr::select(HGNC_symbol)),overlapping_samples_clust]
      METData2[which(METData>0)] <- "Hyper-methylated"  # hyper
      METData2[which(METData<0)] <- "Hypo-methylated"  # hypo
      METData2[which(METData==0)] <- "Not altered" # nothing 
      METData2<-t(METData2)
      Metwidth=ncol(METData2)*0.5
      #number of colors
      Met_col=structure(c("#006992","#d95d39","white"),names=c("Hyper-methylated","Hypo-methylated","Not altered"))
      ha_Met <- Heatmap(METData2, name = "Methylation State", column_title = "Methylation State", cluster_rows=FALSE,cluster_columns=TRUE,show_column_dend=FALSE,show_column_names=TRUE,show_row_names=FALSE,column_names_gp = gpar(fontsize = 6),show_heatmap_legend = TRUE,
                        column_title_gp = gpar(fontsize = 6, fontface = "bold"), col = Met_col, width = unit(Metwidth, "cm"),heatmap_legend_param = list(title_gp = gpar(fontsize = 8),labels_gp = gpar(fontsize = 6)))
      ha_list<- ha_Met + ha_list
    }
  } 
  
  if (!is.null(CNV_matrix)){
    CNVreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$CNV),ModuleRegulators)
    print("CNV regulators will be included when available")
    if (length(CNVreg)>0){
      CNVData2 = CNVData = CNV_matrix[unlist(Alterations %>% dplyr::filter(CNV==1) %>% dplyr::select(HGNC_symbol)),overlapping_samples_clust] 
      CNVData2[which(CNVData>=0.1)] <- "Amplified"  # amplification
      CNVData2[which(CNVData<=(-0.1))] <- "Deleted"  # deletion
      CNVData2[which(CNVData<0.1 & CNVData>(-0.1))] <- "Not altered" # nothing
      CNVData2<-t(CNVData2)
      CNVwidth=ncol(CNVData2)*0.5
      #number of colors
      CNV_col=structure(c("#006992","#d95d39","white"),names=c("Deleted","Amplified","Not altered"))
      ha_CNV <- Heatmap(CNVData2, name = "CNV State", column_title = "CNV State", cluster_rows=FALSE,cluster_columns=TRUE,show_column_dend=FALSE,show_column_names=TRUE,show_row_names=FALSE,column_names_gp = gpar(fontsize = 6),show_heatmap_legend = TRUE,
                        column_title_gp = gpar(fontsize = 6, fontface = "bold"),col = CNV_col,width = unit(CNVwidth, "cm"),heatmap_legend_param = list(title_gp = gpar(fontsize = 8),labels_gp = gpar(fontsize = 6)))
      ha_list<-ha_CNV + ha_list
    }
  }
  
  if (!is.null(SAMPLE_annotation)){
    if (ID %in% colnames(SAMPLE_annotation)){
      SAMPLE_annotation_fil<-as.data.frame(SAMPLE_annotation) %>% dplyr::filter(!!as.name(ID) %in% overlapping_samples_clust)
      suppressWarnings(SAMPLE_annotation_fil<-left_join(as.data.frame(overlapping_samples_clust),SAMPLE_annotation_fil,by=c("overlapping_samples_clust"=ID)))
      SAMPLE_annotation_fil<-column_to_rownames(SAMPLE_annotation_fil,"overlapping_samples_clust")
      cat(nrow(SAMPLE_annotation_fil),"samples have an annotation.\n")
      cat(ncol(SAMPLE_annotation_fil),"annotations are added")
      ha_anot<-Heatmap(SAMPLE_annotation_fil, name="Sample Annotation", column_title = "Sample\nAnnotation", column_title_gp = gpar(fontsize = 6, fontface = "bold"), show_row_names=FALSE,width = unit(4, "mm"),column_names_gp = gpar(fontsize = 6),col=distinctColorPalette(nrow(unique(SAMPLE_annotation_fil))),heatmap_legend_param = list(title_gp = gpar(fontsize = 8),labels_gp = gpar(fontsize = 6)))
      ha_list<-ha_list + ha_anot
    } else {print("The ID is not identified as a column name in the annotation")}
  }
  
  ComplexHeatmap::draw(ha_list,heatmap_legend_side = "bottom",annotation_legend_side="bottom")
}

AMARETTO_CreateModuleData <- function(AMARETTOinit,AMARETTOresults) {
  
    # creating the module data based on the average expression of the genes in the module.
    ModuleData=matrix(0,AMARETTOresults$NrModules,length(colnames(AMARETTOinit$MA_matrix_Var)))
    rownames(ModuleData)=rownames(AMARETTOresults$AutoRegulationReport)
    colnames(ModuleData)=colnames(AMARETTOinit$MA_matrix_Var)
    for (ModuleNr in 1:AMARETTOresults$NrModules) {
        # getting the data
        currentModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership[,1]==ModuleNr,]
        if (length(which(AMARETTOresults$ModuleMembership[,1]==ModuleNr))>1){
            ModuleData[ModuleNr,]=colMeans(currentModuleData)
        } else {
            ModuleData[ModuleNr,]=currentModuleData
        }
    }

    return(ModuleData)
}

AMARETTO_CreateRegulatorPrograms <- function(AMARETTOinit,AMARETTOresults) {
  
    # creating the module data based on the average expression of the genes in the module.
 
    RegulatorProgramData=matrix(0,AMARETTOresults$NrModules,length(colnames(AMARETTOinit$MA_matrix_Var)))
    rownames(RegulatorProgramData)=rownames(AMARETTOresults$AutoRegulationReport)
    colnames(RegulatorProgramData)=colnames(AMARETTOinit$MA_matrix_Var)

    RegulatorNames=rownames(AMARETTOinit$RegulatorData)
    for (ModuleNr in 1:AMARETTOresults$NrModules) {
        currentRegulators = RegulatorNames[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
        weights=AMARETTOresults$RegulatoryPrograms[ModuleNr,currentRegulators]
        RegulatorData=AMARETTOinit$RegulatorData[currentRegulators,]
  
        RegulatorProgramData[ModuleNr,]=weights %*% RegulatorData
    }

    return(RegulatorProgramData)
}

CreateRegulatorData <- function(MA_matrix=MA_matrix,CNV_matrix=NULL,MET_matrix=NULL, Driver_list = NULL,PvalueThreshold=0.001,RsquareThreshold=0.1,method=NULL) {
  
  if(is.null(Driver_list)) DriversList <- NULL
  ### adding null CNV dataset
  if(is.null(CNV_matrix))  CNV_matrix <- matrix(0, nrow = 0, ncol = 0)
  
  # First removing genes with constant CNV status.
  if (nrow(CNV_matrix)>1){
    GeneVariances=rowVars(CNV_matrix)
    CNV_matrix=CNV_matrix[GeneVariances>=0.0001,]
  }
  if (nrow(CNV_matrix)>0){
    CNV_matrix=FindTranscriptionallyPredictive_CNV(MA_matrix,CNV_matrix,PvalueThreshold=PvalueThreshold,RsquareThreshold=RsquareThreshold)
  }
  
  cat('\tFound',length(rownames(CNV_matrix)),'CNV driver genes.\n')
  
  ### adding null MET dataset
  if(is.null(MET_matrix))  MET_matrix <- matrix(0, nrow = 0, ncol = 0)
  # based on MethylMix we have the following regulators
  if (nrow(MET_matrix)>1){
    GeneVariances=rowVars(MET_matrix)
    MET_matrix=MET_matrix[GeneVariances>=0.0001,]
  }
  
  MET_drivers <- intersect(rownames(MET_matrix),rownames(MA_matrix))
  #MET_matrix=FindTranscriptionallyPredictive_MET(MET_matrix,CNV_matrix,PvalueThreshold=PvalueThreshold,RsquareThreshold=RsquareThreshold)
  cat('\tFound',length(MET_drivers),'MethylMix driver genes.\n')
  
  # Driver list data
  if(!is.null(Driver_list))
  {
    #MET_matrix <- matrix(numeric(0),0,0) ;  CNV_matrix <- matrix(numeric(0),0,0) # ignore the rest input matrices
    DriversList <- Driver_list
    DriversList <- intersect(DriversList, rownames(MA_matrix))
    cat('\tFound',length(DriversList),'driver genes from the input list.\n')
    
  }
  
  # set of drivers
  #Drivers <- unique(c(rownames(CNV_matrix),MET_drivers,DriversList))
  if(is.null(Driver_list))  Drivers <- unique(c(rownames(CNV_matrix),MET_drivers,DriversList))
  dataset_drivers <- unique(c(rownames(CNV_matrix),MET_drivers))
  if(!is.null(Driver_list) & method=="union")  Drivers <- union(dataset_drivers,DriversList)
  if(!is.null(Driver_list) & method=="intersect")  Drivers <- intersect(dataset_drivers,DriversList)
  
  cat('\tFound a total of',length(Drivers),'unique drivers with your selected method.\n')
  
  # regulator data
  if (length(Drivers)==0){
    cat("AMARETTO doesn't find any driver genes.")
    return("No driver")
  } else {
    if (length(Drivers)==1){
      RegulatorData_temp <- matrix(0,1,ncol(MA_matrix))
      colnames(RegulatorData_temp) <- colnames(MA_matrix)
      rownames(RegulatorData_temp) <- Drivers
      RegulatorData_temp[1,] <- MA_matrix[Drivers,]
      RegulatorData <- RegulatorData_temp
    } else {
      RegulatorData=MA_matrix[Drivers,]
    }
    RegulatorData=t(scale(t(RegulatorData)))
    
    # alterations of all these drivers
    MET_aberrations <- matrix(0,ncol=3,nrow=length(MET_drivers))
    colnames(MET_aberrations) <- c("Hypo-methylated","No_change","Hyper-methylated")
    rownames(MET_aberrations) <- MET_drivers
    if  (length(MET_drivers)>0){
      MET_aberrations[,1] <- rowSums(MET_matrix[MET_drivers,]>0)/ncol(MET_matrix)
      MET_aberrations[,2] <- rowSums(MET_matrix[MET_drivers,]<0)/ncol(MET_matrix)
      MET_aberrations[,3] <- rowSums(MET_matrix[MET_drivers,]==0)/ncol(MET_matrix)
    }
    
    CNV_alterations <- matrix(0,ncol=3,nrow=nrow(CNV_matrix))
    colnames(CNV_alterations) <- c("Amplification","No_change","Deletion")
    rownames(CNV_alterations) <- rownames(CNV_matrix)
    if (nrow(CNV_alterations)>0){
      CNV_alterations[,1] <- rowSums(CNV_matrix>0)/ncol(CNV_matrix)
      CNV_alterations[,2] <- rowSums(CNV_matrix<0)/ncol(CNV_matrix)
      CNV_alterations[,3] <- rowSums(CNV_matrix==0)/ncol(CNV_matrix)
    }
    
    driverList_alterations <- matrix(1,ncol=1,nrow=length(DriversList))
    #colnames(driverList_alterations) <- Driver_list
    rownames(driverList_alterations) <- DriversList
    #driverList_alterations <-1
    
    Alterations <- matrix(0,nrow=length(Drivers),ncol=3)
    colnames(Alterations) <- c("CNV","MET","Driver List")
    rownames(Alterations) <- Drivers
    
    Alterations[which(Drivers%in% rownames(CNV_matrix)),1] <- rep(1,length(which(Drivers%in% rownames(CNV_matrix))))
    Alterations[which(Drivers%in% rownames(MET_matrix)),2] <- rep(1,length(which(Drivers%in% rownames(MET_matrix))))
    Alterations[which(Drivers%in% rownames(driverList_alterations)),3] <- rep(1,length(which(Drivers %in% rownames(driverList_alterations))))
    
    Alterations <- list(MET=MET_aberrations,CNV=CNV_alterations,Driver_list = driverList_alterations,Summary=Alterations)
    return(list(RegulatorData=RegulatorData,Alterations=Alterations))
  }
}

FindTranscriptionallyPredictive_CNV<- function(MA_matrix,CNV_matrix,PvalueThreshold=0.001,RsquareThreshold=0.1) {
  
    # Overlapping the data sets
    OverlapGenes=Reduce(intersect,list(rownames(MA_matrix),rownames(CNV_matrix)))
    OverlapSamples=Reduce(intersect,list(colnames(MA_matrix),colnames(CNV_matrix)))
    if (length(OverlapGenes)==1){
        CNV_TCGA_temp = MA_matrix_temp <- matrix(0,length(OverlapGenes),length(OverlapSamples))
        rownames(CNV_TCGA_temp) = rownames(MA_matrix_temp) <- OverlapGenes
        colnames(CNV_TCGA_temp) = colnames(MA_matrix_temp) <- OverlapSamples

        CNV_TCGA_temp[1,] <- CNV_matrix[OverlapGenes,OverlapSamples]
        MA_matrix_temp[1,] <- MA_matrix[OverlapGenes,OverlapSamples]
        CNV_matrix = CNV_TCGA_temp
        MA_matrix = MA_matrix_temp
    } else {
        CNV_matrix=CNV_matrix[OverlapGenes,OverlapSamples]
        MA_matrix=MA_matrix[OverlapGenes,OverlapSamples]
    }
    if (length(OverlapGenes)>0 && length(OverlapSamples)>0){
        # Linear modeling
        CNVdrivers=c()
        for (i in 1:length(rownames(CNV_matrix))) {
            res=lm(MA_matrix[i,]~CNV_matrix[i,])
            res.summary=summary(res)
            if (res$coefficients[2]>0 & res.summary$coefficients[2,4]<PvalueThreshold & res.summary$r.squared>RsquareThreshold) {
                CNVdrivers=c(CNVdrivers,rownames(CNV_matrix)[i])
            }
        }

        if (length(CNVdrivers)==1){
            CNV_matrix_temp <- matrix(0,1,ncol(CNV_matrix))
            colnames(CNV_matrix_temp) = colnames(CNV_matrix)
            rownames(CNV_matrix_temp) = CNVdrivers
            CNV_matrix_temp[1,] <- CNV_matrix[CNVdrivers,]
            CNV_matrix = CNV_matrix_temp
        } else {
            CNV_matrix=CNV_matrix[CNVdrivers,]
        }
    }

    # Returning the CNV selected data
    return(CNV_matrix=CNV_matrix)
}

geneFiltering<- function(Type,MAdata,Percentage) {
    switch(Type,
        Variance={
            GeneVariances=rowVars(MAdata)            
            tmpResult=sort(GeneVariances,decreasing=TRUE)
            SortedGenes=names(tmpResult)
            tmpNrGenes=round(length(rownames(MAdata))*Percentage/100)
            MAdata_Filtered=MAdata[SortedGenes[1:tmpNrGenes],]
        },
        MAD={
            GeneVariances=rowMads(MAdata)            
            names(GeneVariances)=rownames(MAdata)
            tmpResult=sort(GeneVariances,decreasing=TRUE)
            SortedGenes=names(tmpResult)
            tmpNrGenes=round(length(rownames(MAdata))*Percentage/100)
            MAdata_Filtered=MAdata[SortedGenes[1:tmpNrGenes],]
        }
    )      
    return(MAdata_Filtered)     
}

printf <- function(...) cat(sprintf(...))
