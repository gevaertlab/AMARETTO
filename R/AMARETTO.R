AMARETTO_Initialize <- function(MA_matrix=MA_matrix,CNV_matrix=NULL,MET_matrix=NULL, Driver_list = NULL,NrModules,
                                VarPercentage,PvalueThreshold=0.001,RsquareThreshold=0.1,pmax=10,NrCores=1,OneRunStop=0){
  
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
  RegulatorInfo=CreateRegulatorData(MA_matrix=MA_matrix,CNV_matrix=CNV_matrix,MET_matrix=MET_matrix,Driver_list=Driver_list,PvalueThreshold=PvalueThreshold,RsquareThreshold=RsquareThreshold)
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

AMARETTO_VisualizeModule <- function(AMARETTOinit,AMARETTOresults,CNV_matrix,MET_matrix,ModuleNr) {
    # getting the data
    if (ModuleNr>AMARETTOresults$NrModules){
        cat('\tCannot plot Module',ModuleNr,'since the total number of modules is',AMARETTOresults$N,'.\n')
    }
    else {
        ModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
        currentRegulators = AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
        RegulatorData=AMARETTOinit$RegulatorData[currentRegulators,]
        ModuleGenes=rownames(ModuleData)
        cat('Module',ModuleNr,'has',length(rownames(ModuleData)),'genes and',length(currentRegulators),'regulators for',length(colnames(ModuleData)),'samples.\n')
  
        # Clustering the module itself
        SampleClustering=hclust(dist(t(ModuleData)), method = "complete", members = NULL)    
        GeneClustering=hclust(dist(ModuleData), method = "complete", members = NULL)
        ClustRegulatorData <- RegulatorData[,SampleClustering$order]
        ClustModuleData <- ModuleData[GeneClustering$order,SampleClustering$order]
        ClustCombinedData <- rbind(ClustModuleData,ClustRegulatorData)
  
        # create annotations 
        Alterations <- rep(0,nrow(ClustCombinedData))
        Alterations[(nrow(ModuleData)+1):nrow(ClustCombinedData)] <- rowSums(cbind(10*AMARETTOinit$RegulatorAlterations$Summary[currentRegulators,1],AMARETTOinit$RegulatorAlterations$Summary[currentRegulators,2]))
        if (length(which(Alterations==11))>0){
            if (length(which(Alterations==1))>0){
                if (length(which(Alterations==10))>0){
                    case = 1 # everything
                } else {
                    case = 2 # both and MET
                }
            } else {
                if (length(which(Alterations==10))>0){
                    case = 3 #both and CNV
                } else {
                    case =4 #Both only
                }
            }
        } else {
            if (length(which(Alterations==1))>0){
                if (length(which(Alterations==10))>0){
                    case = 5 # CNV and MET
                } else {
                    case = 6 # MET only
                }
            } else {
                if (length(which(Alterations==10))>0){
                    case = 7 # CNV only
                }
            }
        }
        Alterations[which(Alterations==1)] <- rep("Methylation aberrations",length(which(Alterations==1)))
        Alterations[which(Alterations=="10")] <- rep("Copy number alterations",length(which(Alterations=="10")))
        Alterations[which(Alterations=="11")] <- rep("Both methylation and copy number alterations",length(which(Alterations=="11")))
        Alterations[which(Alterations=="0")] <- rep(" ",length(which(Alterations=="0")))
        Alterations <- data.frame(Alterations)
  
        if (case==1){
            ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Copy number alterations"="gray","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white","Methylation aberrations"="black","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue"," "= "white","Both methylation and copy number alterations"="bisque4")),
                 which = "row", width = unit(1, "cm"),name="")
        } 
        if (case==2){
            ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Both methylation and copy number alterations"="bisque4"," "="white","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white","Methylation aberrations"="black","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue")),
                 which = "row", width = unit(1, "cm"),name="")
        }
        if (case==3){
            ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Copy number alterations"="gray","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white","Both methylation and copy number alterations"="bisque4"," "="white","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue")),
                 which = "row", width = unit(1, "cm"),name="")
        }
        if (case==4){
            ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Both methylation and copy number alterations"="bisque4"," "="white","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue")),
                 which = "row", width = unit(1, "cm"),name="")
        }
        if (case==5){
            ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Copy number alterations"="gray","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white","Methylation aberrations"="black","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue")),
                 which = "row", width = unit(1, "cm"),name="")
        }
        if (case==6){
            ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Methylation aberrations"="black","Hyper-methylated gene" = "yellow","Hypo-methylated gene" = "cornflowerblue"," "="white")),
                 which = "row", width = unit(1, "cm"),name="")
        }
        if (case==7){
            ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Copy number alterations"="gray","Amplified gene"="sienna1","Deleted gene "="cyan"," "="white")),
                 which = "row", width = unit(1, "cm"),name="")
        }
        CNVreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$CNV),currentRegulators)
        METreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$MET),currentRegulators)
        if (length(CNVreg)>0){
            CNVData <- matrix(0,nrow=length(CNVreg),ncol=ncol(ModuleData))
            colnames(CNVData) <- colnames(ModuleData)[SampleClustering$order]
            rownames(CNVData) <- CNVreg
            CNVData[,colnames(CNV_matrix)] <- CNV_matrix[rownames(CNVData),]  
            CNVData[which(CNVData>0)] <- "Amplified"  # amplification
            CNVData[which(CNVData<0)] <- "Deleted"  # deletion
            CNVData[which(CNVData==0)] <- " " # nothing      

            if (length(METreg)>0){  
                METData <- matrix(0,nrow=length(METreg),ncol=ncol(ModuleData))
                colnames(METData) <- colnames(ModuleData)[SampleClustering$order]
                rownames(METData) <- METreg
                METData[,colnames(MET_matrix)] <- MET_matrix[rownames(METData),]  
                METData[which(METData>0)] <- "Hyper-methylated"  # hyper
                METData[which(METData<0)] <- "Hypo-methylated"  # hypo
                METData[which(METData==0)] <- " " # nothing      
                Genes <- data.frame(t(METData),No=rep(" ",ncol(ModuleData)),t(CNVData))
            } else {
                Genes <- data.frame(No=rep(" ",ncol(ModuleData)),t(CNVData))
            }
        } else {
            METData <- matrix(0,nrow=length(METreg),ncol=ncol(ModuleData))
            colnames(METData) <- colnames(ModuleData)[SampleClustering$order]
            rownames(METData) <- METreg
            METData[,colnames(MET_matrix)] <- MET_matrix[rownames(METData),]  
            METData[which(METData>0)] <- "Hyper-methylated"  # hyper
            METData[which(METData<0)] <- "Hypo-methylated"  # hypo
            METData[which(METData==0)] <- " " # nothing      
            Genes <- data.frame(t(METData),No=rep(" ",ncol(ModuleData)))
        }
        ColAnnotation <- c("Hyper-methylated" = "yellow", "Hypo-methylated" = "cornflowerblue"," "="white","Amplified"="sienna1","Deleted"="cyan")
        ColAnnotation <- rep(list(ColAnnotation),ncol(Genes))
        names(ColAnnotation) <- colnames(Genes)
        haRow = HeatmapAnnotation(df=Genes ,name="test",
             col = ColAnnotation,which="column",show_legend=FALSE)
  
        # plotting
        heatmap <- Heatmap(ClustCombinedData, name = "Gene expression", column_title = paste('Module',ModuleNr), cluster_rows=FALSE,cluster_columns=FALSE,show_column_dend=FALSE,show_column_names=FALSE,row_names_gp=gpar(col=c(rep("white",nrow(ModuleData)),rep("black",nrow(RegulatorData))),fontsize=10),
              column_title_gp = gpar(fontsize = 20, fontface = "bold"),split=c(rep("Module Genes",nrow(ModuleData)),rep(" Regulators",nrow(RegulatorData))),gap = unit(5, "mm"),
              #  col=colorRamp2(c(-max(abs(CombinedData)), 0, max(abs(CombinedData))), c("green", "black", "red")))
              col=colorRamp2(c(-max(abs(ClustCombinedData)), 0, max(abs(ClustCombinedData))), c("green", "black", "red")),heatmap_legend_param = list(color_bar = "continuous"),top_annotation=haRow)
        draw(heatmap+ha)
  
        #    nf <- layout(matrix(c(1,2),2,1, byrow=T), widths=c(6),heights=c(2,3), respect=T)
        #    layout.show(nf)
        #    image(1:length(colnames(ModuleData)),1:length(currentRegulators),t(RegulatorData[,SampleClustering$order]),col=rev(brewer.pal(11,"RdBu")),xlab='',ylab='',main=paste('Module',ModuleNr),axes=FALSE)
        #    axis(LEFT<-2, at=1:length(currentRegulators), labels=currentRegulators, las= 2,cex.axis=0.7)        
        #    image(1:length(colnames(ModuleData)),1:length(ModuleGenes),t(ModuleData[GeneClustering$order,SampleClustering$order]),col=rev(brewer.pal(11,"RdBu")),xlab='Samples',ylab='ModuleGenes',main='',axes=FALSE)
        #axis(LEFT<-2,at=1:length(ModuleGenes), labels=ModuleGenes, las= 2,cex.axis=0.7)    
        if (length(CNVreg)>0){
            if (length(METreg)>0){
                MeanMET <- floor(nrow(METData)/2)+1
                MeanCNV <- floor(nrow(CNVData)/2)+1+nrow(METData)+1
                for(an in colnames(Genes)) {
                    decorate_annotation(an, {
                        if (an=="No"){
                            grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="white"))
                        } else {
                            # annotation names on the right
                            grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="black"))
                        }
                        if (an==colnames(Genes)[MeanMET]){
                            grid.text("MET", unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = c("center","bottom"),rot=90)
                        }
                        if (an==colnames(Genes)[MeanCNV]){
                            grid.text("CNV", unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = c("center","bottom"),rot=90)
                        }
                        # annotation names on the left
                        # grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
                    })
                }
            } else {
                MeanCNV <- floor(nrow(CNVData)/2)+1+1
                for(an in colnames(Genes)) {
                    decorate_annotation(an, {
                    if (an=="No"){
                        grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="white"))
                    } else { 
                        # annotation names on the right
                        grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="black"))
                    }
                    if (an==colnames(Genes)[MeanCNV]){
                        grid.text("CNV", unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = c("center","bottom"),rot=90)
                    }
                    # annotation names on the left
                    # grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
                })
            }
        }
    } else {
        if (length(METreg)>0){
            MeanMET <- floor(nrow(METData)/2)+1+1
                for(an in colnames(Genes)) {
                    decorate_annotation(an, {         
                        if (an=="No"){
                            grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="white"))
                        } else {
                            # annotation names on the right
                            grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left",gp=gpar(fontsize=10,col="black"))
                        }
                        if (an==colnames(Genes)[MeanMET]){
                            grid.text("MET", unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = c("center","bottom"),rot=90)
                        }
                        # annotation names on the left
                        # grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
                    })
                }
            }
        }  
    }
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

CreateRegulatorData <- function(MA_matrix=MA_matrix,CNV_matrix=NULL,MET_matrix=NULL, Driver_list = NULL,PvalueThreshold=0.001,RsquareThreshold=0.1) {
  
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
    DriversList <- Driver_Genes[[Driver_list]]
    DriversList <- intersect(DriversList, rownames(MA_matrix))
    cat('\tFound',length(DriversList),'driver genes from the input list.\n')
    
  }
  
  # set of drivers
  Drivers <- unique(c(rownames(CNV_matrix),MET_drivers,DriversList))
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
    
    driverList_alterations <- matrix(0,ncol=1,nrow=length(DriversList))
    colnames(driverList_alterations) <- Driver_list
    rownames(driverList_alterations) <- DriversList
    driverList_alterations[,Driver_list] <-1
    
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