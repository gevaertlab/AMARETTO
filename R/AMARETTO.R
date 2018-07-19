AMARETTO_Initialize <- function(MA_matrix,CNV_matrix,MET_matrix,NrModules,
                                VarPercentage,PvalueThreshold=0.001,RsquareThreshold=0.1,pmax=10,NrCores=1,OneRunStop=0){
  
  # Fixing default paramters
  AutoRegulation=2
  Lambda2=0.0001
  alpha=1-1e-06
  
  # Creating the parameters structure
  Parameters <- list(AutoRegulation=AutoRegulation,OneRunStop=OneRunStop,Lambda2=Lambda2,Mode='larsen',pmax=pmax,alpha=alpha)
  
  # Create cancer drivers
  RegulatorInfo=CreateRegulatorData(MA_matrix=MA_matrix,CNV_matrix=CNV_matrix,MET_matrix=MET_matrix,PvalueThreshold=PvalueThreshold,RsquareThreshold=RsquareThreshold)
  RegulatorData=RegulatorInfo$RegulatorData
  Alterations = RegulatorInfo$Alterations
  
  # Selecting gene expression data and Clustering
  MA_matrix_Var=geneFiltering('MAD',MA_matrix,VarPercentage)
  MA_matrix_Var=t(scale(t(MA_matrix_Var)))
  if (nrow(MA_matrix_Var)>NrModules){
    KmeansResults=kmeans(MA_matrix_Var,NrModules,iter.max=100)
  } else {
    stop("The number of modules is too large compared to the total number of genes.")
  }
  Clusters=as.numeric(KmeansResults$cluster)
  names(Clusters) <- rownames(MA_matrix_Var)
  
  return(list(MA_matrix_Var=MA_matrix_Var,RegulatorData=RegulatorData,RegulatorAlterations=Alterations,ModuleMembership=Clusters,Parameters=Parameters,NrCores=NrCores))
}

AMARETTO_Run <- function(AMARETTOinit) {
  if (nrow(AMARETTOinit$RegulatorData)==1){
    stop('For cancer ',CancerSite,' only one driver is detected. AMARETTO cannot be run with less than two drivers.\n')
  }
  cat('Running AMARETTO on',length(rownames(AMARETTOinit$MA_matrix_Var)),'genes and',length(colnames(AMARETTOinit$MA_matrix_Var)),'samples.\n')
  cat('\tStopping if less then',0.01*length(rownames(AMARETTOinit$MA_matrix_Var)),'genes reassigned.\n')
  result=AMARETTO_LarsenBased(AMARETTOinit$MA_matrix_Var,AMARETTOinit$ModuleMembership,AMARETTOinit$RegulatorData,AMARETTOinit$Parameters,AMARETTOinit$NrCores)
    
  result$ModuleData=AMARETTO_CreateModuleData(AMARETTOinit,result)
  result$RegulatoryProgramData=AMARETTO_CreateRegulatorPrograms(AMARETTOinit,result)
    
  return(result)
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
  #  result <- list(Beta = regulatoryPrograms$Beta,Error= regulatoryPrograms$error,PreviousClusters = PreviousClusters,AutoRegulationReport = regulatoryPrograms$AutoRegulationReport)
  
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
    
    if (length(which(Clusters==ClusterIDs[i]))==1){
      y = Data[which(Clusters == ClusterIDs[i]),]
    } else {
      y = apply((Data[which(Clusters == ClusterIDs[i]),]),2,mean)
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
    }
    
    # need to do this after the autoregulation, otherwise autoregulation
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
    OldModule = Clusters[i]
    CurrentGeneVector = Data[i,,drop=FALSE]
    Correlations = cor(t(CurrentGeneVector),t(ModuleVectors))
    corr = data.matrix(Correlations,rownames.force = NA)
    MaxCorrelation = max(corr,na.rm=TRUE)
    MaxPosition = which(signif(corr,digits=7) == signif(MaxCorrelation,digits=7))
    MaxPosition = MaxPosition[1] # this is new, to avoid two different reassignements
    
    if (AutoRegulation > 0){ #reassign unrestricted
      if (MaxPosition != OldModule){
        NrReassignGenes = NrReassignGenes + 1
      }
      NewClusters = MaxPosition
      
    }
    else {#only reassign if not a regulator THIS PART IS NOT VALIDATED !!! after update to RegulatorData instead of Regulators referred to as positions
      if (nnzero(rownames(RegulatorData_rownames) %in% GeneNames[i])!=0){
        if (nnzero(which(which(GeneNames %in% rownames(RegulatorData_rownames)) %in% i) %in% which(Beta[MaxPosition,] != 0)) != 0){
          if (MaxPosition != OldModule){
            NrReassignGenes = NrReassignGenes + 1
          }
          NewClusters = MaxPosition
        } else {
          NewClusters = OldModule
        }
      } else {
        if (MaxPosition != OldModule){
          NrReassignGenes = NrReassignGenes + 1
        }
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
        outcome = apply(MA_Data_TestSet[currentClusterGenes,],2,mean)
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
  dimnames(stats) <- list(rownames(stats, do.NULL = FALSE, prefix = "Module_"),
                          c("nrPresReg" ,"nrTotalReg", "nrPresGen", "nrTotGen", "percPresGen",
                            "percWeightPresent", "Rsquare", "RsquareAdjusted","MSE"))
  return(stats)
}

AMARETTO_VisualizeModule <- function(AMARETTOinit,AMARETTOresults,CNV_matrix,MET_matrix,ModuleNr) {
  # getting the data
  if (ModuleNr>AMARETTOresults$NrModules){
    stop('\tCannot plot Module',ModuleNr,'since the total number of modules is',AMARETTOresults$N,'.\n')
  }
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
  Alterations[which(Alterations==1)] <- rep("Methylation aberrations",length(which(Alterations==1)))
  Alterations[which(Alterations=="10")] <- rep("Copy number alterations",length(which(Alterations=="10")))
  Alterations[which(Alterations=="11")] <- rep("Both methylation and copy number alterations",length(which(Alterations=="11")))
  Alterations[which(Alterations=="0")] <- rep(" ",length(which(Alterations=="0")))
  Alterations <- data.frame(Alterations)
    
  ha = HeatmapAnnotation(df = Alterations, col = list(Alterations= c("Copy number alterations"="gray","Methylation aberrations"="black","Both methylation and copy number alterations"="bisque4"," "="white","Amplified gene"="yellow","Hyper-methylated gene" = "yellow","Deleted gene "="cornflowerblue","Hypo-methylated gene" = "cornflowerblue")),
                           which = "row", width = unit(1, "cm"),name="")
    
  CNVreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$CNV),currentRegulators)
  METreg <- intersect(rownames(AMARETTOinit$RegulatorAlterations$MET),currentRegulators)
  if (length(CNVreg)>0){
    CNVData <- matrix(0,nrow=length(CNVreg),ncol=ncol(ModuleData))
    colnames(CNVData) <- colnames(ModuleData)[SampleClustering$order]
    rownames(CNVData) <- CNVreg
    CNVData[,colnames(CNV_matrix)] <- CNV_matrix[rownames(CNVData),]  
    CNVData2 = CNVData
    CNVData2[which(CNVData>0)] <- "Amplified"  # amplification
    CNVData2[which(CNVData<0)] <- "Deleted"  # deletion
    CNVData2[which(CNVData==0)] <- " " # nothing      
 
    if (length(METreg)>0){  
      METData <- matrix(0,nrow=length(METreg),ncol=ncol(ModuleData))
      colnames(METData) <- colnames(ModuleData)[SampleClustering$order]
      rownames(METData) <- METreg
      METData[,colnames(MET_matrix)] <- MET_matrix[rownames(METData),]  
      METData2 = METData
      METData2[which(METData>0)] <- "Hyper-methylated"  # hyper
      METData2[which(METData<0)] <- "Hypo-methylated"  # hypo
      METData2[which(METData==0)] <- " " # nothing      
      Genes <- data.frame(t(METData2),No=rep(" ",ncol(ModuleData)),t(CNVData2))
    } else {
      Genes <- data.frame(No=rep(" ",ncol(ModuleData)),t(CNVData2))
    }
  } else {
    METData <- matrix(0,nrow=length(METreg),ncol=ncol(ModuleData))
    colnames(METData) <- colnames(ModuleData)[SampleClustering$order]
    rownames(METData) <- METreg
    METData[,colnames(MET_matrix)] <- MET_matrix[rownames(METData),]  
    METData2 = METData
    METData2[which(METData>0)] <- "Hyper-methylated"  # hyper
    METData2[which(METData<0)] <- "Hypo-methylated"  # hypo
    METData2[which(METData==0)] <- " " # nothing      
    Genes <- data.frame(t(METData2),No=rep(" ",ncol(ModuleData)))
  }
  ColAnnotation <- c("Hyper-methylated" = "yellow", "Hypo-methylated" = "cornflowerblue"," "="white","Amplified"="yellow","Deleted"="cornflowerblue")
  ColAnnotation <- rep(list(ColAnnotation),ncol(Genes))
  names(ColAnnotation) <- colnames(Genes)
  haRow = HeatmapAnnotation(df=Genes ,name="test",
                              col = ColAnnotation,which="column",show_legend=FALSE)
    
  # plotting
  heatmap <- Heatmap(ClustCombinedData, name = "Gene expression", column_title = paste('Module',ModuleNr), cluster_rows=FALSE,cluster_columns=FALSE,show_column_dend=FALSE,show_column_names=FALSE,row_names_gp=gpar(col=c(rep("white",nrow(ModuleData)),rep("black",nrow(RegulatorData))),fontsize=10),
                    column_title_gp = gpar(fontsize = 20, fontface = "bold"),split=c(rep("Module Genes",nrow(ModuleData)),rep("Regulators",nrow(RegulatorData))),gap = unit(5, "mm"),
                       #  col=colorRamp2(c(-max(abs(CombinedData)), 0, max(abs(CombinedData))), c("green", "black", "red")))
                      col=colorRamp2(c(-max(abs(ClustCombinedData)), 0, max(abs(ClustCombinedData))), c("green", "black", "red")),heatmap_legend_param = list(color_bar = "continuous"),bottom_annotation=haRow)
  draw(heatmap + ha)
    
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

AMARETTO_CreateModuleData <- function(AMARETTOinit,AMARETTOresults) {
  
  # creating the module data based on the average expression of the genes in the module.
  ModuleData=matrix(0,AMARETTOresults$NrModules,length(colnames(AMARETTOinit$MA_matrix_Var)))
  rownames(ModuleData)=rownames(AMARETTOresults$AutoRegulationReport)
  colnames(ModuleData)=colnames(AMARETTOinit$MA_matrix_Var)
  for (ModuleNr in 1:AMARETTOresults$NrModules) {
    # getting the data
    currentModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership[,1]==ModuleNr,]
    ModuleData[ModuleNr,]=colMeans(currentModuleData)
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

CreateRegulatorData <- function(MA_matrix,CNV_matrix,MET_matrix,PvalueThreshold=0.001,RsquareThreshold=0.1) {
  
  # First intersecting data sets
  GeneVariances=rowVars(CNV_matrix)
  CNV_matrix=CNV_matrix[GeneVariances>=0.0001,]
  GeneVariances=rowVars(MET_matrix)
  MET_matrix=MET_matrix[GeneVariances>=0.0001,]
  OverlapGenesCNV = Reduce(intersect,list(rownames(CNV_matrix),rownames(MA_matrix)))
  OverlapGenesMET = Reduce(intersect,list(rownames(MA_matrix),rownames(MET_matrix)))
  if (length(OverlapGenesCNV)==0 && length(OverlapGenesMET)==0){
    stop("AMARETTO cannot be run. Data sets do not share genes.")
  }  
  OverlapSamplesCNV = Reduce(intersect,list(colnames(CNV_matrix),colnames(MA_matrix)))
  OverlapSamplesMET = Reduce(intersect,list(colnames(MA_matrix),colnames(MET_matrix)))
  if (length(OverlapSamplesCNV)==0 && length(OverlapSamplesMET)==0){
    stop("AMARETTO cannot be run. Data sets do not share samples.")
  }  
  
  # Based on copy number, we have the following regulators
  CNV_matrix2 <- mat.or.vec(length(OverlapGenesCNV),length(OverlapSamplesCNV))
  MA_matrix2 <- mat.or.vec(length(OverlapGenesCNV),length(OverlapSamplesCNV))
  rownames(CNV_matrix2) = rownames(MA_matrix2) = OverlapGenesCNV
  colnames(CNV_matrix2) = colnames(MA_matrix2) = OverlapSamplesCNV
  CNV_matrix2[OverlapGenesCNV,OverlapSamplesCNV] = CNV_matrix[OverlapGenesCNV,OverlapSamplesCNV]
  MA_matrix2[OverlapGenesCNV,OverlapSamplesCNV] = MA_matrix[OverlapGenesCNV,OverlapSamplesCNV]
  
  if (length(OverlapGenesCNV)>0 && length(OverlapSamplesCNV)>0){
    CNV_matrix2=FindTranscriptionallyPredictive_CNV(MA_matrix2,CNV_matrix2,PvalueThreshold=PvalueThreshold,RsquareThreshold=RsquareThreshold)
  }
  cat('\tFound',length(rownames(CNV_matrix2)),'CNV driver genes.\n')

  # based on MethylMix we have the following regulators
  MET_matrix = MET_matrix[OverlapGenesMET,OverlapSamplesMET]
  cat('\tFound',length(rownames(MET_matrix)),'MethylMix driver genes.\n')

  # set of drivers
  Drivers <- unique(c(rownames(CNV_matrix2),rownames(MET_matrix)))
  
  # regulator data
  RegulatorData=MA_matrix[Drivers,]
  RegulatorData=t(scale(t(RegulatorData)))
  
  # alterations of all these drivers
  if (nrow(MET_matrix)>0){
    MET_aberrations <- matrix(0,ncol=3,nrow=nrow(MET_matrix))
    colnames(MET_aberrations) <- c("Hypo-methylated","No_change","Hyper-methylated")
    rownames(MET_aberrations) <- rownames(MET_matrix)
    MET_aberrations[,1] <- rowSums(MET_matrix>0)/ncol(MET_matrix)
    MET_aberrations[,2] <- rowSums(MET_matrix<0)/ncol(MET_matrix)
    MET_aberrations[,3] <- rowSums(MET_matrix==0)/ncol(MET_matrix)
  } else {
    MET_aberrations <- NA
  }
  
  if (nrow(CNV_matrix2)>0){
    CNV_alterations <- matrix(0,ncol=3,nrow=nrow(CNV_matrix2))
    colnames(CNV_alterations) <- c("Amplification","No_change","Deletion")
    rownames(CNV_alterations) <- rownames(CNV_matrix2)
    CNV_alterations[,1] <- rowSums(CNV_matrix2>0)/ncol(CNV_matrix2)
    CNV_alterations[,2] <- rowSums(CNV_matrix2<0)/ncol(CNV_matrix2)
    CNV_alterations[,3] <- rowSums(CNV_matrix2==0)/ncol(CNV_matrix2)
  } else {
    CNV_alterations <- NA
  }
  Alterations <- matrix(0,nrow=length(Drivers),ncol=2)
  colnames(Alterations) <- c("CNV","MET")
  rownames(Alterations) <- Drivers
  Alterations[which(Drivers%in% rownames(CNV_matrix2)),1] <- rep(1,length(which(Drivers%in% rownames(CNV_matrix2))))
  Alterations[which(Drivers%in% rownames(MET_matrix)),2] <- rep(1,length(which(Drivers%in% rownames(MET_matrix))))
  
  Alterations <- list(MET=MET_aberrations,CNV=CNV_alterations,Summary=Alterations)
  
  return(list(RegulatorData=RegulatorData,Alterations=Alterations))
}

FindTranscriptionallyPredictive_CNV<- function(MA_matrix,CNV_matrix,PvalueThreshold=0.001,RsquareThreshold=0.1) {
  
      # Linear modeling
      CNVdrivers=c()
      for (i in 1:length(rownames(CNV_matrix))) {
        res=lm(MA_matrix[i,]~CNV_matrix[i,])
        res.summary=summary(res)
        if (res$coefficients[2]>0 & res.summary$coefficients[2,4]<PvalueThreshold & res.summary$r.squared>RsquareThreshold) {
          CNVdrivers=c(CNVdrivers,rownames(CNV_matrix)[i])
        }
      }
      CNV_matrix=CNV_matrix[CNVdrivers,]
  
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
           if (tmpNrGenes>1){
             MAdata_Filtered=MAdata[SortedGenes[1:tmpNrGenes],]
           } else {
             stop('AMARETTO cannot be run. Variable VarPercentage is too small compared to the number of genes.\n')
           }
         }
  )      
  return(MAdata_Filtered)     
}
