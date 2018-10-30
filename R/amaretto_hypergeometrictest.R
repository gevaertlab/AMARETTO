#' Hyper Geometric Geneset Enrichement Test
#' @param gmtfile gmt with reference gene set
#' @param testgmtfile gmt with gene sets to test
#' @param NrCores number of cores for parallelization
#' @param ref.numb.genes The total number of genes teste, standard equal to 45 956 (MSIGDB standard)
#' @import doparallel


HyperGTestGeneEnrichment<-function(gmtfile,testgmtfile,NrCores,ref.numb.genes=45956){
  
  test.gmt<-readGMT(testgmtfile) # our gmt_file_output_from Amaretto
  gmt.path<-readGMT(gmtfile)  # the hallmarks_and_co2...
  
  ###########################  Parallelizing :
  cluster <- makeCluster(c(rep("localhost", NrCores)), type = "SOCK")
  registerDoParallel(cluster,cores=NrCores)
  
  resultloop<-foreach(j=1:length(test.gmt), .combine='rbind') %do% {
    #print(j)
    foreach(i=1:length(gmt.path),.combine='rbind') %dopar% {
      #print(i)
      l<-length(gmt.path[[i]])
      k<-sum(gmt.path[[i]] %in% test.gmt[[j]])
      m<-ref.numb.genes
      n<-length(test.gmt[[j]])
      p1<-phyper(k-1,l,m-l,n,lower.tail=FALSE)
      
      if (k>0){
        overlapping.genes<-gmt.path[[i]][gmt.path[[i]] %in% test.gmt[[j]]]
        overlapping.genes<-paste(overlapping.genes,collapse = ', ')
        c(Geneset=names(gmt.path[i]),Testset=names(test.gmt[j]),p_value=p1,n_Overlapping=k,Overlapping_genes=overlapping.genes)
      }
    }
  }
  
  stopCluster(cluster)
  resultloop<-as.data.frame(resultloop,stringsAsFactors=FALSE)
  resultloop$p_value<-as.numeric(resultloop$p_value)
  resultloop$n_Overlapping<-as.numeric((resultloop$n_Overlapping))
  resultloop[,"padj"]<-p.adjust(resultloop[,"p_value"],method='BH')
  return(resultloop)
}