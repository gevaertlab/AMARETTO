##################################  HTML Report Functions
amaretto_html_report <- function(AMARETTOinit,AMARETTOresults,CNV_matrix,MET_matrix,VarPercentage,hyper_geo_test_bool=TRUE,n_cluster=AMARETTOinit$NrCores,wordcloud_bool=FALSE,hyper_geo_refence_name='H.C2CP.genesets.gmt',output_address='./',output_name='cancer')
{
  
    #You MUST change these two addresses to the address of the corresponding files. 
    hyper_geo_reference_geneset_address=paste('./hyper_geo_test/',hyper_geo_refence_name,sep='')
    all_gene_address="./hyper_geo_test/all_genes.txt"

    set.seed(1234)
   
    
    #required packages 
    list_of_packages<-c("AMARETTO","svglite","R2HTML","GSEABase","rstudioapi","tm","SnowballC","wordcloud","RColorBrewer","foreach","doParallel","tibble","tidyverse")
    lapply(list_of_packages, require, character.only = TRUE)
    ############################################################
    NrModules<-AMARETTOresults$NrModules
    ##################################################################################################################################################################
    #Create necessary folders :
    ##################################################################################################################################################################
    unlink(paste(output_address,"report_htm/*",sep=''))
    unlink(paste(output_address,"report_html/htmls/*",sep=''))
    unlink(paste(output_address,"report_html/htmls/images/*",sep=''))
    unlink(paste(output_address,"report_html/htmls/amaretto_data/*",sep=''))
    unlink(paste(output_address,"report_html/htmls/modules_data/*",sep=''))
    unlink(paste(output_address,"report_html/htmls/hyper_geo_data/*",sep=''))
    unlink(paste(output_address,"report_html/htmls/hyper_geo_data/module_hyper_geo_test/*",sep=''))
    
    dir.create(paste(output_address,"report_html",sep=''))
    dir.create(paste(output_address,"report_html/htmls",sep=''))
    dir.create(paste(output_address,"report_html/htmls/images",sep=''))
    dir.create(paste(output_address,"report_html/htmls/amaretto_data",sep=''))
    dir.create(paste(output_address,"report_html/htmls/modules_data",sep=''))
    dir.create(paste(output_address,"report_html/htmls/hyper_geo_data",sep=''))
    dir.create(paste(output_address,"report_html/htmls/hyper_geo_data/module_hyper_geo_test",sep=''))
    ########################################################
    # Save images of all the modules
    ########################################################   
    
    save(AMARETTOresults, file=paste(output_address,'report_html/',"amarettoResults","_",output_name,".RData",sep=""))
    
    address1=paste(output_address,"report_html",sep="")
    #address2=paste(output_address,"htmls",sep="")
    #address3=paste(output_address,"htmls/images",sep="")
    
    for (ModuleNr in 1:NrModules )
    {
      html_address=paste(output_address,"report_html","/htmls/images","/module",as.character(ModuleNr),".svg",sep="")
      svglite(file =html_address)
      AMARETTO_VisualizeModule(AMARETTOinit, AMARETTOresults=AMARETTOresults, CNV_matrix, MET_matrix, ModuleNr=ModuleNr) 
      dev.off()
    }
    ##############################################################################
    # Create HTMLs for each module
    ##############################################################################   
    if (hyper_geo_test_bool)
    {
        processTCGA_modules(AMARETTOinit,AMARETTOresults)
        b<- HyperGTestGeneEnrichment(hyper_geo_reference_geneset_address, "./TCGA_modules_target_only.gmt", "./output.txt",n_cluster,all_gene_address=all_gene_address,show.overlapping.genes=TRUE)
        df3<-read.delim("./output.genes.txt", header=TRUE, sep="\t")
    }
    
    ############  Removing files
    fn1 <- "./output.genes.txt"
    fn2<-"./TCGA_modules_target_only.gmt"
    if (file.exists(fn1)) file.remove(fn1)
    if (file.exists(fn2)) file.remove(fn2)
    ############

    number_of_significant_gene_overlappings<-c()
    number_of_significant_gene_overlappings_all<-c()
    for (ModuleNr in 1:NrModules )
    {
      module_name=paste("module",as.character(ModuleNr),sep="")
      ModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
      currentRegulators = AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
      RegulatorData=AMARETTOinit$RegulatorData[currentRegulators,]
      module_regulators_weights=AMARETTOresults$RegulatoryPrograms[ModuleNr,][which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
      module_regulators_weights<-data.frame(module_regulators_weights)
      positiveRegulators=AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] > 0)]
      negetiveRegulators=AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] < 0)]
      ModuleGenes=rownames(ModuleData)
      RegulatoryGenes=rownames(RegulatorData)
      module_all_genes_data <- rbind(ModuleData, RegulatorData)
      module_all_genes_data <-module_all_genes_data[order(rownames(module_all_genes_data)),]
      module_all_genes_data <- unique(module_all_genes_data)
      module_annotations<-create_gene_annotations(module_all_genes_data,ModuleGenes,module_regulators_weights)

      if (hyper_geo_test_bool)
      {
        ####################### Hyper Geometric Significance
        module_name2=paste("Module_",as.character(ModuleNr),sep="")
        df3<-df3[order(df3$p.value),]
        filter_indexes_all<-(df3$Testset==module_name2) 
        filter_indexes<-((df3$Testset==module_name2) & (df3$p.value<0.05)) & (df3$n.Overlapping>1)
  
        gene_descriptions<-df3$Description[filter_indexes]
        gene_descriptions_all<-df3$Description[filter_indexes_all]
        
        gene_names<-df3$Geneset[filter_indexes]
        gene_names_all<-df3$Geneset[filter_indexes_all]
        
        overlapping_gene_names<-df3$Overlapping.genes[filter_indexes]
        overlapping_gene_names_all<-df3$Overlapping.genes[filter_indexes_all]
        
        number_overlappings<-df3$n.Overlapping[filter_indexes]
        number_overlappings_all<-df3$n.Overlapping[filter_indexes_all]
        
        p_values<-df3$p.value[filter_indexes]
        p_values_all<-df3$p.value[filter_indexes_all]
        
        q_values<-df3$q.value[filter_indexes]
        q_values_all<-df3$q.value[filter_indexes_all]
        
        number_of_significant_gene_overlappings<-c(number_of_significant_gene_overlappings,length(gene_names))
        number_of_significant_gene_overlappings_all<-c(number_of_significant_gene_overlappings_all,length(gene_names_all))
        
        ####################### WordCloud making
        if (wordcloud_bool)
          {
              mmm<-gene_descriptions
              mm<-as.vector(unique(mmm))
              
              descriptions=""
              for (var in mm) 
                {
                  #descriptions=strsplit(descriptions,",")[[1]]
                  descriptions = paste(descriptions,var,sep=" ")
                  descriptions =gsub(">",",",descriptions)
                }
              if (nchar(descriptions)>0)
                {
                  wordcloud_making(descriptions,module_name2)
                }
           }
      }
      
      address=paste(output_address,"htmls",sep="")
      fname=paste("module",as.character(ModuleNr),sep="")
      tite_page=paste("module",as.character(ModuleNr),sep="")
      graph1=paste("./images","/module",as.character(ModuleNr),".svg",sep = "")
      tmpfic<-HTMLInitFile(paste(output_address,"report_html/htmls/",sep=""),filename=fname,Title = tite_page,CSSFile="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css",useGrid = FALSE,useLaTeX=FALSE,HTMLframe=FALSE)
      ####### CSS ####
      bootstrap1='<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">'
      bootstrap2='<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>'
      bootstrap3='<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>'
      
      HTML(bootstrap1,file=tmpfic)
      HTML(bootstrap2,file=tmpfic)
      HTML(bootstrap3,file=tmpfic)
      
      HTML('<div class="container-fluid">',file=tmpfic)
      
      ################
      HTML('<br />')
      HTML(paste("<h1 class='text-center text-primary'>",
                 '<span class="label label-default">',
                 "AMARETTO Module ",as.character(ModuleNr),' Report',
                 '</span>',
                 '</h1>'),file=tmpfic)
      
      
      HTML('<br /><br /><br />')
      HTML("<h2 class='text-center text-primary'> Regulatory Module Heatmap</h2>",file=tmpfic)
      HTML('<p></p>')
      HTML('<p class="character"></p>')
      HTML('<div class="row">',file=tmpfic)
        HTML('<div style="text-align: center">',file=tmpfic)
              text1<-paste('<a href=',graph1," ",'download> Download Module ',as.character(ModuleNr),' Heatmap.svg</a>',sep='')
              HTML(text1,file=tmpfic)
        HTML('</div>',file=tmpfic)
      HTML('</div>',file=tmpfic)
      
      
      
      
      HTML('<div class="row">',file=tmpfic)
      
        HTML('<div class="col-sm-1">',file=tmpfic)
        HTML('</div>',file=tmpfic)
        
        HTML('<div class="col-sm-10">',file=tmpfic)
            statement<-paste('<p class="character"><img src=',graph1,'border=1 width=1000></p>')
            #statement<-paste('<img class="img-fluid" src=',graph1,'>')
            HTML('<div align="center">',file=tmpfic)
            HTML(statement,file=tmpfic)
            
            
            HTML('</div>',file=tmpfic) 
        HTML('</div>',file=tmpfic)
        
        HTML('<div class="col-sm-1">',file=tmpfic)
        HTML('</div>',file=tmpfic)
        
      HTML('</div>',file=tmpfic)
      ########################################################################################################################################
      ########################################################################################################################################
      HTML('<hr class="col-xs-12">')
      
      ModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
      currentRegulators = AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
      
      positiveRegulators=AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] > 0)]
      positiveRegulators_weights=AMARETTOresults$RegulatoryPrograms[ModuleNr,][which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] > 0)]
      positiveReg_df<-data.frame(positiveRegulators=positiveRegulators, positiveRegulators_weights=positiveRegulators_weights)
      positiveReg_df<-positiveReg_df[order(-positiveReg_df$positiveRegulators_weights),]

      negetiveRegulators=AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] < 0)]
      negetiveRegulators_weights=AMARETTOresults$RegulatoryPrograms[ModuleNr,][which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] < 0)]
      negetiveReg_df<-data.frame(negetiveRegulators=negetiveRegulators, negetiveRegulators_weights=negetiveRegulators_weights)
      negetiveReg_df<-negetiveReg_df[order(negetiveReg_df$negetiveRegulators_weights),]

      module_regulators_data=AMARETTOresults$RegulatoryPrograms[ModuleNr,][which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
  
      HTML("<h2 class='text-center text-primary'>Regulatory Module Gene Expression Data </h2>",file=tmpfic)

      all_gene_expression_file_name_save=paste(module_name,"_","data",".gct",sep="")
      all_gene_expression_file_address=paste(output_address,"report_html/htmls/modules_data",'/',all_gene_expression_file_name_save,sep="")
      
      write_gct(module_all_genes_data,all_gene_expression_file_address)
      ModuleData<-round(ModuleData,2)
     
      annotations_file_name_save=paste(module_name,"_","annotations",".gct",sep="")
      annotations_file_address=paste(output_address,'report_html/htmls/modules_data','/',annotations_file_name_save,sep="")
      write_gct(module_annotations,annotations_file_address)
      #############
      HTML('<div class="row">',file=tmpfic)
            HTML('<div style="text-align: center">',file=tmpfic)
            #HTML('<br />',file=tmpfic)
            HTML(paste(
              
              '<a href=', paste('./modules_data','/',all_gene_expression_file_name_save,sep=""),'  download=',paste('./modules_data','/',all_gene_expression_file_name_save,sep=""),'>',
                       ' Download Module ',as.character(ModuleNr),' Gene Expression Data.gct',
                       '</a>',sep=""),file=tmpfic)

            HTML(paste(
                       '<a href=', paste('./modules_data','/',annotations_file_name_save,sep=""),'  download>',
                       ' Download Module ',as.character(ModuleNr),' Gene Annotations.gct ',
                       '</a>',sep=""),file=tmpfic)

            HTML('</div>',file=tmpfic)
       HTML('</div>',file=tmpfic)

      table_command5=
        '
      <table class="table table-hover .table-striped table-bordered">
      <thead>
      <tr>
      <th  class="text-center">Regulator</th>
      <th scope="col" class="text-center">Weights</th>
      </tr>
      </thead>
      <tbody>
      '
      HTML('<div class="row">',file=tmpfic)
      
        HTML('<div class="col-sm-4">',file=tmpfic)
        HTML('</div>',file=tmpfic)
        
        HTML('<div class="col-sm-4">',file=tmpfic)
        HTML(table_command5,file=tmpfic)
        
        for (kk in 1:length(positiveRegulators))
        {
        HTML(paste('<tr>',
                   '<td align="center" class="text-danger">',positiveReg_df$positiveRegulators[kk],'</td>',
                   '<td align="center">',round(positiveReg_df$positiveRegulators_weights[kk],5),'</td>',
                   '</tr>'),file=tmpfic)
        }
        
        if (length(negetiveRegulators)>0)
        {
            for (kk in 1:length(negetiveRegulators))
            {
              HTML(paste('<tr>',
                         '<td align="center" class="text-primary">',negetiveReg_df$negetiveRegulators[kk],'</td>',
                         '<td align="center">',round(negetiveReg_df$negetiveRegulators_weights[kk],5),'</td>',
                         '</tr>'),file=tmpfic)
            }
        }
        HTML('</tbody></table>',file=tmpfic)
        HTML('</div>',file=tmpfic)
        
        HTML('<div class="col-sm-4">',file=tmpfic)
        HTML('</div>',file=tmpfic)
      
      HTML('</div>',file=tmpfic)
      ########################################################################################################################################
      ########################################################################################################################################
      # Hyper Geometric Test Results
            if ((hyper_geo_test_bool)&(length(gene_names)>0))
            {
                
                HTML('<hr class="col-xs-12">')
                HTML("<h2 class='text-center text-primary'> Functional Enrichment Results </h2>",file=tmpfic)
                HTML('<div class="row">',file=tmpfic)

                ###################################  Download Functional Enrchement table ########################################
                module_hypo_table_header<-c('Gene Set Name','Gene Set Description','Number of Genes in Overlap','Names of Genes in Overlap','P-value','FDR q-value')
                module_hypo_table<-c()
                
                for (kk in 1:length(gene_names_all))
                {
                  rr<-c(as.character(gene_names_all[kk]),gsub(">"," ",gene_descriptions_all[kk]),number_overlappings_all[kk],gsub("  ",' ',as.character(overlapping_gene_names_all[kk])),p_values_all[kk],q_values_all[kk])
                  module_hypo_table<-rbind(module_hypo_table,rr)
                }
                
                colnames(module_hypo_table)<-module_hypo_table_header
                #outfile_hypo_tests=paste('./report_html/htmls/hyper_geo_data/module_hyper_geo_test/Module',ModuleNr,'_hypergeometric_test.tsv',sep='')
                #write.table(module_hypo_table,file=outfile_hypo_tests,sep='\t',quote=F,col.names=T,row.names=F)
                
                outfile_hypo_tests=paste(output_address,'report_html/htmls/hyper_geo_data/module_hyper_geo_test/Module',ModuleNr,'_hypergeometric_test.tsv',sep='')
                write_tsv(module_hypo_table,outfile_hypo_tests,row_name=FALSE)
                
                hypo_geo_file_name_save=  paste('module_hyper_geo_test/Module',ModuleNr,'_hypergeometric_test.tsv',sep='')
                
                HTML('<div style="text-align: center">',file=tmpfic)
                HTML(paste(' <a href=', paste('./hyper_geo_data','/',hypo_geo_file_name_save,sep=""),'  download>',
                           ' Download Module ',as.character(ModuleNr), ' Functional Categories Enrichment Analysis.tsv (All functional categories sorted by P-value)',
                           '</a>',sep=""))
                HTML('</div>',file=tmpfic)
                HTML('</div>',file=tmpfic)
                ######################################################################################################
                
                if (wordcloud_bool){
                        if (nchar(descriptions)==0){
                          HTML("<h4 class='text-center text-danger'> Not enough for wordcloud </h4>",file=tmpfic)
                        } 
                        if (nchar(descriptions)>0)
                        {
                          graph2=paste("./images","/",module_name2,"_WordCloud.svg",sep = "")
                          statement<-paste('<img class="img-fluid" src=',graph2,'>')
                          
                          HTML('<div class="row">',file=tmpfic)
                          
                            HTML('<div class="col-sm-2">',file=tmpfic)
                            HTML('</div>',file=tmpfic)
                          
                            HTML('<div class="col-sm-8">',file=tmpfic)
                          
                              HTML('<div align="center">',file=tmpfic)
                                  HTML(statement,file=tmpfic)
                              HTML('</div>',file=tmpfic) 
                              
                            HTML('</div>',file=tmpfic)
                          
                            HTML('<div class="col-sm-2">',file=tmpfic)
                            HTML('</div>',file=tmpfic)
                          
                          HTML('</div>',file=tmpfic)
                        }
                }
                
                ################  Hyper-geo-table
                if (length(gene_names)>0)
                {
                  HTML('<div class="row">',file=tmpfic)
                    HTML('<div class="col-sm-1">',file=tmpfic)
                    HTML('</div>',file=tmpfic)
                  HTML('<div class="col-sm-10">',file=tmpfic)
                  
                      table_command2=
                        '
                      <table class="table table-hover .table-striped table-bordered">
                      <thead>
                      <tr>
                      <th scope="col" class="text-center">GeneSet Name</th>
                      <th scope="col" class="text-center">GeneSet Description</th>
                      <th scope="col" class="text-center">Number of Genes in Overlap</th>
                      <th scope="col" class="text-center">Name of Genes in Overlap</th>
                      <th scope="col" class="text-center">p-value</th>
                      <th scope="col" class="text-center">FDR q-value</th>
                      
                      </tr>
                      </thead>
                      <tbody>
                      '
                      
                      HTML(table_command2,file=tmpfic)
                      for (kk in 1:length(gene_names))
                            {
                            link_command=paste("<a href=http://software.broadinstitute.org/gsea/msigdb/cards/",as.character(gene_names[kk]),".html>",gene_names[kk],'</a>',sep="")
                              HTML(paste('<tr>',
                                         '<td  style ="word-break:break-all;" align="center">',link_command,'</td>',
                                         '<td align="center"  >',gsub(">"," ",gene_descriptions[kk]),'</td>',
                                         '<td align="center" >',number_overlappings[kk],'</td>',
                                         '<td align="center" >',gsub("  ",' ',as.character(overlapping_gene_names[kk])),'</td>',
                                         '<td align="center" >',round(p_values[kk],4),'</td>',
                                         '<td align="center" >',round(q_values[kk],4),'</td>',
                                         '</tr>'),file=tmpfic)
                            }
                      HTML('</tbody></table>',file=tmpfic)
                      HTML("<p> Threshold : P-value <0.05 and min. 2 overlapping genes </p>",file=tmpfic)
                      ################################################################################################################################################
                    HTML('</div>',file=tmpfic)
                    HTML('<div class="col-sm-1">',file=tmpfic)
                    HTML('</div>',file=tmpfic)
                ##################
                HTML('</div>',file=tmpfic)
                }
          }
      HTML('</div>',file=tmpfic)
    }
    ##############################################################################
    #Create the landing page
    ##############################################################################   
    tmpfic<-HTMLInitFile(paste(output_address,"report_html",sep=""),filename="index",Title = "Amaretto Report",CSSFile="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css",useGrid = FALSE,useLaTeX=FALSE,HTMLframe=FALSE)
    bootstrap1='<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">'
    bootstrap2='<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>'
    bootstrap3='<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>'
    
    HTML(bootstrap1,file=tmpfic)
    HTML(bootstrap2,file=tmpfic)
    HTML(bootstrap3,file=tmpfic)
    
    HTML('<div class="container-fluid">',file=tmpfic)
    #######################################  Create the Title  #########################
    HTML('<br />')
    HTML(' <h1 class="text-primary text-center">
         <span class="label label-primary">
         AMARETTO Report 
         </span>
         </h1>',file=tmpfic)
    HTML('<br /><br /><br /><br />')
    
    #######################################  Saving the Files and creating download links  #########################
    write_gct(AMARETTOresults$ModuleData,paste(output_address,'report_html/htmls/amaretto_data/ModuleData_amaretto.gct',sep=''))
    write_gct(AMARETTOresults$ModuleMembership,paste(output_address,'report_html/htmls/amaretto_data/ModuleMembership_amaretto.gct',sep=''))
    write_gct(AMARETTOresults$RegulatoryProgramData,paste(output_address,'report_html/htmls/amaretto_data/RegulatoryProgramData_amaretto.gct',sep=''))
    write_gct(AMARETTOresults$RegulatoryPrograms,paste(output_address,'report_html/htmls/amaretto_data/RegulatoryPrograms_amaretto.gct',sep=''))
    write_tsv(AMARETTOresults$AllGenes,paste(output_address,'report_html/htmls/amaretto_data/AllGenes_amaretto.tsv',sep=''))
    write_tsv(AMARETTOresults$AllRegulators,paste(output_address,'report_html/htmls/amaretto_data/AllRegulators_amaretto.tsv',sep=''))
    write_tsv(AMARETTOresults$NrModules,paste(output_address,'report_html/htmls/amaretto_data/NrModules_amaretto.tsv',sep=''))

    text1<-'<a href="./htmls/amaretto_data/ModuleData_amaretto.gct" download>Download ModuleData_amaretto.gct</a>'
    text2<-'<a href="./htmls/amaretto_data/ModuleMembership_amaretto.gct" download>Download ModuleMembership_amaretto.gct</a>'
    text3<-'<a href="./htmls/amaretto_data/RegulatoryProgramData_amaretto.gct" download>Download RegulatoryProgramData_amaretto.gct</a>'
    text4<-'<a href="./htmls/amaretto_data/RegulatoryPrograms_amaretto.gct" download>Download RegulatoryPrograms_amaretto.gct</a>'
    text5<-'<a href="./htmls/amaretto_data/AllGenes_amaretto.tsv" download>Download AllGenes_amaretto.tsv</a>'
    text6<-'<a href="./htmls/amaretto_data/AllRegulators_amaretto.tsv" download>Download AllRegulators_amaretto.tsv</a>'
    text7<-'<a href="./htmls/amaretto_data/NrModules_amaretto.tsv" download>Download NrModules_amaretto.tsv</a>'
    
    HTML('<h4 class="text-center">AMARETTO output files : </h4>')
    HTML('<br/>')
    HTML('<div class="text-primary text-center">')
      
      HTML(text1,file=tmpfic)
      HTML(text2,file=tmpfic)
      HTML(text3,file=tmpfic)
      HTML(text4,file=tmpfic)
      HTML(text5,file=tmpfic)
      HTML(text6,file=tmpfic)
      HTML(text7,file=tmpfic)
      
    HTML('</div>')
    HTML('<hr class="col-xs-12">')
    HTML('<br /><br />')
    #######################################  Create the table  #########################
    table_command5=
      '
    <table class="table table-hover.table-striped table-bordered">
    <thead>
    <tr>
    <th scope="col"># of samples</th>
    <th scope="col"># of modules</th>
    <th scope="col"> Var-Percentage</th>
    </tr>
    </thead>
    <tbody>
    '
    
    table_command1=
      '
    
    <table class="table table-hover ">
    <thead ">
    <tr>
    <th scope="col" class="text-center">Regulatory Modules</th>
    <th scope="col" class="text-center">Number of Target Genes</th>
    <th scope="col" class="text-center">Number of Regulator Genes</th>
    <th scope="col" class="text-center">Number of Enriched Functional Categories</th>
    </tr>
    </thead>
    <tbody>
    '
    ModuleNr<-1
    ModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
    number_of_samples_expression=length(colnames(ModuleData))
    number_of_samples_copy_number=dim(CNV_matrix)[2]
    number_of_samples_methylations=dim(MET_matrix)[2]
    
    HTML('<div class="col-sm-1">',file=tmpfic)
    HTML('</div>',file=tmpfic)

    HTML('<div class="col-sm-10">',file=tmpfic)
    
    HTML(paste('<p class=".text-success text-left"> <font size="2"> Number of Samples (Expression) = ',as.character(number_of_samples_expression),'</font></p>'))
    
    
    if (!is.null(CNV_matrix))
    {
          HTML(paste('<p class=".text-success text-left"> <font size="2"> Number of Samples (Copy Number) = ',as.character(number_of_samples_copy_number),'</font></p>'))
    }
    if (!is.null(MET_matrix))
    {
    HTML(paste('<p class=".text-success text-left"> <font size="2"> Number of Samples (Methylation) = ',as.character(number_of_samples_methylations),'</font></p>'))
    }
    HTML(paste('<p class=".text-success text-left"> <font size="2"> Number of Genes ( for ',as.character(VarPercentage),'% most variably expressed genes) = ',as.character(length(AMARETTOresults$AllGenes)),'</font></p>',sep=''))
    HTML(paste('<p class=".text-success text-left"> <font size="2"> Number of Regulatory Modules = ',as.character(NrModules),'</font></p>'))
    
    HTML('<br />')
    
    HTML(table_command1,file=tmpfic)
    
    amaretto_result_table_header<-c('Module_No','number_of_target_genes','number_of_regulator_genes','number_of_significant_gene_overlappings')
    amaretto_result_table<-c()
    
        for (ModuleNr in 1:NrModules )
        {
          module_name=paste("module",as.character(ModuleNr),sep="")
          ###################### find module info
          ModuleData=AMARETTOinit$MA_matrix_Var[AMARETTOresults$ModuleMembership==ModuleNr,]
          currentRegulators = AMARETTOresults$AllRegulators[which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
          RegulatorData=AMARETTOinit$RegulatorData[currentRegulators,]
          module_regulators_weights=AMARETTOresults$RegulatoryPrograms[ModuleNr,][which(AMARETTOresults$RegulatoryPrograms[ModuleNr,] != 0)]
          ModuleGenes=rownames(ModuleData)
          RegulatoryGenes=rownames(RegulatorData)
          
          number_of_genes=length(rownames(ModuleData))
          number_of_regulators=length(currentRegulators)
          number_of_samples=length(colnames(ModuleData))
          
          #########################  creating Link for each module ####################
          address='./htmls'
          htmladdress=paste("'",address,"/module",as.character(ModuleNr),".html","'",sep="")
          
          module_space_name=paste("Module",as.character(ModuleNr),sep=" ")
          
          link_command=paste("<a href=",htmladdress,'>',module_space_name,'</a>',sep="")
          ###########################################################################
          HTML('<tr>',file=tmpfic)
          HTML(paste('<td align="center">',link_command,'</td>'),file=tmpfic)
          HTML(paste('<td align="center">',as.character(number_of_genes),'</td>'),file=tmpfic)
          HTML(paste('<td align="center">',as.character(number_of_regulators),'</td>'),file=tmpfic)
          HTML(paste('<td align="center">',as.character( number_of_significant_gene_overlappings[ModuleNr]),'</td>'),file=tmpfic)
         
          HTML('</tr>',file=tmpfic)
          rr<-c(module_name,as.character(number_of_genes),as.character(number_of_regulators),as.character( number_of_significant_gene_overlappings[ModuleNr]))
          amaretto_result_table<-rbind(amaretto_result_table,rr)
        }
    
    colnames(amaretto_result_table)<-amaretto_result_table_header
    #outfile=paste('./report_html/htmls/hyper_geo_data/amaretto','.tsv',sep='')
    #write.table(amaretto_result_table,file=outfile,sep='\t',quote=F,col.names=T,row.names=F)
    
    outfile=paste(output_address,'report_html/htmls/hyper_geo_data/amaretto','.gct',sep='')
    write_gct(amaretto_result_table,outfile)
    
    
    HTML('</tbody></table>',file=tmpfic)
    HTML('</div>',file=tmpfic)
    HTML('<div class="col-sm-1">',file=tmpfic)
    HTML('</div>',file=tmpfic)
    
    HTML('</div>',file=tmpfic)
    # HTMLEndFile()
    #################################################################################[####
    curr_add<-getwd()
    setwd(output_address)
    zip(zipfile = paste(output_address,'Amaretto_HTML_Report',"_",output_name,sep=''), files = paste('report_html/',sep=''))
    setwd(curr_add)
    ###############################
}

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
create_gene_annotations<-function(module_all_genes_data,Module_Genes_names,Module_regulators_weights)
{
  all_genes_names=rownames(module_all_genes_data)
 
  targets_boolean<-c()
  regulators_boolean<-c()
  regulators_weight<-c()

  for (i in 1:length(all_genes_names))
  {
    
    gene_name=all_genes_names[i]

    
    a=0
    b=0
    c=0
    
    
    if (is.element(gene_name, Module_Genes_names))
    {
      
      a<-1
      
    }
    
    
    
    if (is.element(gene_name,  rownames(Module_regulators_weights)))
    {
      b<-1
      c<-Module_regulators_weights$module_regulators_weights[rownames(Module_regulators_weights)==gene_name]
    }
    
    targets_boolean<-c(targets_boolean,a)
    regulators_boolean<-c(regulators_boolean,b)
    regulators_weight<-c(regulators_weight,c)
  }
  
  df=data.frame(targets_boolean,regulators_boolean,regulators_weight,stringsAsFactors=FALSE)
  
  colnames(df) <- c("targets_boolean", "regulators_boolean","regulators_weight")
  rownames(df)<-all_genes_names

  return(df)
}


write_tsv<-function(data_in,file_address,row_name=TRUE){
  write.table(data_in,file=file_address,sep='\t',quote=F,col.names=T,row.names=row_name)
}


write_gct<-function(data_in,file_address)
{
  
  write('#1.2',file=file_address)
  
  out<-c()
  n_col=ncol(data_in)
  n_row=nrow(data_in)
  second_row<-paste(as.character(n_row),as.character(n_col),sep="\t")
  write(second_row,file=file_address,append = TRUE)
  
  
  third_row=c("Name","Description",colnames(data_in))

  new_df<-data.frame(matrix(,nrow=(n_row+1), ncol=(n_col+2)),stringsAsFactors = FALSE)
  new_df[1,]<-third_row
  
  new_df[2:(n_row+1),1]<-rownames(data_in)
  new_df[2:(n_row+1),2]<-rownames(data_in)
  new_df[2:(n_row+1),3:(n_col+2)]<-data_in
  
  write.table(new_df,file = file_address, sep="\t",quote = FALSE, col.names=FALSE, row.names=FALSE,append = TRUE)

}
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
read_gct<-function(file_address){
  data_fr=read.delim(file_address, skip = 2,sep = '\t',header=TRUE,row.names = 1)
  data_fr <- subset(data_fr, select = -c(Description))
}
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
wordcloud_making<-function(text1,file_name){
  
  save_file_name=paste(file_name,"_WordCloud",sep="")
  # # Install
  # install.packages("tm")  # for text mining
  # install.packages("SnowballC") # for text stemming
  # install.packages("wordcloud") # word-cloud generator 
  # install.packages("RColorBrewer") # color palettes

  #require("tm")
  #require("SnowballC")
  #require("wordcloud")
  #require("RColorBrewer")
  
  text1<-as.character(text1)
  Encoding(text1) <- 'UTF-8'
  # Load the data as a corpus
  doc_ids <- c(1)
  df <- data.frame(doc_id = doc_ids, text = text1, stringsAsFactors = FALSE)
  docs <- Corpus(DataframeSource(df))
  inspect(docs)
  toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  docs <- tm_map(docs, toSpace, "/")
  docs <- tm_map(docs, toSpace, "@")
  docs <- tm_map(docs, toSpace, "\\|")
  
  # Convert the text to lower case
  docs <- tm_map(docs, content_transformer(tolower))
  # Remove numbers
  docs <- tm_map(docs, removeNumbers)
  # Remove english common stopwords
  docs <- tm_map(docs, removeWords, stopwords("english"))
  # Remove your own stop word
  # specify your stopwords as a character vector
  docs <- tm_map(docs, removeWords, c("blabla1", "blabla2","genes","response","involved","events","signaling","encoding","pathway")) 
  # Remove punctuations
  docs <- tm_map(docs, removePunctuation)
  # Eliminate extra white spaces
  docs <- tm_map(docs, stripWhitespace)
  # Text stemminghe
  # docs <- tm_map(docs, stemDocument)
  dtm <- TermDocumentMatrix(docs)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m),decreasing=TRUE)
  d <- data.frame(word = names(v),freq=v)
  head(d, 10)
  
  set.seed(1234)
  
  address=paste(output_address,"report_html/htmls/images/",save_file_name,".svg",sep="")
  svglite(file=address,width = 10, height = 6)
  wordcloud(words = d$word, freq = d$freq, min.freq = 1,
            max.words=200, random.order=FALSE, rot.per=0.35, 
            colors=brewer.pal(8, "Dark2"))
  dev.off()
}
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
processTCGA_modules <- function(AMARETTOinit,AMARETTOresults){
  
  
  all_genes<-unique(c(AMARETTOresults$AllGenes,AMARETTOresults$AllRegulators))
  annotation<-data.frame(Gene=all_genes)
  annotation$Regulator<-annotation$Gene%in%AMARETTOresults$AllRegulators
  
  all_data<-unique(rbind(AMARETTOinit$MA_matrix_Var,AMARETTOinit$RegulatorData))
  
  
  rp<-t(AMARETTOresults$RegulatoryPrograms)
  ma<-AMARETTOresults$ModuleMembership
  
  module_weights<-matrix(0,nrow=length(all_genes),ncol=dim(rp)[2])
  rownames(module_weights)<-all_genes
  colnames(module_weights)<-colnames(rp)
  for(i in 1:dim(module_weights)[1]){
    g<-rownames(module_weights)[i]
    w<-rep(0,dim(rp)[2])
    if(g %in% rownames(rp))
      w<-w+rp[g,]
    if(g %in% rownames(ma))
      w[ma[g,]]<-1+w[ma[g,]]
    module_weights[g,]<-w
  }
  
  annotation<-cbind(annotation,module_weights)
  
  in_module<-module_weights!=0
  colnames(in_module)<-unlist(lapply(colnames(in_module),paste0,"_B"))
  annotation<-cbind(annotation,in_module)
  
  
  gmt_file="./TCGA_modules_target_only.gmt"
  for(i in 1:dim(module_weights)[2]){
    genes<-rownames(in_module)[which(module_weights[,i]==1)]
    if(i==1){
      write(c(colnames(module_weights)[i],colnames(module_weights)[i],genes),file=gmt_file,sep = "\t",ncolumns=length(genes)+2)
    }else{
      write(c(colnames(module_weights)[i],colnames(module_weights)[i],genes),file=gmt_file,sep = "\t",ncolumns=length(genes)+2,append = T)
    }
  }
  
  aa=c(1)
  return(aa)
  
}
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
HyperGTestGeneEnrichment<-function(gmtfile,testgmtfile,outfile,n_cluster,all_gene_address,show.overlapping.genes=FALSE,filter.genes=TRUE,show.unrecognized=FALSE)
  
{
  #library(GSEABase)
  #require(plyr)

  
  
  
  # 
  # ## Help section
  # if("--help" %in% args) {
  #   cat("
  #       The R Script HyperGTestGeneEnrichment is used to perform the hypergeometric test. 
  #  
  #       Arguments:
  #       1) the signature collection
  #       2) the genelist
  #       3) the output file
  #       4) Default FALSE: output the overlapping genes
  #       5) Default TRUE: only count the genes in all_genes.txt (useful for e.g. mouse signatures) 
  #       6) Default TRUE: output genes that are not recognized
  #       --help              - print this text
  #  
  #       Example:
  #       Rscript HyperGTestGeneEnrichment.R h.all.v5.0.symbols.gmt my_genesignatures.gmt output.txt\n\n")
  #   
  #   q(save="no")
  # }
  
  
  #MSigDB 45956
  ref.num<-45956
  
  all_genes<-scan(all_gene_address, what="", sep="\n")
  #ref.num<-length(all_genes)
  
  # gmtfile<-args[1]     #signature collection
  # testgmtfile<-args[2] #genelist
  # outfile<-args[3]     #output file
  # show.overlapping.genes<-args[4]
  if(show.overlapping.genes){
    
  }
  # filter.genes<-as.logical(args[5])
  # show.unrecognized<-as.logical(args[6])
  
  test.gmt<-readGMT(testgmtfile) # our gmt_file_output_from Amaretto
  gmt.path<-readGMT(gmtfile)  # the hallmarks_and_co2...
  out<-c()
  out.genes<-c()
  
  
  if(filter.genes && show.unrecognized){
    genes.in.test.gmt<-unlist(test.gmt$genesets)
    if(sum(!genes.in.test.gmt %in% all_genes)>0) {
      cat("The following genes are not recognized: ",genes.in.test.gmt[!genes.in.test.gmt %in% all_genes],"\n")
    }
  }
  
  
  ###########################  Parallelizing :
  n.cluster <- n_cluster
  cl2 <- makeCluster(c(rep("localhost", n.cluster)), type = "SOCK")
  cluster = cl2
  #library(foreach)
  #library(doParallel)
  registerDoParallel(cores=n.cluster )
  getDoParWorkers()
  
  omidindex=0
  
  resultsss<-c()
  
  resultloop<-foreach(j=1:length(test.gmt$genesets), .combine='rbind') %do%
    
  {
    print(j)
    foreach(i=1:length(gmt.path$genesets),.combine='rbind') %dopar% {
      
      
      set.num<-length(gmt.path$genesets[[i]])
      k<-sum(gmt.path$genesets[[i]] %in% test.gmt$genesets[[j]])
      l<-set.num
      m<-ref.num
      if(filter.genes){
        n<-sum(test.gmt$genesets[[j]] %in% all_genes)
      }else{
        n<-length(test.gmt$genesets[[j]])
      }
      p1<-phyper(k-1,l,m-l,n,lower.tail=FALSE)
     
      
      overlapping.genes<-gmt.path$genesets[[i]][gmt.path$genesets[[i]] %in% test.gmt$genesets[[j]]]
      
      
      overlapping.genes<-gsub('\t',',',as.character(overlapping.genes))
      overlapping.genes<-gsub('  ',',',overlapping.genes)
      overlapping.genes<-gsub(' ',',',overlapping.genes)
      overlapping.genes<-paste(overlapping.genes,collapse = ', ')

      if (k>0)
      {
        c(GENESETNAME=gmt.path$geneset.names[i],GENSETDESCRIPTION=gmt.path$geneset.descriptions[[i]],TESTSETNAME=test.gmt$geneset.names[[j]],p_value=p1,n_Overlapping=k,Overlapping_genes=overlapping.genes)
      }
      
      
    }
    
  }

  resultsss<-resultloop
  
  pp<-resultsss[,4]
  pp.adj<-p.adjust(pp,method='BH')
  resultsss<-cbind(resultsss,pp.adj)
  
  stopCluster(cl2)
 
  outfile.genes<-paste(tools::file_path_sans_ext(outfile),'.genes.',tools::file_ext(outfile),sep="")
  col.names2<-c('Geneset','Description','Testset','p-value','n-Overlapping','Overlapping genes','q-value')
  colnames(resultsss)<-col.names2
  write.table(resultsss,file=outfile.genes,sep='\t',quote=F,col.names=T,row.names=F)
  
  return(c(0))
}
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
HyperGTestGeneEnrichment_serial<-function(gmtfile,testgmtfile,outfile,show.overlapping.genes=FALSE,filter.genes=TRUE,show.unrecognized=FALSE )
  
{
  #library(GSEABase)
  #require(plyr)

  # 
  # ## Help section
  # if("--help" %in% args) {
  #   cat("
  #       The R Script HyperGTestGeneEnrichment is used to perform the hypergeometric test. 
  #  
  #       Arguments:
  #       1) the signature collection
  #       2) the genelist
  #       3) the output file
  #       4) Default FALSE: output the overlapping genes
  #       5) Default TRUE: only count the genes in all_genes.txt (useful for e.g. mouse signatures) 
  #       6) Default TRUE: output genes that are not recognized
  #       --help              - print this text
  #  
  #       Example:
  #       Rscript HyperGTestGeneEnrichment.R h.all.v5.0.symbols.gmt my_genesignatures.gmt output.txt\n\n")
  #   
  #   q(save="no")
  # }
  
  
  #MSigDB 45956
  ref.num<-45956
  
  all_genes<-scan(all_gene_address, what="", sep="\n")

  if(show.overlapping.genes){
    
  }

  test.gmt<-readGMT(testgmtfile)
  gmt.path<-readGMT(gmtfile)
  out<-c()
  out.genes<-c()
  
  
  if(filter.genes && show.unrecognized){
    genes.in.test.gmt<-unlist(test.gmt$genesets)
    if(sum(!genes.in.test.gmt %in% all_genes)>0) {
      cat("The following genes are not recognized: ",genes.in.test.gmt[!genes.in.test.gmt %in% all_genes],"\n")
    }
  }
  

  for(i in 1:length(gmt.path$genesets)){
    
    
    for(j in 1:length(test.gmt$genesets)){
      
      
      
      set.num<-length(gmt.path$genesets[[i]])
      k<-sum(gmt.path$genesets[[i]] %in% test.gmt$genesets[[j]])
      l<-set.num
      m<-ref.num
      if(filter.genes){
        n<-sum(test.gmt$genesets[[j]] %in% all_genes)
      }else{
        n<-length(test.gmt$genesets[[j]])
      }
      p1<-phyper(k-1,l,m-l,n,lower.tail=FALSE)
      r<-c(gmt.path$geneset.names[i],gmt.path$geneset.descriptions[[i]],test.gmt$geneset.names[[j]],ref.num,set.num,n,k,p1)
      out<-rbind(out,r)
      if(show.overlapping.genes){
        overlapping.genes<-gmt.path$genesets[[i]][gmt.path$genesets[[i]] %in% test.gmt$genesets[[j]]]
        
        
        overlapping.genes<-gsub('\t',',',as.character(overlapping.genes))
        overlapping.genes<-gsub('  ',',',overlapping.genes)
        overlapping.genes<-gsub(' ',',',overlapping.genes)
        overlapping.genes<-paste(overlapping.genes,collapse = ', ')
        #if(!identical(overlapping.genes, character(0))){print(overlapping.genes)}
        
        if (as.numeric(k)>0){
          rr<-c(gmt.path$geneset.names[i],gmt.path$geneset.descriptions[[i]],test.gmt$geneset.names[[j]],p1,k,overlapping.genes)
          out.genes<-rbind(out.genes,rr)
        }

      }
      
    }
  }
  
  
  pp<-out.genes[,4]
  pp.adj<-p.adjust(pp,method='BH')
  out.genes<-cbind(out.genes,pp.adj)

  
  outfile.genes<-paste(tools::file_path_sans_ext(outfile),'.genes.',tools::file_ext(outfile),sep="")
  
  
  
  col.names2<-c('Geneset','Description','Testset','p-value','n-Overlapping','Overlapping genes','q-value')
  colnames(out.genes)<-col.names2

  write.table(out.genes,file=outfile.genes,sep='\t',quote=F,col.names=T,row.names=F)
  
  return(c(0))
}
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
readGMT<-function (filename) 
{
  a = scan(filename, what = list("", ""), sep = "\t", quote = NULL, 
           fill = T, flush = T, multi.line = F,strip.white = TRUE)
  geneset.names = a[1][[1]]
  geneset.descriptions = a[2][[1]]
  dd = scan(filename, what = "", sep = "\t", quote = NULL,strip.white = TRUE)
  nn = length(geneset.names)
  n = length(dd)
  ox = rep(NA, nn)
  ii = 1
  for (i in 1:nn) {
    #print(i)
    while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] != 
                                           geneset.descriptions[i])) {
      ii = ii + 1
    }
    ox[i] = ii
    ii = ii + 1
  }
  genesets = vector("list", nn)
  
  for (i in 1:(nn - 1)) {
    # print(i, fill = T)
    if(nn>1){
      i1 = ox[i] + 2
      i2 = ox[i + 1] - 1
      geneset.descriptions[i] = dd[ox[i] + 1]
      current_geneset = dd[i1:i2]
      genesets[[i]] = toupper(current_geneset[current_geneset!=""])
    }else{
      current_geneset = dd[3:n]
      genesets[[1]] = toupper(current_geneset[current_geneset!=""])
    }
  }
  
  geneset.descriptions[nn] = dd[ox[nn] + 1]
  current_geneset = dd[(ox[nn] + 2):n]
  genesets[[nn]] = toupper(current_geneset[current_geneset!=""])
  out = list(genesets = genesets, geneset.names = geneset.names, 
             geneset.descriptions = geneset.descriptions)
  class(out) = "GSA.genesets"
  return(out)
}
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################




