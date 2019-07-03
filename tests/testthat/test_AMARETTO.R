library(AMARETTO)
context("AMARETTO output data objects testing")

data("ProcessedDataLIHC")
NrModules<-10
VarPercentage<-50
AMARETTOinit <- AMARETTO_Initialize(ProcessedDataLIHC,
                                    NrModules = NrModules,
                                    VarPercentage = VarPercentage,
                                    random_seeds = c(42,42))

tryCatch( { AMARETTOresults <- AMARETTO_Run(AMARETTOinit);}
          , warning = function(w) { AMARETTOresults <- AMARETTO_Run(AMARETTOinit) })

testthat::test_that("Checking AMARETTO_init object if it is in decent shape",{
  expect_equal(nrow(AMARETTOinit$MA_matrix_Var)>1,TRUE)
  expect_equal(ncol(AMARETTOinit$MA_matrix_Var)>1,TRUE)
  expect_equal(sum(AMARETTOinit$MA_matrix_Var!=0)>0,TRUE)
  expect_equal(nrow(AMARETTOinit$RegulatorData)>1,TRUE)
  expect_equal(nrow(AMARETTOinit$RegulatorData)>1,TRUE)
  expect_equal(sum(AMARETTOinit$RegulatorData!=0)>0,TRUE)
  expect_equal(length(AMARETTOinit$RegulatorAlterations),4)
  expect_equal(ncol(AMARETTOinit$RegulatorAlterations$MET),3)
  expect_equal(ncol(AMARETTOinit$RegulatorAlterations$CNV),3)
  expect_equal(sum(rowSums(AMARETTOinit$RegulatorAlterations$Summary)>=1) ,nrow(AMARETTOinit$RegulatorAlterations$Summary))
  expect_equal(min(AMARETTOinit$ModuleMembership),1)
  expect_equal(max(AMARETTOinit$ModuleMembership),NrModules)
  expect_equal(length(AMARETTOinit$ModuleMembership)>1,TRUE)
  expect_equal(AMARETTOinit$NrCores!=0,TRUE)
  expect_equal(sum(is.na(AMARETTOinit$MA_matrix_Var)),0)
  expect_equal(sum(is.na(AMARETTOinit$RegulatorData)),0)
  expect_equal(sum(is.na(AMARETTOinit$RegulatorAlterations$MET)),0)
  expect_equal(sum(is.na(AMARETTOinit$RegulatorAlterations$CNV)),0)
})

testthat::test_that("Checking AMARETTOresults object to be in decent shape",{
  expect_equal(length(AMARETTOresults),9)
  expect_equal(AMARETTOresults$NrModules,NrModules)
  expect_equal(nrow(AMARETTOresults$RegulatoryPrograms),NrModules)
  expect_equal(ncol(AMARETTOresults$RegulatoryPrograms)>1,TRUE)
  expect_equal(length(AMARETTOresults$AllRegulators)>1,TRUE)
  expect_equal(ncol(AMARETTOresults$RegulatoryPrograms),length(AMARETTOresults$AllRegulators))
  expect_equal(length(AMARETTOresults$AllGenes)>1,TRUE)
  expect_equal(nrow(AMARETTOresults$ModuleMembership),length(AMARETTOresults$AllGenes))
  expect_equal(ncol(AMARETTOresults$ModuleMembership),1)
  expect_equal(max(AMARETTOresults$ModuleMembership),NrModules)
  expect_equal(min(AMARETTOresults$ModuleMembership),1)
  expect_equal(nrow(AMARETTOresults$AutoRegulationReport),NrModules)
  expect_equal(nrow(AMARETTOresults$ModuleData),NrModules)
  expect_equal(sum(AMARETTOresults$ModuleData!=0)>1,TRUE)
  expect_equal(nrow(AMARETTOresults$RegulatoryProgramData),NrModules)
  expect_equal(sum(AMARETTOresults$RegulatoryProgramData!=0)>1,TRUE)
  expect_equal(ncol(AMARETTOresults$ModuleData),ncol(AMARETTOresults$RegulatoryProgramData))
  expect_equal(sum(is.na(AMARETTOresults$ModuleMembership)),0)
  expect_equal(sum(is.na( AMARETTOresults$RegulatoryPrograms)),0)
  expect_equal(sum(is.na(AMARETTOresults$ModuleData)),0)
  expect_equal(sum(is.na(AMARETTOresults$RegulatoryProgramData)),0)
})



