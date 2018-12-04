context("test-my-test")
set.seed(3.14)
test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("AMARETTO_Initialize", {
  data(ProcessedDataLIHC)
  NrModules <- 50
  VarPercentage <- 50
  expect_that(
    AMARETTOinit <-
      AMARETTO_Initialize(
        MA_matrix = ProcessedDataLIHC$MA_matrix,
        CNV_matrix = ProcessedDataLIHC$CNV_matrix,
        MET_matrix = ProcessedDataLIHC$MET_matrix,
        NrModules = NrModules,
        VarPercentage = VarPercentage
      ),
    is_a('list')
  )
  
})


test_that("AMARETTO_Run", {
  expect_that(AMARETTOresults <- AMARETTO_Run(AMARETTOinit),
              is_a('list'))
  
})


test_that("AMARETTO_EvaluateTestSet", {
  expect_that(
    AMARETTOtestReport <-
      AMARETTO_EvaluateTestSet(
        AMARETTOresults = AMARETTOresults,
        MA_Data_TestSet = AMARETTOinit$MA_matrix_Var,
        RegulatorData_TestSet = AMARETTOinit$RegulatorData
      )
    ,
    is_a('matrix')
  )
  
})