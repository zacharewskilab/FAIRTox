source('../build_dataframes.R', chdir = TRUE)
source('../global.R', chdir = TRUE)
library(testthat)

#' Test buildBatchAnalysisdf function in build_dataframes.R
#' @example buildBatchAnalysisdf(connection)

# Test structure
test_that("test structure", {
  df <- buildBatchAnalysisdf(dbconn)
  expected_colnames <- c("Ensembl_ID", "H_ID", "Project", "AID", "Species", "Strain", "Sex", "Tissue",    
                         "Chemical", "Dose", "Timepoint", "ZT", "Count" )
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Ensembl_ID, "character")
  expect_type(df$H_ID, "integer")
  expect_type(df$Project, "integer")
  expect_type(df$AID, "character")
  expect_type(df$Species, "character")
  expect_type(df$Strain, "character")
  expect_type(df$Sex, "character")
  expect_type(df$Tissue, "character")
  expect_type(df$Chemical, "character")
  expect_type(df$Dose, "double")
  expect_type(df$Timepoint, "integer")
  expect_type(df$ZT, "character")
  expect_type(df$Count, "integer")
})
