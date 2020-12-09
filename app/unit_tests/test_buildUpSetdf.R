source('../build_dataframes.R', chdir = TRUE)
source('../global.R', chdir = TRUE)
library(testthat)

#' Test buildUpSetdf function in build_dataframes.R
#' @example buildUpSetdf(selected, connection)

# Test less than two projects selected
test_that("no projects", {
  # Test no projects
  expect_error(buildUpSetdf(c(), dbconn), "Select at least two projects")
  # Test one project
  expect_error(buildUpSetdf(c("Mouse female 92 days repeated TCDD dose-response"), dbconn), 
               "Select at least two projects")
})

# Test structure of dataframe
test_that("test structure", {
  df <- buildUpSetdf(c("Mouse female 92 days repeated TCDD dose-response",
                       "Mouse female 28 days repeated TCDD dose-response \n"),
                     dbconn)
  expected_colnames <- c("Longname", "P1t", "Dose", "ZT", "FoldChange", "Symbol", "Organ_name")
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Longname, "character")
  expect_type(df$P1t, "double")
  expect_type(df$Dose, "integer")
  expect_type(df$ZT, "integer")
  expect_type(df$FoldChange, "double")
  expect_type(df$Symbol, "character")
  expect_type(df$Organ_name, "character")
})
