source('../build_dataframes.R', chdir = TRUE)
source('../global.R', chdir = TRUE)
library(testthat)

#' Test buildGSEAdf function in build_dataframes.R
#' @example buildGSEAdf(selected, connection)

# Test no genes selected
test_that("no genes selected", {
  expect_error(buildGSEAdf(c(), dbconn), "Select at least one gene")
})

# Test structure
test_that("test structure", {
  df <- buildGSEAdf(c("Cyp1a1", "Cyp1a2"), dbconn)
  expected_colnames <- c("Symbol", "P1t", "Dose", "FoldChange", "Longname", "Species_common", "Sex",
                         "Organ_name", "Chemical_Name", "Assay_Name", "Strain_Name", "Project_ID")
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Symbol, "character")
  expect_type(df$P1t, "double")
  expect_type(df$Dose, "double")
  expect_type(df$FoldChange, "double")
  expect_type(df$Longname, "character")
  expect_type(df$Species_common, "character")
  expect_type(df$Sex, "character")
  expect_type(df$Organ_name, "character")
  expect_type(df$Chemical_Name, "character")
  expect_type(df$Assay_Name, "character")
  expect_type(df$Strain_Name, "character")
  expect_type(df$Project_ID, "integer")
})
