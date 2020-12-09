source('../build_dataframes.R', chdir = TRUE)
source('../global.R', chdir = TRUE)
library(testthat)

#' Test builddfBE function in build_dataframes.R
#' @example builddfBE(annotation, symbol, connection)

# Test with single gene - Symbol
test_that("single gene - symbol", {
  df <- builddfBE("Symbol", "Cyp1a1", dbconn)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Project_ID, "integer")
  expect_type(df$Species_ID, "integer")
  expect_type(df$Sex_ID, "integer")
  expect_type(df$Strain_Name, "character")
  expect_type(df$NormalizedSignal, "double")
  expect_type(df$Species_common, "character")
  expect_type(df$Sex, "character")
})

# Test with multiple genes - Symbol
test_that("multiple genes - symbol", {
  expect_error(builddfBE("Symbol", c("Cyp1a1", "Cyp1a2"), dbconn))
})

# Test with no genes - Symbol
test_that("no genes - symbol", {
  df <- builddfBE("Symbol", "", dbconn)
  expect_equal(nrow(df), 0)
})

# Test with single gene - Ensembl_ID
test_that("single gene - Ensembl_ID", {
  df <- builddfBE("Symbol", "ENSMUSG00000032315", dbconn)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Project_ID, "integer")
  expect_type(df$Species_ID, "integer")
  expect_type(df$Sex_ID, "integer")
  expect_type(df$Strain_Name, "character")
  expect_type(df$NormalizedSignal, "double")
  expect_type(df$Species_common, "character")
  expect_type(df$Sex, "character")
})

# Test with multiple genes - Ensembl_ID
test_that("multiple genes - Ensembl_ID", {
  expect_error(builddfBE("Ensembl_ID", c("ENSMUSG00000032315", "ENSMUSG00000032310"), dbconn))
})

# Test with no genes - Ensembl_ID
test_that("no genes - Ensembl_ID", {
  df <- builddfBE("Ensembl_ID", "", dbconn)
  expect_equal(nrow(df), 0)
})

# Test with single gene - NCBI_ID
test_that("single gene - NCBI_ID", {
  df <- builddfBE("Symbol", "13076", dbconn)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Project_ID, "integer")
  expect_type(df$Species_ID, "integer")
  expect_type(df$Sex_ID, "integer")
  expect_type(df$Strain_Name, "character")
  expect_type(df$NormalizedSignal, "double")
  expect_type(df$Species_common, "character")
  expect_type(df$Sex, "character")
})

# Test with multiple genes - NCBI_ID
test_that("multiple genes - NCBI_ID", {
  expect_error(builddfBE("NCBI_ID", c("13076", "13077"), dbconn),
               "More than one gene selected.")
})

# Test with no genes - NCBI_ID
test_that("no genes - NCBI_ID", {
  df <- builddfBE("NCBI_ID", "", dbconn)
  expect_equal(nrow(df), 0)
})
