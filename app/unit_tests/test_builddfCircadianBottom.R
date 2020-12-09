source('../build_dataframes.R', chdir = TRUE)
source('../global.R', chdir = TRUE)
library(testthat)

#' Test builddfCircadianBottom function in build_dataframes.R
#' @example builddfCircadianBottom(annotation, symbol, connection)

# Test with single gene - Symbol
test_that("single gene - symbol", {
  df <- builddfCircadianBottom("Symbol", "Cyp1a1", dbconn)
  expected_colnames <- c("ZT", "Project_ID", "Dose", "NormalizedSignal")
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$ZT, "integer")
  expect_type(df$Project_ID, "integer")
  expect_type(df$Dose, "integer")
  expect_type(df$NormalizedSignal, "double")
})

# Test with multiple genes - Symbol
test_that("multiple genes - symbol", {
  expect_error(builddfCircadianBottom("Symbol", c("Cyp1a1", "Cyp1a2"), dbconn))
})

# Test with no genes - Symbol
test_that("no genes - symbol", {
  df <- builddfCircadianBottom("Symbol", "", dbconn)
  expect_equal(nrow(df), 0)
})

# Test with single gene - Ensembl_ID
test_that("single gene - Ensembl_ID", {
  df <- builddfCircadianBottom("Symbol", "ENSMUSG00000032315", dbconn)
  expected_colnames <- c("ZT", "Project_ID", "Dose", "NormalizedSignal")
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$ZT, "integer")
  expect_type(df$Project_ID, "integer")
  expect_type(df$Dose, "integer")
  expect_type(df$NormalizedSignal, "double")
})

# Test with multiple genes - Ensembl_ID
test_that("multiple genes - Ensembl_ID", {
  expect_error(builddfCircadianBottom("Ensembl_ID", c("ENSMUSG00000032315", "ENSMUSG00000032310"), dbconn))
})

# Test with no genes - Ensembl_ID
test_that("no genes - Ensembl_ID", {
  df <- builddfCircadianBottom("Ensembl_ID", "", dbconn)
  expect_equal(nrow(df), 0)
})

# Test with single gene - NCBI_ID
test_that("single gene - NCBI_ID", {
  df <- builddfCircadianBottom("Symbol", "13076", dbconn)
  expected_colnames <- c("ZT", "Project_ID", "Dose", "NormalizedSignal")
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$ZT, "integer")
  expect_type(df$Project_ID, "integer")
  expect_type(df$Dose, "integer")
  expect_type(df$NormalizedSignal, "double")
})

# Test with multiple genes - NCBI_ID
test_that("multiple genes - NCBI_ID", {
  expect_error(builddfCircadianBottom("NCBI_ID", c("13076", "13077"), dbconn),
               "More than one gene selected.")
})

# Test with no genes - NCBI_ID
test_that("no genes - NCBI_ID", {
  df <- builddfCircadianBottom("NCBI_ID", "", dbconn)
  expect_equal(nrow(df), 0)
})
