source('../build_dataframes.R', chdir = TRUE)
source('../global.R', chdir = TRUE)
library(testthat)

#' Test builddfCircadianTop function in build_dataframes.R
#' @example builddfCircadianTop(annotation, symbol, connection)

# Test with single gene - Symbol
test_that("single gene - symbol", {
  df <- builddfCircadianTop("Symbol", "Cyp1a1", dbconn)
  expected_colnames <- c("Project_ID", "ZT", "FoldChange", "Longname", "P1t", "Organ_name", "DesignType")
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Project_ID, "character")
  expect_type(df$ZT, "double")
  expect_type(df$FoldChange, "double")
  expect_type(df$Longname, "character")
  expect_type(df$P1t, "double")
  expect_type(df$Organ_name, "character")
  expect_type(df$DesignType, "character")
})

# Test with multiple genes - Symbol
test_that("multiple genes - symbol", {
  expect_error(builddfCircadianTop("Symbol", c("Cyp1a1", "Cyp1a2"), dbconn))
})

# Test with no genes - Symbol
test_that("no genes - symbol", {
  df <- builddfCircadianTop("Symbol", "", dbconn)
  expect_equal(nrow(df), 0)
})

# Test with single gene - Ensembl_ID
test_that("single gene - Ensembl_ID", {
  df <- builddfCircadianTop("Symbol", "ENSMUSG00000032315", dbconn)
  expected_colnames <- c("Project_ID", "ZT", "FoldChange", "Longname", "P1t", "Organ_name", "DesignType")
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Project_ID, "character")
  expect_type(df$ZT, "double")
  expect_type(df$FoldChange, "double")
  expect_type(df$Longname, "character")
  expect_type(df$P1t, "double")
  expect_type(df$Organ_name, "character")
  expect_type(df$DesignType, "character")
})

# Test with multiple genes - Ensembl_ID
test_that("multiple genes - Ensembl_ID", {
  expect_error(builddfCircadianTop("Ensembl_ID", c("ENSMUSG00000032315", "ENSMUSG00000032310"), dbconn))
})

# Test with no genes - Ensembl_ID
test_that("no genes - Ensembl_ID", {
  df <- builddfCircadianTop("Ensembl_ID", "", dbconn)
  expect_equal(nrow(df), 0)
})

# Test with single gene - NCBI_ID
test_that("single gene - NCBI_ID", {
  df <- builddfCircadianTop("Symbol", "13076", dbconn)
  expected_colnames <- c("Project_ID", "ZT", "FoldChange", "Longname", "P1t", "Organ_name", "DesignType")
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Project_ID, "character")
  expect_type(df$ZT, "double")
  expect_type(df$FoldChange, "double")
  expect_type(df$Longname, "character")
  expect_type(df$P1t, "double")
  expect_type(df$Organ_name, "character")
  expect_type(df$DesignType, "character")
})

# Test with multiple genes - NCBI_ID
test_that("multiple genes - NCBI_ID", {
  expect_error(builddfCircadianTop("NCBI_ID", c("13076", "13077"), dbconn),
               "More than one gene selected.")
})

# Test with no genes - NCBI_ID
test_that("no genes - NCBI_ID", {
  df <- builddfCircadianTop("NCBI_ID", "", dbconn)
  expect_equal(nrow(df), 0)
})
