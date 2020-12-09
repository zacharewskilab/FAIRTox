source('../build_dataframes.R', chdir = TRUE)
source('../global.R', chdir = TRUE)
library(testthat)

#' Test builddfTimeCourse function in build_dataframes.R
#' @example builddfTimeCourse(annotation, symbol, connection)

# Test with single gene - Symbol
test_that("single gene - symbol", {
  df <- builddfTimeCourse("Symbol", "Cyp1a1", dbconn)
  expected_colnames <- c("Project_ID", "Longname", "TimePoint", "FoldChange", "Species_common",
                         "Sex", "Organ_name", "P1t", "PMID", "GEO", "DesignType", "Chemical_Name",
                         "Assay_Name", "Strain_Name")
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Project_ID, "character")
  expect_type(df$Longname, "character")
  expect_type(df$TimePoint, "double")
  expect_type(df$FoldChange, "double")
  expect_type(df$Species_common, "character")
  expect_type(df$Sex, "character")
  expect_type(df$Organ_name, "character")
  expect_type(df$P1t, "double")
  expect_type(df$PMID, "integer")
  expect_type(df$GEO, "character")
  expect_type(df$DesignType, "character")
  expect_type(df$Chemical_Name, "character")
  expect_type(df$Assay_Name, "character")
  expect_type(df$Strain_Name, "character")
})

# Test with multiple genes - Symbol
test_that("multiple genes - symbol", {
  expect_error(builddfTimeCourse("Symbol", c("Cyp1a1", "Cyp1a2"), dbconn))
})

# Test with no genes - Symbol
test_that("no genes - symbol", {
  df <- builddfTimeCourse("Symbol", "", dbconn)
  expect_equal(nrow(df), 0)
})

# Test with single gene - Ensembl_ID
test_that("single gene - Ensembl_ID", {
  df <- builddfTimeCourse("Ensembl_ID", "ENSMUSG00000032315", dbconn)
  expected_colnames <- c("Project_ID", "Longname", "TimePoint", "FoldChange", "Species_common",
                         "Sex", "Organ_name", "P1t", "PMID", "GEO", "DesignType", "Chemical_Name",
                         "Assay_Name", "Strain_Name")
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Project_ID, "character")
  expect_type(df$Longname, "character")
  expect_type(df$TimePoint, "double")
  expect_type(df$FoldChange, "double")
  expect_type(df$Species_common, "character")
  expect_type(df$Sex, "character")
  expect_type(df$Organ_name, "character")
  expect_type(df$P1t, "double")
  expect_type(df$PMID, "integer")
  expect_type(df$GEO, "character")
  expect_type(df$DesignType, "character")
  expect_type(df$Chemical_Name, "character")
  expect_type(df$Assay_Name, "character")
  expect_type(df$Strain_Name, "character")
})

# Test with multiple genes - Ensembl_ID
test_that("multiple genes - Ensembl_ID", {
  expect_error(builddfTimeCourse("Ensembl_ID", c("ENSMUSG00000032315", "ENSMUSG00000032310"), dbconn))
})

# Test with no genes - Ensembl_ID
test_that("no genes - Ensembl_ID", {
  df <- builddfTimeCourse("Ensembl_ID", "", dbconn)
  expect_equal(nrow(df), 0)
})

# Test with single gene - NCBI_ID
test_that("single gene - NCBI_ID", {
  df <- builddfTimeCourse("NCBI_ID", "13076", dbconn)
  expected_colnames <- c("Project_ID", "Longname", "TimePoint", "FoldChange", "Species_common",
                         "Sex", "Organ_name", "P1t", "PMID", "GEO", "DesignType", "Chemical_Name",
                         "Assay_Name", "Strain_Name")
  expect_equal(colnames(df), expected_colnames)
  expect_s3_class(df, "data.frame")
  expect_type(df, "list")
  expect_type(df$Project_ID, "character")
  expect_type(df$Longname, "character")
  expect_type(df$TimePoint, "double")
  expect_type(df$FoldChange, "double")
  expect_type(df$Species_common, "character")
  expect_type(df$Sex, "character")
  expect_type(df$Organ_name, "character")
  expect_type(df$P1t, "double")
  expect_type(df$PMID, "integer")
  expect_type(df$GEO, "character")
  expect_type(df$DesignType, "character")
  expect_type(df$Chemical_Name, "character")
  expect_type(df$Assay_Name, "character")
  expect_type(df$Strain_Name, "character")
})

# Test with multiple genes - NCBI_ID
test_that("multiple genes - NCBI_ID", {
  expect_error(builddfTimeCourse("NCBI_ID", c("13076", "13077"), dbconn),
               "More than one gene selected.")
})

# Test with no genes - NCBI_ID
test_that("no genes - NCBI_ID", {
  df <- builddfTimeCourse("NCBI_ID", "", dbconn)
  expect_equal(nrow(df), 0)
})
