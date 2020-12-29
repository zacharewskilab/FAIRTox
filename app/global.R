#' global.R
#' Author: Jack Dodson
#' Michigan State University - Zacharewski Lab

library(shiny) #App framework
library(RSQLite) #SQL database import
library(ggplot2) #Figure generation
library(plotly) #Interactive figures
library(shinyjs) #Tableau implementation
library(shinycssloaders) #Loading animations
library(Seurat) #For single cell data
library(ggm) #Enables powersets
library(reshape2) #Array/vector operations
library(plyr) #Array/vector operations
library(tibble) #Array operations
library(dplyr) #Advanced array operations
library(UpSetR) #UpSet plots
library(ggnewscale) #tzheatmap requirement
library(stringr) #Advanced string ops
library(fgsea) #GSEA analysis
library(data.table) #Advanced dataframe ops
library(readxl) #Read .xl files
library(tidyverse) #Data organization and manipulation
library(factoextra) #Dimensionality reduction for PCA
library(matrixStats) #Required for PCA
library(FactoMineR) #Required for PCA
library(randomcoloR) #Large randomized color palettes
library(mixOmics) #Required for PLS
library(reticulate) #Required for single cell plots
library(ggridges) #Single cell ridge plots

# Load Single Cell Metadata
sc_dataset_meta <- read_excel("C:\\Users\\Jack\\Desktop\\FAIRTox_github\\app\\RData\\SingleCell_Metadata.xlsx")

# Enable server-side bookmarking in /shinybookmarks folder
enableBookmarking(store = "server")

# Load dbZach
dbPath <- "./RData/full-dbzach3.db"
dbconn <- dbConnect(dbDriver("SQLite"), dbname = dbPath)

# Load single cell data
loadSingleCellData = function(directory){
  # directory = the directory path with the single-cell data
  sn.files = list.files(directory, full.names = TRUE)
  umap = read.table(sn.files[grepl('umap', sn.files)], header = TRUE)[-1, ]
  meta = read.table(sn.files[grepl('metadata', sn.files)], header = TRUE, sep = '\t')[-1, ]
  genes.order = read.table(sn.files[grepl('genes.tsv', sn.files)], header = FALSE)$V1
  barcodes.order = read.table(sn.files[grepl('barcodes.tsv', sn.files)], header = FALSE)$V1
  mtxConn = RSQLite::dbConnect(RSQLite::dbDriver("SQLite"), dbname = sn.files[grepl('.sql', sn.files)])
  return(list(umap = umap, meta = meta, gene.order = genes.order, barcodes.order = barcodes.order, conn = mtxConn))
}

# Load gene data
getGeneData = function(gene.list, barcode.list, scData){
  # gene.list = the genes of interest (e.g. c('Cyp1a1', 'Cyp1a2', ...))
  # barcode.list = the barcodes of interest. For all barcodes use {loaded-data}$barcodes.order
  # scData = the list object from the loadSingleCellData function
  genes.ind = match(gene.list, scData$gene.order)
  barcodes.ind = match(barcode.list, scData$barcodes.order)
  query = paste(
    "SELECT * FROM log_norm_umi WHERE gene_row IN (",
    paste(genes.ind, collapse = ", "),
    ") AND barcode_row IN (",
    paste(barcodes.ind, collapse = ", "),
    ");",
    sep = ''
  )
  gene.df = RSQLite::dbGetQuery(scData$conn, query)
  gene.df$gene_row = scData$gene.order[gene.df$gene_row]
  gene.df$barcode_row = scData$barcodes.order[gene.df$barcode_row]
  return(gene.df)
}

# Load Broad Format data
sn_Small <- loadSingleCellData("./RData/BroadFormat/Small")
sn_Large <- loadSingleCellData("./RData/BroadFormat/Large")

# Builds empty dataframe for error catching
empty.frame <- data.frame(matrix(ncol = 1, nrow = 1))
colnames(empty.frame) <- c("No valid studies matching criteria available")

# Creates lists of unique variable names for each factors of interest
species <- (dbGetQuery(dbconn, "SELECT DISTINCT Species_common FROM
                                (SELECT DISTINCT Species_ID FROM ExpressionChange) as A 
                                JOIN Species as S on S.Species_ID = A.Species_ID"))[, 1]
sex <- (dbGetQuery(dbconn, "SELECT DISTINCT Sex FROM Sex"))[, 1]
design <- (dbGetQuery(dbconn, "SELECT DISTINCT DesignType FROM StudyDesigns"))[, 1]
organs <- (dbGetQuery(dbconn, "SELECT DISTINCT Organ_name, Project_ID FROM
                               (SELECT DISTINCT Organ_name_ID, Project_ID FROM ExpressionChange) as A 
                               JOIN Organ_name as O on O.Organ_name_ID = A.Organ_name_ID"))
genes <- (dbGetQuery(dbconn, "SELECT DISTINCT Symbol FROM Annotation"))[, 1]
projects <- (dbGetQuery(dbconn, "SELECT DISTINCT A.Project_ID, Longname, DesignType FROM 
                                (SELECT * FROM ExpressionChange) as A
                                JOIN Project as P ON P.Project_ID = A.Project_ID
                                JOIN StudyDesigns as S ON S.StudyDesign_ID = A.StudyDesign_ID
                                ORDER BY A.Project_ID ASC"))
doses <- (dbGetQuery(dbconn, "SELECT DISTINCT Dose, Project_ID, ZT, DesignType FROM
                               (SELECT DISTINCT Dose, Project_ID, ZT, StudyDesign_ID FROM ExpressionChange) as A 
                              JOIN StudyDesigns as S on S.StudyDesign_ID = A.StudyDesign_ID"))
doses$ZT <- as.double(doses$ZT)
foldChanges <-(dbGetQuery(dbconn, "SELECT DISTINCT FoldChange FROM ExpressionChange ORDER BY FoldChange DESC"))[,1]
longnames <- (dbGetQuery(dbconn, "SELECT DISTINCT Project_ID, Longname FROM Project ORDER BY Project_ID ASC"))[,1:2]
treatments <- (dbGetQuery(dbconn, "SELECT DISTINCT Chemical_Name FROM
                                   (SELECT DISTINCT Chemical_ID FROM ExpressionChange) as A
                                   JOIN Chemical on Chemical.Chemical_ID = A.Chemical_ID"))[,1]
assays <- (dbGetQuery(dbconn, "SELECT DISTINCT Assay_Name FROM Assays"))[,1]
strains <- (dbGetQuery(dbconn, "SELECT DISTINCT Strain_Name FROM Strain"))[,1]
colors <- c("mediumblue", "red3", "darkorange", "seagreen")
linetypes <- c("solid", "dotted", "twodash", "longdash")
metadata <- c("Species", "Sex", "Organ", "Treatment", "Assay", "Strain")

# Read homology ID's
homology_IDs <- (dbGetQuery(dbconn, "SELECT * FROM temp_homology"))

# JavaScript code for app-wide resets
jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

# Read pca and pls files
pca_metadata <- read.delim("countfiles_metadata.txt", header = TRUE, sep = '\t')
pls_metadata <- read.delim("countfiles_metadata.txt", header = TRUE, sep = '\t')

# Create large distinct color palette (used in PCA)
mypal <- distinctColorPalette(100)

# Load single-cell RNA-sequencing data
#load('./10X-Zacharewski-LiverData.RData')
