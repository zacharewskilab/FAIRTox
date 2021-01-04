library(dplyr)
library(tidyr)
library(reshape2)

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
sn_Small <- loadSingleCellData("C:\\Users\\Jack\\Desktop\\FAIRTox_github\\app\\RData\\BroadFormat\\Small")
sn_Large <- loadSingleCellData("C:\\Users\\Jack\\Desktop\\FAIRTox_github\\app\\RData\\BroadFormat\\Large")

sn <- sn_Small

gene_df <- getGeneData(c("Ces3a", "Car3", "Mup7", "C8a", "Bhmt", "Elovl3"), sn$barcodes.order, sn)
colnames(gene_df) <- c("GENE", "NAME", "VALUE")
# Merge on barcode
gene_df <- merge(gene_df, sn$meta, by = "NAME", all.y = TRUE)

data <- spread(data, GENE, VALUE)

data <- data[order(VALUE),]

data <- data[1:50,]
row.names(data) <- data$NAME
data <- data[,3:ncol(data)]
data <- t(data)

saveRDS(data, 'C:\\Users\\Jack\\Desktop\\dittoHeatmapTestData.rds')