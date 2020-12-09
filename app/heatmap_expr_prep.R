#' Prepare a heatmap from bulk RNAseq data
#'
#' Draws a heatmap in the lab style from time-course or dose-response
#' bulk RNAsequencing data important as a tab-delimited text file.
#'
#' @author Rance Nault, \email{naultran@@msu.edu}
#'
#' @param in.path complete path to tab-delimited text file. The file
#' MUST have a column named "Gene". Putative Dioxin Response Element
#' pDRE data column must be titled "DRE" but is \code{not required}.
#' AhR data must be titled "AhR" but is \code{not required}. Count
#' data must be titled "Counts" but is \code{not required}. Rythmic
#' notation data columns must be titled "Veh_Cir" and "TCDD_Cir" but
#' is \code{not required}. Expression data should be a g x n matrix
#' where g are genes and n are treatment groups, represented as fold-
#' changes \code{not in log2 space}. To include asterisks, P1(t)
#' columns should follow the same format as expression data and
#' columns must match exactly expression data beginning with "P1_".
#' The asterisks flag should be set to TRUE.
#' @param expression.order a vector of treatment groups in the order
#' in which they should be plotted.
#' @param exp.lim a vector of minimum and maximum log2 fold-change
#' thresholds
#' @param counts.lim a vector of minimum and maximum count tresholds
#' @param cat.order a vector for the order to plot individual blocks
#' @param decimal set to the character used as decimal for numbers
#' @param asterisks indicates whether to plot asterisks significance
#'
#' @return a ggplot object
#'
#' @examples
#' plot.heatmap('time-course-file.txt', c(2,4,8,12,24,72,168))
#' plot.heatmap('time-course-file.txt', c(2,4,8,12,24,72,168), asterisks = TRUE)
#'
#' @export
heatmap_expr_prep <- function(in.data, title = runif(1, min = 1000, max = 2000), exprColumns, pColumns = c(), countColumn = c(),
                         logFC = TRUE, exp.lim = c(-2,2), count.lim = c(0,3000),
                         fcPal, countPal)
{
  #Check for count columns
  if (colnames(in.data)[1] != 'Gene'){
    stop("ERROR: Gene column missing or incorrectly named")
  }

  #Check for pvalue columns
  if (length(pColumns) == 0){
    sig = FALSE
    print("WARNING: Data does not contain statistical significance")
  } else {
    sig = TRUE
  }

  #Check for count columns
  if (length(countColumn) != 1){
    count = FALSE
    print("WARNING: Data does not contain count data")
  } else {
    count = TRUE
  }

  #Check for count columns
  if (length(exprColumns) != length(pColumns)){
    stop("ERROR: Gene column missing or incorrectly named")
  }

  dataTable = in.data[,c(1,exprColumns)]
  dataTable = reshape2::melt(dataTable, id.vars = 'Gene')
  dataTable$type = 'expression'
  if (logFC) {
    dataTable$value = log(dataTable$value,2)
  }
  dataTable$format = map2color(dataTable$value, fcPal, limits = c(-2,2))
  dataTable$variable = factor(dataTable$variable, levels = colnames(in.data)[exprColumns])


  if (sig) {
    pTable = in.data[,c(1,pColumns)]
    pTable = reshape2::melt(pTable, id.vars = 'Gene')
    pTable$type = 'significance'
    pTable$variable = dataTable$variable
    pTable$format= ''
    pTable[which(dataTable$value < log(1.5,2) & dataTable$value > log((1/1.5),2)), "value"] = 0
    dataTable = rbind(dataTable, pTable)
  }

  if (count) {
    cTable = in.data[,c(1, countColumn)]
    cTable = reshape2::melt(cTable, id.vars = 'Gene')
    cTable$type = 'count'
    cTable$format = '#E5E5FF'
    dataTable = rbind(dataTable, cTable)
  }

  dataTable$title = title
  return(dataTable)
}



















