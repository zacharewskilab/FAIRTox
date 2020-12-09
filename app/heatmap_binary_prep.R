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
heatmap_binary_prep <- function(in.data, title = runif(1, min = 1000, max = 2000), binColumns)
{
  #Check for count columns
  if (colnames(in.data)[1] != 'Gene'){
    stop("ERROR: Gene column missing or incorrectly named")
  }

  dataTable = in.data[,c(1,binColumns)]
  dataTable = reshape2::melt(dataTable, id.vars = 'Gene')
  dataTable$type = 'binary'
  dataTable$format = 'none'
  dataTable$title = title

  return(dataTable)

}
