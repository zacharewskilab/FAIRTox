dittoHeatmap <- function(
  broad_input,
  object = NULL,
  genes = getGenes(object, assay),
  metas = NULL,
  cells.use = NULL,
  annot.by = NULL,
  order.by = .default_order(object, annot.by),
  main = NA,
  cell.names.meta = NULL,
  assay = .default_assay(object),
  slot = .default_slot(object),
  swap.rownames = NULL,
  heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50),
  scaled.to.max = FALSE,
  heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
  annot.colors = c(dittoColors(),dittoColors(1)[seq_len(7)]),
  annotation_col = NULL,
  annotation_colors = NULL,
  data.out=FALSE,
  highlight.features = NULL,
  highlight.genes = NULL, 
  show_colnames = isBulk(object),
  show_rownames = TRUE,
  scale = "row",
  cluster_cols = isBulk(object),
  border_color = NA,
  legend_breaks = NA,
  breaks = NA,
  complex = FALSE,
  ...) {
  
  # If cells.use given as logical, populate as names.
  # cells.use <- .which_cells(cells.use, object)
  # all.cells <- .all_cells(object)
  
  # Handle deprecated argument
  if (!is.null(highlight.genes)) {
    if (!is.null(highlight.features)) {
      stop("you can only specify one of 'highlight.genes' or 'highlight.features'")
    }
    .Deprecated(msg="argument 'highlight.genes' is deprecated, please use 'highlight.features' instead")
    highlight.features <- highlight.genes
  }
  
  ### Obtain all needed data
  object <- .swap_rownames(object, swap.rownames)
  
  # Replace with whatever test data in correct format
  #         sample1   sample2   sample3   ...
  # gene1      .         .         .      ...
  # gene2      .         .         .      ...
  #
  #
  data <- broad_input
  
  print(head(data))
  
  if (!is.null(cell.names.meta)) {
    cell.names <- .var_OR_get_meta_or_gene(cell.names.meta, object)
  } else {
    cell.names <- NULL
  }
  
  if (!is.null(order.by)) {
    order_data <- .var_OR_get_meta_or_gene(order.by, object, assay, slot)
  } else {
    order_data <- NULL
  }
  
  # Make the columns annotations data
  if (is.null(annotation_col)) {
    annotation_col <- data.frame(row.names = cells.use)
  } else if (!all(cells.use %in% rownames(annotation_col))) {
    stop("rows of 'annotation_col' must be cell/sample names of all cells/samples being displayed")
  }
  if (!is.null(annot.by)) {
    annotation_col <- rbind(
      as.data.frame(getMetas(object, names.only = FALSE)[cells.use, annot.by, drop = FALSE]),
      annotation_col[cells.use, , drop=FALSE])
  }
  
  # Set title (off = NA in pheatmap)
  if (is.null(main)) {
    main <- NA
  }
  
  ### Prep inputs for heatmap contructor calls
  args <- .prep_ditto_heatmap(
    data, cells.use, all.cells, cell.names, order_data, main,
    heatmap.colors, scaled.to.max, heatmap.colors.max.scaled, annot.colors,
    annotation_col, annotation_colors, highlight.features, TRUE,
    show_rownames, scale, TRUE, border_color, legend_breaks,
    breaks, ...)
  
  if (data.out) {
    OUT <- args
  } else if (complex) {
    .error_if_no_complexHm()
    OUT <- do.call(ComplexHeatmap::pheatmap, args)
  } else {
    OUT <- do.call(pheatmap::pheatmap, args)
  }
  OUT
}

.prep_ditto_heatmap <- function(
  data,
  cells.use,
  all.cells,
  cell.names,
  order_data,
  main,
  heatmap.colors,
  scaled.to.max,
  heatmap.colors.max.scaled,
  annot.colors,
  annotation_col,
  annotation_colors,
  highlight.features,
  show_colnames,
  show_rownames,
  scale,
  cluster_cols,
  border_color,
  legend_breaks,
  breaks,
  ...) {
  
  # Create the base pheatmap inputs
  args <- list(
    mat = data, main = main, show_colnames = show_colnames,
    show_rownames = show_rownames, color = heatmap.colors,
    cluster_cols = cluster_cols, border_color = border_color,
    scale = scale, breaks = breaks, legend_breaks = legend_breaks, ...)
  
  print(head(args))
  
  # Adjust data
  if (!is.null(order_data)){
    args$mat <- args$mat[,order(order_data[all.cells %in% cells.use])]
    
    if (!identical(annotation_col, NULL)) {
      annotation_col <- annotation_col[colnames(args$mat),, drop = FALSE]
    }
  }
  if (scaled.to.max) {
    args <- .scale_to_max(args, heatmap.colors.max.scaled)
  }
  
  # Add annotation_col / annotation_colors only if needed
  if (ncol(annotation_col)>0 || !is.null(args$annotation_row)) {
    if (ncol(annotation_col)>0) {
      args$annotation_col <- annotation_col
    }
    args$annotation_colors <- annotation_colors
    # Add any missing annotation colors
    args <- .make_heatmap_annotation_colors(args, annot.colors)
  }
  
  # Make a labels_row input for displaying only certain genes/variables if given to 'highlight.features'
  if (!(is.null(highlight.features)) && sum(highlight.features %in% rownames(data))>0) {
    highlight.features <- highlight.features[highlight.features %in% rownames(data)]
    args$labels_row <- rownames(data)
    #Overwrite all non-highlight genes rownames to ""
    args$labels_row[-(match(highlight.features,rownames(data)))] <- ""
    args$show_rownames <- TRUE
  }
  
  # Add cell/sample/row names unless provided separately by user
  if(is.null(args$labels_col) && !is.null(cell.names)) {
    args$labels_col <- as.character(cell.names[colnames(args$mat)])
    args$show_colnames <- TRUE
  }
  
  args
}

# .scale_to_max <- function(args, heatmap.colors.max.scaled) {
#   # Max scales the data matrix,
#   # addss max-sclaing-colors
#   # adjusts the color and legend breaks accordingly (unless set by user)
#   maxs <- apply(args$mat,1,max)
#   args$mat <- args$mat/maxs
#   args$scale <- "none"
#   args$color <- heatmap.colors.max.scaled
#   if (is.na(args$legend_breaks)) {
#     args$legend_breaks <- seq(0, 1, 0.2)
#   }
#   if (is.na(args$breaks)) {
#     args$breaks <- seq(0, 1, length.out = length(args$color)+1)
#   }
#   
#   args
# }

# This next function creates pheatmap annotations_colors dataframe
# list of character vectors, all named.
# vector names = annotation titles
# vector members' (colors') names = annotation identities
# .make_heatmap_annotation_colors <- function(args, annot.colors) {
#   
#   # Extract a default color-set
#   annot.colors.d <- annot.colors
#   annot.colors.n <- rev(annot.colors)
#   
#   # Initiate variables
#   next.color.index.discrete <- 1
#   next.color.index.numeric <- 1
#   col_colors <- NULL
#   row_colors <- NULL
#   user.provided <- names(args$annotation_colors)
#   
#   # Columns First (if there)
#   if (!is.null(args$annotation_col) && ncol(args$annotation_col)>0) {
#     make.these <- !(colnames(args$annotation_col) %in% user.provided)
#     if (any(make.these)) {
#       dfcolors_out <- .pick_colors_for_df(
#         args$annotation_col[,make.these, drop = FALSE],
#         next.color.index.discrete, next.color.index.numeric,
#         annot.colors.d, annot.colors.n)
#       col_colors <- dfcolors_out$df_colors
#       next.color.index.discrete <- dfcolors_out$next.color.index.discrete
#       next.color.index.numeric <- dfcolors_out$next.color.index.numeric
#     }
#   }
#   
#   # Rows Second (if there)
#   if (!is.null(args$annotation_row)) {
#     make.these <- !(colnames(args$annotation_row) %in% user.provided)
#     if (any(make.these)) {
#       dfcolors_out <- .pick_colors_for_df(
#         args$annotation_row[,make.these, drop = FALSE],
#         next.color.index.discrete, next.color.index.numeric,
#         annot.colors.d, annot.colors.n)
#       row_colors <- dfcolors_out$df_colors
#     }
#   }
#   
#   # Combine new with user.provided
#   args$annotation_colors <- c(
#     args$annotation_colors, col_colors, row_colors)
#   
#   args
# }

# Interpret annotations dataframe,
# Pick, name, and add colors.
# .pick_colors_for_df <- function(
#   annotation_df,
#   next.color.index.discrete, next.color.index.numeric,
#   annot.colors.d, annot.colors.n
# ) {
#   
#   df_colors <- list()
#   for (i in seq_len(ncol(annotation_df))){
#     
#     # Make new colors
#     if(!is.numeric(annotation_df[,i])){
#       # Discrete, determine the distinct contents of the annotation first
#       in.this.annot <- levels(as.factor(annotation_df[,i]))
#       
#       # Take colors for each, and name them.
#       new.colors <- annot.colors.d[
#         seq_along(in.this.annot) + next.color.index.discrete - 1
#       ]
#       names(new.colors) <- in.this.annot
#       
#       next.color.index.discrete <-
#         next.color.index.discrete + length(in.this.annot)
#     } else {
#       # Numeric, just need colors for min (white) and max
#       new.colors <-c("white",annot.colors.n[next.color.index.numeric])
#       
#       next.color.index.numeric <- next.color.index.numeric + 1
#     }
#     
#     # Add the new colors as the list
#     df_colors <- c(
#       df_colors,
#       list(new.colors))
#   }
#   names(df_colors) <- names(annotation_df)
#   
#   list(df_colors = df_colors,
#        next.color.index.discrete = next.color.index.discrete,
#        next.color.index.numeric = next.color.index.numeric)
# }

# # Get the heatmap data matrix.
# .get_heatmap_data <- function(object, genes, metas, assay, slot, cells.use) {
#   if (!is.null(genes)) { 
#     data <- as.matrix(.which_data(assay,slot,object)[genes,cells.use]) 
#   } else { 
#     data <- NULL 
#   } 
#   
#   if (!is.null(metas)) { 
#     met.data <- as.matrix(t(getMetas(object, names.only = FALSE)[cells.use, metas])) 
#     data <- rbind(data, met.data) 
#   } 
#   
#   if (is.null(data)) { 
#     stop("No 'genes' or 'metas' requested") 
#   } 
#   
#   if (any(rowSums(data)==0)) { 
#     data <- data[rowSums(data)!=0,] 
#     if (nrow(data)==0) { 
#       stop("No target genes/metadata features have non-zero values in the 'cells.use' subset") 
#     } 
#     warning("Gene(s) or metadata removed due to absence of non-zero values within the 'cells.use' subset") 
#   } 
#   
#   data
# }

dittoHeatmap()