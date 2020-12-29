#' server.R
#' Author: Jack Dodson
#' Michigan State University - Zacharewski Lab

source("./build_dataframes.R")
source("./tzheatmap.R")
source("./tzGSEA.R")

#' The server function. Run once per session.
#' 
#' @param input The server input
#' @param output The server output
#' @param session The server session
shinyServer(function(input, output, session){
  # List possible genes matching input
  output$txt <- renderText({
    # Input must be longer than 2 chars. Dependent on selected annotation type
    if (nchar(input$symbolInput) > 2 && input$annotationInput == "Symbol"){
      sqlcmd <- paste("SELECT DISTINCT Symbol FROM annotation WHERE annotation.Symbol LIKE '%", input$symbolInput, "%'", sep = "")
      gene.list <- unique(toupper(dbGetQuery(dbconn, sqlcmd)[, 1])) 
      gene.list
    }
  })
  
  ### Create the top panel plots (expression data) ###
  
  #' Plot the dose-response data
  #' Initiated by: changing metadata filters, x/Y log scale checkboxes
  getDRPlot <- eventReactive(c(input$enterButton, input$speciesInput, input$sexInput, input$showPlotPreview,
                               input$organInput, input$p1tInput, input$XLogScale, input$YLogScale, 
                               input$plotType, input$treatmentInput, input$assayInput, input$strainInput,
                               input$DRplotExportHeight, input$DRplotExportWidth, input$facetByInput, input$dose_response_toggle_legend),{
    # Builds dataframe
    df <- builddfDoseResponse(input$annotationInput, input$symbolInput, dbconn)
    
    #Filters dataframe based on checkbox inputs
    try(df <- df %>% filter(DesignType %in% c("Dose-response")), TRUE)
    df <- global_filter(df)
    
    # Detects if dataframe is empty
    validate(need(dim(df)[1] != 0, "No data available for selection"))
    
    # Determines P1t threshold
    df$P1t <- ifelse(df$P1t < input$p1tInput, "0", df$Project_ID)
    
    # Determine axis scaling
    XTrans <- "identity"
    YTrans <- "identity"
    if (input$XLogScale){
      XTrans <- "log2"
    }
    if (input$YLogScale){
      YTrans <- "log2"
    }
    
    # Renames columns for easier plot customization
    colnames(df) <- c("Project", "Longname", "Dose", "FoldChange", "Species", "Sex", "Organ", "P1t", "DesignType", "PMID", "GEO", "Chemical", "Assay", "Strain")
    
    # Toggle built in legend
    if(input$dose_response_toggle_legend == TRUE){
      legend_position = "right"
    }
    else{
      legend_position = "none"
    }
    
    # Creates plot(s)
    plot <- ggplot(data = df, mapping = aes(x = Dose, y = FoldChange, text = paste("Project: ", Longname))) +
      geom_line(aes(color = Project, linetype = Organ)) +
      geom_point(aes(color = Project, fill = P1t)) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = c("white", colors)) +
      scale_linetype_manual(values = linetypes) +
      scale_x_continuous(name = "Dose", breaks = unique(df$Dose), trans = XTrans) +
      scale_y_continuous(name = "Fold Change", trans = YTrans) +
      theme_bw() +
      theme(legend.position = legend_position)
    # Faceted plot by selected metadata
    if (input$facetByInput != "None"){
      plot <- ggplot(data = df, mapping = aes(x = Dose, y = FoldChange, text = paste("Project: ", Longname))) +
        geom_line(aes(color = Project, linetype = Organ)) +
        geom_point(aes(color = Project, fill = P1t)) +
        scale_color_manual(values = colors) +
        scale_fill_manual(values = c("white", colors)) +
        scale_linetype_manual(values = linetypes) +
        scale_x_continuous(name = "Dose", breaks = unique(df$Dose), trans = XTrans) +
        scale_y_continuous(name = "Fold Change", trans = YTrans) +
        theme_bw() +
        theme(legend.position = legend_position) +
        facet_wrap(paste0("~", input$facetByInput))
    }
    # Barplot
    if (input$plotType == "Bar"){
      plot <- ggplot(data = df, mapping = aes(x = as.factor(Dose), y = FoldChange, fill = Organ)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_x_discrete(name = "Dose", breaks = unique(df$Dose)) +
        scale_y_continuous(name = "Fold Change") +
        scale_fill_brewer(palette = "Set2") +
        theme_bw() +
        theme(legend.position = legend_position) +
        facet_wrap(~Project)
      p <- ggplotly(plot, tooltip = c("y", "fill"))
    }
    else{
      p <- ggplotly(plot, tooltip = c("y", "shape", "x", "text"))
    }
    
    # Resizes figure
    p <- p %>% layout(autosize = F, width = as.numeric(input$DRplotExportWidth) * 96, height = as.numeric(input$DRplotExportHeight) * 96)
  })
  
  #' Plot the table for hyperlinks to datasets
  #' Initiated by: changing metadata filters
  getDRStudyTable <- eventReactive(c(input$enterButton, input$speciesInput, input$sexInput, input$organInput, 
                                     input$p1tInput, input$treatmentInput, input$assayInput, input$strainInput),{
    # Builds dataframe
    df <- builddfDoseResponse(input$annotationInput, input$symbolInput, dbconn)
    
    # Filters dataframe based on checkbox input
    try(df <- df %>% filter(DesignType %in% c("Dose-response")), TRUE)
    df <- global_filter(df)
    
    # Define which columns to include in completeness score
    QC_cols_to_check <- c('Dose', 'Species_common', 'Sex', 'Organ_name', 'Chemical_Name', 'Assay_Name', 'Strain_Name')
    QC_scores <- c()
    errors <- c()
    # Loop through each project subset of dataframe and calculate completeness
    for(i in unique(df$Project_ID)){
      df_temp <- df %>% filter(Project_ID == i)
      QC_score <- length(QC_cols_to_check)
      for(j in QC_cols_to_check){
        print(paste0('Checking project ', i, ' for ', j))
        temp_col <- df_temp %>% pull(j)
        if(any(is.na(temp_col))){
          QC_score <- QC_score - 1
          error_msg <- paste0('Completeness check failed: Project_ID ', i, ' for data ', j, ' in Dose-response data.')
          errors <- c(errors, error_msg)
          print('failed')
        }
        else{
          print('passed')
        }
      }
      QC_scores <- c(QC_scores, QC_score)
    }
    # Build qc dataframe
    QC_scores_df <- data.frame('Project_ID' = unique(df$Project_ID), 'QC_score' = QC_scores)
    df <- merge(df, QC_scores_df, by.x = 'Project_ID', by.y = 'Project_ID')
    
    # Displays table
    if (dim(df)[1] == 0){
      empty.frame
    }
    else{
      df <- unique(df[, c("Project_ID", "Longname", "PMID", "GEO", "QC_score")])
      df$PMID <- paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed/", df$PMID, "'target='_blank'>", df$PMID, "</a>")
      df$GEO <- paste0("<a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", df$GEO, "'target='_blank'>", df$GEO, "</a>")
      df
    }
  })
  
  #' Create the dose-response plot legend
  #' Initiated by: changing metadata filters
  getDRLegend <- eventReactive(c(input$enterButton, input$speciesInput, input$sexInput, input$organInput, 
                                     input$p1tInput, input$treatmentInput, input$assayInput, input$strainInput),{
    # Builds dataframe
    df <- builddfDoseResponse(input$annotationInput, input$symbolInput, dbconn)
                                     
    # Filters dataframe based on checkbox input
    try(df <- df %>% filter(DesignType %in% c("Dose-response")), TRUE)
    df <- global_filter(df)
    
    # Displays table
    if (dim(df)[1] == 0){
      empty.frame
    }
    else{
      proj_v <- unique(df[, c("Project_ID")])
      organ_v <- unique(df[, c("Organ_name")])
      color_v <- colors[c(1:length(proj_v))]
      linetype_v <- linetypes[c(1:length(organ_v))]
      left <- c(proj_v, organ_v)
      right <- c(color_v, rev(linetype_v))
      df <- data.frame("Data" = left, "Feature" = right)
      df
    }
  })
  
  #' Plot the time-couse data
  #' Initiated by: changing metadata filters, X/Y log scale checkboxes
  getTCPlot <- eventReactive(c(input$enterButton, input$speciesInput, input$sexInput, input$organInput, 
                               input$p1tInput, input$XLogScaleTC, input$YLogScaleTC,
                               input$plotTypeTC, input$treatmentInput, input$assayInput, input$strainInput,
                               input$showPlotPreview, input$TCplotExportHeight, input$TCplotExportWidth,
                               input$facetByInput, input$time_course_toggle_legend),{
    # Builds dataframe
    df <- builddfTimeCourse(input$annotationInput, input$symbolInput, dbconn)
    
    # Filters dataframe based on checkbox inputs
    try(df <- df %>% filter(DesignType %in% c("Time-course", "GWAS")), TRUE)
    df <- global_filter(df)
    
    # Accounts for empty dataframe
    validate(need(dim(df)[1] != 0, "No data available for selection"))
    
    # Determines P1t threshold
    df$P1t <- ifelse(df$P1t < input$p1tInput, "0", df$Project_ID)
    
    # Determine axis scaling
    XTrans <- "identity"
    YTrans <- "identity"
    if (input$XLogScaleTC){
      XTrans <- "log2"
    }
    if (input$YLogScaleTC){
      YTrans <- "log2"
    }
    
    # Rename columns for easier plot customization
    colnames(df) <- c("Project", "Longname", "TimePoint", "FoldChange", "Species", "Sex", "Organ", "P1t", "PMID", "GEO", "DesignType", "Chemical", "Assay", "Strain")
    
    # Toggle built in legend
    if(input$time_course_toggle_legend == TRUE){
      legend_position = "right"
    }
    else{
      legend_position = "none"
    }
    
    # Creates plot
    plot <- ggplot(data = df, mapping = aes(x = TimePoint, y = FoldChange, text = paste("Project: ", Longname))) +
      geom_line(aes(color = Project, linetype = as.factor(Organ))) +
      geom_point(aes(color = Project, shape = as.factor(Species), fill = P1t)) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = c(colors, "white")) +
      scale_linetype_manual(values = linetypes) +
      scale_x_continuous(name = "Time", breaks = unique(df$TimePoint), trans = XTrans) +
      scale_y_continuous(name = "Fold Change", trans = YTrans) +
      theme_bw() +
      theme(legend.position = legend_position)
    # Facetting plot
    if (input$facetByInput != "None"){
      plot <- ggplot(data = df, mapping = aes(x = TimePoint, y = FoldChange, text = paste("Project: ", Longname))) +
        geom_line(aes(color = Project, linetype = as.factor(Organ))) +
        geom_point(aes(color = Project, shape = as.factor(Species), fill = P1t)) +
        scale_color_manual(values = colors) +
        scale_fill_manual(values = c(colors, "white")) +
        scale_linetype_manual(values = linetypes) +
        scale_x_continuous(name = "Time", breaks = unique(df$TimePoint), trans = XTrans) +
        scale_y_continuous(name = "Fold Change", trans = YTrans) +
        theme_bw() +
        theme(legend.position = legend_position) +
        facet_wrap(paste0("~", input$facetByInput))
    }
    # Barplot
    if(input$plotTypeTC == "Bar"){
      plot <- ggplot(data = df, mapping = aes(x = as.factor(TimePoint), y = FoldChange, fill = Organ)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_x_discrete(name = "Dose", breaks = unique(df$TimePoint)) +
        scale_y_continuous(name = "Fold Change") +
        scale_fill_brewer(palette = "Set2") +
        theme_bw() +
        theme(legend.position = legend_position) +
        facet_wrap(~Project)
      p <- ggplotly(plot, tooltip = c("y", "fill"))
    }
    else{
      p <- ggplotly(plot, tooltip = c("y", "shape", "text", "x"))
    }
    
    # Resizes figure
    p <- p %>% layout(autosize = F, width = as.numeric(input$TCplotExportWidth) * 96, height = as.numeric(input$TCplotExportHeight) * 96)
  })
  
  #' Plot the table for hyperlinks to datasets 
  #' Initiated by: changing metadata filters
  getTCStudyTable <- eventReactive(c(input$enterButton, input$speciesInput, input$sexInput, input$organInput, 
                                     input$p1tInput, input$treatmentInput, input$assayInput, input$strainInput),{
    # Builds dataframe
    df <- builddfTimeCourse(input$annotationInput, input$symbolInput, dbconn)
    
    # Filters dataframe based on checkbox input
    try(df <- df %>% filter(DesignType %in% c("Time-course", "GWAS")), TRUE)
    df <- global_filter(df)
    
    # Define which columns to check for completeness
    QC_cols_to_check <- c('TimePoint', 'Species_common', 'Sex', 'Organ_name', 'Chemical_Name', 'Assay_Name', 'Strain_Name')
    QC_scores <- c()
    errors <- c()
    # Loop through each project subset of dataframe and determine completeness for each
    for(i in unique(df$Project_ID)){
      df_temp <- df %>% filter(Project_ID == i)
      QC_score <- length(QC_cols_to_check)
      for(j in QC_cols_to_check){
        print(paste0('Checking project ', i, ' for ', j))
        temp_col <- df_temp %>% pull(j)
        if(any(is.na(temp_col))){
          QC_score <- QC_score - 1
          error_msg <- paste0('Completeness check failed: Project_ID ', i, ' for data ', j, ' in Time-course data.')
          errors <- c(errors, error_msg)
          print('failed')
        }
        else{
          print('passed')
        }
      }
      QC_scores <- c(QC_scores, QC_score)
    }
    # Build completeness dataframe
    QC_scores_df <- data.frame('Project_ID' = unique(df$Project_ID), 'QC_score' = QC_scores)
    df <- merge(df, QC_scores_df, by.x = 'Project_ID', by.y = 'Project_ID')
    
    # Displays table
    if (dim(df)[1] == 0) {
      empty.frame
    } 
    else {
      df <- unique(df[, c("Project_ID", "Longname", "Species_common", "PMID", "GEO", "QC_score")])
      df$PMID <- paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed/", df$PMID, "'target='_blank'>", df$PMID, "</a>")
      df$GEO <- paste0("<a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", df$GEO, "'target='_blank'>", df$GEO, "</a>")
      df
    }
  })
  
  #' Create the dose-response plot legend
  #' Initiated by: changing metadata filters
  getTCLegend <- eventReactive(c(input$enterButton, input$speciesInput, input$sexInput, input$organInput, 
                                 input$p1tInput, input$treatmentInput, input$assayInput, input$strainInput),{
    # Builds dataframe
    df <- builddfTimeCourse(input$annotationInput, input$symbolInput, dbconn)
                                   
    # Filters dataframe based on checkbox input
    try(df <- df %>% filter(DesignType %in% c("Time-course", "GWAS")), TRUE)
    df <- global_filter(df)
                                   
    # Displays table
    if (dim(df)[1] == 0){
      empty.frame
    }
    else{
      proj_v <- unique(df[, c("Project_ID")])
      organ_v <- rev(unique(df[, c("Organ_name")]))
      color_v <- colors[c(1:length(proj_v))]
      linetype_v <- linetypes[c(1:length(organ_v))]
      left <- c(proj_v, organ_v)
      right <- c(color_v, rev(linetype_v))
      df <- data.frame("Data" = left, "Feature" = right)
      df
    }
  })
  
  #' Plot the circadian fold-changes
  #' Initiated by: changing metadata filters
  getCircadianAVG <- eventReactive(c(input$enterButton, input$speciesInput, input$sexInput, input$organInput, input$showPlotPreview, 
                                     input$p1tInput, input$treatmentInput, input$assayInput, input$strainInput,
                                     input$CRTopplotExportWidth, input$CRTopplotExportHeight),{
    # Builds dataframe
    df_orig <- builddfCircadianTop(input$annotationInput, input$symbolInput, dbconn)
    
    # Filters dataframe based on checkbox inputs
    try(df <- df_orig %>% filter(DesignType %in% c("Circadian")), TRUE)
    df <- global_filter(df)
    
    # Accounts for empty dataframe
    validate(need(dim(df)[1] != 0, "No data available for selection"))
    
    # Determines P1t threshold
    df$P1t <- ifelse(df$P1t < input$p1tInput, "0", df$Project_ID)
    
    # Copies datframe for full 48-hour timeframe
    double.ZT <- df$ZT + 24
    double.frame <- df
    double.frame$ZT <- double.ZT
    df <- rbind(df, double.frame)
        
    # Creates plot
    plot <- ggplot(data = df, mapping = aes(x = ZT, y = FoldChange, text = paste("Project: ", Longname))) +
      geom_line(aes(color = Project_ID, linetype = Organ_name, group = Project_ID)) +
      geom_point(aes(color = Project_ID, shape = Organ_name, fill = P1t)) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = c(colors, "white")) +
      scale_x_continuous(name = "ZT", breaks = unique(df$ZT)) +
      scale_y_continuous(name = "Fold Change") +
      theme_bw() +
      theme(legend.position = "none")
    
    
    # Plots figure
    p <- ggplotly(plot, tooltip = c("text", "y", "shape", "x"))
    
    # Resizes figure
    p <- p %>% layout(autosize = F, width = as.numeric(input$CRTopplotExportWidth) * 96, height = as.numeric(input$CRTopplotExportHeight) * 96)
  })
  
  #' Plot the circadian normalized counts
  #' 
  #' Initiated by: changing metadata filters
  getCircadianCounts <- eventReactive(c(input$enterButton, input$speciesInput, input$sexInput, input$organInput, input$p1tInput,
                                        input$showPlotPreview, input$CRBottomplotExportHeight, input$CRBottomplotExportWidth),{
    # Builds dataframe
    df <- builddfCircadianBottom(input$annotationInput, input$symbolInput, dbconn)
    
    # Filters dataframe based on checkbox inputs
    try(df <- df %>% filter(DesignType %in% c("Circadian")), TRUE)
    df <- global_filter(df)
    
    # Accounts for empty dataframe
    validate(need(dim(df)[1] != 0, "No data available for selection"))
    
    # Calculates mean and standard deviation of dataset
    mean <- (tapply(df$NormalizedSignal, list(df$ZT, df$Project_ID, df$Dose), mean))
    sd <- (tapply(df$NormalizedSignal, list(df$ZT, df$Project_ID, df$Dose), sd))
    
    mean <- reshape2::melt(mean)
    colnames(mean) <- c("ZT", "Project_ID", "Dose", "NormalizedSignalMean")
    
    sd <- reshape2::melt(sd)
    colnames(sd) <- c("ZT", "Project_ID", "Dose", "NormalizedSignalSD")
    
    # Merge melted dataframes
    df <- merge(mean, sd, by = c("ZT", "Project_ID", "Dose"))
    
    # Copies datframe for full 48hour timeframe
    double.ZT <- df$ZT + 24
    double.frame <- df
    double.frame$ZT <- double.ZT
    df <- rbind(df, double.frame)
    
    # Convert dose column to a factor
    df$Dose = as.factor(df$Dose)
    
    # Creates plot
    plot <- ggplot(data = df, mapping = aes(x = ZT, y = NormalizedSignalMean, text = paste("Project: ", Project_ID))) +
      geom_line(aes(color = Dose)) +
      geom_point(aes(color = Dose)) +
      geom_errorbar(mapping = aes(ymin = NormalizedSignalMean - NormalizedSignalSD,
                                  ymax = NormalizedSignalMean + NormalizedSignalSD),
                    width = 3,
                    position = position_dodge(0.05)) +
      scale_x_continuous(name = "ZT") +
      scale_y_continuous(name = "Normalized Signal") +
      theme_bw()+
      theme(legend.position = "none")
    
    # Plots figure
    p <- ggplotly(plot, tooltip = c("text", "x", "y"))
    
    # Resizes figure
    p <- p %>% layout(autosize = F, width = as.numeric(input$CRBottomplotExportWidth) * 96, height = as.numeric(input$CRBottomplotExportHeight) * 96)
  })
  
  output$CircadianQCTable <- renderTable({
    df_orig <- builddfCircadianTop(input$annotationInput, input$symbolInput, dbconn)
    
    # Define which columns to check for completeness
    QC_cols_to_check <- c('ZT', 'Organ_name')
    QC_scores <- c()
    errors <- c()
    # Loop through each project subset of the dataframe and determine completeness for each
    for(i in unique(df_orig$Project_ID)){
      df_temp <- df_orig %>% filter(Project_ID == i)
      QC_score <- length(QC_cols_to_check)
      # Check each column for completeness
      for(j in QC_cols_to_check){
        print(paste0('Checking project ', i, ' for ', j))
        temp_col <- df_temp %>% pull(j)
        if(any(is.na(temp_col))){
          QC_score <- QC_score - 1
          error_msg <- paste0('Completeness check failed: Project_ID ', i, ' for data ', j, ' in Circadian regulated data.')
          errors <- c(errors, error_msg)
          print('failed')
        }
        else{
          print('passed')
        }
      }
      QC_scores <- c(QC_scores, QC_score)
    }
    # Write any errors to text file
    fileConn <- file('QC_failures.txt')
    writeLines(errors, fileConn)
    close(fileConn)
    
    # Build final dataframe and display table
    QC_scores_df <- data.frame('Project_ID' = unique(df_orig$Project_ID), 'QC_score' = QC_scores)
    print(QC_scores_df)
  })
  
  #' Plot the normalized expression data
  #' 
  #' Initiated by: changing metadata filters, Y log scale checkbox, groupBy input
  getExprPlot <- eventReactive(c(input$enterButton, input$speciesInput, input$sexInput, input$organInput, input$strainInput, input$p1tInput, 
                                 input$YLogScaleBE, input$groupCountBy, input$showPlotPreview, input$BEplotExportHeight, input$BEplotExportWidth),{
    # Build dataframe
    df <- builddfBE(input$annotationInput, input$symbolInput, dbconn)
    df <- df[, c("Project_ID", "Species_ID", "Sex_ID", "Strain_Name", "NormalizedSignal", "Species_common", "Sex")]
    colnames(df) <- c("Project", "Species_ID", "Sex_ID", "Strain", "NormalizedSignal", "Species", "Sex")
    
    # Filter dataframe
    try(df <- df %>% filter(Species %in% c(input$speciesInput)), TRUE)
    try(df <- df %>% filter(Sex %in% c(input$sexInput)), TRUE)
    try(df <- df %>% filter(Organ_name %in% c(input$organInput)), TRUE)
    try(df <- df %>% filter(Chemical_Name %in% c(input$treatmentInput)), TRUE)
    try(df <- df %>% filter(Assay_Name %in% c(input$assayInput)), TRUE)
    try(df <- df %>% filter(Strain %in% c(input$strainInput)), TRUE)
    
    # Determine axis scaling
    YTrans <- "identity"
    if (input$YLogScaleBE){
      YTrans <- "log"
    }
    
    # Create formula for grouping
    vector <- c(input$groupCountBy)
    form <- as.formula(paste("~interaction(",toString(vector),")"))
    p <- plot_ly(df, x = form, y = ~NormalizedSignal, type = "box")
    p <- layout(p, yaxis = list(type = YTrans))
    
    # Resizes figure
    p <- p %>% layout(autosize = F, width = as.numeric(input$BEplotExportWidth) * 96, height = as.numeric(input$BEplotExportHeight) * 96)
  })
  
  ### Plot UpSet ###
  
  #' Function to create and display venn diagrams
  #' Responsive to plot venn button
  createUpSetPlot <- eventReactive(input$UpSetButton,{
    # Builds the dataframe
    projList <- input$UpSet_selected_projects
    df <- buildUpSetdf(projList, dbconn)
    
    # Require at least 2 projects for comparison
    validate(need(length(projList) >= 2, "Please select at least 2 projects to compare."))
    
    # Apply filters
    significantGenes <- list()
    for (i in 1:length(input$UpSet_selected_projects)){
      subset <- filter(df, Longname == input$UpSet_selected_projects[i])
      if (length(unique(subset$P1t)) == 1){
        filtered <- filter(subset, Dose %in% input[[paste0('DoseFilter', input$UpSet_selected_projects[i])]]
                           & (FoldChange < as.numeric(input[[paste0('FCFilter', input$UpSet_selected_projects[i])]])
                              | FoldChange >= 1/as.numeric(input[[paste0('FCFilter', input$UpSet_selected_projects[i])]]))
                           & Organ_name == input[[paste0('OrganFilter', input$UpSet_selected_projects[i])]])
      }
      else if ((projects %>% filter(projects[,2] == input$UpSet_selected_projects[i]))[,3] == 'Circadian'){
        filtered <- filter(subset, ZT == input[[paste0('DoseFilter', input$UpSet_selected_projects[i])]]
                           & P1t >= input[[paste0('P1tFilter', input$UpSet_selected_projects[i])]]
                           & (FoldChange < as.numeric(input[[paste0('FCFilter', input$UpSet_selected_projects[i])]])
                              | FoldChange >= 1/as.numeric(input[[paste0('FCFilter', input$UpSet_selected_projects[i])]]))
                           & Organ_name == input[[paste0('OrganFilter', input$UpSet_selected_projects[i])]])
      }
      else{
        filtered <- filter(subset, Dose %in% input[[paste0('DoseFilter', input$UpSet_selected_projects[i])]]
                           & P1t >= input[[paste0('P1tFilter', input$UpSet_selected_projects[i])]]
                           & (FoldChange < as.numeric(input[[paste0('FCFilter', input$UpSet_selected_projects[i])]])
                            | FoldChange >= 1/as.numeric(input[[paste0('FCFilter', input$UpSet_selected_projects[i])]]))
                           & Organ_name == input[[paste0('OrganFilter', input$UpSet_selected_projects[i])]])
      }
      significantGenes[[i]] <- sort(unique(filtered$Symbol))
      str <- paste0(substr(projList[i], 1, 20), '...')
      names(significantGenes)[i] <- str
    }
    
    # Plot data
    upset(fromList(significantGenes), order.by = "freq", point.size = 3, line.size = 1, text.scale = 1.5, 
          number.angles = 30, mainbar.y.label = "Project Intersections", sets.x.label = "Genes Per Project",
          set_size.show = FALSE)
  })
  
  # Output each plot in ZED main tab
  output$DRPlot <- renderPlotly({
    getDRPlot()
  })
  output$DRStudyTable <- renderTable({
    getDRStudyTable()
  }, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
  output$TCPlot <- renderPlotly({
    getTCPlot()
  })
  output$TCStudyTable <- renderTable({
    getTCStudyTable()
  }, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
  output$CircadianPlotAVG <- renderPlotly({
    getCircadianAVG()
  })
  output$CircadianCounts <- renderPlotly({
    getCircadianCounts()
  })
  output$seuratPlot <- renderPlot({
    getSeuratPlot()
  })
  output$exprPlot <- renderPlotly({
    getExprPlot()
  })
  output$UpSetPlot <- renderPlot({
    createUpSetPlot()
  })
  
  #' Plot the table for gene references (in ZED)
  #' Initiated by: changing symbol input
  getGeneRefTable <- eventReactive(c(input$enterButton, input$speciesInput),{
    # Convert input to uppercase for SQL query
    gene_symbol_upper <- toupper(input$symbolInput)
    
    # Create list of reference types
    references <- c("GeneCards", "Ensembl - Human", "Ensembl - Mouse", "Ensembl - Rat", "NCBI - Human", "NCBI - Mouse", "NCBI - Rat", "KEGG", "Reactome", "BioGPS")
    
    # Define SQL commands for pulling ENSEMBL and NCBI ID's
    sqlcmd_ens <- paste("SELECT Ensembl_ID FROM Annotation where Annotation.Symbol == '",gene_symbol_upper,"';", sep = "")
    sqlcmd_ncbi <- paste("SELECT NCBI_ID FROM Annotation where Annotation.Symbol == '",gene_symbol_upper,"';", sep = "")
    ensembls <- dbGetQuery(dbconn, sqlcmd_ens)
    ncbis <- dbGetQuery(dbconn, sqlcmd_ncbi)
    
    # Create dataframes of ID's
    ensembls_df <- data.frame(labels = c("Ensembl - Human", "Ensembl - Mouse", "Ensembl - Rat"), ensembls)
    ncbis_df <- data.frame(labels = c("NCBI - Human", "NCBI - Mouse", "NCBI - Rat"), ncbis)
    
    # Build links
    links <- c(paste0("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", input$symbolInput, "'target='_blank'>", input$symbolInput,"</a>"),
               paste0("<a href='https://uswest.ensembl.org/Homo_sapiens/Gene/Summary?g=", ensembls_df$Ensembl_ID[1], "'target='_blank'>", ensembls_df$Ensembl_ID[1],"</a>"),
               paste0("<a href='https://useast.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=", ensembls_df$Ensembl_ID[2], "'target='_blank'>", ensembls_df$Ensembl_ID[2],"</a>"),
               paste0("<a href='https://useast.ensembl.org/Rattus_norvegicus/Gene/Summary?db=core;g=", ensembls_df$Ensembl_ID[3], "'target='_blank'>", ensembls_df$Ensembl_ID[3],"</a>"),
               paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/", ncbis_df$NCBI_ID[1], "'target='_blank'>", ncbis_df$NCBI_ID[1],"</a>"),
               paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/", ncbis_df$NCBI_ID[2], "'target='_blank'>", ncbis_df$NCBI_ID[2],"</a>"),
               paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/", ncbis_df$NCBI_ID[3], "'target='_blank'>", ncbis_df$NCBI_ID[3],"</a>"),
               paste0("<a href='https://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&dbkey=kegg&keywords=", input$symbolInput, "'target='_blank'>", input$symbolInput,"</a>"),
               paste0("<a href='https://reactome.org/content/query?q=", input$symbolInput, "&species=Homo+sapiens&species=Entries+without+species&cluster=true'target='_blank'>", input$symbolInput,"</a>"),
               paste0("<a href='https://biogps.org/#goto=genereport&id=", ncbis_df$NCBI_ID[2], "'target='_blank'>", ncbis_df$NCBI_ID[2],"</a>"))
    
    # Build output dataframe 
    df <- data.frame("Source" = references, "Link" = links)
    
    # Hide data if not selected
    if("human" %in% input$speciesInput == FALSE){
      df <- subset(df, df$Source != "Ensembl - Human" & df$Source != "NCBI - Human")
    }
    if("mouse" %in% input$speciesInput == FALSE){
      df <- subset(df, df$Source != "Ensembl - Mouse" & df$Source != "NCBI - Mouse")
    }
    if("rat" %in% input$speciesInput == FALSE){
      df <- subset(df, df$Source != "Ensembl - Rat" & df$Source != "NCBI - Rat")
    }
    
    # Print dataframe
    df
  })
  
  # Output table populated with gene reference links
  output$geneRefTable <- renderTable({
    getGeneRefTable()
  }, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
  
  # Output legend tables for DR and TC
  output$DRLegend <- renderTable({
    getDRLegend()
  }, hover = TRUE, align = "c", sanitize.text.function = function(x) x)
  output$TCLegend <- renderTable({
    getTCLegend()
  }, hover = TRUE, align = "c", sanitize.text.funtion = function(x) x)
  
  # Show/hide global plot settings in Settings tab
  observeEvent(input$settingsSelected, {
    if("Dose Response" %in% input$settingsSelected){
      shinyjs::show(id = "DRglobalPlotSettings")
    }
    else{
      shinyjs::hide(id = "DRglobalPlotSettings")
    }
    if("Time Course" %in% input$settingsSelected){
      shinyjs::show(id = "TCglobalPlotSettings")
    }
    else{
      shinyjs::hide(id = "TCglobalPlotSettings")
    }
    if("Circadian" %in% input$settingsSelected){
      shinyjs::show(id = "CRglobalPlotSettings")
    }
    else{
      shinyjs::hide(id = "CRglobalPlotSettings")
    }
    if("Basal Expression" %in% input$settingsSelected){
      shinyjs::show(id = "BEglobalPlotSettings")
    }
    else{
      shinyjs::hide(id = "BEglobalPlotSettings")
    }
    if("UpSet" %in% input$settingsSelected){
      shinyjs::show(id = "USglobalPlotSettings")
    }
    else{
      shinyjs::hide(id = "USglobalPlotSettings")
    }
  })
  
  # Show/hide metadata filters in ZED left side panel
  observeEvent(input$filterSelected, {
    for (i in metadata){
      if (i %in% input$filterSelected){
        shinyjs::show(id = i)
      }
      else{
        shinyjs::hide(id = i)
      }
    }
  })
  
  # Show/hide project selected boxes by num selection in UpSet tab
  observeEvent(input$UpSet_selected_projects, {
    for (i in 1:length(projects[,2])){
      if (projects[,2][i] %in% input$UpSet_selected_projects){
        shinyjs::show(id = as.character(projects[,1][i]))
      }
      else{
        shinyjs::hide(id = as.character(projects[,1][i]))
      }
    }
  })
  
  # Hide filter panel if UpSet tab is selected
  observeEvent(input$single_gene_tabs, {
    if (input$single_gene_tabs == 'UpSetTab'){
      shinyjs::hide(id = 'filterPanel')
      shinyjs::hide(id = 'globalP1tSlider')
      shinyjs::hide(id = 'globalGeneSelectPlotButton')
      shinyjs::show(id = 'UpSetPlotButton')
      shinyjs::hide(id = 'globalPlotExportPanel')
      shinyjs::hide(id = 'globalP1tSlider')
      shinyjs::hide(id = 'GeneSelectionPanel')
      shinyjs::hide(id = 'geneRefTable')
      shinyjs::hide(id = 'resetButton')
      shinyjs::hide(id = 'saveLoadBookmarkZED')
    }
    else{
      shinyjs::show(id = 'filterPanel')
      shinyjs::show(id = 'globalP1tSlider')
      shinyjs::show(id = 'globalGeneSelectPlotButton')
      shinyjs::hide(id = 'UpSetPlotButton')
      shinyjs::show(id = 'globalPlotExportPanel')
      shinyjs::show(id = 'globalP1tSlider')
      shinyjs::show(id = 'GeneSelectionPanel')
      shinyjs::show(id = 'geneRefTable')
      shinyjs::show(id = 'resetButton')
      shinyjs::show(id = 'saveLoadBookmarkZED')
    }
    if (input$single_gene_tabs == 'DRtab' || input$single_gene_tabs == 'TCtab'){
      shinyjs::show(id = 'facetOptions')
    }
    else{
      shinyjs::hide(id = 'facetOptions')
    }
  })
  
  # Go to settings tab to edit the plot export size
  observeEvent(input$EditExportSize, {
    updateTabsetPanel(session, "single_gene_tabs", selected = "settingsTab")
  })
  
  # Return to plot tabs after editing settings
  observeEvent(input$return_to_DRTab, {
    updateTabsetPanel(session, 'single_gene_tabs', selected = 'DRtab')
  })
  observeEvent(input$return_to_TCTab, {
    updateTabsetPanel(session, 'single_gene_tabs', selected = 'TCtab')
  })
  observeEvent(input$return_to_CRTab, {
    updateTabsetPanel(session, 'single_gene_tabs', selected = 'CRtab')
  })
  observeEvent(input$return_to_BETab, {
    updateTabsetPanel(session, 'single_gene_tabs', selected = 'BEtab')
  })
  observeEvent(input$return_to_USTab, {
    updateTabsetPanel(session, 'single_gene_tabs', selected = 'UpSetTab')
  })
  
  # Reset plot sizes to default. Dependent on selected tab
  observeEvent(input$ResetToDefault, {
    if(input$single_gene_tabs == 'DRtab'){
      updateTextInput(session, 'DRplotExportWidth', value = '12')
      updateTextInput(session, 'DRplotExportHeight', value = '4.4')
    }
    else if(input$single_gene_tabs == 'TCtab'){
      updateTextInput(session, 'TCplotExportWidth', value = '12')
      updateTextInput(session, 'TCplotExportHeight', value = '4.4')
    }
    else if(input$single_gene_tabs == 'CRtab'){
      updateTextInput(session, 'CRTopplotExportWidth', value = '12')
      updateTextInput(session, 'CRTopplotExportHeight', value = '4')
      updateTextInput(session, 'CRBottomplotExportWidth', value = '12')
      updateTextInput(session, 'CRBottomplotExportHeight', value = '4')
    }
    else if(input$single_gene_tabs == 'BEtab'){
      updateTextInput(session, 'BEplotExportWidth', value = '9')
      updateTextInput(session, 'BEplotExportHeight', value = '4')
    }
  })
 
  ### Gene List Dataviewer (GLD) Tab Panel Code start ###

  # Show/hide heatmaps based on user selection
  observeEvent(input$heatmapsSelected, {
     for (i in projects[,2]){
       if (i %in% input$heatmapsSelected){
         shinyjs::show(id = paste0('heatmap', i))
       }
       else{
         shinyjs::hide(id = paste0('heatmap', i))
       }
     }
  })
  
   # Show/hide heatmap metadata options
   observeEvent(input$heatmap_metadata_selected, {
     for (i in metadata){
       if (i %in% input$heatmap_metadata_selected){
         shinyjs::show(id = paste0('hm', i, 'Filter'))
       }
       else{
         shinyjs::hide(id = paste0('hm', i, 'Filter'))
       }
     }
   })

   # Create file input reactive values
   GLD_fileInput_rv <- reactiveValues(data = NULL)

   # Read in gene list from csv file
   observe({
     req(input$GLD_upload_file_input)
     GLD_fileInput_rv$data <- read.csv(input$GLD_upload_file_input$datapath)
   })

   # Load preset upregulated genes into text input
   observeEvent(input$GLD_upload_preset_up, {
     str <- read.delim("./positive_gene_example.txt")
     print(str)
     str2 <- ""
     for(i in c(1:length(str))){
       if(i == 1){
         str2 <- as.character(colnames(str[i]))
       }
       else{
         str2 <- paste0(str2, "\n", as.character(colnames(str)[i]))
       }
     }
     updateTextInput(session, "GLD_paste_list_input", value = str2)
   })

   # Load preset downregulated genes into text input
   observeEvent(input$GLD_upload_preset_down, {
     str <- read.delim("./negative_gene_example.txt")
     print(str)
     str2 <- ""
     for(i in c(1:length(str))){
       if(i == 1){
         str2 <- as.character(colnames(str[i]))
       }
       else{
         str2 <- paste0(str2, "\n", as.character(colnames(str)[i]))
       }
     }
     updateTextInput(session, "GLD_paste_list_input", value = str2)
   })

   # Clear file upload selection
   observeEvent(input$GLD_reset_file_input, {
     GLD_fileInput_rv$data <- NULL
     reset('GLD_upload_file_input')
   })

   # Create heatmap
   observeEvent(c(input$GLD_submit_list_button, input$heatmapsSelected), {
     if(input$GLDTabs == 'heatmap'){
     # Build gene list from either pasted input or file upload
     if (!is.null(GLD_fileInput_rv$data) & !is.na(input$GLD_paste_list_input)){
       input_list <- read.table(input$GLD_upload_file_input$datapath, sep = '\t')
       gene_list <- vector()
       for (i in input_list){
         gene_list <- c(gene_list, as.character(i))
       }
     }
     else if (!is.null(GLD_fileInput_rv$data)){
       input_list <- read.table(input$GLD_upload_file_input$datapath, sep = '\t')
       gene_list <- vector()
       for (i in input_list){
         gene_list <- c(gene_list, as.character(i))
       }
     }
     else if (!is.na(input$GLD_paste_list_input)){
       gene_list <- as.character(unlist(strsplit(input$GLD_paste_list_input, '\n')))
     }

     gene_list <- unique(gene_list)
     # Pull all data from database
     df <- buildNewHeatmapdf(gene_list, dbconn)
    
     # Create heatmap using lapply for each project
     lapply(1:length(input$heatmapsSelected), function(i) {
       output[[paste0('hm', input$heatmapsSelected[i])]] <- renderPlot({
         # Filter by project (used as check for known missing data)
         df <- filter(df, df$Species_common %in% c(input$heatmap_species_input) &
                          df$Sex %in% c(input$heatmap_sex_input) &
                          df$Organ_name %in% c(input$heatmap_organ_input))

         df <-df[, c('Symbol', 'Dose', 'TimePoint', 'Project_ID', 'FoldChange', 'P1t')]

         df <- filter(df, Project_ID == filter(projects, Longname == input$heatmapsSelected[i])[,1])
         
         validate(need(dim(df)[1] != 0, "No data available for selection"))
         
         # Melt data to be wide instead of long
         if(length(unique(df$Dose)) == 1){
           print(head(df))
           print(unique(df$Dose))
           print(unique(df$ZT))
           fc_wide <- reshape2::dcast(df, Symbol ~ TimePoint, value.var = "FoldChange")
           p1t_wide <- reshape2::dcast(df, Symbol ~ TimePoint, value.var = "P1t")
         }
         else{
           fc_wide <- reshape2::dcast(df, Symbol ~ Dose, value.var = "FoldChange")
           p1t_wide <- reshape2::dcast(df, Symbol ~ Dose, value.var = "P1t")
         }
         
         # Bind FC and P1t data and rename columns for input into heatmap functions
         df_wide <- data.frame(cbind(fc_wide, p1t_wide[2:ncol(p1t_wide)]))
         colnames(df_wide) <- c('Gene', colnames(df_wide)[2:ncol(df_wide)])

         # Determine where FC columns end and P1t columns begin
         cutoff <- (ncol(df_wide) - 1) / 2

         # AhR/DRE data
         ahr.dre.input <- data.frame(Gene = df_wide$Gene, pDRE = 0, AhR = 0)
         AHR <- prepAhR(ahr.dre.input)

         # Prep FC data
         prepped_data <- prepData(gene_list, df_wide, cutoff)

         # Combine data
         combined <- heatmap_combine(gene_list, prepped_data, ahr = AHR)

         # Create heatmap
         hm <- create_heatmap(combined)

         # Output plot
         print(hm)
       })
     })
     }
   })
  
  # Perform GSEA
  observeEvent(input$GLD_submit_list_button, {
    if(input$GLDTabs == 'GSEA'){
    # Build gene list from either pasted input or file upload
    if (!is.null(GLD_fileInput_rv$data) & !is.na(input$GLD_paste_list_input)){
      input_list <- read.table(input$GLD_upload_file_input$datapath, sep = '\t')
      gene_list <- vector()
      for (i in input_list){
        gene_list <- c(gene_list, as.character(i))
      }
    }
    else if (!is.null(GLD_fileInput_rv$data)){
      input_list <- read.table(input$GLD_upload_file_input$datapath, sep = '\t')
      gene_list <- vector()
      for (i in input_list){
        gene_list <- c(gene_list, as.character(i))
      }
    }
    else if (!is.na(input$GLD_paste_list_input)){
      gene_list <- as.character(unlist(strsplit(input$GLD_paste_list_input, '\n')))
    }

    # Build dataframe
    df <- buildGSEAdf(gene_list, dbconn)

    # Run GSEA and sort by NES
    plot_dataframe <- tzGSEA(df, gene_list, input$GSEA_metadata_selected)
    
    # Subset significant pathways for GSEA 4-layer plots
    filtered_path <- subset(plot_dataframe, plot_dataframe$pVal < 0.05)
    filt_p <- as.vector(filtered_path$project_dose)

    # Create GSEA plot objects
    output$plots <- renderUI({
      plot_output_list <- lapply(1:(length(filt_p) + 1), function(i) {
        if (i == 1){
          plotlyOutput("allNES", height = 500, width = 800)
        }
        else{
          plotname <- paste0("plot", filt_p[i - 1])
          plotOutput(plotname, height = 280, width = 500)
        }
      })
      do.call(tagList, plot_output_list)
    })

    # Fill plot objects with GSEA 4-layer plots
    for (i in filt_p) {
      local({
        # Establish local variables
        my_i <- i
        plotname <- paste0("plot", my_i)

        # Split plot title for entry into filter functions below
        project_metadata_vector <- str_split(i, "/", simplify = TRUE)

        # Create renderPlot object
        output[[plotname]] <- renderPlot({
          # Filter by project
          ranked_list <- filter(df, Longname == project_metadata_vector[1])

          # Filter by metadata depending on which is selected
          if('Species' %in% input$GSEA_metadata_selected){
            ranked_list <- filter(ranked_list, Species_common == project_metadata_vector[2])
          }
          else if('Sex' %in% input$GSEA_metadata_selected){
            ranked_list <- filter(ranked_list, Sex == project_metadata_vector[2])
          }
          else if('Organ' %in% input$GSEA_metadata_selected){
            ranked_list <- filter(ranked_list, Organ_name == project_metadata_vector[2])
          }
          else if('Treatment' %in% input$GSEA_metadata_selected){
            ranked_list <- filter(ranked_list, Chemical_Name == project_metadata_vector[2])
          }
          else if('Assay' %in% input$GSEA_metadata_selected){
            ranked_list <- filter(ranked_list, Assay_Name == project_metadata_vector[2])
          }
          else if('Strain' %in% input$GSEA_metadata_selected){
            ranked_list <- filter(ranked_list, Strain_Name == project_metadata_vector[2])
          }
          else if('Dose' %in% input$GSEA_metadata_selected){
            ranked_list <- filter(ranked_list, Dose == project_metadata_vector[2])
          }

          # Sort ranked list, remove duplicates, and resolve ties
          ranked_list <- ranked_list[order(-ranked_list$FoldChange), ]
          ranked_list <- ranked_list[!duplicated(ranked_list$Symbol), ]
          ranked_list <- ranked_list %>% mutate(rank = rank(FoldChange, ties.method = 'random'))

          # Convert ranked list to vector
          ranked_list_vector <- c(ranked_list$FoldChange)
          names(ranked_list_vector) <- ranked_list$Symbol

          # Create gene list from user input
          user_gene_list <- list("mylist" = gene_list)

          # Create unique plot title
          plot_title <- paste0(project_metadata_vector[1], " at condition: ", input$GSEA_metadata_selected, ' ', project_metadata_vector[2])

          # Create plot
          plt <- fgsea::plotEnrichment(user_gene_list[["mylist"]], ranked_list_vector) +
                                       labs(title = plot_title)

          # Display plot
          print(plt)
        })
      })
    }

    # Plot NES barplot
    output[["allNES"]] <- renderPlotly({
      # Create NES plot for all conditions
      plot <- ggplot(plot_dataframe, aes(x = reorder(project_dose, -NES), y = NES)) +
        geom_col(aes(fill=pAdj<0.05)) +
        scale_fill_manual(values = c("#cc2929", "#1b9e35")) +
        labs(x = "", y = "Normalized Enrichment Score") +
        coord_flip() +
        theme_minimal()

      # Convert to plotly
      plot <- ggplotly(plot)

      # Render plot
      plot
    })
    }
  })

  # Reset entire application
  observeEvent(input$reset_button, {js$reset()}) # Call reset method when reset button is pressed
  observeEvent(input$reset_button2, {js$reset()}) # Call reset method when reset button is pressed

  # Save user settings
  # TODO: Devise modular solution so that if new inputs are added, changes to this code are not necessary (priority: LOW)
  observeEvent(input$saveBookmarkZED, {
    inputs_to_save <- c('symbolInput', 'p1tInput', 'filterSelected', 'speciesInput',
                        'sexInput', 'organInput', 'treatmentInput', 'assayInput', 'strainInput',
                        'facetByInput')
    input_types <- c('text', 'slider', 'select', 'check', 'check', 'check', 'check', 'check', 'check', 'select')
    inputs <- list()
    input_value_list <- list()
    for(i in 1:length(inputs_to_save)){
      inputs[[i]] <- input[[inputs_to_save[i]]]
      inputs[[i]] <- c(input_types[i], inputs[[i]])
      input_value_list[[inputs_to_save[i]]] <- inputs[[i]]
    }
    user_file <- paste0("./shiny_bookmarks/", Sys.getenv(c("SHINYPROXY_USERNAME")), ".rds")
    saveRDS(input_value_list, file = user_file)
    print("Settings saved successfully")
  })

  # Load previous user settings
  observeEvent(input$loadBookmarkZED, {
    print("Loading settings...")
    user_file <- paste0("./shiny_bookmarks/", Sys.getenv(c("SHINYPROXY_USERNAME")), ".rds")
    input_value_list <- readRDS(user_file)
    list_names <- list.names(input_value_list)
    for(i in 1:length(input_value_list)){
      print(list_names[i])
      # Depending on input type, update object to match stored settings in rds
      if(input_value_list[[i]][1] == 'text'){
        updateTextInput(session, inputId = list_names[i],
                        value = input_value_list[[i]][2:length(input_value_list[[i]])])
        print(paste0(list_names[i], " loaded successfully"))
      }
      else if(input_value_list[[i]][1] == 'slider'){
        updateSliderInput(session, inputId = list_names[i],
                          value = input_value_list[[i]][2:length(input_value_list[[i]])])
        print(paste0(list_names[i], " loaded successfully"))
      }
      else if(input_value_list[[i]][1] == 'select'){
        updateSelectInput(session, inputId = list_names[i],
                          selected = input_value_list[[i]][2:length(input_value_list[[i]])])
        print(paste0(list_names[i], " loaded successfully"))
      }
      else if(input_value_list[[i]][1] == 'check'){
        updateCheckboxGroupInput(session, inputId = list_names[i],
                          selected = input_value_list[[i]][2:length(input_value_list[[i]])])
        print(paste0(list_names[i], " loaded successfully"))
      }
    }
  })

  ### END GLD CODE ###
  
  ### Begin PCA analysis code ###
  
  #' Build dataset for input into PCA analysis
  #' Activated by 'Update Dataset' button on PCA tab
  observeEvent(input$pca_update_dataset, {
    print('Starting PCA')
    # Read in dataframe
    batchanalysis <- buildBatchAnalysisdf(dbconn)
    print(head(batchanalysis))
    
    # Create unique metadata ID for each row
    batchanalysis$SID = paste("S", batchanalysis$Project, batchanalysis$AID, batchanalysis$Species, batchanalysis$Strain, 
                              batchanalysis$Sex, batchanalysis$Tissue, batchanalysis$Chemical, batchanalysis$Dose, batchanalysis$Timepoint, 
                              batchanalysis$ZT, sep = "_")
    
    # Select only relevant columns and normalize
    batchanalysis_cut = batchanalysis[,c('SID', 'H_ID', 'Count')]
    sumbysample = batchanalysis_cut %>% group_by(SID) %>% summarise(sum = sum(Count))
    
    # dcast long to wide
    batchlong = batchanalysis_cut %>% reshape2::dcast(H_ID ~ SID, sum)
    
    # Remove H_ID column and convert to rownames
    rownames(batchlong) <- batchlong$H_ID
    batchlong$H_ID <- NULL
    
    # Store number of genes to analyze
    num_genes <- input$pca_num_genes_to_analyze
       
    # Ensure that number of genes being analyzed is within bounds
    if(input$pca_num_genes_to_analyze > nrow(batchlong)){
      print("Number of genes selected exceeds those available - defaulting to maximum")
      num_genes <- nrow(batchlong)
    }
    else if(input$pca_num_genes_to_analyze < 10){
      print("Number of genes selected out of range - defaulting to minimum")
      num_genes <- 10
    }
    # Select x most variable genes (as defined by user)
    if(input$pca_select_all_genes_input){
      select = order(rowVars(as.matrix(batchlong)), decreasing = TRUE)
    }
    else{
      select = order(rowVars(as.matrix(batchlong)), decreasing = TRUE)[1:num_genes]
    }
    batchlong <- batchlong[select,]
    
    # Remove irrelevant rows and transpose
    matrix.batchlong <- as.matrix(batchlong)
    row.maximums <- rowMaxs(matrix.batchlong)
    for(i in c(1:length(row.maximums))){
      if(i < 5){
        print(paste0("removing row ", i))
        matrix.batchlong[-c(i),]
      }
    }
    matrix.batchlong <- t(matrix.batchlong)
    print(matrix.batchlong[1:6, 1:6])
    
    data_for_PCA <- matrix.batchlong
    
    # Save files in R space
    saveRDS(data_for_PCA, 'data_for_PCA.rds')
    
    print('PCA Successful')
  })

  #' Create PCA individuals (samples) plot
  #' Responsive to plot button and setting inputs
  createPCAPlot_inds <- eventReactive(c(input$pca_dim1_input, input$pca_dim2_input, input$pca_ellipse_group_by_input, input$pca_show_hide_label_input, input$pca_plot_button, input$pca_show_hide_ellipses_input), {
    if(!file.exists('data_for_PCA.rds')){
      shinyjs::show(id = "pca_no_rds_file_warning")
    }
    else{
      shinyjs::hide(id = "pca_no_rds_file_warning")
    }
    
    # Read in objects
    data_for_PCA <- readRDS('data_for_PCA.rds')
    
    species_v <- c()
    strain_v <- c()
    sex_v <- c()
    tissue_v <- c()
    chemical_v <- c()
    dose_v <- c()
    rows_to_filter <- c()
    for(i in 1:length(rownames(data_for_PCA))){
      str_v <- unlist(strsplit(rownames(data_for_PCA)[i], "_"))
      if(str_v[4] == input$pca_species_input & str_v[5] %in% input$pca_strain_input & str_v[6] %in% input$pca_sex_input &
         str_v[7] %in% input$pca_organ_input & str_v[8] %in% input$pca_chemical_input & str_v[9] %in% input$pca_dose_input){
        species_v <- c(species_v, str_v[4])
        strain_v <- c(strain_v, str_v[5])
        sex_v <- c(sex_v, str_v[6])
        tissue_v <- c(tissue_v, str_v[7])
        chemical_v <- c(chemical_v, str_v[8])
        dose_v <- c(dose_v, str_v[9])
      }
      else{
        rows_to_filter <- c(rows_to_filter, i)
      }
    }
    # Filter data_for_PCA
    data_for_PCA <- data_for_PCA[-rows_to_filter,]
    
    # Perform PCA
    res.pca <- prcomp(as.data.frame(data_for_PCA))
    
    metadata <- data.frame("SPECIES" = species_v,
                           "STRAIN" = strain_v,
                           "SEX" = sex_v,
                           "TISSUE" = tissue_v,
                           "CHEMICAL" = chemical_v,
                           "DOSE" = dose_v)
    
    if(length(input$pca_ellipse_group_by_input) == 1 && length(input$pca_ellipse_group_by_input) != 0){
      grouping_vector <- (metadata[[input$pca_ellipse_group_by_input]])
    }
    else{
      grouping_vector <- metadata[, input$pca_ellipse_group_by_input]
    }
    
    # Set grouping and color palette based on user input and data availability
    if(length(unique(grouping_vector)) == 1 | length(input$pca_ellipse_group_by_input) == 0 | input$pca_show_hide_ellipses_input == FALSE){
      groups = 'black'
      palette = NULL
      ellipse_bool = FALSE
    }
    else{
      if(length(input$pca_ellipse_group_by_input) == 1){
        groups <- as.factor(grouping_vector)
      }
      else{
        table <- grouping_vector
        fac <- table %>% mutate_all(as.character) %>% unite("factor", input$pca_ellipse_group_by_input, remove = TRUE)
        groups = as.factor(fac$factor)
        print(groups)
      }
      palette = mypal[1:length(levels(groups))]
      ellipse_bool = TRUE
    }
    
    # Define which dims to plot based on user
    dimensions <- c(input$pca_dim1_input, input$pca_dim2_input)

    # Show/hide labels based on user
    label = "none"
    if(input$pca_show_hide_label_input){label = "all"}

    # Plot PCA individuals (samples)
    print('plotting pca...')
    plot <- fviz_pca_ind(res.pca,
                axes = dimensions,
                col.ind = groups,
                palette = palette,
                addEllipses = ellipse_bool,
                ellipse.type = "confidence",
                legend.title = "Groups",
                repel = TRUE,
                label = label)
    return(plot)
  })

  #' Create PCA variables (genes) plot
  #' Responsive to plot button and plot settings inputs
  createPCAPlot_vars <- eventReactive(c(input$pca_dim1_input, input$pca_dim2_input, input$pca_ellipse_group_by_input, input$pca_show_hide_label_input, input$pca_plot_button, input$pca_show_hide_ellipses_input), {
    if(!file.exists('data_for_PCA.rds')){
      shinyjs::show(id = "pca_no_rds_file_warning")
    }
    else{
      shinyjs::hide(id = "pca_no_rds_file_warning")
    }
    
    # Read in objects
    data_for_PCA <- readRDS('data_for_PCA.rds')
    
    species_v <- c()
    strain_v <- c()
    sex_v <- c()
    tissue_v <- c()
    chemical_v <- c()
    dose_v <- c()
    rows_to_filter <- c()
    for(i in 1:length(rownames(data_for_PCA))){
      str_v <- unlist(strsplit(rownames(data_for_PCA)[i], "_"))
      if(str_v[4] == input$pca_species_input & str_v[5] %in% input$pca_strain_input & str_v[6] %in% input$pca_sex_input &
         str_v[7] %in% input$pca_organ_input & str_v[8] %in% input$pca_chemical_input & str_v[9] %in% input$pca_dose_input){
        species_v <- c(species_v, str_v[4])
        strain_v <- c(strain_v, str_v[5])
        sex_v <- c(sex_v, str_v[6])
        tissue_v <- c(tissue_v, str_v[7])
        chemical_v <- c(chemical_v, str_v[8])
        dose_v <- c(dose_v, str_v[9])
      }
      else{
        rows_to_filter <- c(rows_to_filter, i)
      }
    }
    # Filter data_for_PCA
    data_for_PCA <- data_for_PCA[-rows_to_filter,]
    
    # Perform PCA
    res.pca <- prcomp(as.data.frame(data_for_PCA))
    
    metadata <- data.frame("SPECIES" = species_v,
                           "STRAIN" = strain_v,
                           "SEX" = sex_v,
                           "TISSUE" = tissue_v,
                           "CHEMICAL" = chemical_v,
                           "DOSE" = dose_v)
    
    if(length(input$pca_ellipse_group_by_input) == 1 && length(input$pca_ellipse_group_by_input) != 0){
      grouping_vector <- (metadata[[input$pca_ellipse_group_by_input]])
    }
    else{
      grouping_vector <- metadata[, input$pca_ellipse_group_by_input]
    }
    
    # Set grouping and color palette based on user input and data availability
    #TODO: fix reactivity and functionality
    if(length(unique(grouping_vector)) == 1 | length(input$pca_ellipse_group_by_input) == 0 | input$pca_show_hide_ellipses_input == FALSE){
      groups = 'black'
      palette = NULL
      ellipse_bool = FALSE
    }
    else{
      if(length(input$pca_ellipse_group_by_input) == 1){
        groups <- as.factor(grouping_vector)
      }
      else{
        table <- grouping_vector
        fac <- table %>% mutate_all(as.character) %>% unite("factor", input$pca_ellipse_group_by_input, remove = TRUE)
        groups = as.factor(fac$factor)
      }
      palette = mypal[1:length(levels(groups))]
      ellipse_bool = TRUE
    }
    
    print(ellipse_bool)
    print(groups)
    
    # Define which dims to plot based on user
    dimensions <- c(input$pca_dim1_input, input$pca_dim2_input)

    # Show/hide labels based on user
    label = "none"
    if(input$pca_show_hide_label_input){label = "all"}

    # Plot PCA variables (genes)
    plot <- fviz_pca_var(res.pca,
                 axes = dimensions,
                 col.ind = groups,
                 palette = palette,
                 addEllipses = ellipse_bool,
                 ellipse.type = "confidence",
                 legend.title = "Groups",
                 repel = TRUE,
                 label = label)
    return(plot)
  })

  #' Create scree plot (barplot of eigenvalues)
  #' Responsive to plot button and plot settings inputs
  createScreePlot <- eventReactive(c(input$pca_dim1_input, input$pca_dim2_input, input$pca_ellipse_group_by_input, input$pca_show_hide_label_input, input$pca_plot_button), {
    # Read in data_for_PCA objec
    data_for_PCA <- readRDS('data_for_PCA.rds')
    
    # Perform PCA
    res.pca <- prcomp(as.data.frame(data_for_PCA))

    # Plot scree
    print(get_eig(res.pca))
    plot <- fviz_eig(res.pca)
    return(plot)
  })

  #' Build table summarizing most variable genes for given dimensions
  #' Responsive to plot button and plot settings inputs
  createPCA_variable_genes_table <- eventReactive(c(input$pca_dim1_input, input$pca_dim2_input, input$pca_ellipse_group_by_input, input$pca_show_hide_label_input, input$pca_plot_button), {
    # Load data from memory and store user selected dimensions
    data_for_PCA <- readRDS('data_for_PCA.rds')
    dimensions <- c(input$pca_dim1_input, input$pca_dim2_input)

    # Perform PCA and extract sorted contribution scores
    res.pca.alt <- PCA(data_for_PCA, scale.unit = FALSE, ncp = 4, graph = FALSE)
    var <- get_pca_var(res.pca.alt)
    dim1_sorted <- res.pca.alt$var$coord[order(-res.pca.alt$var$coord[,dimensions[1]]),]
    dim2_sorted <- res.pca.alt$var$coord[order(-res.pca.alt$var$coord[,dimensions[2]]),]

    # Build genes lists
    dim1_genes <- row.names(dim1_sorted)
    dim2_genes <- row.names(dim2_sorted)

    # Build rank vector
    rank_vector <- c(1:input$pca_top_genes_to_show)
    
    dim1_genes <- homology_IDs[which(homology_IDs$H_ID %in% dim1_genes[1:input$pca_top_genes_to_show]), 2]
    dim2_genes <- homology_IDs[which(homology_IDs$H_ID %in% dim2_genes[1:input$pca_top_genes_to_show]), 2]
    
    # Build final dataframe columns individually
    dim_col <- c(rep(dimensions[1], input$pca_top_genes_to_show), rep(dimensions[2], input$pca_top_genes_to_show))
    genes_col <- c(dim1_genes[1:input$pca_top_genes_to_show], dim2_genes[1:input$pca_top_genes_to_show])
    rank_col <- c(rep(rank_vector, 2))

    # Build and output final dataframe
    table <- data.frame("Dimension" = dim_col, "Gene" = genes_col, "Rank" = rank_col)
    print(table)
  })

  # Output PCA individuals plot
  output$PCA_plot_inds <- renderPlot({
    p <- createPCAPlot_inds()
    p
  })

  # Output PCA variables plot
  output$PCA_plot_vars <- renderPlot({
    p <- createPCAPlot_vars()
    p
  })

  # Output Scree plot
  output$Scree_plot <- renderPlot({
    p <- createScreePlot()
    p
  })

  # Output table summarizing most variable genes in selection (PCA)
  output$PCA_variable_genes_table <- renderTable({
    createPCA_variable_genes_table()
  }, digits = 10, striped = TRUE)

  # Alert user to out of sync plot
  observeEvent(c(input$pca_update_dataset, input$pca_species_input, input$pca_strain_input, input$pca_sex_input, input$pca_organ_input, input$pca_chemical_input, input$pca_dose_input), {
    shinyjs::show(id = "pca_out_of_sync_warning")
  })
  observeEvent(input$pca_plot_button, {
    shinyjs::hide(id = "pca_out_of_sync_warning")
  })
  
  # Inform user if dataset is up to date or not.
  observeEvent(c(input$pca_num_genes_to_analyze, input$pca_select_all_genes_input), {
    shinyjs::show(id = "pca_dataset_not_up_to_date")
    shinyjs::hide(id = "pca_dataset_up_to_date")
  })
  observeEvent(c(input$pca_update_dataset), {
    shinyjs::hide(id = "pca_dataset_not_up_to_date")
    shinyjs::show(id = "pca_dataset_up_to_date")
  })
  
  # Show/hide metadata filters
  observeEvent(input$filterPCASelected, {
    for (i in c("species", "strain", "sex", "organ", "chemical", "dose")){
      if (i %in% input$filterPCASelected){
        shinyjs::show(id = paste0('PCA_', i, '_select'))
      }
      else{
        shinyjs::hide(id = paste0('PCA_', i, '_select'))
      }
    }
  })
  
  ### End PCA analysis code ###

  ### Begin PLS code ###

  #' Build dataset for input into PLS analysis
  #' Activated by 'Update Dataset' button on PLS tab
  loadDataset <- eventReactive(input$pls_update_dataset, {
    print('Starting PLS')
    # Read in dataframe
    batchanalysis <- buildBatchAnalysisdf(dbconn)
    print(head(batchanalysis))
    
    # Create unique metadata ID for each row
    batchanalysis$SID = paste("S", batchanalysis$Project, batchanalysis$AID, batchanalysis$Species, batchanalysis$Strain, 
                              batchanalysis$Sex, batchanalysis$Tissue, batchanalysis$Chemical, batchanalysis$Dose, batchanalysis$Timepoint, 
                              batchanalysis$ZT, sep = "_")
    
    # Select only relevant columns and normalize
    batchanalysis_cut = batchanalysis[,c('SID', 'H_ID', 'Count')]
    sumbysample = batchanalysis_cut %>% group_by(SID) %>% summarise(sum = sum(Count))
    
    # dcast long to wide
    batchlong = batchanalysis_cut %>% reshape2::dcast(H_ID ~ SID, sum)
    
    # Remove H_ID column and convert to rownames
    rownames(batchlong) <- batchlong$H_ID
    batchlong$H_ID <- NULL
    
    # Store number of genes to analyze
    num_genes <- input$pls_num_genes_to_analyze
    
    # Ensure that number of genes being analyzed is within bounds
    if(input$pls_num_genes_to_analyze > nrow(batchlong)){
      print("Number of genes selected exceeds those available - defaulting to maximum")
      num_genes <- nrow(batchlong)
    }
    else if(input$pls_num_genes_to_analyze < 10){
      print("Number of genes selected out of range - defaulting to minimum")
      num_genes <- 10
    }
    
    # Select x most variable genes (as defined by user)
    if(input$pls_select_all_genes_input){
      select = order(rowVars(as.matrix(batchlong)), decreasing = TRUE)
    }
    else{
      select = order(rowVars(as.matrix(batchlong)), decreasing = TRUE)[1:num_genes]
    }
    batchlong <- batchlong[select,]
    
    # Remove irrelevant rows and transpose
    matrix.batchlong <- as.matrix(batchlong)
    row.maximums <- rowMaxs(matrix.batchlong)
    for(i in c(1:length(row.maximums))){
      if(i < 5){
        print(paste0("removing row ", i))
        matrix.batchlong[-c(i),]
      }
    }
    matrix.batchlong <- t(matrix.batchlong)
    print(matrix.batchlong[1:6, 1:6])
    
    # Save data in R space
    saveRDS(matrix.batchlong, 'data_for_PLS.rds')
    
    print('PLS Successful')
  })

  #' Create plot of individuals
  #' Activated by plot button and plot settings inputs
  createPLSPlot_indivs <- eventReactive(c(input$pls_plot_button, input$pls_ellipse_group_by_input, input$pls_point_labels), {
    if(!file.exists('data_for_PLS.rds')){
      shinyjs::show(id = "pls_no_rds_file_warning")
    }
    else{
      shinyjs::hide(id = "pls_no_rds_file_warning")
    }
    
    # Read in dataset and file list from above
    print("loading data for PLS...")
    data_for_PLS <- readRDS('data_for_PLS.rds')
    
    species_v <- c()
    strain_v <- c()
    sex_v <- c()
    tissue_v <- c()
    chemical_v <- c()
    dose_v <- c()
    rows_to_filter <- c()
    for(i in 1:length(rownames(data_for_PLS))){
      str_v <- unlist(strsplit(rownames(data_for_PLS)[i], "_"))
      if(str_v[4] == input$pls_species_input & str_v[5] %in% input$pls_strain_input & str_v[6] %in% input$pls_sex_input &
         str_v[7] %in% input$pls_organ_input & str_v[8] %in% input$pls_chemical_input & str_v[9] %in% input$pls_dose_input){
        species_v <- c(species_v, str_v[4])
        strain_v <- c(strain_v, str_v[5])
        sex_v <- c(sex_v, str_v[6])
        tissue_v <- c(tissue_v, str_v[7])
        chemical_v <- c(chemical_v, str_v[8])
        dose_v <- c(dose_v, str_v[9])
      }
      else{
        rows_to_filter <- c(rows_to_filter, i)
      }
    }
    # Filter data_for_PCA
    data_for_PLS <- data_for_PLS[-rows_to_filter,]
    
    metadata <- list("SPECIES" = species_v,
                     "STRAIN" = strain_v,
                     "SEX" = sex_v,
                     "TISSUE" = tissue_v,
                     "CHEMICAL" = chemical_v,
                     "DOSE" = dose_v)
    
    print(head(metadata))
    
    # Ensure data is compatible with PLS
    validate(need(length(unique(metadata[[input$pls_ellipse_group_by_input]])) >= 2, print("Grouping must have more than 1 level")))
  
    # Run PLS and plot
    print("plotting pls...")
    MyResult.splsda <- splsda(data_for_PLS, metadata[[input$pls_ellipse_group_by_input]], keepX = c(300,300))
    plotIndiv(MyResult.splsda, ind.names = metadata[[input$pls_point_labels]],
              group = as.factor(as.character(metadata[[input$pls_ellipse_group_by_input]])), legend=TRUE, cex = 4,
              ellipse = TRUE, ellipse.level = 0.95, star = FALSE,
              X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

    output$PLS_top_genes <- renderPlot({
      plotLoadings(MyResult.splsda, contrib = 'max', method = 'mean', ndisplay = 50)
    })
  })

  # Render PLS individuals plot
  output$PLS_indivs <- renderPlot({
    loadDataset()
    p <- createPLSPlot_indivs()
    p
  })

  # Alert user to out of sync plot
  observeEvent(input$pls_update_dataset, {
    shinyjs::show(id = "pls_out_of_sync_warning")
  })
  observeEvent(input$pls_plot_button, {
    shinyjs::hide(id = "pls_out_of_sync_warning")
  })
  
  # Show/hide metadata filters
  observeEvent(input$filterPLSSelected, {
    for (i in c("species", "strain", "sex", "organ", "chemical", "dose")){
      if (i %in% input$filterPLSSelected){
        shinyjs::show(id = paste0('PLS_', i, '_select'))
      }
      else{
        shinyjs::hide(id = paste0('PLS_', i, '_select'))
      }
    }
  })
  
  ### End PLS code ###
  
  ### Begin Single Cell ###
  
  # Create file input reactive values
  singlecell_fileInput_rv <- reactiveValues(data = NULL)
  
  # Read in gene list from csv file
  observe({
    req(input$singlecell_upload_file_input)
    singlecell_fileInput_rv$data <- read.csv(input$singlecell_upload_file_input$datapath)
  })
  
  # Load preset upregulated genes into text input
  observeEvent(input$singlecell_upload_preset_up, {
    str <- read.delim("./positive_gene_example.txt")
    for(i in c(1:length(str))){
      if(i == 1){
        str2 <- as.character(colnames(str[i]))
      }
      else{
        str2 <- paste0(str2, "\n", as.character(colnames(str)[i]))
      }
    }
    updateTextInput(session, "singlecell_paste_list_input", value = str2)
  })
  
  # Load preset downregulated genes into text input
  observeEvent(input$singlecell_upload_preset_down, {
    str <- read.delim("./negative_gene_example.txt")
    for(i in c(1:length(str))){
      if(i == 1){
        str2 <- as.character(colnames(str[i]))
      }
      else{
        str2 <- paste0(str2, "\n", as.character(colnames(str)[i]))
      }
    }
    updateTextInput(session, "singlecell_paste_list_input", value = str2)
  })
  
  # Clear file upload selection
  observeEvent(input$singlecell_reset_file_input, {
    singlecell_fileInput_rv$data <- NULL
    reset('singlecell_upload_file_input')
  })
  
  create_singlecell_genelist <- function(){
    # Build gene list from either pasted input or file upload
    if (!is.null(singlecell_fileInput_rv$data) & !is.na(input$singlecell_paste_list_input)){
      input_list <- read.table(input$singlecell_upload_file_input$datapath, sep = '\t')
      gene_list <- vector()
      for (i in input_list){
        gene_list <- c(gene_list, as.character(i))
      }
    }
    else if (!is.null(singlecell_fileInput_rv$data)){
      input_list <- read.table(input$singlecell_upload_file_input$datapath, sep = '\t')
      gene_list <- vector()
      for (i in input_list){
        gene_list <- c(gene_list, as.character(i))
      }
    }
    else if (!is.na(input$singlecell_paste_list_input)){
      gene_list <- as.character(unlist(strsplit(input$singlecell_paste_list_input, '\n')))
    }
    return(gene_list)
  }
  
  # Create UMAP plot
  createUMAP <- eventReactive(input$plot_singlecell, {
    sn <- get(paste0("sn_", input$singlecell_dataset_input))
    # Merge cluster and metadata
    umap_dataframe <- merge(sn$umap, sn$meta, by = "NAME")
    
    if(input$UMAP_show_num_genes == 500){
      umap_dataframe <- umap_dataframe[1:500,]
    }
    
    # Plot
    if(input$UMAP_label_by_input == "None"){
      plot <- ggplot(umap_dataframe, aes(x = as.numeric(as.vector(X)), y = as.numeric(as.vector(Y)), color = as.factor(.data[[input$UMAP_color_by_select]]))) + 
        geom_point(size = 0.1) +
        xlab("UMAP1") +
        ylab("UMAP2") +
        labs(color = as.character(input$UMAP_color_by_select)) +
        theme_bw()
    }
    else{
      plot <- ggplot(umap_dataframe, aes(x = as.numeric(as.vector(X)), y = as.numeric(as.vector(Y)), color = as.factor(.data[[input$UMAP_color_by_select]]), label = as.factor(.data[[input$UMAP_label_by_input]]))) + 
              geom_point(size = 0.1) +
              geom_text(check_overlap = TRUE) +
              xlab("UMAP1") +
              ylab("UMAP2") +
              labs(color = as.character(input$UMAP_color_by_select)) +
              theme_bw()
    }
    
    plot <- ggplotly(plot)
    
    return(plot)
  })
  
  # Create Feature plot
  createFeaturePlot <- eventReactive(input$plot_singlecell, {
    sn <- get(paste0("sn_", input$singlecell_dataset_input))
    # Read in gene list from user
    gene_list <- create_singlecell_genelist()
    
    # Get all barcoded cells for genes selected and rename columns for merge
    gene_df <- getGeneData(c(gene_list), sn$barcodes.order, sn)
    colnames(gene_df) <- c("GENE", "NAME", "VALUE")
    # Merge on barcode
    gene_df <- merge(gene_df, sn$umap, by = "NAME")
    gene_df <- merge(gene_df, sn$meta, by = "NAME")
    
    # Create Plot
    plot <- ggplot(gene_df) + geom_point(aes(x = as.numeric(as.vector(X)), y = as.numeric(as.vector(Y)), color = as.numeric(as.vector(VALUE))), size = 0.1) +
            scale_color_gradient2(low = 'lightgrey', mid = 'grey', high = '#1908ff') +
            facet_wrap(~as.factor(.data[[input$singlecell_metadata_select_feature]])) +
            xlab("UMAP1") +
            ylab("UMAP2") +
            labs(color = "Normalized Expression") +
            theme_bw()
    
    plot <- ggplotly(plot)
    
    return(plot)
  })
  
  # Create Ridge plot
  createRidgePlot <- eventReactive(input$plot_singlecell, {
    sn <- get(paste0("sn_", input$singlecell_dataset_input))
    # Read in gene list from user
    gene_list <- create_singlecell_genelist()
    
    # Get all barcoded cells for genes selected and rename columns for merge
    gene_df <- getGeneData(c(gene_list), sn$barcodes.order, sn)
    colnames(gene_df) <- c("GENE", "NAME", "VALUE")
    # Merge on barcode
    gene_df <- merge(gene_df, sn$meta, by = "NAME")

    # Create plot
    plot <- ggplot(gene_df, aes(x = as.numeric(as.vector(VALUE)), y = as.factor(.data[[input$singlecell_metadata_select]]), fill = as.factor(.data[[input$singlecell_metadata_select]]))) + 
            geom_density_ridges(show.legend = FALSE) +
            facet_wrap(~GENE) + 
            xlab("Normalized Expression") + 
            ylab("Cell Type") +
            theme_bw()
    
    return(plot)
  })
  
  # Create Violin plot
  createViolinPlot <- eventReactive(input$plot_singlecell, {
    sn <- get(paste0("sn_", input$singlecell_dataset_input))
    # Read in gene list from user
    gene_list <- create_singlecell_genelist()
    
    # Get all barcoded cells for genes selected and rename columns for merge
    gene_df <- getGeneData(c(gene_list), sn$barcodes.order, sn)
    colnames(gene_df) <- c("GENE", "NAME", "VALUE")
    # Merge on barcode
    gene_df <- merge(gene_df, sn$meta, by = "NAME")
    
    # Create Plot
    plot <- ggplot(gene_df, aes(x = as.factor(.data[[input$singlecell_metadata_select]]), y = as.numeric(as.vector(VALUE)), color = as.factor(.data[[input$singlecell_metadata_select]]))) +
      geom_violin() + 
      geom_jitter(shape = 16, position = position_jitter(0.2)) +
      facet_wrap(~GENE) +
      xlab("Cell Type") +
      ylab("Normalized Expression") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none")
    
    return(plot)
  })
  
  # Create Dot plot
  createDotPlot <- eventReactive(input$plot_singlecell, {
    sn <- get(paste0("sn_", input$singlecell_dataset_input))
    # Read in gene list from user
    gene_list <- create_singlecell_genelist()
    
    gene_df <- getGeneData(c(gene_list), sn$barcodes.order, sn)
    colnames(gene_df) <- c("GENE", "NAME", "VALUE")
    # Merge on barcode
    gene_df <- merge(gene_df, sn$meta, by = "NAME", all.y = TRUE)
    
    # Convert long to wide
    data_wide <- spread(gene_df[,c("NAME","GENE","VALUE")], NAME, VALUE)
    
    # Convert wide to long and replace NA's with 0
    data_long <- melt(data_wide)
    data_long[is.na(data_long)] <- 0
    
    # Rename columns and merge on metadata
    colnames(data_long) <- c("GENE", "NAME", "VALUE")
    data_long <- as.data.frame(merge(data_long, sn$meta, by = "NAME"))
    
    # Count number of zeros and number of nonzeros grouped by gene and celltype 
    nonzeros <- data_long %>% 
      group_by(GENE, across(all_of(input$singlecell_metadata_select))) %>% 
      summarise(percent_expressed = (1 - (sum(VALUE == 0)/(sum(VALUE > 0) + sum(VALUE == 0)))) * 100, avg_value = mean(VALUE))
    
    nonzeros <- drop_na(nonzeros)
    
    print(head(nonzeros))
    
    saveRDS(nonzeros, 'C:\\Users\\Jack\\Desktop\\nonzeros.rds')
    
    # Create plot
    plot <- ggplot(data = nonzeros, mapping = aes_string(x = 'GENE', y = 'celltype')) +
      geom_point(mapping = aes_string(size = 'percent_expressed', color = "avg_value")) +
      #scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank())
      #guides(size = guide_legend(title = 'Percent Expressed'))
    
    plot <- ggplotly(plot)
    
    return(plot)
  })
  
  # Output UMAP
  output$UMAP <- renderPlotly({
    plot <- createUMAP()
    plot
  })
  
  # Output Feature plot
  output$FeaturePlot <- renderPlotly({
    plot <- createFeaturePlot()
    plot
  })
  
  # Output Ridge plot
  output$RidgePlot <- renderPlot({
    plot <- createRidgePlot()
    plot
  },)
  
  # Output Violin plot
  output$ViolinPlot <- renderPlot({
    plot <- createViolinPlot()
    plot
  })
  
  # Output Dot plot
  output$DotPlot <- renderPlotly({
    plot <- createDotPlot()
    plot
  })
  
  sc_description_maker <- eventReactive(c(input$singlecell_species_select, input$singlecell_sex_select),{
    choices <- sc_dataset_meta %>% filter(Species_common %in% input$singlecell_species_select, Sex %in% input$singlecell_sex_select)
    result <- choices$Project_ID
    return(result)
  })
  
  output$singlecell_dataset_choices <- renderTable({
    result <- sc_description_maker()
    print(result)
    lapply(1:length(result), function(i) {
      output[[paste0('sc_description', i)]] <- renderUI({
        includeHTML(paste0("C:\\Users\\Jack\\Desktop\\FAIRTox_github\\app\\www\\", result[i], ".html"))
      })
    })
  })
  
  # Hide/show divs based on what tab is selected
  observeEvent(input$singlecell_tabs, {
    if(input$singlecell_tabs == "UMAP"){
      shinyjs::hide(id = "singlecell_genelist_left_panel")
      shinyjs::show(id = "singlecell_right_panel_UMAP")
      shinyjs::hide(id = "singlecell_right_panel_all_except_feature")
      shinyjs::hide(id = "singlecell_right_panel_feature")
    }
    else if(input$singlecell_tabs == "Feature Plot"){
      shinyjs::show(id = "singlecell_genelist_left_panel")
      shinyjs::hide(id = "singlecell_right_panel_UMAP")
      shinyjs::hide(id = "singlecell_right_panel_all_except_feature")
      shinyjs::show(id = "singlecell_right_panel_feature")
    }
    else{
      shinyjs::show(id = "singlecell_genelist_left_panel")
      shinyjs::hide(id = "singlecell_right_panel_UMAP")
      shinyjs::show(id = "singlecell_right_panel_all_except_feature")
      shinyjs::hide(id = "singlecell_right_panel_feature")
    }
  })
  
  ### End Single Cell ###
  
  # autoInvalidate <- reactiveTimer(2000)
  # 
  # env <- globalenv()  # can use globalenv(), parent.frame(), etc
  # output$foo <- renderTable({
  #   autoInvalidate()
  #   table <- data.frame(
  #     object = c(ls(env), "total"),
  #     size = c(unlist(lapply(ls(env), function(x) {
  #       object.size(get(x, envir = env, inherits = FALSE))
  #     })), lobstr::mem_used())
  #   )
  # })
  
  ### Integration of Tableau Viz testing ###
  ### Not implemented in current version ###
  # js$init()
  
  #' Filter dataframe based on user input - used in most ZED plots
  #' 
  #' @param df The dataframe to filter
  #' @return The filtered dataframe
  #' @examples
  #' global_filter(dataframe_1)
  global_filter <- function(df){
    try(df <- df %>% filter(Species_common %in% c(input$speciesInput)), TRUE)
    try(df <- df %>% filter(Sex %in% c(input$sexInput)), TRUE)
    try(df <- df %>% filter(Organ_name %in% c(input$organInput)), TRUE)
    try(df <- df %>% filter(Chemical_Name %in% c(input$treatmentInput)), TRUE)
    try(df <- df %>% filter(Assay_Name %in% c(input$assayInput)), TRUE)
    try(df <- df %>% filter(Strain_Name %in% c(input$strainInput)), TRUE)
    return(df)
  }
})
