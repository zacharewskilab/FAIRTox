#' ui.R
#' Author: Jack Dodson
#' Michigan State University - Zacharewski Lab

function(request){
    fluidPage(theme = "./www/bootstrap.css",
    navbarPage("Navigation",
    tabsetPanel(id = "mainTabs", type = "tabs",
        tabPanel("Home", includeHTML("./www/homePage.html")),
        tabPanel("README", includeMarkdown("./README.md")),
        tabPanel("Tutorial", includeHTML("./www/tutorialPage.html")),
        ## Youtube video tutorial (coming soon)
        #tabPanel("Embed Youtube Testing Area", includeHTML("./www/embedYoutubeTest.html")),
        tabPanel("Single Gene Query",
            title = "Single Gene Query", id = "SingleGeneQuery",
            fluidRow
            (
              column(3, 
                wellPanel(
                id = "leftPanel", style = "background: #919191",
                div(id = "GeneSelectionPanel",
                  tags$p("Search Gene:"),
                  div(id = "globalGeneSelectPlotButton",
                      column(12, textInput("symbolInput", label = NULL, value = "Cyp1a1")),
                      fluidRow(column(6, actionButton("enterButton", "Plot", style = "color: white; background-color: #278edb")))
                      ), 
                  br(),
                  "Possible Genes:", br(),
                  textOutput("txt"), br(),
                  radioButtons("annotationInput", "Annotation source:", 
                               choices = c("Symbol", "Ensembl_ID", "NCBI_ID"), selected = "Symbol", inline = TRUE)
                ),
                div(id = "geneRefTable",
                    "Gene Information:",
                    tableOutput("geneRefTable")
                ),
                div(id = "globalP1tSlider", sliderInput("p1tInput", "P1(t) threshold:", 0, 1, 0.8, ticks = FALSE)),
                div(id = "globalPlotExportPanel",
                  HTML('<hr background-color="#333" color = "#333">'),
                  div(id = "plotPreviewPanel",
                      actionButton("EditExportSize", "Edit export size"),
                      actionButton("ResetToDefault", "Reset to default")                      
                  ),
                  HTML('<hr background-color="#333" color = "#333">')
                ),
                div(id = 'UpSetPlotButton', actionButton("UpSetButton", "Plot", style = "color: white; background-color: #278edb"),
                    style = 'position:fixed; width:inherit;'),
                div(id = "filterPanel",
                  useShinyjs(),
                  "Metadata Filters",
                  selectInput("filterSelected", "Show/hide filters:", 
                              choices = metadata, multiple = TRUE, 
                              selected = c("Species")),
                  div(id = "Species", checkboxGroupInput("speciesInput", "Species:", species, selected = species)), br(),
                  div(id = "Sex", checkboxGroupInput("sexInput", "Sexes:", sex, select = sex)), br(),
                  div(id = "Organ", checkboxGroupInput("organInput", "Target Organs:", unique(organs$Organ_name), select = unique(organs$Organ_name))), br(),
                  div(id = "Treatment", selectInput("treatmentInput", "Treatment Chemical:", treatments, selected = "TCDD")), br(),
                  div(id = "Assay", checkboxGroupInput("assayInput", "Assays:", assays, selected = assays)), br(),
                  div(id = "Strain", checkboxGroupInput("strainInput", "Strains", strains, selected = strains))
                ),
                div(id = "facetOptions",
                  "Facet Options",
                  selectInput("facetByInput", "Facet by:",
                              choices = c("None", "Project", "Species", "Sex", "Organ", "Chemical", "Assay", "Strain"), multiple = FALSE)
                ),
                div(id = "resetButton",
                  #useShinyjs(),
                  #extendShinyjs(text = jsResetCode),
                  actionButton("reset_button", "Reset Page")
                ),
                div(id = "saveLoadBookmarkZED",
                  actionButton("saveBookmarkZED", "Save ZED settings"),
                  actionButton("loadBookmarkZED", "load ZED settings")
                )
               )
              ),
              column(9,
                textOutput('single_gene_tabs'),
                tabsetPanel
                (
                  id = "single_gene_tabs", type = "tabs",
                  tabPanel(
                    title = "Dose-response data", id = "DRtab", value = "DRtab", br(), 
                    withSpinner(plotlyOutput("DRPlot")), br(), 
                    column(8,
                      "Please cite the following:", br(), 
                      tableOutput("DRStudyTable"),
                      "Note: Maximum QC score for this dataset is 7."
                    ),
                    column(4,
                        div(id = "doseResponseLegend",
                            "Legend", br(), tableOutput("DRLegend"),
                            checkboxInput(inputId = "dose_response_toggle_legend", label = "Toggle built in legend", value = FALSE)
                        )
                    )
                  ),
                  tabPanel(
                    title = "Time-course data", id = "TCtab", value = "TCtab", br(), 
                    withSpinner(plotlyOutput("TCPlot")), br(), 
                    column(8,
                           "Please cite the following:", br(), 
                           tableOutput("TCStudyTable"),
                           "Note: Maximum QC score for this dataset is 7."
                    ),
                    column(4,
                        div(id = "doseResponseLegend",
                            "Legend", br(), tableOutput("TCLegend"),
                            checkboxInput(inputId = "time_course_toggle_legend", label = "Toggle built in legend", value = FALSE)
                        )
                    )
                  ),
                  tabPanel(title = "Circadian regulated data", id = "CRtab", value = "CRtab", br(),
                           withSpinner(plotlyOutput("CircadianPlotAVG")), br(), 
                           withSpinner(plotlyOutput("CircadianCounts")),
                           withSpinner(tableOutput("CircadianQCTable")),
                           'Note: Maximum QC for dataset is 2 (WIP).'
                           ),
                  tabPanel(title = "Basal expression levels", id = "BEtab", value = "BEtab",
                           column(8, wellPanel(withSpinner(plotlyOutput("exprPlot")))), 
                           column(4, wellPanel(selectInput("groupCountBy", "Group by:", choices = c("Project", "Sex", "Species", "Strain"), selected = c("Project"), multiple = TRUE)))),
                  tabPanel(title = "UpSet Plot", id = "UpSetTab", value = "UpSetTab",
                           useShinyjs(),
                           br(),
                           tags$style("#UpSetPlot {background-color:white;}"),
                           div("UpSetPlot", column(8, wellPanel(withSpinner(plotOutput("UpSetPlot", height = "600px", width = "100%"))))),
                           column(4, wellPanel(
                             checkboxGroupInput(inputId = "UpSet_selected_projects", "Select projects to compare:", choices = projects[,2]), 
                             br(),
                             lapply(1:length(projects[,1]), function(i) {
                               div(id = paste0(projects[,1][i]),
                                   tags$h3(paste0(projects[,2][i])),
                                   checkboxGroupInput(inputId = paste0('DoseFilter', projects[,2][i]), "Select the doses to observe:", 
                                                      ifelse((doses %>% filter(doses[,2] == projects[,1][i]))[,4] == "Circadian",
                                                             sort((doses %>% filter(doses[,2] == projects[,1][i]))[,3]), sort((doses %>% filter(doses[,2] == projects[,1][i]))[,1])), 
                                                      ifelse((doses %>% filter(doses[,2] == projects[,1][i]))[,4] == "Circadian",
                                                             sort((doses %>% filter(doses[,2] == projects[,1][i]))[,3]), sort((doses %>% filter(doses[,2] == projects[,1][i]))[,1]))),
                                   textInput(inputId = paste0('FCFilter', projects[,2][i]), "Enter the fold change threshold to observe:", value = 1.0),
                                   sliderInput(inputId = paste0('P1tFilter', projects[,2][i]), "P1(t) threshold:", 0, 1, 0.8),
                                   checkboxGroupInput(inputId = paste0('OrganFilter', projects[,2][i]), "Select the tissues to observe:", 
                                                      (organs %>% filter(organs[,2] == projects[,1][i]))[,1], (organs %>% filter(organs[,2] == projects[,1][i]))[,1]),
                                   HTML('<hr background-color="#333" color = "#333">')
                                   )
                             })
                             )
                          )
                    ),
                    tabPanel(title = "Settings", id = "settingsTab", value = "settingsTab",
                           useShinyjs(),
                           br(),
                           selectInput("settingsSelected", "Select Settings to Show:", multiple = TRUE,
                                       choices = c("Dose Response", "Time Course", "Circadian", "Basal Expression", "UpSet"),
                                       selected = c("Dose Response")),
                           br(),
                           div(id = "DRglobalPlotSettings",
                               tags$h1("Dose Response Plot Settings"),
                               column(12, selectInput("plotType", "Select plot type:", c("Line" = "Line", "Bar" = "Bar"))),
                               column(12,
                                      column(3, checkboxInput("XLogScale", "Toggle X axis log scale", value = TRUE)),
                                      column(3, checkboxInput("YLogScale", "Toggle Y axis log scale", value = TRUE))),
                               br(),
                               column(12,
                                      tags$h6("Plot export settings"),
                                      column(3, textInput("DRplotExportWidth", "Width (inches): ", value = "12")),
                                      column(3, textInput("DRplotExportHeight", "Height (inches): ", value = "4.4"))),
                               actionButton(inputId = "return_to_DRTab", "Return to dose-response tab"),
                               HTML('<hr background-color="#000000" color = "#000000">')
                               ),
                           div(id = "TCglobalPlotSettings",
                               tags$h1("Time Course Plot Settings"),
                               column(12, selectInput("plotTypeTC", "Select plot type:", c("Line" = "Line", "Bar" = "Bar"))),
                               column(12,
                                      column(3, checkboxInput("XLogScaleTC", "Toggle X axis log scale", value = TRUE)),
                                      column(3, checkboxInput("YLogScaleTC", "Toggle Y axis log scale", value = TRUE))),
                               br(),
                               column(12,
                                      tags$h6("Plot export settings"),
                                      column(3, textInput("TCplotExportWidth", "Width (inches): ", value = "12")),
                                      column(3, textInput("TCplotExportHeight", "Height (inches): ", value = "4.4"))),
                               actionButton(inputId = 'return_to_TCTab', "Return to time-course tab"),
                               HTML('<hr background-color="#000000" color = "#000000">')
                               ),
                           div(id = "CRglobalPlotSettings",
                               tags$h1("Circadian Data Plot Settings"),
                               column(12,
                                      tags$h6("Top plot (avg) export size"),
                                      column(3, textInput("CRTopplotExportWidth", "Width: ", value = "12")),
                                      column(3, textInput("CRTopplotExportHeight", "Height: ", value = "4"))),
                               column(12,
                                      tags$h6("Bottom plot (counts) export size"),
                                      column(3, textInput("CRBottomplotExportWidth", "Width: ", value = "12")),
                                      column(3, textInput("CRBottomplotExportHeight", "Height: ", value = "4"))),
                               actionButton(inputId = "return_to_CRTab", "Return to circadian tab"),
                               HTML('<hr background-color="#000000" color = "#000000">')
                               ),
                           div(id = "BEglobalPlotSettings",
                               tags$h1("Basal Expression Plot Settings"),
                               column(12, checkboxInput("YLogScaleBE", "Toggle Y axis log scale", value = TRUE)),
                               column(12,
                                      tags$h6("Plot export settings"),
                                      column(3, textInput("BEplotExportWidth", "Width: ", value = "9")),
                                      column(3, textInput("BEplotExportHeight", "Height: ", value = "6"))),
                               actionButton(inputId = "return_to_BETab", "Return to basal expression tab"),
                               HTML('<hr background-color="#000000" color = "#000000">')
                               ),
                           div(id = "USglobalPlotSettings",
                               tags$h1("UpSet Plot Settings"),
                               column(12,
                                      tags$h6("Plot export settings"),
                                      column(3, textInput("plotExportWidth", "Width: ", value = "1920")),
                                      column(3, textInput("plotExportHeight", "Height: ", value = "1080"))),
                               actionButton(inputId = "return_to_USTab", "Return to UpSet tab")
                           )
                    )
                  #' Tableau viz import panel
                  #' Modify url in tableauImport.js in order to embed custom viz's
                  #' MUST open in browser in order to get viz to display
                  #' Currently not in use
                  #tabPanel(title = "Tableau Viz Testing", id = "vizContainer", value = "vizContainer",
                  #  br(),
                  #  column(12, wellPanel
                  #    (
                  #    useShinyjs(),
                  #    tags$head(tags$script(src="https://public.tableau.com/javascripts/api/tableau-2.js")),
                  #    extendShinyjs(script = "./www/tableauImport.js", functions = c("init")), 
                  #    tags$div(id = 'vizContainer')
                  #    ))
                  #  )
                  ), br(), br(), br()
                )
              )
            ),
            tabPanel("Gene List Query",
              title = "Gene List Query",
              fluidRow(
                column(2,
                  wellPanel(id = "GLD_left_panel", style = "background: #919191",
                    div(id = "GLD_left_panel_title", 
                        tags$h2("Upload Gene List", style = "color: white"),
                        HTML('<hr background-color="#333" color = "#333">')
                    ),
                    div(id = "GLD_left_panel_enter_list",
                      tags$h4("Enter Gene List:", style = "color: white"),
                      div(id = "GLD_paste_list",
                          textAreaInput(inputId = "GLD_paste_list_input", label = "A. Paste list", 
                                        placeholder = "Cyp1a1 \nCyp1a2 \n...", rows = 10, resize = "vertical"),
                          radioButtons(inputId = "GLD_paste_list_radiobutton", label = NULL, 
                                       choices = c("Symbol", "Ensembl ID"), selected = "Symbol", inline = TRUE)
                      ),
                      tags$h6("or", style = "color: white"),
                      div(id = "GLD_upload_file",
                        fileInput(inputId = "GLD_upload_file_input", label = "B. Choose from a file", 
                                  multiple = FALSE, buttonLabel = "Choose File"),
                        actionButton(inputId = "GLD_upload_preset_up", label = "Load Preset List - Upregulated"),
                        actionButton(inputId = "GLD_upload_preset_down", label = "Load Preset List - Downregulated"),
                        actionButton(inputId = "GLD_reset_file_input", label = "Clear file selection")
                      ),
                      actionButton(inputId = "GLD_submit_list_button", label = "Submit list")
                    ),
                    div(id = "save/load GLD bookmarks",
                      HTML('<hr background-color="#333" color = "#333">'),
                      actionButton("saveBookmarkGLD", "Save GLD Settings"),
                      actionButton("loadBookmarkGLD", "Load GLD Settings"),
                      HTML('<hr background-color="#333" color = "#333">')
                    ),
                    div(id = "resetButton",
                      #useShinyjs(),                                           # Include shinyjs in the UI
                      #extendShinyjs(text = jsResetCode),                      # Add the js code to the page
                      actionButton("reset_button2", "Reset Page")
                    )
                  )
                ),
                column(10,
                  tabsetPanel(
                    id = 'GLDTabs', type = 'tabs',
                    tabPanel(title = 'heatmap', id = 'heatmap', value = 'heatmap',
                      column(9,
                        selectInput("heatmapsSelected", "Select heatmaps to show:", multiple = TRUE,
                                    choices = projects[,2],
                                    selected = projects[,2][1],
                                    width = '100%'),
                        lapply(1:length(projects[,1]), function(i) {
                          div(id = paste0('heatmap', projects[,2][i]),
                            tags$h3(projects[,2][i]),
                            withSpinner(plotOutput(outputId = paste0('hm', projects[,2][i])))
                          )
                        })
                      ),
                      column(3,
                        wellPanel(id = "heatmap_metadata", style = "background: #919191",
                          div(id = "heatmap_metadata_title",
                            tags$h2("Metadata Filters", style = "color: white"),
                            HTML('<hr background-color="#333" color = "#333">')
                          ),
                          div(id = "heatmap_metadata_select",
                            selectInput("heatmap_metadata_selected", "Select metadata options to show:", multiple = TRUE,
                                        choices = metadata, selected = c("Species"))  
                          ),
                          div(id = "hmSpeciesFilter", checkboxGroupInput("heatmap_species_input", "Species:", species, species), br()),
                          div(id = "hmSexFilter", checkboxGroupInput("heatmap_sex_input", "Sex:", sex, sex), br()),
                          div(id = "hmOrganFilter", checkboxGroupInput("heatmap_organ_input", "Organ:", unique(organs$Organ_name), unique(organs$Organ_name)), br()),
                          div(id = "hmTreatmentFilter", checkboxGroupInput("heatmap_treatment_input", "Treatment:", treatments, treatments), br()),
                          div(id = "hmAssayFilter", checkboxGroupInput("heatmap_assay_input", "Assay:", assays, assays), br()),
                          div(id = "hmStrainFilter", checkboxGroupInput("heatmap_strain_input", "Strain:", strains, strains))
                        )
                      )
                    ),
                    tabPanel(title = 'GSEA', id = 'GSEA', value = 'GSEA',
                      column(10,
                        tags$h5("Note: GSEA analysis may take a few minutes."),
                        withSpinner(uiOutput("plots"))
                      ),
                      column(2,
                        wellPanel(id = "GSEA_metadata", style = "background: #919191",
                          div(id = "GSEA_metadata_title",
                              tags$h2("Metadata Filters", style = "color:white"),
                              HTML('<hr background-color="#333" color = "#333">')
                          ),
                          div(id = "GSEA_metadata_select",
                              selectInput("GSEA_metadata_selected", "Select metadata to group by:", multiple = FALSE,
                                          choices = c("Dose", metadata), selected = c("Dose"))
                          )
                        )
                      )
                    )
                  )
                )
              )
            ),
            tabPanel("Batch Analysis",
              tabsetPanel(id = "batch_analysis_tabs", type = "tabs",
              tabPanel("PCA",
                title = "PCA",
                fluidRow(
                  column(2,
                    wellPanel(id = "PCA_left_panel", style = "background: #919191",
                      div(id = "data_to_view_title",
                        tags$h4("Select data to view", style = "color: white")
                      ),
                      div(id = "PCA_num_genes_select",
                          numericInput(inputId = "pca_num_genes_to_analyze", label = "Select number of genes most variable to analyze", value = 500, min = 10, step = 10),
                      ),
                      div(id = "PCA_select_all_genes",
                          checkboxInput(inputId = "pca_select_all_genes_input", label = "Include all genes", value = FALSE)  
                      ),
                      actionButton(inputId = "pca_update_dataset", label = "Update Dataset"),
                      div(id = "pca_dataset_up_to_date",
                          "Dataset up to date."
                      ),
                      div(id = "pca_dataset_not_up_to_date",
                          "Dataset not up to date.
                          "
                      ),
                      selectInput("filterPCASelected", "Show/hide filters:", 
                                  choices = c("species", "strain", "sex", "organ", "chemical", "dose"), multiple = TRUE, 
                                  selected = c("species")),
                      div(id = "PCA_species_select",
                          radioButtons(inputId = "pca_species_input", label = "Select Species Data:", choices = unique(pca_metadata$SPECIES), selected = "MMU")    
                      ),
                      div(id = "PCA_strain_select",
                          checkboxGroupInput(inputId = "pca_strain_input", label = "Select Strain Data:", choices = unique(pca_metadata$STRAIN), selected = unique(pca_metadata$STRAIN))  
                      ),
                      div(id = "PCA_sex_select",
                          checkboxGroupInput(inputId = "pca_sex_input", label = "Select Sex Data:", choices = c("M", "F"), selected = c("M", "F"))  
                      ),
                      div(id = "PCA_organ_select",
                          checkboxGroupInput(inputId = "pca_organ_input", label = "Select Organ Data:", choices = unique(pca_metadata$TISSUE), selected = unique(pca_metadata$TISSUE))
                      ),
                      div(id = "PCA_chemical_select",
                          checkboxGroupInput(inputId = "pca_chemical_input", label = "Select Treatment Data:", choices = unique(pca_metadata$CHEMICAL), selected = unique(pca_metadata$CHEMICAL))  
                      ),
                      div(id = "PCA_dose_select",
                          checkboxGroupInput(inputId = "pca_dose_input", label = "Select Dose Data:", choices = unique(pca_metadata$DOSE), selected = unique(pca_metadata$DOSE))  
                      )
                    )
                  ),
                  column(8,
                    tabsetPanel(id = 'PCA_tabs', type = 'tabs',
                      tabPanel(title = "PCA Individuals Plot", id = "pca_ind_plot", value = "pca_ind_plot",
                        withSpinner(plotOutput(outputId = "PCA_plot_inds"))
                      ),
                      tabPanel(title = "PCA Variables Plot", id = "pca_var_plot", value = "pca_var_plot",
                        withSpinner(plotOutput(outputId = "PCA_plot_vars"))
                      ),
                      tabPanel(title = "Scree plot", id = "scree_plot", value = "scree_plot",
                        withSpinner(plotOutput(outputId = "Scree_plot"))
                      )
                    ),
                    withSpinner(tableOutput(outputId = "PCA_variable_genes_table"))
                  ),
                  column(2,
                    wellPanel(id = "PCA_plot_settings", style = "background: #919191",
                      div(id = "pca_plot_settings_title",
                        tags$h4("Plot Settings", style = "color: white")
                      ),
                      div(id = "pca_dimension_select",
                        "Select dimensions to plot:",
                        numericInput(inputId = "pca_dim1_input", label = NULL, value = 1, min = 1, max = 5, step = 1),
                        numericInput(inputId = "pca_dim2_input", label = NULL, value = 2, min = 1, max = 5, step = 1)
                      ),
                      div(id = "pca_top_genes_select",
                        sliderInput(inputId = "pca_top_genes_to_show", label = "Select number of most variable genes to show:", min = 5, max = 20, value = 5, step = 1)
                      ),
                      div(id = "pca_ellipse_group_by",
                        selectInput(inputId = "pca_ellipse_group_by_input", label = "Group ellipses by:", choices = c("SPECIES", "STRAIN", "SEX", "TISSUE", "CHEMICAL", "DOSE"),
                                    selected = "STRAIN", multiple = TRUE)
                      ),
                      div(id = "pca_show_hide_labels",
                        checkboxInput(inputId = "pca_show_hide_label_input", label = "Show plot labels", value = FALSE)
                      ),
                      div(id = "pca_show_hide_ellipses",
                          checkboxInput(inputId = "pca_show_hide_ellipses_input", label = "Show ellipses", value = TRUE)
                      ),
                      actionButton(inputId = 'pca_plot_button', label = "Plot"),
                      div(id = "pca_out_of_sync_warning",
                        "WARNING: Dataset stored and plot displayed are out of sync - press plot button above to update."
                      ),
                      div(id = "pca_no_rds_file_warning",
                          "WARNING: No Rdata file in memory available for plotting. Please update dataset in left panel and try again."
                      )
                    )
                  )
                )
              ),
              tabPanel("PLS",
                title = "PLS",
                fluidRow(
                  column(2,
                    wellPanel(id = "PLS_left_panel", style = "background: #919191",
                      div(id = "data_to_view_title",
                        tags$h4("Select data to view", style = "color: white")
                      ),
                      div(id = "PLS_num_genes_select",
                          numericInput(inputId = "pls_num_genes_to_analyze", label = "Select number of genes most variable to analyze", value = 500, min = 10, step = 10),
                      ),
                      div(id = "PLS_select_all_genes",
                          checkboxInput(inputId = "pls_select_all_genes_input", label = "Include all genes", value = FALSE)  
                      ),
                      actionButton(inputId = "pls_update_dataset", label = "Update Dataset"),
                      "Update dataset only required when changing number of most variable genes to display.",
                      selectInput("filterPLSSelected", "Show/hide filters:", 
                                  choices = c("species", "strain", "sex", "organ", "chemical", "dose"), multiple = TRUE, 
                                  selected = c("species")),
                      div(id = "PLS_species_select",
                          radioButtons(inputId = "pls_species_input", label = "Select Species Data:", choices = unique(pls_metadata$SPECIES), selected = "MMU")    
                      ),
                      div(id = "PLS_strain_select",
                          checkboxGroupInput(inputId = "pls_strain_input", label = "Select Strain Data:", choices = unique(pls_metadata$STRAIN), selected = unique(pls_metadata$STRAIN))
                      ),
                      div(id = "PLS_sex_select",
                          checkboxGroupInput(inputId = "pls_sex_input", label = "Select Sex Data:", choices = c("M", "F"), selected = c("M", "F"))
                      ),
                      div(id = "PLS_organ_select",
                          checkboxGroupInput(inputId = "pls_organ_input", label = "Select Organ Data:", choices = unique(pls_metadata$TISSUE), selected = unique(pls_metadata$TISSUE))
                      ),
                      div(id = "PLS_chemical_select",
                          checkboxGroupInput(inputId = "pls_chemical_input", label = "Select Treatment Data:", choices = unique(pls_metadata$CHEMICAL), selected = unique(pls_metadata$CHEMICAL))
                      ),
                      div(id = "PLS_dose_select",
                          checkboxGroupInput(inputId = "pls_dose_input", label = "Select Dose Data:", choices = unique(pls_metadata$DOSE), selected = unique(pls_metadata$DOSE))
                      )
                    )
                  ),
                  column(8,
                    withSpinner(plotOutput('PLS_indivs')),
                    withSpinner(plotOutput('PLS_top_genes'))
                  ),
                  column(2,
                    wellPanel(id = "PLS_plot_settings", style = "background: #919191",
                      div(id = "PLS_plot_settings_title",
                        tags$h4("Plot Settings", style = "color: white")
                      ),
                      div(id = "PLS_ellipse_group_by",
                        selectInput(inputId = "pls_ellipse_group_by_input", label = "Group ellipses by:", choices = c("SPECIES", "STRAIN", "SEX", "TISSUE", "CHEMICAL", "DOSE"),
                                    selected = "DOSE", multiple = FALSE)
                      ),
                      div(id = "PLS_labels",
                        selectInput(inputId = "pls_point_labels", label = "Label points by:", choices = c("SPECIES", "STRAIN", "SEX", "TISSUE", "CHEMICAL", "DOSE"),
                                    selected = "TISSUE", multiple = FALSE)
                      ),
                      actionButton(inputId = 'pls_plot_button', label = "Plot"),
                      div(id = "pls_out_of_sync_warning",
                          "WARNING: Dataset stored and plot displayed are out of sync - press plot button above to update."
                      ),
                      div(id = "pls_no_rds_file_warning",
                        "WARNING: No Rdata file in memory available for plotting. Please update dataset in left panel and try again."
                      )
                    )
                  )
                )
              )
            )),
            tabPanel("Single Cell",
              title = "Single Cell Data",
              fluidRow(
              column(2,
                wellPanel(id = "singlecell_left_panel", style = "background: #919191",
                    div(id = "singlecell_select_dataset",
                      selectInput(inputId = "singlecell_dataset_input", label = "Select Dataset:", choices = c("Large", "Small"), selected = "Small", multiple = FALSE)
                    ),
                    div(id = "singlecell_genelist_left_panel",
                      tags$h4("Enter Gene List:", style = "color: white"),
                      div(id = "singlecell_paste_list",
                        textAreaInput(inputId = "singlecell_paste_list_input", label = "A. Paste list",
                                      placeholder = "Cyp1a1 \nCyp1a2 \n...", rows = 10, resize = "vertical"),
                        radioButtons(inputId = "singlecell_paste_list_radiobutton", label = NULL,
                                     choices = c("Symbol", "Ensembl ID"), selected = "Symbol", inline = TRUE)
                      ),
                      tags$h6("or", style = "color: black"),
                      div(id = "singlecell_upload_file",
                        fileInput(inputId = "singlecell_upload_file_input", label = "B. Choose from a file",
                                  multiple = FALSE, buttonLabel = "Choose File"),
                        HTML('<hr background-color="#000000" color = "#000000">'),
                        actionButton(inputId = "singlecell_upload_preset_up", label = "Load Preset List - Upregulated"),
                        actionButton(inputId = "singlecell_upload_preset_down", label = "Load Preset List - Downregulated"),
                        actionButton(inputId = "singlecell_reset_file_input", label = "Clear file selection"),
                        HTML('<hr background-color="#000000" color = "#000000">')
                      )
                    ),
                    div(id = "singlecell_plot_button",
                        br(),
                        actionButton(inputId = "plot_singlecell", label = "Plot")
                    )
                )
              ),
              column(8,
                tabsetPanel(id = "singlecell_tabs", type = "tabs",
                tabPanel("UMAP",
                  withSpinner(plotOutput("UMAP")), br(), br(), br(), br(), br(), br(), br(), br(), br()
                ),
                tabPanel("Feature Plot",
                  withSpinner(plotOutput("FeaturePlot"))
                ),
                tabPanel("Ridge Plot",
                  withSpinner(plotOutput("RidgePlot"))
                ),
                tabPanel("Violin Plot",
                  withSpinner(plotOutput("ViolinPlot"))
                ),
                tabPanel("Dot Plot",
                  withSpinner(plotlyOutput("DotPlot"))
                )
                )
              ),
              column(2,
                wellPanel(id = "singlecell_right_panel",
                  selectInput(inputId = "singlecell_metadata_select", label = "Select metadata to group by:", choices = c("ident", "treatment", "celltype"),
                              selected = "celltype", multiple = FALSE)
                )
              )
              )
            )
        ),
            # Footer section
            div(id = "footer",
              includeHTML("./www/footer.html"),
              tableOutput('foo'), # Monitor global memory usage
            )
        )
    )
}
