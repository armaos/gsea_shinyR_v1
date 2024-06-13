library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)

ui <- dashboardPage(
  dashboardHeader(title = "K.A. RShiny v1.2.2"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data", 
               menuSubItem("Import", tabName = "menu_data_import"),
               menuSubItem("Overview", tabName = "data_overview"),
               icon = icon("th")),
      menuItem("Normalization and DE", tabName = "menu_norm_and_de", icon = icon("th")),
      menuItem("ORA",
               menuSubItem("Configuration And Table", tabName = "configure_ora"),
               menuSubItem("Plots", tabName = "plots_ora"),
               menuSubItem("KEGGscape", tabName = "KEGGscape_ora"),
               icon = icon("th")),
      menuItem("GSEA",
               menuSubItem("Configuration And Table", tabName = "table_gsea"),
               menuSubItem("Plots", tabName = "plots_gsea"),
               menuSubItem("Running Score", tabName = "running_score_gsea"),
               menuSubItem("Custom UpSet", tabName = "custom_upset_gsea"),
               menuSubItem("KEGGscape", tabName = "KEGGscape_gsea"),
               icon = icon("th"))
    )
  ),
  # dashboardBody(
  #   # Boxes need to be put in a row (or column)
  #   tabItems(
  #     # First tab content
  #     tabItem(tabName = "menu_data_import",
  #             h2("Widgets tab content")
  # )) )) 
  
  dashboardBody(
    useShinyjs(),
    # Boxes need to be put in a row (or column)
    tabItems(
      # First tab content
      tabItem(tabName = "data_overview",
              shinyjs::hidden(div(id = "data_overview_box",
                                  
                  box(pickerInput("overview_table_Default_values", "Contrast values:", choices = c("Name"), multiple = TRUE),
                      pickerInput("overview_table_terms", "Annotations:", choices = c(), multiple = TRUE)
                      ))),
              DT::DTOutput('data_overview')),
      tabItem(tabName = "menu_data_import",
              fluidRow(
                box(
                  title = "Select Data source", status = "primary", solidHeader = TRUE,
                  
                  selectInput("choose", "Choose:", choices = c("Pre-compiled", "Custom"), selected = NULL),
                  #actionButton("demo", "Fill in fields with demo files"),
                  conditionalPanel(
                    "input.choose=='Custom'",
                    fileInput("data_file", "Read Counts:", multiple = TRUE, accept = c("text/csv","text/comma-separated-values,text/plain")),
                    actionButton("submit_norm_and_de", label = "Normalise and Combine!")
                    
                  ),
                  conditionalPanel(
                    "input.choose=='Pre-compiled'",
                    h3("Select pre calculated condition:"),
                    radioButtons("Condition", "Choose one:",
                                 c("D_O" = "D_O",
                                   "L_O" = "L_O",
                                   "D_L" = "D_L",
                                   "B_G_in_O" = "B_G_in_O",
                                   "B_G_in_L" = "B_G_in_L",
                                   "B_G_in_D" = "B_G_in_D"
                                 )),
                    actionButton("submit_contrast", label = "Submit Contrast"),
                    shinyjs::hidden(h2(id = "submitted_contrast1", "Contrast submitted Successfully.")),
                    shinyjs::hidden(p(id = "submitted_contrast2", "You can head to ORA and GSEA Analysis now!")),
                    
                  ),
                  
                  
                )
              )
      ),
      
      # Second tab content
      tabItem(tabName = "menu_norm_and_de",
              h2("Here to be some plots and tables for the normalisation/quantiication and contrast analysis"),
              shinyjs::hidden(
                box(
                  id = "menu_norm_and_de_box",
                  title = "DEG table", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                  DT::DTOutput('DEG_table')
                )),
      ),
      tabItem(tabName = "configure_ora",
              box(
                width = 4,
                title = "OverRepresentation Analysis options", status = "warning", solidHeader = TRUE, collapsible = TRUE,
                sliderInput("slider_lfc", label = p("|logFC| Range"), min = 0, step = 0.5,
                            max = 12, value = c(4, 12)),
                
                sliderTextInput("slider_fdr", label = p("FDR treshold"),
                                choices = seq(from = 100, to = 0, by = -1),
                                selected = 10, pre = "1e-"),
                verticalLayout(
                  actionButton("run", label = "Submit Thresholds"))),
              shinyjs::hidden(p(id = "text_configure_ora", "You haven't submitted any contrast yet!")),
              shinyjs::hidden(
                div(id = "box_configure_term_ora",
                    box(
                      width = 4,
                      title = "Term Enrichment", status = "danger", solidHeader = TRUE, collapsible = TRUE,
                      radioButtons("enrichment_term", label = "Select enrichment",
                                   choiceNames = c("GO:Biological Processes", "GO:Molecular Functions", "GO:Cellular Components", "Reactome", "KEGG Pathways"),
                                   choiceValues = c("biological_process", "molecular_function", "cellular_component", "Reactome", "KEGG_Pathways")
                      ),
                      actionButton("submit_enrichment", label = "Submit enrichment")
                    ))),
              shinyjs::hidden(
                div(id = "box_ora_send_to_plot",
                    box(
                      width = 4,
                      title = "Send to plot", status = "success", solidHeader = TRUE, collapsible = TRUE,
                      sliderInput("ORA_max_elements",
                                  label = p("Maximum terms per plot to show"), min = 1,
                                  max = 50, value = 20),
                      actionButton("submit_plot_ora", label = "Plot table")
                    ))),
              shinyjs::hidden(p(id = "text_ora_go_table")),
              DT::DTOutput('ora_go_table')
      ),
      tabItem(tabName = "table_ora",
              p(id = "text_ora_table", "Please submit your configuration first!"),
              
              
      ),
      
      tabItem(tabName = "plots_ora",
              p(id = "text_ora_plots", "Please submit your configuration first and plot your table!"),
              shinyjs::hidden(
                div(id = "ora_plots_tabsetPanel",
                    tabsetPanel(type = "tabs",
                                tabPanel("Bar plot - Gene Counts",
                                         p(id = "ora_barplot1_text", "Processing..."),
                                         plotOutput("barplot1", height = "600px")),
                                tabPanel("Bar plot - qscore",
                                         p(id = "ora_barplot2_text", "Processing..."),                            
                                         plotOutput("barplot2", height = "600px")),
                                tabPanel("Dot Plot",
                                         p(id = "ora_dotplot_text", "Processing..."),                            
                                         plotOutput("dotplot", height = "600px")),
                                tabPanel("Enrichment Map Plot",
                                         p(id = "ora_emapplot_text", "Processing..."),                            
                                         plotOutput("emapplot", height = "600px")),
                                tabPanel("Gene-Concept Network Plot",
                                         p(id = "ora_cnetplot_text", "Processing..."),                            
                                         plotOutput("cnetplot", height = "600px")),
                                tabPanel("HeatMap Plot",
                                         p(id = "ora_heatplot_text", "Processing..."),                            
                                         plotOutput("heatplot", height = "600px")),
                                tabPanel("Tree Plot",
                                         p(id = "ora_treeplot_text", "Processing..."),                            
                                         plotOutput("treeplot", height = "600px")),
                                tabPanel("UpSetPlot",
                                         p(id = "ora_upsetplot_text", "Processing..."),                            
                                         plotOutput("upsetplot")),
                    )))
              
              
      ),
      
      tabItem(tabName = "KEGGscape_ora",
              navbarPage("KEGGscape Visualization", 
                         tabPanel(
                           title = tagList(icon("table"), "Data"),                           
                           shinyjs::hidden(
                             div(id = "ora_kegg_pathway_dropdown",
                                 box(
                                   title = "KEGG pathways", status = "success", solidHeader = TRUE, collapsible = TRUE,
                                   selectInput("choose_map_ora", "Choose:", "",selected = NULL)))),
                           p(id = "text_KEGGscape_ORA",  "Please submit your configuration first!"),        
                           shinyjs::hidden(
                             div(id = "ora_kegg_pathway_table",
                                 box(
                                   title = "KEGG pathway elements", status = "success", solidHeader = TRUE, collapsible = TRUE,
                                   width = 12,
                                   DT::DTOutput('ora_keggscape_table'))))
                           ),
                         tabPanel(
                            title = tagList("Plot"),
                            shinyjs::hidden(
                             div(id = "keggscape_ora_conf",
                                 box( width = 4,
                                    title = "Configure", status = "success", solidHeader = TRUE, collapsible = TRUE,
                                    actionButton("submit_kegg_ora_toggle", label = "Toggle colours")
                                     #verbatimTextOutput("verb")
                                    ))),
                            plotOutput("KEGG_map", height = "600px"))
                         )
              ),
      tabItem(tabName = "table_gsea",
              h2("Gene Set Enrichment Analysis"),
              
              div(id = "run_gsea_box",
                  box(
                    title = "Term Enrichment", status = "success", solidHeader = TRUE, collapsible = TRUE,
                    radioButtons("enrichment_term_gsea", label = "Select enrichment",
                                 choiceNames = c("GO:Biological Processes", "GO:Molecular Functions", "GO:Cellular Components", "KEGG Pathways"),
                                 choiceValues = c("biological_process", "molecular_function", "cellular_component", "KEGG_Pathways")),
                    actionButton("run_gsea", label = "Run GSEA!"))),
              shinyjs::hidden(p(id = "text_configure_gsea", "You haven't submitted any contrast yet!")),
              shinyjs::hidden(p(id = "text_gsea_table")),
              DT::DTOutput('gsea_table')
      ),
      tabItem(tabName = "plots_gsea",
              shinyjs::hidden(
                div(id = "gsea_plots_tabsetPanel",
                    tabsetPanel(type = "tabs",
                                tabPanel("Dot Plot",
                                         p(id = "gsea_dotplot_text", "Processing..."),                                                      
                                         plotOutput("dotplot_gsea", height = "600px") ),
                                tabPanel("Ridge Plot",
                                         p(id = "gsea_ridgeplot_text", "Processing..."),                                                                                
                                         plotOutput("ridgeplot_gsea", height = "600px")))
                ))),
      
      tabItem(tabName = "running_score_gsea",
              shinyjs::hidden(p(id = "text_gsea_plot")),
              shinyjs::hidden(
                div(id = "choose_GO_box",
                    box(
                      title = "Running Score Term Selection", status = "success", solidHeader = TRUE, collapsible = TRUE,
                      selectInput("choose_GO", "Choose:", "", selected = NULL)
                    ))),
              plotOutput("gseaplot", height = "600px")
      ),
      tabItem(tabName = "custom_upset_gsea",
              shinyjs::hidden(
                div(id = "custom_upset_box",
                    box(
                      title = "Custom UpSet Plot", status = "success", solidHeader = TRUE, collapsible = TRUE,
                      pickerInput("custom_upset_terms", "Select terms:", choices = c(), multiple = TRUE),
                      actionButton("submit_custom_upset_plot", "Plot!")
                    ))),
              
              plotOutput("gsea_custom_upsetplot", height = "600px")
      ),
      
      tabItem(tabName = "KEGGscape_gsea",
              navbarPage("KEGGscape Visualization", 
                         tabPanel(
                           title = tagList(icon("table"), "Data"),
                           shinyjs::hidden(p(id = "text_KEGGscape_GSEA")),
                           shinyjs::hidden(
                             div(id = "gsea_kegg_pathway_dropdown",
                                 box(
                                   title = "KEGG pathways", status = "success", solidHeader = TRUE, collapsible = TRUE,
                                   selectInput("choose_map_gsea", "Choose:", "",selected = NULL)
                                 ))),
                           shinyjs::hidden(
                             div(id = "gsea_kegg_pathway_table",
                                 box(
                                   title = "KEGG pathway elements", status = "success", solidHeader = TRUE, collapsible = TRUE,
                                   width = 12,
                                   DT::DTOutput('gsea_keggscape_table')
                                 )))
                         ),
                         tabPanel(
                           title = tagList("Plot"),
                           shinyjs::hidden(
                             div(id = "keggscape_gsea_conf",
                                 box( width = 4,
                                    title = "Configure", status = "success", solidHeader = TRUE, collapsible = TRUE,
                                    actionButton("submit_kegg_gsea_toggle", label = "Toggle colours")
                                     #verbatimTextOutput("verb")
                                    ))),
                           plotOutput("KEGG_map_gsea", height = "600px"))
              )
      )
    )
    
    
  )
)