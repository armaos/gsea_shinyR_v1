fluidRow(
  shinyjs::hidden(
    div(id = "custom_upset_box",
        box(
          width = 12,
          title = "Custom UpSet Plot", status = "success", solidHeader = TRUE, collapsible = TRUE,
          pickerInput("custom_upset_terms", "Select terms:", choices = c(), multiple = TRUE),
          actionButton("submit_custom_upset_plot", "Plot!")
        ))),

  shinyjs::hidden(
    div(id = "custom_upsetplot_box",
        box(
          title = "", status = "success", solidHeader = TRUE, collapsible = TRUE,
          width = 12,
          plotOutput("gsea_custom_upsetplot", height = "600px")
        ))),

  shinyjs::hidden(
    div(id = "custom_upsettable_box",
        box(
          title = "Custom Upset table per Gene", status = "success", solidHeader = TRUE, collapsible = TRUE,
          width = 12,
          downloadButton("download_upset_table_per_gene", "Download"),
          checkboxInput("add_annotations_in_custom_upset_table_per_gene", "Append Annotations and GE data to the table" , FALSE),
          linebreaks(2),
          DT::DTOutput('custom_upset_table_per_gene')          
        ),
        box(
          title = "Custom Upset table per Term Combination", status = "success", solidHeader = TRUE, collapsible = TRUE,
          width = 12,
          downloadButton("download_upset_table", "Download"),
          linebreaks(2),
          DT::DTOutput('custom_upset_table')          
        ),
        shinyjs::disabled(actionButton("generate_search_term", label = "Generate K.A Database Explorer search term")),
        shinyjs::disabled(actionButton("get_custom_upset_genes", label = "Get Gene List"))
    ))

  
  
  
)