fluidRow(
  shinyjs::hidden(
    div(id = "custom_upset_box",
        box(
          title = "Custom UpSet Plot", status = "success", solidHeader = TRUE, collapsible = TRUE,
          pickerInput("custom_upset_terms", "Select terms:", choices = c(), multiple = TRUE),
          actionButton("submit_custom_upset_plot", "Plot!")
        ))),

  plotOutput("gsea_custom_upsetplot", height = "600px"), type = 6
)