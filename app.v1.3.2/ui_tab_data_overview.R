fluidRow(
    shinyjs::hidden(div(id = "data_overview_box",
                        
                        box(pickerInput("overview_table_Default_values", "Contrast values:", choices = c("Name"), multiple = TRUE),
                            pickerInput("overview_table_terms", "Annotations:", choices = c(), multiple = TRUE)
                        ))),
    shinyjs::hidden(div(
    id = "data_overview_table_box",
    box(title = "Data", width = 12, solidHeader = TRUE, status = "primary", 
        collapsible = TRUE,
        DT::DTOutput('data_overview_table'))))
    #DT::DTOutput('data_overview_table')
)