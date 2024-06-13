fluidRow(
    shinyjs::hidden(div(id = "selected_data_overview_box",
                        
                        box(
                            title = "Column Selection", width = 12, solidHeader = TRUE, status = "primary", collapsible = TRUE,
                            h4("Primary columns"),
                            fluidRow(
                                column(width = 4,pickerInput("selected_overview_table_Default_values", "Contrast values:", choices = c("Name"), multiple = TRUE))
                            ),
                            linebreaks(1),
                            h4("Annotations per Software"),
                            fluidRow(
                                column(width = 4, pickerInput("selected_annotatations_ips", "InterProScan:", choices = c(), multiple = TRUE)),
                                column(width = 4, pickerInput("selected_annotatations_em", "eggNOG Mapper:", choices = c(), multiple = TRUE)),
                                column(width = 4, pickerInput("selected_annotatations_kfk", "KFK Koala:", choices = c(), multiple = TRUE))
                            ),
                            fluidRow(    
                                column(width = 4, pickerInput("selected_annotatations_clusters", "Clusters:", choices = c(), multiple = TRUE)),
                                column(width = 4, pickerInput("selected_annotatations_b2g", "Blast2GO:", choices = c(), multiple = TRUE)),
                                column(width = 4, pickerInput("selected_annotatations_dbCAN2", "dbCAN2:", choices = c(), multiple = TRUE))
                            ),
                            linebreaks(1),
                            h4("Rest of Annotations"),
                            fluidRow(
                                column(width = 4, pickerInput("selected_annotatations_others", "Other:", choices = c(), multiple = TRUE))
                            ),
                            linebreaks(1),
                            downloadButton("download_selected_Data", "Download")
                        ))),
    shinyjs::hidden(div(
        id = "selected_data_overview_table_box",
        box(title = "Data", width = 12, solidHeader = TRUE, status = "primary", 
            collapsible = TRUE,
            DT::DTOutput('selected_data_overview_table')))
    )
    #DT::DTOutput('data_overview_table')
)