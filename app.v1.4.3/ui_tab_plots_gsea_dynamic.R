fluidRow(
  shinyjs::hidden(
    div(id = "gsea_plots_dynamic_tabsetPanel",
        tabsetPanel(type = "tabs", id="tabs_plots_gsea_dynamic",
        )
        # tabsetPanel(type = "tabs",
        #             tabPanel("Dot Plot",
        #                      p(id = "gsea_dotplot_text", "Processing..."),                                                      
        #                      plotOutput("dotplot_gsea", height = "600px") ),
        #             tabPanel("Ridge Plot",
        #                      p(id = "gsea_ridgeplot_text", "Processing..."),                                                                                
        #                      plotOutput("ridgeplot_gsea", height = "600px")))
    ))
)