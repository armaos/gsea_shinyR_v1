
#here the plots for GSEA
appendTab(inputId = "tabs_plots_gsea_dynamic",
          tabPanel(gsea_cond,
                   tabsetPanel(type = "tabs",
                               tabPanel("Dot Plot",
                                        p(id = paste0("gsea_dotplot_text_", gsea_cond), "Processing..."),
                                        plotOutput(paste0("dotplot_gsea_", gsea_cond), height = "600px") ),
                               tabPanel("Ridge Plot",
                                        p(id = paste0("gsea_ridgeplot_text_", gsea_cond), "Processing..."),
                                        plotOutput(paste0("ridgeplot_gsea_", gsea_cond), height = "600px"))
                   )
          )
)

# here for the running score
appendTab(inputId = "tabs_running_score_dynamic",
          tabPanel(gsea_cond,
                  # shinyjs::hidden(
                  #   div(id = paste0("choose_GO_box_",gsea_cond ),
                         box(
                           title = "Running Score Term Selection", status = "success", solidHeader = TRUE, collapsible = TRUE, width = 8,
                           selectInput(paste0("choose_GO_", gsea_cond ), "Choose:", "", selected = NULL)
                         ),
                  plotOutput(paste0("gseaplot_", gsea_cond), height = "600px")
                   # ) 
                   #)
          )
)



source(file = file.path(version, "server_GSEA_plots_dynamic.R"),
       local = TRUE,
       encoding = "UTF-8")