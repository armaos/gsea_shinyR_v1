verticalLayout(
  shinyjs::hidden(p(id = "text_gsea_plot")),
  shinyjs::hidden(
    div(id = "choose_GO_box",
        box(
          title = "Running Score Term Selection", status = "success", solidHeader = TRUE, collapsible = TRUE, width = 8,
          selectInput("choose_GO", "Choose:", "", selected = NULL)
        )) 
    ),
  tabsetPanel(type = "tabs", id="tabs_running_score_dynamic")
  
)