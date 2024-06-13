fluidRow(
  shinyjs::hidden(p(id = "text_gsea_plot")),
  shinyjs::hidden(
    div(id = "choose_GO_box",
        box(
          title = "Running Score Term Selection", status = "success", solidHeader = TRUE, collapsible = TRUE,
          selectInput("choose_GO", "Choose:", "", selected = NULL)
        ))),
  plotOutput("gseaplot", height = "600px")
)