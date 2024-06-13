fluidRow(
  p(id = "text_ora_plots_dynamic", "Please submit your configuration first and plot your table!"),
  shinyjs::hidden(
    div(id = "ora_plots_dynamic_tabsetPanel",
        tabsetPanel(type = "tabs", id="tabs_plots_ora_dynamic"
        ))
    )
)