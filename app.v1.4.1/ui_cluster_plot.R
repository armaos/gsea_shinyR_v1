fluidRow(
  div(id = "cluster_heatmap",
      box(
        width = 8,
        title = "Configuration", status = "primary", solidHeader = TRUE,
        sliderInput("number_of_clusters", label = p("How many clusters do you identify?"), 
                    min = 1, step = 1, max = 20, value = 5),
        actionButton("submit_number_of_clusters", label = "Annotate!")),
      plotOutput("all_cluster_heatmap", height = "700px"))
)