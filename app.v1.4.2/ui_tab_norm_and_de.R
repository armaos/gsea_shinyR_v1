fluidRow(
  div(id = "menu_norm_and_de_panel",
      tabsetPanel(type = "tabs", id = "norm_de_plot_tabsetpanel",
                  tabPanel("PCA plot",
                           withSpinner(plotOutput("pca_plot", height = "600px"), type = 6)
                           ),
                  tabPanel("Heat-Map plot (overview)",
                           withSpinner(plotOutput("heatmap_overview_plot", height = "600px"), type = 6)
                  ),
                  # tabPanel("Cluster plot",
                  #          source(file = file.path(version, "ui_cluster_plot.R"),
                  #                 local = TRUE,
                  #                 encoding = "UTF-8")
                  # ),
                  tabPanel("Volcano plot",  
                           tabsetPanel(type = "tabs", id="volcano_norm_de_dynamic")
                           ),
                  tabPanel("Heat-Map plot (pairwise)", 
                           tabsetPanel(type = "tabs", id="heatmap_norm_de_dynamic")
                           )
      ))
)