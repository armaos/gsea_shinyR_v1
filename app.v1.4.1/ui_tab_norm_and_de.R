fluidRow(
  div(id = "menu_norm_and_de_panel",
      # tabsetPanel(type = "tabs", id="tabs_plots_norm_de_dynamic",
      #             tabPanel("PCA plot",
      #                      plotOutput("pca_plot", height = "600px")),
      #                      ))
      tabsetPanel(type = "tabs",
                  tabPanel("PCA plot",
                           plotOutput("pca_plot", height = "600px")
                           ),
                  # tabPanel("Cluster plot",
                  #          source(file = file.path(version, "ui_cluster_plot.R"),
                  #                 local = TRUE,
                  #                 encoding = "UTF-8")
                  # ),
                  tabPanel("Volcano plot",
                           tabsetPanel(type = "tabs", id="volcano_norm_de_dynamic")
                           #plotOutput("volcano0", height = "600px")
                           ),
                  tabPanel("Heat-Map plot",
                           tabsetPanel(type = "tabs", id="heatmap_norm_de_dynamic")
                           #plotOutput("pheatmap", height = "600px")
                           )
      ))
)