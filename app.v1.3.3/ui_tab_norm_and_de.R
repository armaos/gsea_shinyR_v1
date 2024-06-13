fluidRow(
  div(id = "menu_norm_and_de_panel",
      tabsetPanel(type = "tabs",
                  tabPanel("PCA plot",
                           plotOutput("pca_plot", height = "600px")),
                  tabPanel("Volcano plot",
                           plotOutput("volcano0", height = "600px")),
                  tabPanel("Heat-Map plot",
                           plotOutput("pheatmap", height = "600px"))
      ))
)