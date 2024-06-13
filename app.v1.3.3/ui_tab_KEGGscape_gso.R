fluidRow(
  navbarPage("KEGGscape Visualization", 
            tabPanel(
              title = tagList(icon("table"), "Data"),                           
              shinyjs::hidden(
                div(id = "gso_kegg_pathway_dropdown",
                    box(
                      title = "KEGG pathways", status = "success", solidHeader = TRUE, collapsible = TRUE,
                      selectInput("choose_map_gso", "Choose:", "",selected = NULL)))),
              p(id = "text_KEGGscape_gso",  "Please submit your configuration first!"),        
              shinyjs::hidden(
                div(id = "gso_keggscape_table_box",
                    
                    box(
                      title = "KEGG pathway elements", status = "success", solidHeader = TRUE, collapsible = TRUE,
                      width = 12,
                      DT::DTOutput('gso_keggscape_table'))))
            ),
            tabPanel(
              title = tagList("Plot"),
              shinyjs::hidden(
                div(id = "keggscape_gso_conf",
                    box( width = 4,
                          title = "Configure", status = "success", solidHeader = TRUE, collapsible = TRUE,
                          actionButton("submit_kegg_gso_toggle", label = "Toggle colours"),
                          selectInput( "kegg_gso_values", "Choose overlaying data:", 
                                      choices = c("logFC" = "logfc", 
                                                  "KFK annotatation" = "kfk", 
                                                  "eggNOG annotation" = "em") , 
                                      selected = "logfc")
                    ))),
              div(style = "width: 100%; overflow-y: scroll;overflow-x: scroll;",
                  plotOutput("KEGG_map_gso", height = "600px")
              )
              #plotOutput("KEGG_map_gso", height = "600px")
            )
  )
)