fluidRow(
    navbarPage("KEGGscape Visualization",
              tabPanel(
                title = tagList(icon("table"), "Data"),
                shinyjs::hidden(
                  div(id = "ora_kegg_pathway_dropdown",
                      box(
                        title = "KEGG pathways", status = "success", solidHeader = TRUE, collapsible = TRUE,
                        selectInput("choose_map_ora", "Choose:", "",selected = NULL)))),
                p(id = "text_KEGGscape_ORA",  "Please submit your configuration first!"),
                shinyjs::hidden(
                  div(id = "ora_keggscape_table_box",
                      box(
                        title = "KEGG pathway elements", status = "success", solidHeader = TRUE, collapsible = TRUE,
                        width = 12,
                        DT::DTOutput('ora_keggscape_table'))))
              ),
              tabPanel(
                title = tagList("Plot"),
                shinyjs::hidden(
                  div(id = "keggscape_ora_conf",
                      box( width = 4,
                          title = "Configure", status = "success", solidHeader = TRUE, collapsible = TRUE,
                          actionButton("submit_kegg_ora_toggle", label = "Toggle colours")
                          #verbatimTextOutput("verb")
                      ))),
                plotOutput("KEGG_map", height = "600px"))
  )
)