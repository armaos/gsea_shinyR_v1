fluidRow(
    navbarPage("KEGGscape Visualization",
            tabPanel(
              title = tagList(icon("table"), "Data"),
              shinyjs::hidden(p(id = "text_KEGGscape_GSEA")),
              shinyjs::hidden(
                div(id = "gsea_kegg_pathway_dropdown",
                    box(
                      title = "KEGG pathways", status = "success", solidHeader = TRUE, collapsible = TRUE,
                      selectInput("choose_map_gsea", "Choose:", "",selected = NULL)
                    ))),
              shinyjs::hidden(
                div(id = "gsea_kegg_pathway_table",
                    box(
                      title = "KEGG pathway elements", status = "success", solidHeader = TRUE, collapsible = TRUE,
                      width = 12,
                      DT::DTOutput('gsea_keggscape_table')
                    )))
            ),
            tabPanel(
              title = tagList("Plot"),
              shinyjs::hidden(
                div(id = "keggscape_gsea_conf",
                    box( width = 4,
                          title = "Configure", status = "success", solidHeader = TRUE, collapsible = TRUE,
                          actionButton("submit_kegg_gsea_toggle", label = "Toggle colours")
                    ))),
              plotOutput("KEGG_map_gsea", height = "600px"))
  )
)