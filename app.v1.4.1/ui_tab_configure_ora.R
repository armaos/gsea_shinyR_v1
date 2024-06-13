fluidRow(

  shinyjs::hidden(p(id = "text_configure_ora", "You haven't submitted any contrast yet!")),
  shinyjs::hidden(
    div(id = "box_configure_term_ora",
        box(
          width = 4,
          title = "Term Enrichment", status = "danger", solidHeader = TRUE, collapsible = TRUE,
          radioButtons("enrichment_term", label = "Select enrichment",
                      choiceNames = c("GO:Biological Processes", "GO:Molecular Functions", "GO:Cellular Components",  
                                      "GO:Biological Processes (Atha)", "GO:Molecular Functions (Atha)", "GO:Cellular Components (Atha)",  
                                      "KEGG Pathways"),
                      choiceValues = c("biological_process", "molecular_function", "cellular_component", 
                                        "biological_process_arth", "molecular_function_arth", "cellular_component_arth", 
                                        "KEGG_Pathways")
          ),
          selectInput("choose_ora_condition", "Choose Condition:", 
                      choices = c(), multiple = F),
          selectInput("choose_universe", "Choose Universe:", 
                      choices = c("All transcripts" = "all_t", 
                                  #"All annotated transcripts" = "all_annot_t", 
                                  "All DE transcripts" = "all_de") , 
                      selected = "all_t"),
          sliderInput("ORA_minGSSize",label = p("Minimal size of genes annotated for testing"), min = 5, max = 50, value = 10),
          sliderInput("ORA_maxGSSize",label = p("Maximal size of genes annotated for testing"), min = 100, max = 1000, value = 500),
          actionButton("run_ora", label = "Submit enrichment")
        ))),
  shinyjs::hidden(
    div(id = "box_ora_send_to_plot",
        box(
          width = 4,
          title = "Send to plot", status = "success", solidHeader = TRUE, collapsible = TRUE,
          sliderInput("ORA_max_elements",
                      label = p("Maximum terms per plot to show"), min = 1,
                      max = 50, value = 20),
          actionButton("submit_plot_ora", label = "Plot table")
        ))),
  shinyjs::hidden(p(id = "text_ora_table")),
  shinyjs::hidden(div(
    id = "ora_table_box",
    box(title = "Data", width = 12, solidHeader = TRUE, status = "primary", 
        collapsible = TRUE,
        DT::DTOutput('ora_table')
    )
  ))
  #DT::DTOutput('ora_table')
)