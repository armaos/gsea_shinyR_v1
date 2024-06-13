fluidRow(
  h2("Gene Set Enrichment Analysis"),

  div(id = "run_gsea_box",
      box(
        title = "Term Enrichment", status = "success", solidHeader = TRUE, collapsible = TRUE,
        radioButtons("enrichment_term_gsea", label = "Select enrichment",
                    choiceNames = c("GO:Biological Processes", "GO:Molecular Functions", "GO:Cellular Components",  
                                    "GO:Biological Processes (Atha)", "GO:Molecular Functions (Atha)", "GO:Cellular Components (Atha)",  
                                    "KEGG Pathways"),
                    choiceValues = c("biological_process", "molecular_function", "cellular_component", 
                                      "biological_process_arth", "molecular_function_arth", "cellular_component_arth", 
                                      "KEGG_Pathways")),
        selectInput("choose_gsea_condition", "Choose Condition:", 
                    choices = c(), multiple = F),
        
        actionButton("run_gsea", label = "Run GSEA!"))),
  shinyjs::hidden(p(id = "text_configure_gsea", "You haven't submitted any contrast yet!")),
  shinyjs::hidden(p(id = "text_gsea_table")),
  shinyjs::hidden(
    div(id = "gsea_table_box", 
        box(title = "Data", width = 12, solidHeader = TRUE, status = "primary", 
            collapsible = TRUE,
            DT::DTOutput('gsea_table')
        ))
  ),
  #DT::DTOutput('gsea_table')
)