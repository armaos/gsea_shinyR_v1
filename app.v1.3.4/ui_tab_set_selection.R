fluidRow(
  box(
    width = 6,
    title = "Select your query Gene set", status = "primary", solidHeader = TRUE,
    
    selectInput("choose_gso", "Choose:", choices = c("Upon filtering", "Custom Set"), selected = NULL),
    
    
    conditionalPanel(
      "input.choose_gso=='Custom Set'",
      selectInput("choose_custom_gso", "Choose:", choices = c("Text", "File"), selected = NULL),
      conditionalPanel(
        "input.choose_custom_gso=='File'",
        fileInput("gso_file_input", "Uplaod File with Gene names:", multiple = TRUE, accept = c("text/csv","text/comma-separated-values,text/plain"))),
      conditionalPanel(
        "input.choose_custom_gso=='Text'",
        textAreaInput("gso_text_input", "Paste your genes"),
        actionButton("demo", "Fill in fields with demo files")),
      actionButton("submit_custom_gene_set", label = "Submit custom Gene Set")
    ),
    conditionalPanel(
      "input.choose_gso=='Upon filtering'",
      sliderInput("slider_lfc", label = p("|logFC| Range"), min = 0, step = 0.5,
                  max = 12, value = c(4, 12)),
      
      sliderTextInput("slider_fdr", label = p("FDR treshold"),
                      choices = seq(from = 20, to = 0, by = -1),
                      selected = 10, pre = "1e-"),
      verticalLayout(
        actionButton("submit_thresholds", label = "Submit Thresholds"))
      
    ),
  ),
  shinyjs::hidden(
    div(
      id = "summary_gene_set_box",
      box(
        width = 6,
        title = "Summary", status = "info", solidHeader = TRUE,
        htmlOutput("summary_gene_set")
      )))
)