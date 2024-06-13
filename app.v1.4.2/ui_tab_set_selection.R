fluidRow(
  box(
    width = 8,
    title = "Select your query Gene set", status = "primary", solidHeader = TRUE,
    
    selectInput("choose_gso", "Choose:", choices = c("Upon filtering", "Custom Set", "By group"), selected = NULL),
    
    
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
                  max = 12, value = c(1.5, 12)),
      
      #sliderTextInput("slider_fdr", label = p("FDR treshold"),choices = seq(from = 20, to = 0, by = -1),selected = 10, pre = "1e-"),
      #checkboxInput("fdr_by5", "x5: ", TRUE),
      #div(#style="vertical-align: middle",
      fluidRow(
        column(12, 
        div(style="display:inline-block;vertical-align: bottom; width: 70%;", 
                      sliderTextInput("slider_fdr", label = p("FDR treshold"),choices = seq(from = 20, to = 0, by = -1),selected = 2, pre = "1e-")),
        div(style="display:inline-block;vertical-align: bottom; width: 10%;",
                      checkboxInput("fdr_by5", "x5: ", TRUE)),
        div(style="display:inline-block;vertical-align: bottom;margin-bottom:15px; width: 10%;",
                      textOutput("fdr_value", inline=T)),
      ))
       #div(style="display:inline-block;vertical-align: bottom;", sliderTextInput("slider_fdr", label = p("FDR treshold"),choices = seq(from = 20, to = 0, by = -1),selected = 10, pre = "1e-")),
       #div(style="display:inline-block;vertical-align: bottom;",checkboxInput("fdr_by5", "x5: ", TRUE)),
      # div(style="display:inline-block;vertical-align: bottom; width: 40px;",verbatimTextOutput("fdr_value"))
      #)
      ,
      
      verticalLayout(
        actionButton("submit_thresholds", label = "Submit Thresholds"))
      
    ),
    conditionalPanel(
      "input.choose_gso=='By group'",
      fluidRow(
        column(width = 6, pickerInput("gene_select_by_group", "Grouping by:", choices = c(), multiple = F))
      ),
      fluidRow(
        column(width = 12, pickerInput("gene_select_by_group_values", "Value:", choices = c(), multiple = TRUE))
      ),
      linebreaks(1),
      shinyjs::disabled(actionButton("submit_gene_select_by_group", label = "Submit selection by group"))
      
      
    ),
  ),
  shinyjs::hidden(
    div(
      id = "summary_gene_set_box",
      box(
        width = 4,
        title = "Summary", status = "info", solidHeader = TRUE,
        htmlOutput("summary_gene_set")
      ))),
  
  ### a table like the one in the data overview with the selected data overview
  source(file = file.path(version, "ui_tab_selected_data_overview.R"),
         local = TRUE,
         encoding = "UTF-8")
  
)