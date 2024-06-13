fluidRow(
  shinyjs::hidden(div(
    id = "data_import_box",
    box(
      width = 12,
      title = "Select Data source", status = "primary", solidHeader = TRUE,
      #column(3, p("Choose: ")),
      column(3, style = "margin-top: 7.5px;", p("Choose")),
      column(9, selectInput("choose", NULL, choices = c("Pre-compiled", "Custom"), selected = NULL )),
      
      #actionButton("demo", "Fill in fields with demo files"),
      conditionalPanel(
        "input.choose=='Custom'",
        #fileInput("data_file", "Read Counts:", multiple = TRUE, accept = c("text/csv","text/comma-separated-values,text/plain")),
        column(3, style = "margin-top: 7.5px;", p("Choose")),
        column(9, selectInput("choose_custom_tcmp", NULL, choices = c("MSI" ="testdata", "Galaxy run" = "galaxy"), selected = "testdata")), 
        
        conditionalPanel(
          "input.choose_custom_tcmp=='galaxy'",
          box(
            id = "account_box", width = 12, title = "Account", status = "primary", solidHeader = TRUE, collapsible = TRUE,
            fluidRow(
              column(3, style = "margin-top: 7.5px;", p("Username:")),
              column(6, textAreaInput("username", NULL, "Marton", rows = 1)),
              #column(3, style = "margin-top: 25px;", actionButton("submit_username", label = "Submit"))
              column(2, actionButton("submit_username", label = "Submit"))
            )),
          shinyjs::hidden(
            div(id = "galaxy_submissions_box",
                box(
                  width = 12, title = "Galaxy submissions", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                  fluidRow(
                    column(3, style = "margin-top: 7.5px;", p("Galaxy runs:")),
                    column(5, selectInput("galaxy_runs", NULL ,choices = c(), multiple = T)),
                    column(2,  shinyjs::hidden(actionButton("submit_galaxy_runs", label = "Fetch!")))
                  )))
          )
          ,
          shinyjs::hidden(
            div(id = "galaxy_outputs_box",
                box(
                  width = 12, title = "Galaxy outputs", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                  fluidRow(
                    column(3, style = "margin-top: 7.5px;", p("Select your count data")),
                    column(9, 
                           actionButton("selectall_galaxy", label="Select/Deselect all"),
                           #selectInput("galaxy_outputs", NULL, choices = c("SRR1207056.tabular","SRR7637235.tabular", "SRR7637208.tabular", "SRR7637249.tabular", "SRR7637250.tabular", "SRR7637252.tabular","SRR7637253.tabular" ), multiple = T)
                           selectInput("galaxy_outputs", NULL, choices = c( ), multiple = T)
                    ),
                    #column(9, selectInput("galaxy_outputs", NULL, choices = c(), multiple = T)),
                    #column(3, style = "margin-top: 25px;", actionButton("submit_galaxy_outputs", label = "Combine!"))
                  ),
                  radioButtons("metadata_input_type", "Define Metadata:",
                               c("Manual" = "manual",
                                 "File input" = "file",
                                 "Galaxy pre-defined" = "galaxy",
                                 "Text input" = "text"
                               )),
                  conditionalPanel(
                    "input.metadata_input_type=='manual'",
                    fluidRow(
                      column(6, uiOutput("textbox_ui")),
                      column(3, uiOutput("textbox_donor_ui")),
                      column(3, uiOutput("textbox_group_ui"))
                      #column(7, tags$div(id="inline1", class="inline_textinput", uiOutput("textbox_donor_ui"))),
                      #column(5, tags$div(id="inline2", class="inline_textinput", uiOutput("textbox_group_ui")))
                    )),
                  conditionalPanel(
                    "input.metadata_input_type=='file'",
                    fileInput("metadata_file", "Choose File",
                              multiple = FALSE,
                              accept = c("text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")),
                    
                    # Input: Checkbox if file has header ----
                    checkboxInput("metadata_header", "Is there a header in the file?", TRUE),
                    
                  ),
                  conditionalPanel(
                    "input.metadata_input_type=='text'",
                    textAreaInput("metadata_text", "Define meetadata (3 columns per line)", "SampleID Subgroup Treatment", height = "100px"),
                  ),
                  conditionalPanel(
                    "input.metadata_input_type=='galaxy'",
                    
                    p("Under construction!")
                  ),
                  linebreaks(1),
                  shinyjs::hidden(div( id = "metadata_error_p", p(id = "metadata_error_message"))),
                  shinyjs::hidden(div(
                    id = "metadata_table_box",
                    box(width = 12, title = "Metadata table", status = "primary", solidHeader = TRUE,
                        DT::DTOutput('metadata_table')))),
                  shinyjs::disabled(actionButton("submit_metadata", label = "Confirm metadata")),
                  shinyjs::disabled(actionButton("submit_galaxy_outputs", label = "Submit"))
                  
                )))
        )
        ,
        conditionalPanel(
          "input.choose_custom_tcmp=='testdata'",
          actionButton("submit_custom_data_preprocessing", label = "Submit Data Pre-processing"),
        ),
        
        shinyjs::hidden(
          div(id = "menu_de_box",
              selectInput("choose_conditions", "Choose contrast:",choices = c(),  multiple = T),
              actionButton("run_DE", label = "Submit DE analysis")
          ))
        # DT::DTOutput('DE_table')
      ),
      conditionalPanel(
        "input.choose=='Pre-compiled'",
        h3("Select pre calculated condition:"),
        selectInput("Condition", "Choose one:",
                    choices = c("Low vs Medium" = "D_O",
                                "High vs Medium" = "L_O",
                                "Low vs High" = "D_L",
                                "Brown vs Green in Medium" = "B_G_in_O",
                                "Brown vs Green in High" = "B_G_in_L",
                                "Brown vs Green in Low" = "B_G_in_D"
                    ),
                    multiple = T),
        shinyjs::disabled(actionButton("submit_contrast", label = "Submit Contrast")),
        shinyjs::hidden(h2(id = "submitted_contrast1", "Contrast submitted Successfully.")),
        #shinyjs::hidden(p(id = "submitted_contrast2", "You can head to ORA and GSEA Analysis now!")),
        
      )
    ))),
  shinyjs::hidden(p(id = "overview_text", "Initializing app..."))
)
