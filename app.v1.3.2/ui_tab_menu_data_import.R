fluidRow(
  shinyjs::hidden(div(
    id = "data_import_box",
    box(
      title = "Select Data source", status = "primary", solidHeader = TRUE,
      
      selectInput("choose", "Choose:", choices = c("Pre-compiled", "Custom"), selected = NULL),
      #actionButton("demo", "Fill in fields with demo files"),
      conditionalPanel(
        "input.choose=='Custom'",
        #fileInput("data_file", "Read Counts:", multiple = TRUE, accept = c("text/csv","text/comma-separated-values,text/plain")),
        selectInput("choose_custom_tcmp", "Choose:", choices = c("MSI" ="testdata", "Galaxy run" = "galaxy"), selected = "testdata"),
        conditionalPanel(
          "input.choose_custom_tcmp=='galaxy'",
          box(
            id = "account_box", width = 12, title = "Account", status = "primary", solidHeader = TRUE, collapsible = TRUE,
            fluidRow(
              column(6, textAreaInput("username", "Give your username", "Marton", rows = 1)),
              column(3, style = "margin-top: 25px;", actionButton("submit_username", label = "Submit"))
            )),
          shinyjs::hidden(
            div(id = "galaxy_submissions_box",
                box(
                  width = 12, title = "Galaxy submissions", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                  fluidRow(
                    column(6, selectInput("galaxy_runs", "Your Galaxy runs:",choices = c(), multiple = T)),
                    column(3, style = "margin-top: 25px;", actionButton("submit_galaxy_runs", label = "Fetch!"))
                  )))),
          shinyjs::hidden(
            div(id = "galaxy_outputs_box",
                box(
                  width = 12, title = "Galaxy outputs", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                  fluidRow(
                    column(6, selectInput("galaxy_outputs", "Select your count data", choices = c(), multiple = T)),
                    column(3, style = "margin-top: 25px;", actionButton("submit_galaxy_outputs", label = "Combine!"))
                  ))))
        ),
        conditionalPanel(
          "input.choose_custom_tcmp=='testdata'",
          actionButton("submit_norm_and_de", label = "Submit Data Pre-processing"),
        ),
        
        shinyjs::hidden(
          div(id = "menu_de_box",
              selectInput("choose_conditions", "Choose conditions:",choices = c()),
              actionButton("run_DE", label = "Submit DE analysis")
          ))
        # DT::DTOutput('DE_table')
      ),
      conditionalPanel(
        "input.choose=='Pre-compiled'",
        h3("Select pre calculated condition:"),
        radioButtons("Condition", "Choose one:",
                     c("D_O" = "D_O",
                       "L_O" = "L_O",
                       "D_L" = "D_L",
                       "B_G_in_O" = "B_G_in_O",
                       "B_G_in_L" = "B_G_in_L",
                       "B_G_in_D" = "B_G_in_D"
                     )),
        actionButton("submit_contrast", label = "Submit Contrast"),
        shinyjs::hidden(h2(id = "submitted_contrast1", "Contrast submitted Successfully.")),
        #shinyjs::hidden(p(id = "submitted_contrast2", "You can head to ORA and GSEA Analysis now!")),
        
      )
    ))),
  p(id = "overview_text", "Initializing app...")
)
