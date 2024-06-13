linebreaks <- function(n){HTML(strrep(br(), n))}

ui <- dashboardPage(
  dashboardHeader(title = "K.A. RShiny v1.3.4"),
  dashboardSidebar(
    source(
      file = file.path(version, "ui_dashboardSidebar.R"),
      local = TRUE,
      encoding = "UTF-8"
    )$value

  ),
  # dashboardBody(
  #   # Boxes need to be put in a row (or column)
  #   tabItems(
  #     # First tab content
  #     tabItem(tabName = "menu_data_import",
  #             h2("Widgets tab content")
  # )) )) 
  
  dashboardBody(
    tags$head(tags$style(HTML(
      '.myClass {
            font-size: 20px;
            line-height: 50px;
            text-align: left;
            font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
            padding: 0 15px;
            overflow: hidden;
            color: white;
            }
            '))),
            
    tags$head(
      tags$style(type="text/css", 
      "
      .inline_textinput { display: table-row; } 
      .inline_text { 
          height : 34px;
          line-height:34px;
          display: table-row; 
          margin: 0;  
          }"
    )),
    
    tags$script('$(document).on("shiny:sessioninitialized",function(){$.get("https://api.ipify.org", function(response) {Shiny.setInputValue("getIP", response);});})'),
    tags$script(HTML(header_title_javascript)),
    useShinyjs(),
    #extendShinyjs(text = collapsa_box_jscode, functions = c("collapse") ),
    # Boxes need to be put in a row (or column)
    tabItems(
      # First tab content
      tabItem(tabName = "data_overview",
              source(
                file = file.path(version, "ui_tab_data_overview.R"),
                local = TRUE,
                encoding = "UTF-8"
              )$value
              
      ),
      tabItem(tabName = "menu_data_import",
              source(
                file = file.path(version, "ui_tab_menu_data_import.R"),
                local = TRUE,
                encoding = "UTF-8"
              )$value
      ),
      # Second tab content
      tabItem(tabName = "menu_norm_and_de",
              source(
                file = file.path(version, "ui_tab_norm_and_de.R"),
                local = TRUE,
                encoding = "UTF-8"
              )$value
      ),
      # Third tab content
      tabItem(tabName = "set_selection",
              source(
                file = file.path(version, "ui_tab_set_selection.R"),
                local = TRUE,
                encoding = "UTF-8"
              )$value
      ),
      # GSO
      tabItem(tabName = "gso_plots",
              p(id = "text_gso_plots", "Something here")),
      tabItem(tabName = "KEGGscape_gso",
              source(
                file = file.path(version, "ui_tab_KEGGscape_gso.R"),
                local = TRUE,
                encoding = "UTF-8"
              )$value
      ),
      
      
      tabItem(tabName = "configure_ora",
              source(
                file = file.path(version, "ui_tab_configure_ora.R"),
                local = TRUE,
                encoding = "UTF-8"
              )$value
      ),
      tabItem(tabName = "table_ora",
              p(id = "text_ora_table", "Please submit your configuration first!"),

      ),
      
      tabItem(tabName = "plots_ora",
              source(
                file = file.path(version, "ui_tab_plots_ora.R"),
                local = TRUE,
                encoding = "UTF-8"
              )$value

              
      ),
      # tabItem(tabName = "KEGGscape_ora",
                # source(
                #   file = file.path(version, "ui_tab_KEGGscape_ora.R"),
                #   local = TRUE,
                #   encoding = "UTF-8"
                # )$value
      # ),
      tabItem(tabName = "table_gsea",
              source(
                file = file.path(version, "ui_tab_table_gsea.R"),
                local = TRUE,
                encoding = "UTF-8"
              )$value
             
      ),
      tabItem(tabName = "plots_gsea",
              source(
                file = file.path(version, "ui_tab_plots_gsea.R"),
                local = TRUE,
                encoding = "UTF-8"
              )$value
  
      ),
      tabItem(tabName = "running_score_gsea",
              source(
                file = file.path(version, "ui_tab_running_score_gsea.R"),
                local = TRUE,
                encoding = "UTF-8"
              )$value
    
      ),
      tabItem(tabName = "custom_upset_gsea",
              source(
                file = file.path(version, "ui_tab_custom_upset_gsea.R"),
                local = TRUE,
                encoding = "UTF-8"
              )$value
      )
      # tabItem(tabName = "KEGGscape_gsea",
      #         source(
      #           file = file.path(version, "ui_tab_KEGGscape_gsea.R"),
      #           local = TRUE,
      #           encoding = "UTF-8"
      #         )$value
      # )
    )
    
    
  )
)
