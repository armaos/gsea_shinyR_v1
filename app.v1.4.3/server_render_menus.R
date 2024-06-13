output$menu_fastqc_mout <- renderMenu({
  if(!is.null(rv$show_menu_fastqc ) & rv$show_menu_fastqc  == T)
  {
    menuItem("Quality Control", tabName = "fastqc_menu",
             menuSubItem("Plots", tabName = "fastqc_tab"), 
             selectInput("fastqc_input", "Select a Run acc for the individual plots", choices = c(), selected = NULL),
             icon = icon("check-square-o"))
  }else{
    shinyjs::hide(selector = "a[data-value='fastqc_menu']" )
    #menuItem(NULL)
  }
})


output$menu_norm_and_de_mout <- renderMenu({
  if(!is.null(rv$show_norm_and_de_mout ))
  {
    menuItem("Normalization and DE", 
             menuSubItem("Plots", tabName = "menu_norm_and_de"), 
             icon = icon("chart-bar"))
  }
})


output$gso_mout <- renderMenu({
  if(!is.null(rv$show_gso_mout))
  {
    menuItem("Gene set Analysis",
             menuSubItem("Gene set selection", tabName = "set_selection"),
             # menuItem("Gene set Analytics",
             #          menuSubItem("Plots", tabName = "gso_plots"),
             #          #menuSubItem("KEGGscape", tabName = "KEGGscape_gso"),
             #          icon = icon("th")),
             menuItemOutput("ora_mout"),
             icon = icon("filter"))
  }
})

output$ora_mout <- renderMenu({
  if(!is.null(rv$show_ora_mout ))
  {
    menuItem("Enrichment Analysis (ORA)",
             menuSubItem("Configuration And Table", tabName = "configure_ora"),
             #menuSubItem("Plots", tabName = "plots_ora"),
             menuSubItem("Plots", tabName = "plots_ora_dynamic"),
             icon = icon("chart-bar"))
  }
})


output$gsea_mout <- renderMenu({
  if(!is.null(rv$show_gsea_mout ))
  {
    menuItem("GSEA",
             menuSubItem("Configuration And Table", tabName = "table_gsea"),
             menuItemOutput("gsea_plots_mout"),
             icon = icon("chart-bar")
    )
  }
})

output$gsea_plots_mout <- renderMenu({
  if(!is.null(rv$show_gsea_plots_mout))
  {
    menuItem("GSEA Plots",
             #menuSubItem("Plots", tabName = "plots_gsea"),
             menuSubItem("Plots", tabName = "plots_gsea_dynamic"),
             #menuSubItem("Running Score", tabName = "running_score_gsea"),
             menuSubItem("Running Score", tabName = "running_score_gsea_dynamic"),
             menuSubItem("Custom UpSet", tabName = "custom_upset_gsea"),
             icon = icon("chart-pie"))
  }
})

output$KeggScape_mout <- renderMenu({
  if(!is.null(rv$show_KeggScape_mout))
  {
    menuItem("KEGGscape", tabName = "KEGGscape_gso", icon = icon("circle-nodes"))
  }
  
})


output$Clustering_mout <- renderMenu({
  if(!is.null(rv$show_clustering_mout))
  {
    menuSubItem("Clustering", tabName = "clustering_tab", icon = icon("delicious"))
  }
})
