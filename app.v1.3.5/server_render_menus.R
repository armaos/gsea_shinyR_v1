output$menu_norm_and_de_mout <- renderMenu({
  if(!is.null(rv$show_norm_and_de_mout ))
    menuItem("Normalization and DE", 
             menuSubItem("Plots", tabName = "menu_norm_and_de"), 
             icon = icon("th"))
})

output$gso_mout <- renderMenu({
  if(!is.null(rv$show_gso_mout))
    menuItem("Gene set Analysis",
             menuSubItem("Gene set selection", tabName = "set_selection"),
             # menuItem("Gene set Analytics",
             #          menuSubItem("Plots", tabName = "gso_plots"),
             #          #menuSubItem("KEGGscape", tabName = "KEGGscape_gso"),
             #          icon = icon("th")),
             menuItemOutput("ora_mout"),
             icon = icon("th"))
})

output$ora_mout <- renderMenu({
  if(!is.null(rv$show_ora_mout )){
    menuItem("Enrichment Analysis (ORA)",
             menuSubItem("Configuration And Table", tabName = "configure_ora"),
             #menuSubItem("Plots", tabName = "plots_ora"),
             menuSubItem("Plots", tabName = "plots_ora_dynamic"),
             icon = icon("th"))
  }
})


output$gsea_mout <- renderMenu({
  if(!is.null(rv$show_gsea_mout )){
    menuItem("GSEA",
             menuSubItem("Configuration And Table", tabName = "table_gsea"),
             menuItemOutput("gsea_plots_mout"),
             icon = icon("th")
    )
  }
})

output$gsea_plots_mout <- renderMenu({
  if(!is.null(rv$show_gsea_plots_mout)){
    menuItem("GSEA Plots",
             #menuSubItem("Plots", tabName = "plots_gsea"),
             menuSubItem("Plots", tabName = "plots_gsea_dynamic"),
             #menuSubItem("Running Score", tabName = "running_score_gsea"),
             menuSubItem("Running Score", tabName = "running_score_gsea_dynamic"),
             menuSubItem("Custom UpSet", tabName = "custom_upset_gsea"),
             icon = icon("th"))
  }
})

output$KeggScape_mout <- renderMenu({
  if(!is.null(rv$show_KeggScape_mout))
    menuItem("KEGGscape", tabName = "KEGGscape_gso", icon = icon("th"))
})
