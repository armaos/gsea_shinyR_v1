
sidebarMenu(
  menuItem("Data", 
           menuSubItem("Import", tabName = "menu_data_import"),
           menuSubItem("Overview", tabName = "data_overview"),
           icon = icon("th")),
  menuItemOutput("menu_fastqc_mout"),
  
  menuItemOutput("menu_norm_and_de_mout"),
  
  menuItemOutput( "gso_mout"),
  menuItemOutput( "Clustering_mout"),
  menuItemOutput( "KeggScape_mout"),
  
  menuItemOutput( "gsea_mout")
  
)