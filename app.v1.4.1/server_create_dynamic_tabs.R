ora_tab_list %>%
  walk(~removeTab("tabs_plots_ora_dynamic", .x))
ora_tab_list <<- NULL

for(id in rv$condition){
  #ora_tab_list <<- c(ora_tab_list, id)
  
  # here the plots for norm and de
  appendTab(inputId = "volcano_norm_de_dynamic",
            tabPanel(id,
                     plotOutput(paste0("volcano0_", id), height = "600px"))
  )
  
  appendTab(inputId = "heatmap_norm_de_dynamic",
            tabPanel(id,
                     plotOutput(paste0("pheatmap", id), height = "600px"))
  )
  

  
  source(file = file.path(version, "server_NORM_plots_dynamic.R"),
         local = TRUE,
         encoding = "UTF-8")
  



}
