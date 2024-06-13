ora_tab_list %>%
  walk(~removeTab("tabs_plots_ora_dynamic", .x))
ora_tab_list <<- NULL

for(id in rv$condition){
  ora_tab_list <<- c(ora_tab_list, id)
  
  # here the plots for ora
  appendTab(inputId = "tabs_plots_ora_dynamic",
            tabPanel(id,
                     tabsetPanel(type = "tabs",
                                 tabPanel("Bar plot - Gene Counts",
                                          p(id = paste0("ora_barplot1_text_", id), "Processing..."),
                                          plotOutput(paste0("barplot1_", id), height = "600px")),
                                 tabPanel("Bar plot - qscore",
                                          p(id = paste0("ora_barplot2_text_", id), "Processing..."),                            
                                          plotOutput(paste0("barplot2_", id), height = "600px")),
                                 tabPanel("Dot Plot",
                                          p(id = paste0("ora_dotplot_text_", id), "Processing..."),                            
                                          plotOutput(paste0("dotplot_", id), height = "600px")),
                                 tabPanel("Enrichment Map Plot",
                                          p(id = paste0("ora_emapplot_text_", id), "Processing..."),                            
                                          plotOutput(paste0("emapplot_", id), height = "600px")),
                                 tabPanel("Gene-Concept Network Plot",
                                          p(id = paste0("ora_cnetplot_text_", id), "Processing..."),                            
                                          plotOutput(paste0("cnetplot_", id), height = "600px")),
                                 tabPanel("HeatMap Plot",
                                          p(id = paste0("ora_heatplot_text_", id), "Processing..."),                            
                                          plotOutput(paste0("heatplot_", id), height = "600px")),
                                 tabPanel("Tree Plot",
                                          p(id = paste0("ora_treeplot_text_", id), "Processing..."),                            
                                          plotOutput(paste0("treeplot_", id), height = "600px")),
                                 tabPanel("UpSetPlot",
                                          p(id = paste0("ora_upsetplot_text_", id), "Processing..."),                            
                                          plotOutput(paste0("upsetplot_", id)))  
                     )
            )
            
  )
  
  #here the plots for GSEA
  appendTab(inputId = "tabs_plots_gsea_dynamic",
            tabPanel(id,
                     tabsetPanel(type = "tabs",
                                 tabPanel("Dot Plot",
                                          p(id = paste0("gsea_dotplot_text_", id), "Processing..."),                                                      
                                          plotOutput(paste0("dotplot_gsea_", id), height = "600px") ),
                                 tabPanel("Ridge Plot",
                                          p(id = paste0("gsea_ridgeplot_text_", id), "Processing..."),                                                                                
                                          plotOutput(paste0("ridgeplot_gsea_", id), height = "600px"))
                     )
            )
  )
                
  # here for the running score
  appendTab(inputId = "tabs_running_score_dynamic",
            tabPanel(id,
                     plotOutput(paste0("gseaplot_", id), height = "600px")
            )
  )
  
  
  for(new_rv in c("barplot1", "barplot2", "dotplot", "emapplot", "cnetplot", "heatplot", "treeplot", "upsetplot", "dotplot_gsea", "ridgeplot_gsea", "gseaplot")){
    rv[[paste0(new_rv,"_",id)]] <- NULL
  }
  
  source(file = file.path(version, "server_ORA_plots_dynamic.R"),
         local = TRUE,
         encoding = "UTF-8")

  source(file = file.path(version, "server_GSEA_plots_dynamic.R"),
         local = TRUE,
         encoding = "UTF-8")
  


}
