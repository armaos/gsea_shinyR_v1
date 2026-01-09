if(!ora_cond %in% ora_tab_list){
  ora_tab_list <<- c(ora_tab_list, ora_cond)
  # here the plots for ora
  appendTab(inputId = "tabs_plots_ora_dynamic",
            tabPanel(ora_cond,
                     p(id = "text_tabs_plots_ora_dynamic"),
                     tabsetPanel(type = "tabs", id = paste0("ora_plot_tabsetpanel_",ora_cond),
                                 tabPanel("Bar plot - Gene Counts",
                                          #p(id = paste0("ora_barplot1_text_", ora_cond), "Processing..."),
                                          withSpinner(plotOutput(paste0("barplot1_", ora_cond), height = "600px"), type = 6)),
                                 tabPanel("Bar plot - qscore",
                                          #p(id = paste0("ora_barplot2_text_", ora_cond), "Processing..."),                            
                                          withSpinner(plotOutput(paste0("barplot2_", ora_cond), height = "600px"), type = 6)),
                                 tabPanel("Dot Plot",
                                          #p(id = paste0("ora_dotplot_text_", ora_cond), "Processing..."),                            
                                          withSpinner(plotOutput(paste0("dotplot_", ora_cond), height = "600px"), type = 6)),
                                 tabPanel("Enrichment Map Plot",
                                          #p(id = paste0("ora_emapplot_text_", ora_cond), "Processing..."),                            
                                          withSpinner(plotOutput(paste0("emapplot_", ora_cond), height = "600px"), type = 6)),
                                 tabPanel("Gene-Concept Network Plot",
                                          #p(id = paste0("ora_cnetplot_text_", ora_cond), "Processing..."),                            
                                          withSpinner(plotOutput(paste0("cnetplot_", ora_cond), height = "600px"), type = 6)),
                                 tabPanel("HeatMap Plot",
                                          #p(id = paste0("ora_heatplot_text_", ora_cond), "Processing..."),                            
                                          withSpinner(plotOutput(paste0("heatplot_", ora_cond), height = "600px"), type = 6)),
                                 tabPanel("Tree Plot",
                                          p(id = paste0("ora_treeplot_text_", ora_cond), "Processing..."),                            
                                          withSpinner(plotOutput(paste0("treeplot_", ora_cond), height = "600px"), type = 6)),
                                 tabPanel("UpSetPlot",
                                          p(id = paste0("ora_upsetplot_text_", ora_cond), "Processing..."),                            
                                          withSpinner(plotOutput(paste0("upsetplot_", ora_cond)), type = 6))  
                     )
            )
            
  )
  
  
  for(new_rv in c("barplot1", "barplot2", "dotplot", "emapplot", "cnetplot", "heatplot", "treeplot", "upsetplot", "dotplot_gsea", "ridgeplot_gsea", "gseaplot")){
    rv[[paste0(new_rv,"_",ora_cond)]] <- NULL
  }
  
  source(file = file.path(version, "server_ORA_plots_dynamic.R"),
         local = TRUE,
         encoding = "UTF-8")  
}
