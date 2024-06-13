output[[paste0("gseaplot","_",id)]] <- renderPlot({
  if (!is.null(rv[[paste0("gseaplot","_",id)]]) & !is.null(rv$gsea)){
    shinyjs::enable("run")
    #rv$gsea_plots[rv$go_term]
    shinyjs::hide(paste("gsea", "gsea_plot" ,"text", id, sep = "_"))
    
    rv[[paste0("gseaplot","_",id)]]
    
  }else{
    empty_image()
  }
})


output[[paste0("dotplot_gsea","_",id)]] <- renderPlot({
  if (!is.null(rv[[paste0("dotplot_gsea","_",id)]]) & !is.null(rv$gsea)){
    shinyjs::hide("text_gsea_plot_overview_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("gsea", "dotplot" ,"text", id, sep = "_"))
    
    rv[[paste0("dotplot_gsea","_",id)]]
  }else{
    empty_image()
  }
})

output[[paste0("ridgeplot_gsea","_",id)]] <- renderPlot({
  if (!is.null(rv[[paste0("ridgeplot_gsea","_",id)]]) & !is.null(rv$gsea)){
    shinyjs::hide("text_gsea_plot_overview_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("gsea", "ridgeplot" ,"text", id, sep = "_"))
    
    rv[[paste0("ridgeplot_gsea","_",id)]]
    
  }else{
    empty_image()
  }
})



#### _custom_upset_plot
my_custom_upset_plot <- eventReactive(input$submit_custom_upset_plot,{
  as.data.frame(isolate(rv$gsea)) %>% 
    filter(Description %in% isolate(input$custom_upset_terms)) %>%
    select(Description, core_enrichment) %>% 
    separate_rows(core_enrichment, sep = "/") %>% 
    mutate(Description = ifelse(str_length(Description) > 30, paste0(str_sub(Description,1,30), "..."), Description)) %>%
    group_by(core_enrichment) %>% 
    summarise(terms = list(Description)) %>% 
    ggplot(aes(x=terms)) +
    geom_bar() +
    scale_x_upset()
})

output$gsea_custom_upsetplot <- renderPlot({
  my_custom_upset_plot()
})

#### 