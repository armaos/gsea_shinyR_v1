output$gseaplot <- renderPlot({
  if(rv$okplot_gseaplot){
    shinyjs::enable("run")
    #rv$gsea_plots[rv$go_term]
    shinyjs::hide("text_gsea_plot")
    shinyjs::show("choose_GO_box")
    # if ((!is.null(rv$data ) && rv$data != "" ) && (!is.null(rv$go_term ) && rv$go_term != ""  && rv$go_term != "NA" && !is.na(rv$go_term ) ) ){
    #   ggarrange(plotlist=gsea_plots, col = 1)
    #   print("in gseaplot")
    #   shinyjs::hide("text_gsea_plot")
    # }
    if( rv$go_term != "" && length(rv$go_term) > 0 && !is.na(rv$go_term )){  
      print(paste("output$gseaplot" ,rv$go_term ))
      print(rv$gsea[1,1])
      print(which(rv$gsea$Description == rv$go_term))
      gseaplot2(rv$gsea, geneSetID = which(rv$gsea$Description == rv$go_term) , title = rv$go_term)
    }
  }
})


output$dotplot_gsea <- renderPlot({
  if (!is.null(rv$dotplot_gsea) & !is.null(rv$gsea)){
    shinyjs::hide("text_gsea_plot_overview")
    shinyjs::enable("run")
    shinyjs::hide("gsea_dotplot_text")
    
    rv$dotplot_gsea
  }else{
    empty_image()
  }
})

output$ridgeplot_gsea <- renderPlot({
  if (!is.null(rv$ridgeplot_gsea) & !is.null(rv$gsea)){
    shinyjs::hide("text_gsea_plot_overview")
    shinyjs::enable("run")
    shinyjs::hide("gsea_ridgeplot_text")
    
    rv$ridgeplot_gsea
    
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