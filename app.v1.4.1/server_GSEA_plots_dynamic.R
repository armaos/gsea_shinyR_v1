observeEvent(input[[paste0("choose_GO_", gsea_cond )]], {
  # rv$okplot_gseaplot <- TRUE
  #rv$go_term <- input$choose_GO
  if( input[[paste0("choose_GO_", gsea_cond )]] != "" && length(input[[paste0("choose_GO_", gsea_cond )]]) > 0 && !is.na(input[[paste0("choose_GO_", gsea_cond )]])){  
    
    shinyjs::hide("text_gsea_plot")
    print("input$choose_GO   HERE I SHOULD THINK OF SOMETHING")
    
    rv[[paste0("gseaplot","_", gsea_cond)]] = gseaplot2(rv$gsea_list[[gsea_cond]], 
                                                                          geneSetID = which(rv$gsea_list[[gsea_cond]]$Description == input[[paste0("choose_GO_", gsea_cond )]]) , 
                                                                          title = input[[paste0("choose_GO_", gsea_cond )]])
    
    print("input$choose_GO  OK")
  }
})


output[[paste0("gseaplot","_",gsea_cond)]] <- renderPlot({
  if (!is.null(rv[[paste0("gseaplot","_",gsea_cond)]]) & !is.null(rv$gsea_list[[gsea_cond]])){
    shinyjs::enable("run")
    #rv$gsea_plots[rv$go_term]
    shinyjs::hide(paste("gsea", "gsea_plot" ,"text", gsea_cond, sep = "_"))

    rv[[paste0("gseaplot","_",gsea_cond)]]

  }else{
    empty_image()
  }
})


output[[paste0("dotplot_gsea","_",gsea_cond)]] <- renderPlot({
  if (!is.null(rv[[paste0("dotplot_gsea","_",gsea_cond)]]) & !is.null(rv$gsea_list[[gsea_cond]])){
    shinyjs::hide("text_gsea_plot_overview_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("gsea", "dotplot" ,"text", gsea_cond, sep = "_"))

    rv[[paste0("dotplot_gsea","_",gsea_cond)]]
  }else{
    empty_image()
  }
})

output[[paste0("ridgeplot_gsea","_",gsea_cond)]] <- renderPlot({
  if (!is.null(rv[[paste0("ridgeplot_gsea","_",gsea_cond)]]) & !is.null(rv$gsea_list[[gsea_cond]])){
    shinyjs::hide("text_gsea_plot_overview_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("gsea", "ridgeplot" ,"text", gsea_cond, sep = "_"))

    rv[[paste0("ridgeplot_gsea","_",gsea_cond)]]

  }else{
    empty_image()
  }
})




output$gsea_custom_upsetplot <- renderPlot({
  if (!is.null(rv$gsea_custom_upsetplot)){
    rv$gsea_custom_upsetplot
  }else{
    empty_image()
  }
})

