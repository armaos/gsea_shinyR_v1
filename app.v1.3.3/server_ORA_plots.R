output$barplot1 <- renderPlot({
  if (!is.null(rv$barplot1) & !is.null(rv$ora_desc)){
    shinyjs::hide("text_ora_plots")
    shinyjs::enable("run")
    shinyjs::hide("ora_barplot1_text")
    rv$barplot1
  }else{
    empty_image()
  }
})

output$barplot2 <- renderPlot({
  if (!is.null(rv$barplot2) & !is.null(rv$ora_desc)){
    shinyjs::hide("text_ora_plots")
    shinyjs::enable("run")
    shinyjs::hide("ora_barplot2_text")      
    rv$barplot2
  }else{
    empty_image()
  }
})

output$dotplot <- renderPlot({
  if(!is.null(rv$dotplot) & !is.null(rv$ora_desc)){
    shinyjs::enable("run")
    shinyjs::hide("ora_dotplot_text")
    rv$dotplot
  }else{
    shinyjs::enable("run")
    html("text_ora_plots", "<strong>Make your configurarion for plotting from the previous tab panel first!</strong>")
    empty_image()
  }
})

output$emapplot <- renderPlot({
  if (!is.null(rv$emapplot) & !is.null(rv$ora_desc)){
    shinyjs::hide("text_ora_plots")
    shinyjs::enable("run")
    shinyjs::hide("ora_emapplot_text")
    
    rv$emapplot
  }else{
    empty_image()
  }
})

output$cnetplot <- renderPlot({
  if (!is.null(rv$cnetplot) & !is.null(rv$ora_desc)){
    shinyjs::hide("text_ora_plots")
    shinyjs::enable("run")
    shinyjs::hide("ora_cnetplot_text")
    
    rv$cnetplot
  }else{
    empty_image()
  }
})
output$heatplot <- renderPlot({
  if (!is.null(rv$heatplot) & !is.null(rv$ora_desc)){
    shinyjs::hide("text_ora_plots")
    shinyjs::enable("run")
    shinyjs::hide("ora_heatplot_text")
    
    rv$heatplot
  }else{
    empty_image()
  }
})

output$treeplot <- renderPlot({
  if (!is.null(rv$treeplot) & !is.null(rv$ora_desc)){
    shinyjs::hide("text_ora_plots")
    shinyjs::enable("run")
    shinyjs::hide("ora_treeplot_text")
    
    rv$treeplot
  }else{
    empty_image()
  }
})


output$upsetplot <- renderPlot({
  if (!is.null(rv$upsetplot) & !is.null(rv$ora_desc)){
    shinyjs::hide("text_ora_plots")
    shinyjs::enable("run")
    shinyjs::hide("ora_upsetplot_text")
    
    rv$upsetplot
  }else{
    empty_image()
  }
})
