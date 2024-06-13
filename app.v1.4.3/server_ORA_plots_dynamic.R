output[[paste0("barplot1","_",ora_cond)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("barplot1","_",ora_cond)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "barplot1" ,"text", ora_cond, sep = "_"))
    rv[[paste0("barplot1","_",ora_cond)]]
  }else{
    empty_image()
  }
})


output[[paste0("barplot2","_",ora_cond)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("barplot2","_",ora_cond)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "barplot2" ,"text", ora_cond, sep = "_"))
    rv[[paste0("barplot2","_",ora_cond)]]
  }else{
    empty_image()
  }
})


output[[paste0("dotplot","_",ora_cond)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("dotplot","_",ora_cond)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "dotplot" ,"text", ora_cond, sep = "_"))
    rv[[paste0("dotplot","_",ora_cond)]]
  }else{
    empty_image()
  }
})

output[[paste0("emapplot","_",ora_cond)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("emapplot","_",ora_cond)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "emapplot" ,"text", ora_cond, sep = "_"))
    rv[[paste0("emapplot","_",ora_cond)]]
  }else{
    empty_image()
  }
})


output[[paste0("cnetplot","_",ora_cond)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("cnetplot","_",ora_cond)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "cnetplot" ,"text", ora_cond, sep = "_"))
    rv[[paste0("cnetplot","_",ora_cond)]]
  }else{
    empty_image()
  }
})


output[[paste0("heatplot","_",ora_cond)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("heatplot","_",ora_cond)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "heatplot" ,"text", ora_cond, sep = "_"))
    rv[[paste0("heatplot","_",ora_cond)]]
  }else{
    empty_image()
  }
})


output[[paste0("treeplot","_",ora_cond)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("treeplot","_",ora_cond)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "treeplot" ,"text", ora_cond, sep = "_"))
    rv[[paste0("treeplot","_",ora_cond)]]
  }else{
    empty_image()
  }
})


output[[paste0("upsetplot","_",ora_cond)]] <- renderPlot({
 
  if (!is.null(rv[[paste0("upsetplot","_",ora_cond)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "upsetplot" ,"text", ora_cond, sep = "_"))
    rv[[paste0("upsetplot","_",ora_cond)]]
  }else{
    empty_image()
  }
})
