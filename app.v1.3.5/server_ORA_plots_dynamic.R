output[[paste0("barplot1","_",id)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("barplot1","_",id)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "barplot1" ,"text", id, sep = "_"))
    rv[[paste0("barplot1","_",id)]]
  }else{
    empty_image()
  }
})


output[[paste0("barplot2","_",id)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("barplot2","_",id)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "barplot2" ,"text", id, sep = "_"))
    rv[[paste0("barplot2","_",id)]]
  }else{
    empty_image()
  }
})


output[[paste0("dotplot","_",id)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("dotplot","_",id)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "dotplot" ,"text", id, sep = "_"))
    rv[[paste0("dotplot","_",id)]]
  }else{
    empty_image()
  }
})

output[[paste0("emapplot","_",id)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("emapplot","_",id)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "emapplot" ,"text", id, sep = "_"))
    rv[[paste0("emapplot","_",id)]]
  }else{
    empty_image()
  }
})


output[[paste0("cnetplot","_",id)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("cnetplot","_",id)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "cnetplot" ,"text", id, sep = "_"))
    rv[[paste0("cnetplot","_",id)]]
  }else{
    empty_image()
  }
})


output[[paste0("heatplot","_",id)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("heatplot","_",id)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "heatplot" ,"text", id, sep = "_"))
    rv[[paste0("heatplot","_",id)]]
  }else{
    empty_image()
  }
})


output[[paste0("treeplot","_",id)]] <- renderPlot({
  
  if (!is.null(rv[[paste0("treeplot","_",id)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "treeplot" ,"text", id, sep = "_"))
    rv[[paste0("treeplot","_",id)]]
  }else{
    empty_image()
  }
})


output[[paste0("upsetplot","_",id)]] <- renderPlot({
 
  if (!is.null(rv[[paste0("upsetplot","_",id)]]) ){
    shinyjs::hide("text_ora_plots_dynamic")
    shinyjs::enable("run")
    shinyjs::hide(paste("ora", "upsetplot" ,"text", id, sep = "_"))
    rv[[paste0("upsetplot","_",id)]]
  }else{
    empty_image()
  }
})
