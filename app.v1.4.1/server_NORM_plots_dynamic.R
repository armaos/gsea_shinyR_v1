output[[paste0("pheatmap","_",id)]] <- renderPlot({
  if (!is.null(rv[[paste0("pheatmap","_",id)]]) ){
    rv[[paste0("pheatmap","_",id)]]
  }else{
    empty_image()
  }
})


output[[paste0("volcano0","_",id)]] <- renderPlot({
  if (!is.null(rv[[paste0("volcano0","_",id)]]) ){
    rv[[paste0("volcano0","_",id)]]
  }else{
    empty_image()
  }
})

