
output$pheatmap <- renderPlot({
  if (!is.null(rv$pheatmap) ){
    rv$pheatmap
  }else{
    empty_image()
  }
})

output$volcano0 <- renderPlot({
  if (!is.null(rv$volcano0) ){
    rv$volcano0
  }else{
    empty_image()
  }
})

output$pca_plot <- renderPlot({
  if (!is.null(rv$pca_plot) ){
    rv$pca_plot
  }else{
    empty_image()
  }
})
