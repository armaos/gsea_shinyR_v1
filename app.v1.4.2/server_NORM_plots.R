output$pca_plot <- renderPlot({
  if (!is.null(rv$pca_plot) ){
    rv$pca_plot
  }else{
    empty_image()
  }
})


output$heatmap_overview_plot <- renderPlot({
  if (!is.null(rv$heatmap_overview_plot) ){
    rv$heatmap_overview_plot
  }else{
    empty_image()
  }
})


output$all_cluster_heatmap <- renderPlot({
  if (!is.null(rv$all_cluster_heatmap) ){
    rv$all_cluster_heatmap
  }else{
    empty_image()
  }
})