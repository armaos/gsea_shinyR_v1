

output[[paste0("volcano0_", cond)]] <- renderPlot({
  lfcT=1.5
  fdrT=0.05
  
  
  v0 <- volcano_plot_data_list[[cond]]
  
  ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
    geom_point(alpha=0.6, size=1) + 
    theme_bw(base_size=15) + 
    theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
    ggtitle(paste0("Volcanoplot of ", cond )) + 
    scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
    scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
    geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
})


output[[paste0("pca_", cond)]] <- renderPlot({
  if (!is.null(rv[[paste0("pca_", cond)]]) ){
    rv[[paste0("pca_", cond)]]
  }else{
    empty_image()
  }
})


output[[paste0("pheatmap_", cond)]] <- renderPlot({
  lfcT=1.5
  fdrT=0.05
  if (!is.null(rv[[paste0("pheatmap_", cond)]]) ){
    rv[[paste0("pheatmap_", cond)]]
  }else{
    empty_image()
  }
  # id2 <- heatmap_plot_data_list[[cond]]
  # pheatmap(id2,
  #          angle_col=90,
  #          fontsize=8,
  #          
  #          cutree_rows=2,
  #          cutree_col=2,
  #          scale="row",
  #          cluster_cols=F,
  #          treeheight_col=5,
  #          treeheight_row=15,
  #          color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  #          annotation=annotation,
  #          width=6,
  #          height=12,
  #          show_rownames=F,
  #          main = paste0("Expression of all DEGs ( with FDR < " ,fdrT ," and |logFC| > ", lfcT ," )"))
  # 
    
  })