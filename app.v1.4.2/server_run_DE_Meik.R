# shinyjs::show("overview_text")
withProgress(message = 'Running DE analysis', value = 0, {
  
  
  filtered = rv$filtered
  # in case there are other conditions in the data just select those that are required
  filtered = as.data.frame(filtered)  %>% 
    select(contains(unique(unlist(strsplit(rv$condition, "-")))))
  
  shinyjs::html("pageHeader", paste(rv$condition, collapse = ", "))
  
  
  
  donor <- factor(gsub("\\.\\w+$","",colnames(filtered),perl=T))
  group <- factor(gsub("^\\w+\\.","",colnames(filtered),perl=T))
  
  x<-group
  #input
  design <- model.matrix(~0+group)             #when paired analysis required, this need to be ~0+group+donor
  rownames(design)<-colnames(filtered)
  colnames(design)<-levels(factor(group))
  
  ########################################################################
  html("overview_text", "Processing normalization...")
  incProgress(0.1, detail = paste("Processing PCA..."))
 
  source(file = file.path(version, "server_PCA_plot.R"),
         local = TRUE,
         encoding = "UTF-8")
  print("PCA processed")
  ########################################################################
  html("overview_text", "Processing Constrasts...")
  incProgress(0.1, detail = paste("Processing Constrasts..."))
  
  ####
  showTab(inputId = "norm_de_plot_tabsetpanel", target = "Volcano plot")
  showTab(inputId = "norm_de_plot_tabsetpanel", target = "Heat-Map plot")
  
  
  
  # count replicates per group
  single_replicates = as.data.frame(design) %>% 
    gather(k, v) %>% 
    group_by(k) %>% 
    summarise(n = sum(v==1)) %>% 
    filter(n == 1) 
  
  
  ## overview Heatmap plot
  html("overview_text", "Processing overview Heat-Map plot...")
  incProgress(0.1, detail = paste("Processing overview Heat-Map plot..."))
  
  #id1 <- head(y, 1000)
  id1 <- y
  # this to put the column sin order in the heatmap plot
  selectedCols <- c()
  for(g in unique(group)){
    selectedCols <-c(selectedCols, grep(paste0("\\.", g) ,colnames(id1)))
  }
  
  id2 <- id1[,selectedCols] 
  annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
  rownames(annotation) <- colnames(id2)              # check out the row names of annotation
  rv$heatmap_overview_plot <- pheatmap(id2,
                                       fontsize=8, 
                                       angle_col=90, 
                                       cutree_rows=2, 
                                       cutree_col=2, scale="row", 
                                       cluster_cols=F, 
                                       treeheight_col=5,
                                       treeheight_row=15, 
                                       color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
                                       annotation=annotation,
                                       width=6, height=12, 
                                       show_rownames=F)
  
  
  incProgress(0.1, detail = paste("Fitting glmFit..."))
  # if single_replicates >0 , this means that under the group we have only one replicate
  if(dim(single_replicates)[1] > 1){
    y$common.dispersion <- 0.1
    fit <- glmFit(y, design, robust=TRUE)
    warning_message = "WARNING: to few replicates per treatment/subgroup."
    html("overview_text", warning_message)
    html("overview_text", paste("Processing Constrasts: Fitting glmFit...", "</br>","</br>", warning_message))
  }else{
    y.Disp <- estimateDisp(y, design, robust = TRUE)
    warning_message = ""
    html("overview_text", "Processing Constrasts: Fitting glmFit...")
    fit <- glmFit(y.Disp, design, robust=TRUE)
  }
  
  
  
  
  
  
  
  
  lfcT=1.5
  fdrT=0.05
  
  html("overview_text", paste("Processing Constrasts: Fitting glmLRT...", "</br>","</br>", warning_message))
  incProgress(0.1, detail = paste("Fitting glmLRT..."))
  
  for(cond in rv$condition ){
    setProgress(detail = paste("Processing Constrast:", cond))
    
    #cond = str_replace(cond, "_", "-")
    
    cmd <- paste("my.contrasts <- makeContrasts(", cond, ", levels = design)", sep ='"')
    #cmd <- paste("my.contrasts <- makeContrasts(", "G1-G2", ", levels = design)", sep ='"')
    #my.contrasts = makeContrasts(G1_G2="G1-G2", levels=design)              # need to rewrite acccording to the task
    
    eval(parse(text = cmd))
    
    
    #my.contrasts = makeContrasts(cond, levels=design)              # need to rewrite acccording to the task
    print("makeContrasts OK")
  
    lrt.D_O <- glmLRT(fit, contrast=my.contrasts)
    # lrt.aDTP_WT <- glmQLFTest(fit, contrast=my.contrasts[,id])  # quasi-likelihood (QL)
    
    # add FDR
    deg.edger <- lrt.D_O$table[p.adjust(lrt.D_O$table$PValue, method = "BH") < fdrT, ]
    #dim(deg.edger)
    
    #################################################################
    # final volcano plot           ##################################
    print("Proc volcano Data")
    html("overview_text", paste("Processing Constrasts: Processing volcano data...", "</br>","</br>", warning_message))
    setProgress(detail = paste("Processing Constrast:", cond, ":`volcano plot"))
    
    v0<-lrt.D_O$table
    
    v0$FDR <-p.adjust(v0$PValue, method = "BH") 
    deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1
    
    
    deg.edger1 <- deg.edger0[order(abs(deg.edger0$logFC)*(-log10(deg.edger0$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger0)>150,150,nrow(deg.edger0))],]  #select top150
    degnames1<-rownames(deg.edger1)
    degnames1<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames1,invert=T,value=T)
    #nrow(deg.edger1)
    v0$genelabels<-NA
    v0[rownames(v0) %in% degnames1,]$genelabels<-T  
    options(ggrepel.max.overlaps = Inf)                            # show only DEGs
    v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
   
    print("Proc Volcano PLot")
    html("overview_text", paste("Processing Constrasts: Processing volcano plot...", "</br>","</br>", warning_message)) 
    volcano_plot_data_list[[cond]] <<- v0
    
    # output[[paste0("volcano0_", cond)]] <- renderPlot({
    #   ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
    #   geom_point(alpha=0.6, size=1) + 
    #   theme_bw(base_size=15) + 
    #   theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
    #   ggtitle(paste0("Volcanoplot of ", cond )) + 
    #   scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
    #   scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
    #   geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
    # })
    print(paste0("Volcanoplot of ", cond ))
    #output[[paste0("volcano0_", cond)]] <- renderPlot({rv[[paste0("volcano0_", cond)]] })
    #################################################################
    
    html("overview_text", paste("Processing Constrasts: Processing heatmap plot...", "</br>","</br>", warning_message)) 
    setProgress(detail = paste("Processing Constrast:", cond, ":`heatmap plot"))
    print("Heatmap plot")
    id1 <- log2(cpm(y)+1)
    id1 <- y
    group_1 = str_split(cond, "-")[[1]][1]
    group_2 = str_split(cond, "-")[[1]][2]
    print(group_1)
    print(group_2)
    selectedCols <- c(grep(paste0("\\.",group_1),colnames(y)), grep(paste0("\\.",group_2),colnames(y)))  #select some columns
    print(selectedCols)
    id2 <- id1[rownames(deg.edger0),selectedCols]
    annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
    rownames(annotation) <- colnames(id2)              # check out the row names of annotation
    
    #pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
    
    #pheatmap(id2,fontsize=3,angle_col=0,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_Heatmap.png"), annotation=annotation, width=6,height=12, show_rownames=F)
    
  
    
    heatmap_plot_data_list[[cond]] <- id2
    rv[[paste0("pheatmap_", cond)]] <- pheatmap(id2,
                                                angle_col=90,
                                                fontsize=8,
                                                
                                                cutree_rows=2,
                                                cutree_col=2,
                                                scale="row",
                                                cluster_cols=F,
                                                treeheight_col=5,
                                                treeheight_row=15,
                                                color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                                                annotation=annotation,
                                                width=6,
                                                height=12,
                                                show_rownames=F,
                                                main = paste0("Expression of all DEGs ( with FDR < " ,fdrT ," and |logFC| > ", lfcT ," ) for each sample"))
               

    
    #output[[paste0("pheatmap_", cond)]] <- renderPlot({rv[[paste0("pheatmap_", cond)]] })
    #rv$pheatmap <- pheatmap(id2,fontsize=8, angle_col=0, cutree_rows=2, cutree_col=2, scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), annotation=annotation, width=6, height=12, show_rownames=F)
    
    v0 <- v0 %>%
      mutate(Name = rownames(v0)) %>%
      select(Name, everything())  %>% 
      select(-c(genelabels, change))
    rownames(v0) <- NULL
    
    counts = as.data.frame(y$counts)
    counts <- counts %>% mutate(Name = rownames(counts))
    rownames(counts) <- NULL
    
    v0 <- v0 %>% 
      left_join(counts)
    
    #this to happen only when i load new. If on the precompiled , there is no need to do
    if(create_msi_df){
      for(col in c("logCPM", "logFC", "LR", "PValue", "PValue", "FDR", "genelabels", "change")){
        colnames(v0) = rename_msi_df_columns(colnames(v0), col, paste(cond, col, sep = "."))
      }
      if(dim(rv$msi_df)[1] == 0){
        rv$msi_df = v0
      }else{
        rv$msi_df <- rv$msi_df %>%
          full_join(v0)
      }
      
    }
  }
    
  setProgress(0.9, detail = paste("Merging annotations"))
  
  names(volcano_plot_data_list) <- rv$condition
  #this to happen only when i load new. If on the precompiled , there is no need to do
  if(create_msi_df){
    print("server_run_DE:")
    print(paste( colnames(rv$msi_df), collapse = ", "))
    rv$update_annotations =  rv$update_annotations + 1
    rv$update_msi_df =  rv$update_msi_df + 1
    rv$update_msi_clustering =  rv$update_msi_clustering + 1
  }
  shinyjs::show("menu_norm_and_de_panel")
  print("Done proc DE")
  html("overview_text", paste("DE analysis: DONE...", "</br>","</br>", warning_message)) 
  setProgress(1, detail = paste("DE analysis: DONE.."))
  
  rv$show_norm_and_de_mout = T  
  
})