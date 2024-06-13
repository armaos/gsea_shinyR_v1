
get_ora_vector <- function(ora_cond){
  
  print("get_ora_vector")
  print(ora_cond)
  
  
  if(rv$submission_by == "By group"){
    group = input$gene_select_by_group
    values = input$gene_select_by_group_values
    
    ora_input_data_genes = rv$msi_df %>% 
      select(primary_name, all_of(group)) %>%
      filter(!is.na(!!as.symbol(group))) %>% 
      separate_rows(!!as.symbol(group), sep = "[,]") %>%  
      filter(!!as.symbol(group) %in% values ) %>%
      unique() %>%
      #separate_rows(primary_name, sep = ",") %>%
      distinct(primary_name)
    
    
    ora_input_data <- rv$msi_df %>% 
      select(primary_name, contains("logFC"), contains("FDR")) %>%
      #separate_rows(Name, sep = ",") %>%
      filter(primary_name %in% ora_input_data_genes$primary_name ) %>%
      gather(k, v , -primary_name) %>% 
      separate(k, c("cond", "value_type"), sep = "[.]") %>% 
      spread(value_type, v)
    
    
  }else if(rv$submission_by == "Upon filtering"){
    
    fc_filter_low <<- input$slider_lfc[1] 
    fc_filter_high <<- input$slider_lfc[2] 
    
    ora_input_data = get_filtered_msi_df(rv$msi_df %>% 
                                       select(primary_name, contains("logFC"), contains("FDR")), ora_cond, fdr_filter, fc_filter_low )
      
    
  }else if(rv$submission_by == "Custom Set"){
    ora_input_data_genes <- rv$msi_df %>%  
      separate_rows(Name, sep = ",") %>%
      filter(Name %in% rv$input_genes) %>%
      distinct(primary_name)
    
    ora_input_data <- rv$msi_df %>% 
      select(primary_name, contains("logFC"), contains("FDR")) %>%
      #separate_rows(Name, sep = ",") %>%
      filter(primary_name %in% ora_input_data_genes$primary_name ) %>%
      gather(k, v , -primary_name) %>% 
      separate(k, c("cond", "value_type"), sep = "[.]") %>% 
      spread(value_type, v) 
    
  }
  
  
  
  if(ora_cond != ""){
    ora_input_data <- ora_input_data %>% 
      filter(cond == ora_cond)
  }
  return(ora_input_data)
  
}

rename_msi_df_columns = function(colnames_all, old_name, new_name){
  colnames_all[colnames_all %in% c(old_name)]= new_name
  return(colnames_all)
}

get_filtered_msi_df <- function(df, condition, fdr_filter, logfc_filter){
  print(paste(colnames(df)))
  return(df %>%
           gather(k, v , -primary_name) %>% 
           separate(k, c("cond", "value_type"), sep = "[.]") %>% 
           spread(value_type, v) %>% 
           filter(cond == condition) %>%
           filter(abs(logFC) >= logfc_filter & FDR <= fdr_filter)
  )
}



see_pathview <- function(..., save_image = FALSE)
{
  arguments <- list(...)
  pathview(gene.data = arguments$gene.data,
           gene.idtype = arguments$gene.idtype,
           pathway.id = arguments$pathway.id,
           species = arguments$species,
           bins = arguments$bins,
           limit = arguments$limit,
           low = arguments$low,
           mid = arguments$mid,
           high = arguments$high,
           na.col = arguments$na.col,
           kegg.native = arguments$kegg.native,
           kegg.dir = arguments$kegg.dir,
           same.layer = arguments$same.layer,
           sign.pos = arguments$sign.pos,
           new.signature =  arguments$new.signature ,
           out.suffix	= arguments$out.suffix
           
           
  )
  
  msg <- capture.output(pathview::pathview(...), type = "message")
  
  msg <- grep("image file", msg, value = T)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  file.remove(filename)
  img <- png::readPNG(filename)
  grid::grid.raster(img)
  
}

empty_image <- function(..., save_image = FALSE)
{
  print("empty image!")
  grid::grid.raster(png::readPNG("empty.png"))
}



get_manual_metadata <- function(){
  x <- reactiveValuesToList(input)
  x <- c(x[startsWith(names(x), "Subgroup_")],  x[startsWith(names(x), "Treatment_") ] )
  
  metadata_df = data.frame(
    names = names(x),
    values = unlist(x, use.names = FALSE)
  ) %>% 
    filter(grepl("Treatment_", names) | grepl("Subgroup_", names) ) %>%
    mutate(
      names = str_remove(names, "Treatment_"),
      names = str_remove(names, "Subgroup_")
    ) %>%
    group_by(names) %>%
    summarise(values = paste(values, collapse = "."))

  return(metadata_df)
}

rename_dol <- function(dol_df){
  dol_df <- dol_df %>% 
    rename_at(vars(ends_with('.O')), funs(str_replace(., ".O$", ".Medium"))) %>%
    rename_at(vars(ends_with('.D')), funs(str_replace(., ".D$", ".Low"))) %>%
    rename_at(vars(ends_with('.L')), funs(str_replace(., ".L$", ".High")))  %>%
    rename_at(vars(matches('[B][1-9]')), funs(str_replace(., "^B", "Brown")))  %>%
    rename_at(vars(matches('[G][1-9]')), funs(str_replace(., "^G", "Green")))
  return(dol_df)
}

rename_condition_dol <- function(cond){
  cond_label = switch (cond,
                       "D_O" = "Low_Medium",
                       "L_O" = "High_Medium", 
                       "D_L" = "Low_High" ,
                       "B_G_in_O" = "Brown_Green_in_Medium",
                       "B_G_in_L" = "Brown_Green_in_High",
                       "B_G_in_D" = "Brown_Green_in_Low",
  )
  return(cond_label)
}

clean_RunAcc <- function(runacc){
  return( str_remove(runacc, "_001.fastq.tabular|_002.fastq.tabular|_fastq.tabular|.fastq.tabular|.tabular"))
}


get_choices_for_contrasts <- function(colnames_filtered, subgroup, subgroup_class, treatment){
  
  # extract inter contrasts if treatments
  if(length(unique(treatment)) > 1){
    inter_choices = subset(as.data.frame(t(combn(levels(factor(treatment)), 2))) %>% 
                             mutate(contrast = paste0(V1, "-", V2)))$contrast
  }else{
    inter_choices = c()
  }
  
  
  
  
  n <- 1
  regex <- paste0("\\w+(?:-\\w+){0,", n, "}")
  # extract intra contrasts per treatment
  intra_choices = data.frame(cols =  colnames_filtered, 
                             subgroup = subgroup, 
                             treatment = treatment , 
                             subgroup_class = subgroup_class )  %>% 
    group_by(treatment) %>% 
    filter(n_distinct(subgroup_class) > 1) 
  
  if(dim(intra_choices)[1]>0){
    intra_choices <- intra_choices %>%
      summarise(intra_groups = paste(combn(unique(subgroup_class), 2), collapse = "-", sep = ",")) %>% 
      ungroup() %>%
      rowwise() %>%
      mutate(intra_groups =  paste(str_extract_all(intra_groups, regex, simplify = T), collapse= ",")) %>%
      separate_rows(intra_groups, sep = ",") %>%
      mutate(intra_contrasts = paste(intra_groups, treatment, sep = "_in_"))
      
    if(length(inter_choices) > 0){
      choices = rbind(data.frame(contrasts = inter_choices , type = "inter"),
                      data.frame(contrasts = intra_choices$intra_contrasts,  type = "intra"))
    }else{
      choices = data.frame(contrasts = intra_choices$intra_contrasts,  type = "intra")
    }
    
  }else{
    if(length(inter_choices) > 0){
      choices = data.frame(contrasts = inter_choices , type = "inter")
    }else{
      choices = c()
    }
              
  }
    
  return(choices)
}



make_overview_plots <- function(group, filtered, type){
  print("make_overview_plots")
  # print("colnames(filtered)")
  # print(colnames(filtered))
  
  #input
  design <- model.matrix(~0+group)             #when paired analysis required, this need to be ~0+group+donor
  rownames(design)<-colnames(filtered)
  colnames(design)<-levels(factor(group))
  
  
  ########################################################################
  html("overview_text", paste("Processing normalization..."))
  incProgress(0.1, detail = paste("Processing PCA..."))
  
  y <- calculate_y(filtered, group)
  # print(paste("Done y"))  
  # print("doing PCA")
  source(file = file.path(version, "server_PCA_plot.R"),
         local = TRUE,
         encoding = "UTF-8")
  rv$pca_plot = pca_plot
  # print("PCA processed")
  
  ########################################################################
  html("overview_text", paste("Processing Constrasts..."))
  incProgress(0.05, detail = paste("Processing Constrasts..."))
  
  ####
  

  
  
  ## overview Heatmap plot
  html("overview_text", "Processing overview Heat-Map plot...")
  incProgress(0.1, detail = paste("Processing overview Heat-Map plot..."))
  
  # id1 <- head(y, 1000)
  
  id1 <- y
  # this to put the columns in order in the heatmap plot
  
  
  selectedCols <- c()
  if(type == "inter"){
    for(g in unique(group)){
      selectedCols <-c(selectedCols, grep(paste0("\\.", g) ,colnames(id1)))
    }
  }else{
    for(g in unique(group)){
      selectedCols <-c(selectedCols, grep(paste0(g,"[1-9,a-z,A-Z]*\\.") ,colnames(id1)))
    }
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
  showTab(inputId = "norm_de_plot_tabsetpanel", target = "Heat-Map plot (overview)")
  # print("Done heatmap overview")

  print("make_overview_plots END")
}


make_fit_y <- function(single_replicates, y, design){
  
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
  
  return(list(y , fit))
}

calculate_y <- function(filtered, group){
  y <- DGEList(counts=filtered, group=group)         	# 
  keep <- rowSums(cpm(y) > 1 ) >= 2
  y <- y[keep,keep.lib.sizes=FALSE]    
  
  y <- calcNormFactors(y,method="none") #tximport norm before
  return(y)
}


get_group_subgroup_treatment <- function(colnames_filtered){

  if(input$choose_custom_tcmp == "testdata"){
    
    subgroup <- factor(gsub("\\.\\w+$","", colnames_filtered, perl=T))  # this is the subgroup
    subgroup_class <- factor(gsub("\\w\\.\\w+$","", colnames_filtered, perl=T))  # this is the subgroup without numbering
    treatment <- factor(gsub("^\\w+\\.","", colnames_filtered, perl=T))  # this is the treatment
  }else if(input$choose_custom_tcmp == "galaxy"){
    sst = data.frame(colnames = colnames_filtered) %>%
      separate(colnames, c("subgroup", "treatment"), sep = "[.]", remove = F) %>%
      mutate(subgroup_class = subgroup)
      
    treatment = sst$treatment
    subgroup = sst$subgroup
    subgroup_class = sst$subgroup_class
  }
  return(list(subgroup, subgroup_class, treatment))
}
  
  
  

  
update_notifications <- function(message_vector){
  
  nots <- apply(matrix(message_vector), 1,function(row) {
    notificationItem(text = row, 
                     icon("check"), 
                     status = "success")
  })
  
  output$notifications.type <- renderMenu(
    dropdownMenu(type = "notifications",  .list = nots)
  )
}

  
