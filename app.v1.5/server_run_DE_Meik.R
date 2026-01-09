# shinyjs::show("overview_text")
withProgress(message = 'Running DE analysis', value = 0, {
  print("starting server_run_DE_Meik.R")  
  
  filtered = rv$filtered
 
  # in case there are other conditions in the data just select those that are required
  select_filtered_columns = data.frame(cond = rv$condition) %>% mutate(cond = str_split(cond, "_in_", simplif =T)[,1]) 
  select_filtered_columns %>% mutate(cond = str_split(cond, "_in_", simplif =T)[,1])
  filtered = as.data.frame(filtered)  %>% 
    select(contains(unique(unlist(strsplit(select_filtered_columns$cond, "-")))))
  
  update_notifications(rv$condition)
  #shinyjs::html("pageHeader", paste(rv$condition, collapse = ", "))
 
  
  
  sst = get_group_subgroup_treatment(colnames(filtered))
  subgroup = sst[[1]]
  subgroup_class = sst[[2]]
  treatment = sst[[3]]
  # print( treatment)
  # print( subgroup_class)
  # print( subgroup)
  
  choices = get_choices_for_contrasts(colnames(filtered), subgroup, subgroup_class, treatment)
  
  contrast_types_submitted = unique(subset(choices %>% filter(contrasts %in% rv$condition))$type)
  # print(contrast_types_submitted)
  
  
  lfcT=1.5
  fdrT=0.05
  
  # make overview plots by treatment
  if(length(unique(treatment)) > 1){
    make_overview_plots(treatment, filtered, "inter")
  }else if(length(unique(subgroup_class)) > 1 ){
    make_overview_plots(subgroup_class, filtered, "intra")
  }
  
  
  showTab(inputId = "norm_de_plot_tabsetpanel", target = "Volcano plot (pairwise)")
  showTab(inputId = "norm_de_plot_tabsetpanel", target = "Heat-Map plot (pairwise)")
  showTab(inputId = "norm_de_plot_tabsetpanel", target = "PCA plot (pairwise)")
  
  
  
  for(type_constrast in c("inter", "intra")){
    if(type_constrast %in% contrast_types_submitted){
      # print(paste("doing type_constrast", type_constrast))  

      incProgress(0.1, detail = paste("Fitting glmLRT..."))
      # print(paste("Ready to start per condition"))  
     
      conditions_in_to_run = choices %>% 
        filter(contrasts %in% rv$condition & type == type_constrast) 
       
      for(cond in conditions_in_to_run$contrasts ){
        warning_message = ""
        setProgress(detail = paste("Processing Constrast:", cond))
        source(file = file.path(version, "server_run_condition_contrast.R"),
               local = TRUE,
               encoding = "UTF-8") 
        
      }
    }
  }


 
    
  setProgress(0.9, detail = paste("Merging annotations"))
  
  names(volcano_plot_data_list) <- rv$condition
  #this to happen only when i load new. If on the precompiled , there is no need to do
  if(create_msi_df){
    print("server_run_DE:")
    # print(paste( colnames(rv$msi_df), collapse = ", "))
    rv$update_annotations =  rv$update_annotations + 1
    rv$update_msi_df =  rv$update_msi_df + 1
    rv$update_msi_clustering =  rv$update_msi_clustering + 1
  }
  shinyjs::show("menu_norm_and_de_panel")
  # print("Done proc DE")
  html("overview_text", paste("DE analysis: DONE...", "</br>","</br>", warning_message)) 
  setProgress(1, detail = paste("DE analysis: DONE.."))
  
  rv$show_norm_and_de_mout = T  
  
})