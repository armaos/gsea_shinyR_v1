### submitting the username
observeEvent(input$submit_username, {
  
  username = input$username
  #run_id = "example"
  #count_file = "tcpm.csv"
  
  s3 <- paws::s3(config = list(region = region))
  
  #get file names
  l = s3$list_objects_v2(
    Bucket = bucket_name, 
    Prefix = paste0(username, "/")
    #Prefix = "",
    #Delimiter = "zip"
  )$Contents
  
  if(length(l) > 0){
    # get data frame with all the files of the user
    rv$runs_galaxy_df <- data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))
    # colnames(rv$runs_galaxy_df) = names(l)
    rv$runs_galaxy_df <- rv$runs_galaxy_df %>% 
      select(X1) %>% 
      rename(Key = X1)  %>% 
      mutate(path = Key)  %>% 
      separate(Key, c("user", "run", "file"), "/") %>%
      group_by( run) %>%
      filter(grepl("output.zip", paste(file, collapse = ","))) %>%
      ungroup
    
    shinyjs::show("galaxy_submissions_box")
    #js$collapse("account_box")
    updateSelectInput(session, "galaxy_runs", choices = unique(rv$runs_galaxy_df$run) , selected = ifelse(unique(rv$runs_galaxy_df$run)[1], ""))

  }else{
    showNotification(paste("No User with that username"), duration = 10, type = "error")
  }
})



### check for submitted run
observeEvent(input$galaxy_runs, {
  if(is.null(input$galaxy_runs)){
    shinyjs::hide("submit_galaxy_runs")
  }else if(input$galaxy_runs == "" ){
    shinyjs::hide("submit_galaxy_runs")
  }else{
    shinyjs::show("submit_galaxy_runs")
  }
  
})



### submitting the different runs of the usernmae
observeEvent(input$submit_galaxy_runs, {
  rv$aws_files <- c()
  rv$galaxy_outputs <- data.frame()
  rv$metadata_df  <- data.frame()
  source(file = file.path(version, "server_download_counts.R"),
         local = TRUE,
         encoding = "UTF-8")
})

### submitting the selected output of the usernmae
observeEvent(input$submit_galaxy_outputs, {
  
  shinyjs::disable("submit_galaxy_outputs")
  shinyjs::disable("submit_metadata")
  shinyjs::hide("menu_de_box")
  shinyjs::hide("metadata_table_box")
  
  rv$metadata_df <- get_manual_metadata()


  # rv$galaxy_outputs = data.frame(file = c("SRR7637235.tabular" ,"SRR1207056.tabular","SRR7637208.tabular" , "SRR7637249.tabular" , "SRR7637250.tabular", "SRR7637252.tabular", "SRR7637253.tabular"),
  #                                path = c("tmp/SRR7637235.tabular", "tmp/SRR1207056.tabular", "tmp/SRR7637208.tabular", "tmp/SRR7637249.tabular" , "tmp/SRR7637250.tabular","tmp/SRR7637252.tabular","tmp/SRR7637253.tabular"))

 
  fastqc_runs <- rv$galaxy_outputs$file
  source(file = file.path(version, "server_fastqc_reports.R"),
         local = TRUE,
         encoding = "UTF-8")
  
  
  print("calling server_combine_galaxy_output.R")
  source(file = file.path(version, "server_combine_galaxy_output.R"),
         local = TRUE,
         encoding = "UTF-8")

  shinyjs::enable("submit_galaxy_outputs")
  shinyjs::enable("submit_metadata")
})



observeEvent(input$galaxy_outputs, 
             ignoreNULL = FALSE, {

  if(is.null(input$galaxy_outputs)){
    rv$aws_files <- c()
    shinyjs::disable("submit_metadata")
    #shinyjs::disable("submit_galaxy_outputs")
  }else if(length(input$galaxy_outputs) > 1 ){
    rv$aws_files <- input$galaxy_outputs
    shinyjs::enable("submit_metadata")
    #shinyjs::enable("submit_galaxy_outputs")
  }else if(length(input$galaxy_outputs) > 0 ){
    rv$aws_files <- input$galaxy_outputs
    shinyjs::disable("submit_metadata")
    #shinyjs::disable("submit_galaxy_outputs")
  }
})


texts <- reactive({
  
  files <- rv$aws_files
  n_files = length(files)
  if (n_files > 0) {
    lapply(seq_len(n_files), function(i) {
      p(paste0( files[i], ": ")) %>% 
        tagAppendAttributes(class = 'inline_text')
      
    })
  }
  
})

textboxes_donor <- reactive({
  files <- rv$aws_files
  
  n_files = length(files)
  if (n_files > 0) {
    metadata <-  data.frame(file =  rv$aws_files) %>%
      mutate(Run = clean_RunAcc(file)) %>%
      inner_join(precomp_metadata, by = "Run")
    
    
    
    lapply(seq_len(n_files), function(i) {
      my_value = subset(metadata %>% filter(file == files[i]))$SampleName
      my_value = ifelse(length(my_value) == 0, paste0("Subgroup",i), my_value)
      my_value = str_replace(my_value, "$ ", "")
      my_value = str_replace_all(my_value, " ", "_")
      my_value = str_replace_all(my_value, "[.]", "_")
      
      
      textInput(inputId = paste0("Subgroup_", files[i]),
                label = NULL,
                value = my_value,
                width = '100px') %>% 
        tagAppendAttributes(class = 'inline_textinput')
    })
  }
  
})


textboxes_group <- reactive({
  
  files <- rv$aws_files
  n_files = length(files)
  if (n_files > 0) {
    metadata <-  data.frame(file =  rv$aws_files) %>%
      mutate(Run = clean_RunAcc(file)) %>%
      inner_join(precomp_metadata, by = "Run")
    
    lapply(seq_len(n_files), function(i) {
      my_value = subset(metadata %>% filter(file == files[i]))$BioProject
      my_value = ifelse(length(my_value) == 0 , paste0("Treatment",i), my_value)
      my_value = str_replace(my_value, "$ ", "")
      my_value = str_replace_all(my_value, " ", "_")
      my_value = str_replace_all(my_value, "[.]", "_")
      
      textInput(inputId = paste0("Treatment_", files[i]),
                label = NULL,
                value = my_value,
                width = '100px') %>% 
        tagAppendAttributes(class = 'inline_textinput')
    })
  }
  
})

output$textbox_ui <- renderUI({ texts() })
output$textbox_donor_ui <- renderUI({ textboxes_donor() })
output$textbox_group_ui <- renderUI({ textboxes_group() })
