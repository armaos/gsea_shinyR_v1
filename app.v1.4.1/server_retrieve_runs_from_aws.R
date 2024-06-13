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
  
  # get data frame with all the files of the user
  rv$runs_galaxy_df <- data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE))
  # colnames(rv$runs_galaxy_df) = names(l)
  rv$runs_galaxy_df <- rv$runs_galaxy_df %>% select(X1) %>% rename(Key = X1)  %>% mutate(path = Key)  %>% separate(Key, c("user", "run", "file"), "/")
  
  shinyjs::show("galaxy_submissions_box")
  #js$collapse("account_box")
  updateSelectInput(session, "galaxy_runs", choices = rv$runs_galaxy_df$run , selected = ifelse(rv$runs_galaxy_df$run[1], ""))

  
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
  source(file = file.path(version, "server_download_counts.R"),
         local = TRUE,
         encoding = "UTF-8")
})

### submitting the selected output of the usernmae
observeEvent(input$submit_galaxy_outputs, {
  shinyjs::hide("menu_de_box")
  
  print("submitted galaxy outputs")
  x <- reactiveValuesToList(input)
  x <- c(x[startsWith(names(x), "donor_")],  x[startsWith(names(x), "group_") ] )

  rv$custom_conditions_names = data.frame(
    names = names(x),
    values = unlist(x, use.names = FALSE)
  ) %>% 
    filter(grepl("group_", names) | grepl("donor_", names) ) %>%
    mutate(
      names = str_remove(names, "group_"),
      names = str_remove(names, "donor_")
    ) %>%
    group_by(names) %>%
    summarise(values = paste(values, collapse = "."))
  
  #print(dim(rv$custom_conditions_names)[1])
  
  #print(rv$custom_conditions_names$names)
  
  #print(rv$custom_conditions_names$values)
  
  print("calling server_combine_galaxy_output.R")
  source(file = file.path(version, "server_combine_galaxy_output.R"),
         local = TRUE,
         encoding = "UTF-8")
})



observeEvent(input$galaxy_outputs, {
  rv$aws_files <- input$galaxy_outputs
  
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
    lapply(seq_len(n_files), function(i) {
      textInput(inputId = paste0("donor_", files[i]),
                label = NULL,
                value = paste0("Subgroup",i),
                width = '100px') %>% 
        tagAppendAttributes(class = 'inline_textinput')
    })
  }
  
})


textboxes_group <- reactive({
  
  files <- rv$aws_files
  n_files = length(files)
  if (n_files > 0) {
    lapply(seq_len(n_files), function(i) {
      textInput(inputId = paste0("group_", files[i]),
                label = NULL,
                value = paste0("Treatment",i),
                width = '100px') %>% 
        tagAppendAttributes(class = 'inline_textinput')
    })
  }
  
})

output$textbox_ui <- renderUI({ texts() })
output$textbox_donor_ui <- renderUI({ textboxes_donor() })
output$textbox_group_ui <- renderUI({ textboxes_group() })
