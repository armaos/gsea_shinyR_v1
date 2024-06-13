# download the folder requested
s3 <- paws::s3(config = list(region = region))
###
for(run_id in input$galaxy_runs){
  #run_id = "example_with_metadata_and_fastqc"
  requested_paths = unique(subset(rv$runs_galaxy_df %>% filter(run == run_id)))
  for(i in 1:dim(requested_paths)[1]){
    
    requested_path = requested_paths[i,]$path
    requested_file = requested_paths[i,]$file
    print(paste(run_id,requested_file))
    
    s3_download <- s3$get_object(
      Bucket = bucket_name,
      Key = requested_path,
    )
    #download_zip <- file.path(tempdir(), paste(input$username, run_id, "output.zip", sep = "_"))
    #download_zip_folder <- str_remove_all(download_zip, ".zip")
    download_zip <- file.path(tempdir(), paste(input$username, run_id, requested_file, sep = "_"))
    download_zip_folder  <- file.path(tempdir(), paste(input$username, run_id, sep = "_"))
    writeBin(s3_download$Body, con = download_zip)
    if(grepl(".zip", requested_file)){
      unzip(download_zip, exdir = download_zip_folder)  
      file.remove(download_zip)
    }else{
      file.rename(from = download_zip,  to = file.path(download_zip_folder, requested_file))
    }
    
    
  }
  
  rv$galaxy_outputs = rbind(rv$galaxy_outputs, 
                         data.frame(file = list.files(file.path(download_zip_folder, "pipeline_output/")),
                                    path = list.files(file.path(download_zip_folder, "pipeline_output/"), full.names=T),
                                    run_id = run_id) %>%
                           mutate(fastq_path = file.path(download_zip_folder, "fastqc_output", str_remove(file, ".tabular") ))
  ) 
  
}


###
# for(run_id in input$galaxy_runs){
#   requested_path = subset(rv$runs_galaxy_df %>% filter(run == run_id))$path
#   s3_download <- s3$get_object(
#     Bucket = bucket_name,
#     Key = requested_path,
#   )
#   download_zip <- file.path(tempdir(), paste(input$username, run_id, "output.zip", sep = "_"))
#   download_zip_folder <- str_remove_all(download_zip, ".zip")
#   writeBin(s3_download$Body, con = download_zip)
#   unzip(download_zip, exdir = download_zip_folder)
#   rv$galaxy_outputs = rbind(rv$galaxy_outputs, data.frame(file = list.files(file.path(download_zip_folder, "pipeline_output/")),
#                                                           path = list.files(file.path(download_zip_folder, "pipeline_output/"), full.names=T)
#   ))
#   file.remove(download_zip)

# }



shinyjs::show("galaxy_outputs_box")
#js$collapse("galaxy_runs_box")
#js$collapse("account_box")

rv$galaxy_outputs <- rv$galaxy_outputs  %>%
  group_by(file) %>%
  mutate(r = 1:n(), n=n_distinct(run_id)) %>% 
  ungroup() %>%
  mutate(file2 = ifelse(n>1, paste(file, r,sep = " - " ), file)) 

g_list <- list()
for(i in unique(rv$galaxy_outputs$run_id)){
  lnew = as.list(subset(rv$galaxy_outputs %>% filter(run_id == i))$file2)
  lnew = list(lnew)
  names(lnew) = i
  g_list = append(g_list, lnew)
}
galaxy_outputs_choices_list <<- g_list
#rv$galaxy_outputs = data.frame(path = list.files("~/pipeline_output/", full.names = T), file = list.files("~/pipeline_output/")) 
#updateSelectInput(session, "galaxy_outputs", choices = rv$galaxy_outputs$file)
updateSelectInput(session, "galaxy_outputs", choices = galaxy_outputs_choices_list)
