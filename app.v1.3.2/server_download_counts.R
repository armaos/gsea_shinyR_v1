# download the folder requested
s3 <- paws::s3(config = list(region = region))


for(run_id in input$galaxy_runs){
  requested_path = subset(rv$runs_galaxy_df %>% filter(run == run_id))$path
  s3_download <- s3$get_object(
    Bucket = bucket_name,
    Key = requested_path,
  )
  download_zip <- file.path(tempdir(), paste(input$username, run_id, "output.zip", sep = "_"))
  download_zip_folder <- str_remove_all(download_zip, ".zip")
  writeBin(s3_download$Body, con = download_zip)
  unzip(download_zip, exdir = download_zip_folder)
  rv$galaxy_outputs = rbind(rv$galaxy_outputs, data.frame(file = list.files(file.path(download_zip_folder, "pipeline_output/")),
                                                    path = list.files(file.path(download_zip_folder, "pipeline_output/"), full.names=T)
                                                    ))
  file.remove(download_zip)
}

shinyjs::show("galaxy_outputs_box")
#js$collapse("galaxy_runs_box")
#js$collapse("account_box")
updateSelectInput(session, "galaxy_outputs", choices = rv$galaxy_outputs$file)
