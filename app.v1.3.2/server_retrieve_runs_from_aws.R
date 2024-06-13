username = input$username
#run_id = "example"
#count_file = "tcpm.csv"

s3 <- paws::s3(config = list(region = region))

#get file names
l = unlist(s3$list_objects_v2(
  Bucket = bucket_name, 
  Prefix = paste0(username, "/")
  #Prefix = "",
  #Delimiter = "zip"
)$Contents)

# get data frame with all the files of the user
rv$runs_galaxy_df <- data.frame(matrix(unlist(l), ncol=length(l), byrow=TRUE, ))
colnames(rv$runs_galaxy_df) = names(l)
rv$runs_galaxy_df <- rv$runs_galaxy_df %>% select(Key) %>% mutate(path = Key)  %>% separate(Key, c("user", "run", "file"), "/")

shinyjs::show("galaxy_submissions_box")
#js$collapse("account_box")
updateSelectInput(session, "galaxy_runs", choices = rv$runs_galaxy_df$run , selected = ifelse(rv$runs_galaxy_df$run[1], ""))


