#install.packages("paws")
#install.packages("aws.s3")
library("aws.s3")
#library(paws)



# Create an S3 client, then list the objects within my bucket `my-bucket`.
b = "forjazul.galaxy"
obj <- get_bucket(b)

get_bucket("s3://forjazul.galaxy")
my_Bucket <- get_bucket(b)
files <- get_bucket(my_Bucket, prefix="demo/input/test")
fileList <- list.files(path = "files", pattern=".csv", all.files = TRUE, recursive=FALSE)
fileList

obj <-aws.s3::get_object("s3://forjazul.galaxy/alex/my_history/tcpm.csv")  

get_bucket("forjazul.galaxy", prefix="alex")

data <- 
  aws.s3::s3read_using(read.csv, object = "s3://forjazul.galaxy/alex/my_history/tcpm.csv")


data <- 
  aws.s3::s3read_using(fread, object = "s3://forjazul.galaxy/alex/my_history/tcpm.csv") 

csvcharobj <- rawToChar(obj)  
con <- textConnection(csvcharobj)  
data <- read.csv(file = con)



#### this is what is the current that i want
bucket_name= "forjazul.galaxy"
username = "Marton"
run_id = "example"
#count_file = "tcpm.csv"
region = "us-east-1"

s3 <- paws::s3(config = list(region = region))

#s3$s3_list_objects(bucket_name)


#get file names
l = unlist(s3$list_objects_v2(
  Bucket = bucket_name, 
  Prefix = paste0(username, "/")
  #Prefix = "",
  #Delimiter = "zip"
)$Contents)

# get data frame with all the files of the user
df <- data.frame(matrix(unlist(l), ncol=length(l), byrow=TRUE, ))
colnames(df) = names(l)
df <- df %>% select(Key) %>% mutate(path = Key)  %>% separate(Key, c("user", "run", "file"), "/")

###
galaxy_outputs = data.frame()
for(run_id in c("example_with_metadata_and_fastqc", "metadata_test", "full_example")){
  #run_id = "example_with_metadata_and_fastqc"
  requested_paths = unique(subset(runs_galaxy_df %>% filter(run == run_id)))
  for(i in 1:dim(requested_paths)[1]){
    
    requested_path = requested_paths[i,]$path
    requested_file = requested_paths[i,]$file
    print(paste(run_id,requested_file))
    
    s3_download <- s3$get_object(
      Bucket = bucket_name,
      Key = requested_path,
    )
    download_zip <- file.path(tempdir(), paste("alex1", run_id, requested_file, sep = "_"))
    download_zip_folder  <- file.path(tempdir(), paste("alex1", run_id, sep = "_"))
    writeBin(s3_download$Body, con = download_zip)
    if(grepl(".zip", requested_file)){
      unzip(download_zip, exdir = download_zip_folder)  
      file.remove(download_zip)
    }else{
      file.rename(from = download_zip,  to = file.path(download_zip_folder, requested_file))
    }
  }
  
  galaxy_outputs = rbind(galaxy_outputs, 
                            data.frame(file = list.files(file.path(download_zip_folder, "pipeline_output/")),
                                       path = list.files(file.path(download_zip_folder, "pipeline_output/"), full.names=T)) %>%
                              mutate(fastq_path = file.path(download_zip_folder, "fastqc_output", str_remove(file, ".tabular") ))
                            ) 

}

###

# download the folder requested
requested_path = subset(df %>% filter(run == run_id))$path
s3_download <- s3$get_object(
  Bucket = bucket_name,
  Key = requested_path,
)
download_zip <- file.path(tempdir(), paste(username, run_id, "output.zip", sep = "_"))
download_zip_folder <- str_remove_all(download_zip, ".zip")
writeBin(s3_download$Body, con = download_zip)
unzip(download_zip, exdir = download_zip_folder)
galaxy_outputs = list.files(file.path(download_zip_folder, "pipeline_output/"))
file.remove(download_zip)
      
# get_bucket_df(
#   bucket = "s3://forjazul.galaxy/", 
#   region = region
#   )
# 
# save_object(
#   object = "dwca-nsw_avh-v1.0.zip",
#   bucket = bucket_name, 
#   region = region,
#   file = requested_path
# ) 
#save_object("s3://mybucket/input/test.zip", file = "/home/test.zip", bucket = "mybucket")


writeBin(s3_download$Body, con = file_name2)
count_data = as.data.frame(read.table(file_name2, sep = ",", header = T)) %>% head



#### until here

file_name1 = file.path(username, run_id, count_file)

s3$s3_list_objects(b)


s3_download <- s3$get_object(
  Bucket = bucket_name,
  Key = file_name1,

)

#if i want to write in the temdir folder
#file_name2 <- file.path(tempdir(), paste(username, run_id, count_file, sep = "_"))
#writeBin(s3_download$Body, con = file_name2)
#count_data = as.data.frame(read.table(file_name2, sep = ",", header = T)) %>% head

#or direclty load:
count_data <- s3_download$Body %>% rawToChar %>% read.csv(text = .)



###


username = "alex"
s3$list_objects(Bucket =  file.path("forjazul.galaxy", username , "my_history"))

s3$list_objects(Bucket = "forjazul.galaxy")

s3$get_bucket_location(
  Bucket = "forjazul.galaxy"
)

s3$get_bucket(
  Bucket = "forjazul.galaxy"
)


s3_other <- paws::s3(
  config = list(
    credentials = list(
      creds = list(
        access_key_id = my_access_key,
        secret_access_key = my_secret_key,
        session_token = my_token
      )
    ),
    region = my_region
  )
)


s3$list_objects(Bucket = "forjazul.galaxy", regio )

s3://forjazul.galaxy/<username>/<historyname>/

s3://forjazul.galaxy/alex/my_history/tcpm.csv

Bucket name: forjazul.galaxy


