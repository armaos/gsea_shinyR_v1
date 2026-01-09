withProgress(message = 'Data preprocessing', value = 0.1, {
  print("combine galaxy output Start")
  
  incProgress(0.1, detail = paste("Retrieving Count Data..."))
  
  files <- c()
  txi <- c()
  cts <- c()
  normMat <- c()
  eff.lib <- c()
  y <- c()
  data <- data.frame()
  
  print("rv$aws_files")
  print(rv$aws_files)
  print("input$galaxy_outputs")
  print(input$galaxy_outputs)
  print("rv$galaxy_outputs")
  print(rv$galaxy_outputs)
  files <- data.frame(file2 =  rv$aws_files) %>%
    left_join(rv$galaxy_outputs)
  
  #files <- files$path
  print("files")
  print(files)
  #files <- rv$galaxy_outputs$path
  #files <- rv$aws_files
  #files  <- list.files("~/pipeline_output/")
  #files <- list.files("~/OneDrive - Fondazione Istituto Italiano Tecnologia/B_CRO_projects/re_algae_c026/Meiks_analysis/", ignore.case = T, pattern = ".tabular")
  sampleNames = data.frame(names = files$file2) %>% 
    left_join(rv$metadata_df)
  
  
  #sampleName<-gsub("^.*kallisto\\/","",files)
  #sampleName<-gsub("\\/abundance.*$","",sampleName)
  
  #sampleName <- rv$metadata_df$values
  
  incProgress(0.1, detail = paste("Reading Count Data..."))
  txi <- tximport(files$path, type = "kallisto", ignoreAfterBar = TRUE, txOut=T )     #tx2gene = tx2gene,
  #head(txi)
  #tx2gene = tx2gene we can ignore
  
  cts <- txi$counts
  colnames(cts)<-sampleNames$values
  normMat <- txi$length
  colnames(normMat)<-sampleNames$values
  # Obtaining per-observation scaling factors for length, adjusted to avoid
  # changing the magnitude of the counts.
  normMat <- normMat/exp(rowMeans(log(normMat)))
  normCts <- cts/normMat
  
  incProgress(0.1, detail = paste("Estimating DE..."))
  
  # Computing effective library sizes from scaled counts, to account for
  # composition biases between samples.
  #library(edgeR)
  eff.lib <- calcNormFactors(normCts) * colSums(normCts)
  # Combining effective library sizes with the length factors, and calculating
  # offsets for a log-link GLM.
  normMat <- sweep(normMat, 2, eff.lib, "*")
  incProgress(0.1, detail = paste("Estimating DE..."))
  normMat <- log(normMat)
  # Creating a DGEList object for use in edgeR.
  y <- DGEList(cts)
  incProgress(0.1, detail = paste("Estimating DE..."))
  y <- scaleOffset(y, normMat)
  #head(y)
  # filtering
  keep <- filterByExpr(y)
  ## Warning in filterByExpr.DGEList(y): All samples appear to belong to the same
  ## group.
  y <- y[keep, ]
  data = as.data.frame(cpm(y))
  
  data$Geneid = row.names(data)
  
  incProgress(0.1, detail = paste("Estimating DE..."))
  row.names(data) = NULL
  
  
  data = data %>% select(Geneid, everything())
  
})

print("combine galaxy output END. Heading to server_data_preprocessing.R")
source(file = file.path(version, "server_data_preprocessing.R"),
       local = TRUE,
       encoding = "UTF-8")

#write.table(cpm(y),"cpm_values_redalgar.csv",sep="\t",col.names=NA)
# y is now ready for estimate dispersion functions (edgeR)

