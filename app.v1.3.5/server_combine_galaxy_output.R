files = c()
# rv$galaxy_outputs = data.frame(file = c("B1L_S8_001.fastq.tabular" , "G2_S5_001.fastq.tabular"),
#                                path = c("B1L_S8_001.fastq.tabular" , "G2_S5_001.fastq.tabular"))
                               
files <- data.frame(file =  rv$aws_files) %>%
  left_join(rv$galaxy_outputs)

files <- files$path
#files <- rv$galaxy_outputs$path
#files <- rv$aws_files
#files  <- list.files("~/pipeline_output/")
#files <- list.files("~/OneDrive - Fondazione Istituto Italiano Tecnologia/B_CRO_projects/re_algae_c026/Meiks_analysis/", ignore.case = T, pattern = ".tabular")
sampleNames = data.frame(names = basename(files)) %>% 
  left_join(rv$custom_conditions_names)

# print("files")
# print(files)
#print(sampleNames)

#sampleName<-gsub("^.*kallisto\\/","",files)
#sampleName<-gsub("\\/abundance.*$","",sampleName)

#sampleName <- rv$custom_conditions_names$values


txi <- tximport(files, type = "kallisto", ignoreAfterBar = TRUE, txOut=T, )     #tx2gene = tx2gene,
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


# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
#library(edgeR)
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)
#head(y)
# filtering
keep <- filterByExpr(y)
## Warning in filterByExpr.DGEList(y): All samples appear to belong to the same
## group.
y <- y[keep, ]
data = as.data.frame(cpm(y))

data$Geneid = row.names(data)

row.names(data) = NULL


data = data %>% select(Geneid, everything())




source(file = file.path(version, "server_data_preprocessing.R"),
       local = TRUE,
       encoding = "UTF-8")

#write.table(cpm(y),"cpm_values_redalgar.csv",sep="\t",col.names=NA)
# y is now ready for estimate dispersion functions (edgeR)

