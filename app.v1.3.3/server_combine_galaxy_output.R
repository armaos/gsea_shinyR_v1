files = c()


files <- rv$galaxy_outputs$path

sampleName<-gsub("^.*kallisto\\/","",files)
sampleName<-gsub("\\/abundance.*$","",sampleName)


txi.kallisto.tsv<- tximport(files, type = "kallisto", ignoreAfterBar = TRUE, txOut=T, )     #tx2gene = tx2gene,
txi<-txi.kallisto.tsv
#head(txi)
#tx2gene = tx2gene we can ignore

cts <- txi$counts
colnames(cts)<-sampleName
normMat <- txi$length
colnames(normMat)<-sampleName
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

