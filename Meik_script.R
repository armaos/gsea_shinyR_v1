###   
###
# R pipeline for kallisto normalization and DEG obj
###   
###   
library(dplyr)
library(tidyr)


path = "~/OneDrive - Fondazione Istituto Italiano Tecnologia/B_CRO_projects/re_algae_c026/Meiks_analysis/"
setwd(path)
files <- list.files("~/OneDrive - Fondazione Istituto Italiano Tecnologia/B_CRO_projects/re_algae_c026/Meiks_analysis/", ignore.case = T, pattern = ".tabular")


path = "~/OneDrive - Fondazione Istituto Italiano Tecnologia/B_CRO_projects/re_algae_c026/gsea_shinyR/tmp/tmp1/"
setwd(path)
files <- list.files("~/OneDrive - Fondazione Istituto Italiano Tecnologia/B_CRO_projects/re_algae_c026/gsea_shinyR/tmp/tmp1", ignore.case = T, pattern = ".tabular")
#files <- files[1:3]
#files <- c(files[1], files[1])
files
sampleName<-gsub("^.*kallisto\\/","",files)
sampleName<-gsub("\\/abundance.*$","",sampleName)
sampleName
# file.path(dir, "kallisto", "abundance.tsv.gz")
# file.path(dir, "kallisto", samples$run, "abundance.tsv.gz")
#names(files) <- paste0("sample", 1:6)
library("tximport")
txi.kallisto.tsv<- tximport(files, type = "kallisto", ignoreAfterBar = TRUE, txOut=T, )     #tx2gene = tx2gene,
txi<-txi.kallisto.tsv
#head(txi)

#tx2gene = tx2gene we can ignore
library(edgeR)
cts <- txi$counts
as.data.frame( txi$counts) %>% gather(k,v) %>% filter(v!=0 ) %>% arrange(desc(v)) %>% head
colnames(cts)<-sampleName
normMat <- txi$length
colnames(normMat)<-sampleName
# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat


# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
library(edgeR)
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)
head(y,5)
head(y$counts[y$counts !=0], 4)
as.data.frame( y$counts) %>% gather(k,v) %>% filter(v!=0 ) %>% arrange(desc(v)) %>% head

# filtering
keep <- filterByExpr(y)
## Warning in filterByExpr.DGEList(y): All samples appear to belong to the same
## group.
y <- y[keep, ]
data = as.data.frame(cpm(y))
dim(data)
data$Geneid = row.names(data)
row.names(data) = NULL
head(data)

write.table(cpm(y),"cpm_values_redalgar.csv",sep="\t",col.names=NA)
# y is now ready for estimate dispersion functions (edgeR)



head(y)

