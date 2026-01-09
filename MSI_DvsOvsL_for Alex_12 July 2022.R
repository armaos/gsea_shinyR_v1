#project MSI, Kal_t, Kallisto quant, tximport norm before 
#D vs. O vs. L DEG, GSEA


###### not obligatory, only generate special colors
library(randomcoloR)
colors <- distinctColorPalette(20)
lighten <- function(color, factor = 0.4) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
colors <- c('purple','orange','darkred','blue','chartreuse','darkgoldenrod4','tan2','darkgreen','red3','darkmagenta','deeppink','violet','navy','red','dodgerblue','darkcyan')
colors<-sapply(colors,lighten)
###### not obligatory, only generate special colors

###### star edgeR ####################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library("edgeR")
getwd()
setwd("C:/Users/biokurs/Documents")
data0<-read.table(file = 'tcpm.csv', sep = ',', header = TRUE)
data0<-data0[,-ncol(data0)]
head(data0)
dim(data0)
colnames(data0)[colnames(data0)=="gene_name"]="GeneSymbol"
#################################
dim(data0)
data<-data0
install.packages("ggplot2")
library(ggplot2)
colnames(data)[colnames(data)=="Geneid"]="Geneid"
rowMeans <- apply(data,1, function(x) mean(as.numeric(x),na.rm=T))
c <- data[!duplicated(data$Geneid),]
rownames(c) <- c$Geneid
########################################
dim(c)
head(c)
countData.all<-as.matrix(c[2:ncol(c)])
rownames(countData.all)<-c$Geneid
dim(countData.all)
head(countData.all)

countData<-countData.all
dim(countData)
head(countData)
names(countData)
names(data0)
countData

#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RUVSeq")
#
library(RUVSeq)
plotRLE(countData,outline=FALSE,las=2)    # must appear bad due to "0"
########################################
# filter out "0"
filter<-apply(countData,1,function(x) length(x[x>10])>=ncol(countData)/5)   # 20967 gene features, not always 5, should be (ncol/4) or /5
filtered<-countData[filter, ]
dim(filtered) # 12491
head(filtered)
dev.off()
plotRLE(filtered,outline=FALSE,las=2)
plotPCA(filtered)
########################################
############## PCA plot before normalization ############## 

data.matrix <- as.matrix(log2(filtered+1))
data.matrix[is.na(data.matrix)]<-0                               # remove na
pca <- prcomp(t(data.matrix), scale=F)
summary(pca)                                                     # check the result here
plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
pca.data <- data.frame(Sample=rownames(pca$x), 
                       Donor=gsub('^.+\\.','',rownames(pca$x),perl=T), 
                       Group=gsub('\\..+$','',rownames(pca$x),perl=T),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data


Donor=gsub('^.+\\.','',rownames(pca$x),perl=T), 
Group=gsub('\\..+$','',rownames(pca$x),perl=T),
X=pca$x[,1],
Y=pca$x[,2])

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, col=Donor)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+ 
  stat_ellipse(data=pca.data, aes(fill=Group), type = "t", level=0.9, geom="polygon",alpha=0.2,color=NA)+guides(fill=F) +
  theme_bw()+theme(legend.position = "right")+ggtitle("Samples")

############################ 
ggsave("NEWPCA_processed_data.both.png",width=9,height=9,bg="white")

donor <- factor(gsub("\\.\\w+$","",colnames(filtered),perl=T))
group <- factor(gsub("^\\w+\\.","",colnames(filtered),perl=T))

x<-group
x

########################################################################
############# edgeR #############
library("edgeR")

#        Sample Donor Group        
#        B1.D     D    B1

########################################################################
design <- model.matrix(~0+group)             #when paired analysis required, this need to be ~0+group+donor
rownames(design)<-colnames(filtered)
colnames(design)<-levels(factor(group))
########################################################################

y <- DGEList(counts=filtered, group=group)         	# 
keep <- rowSums(cpm(y) > 1 ) >= 2
y <- y[keep,,keep.lib.sizes=FALSE]    	 

nrow(y)  				
nrow(data)

y <- calcNormFactors(y,method="none") #tximport norm before

plotRLE(log2(cpm(y)+1), outline=FALSE, col=colors[x], las=2)
plotPCA(log2(cpm(y)+1), col=colors[x]) 

############## PCA plot after normalization ############## 
data.matrix <- as.matrix(log2(cpm(y)+1))
data.matrix[is.na(data.matrix)]<-0                  # remove na
pca <- prcomp(t(data.matrix), scale.=F)
summary(pca)                                        # check the result here
plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
pca.data <- data.frame(Sample=rownames(pca$x), 
                       Group=x, 
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

library(ggplot2)    
library(ggrepel) 
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, col=Group)) +
  geom_point(size=4, alpha=0.6 ) + geom_label_repel(size = 4, segment.color='lightgrey')+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+ 
  stat_ellipse(data=pca.data, aes(fill=Group), type = "t", level=0.95, geom="polygon",alpha=0.2,color=NA)+guides(fill=F) +
  theme_bw()+theme(legend.position = "right")+ggtitle("Samples")
############## PCA plot after normalization ############## 

ggsave("fPCA_Group_processed_data.both.png",width=9,height=9, bg="white") 

########################################################################

my.contrasts = makeContrasts(D_O="D-O", L_O="L-O", D_L="D-L", levels=design)              # need to rewrite acccording to the task

########################################################################
#        Sample Donor Group        
#        B1.D     D    B1
#y <- filtered
head(y)
dim(y)
# check dispersion
install.packages("statmod")
library(statmod)
y.Disp <- estimateDisp(y, design, robust = TRUE)
plotBCV(y.Disp)
plotMeanVar(y.Disp, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE, ylim=c(1e+1,1e+12))

fit <- glmFit(y.Disp, design, robust=TRUE)
head(fit$coefficients)

# fit <- glmQLFit(y.Disp, design, robust=TRUE)       # quasi-likelihood (QL), more stringent, only with high DEGs
#################################################################
########################################


id="D_O"                        # comparison ID

lfcT=1.5
fdrT=0.01

lrt.D_O <- glmLRT(fit, contrast=my.contrasts[,id])
# lrt.aDTP_WT <- glmQLFTest(fit, contrast=my.contrasts[,id])  # quasi-likelihood (QL)

# add FDR
deg.edger <- lrt.D_O$table[p.adjust(lrt.D_O$table$PValue, method = "BH") < fdrT, ]
dim(deg.edger)
#################################################################


#################################################################
# final volcano plot           ##################################

library("pheatmap")
library("ggplot2")
library("ggrepel")
v0<-lrt.D_O$table

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1
nrow(deg.edger0)
# 2338
table(deg.edger0$logFC>0)
# 813 UP 1525 DOWN

deg.edger1 <- deg.edger0[order(abs(deg.edger0$logFC)*(-log10(deg.edger0$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger0)>150,150,nrow(deg.edger0))],]  #select top150
degnames1<-rownames(deg.edger1)
degnames1<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames1,invert=T,value=T)
nrow(deg.edger1)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames1,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0

volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0




#################################################################
ggsave(paste0(id,".lfc",lfcT,"fdr",fdrT,"_Volcanoplot.png"),width=11,height=11,bg="white")
#################################################################

id1 <- log2(cpm(y)+1)
id1 <- y
selectedCols <-c(grep("\\.O",colnames(y)),grep("\\.D",colnames(y)))  #select some columns
id2 <- id1[rownames(deg.edger0),selectedCols]
annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

pheatmap(id2,fontsize=3,angle_col=0,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_Heatmap.png"), annotation=annotation, width=6,height=12, show_rownames=F)

dev.off()
dev.new() 
########################################################################

write.table(cbind(deg.edger0,id2),paste0(id,"_Heatmap.csv"),sep="\t", col.names=NA)


# write all records
write.table(cbind(v0,id1[rownames(v0),selectedCols]),paste0(id,"_all_values.csv"),sep="\t", col.names=NA)

## check only meaning genenames

#################################################################
# final volcano plot           ##################################

library("pheatmap")
library("ggplot2")
library("ggrepel")
v0<-lrt.D_O$table
#deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5 & p.adjust(v0$PValue, method = "BH") < 0.1, ]     # FDR 0.1

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1


nrow(deg.edger0)
# 2338
table(deg.edger0$logFC>0)
# 813 UP 1525 DOWN

deg.edger.onlymeaningful <- deg.edger0[grep("^Gm\\d+|Rik$|^RP\\d+|^\\d_",rownames(deg.edger0),invert=T),]


deg.edger.onlymeaningful <- deg.edger.onlymeaningful[order(abs(deg.edger.onlymeaningful$logFC)*(-log10(deg.edger.onlymeaningful$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger.onlymeaningful)>150,150,nrow(deg.edger.onlymeaningful))],]  #select top150
degnames2<-rownames(deg.edger.onlymeaningful)
degnames2<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames2,invert=T,value=T)
#nrow(deg.edger2)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames2,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0
#################################################################
ggsave(paste0(id,".lfc",lfcT,"fdr",fdrT,"_onlymeaningful_Volcanoplot.png"),width=11,height=11,bg="white")
#################################################################
id1 <- y
#id1 <- log2(cpm(y)+1)
selectedCols <-c(grep("\\.O",colnames(y)),grep("\\.D",colnames(y)))  #select some columns
id2 <- id1[rownames(deg.edger.onlymeaningful),selectedCols]
annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_onlymeaningful_Heatmap.png"), annotation=annotation, width=6,height=9, show_rownames=F)

dev.off()
dev.new() 
########################################################################


#############################################################################
# fgsea  version 2020-11-03 #################################################
#############################################################################

library(fgsea)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)

l0<-lrt.D_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/MSI_pathways.gmt") 
#pathways3
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=50, maxSize=500) #5
dim(fgseaRes)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =2.45
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.00001 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.00001  & abs(fgseaResTidy0$NES) >=nesThres-0.01,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.00001 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEA.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAsig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEA.csv"),sep="\t",col.names=NA)

###################
###################

l0<-lrt.D_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/koala_KEGG_ko_pw.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.5
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAkoala.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAkoalasig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAkoala.csv"),sep="\t",col.names=NA)


###################
###################
###################

l0<-lrt.D_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/em_KEGG_Pathway_pw1.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.5
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAkegg_pw.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAkeggsig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAkegg.csv"),sep="\t",col.names=NA)


###################
###################
###################
dev.off()
dev.new()
l0<-lrt.D_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/em_KEGG_ko_pw.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.5
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAkegg_ko.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAkegg_kosig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAkegg_ko.csv"),sep="\t",col.names=NA)


###################
###################
###################

l0<-lrt.D_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/ips_TIGRFAM_pw.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.0
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=2)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAtigrfam.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAtigrfam_sig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAtigrfam.csv"),sep="\t",col.names=NA)


###################
###################
###################
dev.off()
dev.new()
l0<-lrt.D_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/metacyc.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=5, maxSize=500) #5
dim(fgseaRes)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.5
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(pval)<=0.05 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$pval<0.05  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$pval<0.05 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAmetacyc.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAmetacycsig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAmetacyc.csv"),sep="\t",col.names=NA)


###################
###################
###################
dev.off()
dev.new()
l0<-lrt.D_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/reactome.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=5, maxSize=500) #5
dim(fgseaRes)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =2.0
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAreactome.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAreactomesig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAreactome.csv"),sep="\t",col.names=NA)


id="L_O"                        # comparison ID

lfcT=1.5
fdrT=0.01

lrt.L_O<- glmLRT(fit, contrast=my.contrasts[,id])
# lrt.aDTP_WT <- glmQLFTest(fit, contrast=my.contrasts[,id])  # quasi-likelihood (QL)

# add FDR
deg.edger <- lrt.L_O$table[p.adjust(lrt.L_O$table$PValue, method = "BH") < fdrT, ]
dim(deg.edger)
#################################################################


#################################################################
# final volcano plot           ##################################

library("pheatmap")
library("ggplot2")
library("ggrepel")
v0<-lrt.L_O$table

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1
nrow(deg.edger0)
# 2024
table(deg.edger0$logFC>0)
# 1030 UP 994 DOWN

deg.edger1 <- deg.edger0[order(abs(deg.edger0$logFC)*(-log10(deg.edger0$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger0)>150,150,nrow(deg.edger0))],]  #select top150
degnames1<-rownames(deg.edger1)
degnames1<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames1,invert=T,value=T)
nrow(deg.edger1)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames1,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0
#################################################################
ggsave(paste0(id,".lfc",lfcT,"fdr",fdrT,"_Volcanoplot.png"),width=11,height=11,bg="white")
#################################################################

id1 <- log2(cpm(y)+1)
id1 <- y
selectedCols <-c(grep("\\.O",colnames(y)),grep("\\.L",colnames(y)))  #select some columns
id2 <- id1[rownames(deg.edger0),selectedCols]
annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

pheatmap(id2,fontsize=3,angle_col=0,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_Heatmap.png"), annotation=annotation, width=6,height=12, show_rownames=F)

dev.off()
dev.new() 
########################################################################

write.table(cbind(deg.edger0,id2),paste0(id,"_Heatmap.csv"),sep="\t", col.names=NA)


# write all records
write.table(cbind(v0,id1[rownames(v0),selectedCols]),paste0(id,"_all_values.csv"),sep="\t", col.names=NA)

## check only meaning genenames

#################################################################
# final volcano plot           ##################################

library("pheatmap")
library("ggplot2")
library("ggrepel")
v0<-lrt.L_O$table
#deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5 & p.adjust(v0$PValue, method = "BH") < 0.1, ]     # FDR 0.1

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1


nrow(deg.edger0)
# 2024
table(deg.edger0$logFC>0)
# 1030 UP 994 DOWN

deg.edger.onlymeaningful <- deg.edger0[grep("^Gm\\d+|Rik$|^RP\\d+|^\\d_",rownames(deg.edger0),invert=T),]


deg.edger.onlymeaningful <- deg.edger.onlymeaningful[order(abs(deg.edger.onlymeaningful$logFC)*(-log10(deg.edger.onlymeaningful$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger.onlymeaningful)>150,150,nrow(deg.edger.onlymeaningful))],]  #select top150
degnames2<-rownames(deg.edger.onlymeaningful)
degnames2<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames2,invert=T,value=T)
#nrow(deg.edger2)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames2,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0
#################################################################
ggsave(paste0(id,".lfc",lfcT,"fdr",fdrT,"_onlymeaningful_Volcanoplot.png"),width=11,height=11,bg="white")
#################################################################


#id1 <- log2(cpm(y)+1)
id1 <- y
selectedCols <-c(grep("\\.O",colnames(y)),grep("\\.L",colnames(y)))  #select some columns
id2 <- id1[rownames(deg.edger.onlymeaningful),selectedCols]
annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_onlymeaningful_Heatmap.png"), annotation=annotation, width=6,height=9, show_rownames=F)

dev.off()
dev.new() 
########################################################################



#############################################################################
# fgsea  version 2020-11-03 #################################################
#############################################################################

library(fgsea)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)

l0<-lrt.L_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/MSI_pathways.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=50, maxSize=500) #5
dim(fgseaRes)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =2.3
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.00001 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.00001  & abs(fgseaResTidy0$NES) >=nesThres-0.01,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.00001 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEA.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAsig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEA.csv"),sep="\t",col.names=NA)



###################
###################
###################

l0<-lrt.L_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/koala_KEGG_ko_pw.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.5
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAkoala.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAkoalasig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAkoala.csv"),sep="\t",col.names=NA)


###################
###################
###################

l0<-lrt.L_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/em_KEGG_Pathway_pw1.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =2.0
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAkegg_pw.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAkeggsig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAkegg.csv"),sep="\t",col.names=NA)


###################
###################
###################
dev.off()
dev.new()
l0<-lrt.L_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/em_KEGG_ko_pw.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.5
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAkegg_ko.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAkegg_kosig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAkegg_ko.csv"),sep="\t",col.names=NA)


###################
###################
###################

l0<-lrt.L_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/ips_TIGRFAM_pw.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.0
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=2)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAtigrfam.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAtigrfam_sig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAtigrfam.csv"),sep="\t",col.names=NA)


###################
###################
###################
dev.off()
dev.new()
l0<-lrt.L_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/metacyc.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=5, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.0
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(pval)<=0.05 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$pval<0.05  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$pval<0.05 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAmetacyc.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAmetacycsig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAmetacyc.csv"),sep="\t",col.names=NA)


###################
###################
###################
dev.off()
dev.new()
l0<-lrt.L_O$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/reactome.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=5, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.75
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.05 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.05  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.05 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=2)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAreactome.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAreactomesig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAreactome.csv"),sep="\t",col.names=NA)





id="D_L"                        # comparison ID

lfcT=1.5
fdrT=0.01

lrt.D_L<- glmLRT(fit, contrast=my.contrasts[,id])
# lrt.aDTP_WT <- glmQLFTest(fit, contrast=my.contrasts[,id])  # quasi-likelihood (QL)

# add FDR
deg.edger <- lrt.D_L$table[p.adjust(lrt.D_L$table$PValue, method = "BH") < fdrT, ]
dim(deg.edger)
# 3896
#################################################################


#################################################################
# final volcano plot           ##################################

library("pheatmap")
library("ggplot2")
library("ggrepel")
v0<-lrt.D_L$table
#deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5 & p.adjust(v0$PValue, method = "BH") < 0.1, ]     # FDR 0.1

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1
nrow(deg.edger0)
# 1163
table(deg.edger0$logFC>0)
# 205 UP 958 DOWN

deg.edger1 <- deg.edger0[order(abs(deg.edger0$logFC)*(-log10(deg.edger0$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger0)>150,150,nrow(deg.edger0))],]  #select top150
degnames1<-rownames(deg.edger1)
degnames1<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames1,invert=T,value=T)
nrow(deg.edger1)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames1,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0
#################################################################
ggsave(paste0(id,".lfc",lfcT,"fdr",fdrT,"_Volcanoplot.png"),width=11,height=11,bg="white")
#################################################################

id1 <- log2(cpm(y)+1)
id1 <- y
selectedCols <-c(grep("\\.L",colnames(y)),grep("\\.D",colnames(y)))  #select some columns
id2 <- id1[rownames(deg.edger0),selectedCols]
annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

#pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

pheatmap(id2,fontsize=5,angle_col=0,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_Heatmap.png"), annotation=annotation, width=6,height=12, show_rownames=F)

dev.off()
dev.new() 
########################################################################

write.table(cbind(deg.edger0,id2),paste0(id,"_Heatmap.csv"),sep="\t", col.names=NA)


# write all records
write.table(cbind(v0,id1[rownames(v0),selectedCols]),paste0(id,"_all_values.csv"),sep="\t", col.names=NA)





## check only meaning genenames

#################################################################
# final volcano plot           ##################################

library("pheatmap")
library("ggplot2")
library("ggrepel")
v0<-lrt.D_L$table
#deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5 & p.adjust(v0$PValue, method = "BH") < 0.1, ]     # FDR 0.1

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1


nrow(deg.edger0)
# 1163
table(deg.edger0$logFC>0)
# 205 UP 958 DOWN

deg.edger.onlymeaningful <- deg.edger0[grep("^Gm\\d+|Rik$|^RP\\d+|^\\d_",rownames(deg.edger0),invert=T),]


deg.edger.onlymeaningful <- deg.edger.onlymeaningful[order(abs(deg.edger.onlymeaningful$logFC)*(-log10(deg.edger.onlymeaningful$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger.onlymeaningful)>150,150,nrow(deg.edger.onlymeaningful))],]  #select top150
degnames2<-rownames(deg.edger.onlymeaningful)
degnames2<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames2,invert=T,value=T)
#nrow(deg.edger2)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames2,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0
#################################################################
ggsave(paste0(id,".lfc",lfcT,"fdr",fdrT,"_onlymeaningful_Volcanoplot.png"),width=11,height=11,bg="white")
#################################################################

#id1 <- log2(cpm(y)+1)
id1 <- y
selectedCols <-c(grep("\\.L",colnames(y)),grep("\\.D",colnames(y)))  #select some columns
id2 <- id1[rownames(deg.edger.onlymeaningful),selectedCols]
annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_onlymeaningful_Heatmap.png"), annotation=annotation, width=6,height=9, show_rownames=F)

dev.off()
dev.new() 
########################################################################
#################################################################

#############################################################################
# fgsea  version 2020-11-03 #################################################
#############################################################################

library(fgsea)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(patchwork)

l0<-lrt.D_L$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/MSI_pathways.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=50, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =2.44
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.00001 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.00001  & abs(fgseaResTidy0$NES) >=nesThres-0.01,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.00001 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEA.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAsig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEA.csv"),sep="\t",col.names=NA)

#############
#############
#############


l0<-lrt.D_L$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/koala_KEGG_ko_pw.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.5
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAkoala.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAkoalasig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAkoala.csv"),sep="\t",col.names=NA)


###################
###################
###################

l0<-lrt.D_L$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/em_KEGG_Pathway_pw1.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =2.0
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAkegg_pw.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAkeggsig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAkegg.csv"),sep="\t",col.names=NA)


###################
###################
###################
dev.off()
dev.new()
l0<-lrt.D_L$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/em_KEGG_ko_pw.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.5
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAkegg_ko.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAkegg_kosig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAkegg_ko.csv"),sep="\t",col.names=NA)


###################
###################
###################

l0<-lrt.D_L$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/ips_TIGRFAM_pw.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=10, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.0
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=2)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAtigrfam.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAtigrfam_sig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAtigrfam.csv"),sep="\t",col.names=NA)


###################
###################
###################
dev.off()
dev.new()
l0<-lrt.D_L$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/metacyc.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=5, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =1.0
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(pval)<=0.05 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$pval<0.05  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$pval<0.05 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=8)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAmetacyc.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAmetacycsig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAmetacyc.csv"),sep="\t",col.names=NA)


###################
###################
###################
dev.off()
dev.new()
l0<-lrt.D_L$table
l0$FDR <-p.adjust(l0$PValue, method = "BH")
l0$nameUp<-rownames(l0)   #toupper(rownames(l0))
l<-l0[!duplicated(l0$nameUp),]
#rownames(l)<-toupper(rownames(l))
pgenelist<- as.data.frame(l) %>% rownames_to_column(var="hgnc_symbol")
res2 <- pgenelist %>% dplyr::select(hgnc_symbol, logFC, FDR) %>%  na.omit() %>% distinct() %>%  group_by(hgnc_symbol) %>%  summarize(stat=mean(logFC))  #mean(-logFC*log10(FDR))
ranks <- deframe(res2)  # head(ranks, 20)
head(ranks)
pathways3 <- gmtPathways("~/reactome.gmt") 
names(pathways3)<-gsub('^test_','',names(pathways3),perl=T)
pathways.hallmark<-c(pathways3)
pathways.hallmark %>%  head() %>%  lapply(head)
set.seed(123)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, minSize=5, maxSize=500) #5
dim(fgseaRes)
#fgseaRes$pathway = stringr::str_replace(fgseaRes$pathway, "Immune_" , "")
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES)) 
fgseaResTidy
nesThres =2.0
selectedSig <-fgseaResTidy %>% dplyr::select(-leadingEdge, -ES,  leadingEdge) %>% filter(abs(padj)<=0.01 & abs(NES)>=nesThres) %>% arrange(desc(NES)) %>% print(n=Inf)
fgseaResTidy0<-fgseaResTidy
fgseaResTidy0<-fgseaResTidy0[fgseaResTidy0$padj<0.01  & abs(fgseaResTidy0$NES) >=nesThres-0.1,]
fgseaResTidy0$Significant<-'Not'
fgseaResTidy0$Significant[fgseaResTidy0$padj<0.01 & abs(fgseaResTidy0$NES) >=nesThres]<-'Sig.'
t<-theme(title=element_text(size=7, face='bold')) 
p0<-ggplot(fgseaResTidy0, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Significant)) + coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title=id)+t+  scale_fill_brewer(palette = "Paired") 
p0
pp<- plotEnrichment(pathways.hallmark[[selectedSig$pathway[1]]], ranks)+ labs(title=selectedSig$pathway[1], subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[1],]$NES)) +t
for(i in 2:length(selectedSig$pathway)) { pp<-pp+ (plotEnrichment(pathways.hallmark[[selectedSig$pathway[i]]], ranks)+ labs(title=selectedSig$pathway[i],subtitle=paste0("NES:",fgseaResTidy0[fgseaResTidy0$pathway==selectedSig$pathway[i],]$NES)) +t  ) }
pp<- pp+ plot_layout(ncol=2)
pp<-p0+ (pp) + plot_layout(design = "ABBBBB") 
pp

ggsave(paste0(id,"_Scott_GSEAreactome.png"),width=30,height=18,bg="white")
write.table(fgseaResTidy0[,-8],paste(id,"_Scott_GSEAreactomesig.csv"),sep="\t",col.names=NA)
write.table(fgseaResTidy[,-8],paste(id,"_Scott_GSEAreactome.csv"),sep="\t",col.names=NA)


###############
###############
###############
###############
###############
###############
###############

id="B_G"                        # comparison ID

lfcT=1.5
fdrT=0.01

lrt.B_G <- glmLRT(fit, contrast=my.contrasts[,id])
# lrt.aDTP_WT <- glmQLFTest(fit, contrast=my.contrasts[,id])  # quasi-likelihood (QL)

# add FDR
deg.edger <- lrt.B_G$table[p.adjust(lrt.B_G$table$PValue, method = "BH") < fdrT, ]
dim(deg.edger)
# 1059
#################################################################


#################################################################
# final volcano plot           ##################################

library("pheatmap")
library("ggplot2")
library("ggrepel")
v0<-lrt.B_G$table
#deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5 & p.adjust(v0$PValue, method = "BH") < 0.1, ]     # FDR 0.1

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1
nrow(deg.edger0)
# 636
table(deg.edger0$logFC>0)
# 485 UP 151 DOWN

deg.edger1 <- deg.edger0[order(abs(deg.edger0$logFC)*(-log10(deg.edger0$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger0)>150,150,nrow(deg.edger0))],]  #select top150
degnames1<-rownames(deg.edger1)
degnames1<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames1,invert=T,value=T)
nrow(deg.edger1)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames1,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0
#################################################################
ggsave(paste0(id,".lfc",lfcT,"fdr",fdrT,"_Volcanoplot.png"),width=11,height=11,bg="white")
#################################################################

id1 <- log2(cpm(y)+1)
#selectedCols <-c(grep("\\.B",colnames(y)),grep("\\.G",colnames(y)))  #select some columns
#id2 <- id1[rownames(deg.edger0),selectedCols]
id2 <- id1[rownames(deg.edger0),]
#annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
annotation <- data.frame(Group=gsub("\\d\\.L$","",colnames(id2)))


rownames(annotation) <- colnames(id2)              # check out the row names of annotation

#pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

pheatmap(id2,fontsize=3,angle_col=0,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_Heatmap.png"), annotation=annotation, width=6,height=12, show_rownames=F)

dev.off()
dev.new() 
########################################################################

write.table(cbind(deg.edger0,id2),paste0(id,"_Heatmap.csv"),sep="\t", col.names=NA)
# write all records
write.table(cbind(v0,id1[rownames(v0),]),paste0(id,"_all_values.csv"),sep="\t", col.names=NA)



#################################################################
##check only meaning genenames###################################
library("pheatmap")
library("ggplot2")
library("ggrepel")
v0<-lrt.B_G$table
#deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5 & p.adjust(v0$PValue, method = "BH") < 0.1, ]     # FDR 0.1

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1


nrow(deg.edger0)
# 636
table(deg.edger0$logFC>0)
# 485 UP 151 DOWN

deg.edger.onlymeaningful <- deg.edger0[grep("^Gm\\d+|Rik$|^RP\\d+|^\\d_",rownames(deg.edger0),invert=T),]


deg.edger.onlymeaningful <- deg.edger.onlymeaningful[order(abs(deg.edger.onlymeaningful$logFC)*(-log10(deg.edger.onlymeaningful$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger.onlymeaningful)>150,150,nrow(deg.edger.onlymeaningful))],]  #select top150
degnames2<-rownames(deg.edger.onlymeaningful)
degnames2<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames2,invert=T,value=T)
#nrow(deg.edger2)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames2,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0
#################################################################
ggsave(paste0(id,".lfc",lfcT,"fdr",fdrT,"_onlymeaningful_Volcanoplot.png"),width=11,height=11,bg="white")
#################################################################

id1 <- log2(cpm(y)+1)
#selectedCols <-c(grep("\\.B",colnames(y)),grep("\\.G",colnames(y)),grep("\\.L",colnames(y)))  #select some columns
#id2 <- id1[rownames(deg.edger.onlymeaningful),selectedCols]
id2 <- id1[rownames(deg.edger.onlymeaningful),]
#annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
annotation <- data.frame(Group=gsub("\\d\\.L$","",colnames(id2)))
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

#pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_onlymeaningful_Heatmap.png"), annotation=annotation, width=6,height=9, show_rownames=T)


dev.off()
dev.new()

########################################################################
########################################################################
########################################################################

donor <- factor(gsub("\\.\\w+$","",colnames(filtered),perl=T))
group <- factor(gsub("^\\w+\\.","",colnames(filtered),perl=T))

x<-group


########################################################################
############# edgeR #############
library("edgeR")

#        Sample Donor Group        
#        B1.D     D    B1

########################################################################
design <- model.matrix(~0+group)             #when paired analysis required, this need to be ~0+group+donor
rownames(design)<-colnames(filtered)
colnames(design)<-levels(factor(group))
########################################################################

y <- DGEList(counts=filtered, group=group)         	# 
keep <- rowSums(cpm(y) > 1 ) >= 2
y <- y[keep,,keep.lib.sizes=FALSE]    	 

nrow(y)  						#    12143


y <- calcNormFactors(y,method="TMM")

plotRLE(log2(cpm(y)+1), outline=FALSE, col=colors[x], las=2)
plotPCA(log2(cpm(y)+1), col=colors[x]) 

############## PCA plot after normalization ############## 
data.matrix <- as.matrix(log2(cpm(y)+1))
data.matrix[is.na(data.matrix)]<-0                  # remove na
pca <- prcomp(t(data.matrix), scale.=F)
summary(pca)                                        # check the result here
plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
pca.data <- data.frame(Sample=rownames(pca$x), 
                       Group=x, 
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

library(ggplot2)    
library(ggrepel) 
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, col=Group)) +
  geom_point(size=4, alpha=0.6 ) + geom_label_repel(size = 4, segment.color='lightgrey')+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+ 
  stat_ellipse(data=pca.data, aes(fill=Group), type = "t", level=0.95, geom="polygon",alpha=0.2,color=NA)+guides(fill=F) +
  theme_bw()+theme(legend.position = "right")+ggtitle("Samples")
############## PCA plot after normalization ############## 

ggsave("PCA_Group_processed_data.both.png",width=9,height=9, bg="white") 

#loading_scores <- pca$rotation[,1]
#gene_scores <- abs(loading_scores)
#gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
#top_80_genes <- names(gene_score_ranked[1:80])                  
#top_80_genesymbol<-rownames(data.matrix[top_80_genes,])         

########################################################################

my.contrasts = makeContrasts(D_O="D-O", L_O="L-O", D_L="D-L", A_O="A-O", levels=design)              # need to rewrite acccording to the task

########################################################################
#        Sample Donor Group        
#        B1.D     D    B1

# check dispersion
install.packages("statmod")
library(statmod)
y.Disp <- estimateDisp(y, design, robust = TRUE)
plotBCV(y.Disp)
plotMeanVar(y.Disp, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE, ylim=c(1e+1,1e+12))


fit <- glmFit(y.Disp, design, robust=TRUE)
head(fit$coefficients)

# fit <- glmQLFit(y.Disp, design, robust=TRUE)       # quasi-likelihood (QL), more stringent, only with high DEGs
#################################################################

id="A_O"                        # comparison ID

lfcT=1.5
fdrT=0.01

lrt.A_O <- glmLRT(fit, contrast=my.contrasts[,id])
# lrt.aDTP_WT <- glmQLFTest(fit, contrast=my.contrasts[,id])  # quasi-likelihood (QL)

# add FDR
deg.edger <- lrt.A_O$table[p.adjust(lrt.A_O$table$PValue, method = "BH") < fdrT, ]
dim(deg.edger)
# 7898
#################################################################


#################################################################
# final volcano plot           ##################################

library("pheatmap")
library("ggplot2")
library("ggrepel")
v0<-lrt.A_O$table
#deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5 & p.adjust(v0$PValue, method = "BH") < 0.1, ]     # FDR 0.1

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1
nrow(deg.edger0)
# 7818
table(deg.edger0$logFC>0)
# 238 UP 7580 DOWN

deg.edger1 <- deg.edger0[order(abs(deg.edger0$logFC)*(-log10(deg.edger0$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger0)>150,150,nrow(deg.edger0))],]  #select top150
degnames1<-rownames(deg.edger1)
degnames1<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames1,invert=T,value=T)
nrow(deg.edger1)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames1,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0
#################################################################
ggsave(paste0(id,".lfc",lfcT,"fdr",fdrT,"_Volcanoplot.png"),width=11,height=11,bg="white")
#################################################################

id1 <- log2(cpm(y)+1)
selectedCols <-c(grep("\\.O",colnames(y)),grep("\\.A",colnames(y)))  #select some columns
id2 <- id1[rownames(deg.edger0),selectedCols]
annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

#pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

pheatmap(id2,fontsize=3,angle_col=0,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_Heatmap.png"), annotation=annotation, width=6,height=12, show_rownames=F)

dev.off()
dev.new() 
########################################################################

write.table(cbind(deg.edger0,id2),paste0(id,"_Heatmap.csv"),sep="\t", col.names=NA)
# write all records
write.table(cbind(v0,id1[rownames(v0),selectedCols]),paste0(id,"_all_values.csv"),sep="\t", col.names=NA)




MEIK_NEW###################### MEIK_NEW ################################
##############################################################
##############################################################

######################## control vs. treatment average ##################
my.contrasts = makeContrasts(DL_O="(D+L)/2-O", levels=design)              # need to rewrite acccording to the task
#makeContrasts((treatmentI+treatmentII+treatmentIII)/3-treatmentCTL,
#              levels=colnames(design))
#https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html

install.packages("statmod")
library(statmod)
y.Disp <- estimateDisp(y, design, robust = TRUE)
plotBCV(y.Disp)
plotMeanVar(y.Disp, show.raw=TRUE, show.tagwise=TRUE, show.binned=TRUE, ylim=c(1e+1,1e+12))


fit <- glmFit(y.Disp, design, robust=TRUE)
head(fit$coefficients)

# fit <- glmQLFit(y.Disp, design, robust=TRUE)       # quasi-likelihood (QL)
#################################################################

id="DL_O"                        # comparison ID

lfcT=1.5
fdrT=0.01

lrt.DL_O <- glmLRT(fit, contrast=my.contrasts[,id])
# lrt.aDTP_WT <- glmQLFTest(fit, contrast=my.contrasts[,id])  # quasi-likelihood (QL)

# add FDR
deg.edger <- lrt.DL_O$table[p.adjust(lrt.DL_O$table$PValue, method = "BH") < fdrT, ]
dim(deg.edger)
# 6211
#################################################################
# final volcano plot           ##################################

library("pheatmap")
library("ggplot2")
library("ggrepel")
v0<-lrt.DL_O$table
#deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5 & p.adjust(v0$PValue, method = "BH") < 0.1, ]     # FDR 0.1

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.01 and lfc 1.5
nrow(deg.edger0)
# 1965
table(deg.edger0$logFC>0)
# 775 UP 1190 DOWN

deg.edger1 <- deg.edger0[order(abs(deg.edger0$logFC)*(-log10(deg.edger0$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger0)>150,150,nrow(deg.edger0))],]  #select top150
degnames1<-rownames(deg.edger1)
degnames1<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames1,invert=T,value=T)
nrow(deg.edger1)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames1,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0
#################################################################
ggsave(paste0(id,".lfc",lfcT,"fdr",fdrT,"_Volcanoplot.png"),width=11,height=11,bg="white")
#################################################################

id1 <- log2(cpm(y)+1)
selectedCols <-c(grep("\\.O",colnames(y)),  grep("\\.D",colnames(y)),grep("\\.L",colnames(y)))  #select some columns
id2 <- id1[rownames(deg.edger0),selectedCols]
annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

pheatmap(id2,fontsize=3,angle_col=0,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_Heatmap.png"), annotation=annotation, width=6,height=12, show_rownames=F)

dev.off()
dev.new() 
########################################################################

write.table(cbind(deg.edger0,id2),paste0(id,"_Heatmap.csv"),sep="\t", col.names=NA)
# write all records
write.table(cbind(v0,id1[rownames(v0),selectedCols]),paste0(id,"_all_values.csv"),sep="\t", col.names=NA)


#################################################################
##check only meaning genenames###################################
library("pheatmap")
library("ggplot2")
library("ggrepel")
v0<-lrt.DL_O$table
#deg.edger0 <- v0[v0$PValue<=0.01 & abs(v0$logFC)>=1.5 & p.adjust(v0$PValue, method = "BH") < 0.1, ]     # FDR 0.1

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1


nrow(deg.edger0)
# 1965
table(deg.edger0$logFC>0)
# 775 UP 1190 DOWN

deg.edger.onlymeaningful <- deg.edger0[grep("^Gm\\d+|Rik$|^RP\\d+|^\\d_",rownames(deg.edger0),invert=T),]


deg.edger.onlymeaningful <- deg.edger.onlymeaningful[order(abs(deg.edger.onlymeaningful$logFC)*(-log10(deg.edger.onlymeaningful$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger.onlymeaningful)>150,150,nrow(deg.edger.onlymeaningful))],]  #select top150
degnames2<-rownames(deg.edger.onlymeaningful)
degnames2<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames2,invert=T,value=T)
#nrow(deg.edger2)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames2,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
  theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
volcano0
#################################################################
ggsave(paste0(id,".lfc",lfcT,"fdr",fdrT,"_onlymeaningful_Volcanoplot.png"),width=11,height=11,bg="white")
#################################################################


id1 <- log2(cpm(y)+1)
selectedCols <-c(grep("\\.O",colnames(y)),grep("\\.D",colnames(y)),grep("\\.L",colnames(y)))  #select some columns
id2 <- id1[rownames(deg.edger.onlymeaningful),selectedCols]
annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

#pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_onlymeaningful_Heatmap.png"), annotation=annotation, width=6,height=9, show_rownames=T)

dev.off()
dev.new() 

######################## END D vs. O vs. L ##################
#############################################################################