############## PCA plot after normalization ############## 
html("overview_text", "Processing Data for PCA...")
print("Proc PCA")
data.matrix <- as.matrix(log2(cpm(y)+1))
data.matrix[is.na(data.matrix)]<-0                  # remove na
pca <- prcomp(t(data.matrix), scale.=F)
#summary(pca)                                        # check the result here
#plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
#barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
pca.data <- data.frame(Sample=rownames(pca$x), 
                       Group=group, 
                       X=pca$x[,1],
                       Y=pca$x[,2])
print("Proc PCA plot")
html("overview_text", "Processing PCA plot...")
pca_plot = ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, col=Group)) +
  geom_point(size=4, alpha=0.6 ) + geom_label_repel(size = 4, segment.color='lightgrey')+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+ 
  stat_ellipse(data=pca.data, aes(fill=Group), type = "t", level=0.95, geom="polygon",alpha=0.2,color=NA)+guides(fill=F) +
  theme_bw()+theme(legend.position = "right")+ggtitle("Samples")
html("overview_text", "Processing PCA plot...OK")