# shinyjs::show("overview_text")
filtered = rv$filtered

rv$condition  <- input$choose_conditions
shinyjs::html("pageHeader", rv$condition)


print("start DE")
donor <- factor(gsub("\\.\\w+$","",colnames(filtered),perl=T))
group <- factor(gsub("^\\w+\\.","",colnames(filtered),perl=T))
x<-group
#input
design <- model.matrix(~0+group)             #when paired analysis required, this need to be ~0+group+donor
rownames(design)<-colnames(filtered)
colnames(design)<-levels(factor(group))
print("design OK")
########################################################################
html("overview_text", "Processing normalization...")
print("Proc Norm")
y <- DGEList(counts=filtered, group=group)         	# 
keep <- rowSums(cpm(y) > 1 ) >= 2
y <- y[keep,,keep.lib.sizes=FALSE]    

y <- calcNormFactors(y,method="none") #tximport norm before

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
                       Group=x, 
                       X=pca$x[,1],
                       Y=pca$x[,2])
print("Proc PCA plot")
html("overview_text", "Processing PCA plot...")
rv$pca_plot = ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, col=Group)) +
  geom_point(size=4, alpha=0.6 ) + geom_label_repel(size = 4, segment.color='lightgrey')+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+ 
  stat_ellipse(data=pca.data, aes(fill=Group), type = "t", level=0.95, geom="polygon",alpha=0.2,color=NA)+guides(fill=F) +
  theme_bw()+theme(legend.position = "right")+ggtitle("Samples")

########################################################################
html("overview_text", "Processing Constrasts...")

print(paste("Proc Constrast", input$choose_conditions))
cmd <- paste("my.contrasts <- makeContrasts(", input$choose_conditions, ", levels = design)", sep ='"')
#cmd <- paste("my.contrasts <- makeContrasts(", "G1-G2", ", levels = design)", sep ='"')
#my.contrasts = makeContrasts(G1_G2="G1-G2", levels=design)              # need to rewrite acccording to the task

eval(parse(text = cmd))


#my.contrasts = makeContrasts(cond, levels=design)              # need to rewrite acccording to the task
print("makeContrasts OK")


# count replicates per group
single_replicates = as.data.frame(design) %>% 
  gather(k, v) %>% 
  group_by(k) %>% 
  summarise(n = sum(v==1)) %>% 
  filter(n == 1) 


# if single_replicates >0 , this means that iunder the group we have only one replicate
if(dim(single_replicates)[1] > 1){
  y$common.dispersion <- 0.1
  fit <- glmFit(y, design, robust=TRUE)
  warning_message = "WARNING: to few replicates per treatment/subgroup."
  html("overview_text", warning_message)
  html("overview_text", paste("Processing Constrasts: Fitting glmFit...", "</br>","</br>", warning_message))
}else{
  y.Disp <- estimateDisp(y, design, robust = TRUE)
  warning_message = ""
  html("overview_text", "Processing Constrasts: Fitting glmFit...")
  fit <- glmFit(y.Disp, design, robust=TRUE)
}


lfcT=1.5
fdrT=0.01

html("overview_text", paste("Processing Constrasts: Fitting glmLRT...", "</br>","</br>", warning_message))
print("glmLRT")
lrt.D_O <- glmLRT(fit, contrast=my.contrasts)
# lrt.aDTP_WT <- glmQLFTest(fit, contrast=my.contrasts[,id])  # quasi-likelihood (QL)

# add FDR
deg.edger <- lrt.D_O$table[p.adjust(lrt.D_O$table$PValue, method = "BH") < fdrT, ]
#dim(deg.edger)

#################################################################
# final volcano plot           ##################################
print("Proc volcano Data")
html("overview_text", paste("Processing Constrasts: Processing volcano data...", "</br>","</br>", warning_message))

v0<-lrt.D_O$table

v0$FDR <-p.adjust(v0$PValue, method = "BH") 
deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1


deg.edger1 <- deg.edger0[order(abs(deg.edger0$logFC)*(-log10(deg.edger0$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger0)>150,150,nrow(deg.edger0))],]  #select top150
degnames1<-rownames(deg.edger1)
degnames1<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames1,invert=T,value=T)
#nrow(deg.edger1)
v0$genelabels<-NA
v0[rownames(v0) %in% degnames1,]$genelabels<-T  
options(ggrepel.max.overlaps = Inf)                            # show only DEGs
v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))
# volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
#   geom_point(alpha=0.6, size=1) + 
#   geom_label_repel(size = 3.5, segment.color='lightgrey', color='white',fontface='bold', force=0.4, aes(x =logFC, y =-log10(FDR), label = ifelse(genelabels == T, rownames(v0),""))) + 	
#   theme_bw(base_size=15) + theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
#   ggtitle(paste0("Volcanoplot of ",id)) + scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
#   scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
#   geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
# volcano0
print("Proc Volcano PLot")
html("overview_text", paste("Processing Constrasts: Processing volcano plot...", "</br>","</br>", warning_message)) 

rv$volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
  geom_point(alpha=0.6, size=1) + 
  theme_bw(base_size=15) + 
  theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  ggtitle(paste0("Volcanoplot of ",input$choose_conditions)) + 
  scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)


#################################################################

html("overview_text", paste("Processing Constrasts: Processing heatmap plot...", "</br>","</br>", warning_message)) 
print("Heatmap plot")
id1 <- log2(cpm(y)+1)
id1 <- y
group_1 = str_split(input$choose_conditions, "-")[[1]][1]
group_2 = str_split(input$choose_conditions, "-")[[1]][2]
print(group_1)
print(group_2)
selectedCols <- c(grep(paste0("\\.",group_1),colnames(y)), grep(paste0("\\.",group_2),colnames(y)))  #select some columns
print(selectedCols)
id2 <- id1[rownames(deg.edger0),selectedCols]
annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

#pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#pheatmap(id2,fontsize=3,angle_col=0,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_Heatmap.png"), annotation=annotation, width=6,height=12, show_rownames=F)

rv$pheatmap <- pheatmap(id2,fontsize=8, angle_col=0, cutree_rows=2, cutree_col=2, scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), annotation=annotation, width=6, height=12, show_rownames=F)

rv$msi_df <- v0 %>%
  mutate(Name = rownames(v0)) %>%
  select(Name, everything())

rownames(rv$msi_df) <- NULL

rv$update_annotations =  rv$update_annotations + 1
rv$update_msi_df =  rv$update_msi_df + 1
shinyjs::show("menu_norm_and_de_panel")
print("Done proc DE")
html("overview_text", paste("DE analysis: DONE...", "</br>","</br>", warning_message)) 

rv$show_norm_and_de_mout = T  

