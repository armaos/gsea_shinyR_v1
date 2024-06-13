cond_clean = ifelse(grepl("_in_", cond ), str_split(cond, "_in_", simplify = T)[,1] ,cond ) 
conditional_treatment = ifelse(grepl("_in_", cond ), str_split(cond, "_in_", simplify = T)[,2] , "" ) 
type_constrast = type_constrast
lfcT=1.5
fdrT=0.05

group_1 = str_split(cond_clean, "-")[[1]][1]
group_2 = str_split(cond_clean, "-")[[1]][2]
print(group_1)
print(group_2)


print(paste("Runnign condition contrast", cond_clean))

if(conditional_treatment == ""){
  filtered_to_contrast = filtered[, c(grep(paste0("\\.",group_1),colnames(filtered)), grep(paste0("\\.",group_2),colnames(filtered)))]   
  
  group_by = factor(gsub("^\\w+\\.","",colnames(filtered_to_contrast),perl=T))  # this is the treatment
}else{
  filtered_to_contrast = filtered[,grep(paste0("[.]", conditional_treatment ,"$"),colnames(filtered),perl=T)]   
  filtered_to_contrast = filtered_to_contrast[, c(grep(paste0("^",group_1,"\\d*[.]"), colnames(filtered_to_contrast), perl = T), 
                                                  grep(paste0("^",group_2,"\\d*[.]"), colnames(filtered_to_contrast), perl = T))]
                                                  
  sst = get_group_subgroup_treatment(colnames(filtered_to_contrast))
  
  subgroup = sst[[1]]
  subgroup_class = sst[[2]]
  treatment = sst[[3]]
  
  group_by = subgroup_class
}


## make design of the contrast
design <- model.matrix(~0+group_by)             #when paired analysis required, this need to be ~0+group+donor
rownames(design)<-colnames(filtered_to_contrast)
colnames(design)<-levels(factor(group_by))

# count replicates per group
single_replicates = as.data.frame(design) %>% 
  gather(k, v) %>% 
  group_by(k) %>% 
  summarise(n = sum(v==1)) %>% 
  filter(n == 1) 

print(colnames(filtered_to_contrast))

y <- calculate_y(filtered_to_contrast, group_by)
print(paste("Done y"))  
incProgress(0.1, detail = paste("Fitting glmFit..."))



y_fit_list  = make_fit_y(single_replicates, y, design)
y = y_fit_list[[1]]
fit = y_fit_list[[2]]

group = group_by
source(file = file.path(version, "server_PCA_plot.R"),
       local = TRUE,
       encoding = "UTF-8")
rv[[paste0("pca_", cond)]] = pca_plot
print("PCA processed")


cmd <- paste("my.contrasts <- makeContrasts(", cond_clean, ", levels = design)", sep ='"')
#cmd <- paste("my.contrasts <- makeContrasts(", "G1-G2", ", levels = design)", sep ='"')
#my.contrasts = makeContrasts(G1_G2="G1-G2", levels=design)              # need to rewrite acccording to the task

eval(parse(text = cmd))


#my.contrasts = makeContrasts(cond, levels=design)              # need to rewrite acccording to the task
print("makeContrasts OK")

lrt.D_O <- glmLRT(fit, contrast=my.contrasts)
# lrt.aDTP_WT <- glmQLFTest(fit, contrast=my.contrasts[,id])  # quasi-likelihood (QL)

# add FDR
deg.edger <- lrt.D_O$table[p.adjust(lrt.D_O$table$PValue, method = "BH") < fdrT, ]
#dim(deg.edger)

#################################################################
# final volcano plot           ##################################
print("Proc volcano Data")
html("overview_text", paste("Processing Constrasts: Processing volcano data...", "</br>","</br>", warning_message))
setProgress(detail = paste("Processing Constrast:", cond, ":`volcano plot"))

v0<-lrt.D_O$table


v0$FDR <-p.adjust(v0$PValue, method = "BH") 

deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1
if(dim(deg.edger0)[1] == 0 ){
  lfcT = 0
  fdrT = 1
  deg.edger0 <- v0[abs(v0$logFC)>=lfcT & v0$FDR <= fdrT, ]  
  showNotification(paste("Not enough DEG.", cond), duration = 5, type = "error")
  showNotification(paste("LFC threshold set to 0"), duration = 5, type = "error")
  showNotification(paste("FDR threshold set to 1"), duration = 5, type = "error")
}
print("dim(deg.edger0)[1]")
print(dim(deg.edger0)[1])
print(paste(lfcT, fdrT))
#deg.edger0 <- v0[abs(v0$logFC)>=500 & v0$FDR <= fdrT, ]     # FDR 0.05 and pV 1

deg.edger1 <- deg.edger0[order(abs(deg.edger0$logFC)*(-log10(deg.edger0$FDR)),decreasing=T)[1:ifelse(nrow(deg.edger0)>150,150,nrow(deg.edger0))],]  #select top150
degnames1<-rownames(deg.edger1)

degnames1<- grep("^Gm\\d+|Rik$|^RP\\d+",degnames1,invert=T,value=T)
#nrow(deg.edger1)

v0$genelabels<-NA

if(dim(deg.edger1 %>% filter(!is.na(PValue)))[1] > 0){
  v0[rownames(v0) %in% degnames1,]$genelabels <- T  
}


options(ggrepel.max.overlaps = Inf)                            # show only DEGs

v0$change<-as.factor(ifelse(v0$FDR<=fdrT & abs(v0$logFC)>=lfcT, ifelse(v0$logFC>=lfcT,"Up","Down"),"NotSig"))

print("Proc Volcano PLot")
html("overview_text", paste("Processing Constrasts: Processing volcano plot...", "</br>","</br>", warning_message)) 
volcano_plot_data_list[[cond]] <<- v0

# output[[paste0("volcano0_", cond)]] <- renderPlot({
#   ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
#   geom_point(alpha=0.6, size=1) + 
#   theme_bw(base_size=15) + 
#   theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
#   ggtitle(paste0("Volcanoplot of ", cond )) + 
#   scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
#   scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
#   geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
# })
print(paste0("Volcanoplot of ", cond ))
#output[[paste0("volcano0_", cond)]] <- renderPlot({rv[[paste0("volcano0_", cond)]] })
#################################################################

html("overview_text", paste("Processing Constrasts: Processing heatmap plot...", "</br>","</br>", warning_message)) 
setProgress(detail = paste("Processing Constrast:", cond, ":`heatmap plot"))
print("Heatmap plot")
id1 <- log2(cpm(y)+1)
id1 <- y


if(type_constrast == "inter"){
  selectedCols <- c(grep(paste0("\\.",group_1),colnames(y)), grep(paste0("\\.",group_2),colnames(y)))  #select some columns 
}else{
  selectedCols <- c(grep(paste0(group_1,"[1-9,a-z,A-Z]*\\."),colnames(y)), grep(paste0(group_2,"[1-9,a-z,A-Z]*\\."),colnames(y)))  #select some columns 
}
print(selectedCols)

id2 <- id1[rownames(deg.edger0),selectedCols]
annotation <- data.frame(Group=gsub("^.*\\.","",colnames(id2)))
print(annotation)
rownames(annotation) <- colnames(id2)              # check out the row names of annotation

#pheatmap(id2,fontsize=8,angle_col=45,cutree_rows=2,cutree_col=1,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, annotation=annotation, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#pheatmap(id2,fontsize=3,angle_col=0,cutree_rows=2,cutree_col=2,scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), file=paste0(id,"_Heatmap.png"), annotation=annotation, width=6,height=12, show_rownames=F)


print("doing heatmap plot")
print(dim(id2))
heatmap_plot_data_list[[cond]] <- id2

rv[[paste0("pheatmap_", cond)]] <- tryCatch(
  {
    pheatmap(
              id2,
              angle_col=90,
              fontsize=8,

              cutree_rows=2,
              cutree_col=2,
              scale="row",
              cluster_cols=F,
              treeheight_col=5,
              treeheight_row=15,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
              annotation=annotation,
              width=6,
              height=12,
              show_rownames=F,
              main = paste0("Expression of all DEGs ( with FDR < " ,fdrT ," and |logFC| > ", lfcT ," ) for each sample"))
  },
  error=function(cond) {
    empty_image()
  }
)


#output[[paste0("pheatmap_", cond)]] <- renderPlot({rv[[paste0("pheatmap_", cond)]] })
#rv$pheatmap <- pheatmap(id2,fontsize=8, angle_col=0, cutree_rows=2, cutree_col=2, scale="row", cluster_cols=F, treeheight_col=5, treeheight_row=15, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), annotation=annotation, width=6, height=12, show_rownames=F)

v0 <- v0 %>%
  mutate(Name = rownames(v0)) %>%
  select(Name, everything())  %>% 
  select(-c(genelabels, change))
rownames(v0) <- NULL

counts = as.data.frame(y$counts)
counts <- counts %>% mutate(Name = rownames(counts))
rownames(counts) <- NULL

v0 <- v0 %>% 
  left_join(counts)

#this to happen only when i load new. If on the precompiled , there is no need to do
if(create_msi_df){
  for(col in c("logCPM", "logFC", "LR", "PValue", "PValue", "FDR", "genelabels", "change")){
    colnames(v0) = rename_msi_df_columns(colnames(v0), col, paste(cond, col, sep = "."))
  }
  if(dim(rv$msi_df)[1] == 0){
    rv$msi_df = v0
  }else{
    rv$msi_df <- rv$msi_df %>%
      full_join(v0)
  }
  
}
print(paste("Runnign condition contrast END", cond))