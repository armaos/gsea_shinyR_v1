library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(DT)
#library(plotly)
library(shinyjs)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(ggupset )
library(stringr)
#
library("edgeR")
library(ggrepel) 
library(statmod)
library(pheatmap)
#library(RUVSeq)

if (!interactive()) sink(stderr(), type = "output")
wd = "/srv/shiny-server/gsea_shinyR/"



see_pathview <- function(..., save_image = FALSE)
{
  arguments <- list(...)
  pathview(gene.data = arguments$gene.data,
           gene.idtype = arguments$gene.idtype,
           pathway.id = arguments$pathway.id,
           species = arguments$species,
           bins = arguments$bins,
           limit = arguments$limit,
           low = arguments$low,
           mid = arguments$mid,
           high = arguments$high,
           na.col = arguments$na.col,
           kegg.native = arguments$kegg.native,
           kegg.dir = arguments$kegg.dir,
           same.layer = arguments$same.layer,
           sign.pos = arguments$sign.pos,
           new.signature =  arguments$new.signature ,
           out.suffix	= arguments$out.suffix
           
               
  )
  
  msg <- capture.output(pathview::pathview(...), type = "message")
  
  msg <- grep("image file", msg, value = T)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  file.remove(filename)
  img <- png::readPNG(filename)
  grid::grid.raster(img)
  
}

empty_image <- function(..., save_image = FALSE)
{
  print("empty")
  grid::grid.raster(png::readPNG("empty.png"))
}




server <- function(input, output, session) {
  
  # the annotation table 
  setwd(wd)
  #annotations_default = as.data.frame(read.delim2("Ka.3.1_annotation_table_test.tsv", header = T, sep = "\t", fill = T))
  annotations_default = as.data.frame(read.delim2("Ka.2_plus_3.2_annotation_table.tsv", header = T, sep = "\t", fill = T,  quote = ""))
  #annotations_default = as.data.frame(read.delim2("Ka.2_plus_3.2_annotation_table_test.tsv", header = T, sep = "\t", fill = T,  quote = ""))
  #annotations_default = as.data.frame(read.delim2("Ka.2_plus_3.1_annotation_table.tsv", header = T, sep = "\t", fill = T))
  
  representatives = annotations_default %>% select(Name)
  annotations_default <- annotations_default %>% 
    select(colnames(annotations_default)[!grepl("exp_", colnames(annotations_default))]) %>%
    select(-any_of(c("protein_product", "seqid", "start", "end", "strand"))) %>%
    mutate(Similar_to = ifelse(Similar_to == "-", Name,  paste(Similar_to, Name, sep = ","))) %>% 
    separate_rows(Similar_to, sep = "[,]") %>% 
    mutate(Name = ifelse(Similar_to == "-", Name, Similar_to)) %>%
    select(-Similar_to) %>%
    select_if(function(x) !(all(is.na(x)) | all(x=="-")))
  
  annotations_columns = sort(colnames(annotations_default))
  print(paste("dim annotations_default", dim(annotations_default), collapse = " "))
  
  
  
  
  
  # here AnnotationDbi::Ontology("GO:0006486") for go ontology grouping
  ##go2gene_all <- as.data.frame(read.table("kal_v3_all_GO.txt", header = T))
  # go2gene_all <- annotations_default %>%
  #   select(Name, Ontology_term , em_GOs, b2go_GO_id, ips_GO, arabThal_GO_ID) %>%
  #   gather(key, val, c("Ontology_term" , "em_GOs", "b2go_GO_id", "ips_GO", "arabThal_GO_ID")) %>%
  #   filter(val != "-") %>%
  #   select(val, Name) %>%
  #   separate_rows(val, sep = ",") %>%
  #   unique
  # go2gene_all %>% write.table("Ka.2_plus_3.1_GO_terms_manipulate.txt")
  go2gene_all <- as.data.frame(read.table("Ka.2_plus_3.1_GO_terms_manipulate.txt", header = T))
  
  go_terms =  rbind(
    as.data.frame(read.delim2("GOIDs_GOterms.txt", header = T, sep = "\t", stringsAsFactors = F, fill = T)) %>%
      mutate(ArabT = "0"),
    as.data.frame(read.delim2("GOIDs_GOterms_ArabT.txt", header = F, sep = " ", stringsAsFactors = F, fill = T)) %>%
      mutate(ArabT = "1") %>%
      rename(GO_ID = V2) %>%
      separate(V1, c("GO_namespace", "GO_term"), sep = "_AT_") %>%
      mutate(GO_namespace = ifelse(GO_namespace == "GO_MF", "molecular_function", ifelse(
        GO_namespace == "GO_BP", "biological_process", "cellular_component" ) ) ) %>%
      mutate(GO_term = str_replace_all(GO_term, "-", " ")) %>%
      mutate(GO_term = str_replace_all(GO_term, "_", " ")) %>%
      mutate(GO_term = str_to_lower(GO_term))
  )
  
  
  demo_genes <- as.data.frame(read.table("genes_to_test.txt", header = T))
  
  # these are the KEGG thresholds per node coming from the prepare_graph_nodes_table.py scirpt that reads all the nodes from the paths.
  # koala_thresholds = as.data.frame(read.table("koala_genes_vector.txt"))
  # koala_thresholds <- koala_thresholds %>%
  #   group_by(V1 ) %>% 
  #   summarise(V2 = max(V2))
  
  #### read KO annotations (i have manipulated the data to load faster)
  ## ko_annots <- as.data.frame(read.table("Ka.2_plus_3.1_KO_terms.txt", header = T, sep = "\t"))
  # ko_annots <- annotations_default %>% select(Name, koala_KEGG_ko , em_KEGG_ko) # this is better since i get the annotations directly from the table
  # ko_annots = full_join(
  #   ko_annots %>%
  #     select(Name, koala_KEGG_ko) %>%
  #     filter(koala_KEGG_ko != "-") %>%
  #     separate_rows(koala_KEGG_ko, sep = "[,]") %>%
  #     separate(koala_KEGG_ko, c("KO", "koala_threshold"), sep = " ") %>%
  #     mutate(koala_threshold = str_remove(koala_threshold, "[(]"),
  #            koala_threshold = str_remove(koala_threshold, "[)]"),
  #            koala_threshold = str_remove(koala_threshold, "[%]"),
  #            koala_threshold = as.numeric(koala_threshold)),
  #   ko_annots %>%
  #     select(Name, em_KEGG_ko) %>%
  #     filter(em_KEGG_ko != "-") %>%
  #     separate_rows(em_KEGG_ko, sep = "[,]") %>%
  #     mutate(em_KEGG_ko = str_remove(em_KEGG_ko, "ko:"),
  #            em = 1),
  #   by = c("Name" = "Name", "KO" = "em_KEGG_ko")
  # ) %>% unique
  # ko_annots %>% write.table("Ka.2_plus_3.1_KO_terms_manipulate.txt", quote = F, row.names = F)

  ko_annots <- as.data.frame(read.table("Ka.2_plus_3.1_KO_terms_manipulate.txt", header = T)) 
  koala_thresholds <- suppressWarnings(ko_annots %>%
    group_by(KO) %>%
    summarise(koala_threshold = max(koala_threshold, na.rm = T),
              em = max(em)) %>%
    mutate(koala_threshold = ifelse(is.infinite(koala_threshold), NA, koala_threshold)))
    
  
  ko_symbols <- as.data.frame(read.delim2("KO_titles_new.txt", header = F, sep = ";", quote = "")) %>%
    mutate(V3 = str_remove(V3, " ")) %>%
    rename(Symbol = V2, Node_Description = V3)
  
  
  #### read kegg mapp annotations (i have manipulated the data to load faster
  ## ko_Pathway = "Ka.2_plus_3.1_KO_Pathways.txt"
  ## ko_Pathway <- as.data.frame(read.table(ko_Pathway, header = T))
  # ko_Pathway <- annotations_default %>% select(Name, em_KEGG_Pathway , koala_KEGG_Pathway) # this is better since i get the annotations directly from the table
  # ko_Pathway <- ko_Pathway %>%
  #   filter(em_KEGG_Pathway != "-" | koala_KEGG_Pathway != "-" ) %>%
  #   gather(key, val, c(em_KEGG_Pathway, koala_KEGG_Pathway)) %>%
  #   separate_rows(val, sep = "[,]") %>%
  #   filter(val != "-" & grepl("ko", val)) %>%
  #   select(-key) %>%
  #   unique
  # ko_Pathway %>% mutate(val = str_replace(val ,"ko", "")) %>% write.table("Ka.2_plus_3.1_KO_Pathways_manipulate.txt", quote = F, row.names = F)
  ko_Pathway <- as.data.frame(read.table("Ka.2_plus_3.1_KO_Pathways_manipulate.txt", header = T, colClasses = c("character", "character"))) %>%
    dplyr::select(val, Name)  %>% 
    mutate(val = str_replace(val ,"ko", "")) 
  
  ######
  
  ko_Pathway_terms <- as.data.frame(read.table("Ka_KEGGpathways.txt", header = F, sep = "\t", stringsAsFactors = F, 
                                               colClasses = c("character", "character")))
  
  ko_nodes <- as.data.frame(read.table("KEGGmaps_and_nodes.txt", header = F, sep = "\t") ) %>% 
    rename(map = V1, nodes = V2) %>%
    mutate(map = str_replace(map, "ko", ""))
  
  
  
  rv <- reactiveValues(gsea_plots = NULL,
                       IP = NULL, 
                       show_norm_and_de_mout = NULL,
                       show_gso_mout = NULL,
                       show_KeggScape_mout = NULL,
                       show_gsea_mout = NULL,
                       show_gsea_plots_mout = NULL,
                       show_ora_mout = NULL,
                       
                       condition = "",
                       
                       update_annotations = 0,
                       update_msi_df = 0,
                       
                       norm_and_de_signal =0,
                       data = NULL, 
                       text = NULL, 
                       okplot = FALSE, 
                       #okplot_KEGGmap_ora = FALSE,
                       okplot_KEGGmap_gso = FALSE,
                       okplot_KEGGmap_gsea = FALSE,
                       okplot_gseaplot = FALSE,
                       map = "",
                       map_gsea = "",
                       go_term = "",
                       run_counter = 0, 
                       
                       gse = NULL,
                       gse_desc = NULL,
                       
                       barplot1 = NULL,
                       barplot2 = NULL,
                       dotplot = NULL,
                       emapplot = NULL,
                       cnetplot = NULL,
                       heatplot = NULL,
                       treeplot = NULL,
                       upsetplot = NULL,
                       
                       dotplot_gsea = NULL,
                       ridgeplot_gsea = NULL,
                       
                       filtered = data.frame(), 
                       pheatmap = NULL,
                       volcano0 = NULL,
                       pca_plot = NULL,
                       
                       gsea = NULL,
                       gseaplot = NULL,
                       
                       enriched_kegg_pathways = NULL,
                       #enriched_kegg_pathways_desc = NULL,
                       #enriched_kegg_pathways_id = NULL,
                       
                       input_data = NULL,
                       input_genes_with_annotations = c(),
                       
                       kegg_ora_vector = NULL,
                       kegg_gsea_vector = NULL,
                       go_term = NULL,
                       
                       ora = data.frame(), 
                       msi_df = data.frame(), 
                       gsea_input = NULL,
                       gsea = data.frame(),
                       gsea_kegg = data.frame(),
                       go = NULL,
                       enricher_plots = NULL,
                       
                       fc = c(),
                       
                       dim_msi_df = NULL ,
                       
                       kegg_gso_col_low = "red",
                       kegg_gso_col_high = "green",
                       kegg_gsea_col_low =  "red",
                       kegg_gsea_col_high = "green"
                       
  )
  
  html("overview_text", "Initialization complete.")
  shinyjs::show("data_import_box")
  # button to fill in the demo genes
  observeEvent(input$demo, {
    updateTextAreaInput(session, "gso_text_input",
                        value = paste(head(unique(demo_genes$Name) , 40), collapse = "\n") )
    updateSelectInput(session, "choose_gso", selected = "Custom Set")

  })

  #dropdown for selecting text or file for custom genes
  observeEvent(input$choose_custom_gso, {
    if(input$choose_custom_gso == "File"){
      if(is.null(input$gso_file_input$datapath) || input$gso_file_input$datapath == ""  ){
        shinyjs::disable("submit_custom_gene_set")
      }else{
        shinyjs::enable("submit_custom_gene_set")
      }
    }else if(input$choose_custom_gso == "text"){
      if(is.null(input$gso_text_input) || input$gso_text_input == ""){
        shinyjs::disable("submit_custom_gene_set")
      }else{
        shinyjs::enable("submit_custom_gene_set")
      }
    }
  })

  #enabling the submition with text for custom gene sets
  observeEvent(input$gso_text_input, {
    if(input$choose_custom_gso == "Text"){
      if(!is.null(input$gso_text_input) && input$gso_text_input != ""){
        shinyjs::enable("submit_custom_gene_set")
      }else{
        shinyjs::disable("submit_custom_gene_set")
      }
    }
  })
  
  #loading the count data for custom contrast
  # observeEvent(input$data_file$datapath, {
  #   if(input$choose == "file"){
  #     if(!is.null(input$data_file$datapath) && input$data_file$datapath != "" ){
  #       print("file loaded")
  #       shinyjs::enable("submit_norm_and_de")
  #     }else{
  #       shinyjs::disable("submit_norm_and_de")
  #       
  #     }
  #   }
  # })
  
  # this returns the filtered table from the tcpm input table
  observeEvent(input$submit_norm_and_de, {
    print("eventReactive(submit_norm_and_de")
    shinyjs::show("overview_text")
    html("overview_text", "Reading Count Data...")
    
    print("Reading...")
    data0 <- switch(input$choose_custom_tcmp,
                    testdata = read.table(file = '../Meiks_analysis/tcpm.csv', sep = ',', header = TRUE) %>% select(-"Undetermined_S0"),
                    testdata = read.table(file = '../Meiks_analysis/tcpm.csv', sep = ',', header = TRUE) %>% select(-"Undetermined_S0")
    )
    colnames(data0)[colnames(data0)=="gene_name"]="GeneSymbol"
    
    #data<- head(data0,1000)
    data <- data0
    print(paste("data", head(data, 1)))
    colnames(data)[colnames(data)=="Geneid"]="Geneid"
    
    
    suppressWarnings(rowMeans <- apply(data,1, function(x) mean(as.numeric(x),na.rm=T)))
    c <- data[!duplicated(data$Geneid),]
    rownames(c) <- c$Geneid
    ########################################
    html("overview_text", "Processing Count Data...")
    countData.all<-as.matrix(c[2:ncol(c)])
    
    rownames(countData.all)<-c$Geneid
    
    countData<-countData.all

    html("overview_text", "Filtering Count Data...")
    print("Filtering Data")
    filter<-apply(countData,1,function(x) length(x[x>10])>=ncol(countData)/5)   # 20967 gene features, not always 5, should be (ncol/4) or /5
    filtered<-countData[filter, ]
    
    donor <- factor(gsub("\\.\\w+$","",colnames(filtered),perl=T))
    group <- factor(gsub("^\\w+\\.","",colnames(filtered),perl=T))
    choices = subset(as.data.frame(t(combn(levels(factor(group)), 2))) %>% 
                       mutate(contrast = paste0(V1, "-", V2)))$contrast
    print(choices)
    updateSelectInput(session, "choose_conditions", choices = choices, selected = choices[1])
    
    
    html("overview_text", "Extracting conditions : DONE")
    print("Extracting conditions")
    shinyjs::show("menu_de_box")
    rv$filtered = filtered
    print("eventReactive(submit_norm_and_de END")
  })
  
  
  observeEvent(input$run_DE, {

    # shinyjs::show("overview_text")
    filtered = rv$filtered
    
    rv$condition  <- input$choose_conditions
    shinyjs::html("pageHeader", rv$condition)
    
    
    print("start DE")
    donor <- factor(gsub("\\.\\w+$","",colnames(filtered),perl=T))
    group <- factor(gsub("^\\w+\\.","",colnames(filtered),perl=T))
    x<-group
    input
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
    eval(parse(text = cmd))
    
    
    #my.contrasts = makeContrasts(cond, levels=design)              # need to rewrite acccording to the task
    print("makeContrasts OK")
    
    y.Disp <- estimateDisp(y, design, robust = TRUE)
    
    print("estimateDisp OK")
    html("overview_text", "Processing Constrasts: Fitting glmFit...")
    fit <- glmFit(y.Disp, design, robust=TRUE)
    print("Fitting glmfit")
    
    lfcT=1.5
    fdrT=0.01
    
    html("overview_text", "Processing Constrasts: Fitting glmLRT...")
    print("glmLRT")
    lrt.D_O <- glmLRT(fit, contrast=my.contrasts)
    # lrt.aDTP_WT <- glmQLFTest(fit, contrast=my.contrasts[,id])  # quasi-likelihood (QL)
    
    # add FDR
    deg.edger <- lrt.D_O$table[p.adjust(lrt.D_O$table$PValue, method = "BH") < fdrT, ]
    #dim(deg.edger)
    
    #################################################################
    # final volcano plot           ##################################
    print("Proc volcano Data")
    html("overview_text", "Processing Constrasts: Processing volcano data...")
    
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
    html("overview_text", "Processing Constrasts: Processing volcano plot...") 
    
    rv$volcano0 <- ggplot(data=v0, aes(x=logFC, y=-log10(FDR), color=change, fill=change), fontface='bold') +
      geom_point(alpha=0.6, size=1) + 
      theme_bw(base_size=15) + 
      theme(legend.position='none',panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
      ggtitle(paste0("Volcanoplot of ",input$choose_conditions)) + 
      scale_fill_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
      scale_color_manual(name="", values=c("red", "darkorange", "darkgrey"), limits=c("Up", "Down", "NotSig")) + 
      geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + geom_hline(yintercept=-log10(fdrT), lty=2, col="gray", lwd=0.5)
    
    
    #################################################################
    
    html("overview_text", "Processing Constrasts: Processing heatmap plot...") 
    print("Heatmap plot")
    id1 <- log2(cpm(y)+1)
    id1 <- y
    selectedCols <-c(grep("\\.O",colnames(y)),grep("\\.D",colnames(y)))  #select some columns
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
    html("overview_text", "DE analysis: DONE...") 
    
    rv$show_norm_and_de_mout = T
  })
  
  
  # add DE values to the annotations table
  annotations <- eventReactive(
    rv$update_annotations, {
      if(rv$update_annotations == 0){
        data.frame()
      }else{
        
        annotations_default %>% 
          select(all_of(annotations_columns)) %>% 
          full_join(rv$msi_df) %>% 
          mutate(In_DE_data = ifelse(!is.na(logFC), "Yes", "No"))
      }
      
  })
  
  
  observeEvent(input$submit_contrast, {
    #shinyjs::disable("submit_contrast")
    #shinyjs::disable("submit_norm_and_de")
    #print("submit_norm_and_de")
    shinyjs::show("submitted_contrast1")
    #shinyjs::show("submitted_contrast2")
    rv$condition  <- input$Condition
    
    shinyjs::html("pageHeader", rv$condition)
    
    rv$msi_df = as.data.frame(read.table(paste0("MSI_data/", rv$condition, "_all_values.csv"), header = T , sep = "\t"))
    colnames(rv$msi_df)[1] = "Name"
    rv$msi_df <- rv$msi_df %>% select(colnames(rv$msi_df)[!grepl("X", colnames(rv$msi_df))])
    
    # add DE values to the annotations table
    rv$update_annotations =  rv$update_annotations + 1
    rv$update_msi_df =  rv$update_msi_df + 1
  })
  
  observeEvent(
    rv$update_msi_df, {
      if(rv$update_msi_df > 0){
        
        rv$msi_df <- rv$msi_df %>% arrange(desc(logFC))
      # annotations <- annotations %>% select(all_of(annotations_columns))
      # annotations <- annotations %>% full_join(rv$msi_df) %>% mutate(In_DE_data = ifelse(!is.na(logFC), "Yes", "No"))
      # print(paste("dim annotations in submit contrast", dim(annotations), collapse = " "))
      # 
      # prepare gsea input
      rv$gsea_input = rv$msi_df$logFC
      rv$gsea_kegg = data.frame()
      names(rv$gsea_input) <- rv$msi_df$Name
      
      # summarise the DE table to reduce redundancy
      print(paste("msi df dim before group.", dim(rv$msi_df)))
      rv$msi_df <- rv$msi_df %>% group_by(across(c(-Name))) %>% summarise(Name = paste(Name, sep = ",", collapse = ",")) %>% ungroup() 
      rv$msi_df <- rv$msi_df %>% select(Name, everything())
      print(paste("msi df dim After group.", dim(rv$msi_df)))
      # add annotations to the DE table
      rv$msi_df = rv$msi_df %>% mutate(primary_name = Name) %>% separate(primary_name, c("primary_name"), ",") %>% select(Name, everything() ) #%>% select(-condition, )
      rv$msi_df <- rv$msi_df %>% left_join(annotations_default, by = c("primary_name" = "Name"))
      
      #rv$msi_df = as.data.frame(read.table(paste0("MSI_Kal_t_18_june22_annotated/", rv$condition, "_all_values.csv"), header = T , sep = "\t", quote = ""))
      #rv$msi_df = rv$msi_df %>% mutate(primary_name = Name) %>% separate(primary_name, c("primary_name"), ",") %>% select(Name, everything() ) #%>% select(-condition, )
      #
      print("OK msi df")
     
      
      default_columns = c("In_DE_data", colnames(rv$msi_df)) #[1:which(colnames(rv$msi_df) == "Ontology_term") - 1]
      #annotations_columns = colnames(rv$msi_df)[which(colnames(rv$msi_df) == "Ontology_term" ) : ( length(colnames(rv$msi_df) ) -1 )] 
      updatePickerInput(session, "overview_table_Default_values", choices = default_columns, selected = c("Name", "In_DE_data", "logFC", "logCPM", "FDR"))
      #rv$msi_df <- rv$msi_df %>% left_join(annotations)
      updatePickerInput(session, "overview_table_terms", choices = annotations_columns[annotations_columns != "Name"])
      
      
      shinyjs::show("run_gsea_box")
      shinyjs::show("data_overview_box")
      shinyjs::show("data_overview_table_box")
      
      print("msi_df and contrast load succesfully!")
      p(id = "overview_text", "Data loaded succesfully. You can head on to the downstream analysis on the left panel.")
      rv$show_gso_mout = T
      rv$show_gsea_mout = T
    }
  })
  
  
  ##### render dynamic menus
  output$menu_norm_and_de_mout <- renderMenu({
    if(!is.null(rv$show_norm_and_de_mout ))
      menuItem("Normalization and DE", 
               menuSubItem("Plots", tabName = "menu_norm_and_de"), 
               icon = icon("th"))
  })
 
  output$gso_mout <- renderMenu({
    if(!is.null(rv$show_gso_mout))
      menuItem("Gene set Analysis",
               menuSubItem("Gene set selection", tabName = "set_selection"),
               # menuItem("Gene set Analytics",
               #          menuSubItem("Plots", tabName = "gso_plots"),
               #          #menuSubItem("KEGGscape", tabName = "KEGGscape_gso"),
               #          icon = icon("th")),
               menuItemOutput("ora_mout"),
               icon = icon("th"))
  })
  
  output$ora_mout <- renderMenu({
    if(!is.null(rv$show_ora_mout )){
      menuItem("Enrichment Analysis (ORA)",
               menuSubItem("Configuration And Table", tabName = "configure_ora"),
               menuSubItem("Plots", tabName = "plots_ora"),
               icon = icon("th"))
      }
    })
  
  
  output$gsea_mout <- renderMenu({
    if(!is.null(rv$show_gsea_mout )){
      menuItem("GSEA",
               menuSubItem("Configuration And Table", tabName = "table_gsea"),
               menuItemOutput("gsea_plots_mout"),
               icon = icon("th")
               )
    }
  })
  
  output$gsea_plots_mout <- renderMenu({
    if(!is.null(rv$show_gsea_plots_mout)){
      menuItem("GSEA Plots",
               menuSubItem("Plots", tabName = "plots_gsea"),
               menuSubItem("Running Score", tabName = "running_score_gsea"),
               menuSubItem("Custom UpSet", tabName = "custom_upset_gsea"),
               icon = icon("th"))
    }
  })
  
  output$KeggScape_mout <- renderMenu({
    if(!is.null(rv$show_KeggScape_mout))
      menuItem("KEGGscape", tabName = "KEGGscape_gso", icon = icon("th"))
  })
  
  ######## end render dynamic menus
  
  # observeEvent(input$submit_norm_and_de, {
  #   rv$condition  <- input$Condition
  #   rv$norm_and_de_signal =  rv$norm_and_de_signal + 1 
  #   print("observeEvent(input$submit_norm_and_de")
  #   #shinyjs::disable("submit_contrast")
  #   #shinyjs::disable("submit_norm_and_de")
  # })
  
  
  #picker for KEGGscape in GSO
  observeEvent(input$choose_map_gso, {
    rv$map <- subset(ko_Pathway_terms %>% filter(V2 == input$choose_map_gso))$V1 #str_replace(input$choose_map_ora, "ko", "")
    #print(paste0("observeEvent choose_map_ora rv$map-", rv$map,"-", is.null(rv$map),"-", is.na(rv$map), "-",length(rv$map)))
  })
  
  #picker for KEGGscape in ORA
  observeEvent(input$choose_map_ora, {
    rv$map <- subset(ko_Pathway_terms %>% filter(V2 == input$choose_map_ora))$V1 #str_replace(input$choose_map_ora, "ko", "")
    #print(paste0("observeEvent choose_map_ora rv$map-", rv$map,"-", is.null(rv$map),"-", is.na(rv$map), "-",length(rv$map)))
  })
  
  #picker for KEGGscape in GSEA
  observeEvent(input$choose_map_gsea, {
    rv$map_gsea <- subset(ko_Pathway_terms %>% filter(V2 == input$choose_map_gsea))$V1 
    #print(paste0("observeEvent choose_map_gsea rv$map_gsea-", rv$map_gsea,"-", is.null(rv$map_gsea),"-", is.na(rv$map_gsea), "-",length(rv$map_gsea)))
  })
  
  
  observeEvent(input$choose_GO, {
    # rv$okplot_gseaplot <- TRUE
    rv$go_term <- input$choose_GO
    #rv$gseaplot = gseaplot(gsea, by = "all", title = rv$go_term, geneSetID = 1)
    #print(paste("observeEvent choose_GO rv$go_term",rv$okplot_gseaplot , rv$go_term))
  })
  
  
  
  #if i define gene set by filters
  observeEvent(
    eventExpr = input$run, {
      # this is the reading of MSI data
      fc_filter_low = input$slider_lfc[1] 
      fc_filter_high = input$slider_lfc[2] 
      
      #pval_filter = input$slider_pval[1] 
      fdr_filter = input$slider_fdr[1] 
      
      msi_df <- rv$msi_df %>%  
        filter(abs(logFC) >= fc_filter_low & abs(logFC) <= fc_filter_high & FDR <= 10^(-fdr_filter)) %>%
        select(Name, logFC) %>%
        separate_rows(Name, sep = ",")
      
      rv$input_data <- msi_df$Name
      rv$fc = msi_df$logFC
      names(rv$fc)= rv$input_data
      print(paste("filtered genes load succesfully", length(rv$input_data)[1]))
    })
  
  #if i define gene set by input
  observeEvent(
    eventExpr = input$submit_custom_gene_set, {
      
      input_genes = switch(input$choose_custom_gso,
                       File = read_file(input$gso_file_input$datapath),
                       Text = input$gso_text_input
      )
      
      input_genes = str_split(input_genes, "\n")[[1]] 
      msi_df <- rv$msi_df %>%  
        separate_rows(Name, sep = ",") %>%
        filter(Name %in% input_genes) %>%
        select(Name, logFC) 
      
      rv$input_data <- input_genes
      
      rv$fc = msi_df$logFC
      names(rv$fc)= msi_df$Name
      print(paste("custom genes load succesfully", length(rv$input_data)))
    })
  
  #when i have defined the input data either by filters or by input
  observeEvent(
    eventExpr = rv$input_data, {
      
      if (rv$condition == ""){
        shinyjs::show("text_configure_ora")
      }else{
        shinyjs::hide("text_configure_ora")
        #shinyjs::disable("run")
        
        html("text_ora_plots", "You should submit your ploting configuration first and Plot your table!!")
        shinyjs::show("text_ora_plots")
        
        # html("text_KEGGscape_ORA", "Processing...")
        # shinyjs::show("text_KEGGscape_ORA")
        html("text_KEGGscape_gso", "Processing...")
        shinyjs::show("text_KEGGscape_gso")
        
        # rv$data = switch(input$choose,
        #                  file = read_file(input$data_file$datapath),
        #                  text = input$input_text
        # ) 
        print(paste("dim annotations() in rv$input_data", dim(annotations()), collapse = " "))
        
        shinyjs::show("summary_gene_set_box")
        # rv$input_genes_with_annotations = annotations() %>%
        #   filter(Name %in% rv$input_data) %>%
        #   filter_at(annotations_columns[annotations_columns != "Name" & annotations_columns != "mcl_cluster"], any_vars(. != "-")) %>%
        #   distinct(Name)
        
        # print(paste0("N genes: ", length(rv$input_data), 
        #              "<br/>",
        #              "N annotated genes: ", dim(rv$input_genes_with_annotations)[1]))
        
        output$summary_gene_set <- renderUI({ 
          HTML(paste0("N genes: ", length(rv$input_data)
                      #"<br/>",
                      #"N annotated genes: ", dim(rv$input_genes_with_annotations)[1]
                      ))
          })
        
        
        rv$IP = input$getIP
        rv$run_counter = rv$run_counter +1
        

        
        rv$map <- ""
        rv$go_term <- ""
        
        
        
        #input_data = str_split(rv$data, "\n")[[1]] 
        #print(paste("demo inut data", length(input_data)))
        
        ### this is the diana process
        #genes_to_GO = go2gene_all %>% select(Name)
        # rv$go <- enrich_go_custom(
        #   input_data,
        #   universe = genes_to_GO[, 1],
        #   genes_to_GO,
        #   qvalue = 0.1,
        #   pvalue = 0.05,
        #   GO_type = "BP"
        # )
        
       
        print(paste("MSI input data", length(rv$input_data)))
        rv$dim_msi_df <- length(rv$input_data)
       
       
        
        ## UPDATE the GSO kegg map picker
        gso_kegg_maps = ko_Pathway %>% 
          filter(Name %in% rv$input_data) %>% 
          distinct(val) %>% 
          left_join(ko_Pathway_terms, by = c("val" = "V1")) %>%
          arrange(V2) %>%
          filter(!is.na(V2))
        updateSelectInput(session, "choose_map_gso", choices = gso_kegg_maps$V2)
        rv$map <- ifelse(length(gso_kegg_maps$V2) > 0 , gso_kegg_maps[1,2], rv$map)
        
        shinyjs::show("gso_kegg_pathway_dropdown")
        shinyjs::show("gso_keggscape_table_box")
        
        
        # find the ko node id terms
        # print("FIND ko_annots")
        # gso_kegg <- ko_annots %>%
        #   filter(Name %in% rv$input_data) %>%
        #   left_join(rv$msi_df %>% select(Name, logFC) %>% separate_rows(Name, sep = ",") ) %>%
        #   group_by(val) %>%
        #   summarise(logFC = max(logFC, na.rm = T)) %>%
        #   ungroup()
        # 
        # rv$kegg_gso_vector <- gso_kegg$logFC
        # names(rv$kegg_gso_vector) <- gso_kegg$val
        # print(paste("kegg_gso_vector",dim(rv$kegg_gso_vector)[1]))
        
        #rv$okplot <- TRUE
        #rv$okplot_KEGGmap_ora <- TRUE 
        rv$okplot_KEGGmap_gso <- TRUE 
        
        print("and now show!")
        shinyjs::show("box_configure_term_ora")
        shinyjs::show("enrichment_term")
        shinyjs::show("submit_enrichment")
        #shinyjs::hide("text_KEGGscape_ORA")
        shinyjs::hide("text_KEGGscape_gso")
        #shinyjs::show("ora_keggscape_table_box")
        
        rv$show_KeggScape_mout = T
        rv$show_ora_mout = T
      }
    })
  
  
  observeEvent(
    eventExpr = input$run_gsea, {
      if (rv$condition == ""){
        shinyjs::show("text_configure_gsea")
      }else{
        shinyjs::hide("text_configure_gsea")
        
        
        html("text_gsea_plot_overview", "Processing...")
        shinyjs::show("text_gsea_plot_overview")
        
        html("text_gsea_plot", "Processing...")
        shinyjs::show("text_gsea_plot")
        
        html("text_gsea_table", "Processing...")
        shinyjs::show("text_gsea_table")
        
        html("text_KEGGscape_GSEA", "Processing...")
        shinyjs::show("text_KEGGscape_GSEA")
        
        ######
        
        rv$ridgeplot_gsea <- NULL
        rv$dotplot_gsea <- NULL
        rv$gsea <- data.frame()
        rv$okplot_KEGGmap_gsea <- FALSE
        rv$map_gsea <- ""
        rv$okplot_gseaplot <-  FALSE
        
        #####
        gsea_kegg_v = ko_annots %>% 
          inner_join(rv$msi_df %>% select(Name, logFC) %>% separate_rows(Name,sep = ",")) %>%
          #inner_join(msi_df %>% select(Name, logFC) %>% separate_rows(Name,sep = ",")) %>%
          group_by(KO) %>%
          summarise(logFC = max(logFC, na.rm = T)) %>%
          ungroup()
        rv$kegg_gsea_vector = gsea_kegg_v$logFC
        names(rv$kegg_gsea_vector) <- gsea_kegg_v$KO
        
        print("GSEA")
        # first run the KEGG gsea just to have it stored so that the maps can be loaded
        if(dim(rv$gsea_kegg)[1] == 0){
          rv$gsea_kegg = GSEA(rv$gsea_input, TERM2GENE=ko_Pathway, TERM2NAME = ko_Pathway_terms)
          
          updateSelectInput(session, "choose_map_gsea", choices = rv$gsea_kegg$Description[!is.na(rv$gsea_kegg$Description)])
          rv$map_gsea = ifelse(length(rv$gsea_kegg[, 1]) > 0 , rv$gsea_kegg[1, 1], rv$map_gsea)
          shinyjs::show("gsea_kegg_pathway_dropdown")
          print(paste("rv$map_gsea ", rv$map_gsea ))
        }
        
        print(paste("length rv$gsea_kegg" , length(rv$gsea_kegg)))
        
        if (input$enrichment_term_gsea %in% c("biological_process", "molecular_function", "cellular_component")){
          term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == input$enrichment_term_gsea & ArabT == "0"))$GO_ID)
          term2name = go_terms %>% filter(GO_namespace == input$enrichment_term_gsea & ArabT == "0")
          
          rv$gsea <- GSEA(rv$gsea_input, 
                          TERM2GENE = term2gene,
                          TERM2NAME = term2name)
          
        }else if (input$enrichment_term_gsea %in% c("biological_process_arth", "molecular_function_arth", "cellular_component_arth")){
          term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == str_remove(input$enrichment_term_gsea, "_arth") & ArabT == "1"))$GO_ID)
          term2name = go_terms %>% filter(GO_namespace == str_remove(input$enrichment_term_gsea, "_arth") & ArabT == "1")
          
          rv$gsea <- GSEA(rv$gsea_input, 
                          TERM2GENE = term2gene,
                          TERM2NAME = term2name)
          
        }else if(input$enrichment_term_gsea == "KEGG_Pathways"){
          term2gene = ko_Pathway
          term2name = ko_Pathway_terms
          rv$gsea = rv$gsea_kegg
        }
        
        
        
        
        
        rv$dotplot_gsea = dotplot(rv$gsea)
        print("doing ridgeplot_gsea")
        rv$ridgeplot_gsea = ridgeplot(rv$gsea )
        print("ridgeplot_gsea DONE")
        updateSelectInput(session, "choose_GO", choices = rv$gsea$Description[!is.na(rv$gsea$Description)])
        #updateCheckboxGroupInput(session, "custom_upset_terms", choices = rv$gsea$Description[!is.na(rv$gsea$Description)])
        updatePickerInput(session, "custom_upset_terms", choices = rv$gsea$Description[!is.na(rv$gsea$Description)])
        
        rv$show_gsea_plots_mout = T
        
        shinyjs::show("gsea_plots_tabsetPanel")
        
        shinyjs::show("custom_upset_box")
        
        shinyjs::hide("text_KEGGscape_GSEA")
        shinyjs::show("gsea_kegg_pathway_table")
        
        
        
        
        html("text_gsea_table", "GSEA: DONE")
        print("GSEA: DONE")
        
        shinyjs::show("gsea_table_box")
        rv$okplot_KEGGmap_gsea <- TRUE
        rv$okplot_gseaplot <- TRUE
      }
    })
  
  
  
  observeEvent(
    eventExpr = input$submit_enrichment, {
      print(" input$submit_enrichment")
      shinyjs::hide("box_ora_send_to_plot")
      shinyjs::hide("ora_table_box")
      #rv$okplot_KEGGmap_ora <- FALSE
      rv$barplot1 <- NULL
      rv$barplot2 <- NULL
      rv$dotplot <- NULL
      rv$emapplot <- NULL
      rv$cnetplot <- NULL
      rv$heatplot <- NULL
      rv$treeplot <- NULL
      rv$upsetplot <- NULL
      
      rv$ora <- data.frame()
      rv$ora_desc <- NULL
      
      html("text_ora_table", "Processing...")
      shinyjs::show("text_ora_table")
      print("ORA with enricher")
      
      universe = switch(input$choose_universe,
                       all_t = annotations() %>%
                         distinct(Name),
                       all_annot_t = annotations() %>%
                         filter_at(annotations_columns[annotations_columns != "Name" & annotations_columns != "mcl_cluster"], any_vars(. != "-")) %>%
                         distinct(Name),
                       all_de = rv$msi_df %>% select(Name) %>% separate_rows(Name, sep = ",") 
                       
      )
      # we should also add the input data to the universe
      universe <- rbind(universe, data.frame(Name = rv$input_data)) %>% distinct(Name)
      print(paste("universe:", dim(universe)[1]))
      
      if (input$enrichment_term %in% c("biological_process", "molecular_function", "cellular_component")){
        term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == input$enrichment_term & ArabT == "0"))$GO_ID) 
        term2name = go_terms %>% filter(GO_namespace == input$enrichment_term & ArabT == "0")
      
      }else if (input$enrichment_term %in% c("biological_process_arth", "molecular_function_arth", "cellular_component_arth")){
        term2gene = go2gene_all %>% filter(val %in% subset(go_terms %>% filter(GO_namespace == str_remove(input$enrichment_term, "_arth") & ArabT == "1"))$GO_ID)
        term2name = go_terms %>% filter(GO_namespace == str_remove(input$enrichment_term, "_arth") & ArabT == "1")

      }else if(input$enrichment_term == "KEGG_Pathways"){
        term2gene = ko_Pathway
        term2name = ko_Pathway_terms
      }
      
      rv$ora <- enricher(rv$input_data, 
                         universe = universe$Name,
                         TERM2GENE = term2gene,
                         TERM2NAME = term2name, 
                         minGSSize = input$ORA_minGSSize,
                         maxGSSize = input$ORA_maxGSSize)
      
      print("ORA with enricher :DONE")
      if(dim(data.frame(rv$ora))[1] > 0 ){
        shinyjs::show("box_ora_send_to_plot")
        shinyjs::show("ora_table_box")
        html("text_ora_table", "ORA completed!")
      }else{
        html("text_ora_table", paste0("No gene sets have size between ",input$ORA_minGSSize , " and ",  input$ORA_maxGSSize," 500 ..."))
      }

    })
  
  
  observeEvent(
    eventExpr = input$submit_plot_ora, {
      html("text_ora_plots", "Processing...!")
      
      print(paste("input$submit_plot_ora", dim(as.data.frame(rv$ora))[1]))
      if(length(rv$ora ) > 0){
        #if(!is.null(rv$ora) & rv$ora != ""){
        
        max_to_show = input$ORA_max_elements
        print(max_to_show)
        p.adjust_low = 0
        p.adjust_high = 1
        if(!is.null(input$ORA_p.adjust)){
          p.adjust_low = as.numeric(str_split(input$ORA_p.adjust, ' ')[[1]][1]) 
          p.adjust_low = ifelse(is.numeric(p.adjust_low), p.adjust_low, 0)
          
          p.adjust_high = as.numeric(str_split(input$ORA_p.adjust, ' ')[[1]][3])
          p.adjust_high = ifelse(is.numeric(p.adjust_high), p.adjust_high, 1)
          
        }
        
        print(paste(p.adjust_low , p.adjust_high))
        
        rv$ora_desc = head(
          subset(
            as.data.frame(rv$ora) %>%
              filter(p.adjust >= p.adjust_low & p.adjust <= p.adjust_high )
          )$Description,
          max_to_show )
        
        print(paste("length of ora desc" , length(rv$ora_desc)))
        if(length(rv$ora_desc) > 0 ){
          print("barplot1")
          rv$barplot1 = barplot(rv$ora, showCategory=rv$ora_desc, main = "Barplot (count)" )
          rv$barplot2 = mutate(rv$ora,
                               qscore = -log(p.adjust, base=10)) %>%
            barplot(x="qscore", showCategory=rv$ora_desc, main ="Barplot (qscore)")
          print("dotplot")
          rv$dotplot = dotplot(rv$ora, showCategory=rv$ora_desc)
          print("pairwise_termsim")
          p_t = pairwise_termsim(rv$ora)
          print("emapplot")
          rv$emapplot = emapplot(p_t, showCategory = rv$ora_desc)
          print("cnetplot")
          rv$cnetplot = cnetplot(rv$ora, categorySize="pvalue", showCategory = rv$ora_desc, foldChange = rv$fc )
          print("heatplot")
          rv$heatplot = heatplot(rv$ora, showCategory = rv$ora_desc, foldChange = rv$fc )
          print("upsetplot")
          rv$upsetplot = upsetplot(rv$ora, showCategory = rv$ora_desc)
          print("treeplot")
          #rv$treeplot = treeplot(p_t)
          print("DONE plots")
          shinyjs::hide("text_ora_plots")
          shinyjs::show("ora_plots_tabsetPanel")
          
          
        }else{
          html("text_ora_plots", "<strong>No enriched terms, try with a different threshold!</strong>")
        }
        
      }
    })
  
  
  # output$enriched_go <- renderPlotly({
  #   if(rv$okplot){
  #     shinyjs::enable("run")
  #     if (!is.null(rv$data ) && rv$data != ""){
  #       draw_enrich_go(rv$go, max_go = 20)
  #       shinyjs::hide("text1")
  #     }else{
  #       html("text1", "<strong>bold</strong> WRONG input data")
  #     }
  #   }
  # })
  
  
  # output$enriched_map <- renderPlot({
  #   if(rv$okplot){
  #     shinyjs::enable("run")
  #     if (!is.null(rv$data ) && rv$data != ""){
  #       draw_enrich_go_map(rv$go)
  #       shinyjs::hide("text1")
  #     }else{
  #       html("text1", "<strong>bold</strong> WRONG input data")
  #     }
  #   }
  # })
  
  
  # output$enriched_map <- renderPlot({
  #   if(rv$okplot){
  #     shinyjs::enable("run")
  #     if (!is.null(rv$data ) && rv$data != ""){
  #       
  #       hist(rnorm(100, 4, 1),breaks = 50)
  #       shinyjs::hide("text1")
  #       shinyjs::hide("text_gsea_diana")
  #     }else{
  #       html("text1", "<strong>bold</strong> WRONG input data")
  #       html("text_gsea_diana", "<strong>bold</strong> WRONG input data")
  #     }
  #   }
  # })
  
  #TABLES
  # output$DE_table <- DT::renderDT({
  #   as.data.frame(rv$y) %>% head(10)
  # },
  # filter="top",
  # options = list(
  #   scrollX = TRUE,
  #   dom = '<"top" lpif>',
  #   pageLength = 10)
  # )
  
  
  output$data_overview_table <- DT::renderDT({
    colnames_to_show = c(input$overview_table_Default_values,  input$overview_table_terms)
    print(paste("dim annotations() - render data_overview_table", dim(annotations()), collapse = " "))
    
    
    annotations() %>% 
      #filter(Name %in% rv$input_data) %>%
      select(colnames_to_show)
  },
  filter="top",
  options = list(
    scrollX = TRUE, 
    scrollY = "400px",
    dom = '<"top" lpift>')
  )
  
  
  output$ora_table <- DT::renderDT({
    if(length(rv$ora ) > 0){
      #if(!is.null(rv$ora) & rv$ora != ""){
      shinyjs::enable("run")
      shinyjs::hide("text_ora_table")
      shinyjs::show("submit_plot_ora")
      shinyjs::show("ORA_max_elements")

      
      as.data.frame(rv$ora) %>% select(-geneID) 
      #as.data.frame(rv$ora[rv$ora$pvalue <= 10^(-isolate(input$go_pval))])  %>% select(-geneID)
    }
    
  },
  callback = JS(
    # "table.on( 'search.dt', function () { Shiny.setInputValue( 'search', table.search() ); } );",
    "$( 'input').on( 'input', function(){",
    "var value = this.value;",
    "var clicked_td = $(this).closest('td');",
    "var td_index = clicked_td.index(); ",
    "var tr_index = clicked_td.parent().index();",
    #"console.log( 'search11', td_index , tr_index, table.column(td_index).header().textContent);",
    "Shiny.setInputValue( 'ORA_' + table.column( td_index ).header().textContent , value);",
    "});"
  ),
  filter="top",
  options = list(
    scrollX = TRUE, 
    scrollY = "400px",
    dom = '<"top" lpift>')
  )
  
  
  # output$ora_kegg_table <- DT::renderDT({
  #   if(!is.null(rv$enriched_kegg_pathways)){
  #     shinyjs::enable("run")
  #     shinyjs::hide("text_ora_kegg_table")
  #     shinyjs::show("kegg_pval")
  #     as.data.frame(rv$enriched_kegg_pathways[rv$enriched_kegg_pathways$pvalue <= 10^(-input$kegg_pval)])  %>% select(-geneID)}
  # },
  # filter="top",
  # options = list(
  #   scrollX = TRUE, 
  #   scrollY = "400px",
  #   dom = '<"top" lpift>')
  # )
  
  
  output$gsea_table <- DT::renderDT({
    if(length(rv$gsea ) > 0 ){
      shinyjs::enable("run")
      shinyjs::hide("text_gsea_table")
      
      
      as.data.frame(rv$gsea) %>% select(-core_enrichment)}
  },
  filter="top",
  options = list(
    scrollX = TRUE, 
    scrollY = "400px",
    dom = '<"top" lpift>')
  )
  
  
  # ORA plots  
  output$barplot1 <- renderPlot({
    if (!is.null(rv$barplot1) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_barplot1_text")
      rv$barplot1
    }else{
      empty_image()
    }
  })
  
  output$barplot2 <- renderPlot({
    if (!is.null(rv$barplot2) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_barplot2_text")      
      rv$barplot2
    }else{
      empty_image()
    }
  })
  
  output$dotplot <- renderPlot({
    if(!is.null(rv$dotplot) & !is.null(rv$ora_desc)){
      shinyjs::enable("run")
      shinyjs::hide("ora_dotplot_text")
      rv$dotplot
    }else{
      shinyjs::enable("run")
      html("text_ora_plots", "<strong>Make your configurarion for plotting from the previous tab panel first!</strong>")
      empty_image()
    }
  })
  
  output$emapplot <- renderPlot({
    if (!is.null(rv$emapplot) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_emapplot_text")
      
      rv$emapplot
    }else{
      empty_image()
    }
  })
  
  output$cnetplot <- renderPlot({
    if (!is.null(rv$cnetplot) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_cnetplot_text")
      
      rv$cnetplot
    }else{
      empty_image()
    }
  })
  output$heatplot <- renderPlot({
    if (!is.null(rv$heatplot) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_heatplot_text")
      
      rv$heatplot
    }else{
      empty_image()
    }
  })
  
  output$treeplot <- renderPlot({
    if (!is.null(rv$treeplot) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_treeplot_text")
      
      rv$treeplot
    }else{
      empty_image()
    }
  })
  
  
  output$upsetplot <- renderPlot({
    if (!is.null(rv$upsetplot) & !is.null(rv$ora_desc)){
      shinyjs::hide("text_ora_plots")
      shinyjs::enable("run")
      shinyjs::hide("ora_upsetplot_text")
      
      rv$upsetplot
    }else{
      empty_image()
    }
  })
  
  ### NORM plots
  output$pheatmap <- renderPlot({
    if (!is.null(rv$pheatmap) ){
      rv$pheatmap
    }else{
      empty_image()
    }
  })
  
  output$volcano0 <- renderPlot({
    if (!is.null(rv$volcano0) ){
      rv$volcano0
    }else{
      empty_image()
    }
  })
  
  output$pca_plot <- renderPlot({
    if (!is.null(rv$pca_plot) ){
      rv$pca_plot
    }else{
      empty_image()
    }
  })
  
  
  #GSEA overview
  output$dotplot_gsea <- renderPlot({
    if (!is.null(rv$dotplot_gsea) & !is.null(rv$gsea)){
      shinyjs::hide("text_gsea_plot_overview")
      shinyjs::enable("run")
      shinyjs::hide("gsea_dotplot_text")
      
      rv$dotplot_gsea
    }else{
      empty_image()
    }
  })
  
  output$ridgeplot_gsea <- renderPlot({
    if (!is.null(rv$ridgeplot_gsea) & !is.null(rv$gsea)){
      shinyjs::hide("text_gsea_plot_overview")
      shinyjs::enable("run")
      shinyjs::hide("gsea_ridgeplot_text")
      
      rv$ridgeplot_gsea
      
    }else{
      empty_image()
    }
  })
  # output$treeplot <- renderPlot({
  #   if (!is.null(rv$treeplot)){
  #     shinyjs::hide("text_ora_plots")
  #     rv$treeplot
  #   }
  # })
  
  #GSEA PLOT
  output$gseaplot <- renderPlot({
    if(rv$okplot_gseaplot){
      shinyjs::enable("run")
      #rv$gsea_plots[rv$go_term]
      shinyjs::hide("text_gsea_plot")
      shinyjs::show("choose_GO_box")
      # if ((!is.null(rv$data ) && rv$data != "" ) && (!is.null(rv$go_term ) && rv$go_term != ""  && rv$go_term != "NA" && !is.na(rv$go_term ) ) ){
      #   ggarrange(plotlist=gsea_plots, col = 1)
      #   print("in gseaplot")
      #   shinyjs::hide("text_gsea_plot")
      # }
      if( rv$go_term != "" && length(rv$go_term) > 0 && !is.na(rv$go_term )){  
        print(paste("output$gseaplot" ,rv$go_term ))
        print(rv$gsea[1,1])
        print(which(rv$gsea$Description == rv$go_term))
        gseaplot2(rv$gsea, geneSetID = which(rv$gsea$Description == rv$go_term) , title = rv$go_term)
      }
    }
  })
  
  
  
  my_custom_upset_plot <- eventReactive(input$submit_custom_upset_plot,{
    as.data.frame(isolate(rv$gsea)) %>% 
      filter(Description %in% isolate(input$custom_upset_terms)) %>%
      select(Description, core_enrichment) %>% 
      separate_rows(core_enrichment, sep = "/") %>% 
      mutate(Description = ifelse(str_length(Description) > 30, paste0(str_sub(Description,1,30), "..."), Description)) %>%
      group_by(core_enrichment) %>% 
      summarise(terms = list(Description)) %>% 
      ggplot(aes(x=terms)) +
      geom_bar() +
      scale_x_upset()
  })
  output$gsea_custom_upsetplot <- renderPlot({
    my_custom_upset_plot()
  })
  
  
  ## ENRICHER kegg table and plot
  output$gso_keggscape_table <- DT::renderDT({
    #shinyjs::hide("text_KEGGscape_ORA")
    #shinyjs::show("ora_keggscape_table_box")
    if( rv$map != "" && length(rv$map) > 0 && !is.na(rv$map )){  
      ko_nodes %>% 
        filter(map == rv$map) %>%
        #filter(map == "00630") %>%
        separate_rows(nodes, sep = ",") %>%
        left_join(ko_symbols, by = c("nodes" = "V1")) %>%
        inner_join(ko_annots, by=c("nodes" = "KO")) %>%
        filter(Name %in% representatives$Name) %>%
        mutate(Name = paste0(Name, " (KFK:", koala_threshold, ", em:", em, ")" )) %>%
        #filter(nodes %in% names(rv$kegg_gso_vector)) %>%
        #filter(Name %in% rv$input_data) %>%
        group_by(nodes, Symbol, Node_Description) %>% 
        summarise(Transcripts = paste(unique(Name), collapse= "; ")) %>%
        ungroup() %>%
        rename(KEGG_ko = nodes) 
    }
  },
  filter="top",
  options = list(
    scrollX = TRUE, 
    scrollY = "400px",
    dom = '<"top" lpift>')
  )
  
  
  observeEvent(
    eventExpr = input$submit_kegg_gso_toggle, {
      rv$kegg_gso_col_low = ifelse(rv$kegg_gso_col_low == "red", "green", "red")
      rv$kegg_gso_col_high = ifelse(rv$kegg_gso_col_high == "red", "green", "red")
    })
  
  
  observeEvent(
    eventExpr = input$submit_kegg_gsea_toggle, {
      rv$kegg_gsea_col_low = ifelse(rv$kegg_gsea_col_low == "red", "green", "red")
      rv$kegg_gsea_col_high = ifelse(rv$kegg_gsea_col_high == "red", "green", "red")
    })
  
  output$KEGG_map_gso <- renderImage({
    outfile <- tempfile(fileext='.png')
    png(outfile, width=400, height=300)
    frame()
    dev.off()
    blank <- list(src = outfile,
                  contentType = 'image/png',
                  width = 400,
                  height = 300,
                  alt = " ")
    
    if(rv$okplot_KEGGmap_gso){
      
      
      if( rv$map != "" && length(rv$map) > 0 && !is.na(rv$map )){  
        tempFolder <- tempdir() # tempFolder = "temp";
        #tempFolder <- "~/OneDrive - Fondazione Istituto Italiano Tecnologia/B_CRO_projects/re_algae_c026/gsea_shinyR"
        # this is the previous implementation  
        # print("FIND ko_annots")
        # gso_kegg <- ko_annots %>%
        #   filter(Name %in% rv$input_data) %>%
        #   left_join(rv$msi_df %>% select(Name, logFC) %>% separate_rows(Name, sep = ",") ) %>%
        #   group_by(val) %>%
        #   summarise(logFC = max(logFC, na.rm = T)) %>%
        #   ungroup()
        # 
        # rv$kegg_gso_vector <- gso_kegg$logFC
        # names(rv$kegg_gso_vector) <- gso_kegg$val
        # print(paste("kegg_gso_vector",dim(rv$kegg_gso_vector)[1]))
        # 
        # 
        # kegg_gso_vector = rv$kegg_gso_vector[names(rv$kegg_gso_vector) %in% str_split(subset(ko_nodes %>% filter(map == rv$map))$nodes, pattern = "[,]")[[1]]]
        # 
        # max_range = ceiling(max(abs(min(kegg_gso_vector)), abs( max(kegg_gso_vector)) ))
        # ##
        #
        kegg_gso_vector <- ko_nodes %>% 
          filter(map == rv$map) %>% 
          #filter(map == "05221") %>% 
          select(nodes) %>% 
          separate_rows(nodes) %>%
          distinct(nodes) %>%
          left_join(ko_annots, by = c("nodes" = "KO")) %>%
          left_join(rv$msi_df %>% select(Name, logFC) %>% separate_rows(Name, sep = ",") ) %>%
          #left_join(msi_df %>% select(Name, logFC) %>% separate_rows(Name, sep = ",") ) %>%
          group_by(nodes) %>%
          summarise(logFC = max(logFC, na.rm = T)) %>%
          ungroup() %>%
          mutate(logFC = ifelse(is.finite(logFC), logFC, NA)) %>%
          left_join(koala_thresholds , by = c("nodes" = "KO")) 
        
        kegg_gso_vector$koala_threshold = kegg_gso_vector$koala_threshold / 100
        print(paste("length kegg_gso_vector", length(kegg_gso_vector)))
        print(paste("input$kegg_gso_values", input$kegg_gso_values))
        #select the values to overlay in the kegg pathway
        if (input$kegg_gso_values == "logfc"){
          gene_vector =  kegg_gso_vector$logFC
          names(gene_vector) = kegg_gso_vector$nodes
          max_range = ceiling(max(abs(min(kegg_gso_vector$logFC, na.rm = T)), abs( max(kegg_gso_vector$logFC, na.rm = T)) ))
          min_range = -max_range
        }else if (input$kegg_gso_values == "em"){
          gene_vector =  kegg_gso_vector$em
          names(gene_vector) = kegg_gso_vector$nodes
          max_range = 1
          min_range = 0
        }else if (input$kegg_gso_values == "kfk"){
          gene_vector =  kegg_gso_vector$koala_threshold
          names(gene_vector) = kegg_gso_vector$nodes
          max_range = 1
          min_range = 0
        }else if (input$kegg_gso_values == "all"){
          kegg_gso_vector$logFC = kegg_gso_vector$logFC / max(abs(kegg_gso_vector$logFC), na.rm = T)
          gene_vector =  cbind(matrix(kegg_gso_vector$logFC), kegg_gso_vector$koala_threshold, kegg_gso_vector$em)
          rownames(gene_vector) =  kegg_gso_vector$nodes
          
          max_range = 1
          min_range = -1
        }
        
        print(paste("gene_vector", length(gene_vector)))
        
        
        #print(gene_vector)
        
        
        outfile <- ifelse(input$kegg_gso_values == "all",
                          paste( tempFolder,"/ko", rv$map, ".", rv$IP, ".gso.multi.png",sep=""),
                          paste( tempFolder,"/ko", rv$map, ".", rv$IP, ".gso.png",sep=""))
        #outfile <- paste( tempFolder,"/ko", rv$map, ".", rv$IP, ".gso.png",sep="")
        wd = getwd()
        setwd(tempFolder)
        
        pvout = pathview(gene.data = gene_vector,
                         gene.idtype = "KEGG",
                         pathway.id = rv$map,
                         species = "ko",
                         bins = list(gene = 8),
                         #limit = list(gene = c(0,1)),
                         limit = list(gene = c(min_range, max_range)),
                         low = list(gene = rv$kegg_gso_col_low  ),
                         #low = list(gene = "red"  ),
                         mid = list(gene = "yellow"),
                         #high = list(gene =  "green"),
                         high = list(gene =  rv$kegg_gso_col_high),
                         na.col = "transparent",
                         kegg.dir = file.path(wd, "kegg_data"),
                         #out.suffix = suffix,
                         kegg.native = T,
                         same.layer = T,
                         sign.pos = "bottomleft",
                         new.signature = F,
                         #key.pos = "topleft",
                         out.suffix = paste0(rv$IP,".gso"))
        setwd(wd)
        shinyjs::show("keggscape_gso_conf")        
        
        print(outfile)
        return( list(src = outfile,
                     contentType = 'image/png',
                     width = "150%",
                     height = "150%",
                     alt = "This is alternate text") )
        
        
      }else{
        print("empty images KEGG")
        html("text_KEGGscape_gso", "<strong>No KEGG terms! Try again with different filters</strong>")
        return (blank)
        #empty_image()
      }
    }else{
      return (blank)
    }
  }, deleteFile = F)
  
  
  ### and the gsea kegg table and kegg scape
  # output$gsea_keggscape_table <- DT::renderDT({
  #   if( rv$map_gsea != "" && length(rv$map_gsea) > 0 && !is.na(rv$map_gsea )){  
  #     ko_nodes %>% 
  #       filter(map == rv$map_gsea) %>%
  #       separate_rows(nodes, sep = ",") %>%
  #       left_join(ko_symbols, by = c("nodes" = "V1")) %>%
  #       inner_join(ko_annots, by=c("nodes" = "KO")) %>%
  #       filter(nodes %in% names(rv$kegg_gsea_vector)) %>%
  #       mutate(Name = paste0(Name, " (KFK:", koala_threshold, ", em:", em, ")" )) %>%
  #       group_by(nodes, Symbol, Node_Description) %>% 
  #       summarise(Transcripts = paste(unique(Name), collapse= "; ")) %>%
  #       rename(KEGG_ko = nodes) 
  #   }
  # },
  # filter="top",
  # options = list(
  #   scrollX = TRUE, 
  #   scrollY = "400px",
  #   dom = '<"top" lpift>')
  # )
  
  #the pathview GSEA
  output$KEGG_map_gsea <- renderImage({
    outfile <- tempfile(fileext='.png')
    png(outfile, width=400, height=300)
    frame()
    dev.off()
    blank <- list(src = outfile,
                  contentType = 'image/png',
                  width = 400,
                  height = 300,
                  alt = " ")
    
    if(rv$okplot_KEGGmap_gsea){
      shinyjs::enable("run")
      print(paste0("output$KEGG_map_gsea  rv$map_gsea:",rv$map_gsea,"-",is.null(rv$map_gsea)))
      if( rv$map_gsea != "" && length(rv$map_gsea) > 0 && !is.na(rv$map_gsea )){  
        # hist(rnorm(100, 4, 1),breaks = 50, main = rv$run_counter)
        
        tempFolder <- tempdir() # tempFolder = "temp";
        #pathID = rv$map_gsea
        #myIP = rv$IP
        #kegg_gsea_vector = rv$kegg_gsea_vector
        
        kegg_gsea_vector = rv$kegg_gsea_vector[names(rv$kegg_gsea_vector) %in% str_split(subset(ko_nodes %>% filter(map == rv$map_gsea))$nodes, pattern = "[,]")[[1]]]
        
        max_range = ceiling(max(abs(min(kegg_gsea_vector)), abs( max(kegg_gsea_vector)) ))
        
        outfile <- paste( tempFolder,"/ko", rv$map_gsea, ".", rv$IP, ".gsea.png",sep="")
        setwd(tempFolder)
        
        pvout = pathview(gene.data = kegg_gsea_vector,
                         gene.idtype = "KEGG",
                         pathway.id = rv$map_gsea,
                         species = "ko",
                         bins = list(gene = 8),
                         #limit = list(gene = c(0,1)),
                         limit = list(gene = c(-max_range, max_range)),
                         low = list(gene = rv$kegg_gsea_col_low  ),
                         mid = list(gene = "yellow"),
                         high = list(gene =  rv$kegg_gsea_col_high),
                         na.col = "transparent",
                         kegg.dir = file.path(wd, "kegg_data"),
                         #out.suffix = suffix,
                         kegg.native = T,
                         same.layer = T,
                         sign.pos = "bottomleft",
                         new.signature = F,
                         #key.pos = "topleft",
                         out.suffix = paste0(rv$IP,".gsea"))
        setwd(wd)
        shinyjs::show("keggscape_gsea_conf")
        
        print(outfile)
        return( list(src = outfile,
                     contentType = 'image/png',
                     width = "150%",
                     height = "150%",
                     alt = "This is alternate text") )
        
        
      }else{
        print("empty images KEGG")
        html("text_KEGGscape_GSEA", "<strong>No enriched KEGG terms! Try again with different filters</strong>")
        return (blank)
        #empty_image()
      }
    }else{
      return (blank)
    }
  }, deleteFile = T)
  
}