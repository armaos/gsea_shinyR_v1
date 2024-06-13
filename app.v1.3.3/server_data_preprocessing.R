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