# this is loaded either from the button listener when tcpm default custom is selected or 
# from the server_combine_galaxy when the data has been read.
withProgress(message = 'Data preprocessing', value = 0.4, {
  suppressWarnings(rowMeans <- apply(data,1, function(x) mean(as.numeric(x),na.rm=T)))
  
  
  c <- data[!duplicated(data$Geneid),]
  rownames(c) <- c$Geneid
  
  ########################################
  html("overview_text", "Processing Count Data...")
  incProgress(0.1, detail = paste("Processing Count Data..."))
  countData.all<-as.matrix(c[2:ncol(c)])
  
  rownames(countData.all)<-c$Geneid
  
  countData<-countData.all
  
  html("overview_text", "Filtering Count Data...")
  incProgress(0.1, detail = paste("Filtering Count Data..."))
  
  filter<-apply(countData,1,function(x) length(x[x>10])>=ncol(countData)/5)   # 20967 gene features, not always 5, should be (ncol/4) or /5
  filtered<-countData[filter, ]
  
  incProgress(0.1, detail = paste("Extracting conditions"))
  print("colnames(filtered)")
  print(colnames(filtered))
 
  
  # subgroup <- factor(gsub("\\.\\w+$","",colnames(filtered),perl=T))  # this is the subgroup
  # subgroup_class <- factor(gsub("\\w\\.\\w+$","",colnames(filtered),perl=T))  # this is the subgroup without numbering
  # treatment <- factor(gsub("^\\w+\\.","",colnames(filtered),perl=T))  # this is the treatment
  sst = get_group_subgroup_treatment(colnames(filtered))
  subgroup = sst[[1]]
  subgroup_class = sst[[2]]
  treatment = sst[[3]]
  print(subgroup_class)
  print(subgroup)
  print(treatment)
  
  
  choices = get_choices_for_contrasts(colnames(filtered), subgroup, subgroup_class, treatment)
 
  print(choices)
  
  rv$conditions = choices$contrasts
  updateSelectInput(session, "choose_conditions", choices = choices$contrasts, selected = choices$contrasts[1])
  
  
  html("overview_text", "Extracting conditions : DONE")
  setProgress(detail = paste("Extracting conditions : DONE"))
  
  shinyjs::show("menu_de_box")
  rv$filtered = filtered
  print("eventReactive submit_norm_and_de END")
  setProgress(1, detail = paste("Data preprocessing: DONE"))
})