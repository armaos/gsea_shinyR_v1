#fileDir <- system.file("fastqc/", package = "ngsReports")
print("Starting server FastQC")
fastqc_files <<- data.frame(files = list.files("fastqc", pattern = "fastqc.zip$", full.names = TRUE)) %>%
  mutate(Runacc = basename(files)) %>%
  separate(Runacc, c("Runacc_num"), sep = "_fastqc", remove = F) %>%
  separate(Runacc, c("Runacc"), sep = "_") %>%
  filter(Runacc %in% clean_RunAcc(rv$galaxy_outputs$file)) 

  print(fastqc_files)
  

if(length(fastqc_files$files) > 0)
{
  updateSelectInput(session, "fastqc_input", choices = fastqc_files$Runacc_num, selected = NULL)
  fdl <- FastqcDataList(fastqc_files$files)
  reads <- readTotals(fdl)
  rv$show_menu_fastqc = TRUE
}else{
  showNotification(paste("No available FastQC data"), duration = 10, type = "error")
  rv$show_menu_fastqc = FALSE
}



#updateSelectInput(session, "fastqc_input", choices = c("illumina", "iontorrent", "SRR1553606_1"), selected = NULL)

#maybe as table
# reads %>%
#   pander::pander(
#     big.mark = ",",
#     caption = "Read totals from R1 libraries", 
#     justify = "lr"
#   )


#
output$fastqc_summary = renderPlot({plotSummary(fdl)})

output$fastqc_readstotal = renderPlot({
  plotReadTotals(fdl) +
    theme(
      legend.position = c(1, 1), 
      legend.justification = c(1, 1),
      legend.background = element_rect(colour = "black")
    )
})


observeEvent(input$fastqc_input, {
  if(!is.null(input$fastqc_input)) {
    rv$fastqc_run_acc = subset(fastqc_files %>%
                                 mutate(r = 1:n() ) %>%
                                 filter(Runacc_num == input$fastqc_input))$r
  }else{rv$fastqc_run_acc = 0}
  
})



output$fastqc_perBaseQ = renderPlot({
  if(rv$fastqc_run_acc != 0){
    plotBaseQuals(fdl[[rv$fastqc_run_acc]])
  }else{
    empty_image()
  }
         
    
})
output$fastqc_fastqc_perBaseQ_all = renderPlot({plotBaseQuals(fdl)})

###
output$fastqc_SeqQ = renderPlot({
  if(rv$fastqc_run_acc != 0){
    plotSeqQuals(fdl[[rv$fastqc_run_acc]])
  }else{
    empty_image()
  }
})

output$fastqc_SeqQ_all_heat = renderPlot({plotSeqQuals(fdl)})
output$fastqc_SeqQ_all_line = renderPlot({plotSeqQuals(fdl, plotType = "line")})

####
output$fastqc_SeqCo = renderPlot({
  if(rv$fastqc_run_acc != 0){
    plotSeqContent(fdl[[rv$fastqc_run_acc]])
  }else{
    empty_image()
  }
  
})
output$fastqc_SeqCo_all_heat = renderPlot({plotSeqContent(fdl)})
output$fastqc_SeqCo_all_line = renderPlot({plotSeqContent(fdl, plotType = "line", nc = 1) })

####
output$fastqc_AdaptCo = renderPlot({
  if(rv$fastqc_run_acc != 0){
    plotAdapterContent(fdl[[rv$fastqc_run_acc]])
  }else{
    empty_image()
  }
 
})
output$fastqc_AdaptCo_all_heat = renderPlot({plotAdapterContent(fdl)})
output$fastqc_AdaptCo_all_line = renderPlot({plotAdapterContent(fdl, plotType = "line") })

###
output$fastqc_DupL  = renderPlot({
  if(rv$fastqc_run_acc != 0){
    plotDupLevels(fdl[[rv$fastqc_run_acc]])
  }else{
    empty_image()
  }
})
output$fastqc_DupL_all_heat  = renderPlot({plotDupLevels(fdl)})

###
output$fastqc_GC = renderPlot({
  if(rv$fastqc_run_acc != 0){
    plotGcContent(fdl[[rv$fastqc_run_acc]])
  }else{
    empty_image()
  }
  
})


###
output$fastqc_GC_all_heat = renderPlot({
  plotGcContent(fdl, species = "Athaliana", gcType = "Transcriptome")
})
output$fastqc_GC_all_line = renderPlot({
  plotGcContent(fdl, plotType = "line", species = "Athaliana", gcType = "Transcriptome")
})

####
output$fastqc_OverRep = renderPlot({
  if(rv$fastqc_run_acc != 0){
    plotOverrep(fdl[[rv$fastqc_run_acc]])
  }else{
    empty_image()
  }
 
})
output$fastqc_OverRep_all  = renderPlot({plotOverrep(fdl)})

