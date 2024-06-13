output$data_overview_table <- DT::renderDT({
  colnames_to_show = c(input$overview_table_Default_values,  
                       input$annotatations_ips,
                       input$annotatations_em,
                       input$annotatations_kfk,
                       input$annotatations_clusters,
                       input$annotatations_b2g,
                       input$annotatations_others
                       )
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
