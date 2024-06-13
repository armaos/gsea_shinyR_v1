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