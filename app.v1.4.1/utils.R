rename_msi_df_columns = function(colnames_all, old_name, new_name){
  colnames_all[colnames_all %in% c(old_name)]= new_name
  return(colnames_all)
}

get_filtered_msi_df <- function(df, condition, fdr_filter, logfc_filter){
  return(df %>%
           gather(k, v , -Name) %>% 
           separate(k, c("cond", "value_type"), sep = "[.]") %>% 
           spread(value_type, v) %>% 
           filter(cond == condition) %>%
           filter(abs(logFC) >= logfc_filter & FDR <= fdr_filter)
  )
}


# get_fdr_threshold = function(input_fdr){
#   fdr_filter <<- ifelse(input$fdr_by5 == T,  10^(-input_fdr) *5,  10^(-input_fdr) )
#   output$fdr_value <- renderText({ fdr_filter })
#   return(fdr_filter)
# }




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