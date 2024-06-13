

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