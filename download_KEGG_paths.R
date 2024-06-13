library(pathview)

ko_Pathway_terms <- as.data.frame(read.table("Ka_KEGGpathways.txt", header = F, sep = "\t", stringsAsFactors = F, 
                                             colClasses = c("character", "character")))

for(kegg in ko_Pathway_terms$V1){
  print(kegg)
  if(file.exists(paste0("kegg_data/ko", kegg, ".png"))){
    next
  }
  download.kegg(kegg, species = "ko", kegg.dir = "kegg_data/", file.type =c("xml", "png"))
  # pathview(gene.data = c("00062", "01110"),
  #          gene.idtype = "KEGG",
  #          pathway.id = kegg,
  #          species = "ko",
  #          bins = list(gene = 20),
  #          limit = list(gene = c(0,1)),
  #          low = list(gene = "red"),
  #          mid = list(gene = "yellow"),
  #          high = list(gene = "green"),
  #          na.col = "transparent",
  #          kegg.dir = "kegg_data",
  #          #out.suffix = suffix,
  #          kegg.native = T,
  #          same.layer = T,
  #          sign.pos = "bottomleft",
  #          new.signature = F
  # )
  file.remove(paste0("ko", kegg , ".pathview.png"))
  Sys.sleep(2)
}
