
library(RColorBrewer)
library(purrr)

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
library(readr)
library(DT)
#library(plotly)

library(clusterProfiler)
library(enrichplot)
library(pathview)
library(ggupset )
library(paws)
#
library("edgeR")
library(ggrepel) 
library(statmod)
library(pheatmap)
library("tximport")
#library(RUVSeq)

library(ngsReports)


if (!interactive()) sink(stderr(), type = "output")
wd = "/srv/shiny-server/gsea_shinyR_v1/"
version = "app.v1.4.3"

#print(version)
#setwd(wd)

source(file = file.path(version, "server_load_files.R"),
         local = TRUE,
         encoding = "UTF-8")

server <- function(input, output, session) {
  
  source(file = file.path(version, "server_app.R"),
         local = TRUE,
         encoding = "UTF-8")
}


